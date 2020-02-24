#!/usr/bin/env python3
"""
(c) MGH Center for Integrated Diagnostics

MGH CID VarGroup
Merge variant calls within a certain distance to one call
Started by Ryan Schmidt in R - https://bitbucket.org/mghcid/cidtools/branch/rs_variant_aggregator
Pythonized by Allison MacLeay

"""

import pandas as pd
from Bio import SeqIO
import click
import vcf
# from urllib2 import urlopen, URLError  # python2
from urllib.request import urlopen
from urllib.error import URLError
from xml.parsers.expat import ExpatError
import xmltodict
from collections import OrderedDict, defaultdict
import pysam
import numpy as np
from celltics.lib import chunks, get_indel_from_cigar
import time
import multiprocessing as mp


MAX_GROUPED = 200


class PickleMe(object):
    """ Errors with the pickle library resulted in the need to remove the CallData objects from vcf.model.Record"""
    def __init__(self, records, var_dict):
        """ Clean data """
        self.records = [self.flatten(rec) for rec in records]
        self.var_dict = self.get_reverse_dict(var_dict)

    @staticmethod
    def get_reverse_dict(variant_groups):
        """
        :param variant_groups:
        :return: Reversed dictionary from variant to the group they were in
        """
        rev_dict = defaultdict(list)
        for group in variant_groups:
            for variant in variant_groups[group]:
                variant_id = get_id(variant)
                rev_dict[variant_id].append(group)
        return rev_dict

    def flatten(self, variant):
        """ remove CallData from variant group lists
            causes sample data string (ex.  ./.:0/1:203,16:8.0:219:0.073:.:.:.:.:.:.:.:.) omission from grouped variants
        """
        for sample in variant.samples:
            sample.data = None
        return variant

    def get_fat(self):
        """ convert data flattened for the pickle library to original version (or close to it) """
        if len(self.records) == 0 or self.records[0].FORMAT is None:
            return self.records, self.var_dict
        fields = [fld for fld in self.records[0].FORMAT.split(':')]
        calldata = vcf.model.make_calldata_tuple(fields)
        args = ['./.']
        args.extend([None for _ in fields[1:]])
        cd_obj = calldata(*args)
        for rec in self.records:
            rec.samples[0].data = cd_obj
        return self.records, self.var_dict


class VariantGroup(object):
    """
    A group of PyVcf Record objects
    """
    def __init__(self, var_id, var_list):
        arlen = len(var_list)
        self.variant_id = var_id
        self.variant_list = var_list
        self.valid = None
        self.chrom = var_list[0].CHROM
        self.pos = min([var.start for var in var_list])
        self.end = max([var.end for var in var_list])
        self.exists = False
        self.coverage_array = np.zeros((arlen, arlen))
        self.existence_array = np.zeros((arlen, arlen))
        self.a_not_b_array = np.zeros((arlen, arlen))
        self.b_not_a_array = np.zeros((arlen, arlen))
        self.filter = np.zeros((arlen, arlen)) == 0

    def get_chr(self, append_chr=False):
        """\
        Gets the formatted chromosome name
        :param append_chr: whether to prefix with 'chr'
        :return: string of chromosome name
        """
        if append_chr:
            return 'chr{}'.format(self.chrom)
        return self.chrom

    # pylint: disable=too-many-nested-blocks
    def _get_subgroups(self):
        """ Returns array of unrelated groups using self.filter boolean matrix """
        groups = []  # array of arrays
        for i in range(self.filter.shape[0]):
            for j in range(i):
                if self.filter[i][j]:
                    if len(groups) < 1:
                        groups.append([j, i])
                        continue
                    found = False
                    for group_i, _ in enumerate(groups):
                        if i in groups[group_i]:
                            if j not in groups[group_i]:
                                groups[group_i].append(j)
                            found = True
                        elif j in groups[group_i]:
                            if i not in groups[group_i]:
                                groups[group_i].append(i)
                            found = True
                    if not found:
                        groups.append([i, j])
        return groups

    def split_and_trim(self):
        """
        Operation uses self.filter boolean matrix for variant relationship status
        Trim out variants that do not co-occur with other variants
        Split into subgroups if this group contains unrelated groups
        :return variant_subgroups (array), rejected_variants (array)
        """
        indeces_grouped = []
        rejected_variants = []
        groups = self._get_subgroups()
        variant_groups = []
        for group in groups:
            variant_list = []
            for i in group:
                indeces_grouped.append(i)
                variant_list.append(self.variant_list[i])
            var_id = '{}_{}'.format(variant_list[0].CHROM, min([var.start for var in variant_list]))
            variant_groups.append(VariantGroup(var_id, variant_list))
        for v, variant in enumerate(self.variant_list):
            if v not in indeces_grouped:
                rejected_variants.append(variant)
        return variant_groups, rejected_variants

    def unsplit(self, variant_groups):
        """ Merge variant groups into this group """
        for vargroup in variant_groups:
            self.variant_list.extend(vargroup.variant_list)
        self.pos = min([var.start for var in self.variant_list])
        self.end = max([var.end for var in self.variant_list])

    def parse_sam(self, sam_handle, append_chr=False):
        """ Set counts of dual coverage and variant presence """
        vargroup_reads = np.asarray([read
                                     for read in sam_handle.fetch(self.get_chr(append_chr), self.pos, self.end)
                                     if not read.is_duplicate])

        # Convert some key read information into a dataframe to speed up filtering
        read_df = pd.DataFrame(columns=['rn', 'start', 'end', 'read', 'indels'],
                               data=[(rn, read.reference_start, read.aend, read, get_indel_from_cigar(read.cigar))
                                     for rn, read in enumerate(vargroup_reads)])

        reads_coverage = np.zeros((len(vargroup_reads), len(self.variant_list)))
        reads_existence = np.zeros((len(vargroup_reads), len(self.variant_list)))

        if len(vargroup_reads) == 0:
            print('Warning: No reads found at {}:{}-{}'.format(self.chrom, self.pos, self.end))
            return self._build_existence_matrix(reads_existence, reads_coverage)

        # pylint: disable=invalid-name
        for vn, variant in enumerate(self.variant_list):
            # Cache variant properties: those lookups are expensive in PySam
            var_type = variant.var_type
            is_indel = variant.is_indel
            is_deletion = variant.is_deletion

            read_overlap_mask = (read_df['start'] <= variant.POS) & (read_df['end'] >= variant.POS)

            # Coverage is easy: all reads which overlap this variant get a coverage of 1
            reads_coverage[read_overlap_mask, vn] = 1

            # SNPs
            if var_type == 'snp':
                # for rn, read, indels in itertools.izip(read_df[read_overlap_mask]['rn'],  # python2
                for rn, read, indels in zip(read_df[read_overlap_mask]['rn'], read_df[read_overlap_mask]['read'],
                                            read_df[read_overlap_mask]['indels']):
                    # get start position using the cigar string to find the offset
                    variant_start = self._get_start(variant, read.reference_start, read.cigar, ignore_softclip=True)
                    # If the base matches the alternate read add it to the existence array
                    read_alt = read.query[variant_start: variant.end - variant.POS + variant_start + 1]
                    if read_alt == variant.ALT[0].sequence:
                        reads_existence[rn, vn] = 1

            # Insertions/Deletions
            elif is_indel:
                # for rn, read, indels in itertools.izip(read_df[read_overlap_mask]['rn'],  # python2
                for rn, read, indels in zip(read_df[read_overlap_mask]['rn'], read_df[read_overlap_mask]['read'],
                                            read_df[read_overlap_mask]['indels']):
                    iloc = self._get_indel_pos(variant.POS, read)
                    # If the insertion/deletion exist in the cigar string add it to the existence array
                    if is_deletion and iloc in indels and indels[iloc][0] == 'D':  # Deletions
                        reads_existence[rn, vn] = 1
                    elif not is_deletion and iloc in indels and indels[iloc][0] == 'I':  # Insertions
                        if variant.ALT[0] == read.seq[iloc:iloc + 1 + indels[iloc][1]]:
                            reads_existence[rn, vn] = 1
            else:
                print('Warning: Unknown type found: {}'.format(variant.var_type))

        return self._build_existence_matrix(reads_existence, reads_coverage)

    def set_filter_fq_pab(self, threshold):
        """ Set a filter on the frequency of observing 2 variants together
            Probability of A and B -> P(AB)
        """
        frequency_table = self._get_existence_frequency()
        self.filter = frequency_table > threshold

    def set_filter_fq_pagb(self, threshold, take_max):
        """ Set a filter on the frequency of observing 2 variants together
            Probability of A given B -> P(A|B)
        """
        with np.errstate(invalid='ignore'):
            prob_a_given_b = np.nan_to_num(self.existence_array / (self.existence_array + self.b_not_a_array))
            prob_b_given_a = np.nan_to_num(self.existence_array / (self.existence_array + self.a_not_b_array))
            self.filter = self._matrix_wise_max_or_min(prob_a_given_b, prob_b_given_a, take_max) > threshold

    @staticmethod
    # pylint: disable=invalid-name
    def _matrix_wise_max_or_min(X, Y, is_max):
        arlen = len(X)
        matval = np.zeros((arlen, arlen))
        for i in range(arlen):
            for j in range(i):
                if is_max:
                    matval[i][j] = np.max([X[i][j], Y[i][j]])
                else:
                    matval[i][j] = np.min([X[i][j], Y[i][j]])
        return matval * 100

    def add_filter_min_reads(self, min_reads):
        """
        Add a filter for minimum occurrences of variants existing on same read
        :param min_reads: default 3
        :return: N/A
        """
        mr_filter = self.existence_array > min_reads
        self.filter &= mr_filter

    def reset_filter(self):
        """ Reset filter to matrix of True """
        arlen = len(self.variant_list)
        self.filter = np.zeros((arlen, arlen)) == 0

    def _get_start(self, variant, reference_start, cigar, ignore_softclip=False):
        """
        Return the read start from reference start position and cigar string
        :param variant:
        :param reference_start:
        :param cigar:
        :param ignore_softclip: does not add softclip bases to position.  Needed for SNV start position
        :return:
        """
        indels = get_indel_from_cigar(cigar, ignore_softclip)
        start = variant.POS - reference_start - 1
        # for pos, val in indels.iteritems(): # python2
        for pos, val in indels.items():
            if pos > start:
                break
            if val[0] == 'I':
                start += val[1]
            elif val[0] == 'D':
                start -= val[1]
        return start

    def _get_indel_pos(self, variant_pos, read):
        """ Get position of indel """
        hardclipped = 0 if read.cigartuples[0][0] != 5 else read.cigartuples[0][1]  # read location must be adjusted for
        # number of hardclipped bases represented in cigar but not in read_seq  https://www.biostars.org/p/119537/
        iloc = variant_pos - read.reference_start + read.query_alignment_start - 1 + hardclipped
        return iloc

    def _build_existence_matrix(self, reads_existence, reads_coverage):
        self.exists = False
        if reads_existence.shape[0] == 0:
            return

        n_variants = len(self.variant_list)
        self.coverage_array = np.zeros((n_variants, n_variants))
        self.existence_array = np.zeros((n_variants, n_variants))
        self.a_not_b_array = np.zeros((n_variants, n_variants))
        self.b_not_a_array = np.zeros((n_variants, n_variants))

        for i in range(n_variants):  # for each variant
            for j in range(i):  # create a relationship matrix
                covered = (reads_coverage[:, i] == 1) & (reads_coverage[:, j] == 1)
                self.coverage_array[i, j] = np.sum(covered)
                self.existence_array[i, j] = np.sum((reads_existence[:, i] == 1) & (reads_existence[:, j] == 1))
                self.a_not_b_array[i, j] = np.sum((reads_existence[:, i] == 1) & (reads_existence[:, j] == 0) & covered)
                self.b_not_a_array[i, j] = np.sum((reads_existence[:, i] == 0) & (reads_existence[:, j] == 1) & covered)

        self.exists = np.any(self.existence_array > 0)

    def _get_existence_frequency(self):
        """ Calculate Frequency of 2 variants occurring together """
        with np.errstate(invalid='ignore'):
            return np.nan_to_num(self.existence_array / self.coverage_array * 100)


def check_for_chr(sam):
    """ Check sam file to see if 'chr' needs to be prepended to chromosome """
    if 'chr' in sam.references[0]:
        return True
    return False


def inspect_bam(bam_file, variant_dict, threshold, min_reads, filter_type='pagb'):
    """
    Determine whether variant groups are observed on the same read
    :param bam_file:
    :param variant_dict:
    :return: {variant_groups}, [rejected_variant_groups]
    """
    sam = pysam.AlignmentFile(bam_file, 'rb')  # pylint: disable=no-member
    append_chr = check_for_chr(sam)
    valid_variants = {}
    rejected_variants = []
    for key in variant_dict:
        vargroup = VariantGroup(key, variant_dict[key])
        if len(vargroup.variant_list) > MAX_GROUPED:  # TODO: split into smaller groups
            print('Skipping variant group {} ({} of {})\nSize: {}  Chrom: {}  Start: {}  End: {}'
                  ''.format(key, variant_dict.keys().index(key), len(variant_dict.keys()), len(vargroup.variant_list),
                            vargroup.chrom, vargroup.pos, vargroup.end))
            rejected_variants.extend(vargroup.variant_list)
            continue
        vargroup.parse_sam(sam, append_chr)
        if filter_type == 'pab':  # Probability of A and B
            vargroup.set_filter_fq_pab(threshold)
        if filter_type in ['pagb', 'max_pagb']:  # Probability of A given B.
            vargroup.set_filter_fq_pagb(threshold, 'm' in filter_type)
        vargroup.add_filter_min_reads(min_reads)

        if vargroup.exists:
            split_vargroups, reject = vargroup.split_and_trim()
            if reject:
                # print('Trimming {} from group {}'.format(len(reject), key))
                rejected_variants.extend(reject)
            for k, s_vargroup in enumerate(split_vargroups):
                s_key = '{}_{}'.format(key, k)
                if not s_vargroup.variant_list:
                    pass
                    # print('No variant groups from {} passed the threshold of {}'.format(key, threshold))
                else:
                    valid_variants[s_key] = s_vargroup.variant_list
        else:
            rejected_variants.extend(vargroup.variant_list)

    return valid_variants, rejected_variants


def get_id(variant):
    """
    :param variant:
    :return: chrom_pos_ref_alt
    """
    return '{}_{}_{}_{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])


def merge_records(variants, group_id, seq_dict=None):
    """
    merge one list of variants
    """
    # pylint: disable=too-many-locals
    nvariants = len(variants)
    chrom = variants[0].CHROM
    start = min([variant.POS for variant in variants])
    info = variants[0].INFO
    end = max([variant.end for variant in variants])
    ref = get_reference_seq(chrom, start, end, seq_dict)
    alt = ref

    agg_qual, agg_alt_af, agg_alt_dp, agg_ref_dp = 0, 0, 0, 0
    shift = 0
    for variant in sorted(variants, key=lambda var: int(var.POS)):
        var_len = variant.end - variant.POS + 1
        var_alt = [sub.sequence for sub in variant.ALT]
        if variant.is_deletion:
            del_length = variant.end - variant.POS
            alt = alt[:shift + variant.POS - start] + "".join(var_alt) + \
                  alt[shift + variant.end - start + len(var_alt[0]):]
            shift -= del_length
        else:
            alt = alt[:shift + variant.POS - start] + "".join(var_alt) + alt[shift + variant.POS - start + var_len:]
            shift += len(variant.ALT[0]) - len(variant.REF)  # for insertions
        if variant.QUAL is not None:
            agg_qual += variant.QUAL
        agg_alt_af += variant.INFO.get('ALT_AF', 0) or 0
        agg_alt_dp += variant.INFO.get('ALT_DP', 0) or 0
        agg_ref_dp += variant.INFO.get('REF_DP', 0) or 0

    info['ALT_AF'] = agg_alt_af / len(variants)
    info['ALT_DP'] = agg_alt_dp / len(variants)
    info['REF_DP'] = agg_ref_dp / len(variants)
    info['PRIMARY_CALLER'] = ['vargrouper']
    info['GROUP_ID'] = group_id
    info['IN_GROUP'] = []
    alt_obj = vcf.model._Substitution(alt)  # pylint: disable=protected-access
    # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes, samples=None
    # pylint: disable=protected-access
    return vcf.model._Record(chrom, start, variants[0].ID, ref, [alt_obj], agg_qual / nvariants, None, info,
                             variants[0].FORMAT, variants[0]._sample_indexes, variants[0].samples)


def append_group_info(record, groups):
    """ Add new info fields """
    record.INFO['IN_GROUP'] = groups
    record.INFO['GROUP_ID'] = ''
    return record


def get_reference_seq(chrom, start, end, seq_dict=None):
    """
       Lookup in local full genome fasta file to return reference sequence
       :param chrom:
       :param start:
       :param end:
       :return: dna reference sequence
    """
    if not seq_dict:
        return get_reference_seq_ucsc(chrom, start, end)
    # ex. start = 1, end = 3, return [0,1,2]
    # because computer scientists are exclusive starting from 0 but biologists are inclusive starting from 1
    start = int(start) - 1
    try:
        dna = str(seq_dict[chrom][start:end].seq)
    except IndexError as e:
        raise Exception("Error: Could not find that sequence: %s" % str(e))
    except KeyError as e:
        print("No chromosome named: %s\nTrying UCSC..." % str(e))
        dna = get_reference_seq(chrom, start, end)
    return dna.lower()


def get_reference_seq_ucsc(chrom, start, end):
    """
    UCSC http request to return reference sequence
    :param chrom:
    :param start:
    :param end:
    :return: dna reference sequence
    """
    if chrom.startswith('chr'):
        chrom = chrom.replace('chr', '')
    request = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr{}:{},{}'.format(chrom, start, end)
    try:
        dna = xmltodict.parse(urlopen(request).read())['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n', '')
    except (URLError, ExpatError) as e:
        print('Could not open UCSC url.  Please check your internet connection.\n{}\n{}'.format(request, e.message))
        dna = "n" * (start - end)
    return dna


def parse_vcf(reader, merge_distance, skip_overlap):
    """ Read through a sorted VCF file and aggregate candidates for variant merging """
    empty_record = vcf.model._Record(None, 0, None, '', '', '', '', '', '', '')  # pylint: disable=protected-access
    last_record = empty_record
    vars_to_group = OrderedDict()
    not_aggregated = []
    current_group = None
    for record in reader:
        if last_record is not None and last_record.CHROM != record.CHROM:
            last_record = empty_record
        elif last_record.start > record.start:
            raise ValueError('VCF must be sorted.')
        # if same chromosome, within distance
        if record.CHROM == last_record.CHROM and record.POS - last_record.POS < merge_distance:
            # add if overlapping, because it will get detangled after.  Skip overlap is True when bam_file is None
            if not last_record.end < record.POS and skip_overlap:
                not_aggregated.append(record)
                last_record = record
                continue
            # Add to data structure
            if current_group is None:
                current_group = '_'.join([last_record.CHROM, str(last_record.POS)])
                vars_to_group[current_group] = [last_record, record]
            else:
                vars_to_group[current_group].append(record)
        else:
            current_group = None
            not_aggregated.append(record)
        last_record = record
    return not_aggregated, vars_to_group


def dict_chunks(udict, nchunks):
    """ Chunk dictionary """
    cbr = []  # candidate breaks (at chromosome splits)
    sizes = {}
    chroms = {}
    ukeys = sorted(udict.keys())
    last = None
    tsize = 0
    for i, k in enumerate(ukeys):
        tsize += 1
        chrom = k.split('_')[0]
        if chrom != last:
            cbr.append(i)
            sizes[last] = tsize  # TODO - balance split on size
            chroms[i] = chrom
            tsize = 0
        last = chrom
    chrom_chunks = chunks(nchunks, list(chroms.values()))
    cloc = {}
    for i, chrom_arr in enumerate(chrom_chunks):
        for chrom in chrom_arr:
            cloc[chrom] = i
    cdict = [OrderedDict() for _ in chrom_chunks]
    for dkey in udict:
        chrom = dkey.split('_')[0]
        cdict[cloc[chrom]][dkey] = udict[dkey]
    return cdict, cloc


def get_ref_seq_dict(ref_seq):
    """ get ref seq dict """
    return SeqIO.to_dict(SeqIO.parse(ref_seq, 'fasta')) if ref_seq else None


def get_merged_records(variant_dict, fdict):
    """
        Return merged variant calls
    """
    records = []
    for var in variant_dict:
        record = merge_records(variant_dict[var], var, seq_dict=fdict)
        records.append(record)
    return records


def bam_and_merge(bam_file, vars_to_group, fq_threshold, min_reads, filter_type, fdict):
    """ inspect bam and merge records """
    if bam_file is not None:
        print('Processing bam: {}'.format(bam_file))
        vars_to_group, _ = inspect_bam(bam_file, vars_to_group, fq_threshold, min_reads, filter_type)
    return PickleMe(get_merged_records(vars_to_group, fdict), vars_to_group)


def write_vcf(records, var_dict, output_file, reader, original_variants, write_mode):
    """
    Write merged variant calls to vcf
    """
    # pylint: disable=protected-access
    reader.infos['GROUP_ID'] = vcf.parser._Info('GROUP_ID', 0, 'String', 'Group ID given by VarGrouper', 'VarGrouper',
                                                '2.0')
    # pylint: disable=protected-access
    reader.infos['IN_GROUP'] = vcf.parser._Info('IN_GROUP', 0, 'String', 'VarGrouper group ID this variant was put in',
                                                'VarGrouper', '2.0')
    with open(output_file, 'w+') as ovcf:
        writer = vcf.Writer(ovcf, template=reader)
        if write_mode != 'merged_only':
            for record in original_variants:
                rec_id = get_id(record)
                rec_status = var_dict.get(rec_id, '')
                if write_mode != 'intersect' or rec_status == '':
                    append_group_info(record, rec_status)
                    writer.write_record(record)
        for record in records:
            writer.write_record(record)


def split_ref_seq(fdict, cloc, nthreads):
    """ Fix shared memory error with multiprocessing by splitting reference dictionary """
    if fdict is None:
        return [None for _ in range(nthreads)]
    ref_array = [{} for _ in range(nthreads)]
    for chrom in cloc:
        ref_array[cloc[chrom]][chrom] = fdict[chrom]
    return ref_array


# pylint: disable=too-many-locals
def bam_and_merge_multiprocess(bam_file, vars_to_group, fq_threshold, min_reads, bam_filter_mode, ref_seq, nthreads,
                               debug=False):
    """ Multiprocess inspect bam and merge step together """
    fdict = None
    if bam_file is not None:
        if ref_seq is not None:
            print('Opening reference sequence: {}'.format(ref_seq))
        fdict = get_ref_seq_dict(ref_seq)
    start = time.time()
    var_chunks, chrom_pos = dict_chunks(vars_to_group, nthreads)
    refs = split_ref_seq(fdict, chrom_pos, len(var_chunks))
    del fdict
    if not debug and nthreads > 1:
        pool = mp.Pool(processes=nthreads)
        res = [pool.apply_async(bam_and_merge, args=(bam_file, achunk, fq_threshold, min_reads, bam_filter_mode,
                                                     ref)) for achunk, ref in zip(var_chunks, refs)]
        pool.close()
    else:
        res = [bam_and_merge(bam_file, achunk, fq_threshold, min_reads, bam_filter_mode, ref)
               for achunk, ref in zip(var_chunks, refs)]
    records = []
    var_dict = {}
    for r in res:
        recs, var_dict_part = r.get().get_fat() if not debug and nthreads > 1 else r.get_fat()
        records.extend(recs)
        var_dict.update(var_dict_part)
    how_long = time.time() - start
    print("%d seconds to process with %d threads" % (how_long, nthreads))
    return records, var_dict


# pylint: disable=too-many-arguments,too-many-locals
def main(input_file=None, output_file=None, bam_file=None, merge_distance=9, fq_threshold=0, min_reads=3,
         bam_filter_mode='pagb', write_mode='append', ref_seq=None, threads=None, debug=False):
    """the main function"""
    nthreads = mp.cpu_count() if threads is None else min(int(threads), mp.cpu_count())
    print('Grouping file: {}'.format(input_file))
    merge_distance = int(merge_distance)
    with open(input_file, 'r') as openfile:
        reader = vcf.Reader(openfile)
        original_variants = list(reader)
        _, vars_to_group = parse_vcf(original_variants, merge_distance, bam_file is None)
        n_candidates = len(vars_to_group)
        records, var_dict = bam_and_merge_multiprocess(bam_file, vars_to_group, fq_threshold, min_reads,
                                                       bam_filter_mode, ref_seq, nthreads, debug)
        write_vcf(records, var_dict, output_file, reader, original_variants, write_mode)
    print("Found %s variants to merge\nRejected groups: %s" % (len(records), n_candidates - len(records)))


@click.group(invoke_without_command=True)
@click.option('--input-file', '-i', required=True, type=click.Path(exists=True),
              help='Path to input VarVetter input file (required)')
@click.option('--output-file', '-o', required=True, help='Path to output VarVetter input file (required)')
@click.option('--bam-file', '-b', help='Path to bam file')
@click.option('--merge-distance', '-m', default=9, help='Find all variants within X bases. (default=9)')
@click.option('--fq-threshold', '-ft', default=0, help='Minimim frequency for grouping. (default=0)')
@click.option('--min-reads', '-r', default=3, help='Minimum supporting reads (default=3)')
@click.option('--bam-filter-mode', '-f', default='pab', type=click.Choice(['pab', 'pagb', 'max_pagb']),
              help='pagb - (default) Minimum probability of A given B.  Value is min of P(A|B) and P(B|A).\n'
                   'max_pagb - Maximum probability of A given B.  Value is max of P(A|B) and P(B|A).'
                   'pab - Probability of A and B\n')
@click.option('--write-mode', '-w', default='append', type=click.Choice(['append', 'intersect', 'merged_only']),
              help='append - (default) adds merged calls to vcf\nintersect - removes single calls that were merged'
                   '\nmerged_only - outputs only merged calls')
@click.option('--ref-seq', '-rf', default=None, type=click.Path(exists=True),
              help='Path to reference genome file (required)')
@click.option('--threads', '-t', default=None, help='Number of threads')
@click.option('--debug', '-db', is_flag=True, default=False, help='Run in debug mode (No multiprocessing)')
# pylint: disable=too-many-arguments
def cli(input_file=None, output_file=None, bam_file=None, merge_distance=9, fq_threshold=0, min_reads=3,
        bam_filter_mode='min_pagb', write_mode='append', ref_seq=None, threads=None, debug=False):
    """
    :param input_file: VCF file with variants to merge
    :param merge_distance: merge variants within this distance
    :return: merged VCF file
    """
    main(input_file=input_file, output_file=output_file, bam_file=bam_file, merge_distance=merge_distance,
         fq_threshold=fq_threshold, min_reads=min_reads, bam_filter_mode=bam_filter_mode, write_mode=write_mode,
         ref_seq=ref_seq, threads=threads, debug=debug)


if __name__ == '__main__':
    cli()
elif not __name__.startswith('celltics'):
    print('WARNING:  {}: This module cannot be called using the cli library.\n'
          'Multiprocessing does not work because of dependency on the pickle library.  Classes created at the '
          'interpreter level cannot be pickled.\n'.format(__name__))
