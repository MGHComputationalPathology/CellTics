"""
(c) MGH Center for Integrated Diagnostics

MGH CID VarGroup
Merge variant calls within a certain distance to one call
Started by Ryan Schmidt in R - https://bitbucket.org/mghcid/cidtools/branch/rs_variant_aggregator
Pythonized by Allison MacLeay
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
import click
import vcf
from urllib2 import urlopen, URLError
from xml.parsers.expat import ExpatError
import xmltodict
import pysam
import numpy as np
from celltics.lib import CIGAR_CODES


MAX_GROUPED = 200


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
        reads_existence, reads_coverage = [], []
        for read in sam_handle.fetch(self.get_chr(append_chr), self.pos, self.end):
            if read.is_duplicate:
                continue
            existence, coverage = self._build_existence_array(read)
            reads_existence.append(existence)
            reads_coverage.append(coverage)
        if len(reads_existence) < 1:
            print('Warning: No reads found at {}:{}-{}'.format(self.chrom, self.pos, self.end))

        self._build_existence_matrix(reads_existence, reads_coverage)

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

    @staticmethod
    def _get_indel_from_cigar(cigar, ignore_softclip=False):
        """
        Get a dictionary of positions and cigar codes
        :param cigar:
        :param ignore_softclip: does not add softclip bases to position.  Needed for SNV start position
        :return: {pos: cigar_code, pos: cigar_code}
        """
        pos = 0
        indels = {}
        for indel, num in cigar:
            if CIGAR_CODES[indel] in 'DI':
                indels[pos - 1] = (CIGAR_CODES[indel], num)
            if CIGAR_CODES[indel] != 'S' or not ignore_softclip:
                pos += num
        return indels

    def _get_start(self, variant, reference_start, cigar, ignore_softclip=False):
        """
        Return the read start from reference start position and cigar string
        :param variant:
        :param reference_start:
        :param cigar:
        :param ignore_softclip: does not add softclip bases to position.  Needed for SNV start position
        :return:
        """
        indels = self._get_indel_from_cigar(cigar, ignore_softclip)
        start = variant.POS - reference_start - 1
        for pos in indels:
            if pos > start:
                break
            if indels[pos][0] == 'I':
                start += indels[pos][1]
            elif indels[pos][0] == 'D':
                start -= indels[pos][1]
        return start

    def _build_existence_array(self, read):
        """
        Get array of 0 and 1 signifying whether the nth variant was observed in this read
        :param read: PySam read object
        :param vargroup: VariantGroup object
        :return: array of 0 and 1 [0,0,0,1]
        """
        existence = np.zeros(len(self.variant_list))
        coverage = np.zeros(len(self.variant_list))
        # pylint: disable=invalid-name
        for vn, variant in enumerate(self.variant_list):
            if variant.POS < read.reference_start or variant.POS > read.aend:
                continue
            coverage[vn] = 1
            # SNPs
            if variant.var_type == 'snp':
                # get start position using the cigar string to find the offset
                variant_start = self._get_start(variant, read.reference_start, read.cigar, ignore_softclip=True)
                # If the base matches the alternate read add it to the existence array
                if read.query[variant_start: variant.end - variant.POS + variant_start + 1] == variant.ALT[0].sequence:
                    existence[vn] = 1
            # Insertions/Deletions
            elif variant.is_indel:
                indels = self._get_indel_from_cigar(read.cigar)
                iloc = variant.POS - read.reference_start + read.query_alignment_start - 1
                # If the insertion/deletion exist in the cigar string add it to the existence array
                if variant.is_deletion and iloc in indels and indels[iloc][0] == 'D':  # Deletions
                    existence[vn] = 1
                elif not variant.is_deletion and iloc in indels and indels[iloc][0] == 'I':  # Insertions
                    if variant.ALT[0] == read.seq[iloc:iloc + 1 + indels[iloc][1]]:
                        existence[vn] = 1
            else:
                print('Warning: Unknown type found: {}'.format(variant.var_type))
        return existence, coverage

    def _build_existence_matrix(self, reads_existence, reads_coverage):
        self.exists = False
        if len(reads_existence) < 1:
            pass
        for existence, coverage in zip(reads_existence, reads_coverage):  # for each read
            # pylint: disable=consider-using-enumerate
            for i in range(len(existence)):  # for each variant
                for j in range(i):  # create a relationship matrix
                    if existence[i] == 1 and existence[j] == 1:  # if both variants exist
                        self.existence_array[i][j] += 1
                        self.exists = True
                    if existence[i] == 1 and existence[j] == 0:  # if a not b
                        self.a_not_b_array[i][j] += 1
                    if existence[i] == 0 and existence[j] == 1:  # if b not a
                        self.b_not_a_array[i][j] += 1
                    if coverage[i] == 1 and coverage[j] == 1:  # if there was coverage for both variant locations
                        self.coverage_array[i][j] += 1

    def _get_existence_frequency(self):
        """ Calculate Frequency of 2 variants occurring together """
        with np.errstate(invalid='ignore'):
            return np.nan_to_num(self.existence_array / self.coverage_array * 100)


def check_for_chr(sam):
    """ Check sam file to see if 'chr' needs to be prepended to chromosome """
    if 'chr' in sam.references[0]:
        return True
    return False


def inspect_bam(bam_file, variant_dict, threshold, min_reads, filter_type='min_pagb'):
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
        if filter_type in ['min_pagb', 'max_pagb']:  # Probability of A given B.
            vargroup.set_filter_fq_pagb(threshold, 'max' in filter_type)
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


def get_reverse_dict(variant_groups):
    """
    :param variant_groups:
    :return: Reversed dictionary from variant to the group they were in
    """
    rev_dict = {}
    for group in variant_groups:
        for variant in variant_groups[group]:
            variant_id = get_id(variant)
            if variant_id not in rev_dict:
                rev_dict[variant_id] = []
            rev_dict[variant_id].append(group)
    return rev_dict


def merge(variant_dict, unmerged, output_file, reader, ref_seq, mode='append'):
    """
    Return merged variant calls
    """
    # pylint: disable=protected-access
    reader.infos['GROUP_ID'] = vcf.parser._Info('GROUP_ID', 0, 'String', 'Group ID given by VarGrouper', '', '')
    # pylint: disable=protected-access
    reader.infos['IN_GROUP'] = vcf.parser._Info('IN_GROUP', 0, 'String', 'VarGrouper group ID this variant was put in',
                                                '', '')
    writer = vcf.Writer(open(output_file, 'w+'), template=reader)
    if mode != 'merged_only':
        var_map = get_reverse_dict(variant_dict)
        for record in unmerged:
            record_id = get_id(record)
            rec_status = var_map.get(record_id, '')
            if mode != 'intersect' or rec_status == '':
                record = append_group_info(record, rec_status)
                writer.write_record(record)
    if ref_seq:
        fdict = SeqIO.to_dict(SeqIO.parse(ref_seq, 'fasta'))
    else:
        fdict = None
    for var in variant_dict:
        record = merge_records(variant_dict[var], var, seq_dict=fdict)
        writer.write_record(record)
    writer.close()


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
    return vcf.model._Record(chrom, start, variants[0].ID, ref, [alt_obj], agg_qual / nvariants, '', info,
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
    request = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr{}:{},{}'.format(chrom, start, end)
    try:
        dna = xmltodict.parse(urlopen(request).read())['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n', '')
    except (URLError, ExpatError) as e:
        print('Could not open UCSC url.  Please check your internet connection.\n{}\n{}'.format(request, e.message))
        dna = "n" * (start - end)
    return dna


def parse_vcf(reader, merge_distance):
    """ Read through a sorted VCF file and aggregate candidates for variant merging """
    empty_record = vcf.model._Record(None, 0, None, '', '', '', '', '', '', '')  # pylint: disable=protected-access
    last_record = empty_record
    vars_to_group = {}
    not_aggregated = []
    current_group = None
    for record in reader:
        if last_record is not None and last_record.CHROM != record.CHROM:
            last_record = empty_record
        elif last_record.start > record.start:
            raise ValueError('VCF must be sorted.')
        # if same chromosome, within distance
        if record.CHROM == last_record.CHROM and record.POS - last_record.POS < merge_distance:
            # if overlapping, skip
            if not last_record.end < record.POS:
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


# pylint: disable=too-many-arguments
def main(input_file=None, output_file=None, bam_file=None, merge_distance=9, fq_threshold=0, min_reads=3,
         bam_filter_mode='pagb', write_mode='append', ref_seq=None):
    """The main function"""
    print('Grouping file: {}'.format(input_file))
    rejected_groups = []
    merge_distance = int(merge_distance)
    reader = vcf.Reader(open(input_file, 'r'))
    original_variants = list(reader)
    _, vars_to_group = parse_vcf(original_variants, merge_distance)
    if bam_file is not None:
        print("Inspecting bam file.")
        vars_to_group, rejected_groups = inspect_bam(bam_file, vars_to_group, fq_threshold, min_reads,
                                                     filter_type=bam_filter_mode)
    print("Merging records.")
    merge(vars_to_group, original_variants, output_file, reader, ref_seq, write_mode)
    print("Found %s variants to merge\nRejected groups: %s" % (len(vars_to_group.keys()), len(rejected_groups)))


@click.group(invoke_without_command=True)
@click.option('--input-file', '-i', required=True, type=click.Path(exists=True), help='Path to input vcf file')
@click.option('--output-file', '-o', required=True, help='Path to output vcf file')
@click.option('--bam-file', '-b', help='Path to bam file')
@click.option('--merge-distance', '-m', default=9, help='Find all variants within X bases. (default=9)')
@click.option('--fq-threshold', '-ft', default=0, help='Minimim frequency for grouping. (default=0)')
@click.option('--min-reads', '-r', default=3, help='Minimum supporting reads (default=3)')
@click.option('--bam-filter-mode', '-f', default='pab', type=click.Choice(['pab', 'min_pagb', 'max_pagb']),
              help='min_pagb - (default) Minimum probability of A given B.  Value is min of P(A|B) and P(B|A).\n'
                   'max_pagb - Maximum probability of A given B.  Value is max of P(A|B) and P(B|A).'
                   'pab - Probability of A and B\n')
@click.option('--write-mode', '-w', default='append', type=click.Choice(['append', 'intersect', 'merged_only']),
              help='append - adds merged calls to vcf\nintersect - removes single calls that were merged'
                   '\nmerged_only - outputs only merged calls')
@click.option('--ref-seq', '-rf', default=None, type=click.Path(exists=True),
              help='Path to reference genome file (required)')
# pylint: disable=too-many-arguments
def cli(input_file=None, output_file=None, bam_file=None, merge_distance=9, fq_threshold=0, min_reads=3,
        bam_filter_mode='min_pagb', write_mode='append', ref_seq=None):
    """
    Group variants in vcf
    """
    main(input_file=input_file, output_file=output_file, bam_file=bam_file, merge_distance=merge_distance,
         fq_threshold=fq_threshold, min_reads=min_reads, bam_filter_mode=bam_filter_mode, write_mode=write_mode,
         ref_seq=ref_seq)


if __name__ == '__main__':
    cli()
