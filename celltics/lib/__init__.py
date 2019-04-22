"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import

import numpy as np

__author__ = 'Allison MacLeay'

CIGAR_CODES = {0: 'M',
               1: 'I',
               2: 'D',
               3: '',
               4: 'S',
               5: 'H'}


def get_read_location(read, position, is_indel=False):
    """
    Return corrected location in read
    pysam read
    position - 1 indexed variant.POS from pysam
    is_indel: if True return the inserted sequence if within an insert
              if False ignore inserted sequence and return base at the genomic location
    :returns genomic location, cigar code, cigar size, distance from cigar break, subsequent deletion size if applicable
             ndel_sz is the size of the insertion or deletion at the next base position because indels are annotated
                     with the location of the base before the event, but this is usually at a M (match) positiondit

    """
    query_pos = position - read.pos - 1
    code = None
    curpos = 0
    shifted = 0
    ndel_sz = None
    size = 0
    for icount, (code, size) in enumerate(read.cigar):
        ndel_sz = None
        if CIGAR_CODES[code] in ['I', 'M']:
            curpos += size
        if CIGAR_CODES[code] == 'I':
            if not is_indel or (is_indel and curpos < query_pos):
                # ignore inserted sequence
                query_pos += size
            else:
                # return location inside inserted sequence
                shifted = query_pos - size
                break
        if CIGAR_CODES[code] == 'D':
            if curpos > query_pos - size:
                break
            else:
                query_pos -= size
        shifted = curpos - query_pos
        # shifted == 1 when base is the last before next cigar
        if shifted == 1 and len(read.cigar) > icount + 1 and CIGAR_CODES[read.cigar[icount + 1][0]] in ['D', 'I']:
            ndel_sz = read.cigar[icount + 1][1]
            if CIGAR_CODES[read.cigar[icount + 1][0]] == 'D':
                ndel_sz *= -1
        if curpos > query_pos:
            break
    return query_pos, code, size, shifted, ndel_sz


def get_pos(read, position, show_del=True):
    """Return the base at a read position from a genomic coordinate """
    query_pos, code, _, _, ndel_sz = get_read_location(read, position)
    base = read.query[query_pos:query_pos + 1]
    if show_del and (CIGAR_CODES[code] == 'D' or ndel_sz is not None):
        base = '-'
    return base


def has_variant(read, variant, strict=True):
    """ Does pileup read contain the variant?
        strict=True: guarantee correct size of indel (Needed for repeated sequence deletions or insertions)
    """
    query_pos, code, size, shift, ndel_sz = get_read_location(read, variant.POS, variant.is_indel)
    if read.query[query_pos:query_pos + len(variant.ALT[0])] == variant.ALT[0]:
        if variant.var_subtype == 'ins' and code is not None:
            if CIGAR_CODES[code] == 'M' and ndel_sz is not None:
                return not strict or len(variant.ALT[0]) == ndel_sz + 1
            #  add 1 to size because pyvcf returns the base before the insertion
            return CIGAR_CODES[code] == 'I' and (not strict or (len(variant.ALT[0]) == size + 1 and shift == 0))
        if variant.var_subtype == 'del' and code is not None:
            # CIGAR could be 'M' or 'D' since it uses the base before the deletion
            if CIGAR_CODES[code] == 'M' and ndel_sz is not None:
                return not strict or len(variant.REF) == -1 * ndel_sz + 1
            if CIGAR_CODES[code] == 'D':
                return not strict or (len(variant.REF) == size + 1) and shift == 0
            return False
        return code == 0  # False if not SNP
    return False


def bases_counter(bam, chrom_name, start, end, min_mapping_quality):
    """
    For a given position return a dictionary of reads supporting A,C,G,T or N (raw calls and >mmq).
    NOTE: This method supports SNV only (1bp size) not indels.
    :param chrom: String, chromosome id
    :param start: Integer , variant start genomic coordinate
    :param end: Integer , variant end genomic coordinate (should the same as start for SNV).
    :return: Dictionary of counts
    """
    nucleotides = 'ACGTN'
    start -= 1  # since pysam is 0-based index
    end -= 1  # since pysam is 0-based index
    is_snv = start == end - 1
    if is_snv:  # is SNV
        holder = 0
    else:  # indels , cannot handle them yet TODO
        holder = 'NA'

    counts = {'raw': dict(zip(list(nucleotides), [holder for _ in nucleotides])),  # All (unfiltered) reads
              'hq+': dict(zip(list(nucleotides), [holder for _ in nucleotides])),  # high quality reads on plus
              'hq-': dict(zip(list(nucleotides), [holder for _ in nucleotides]))}  # high quality reads on minus

    if not is_snv:
        return counts  # return counts with NA values since we can't count bases for INDELs
    for column in bam.pileup(chrom_name, start, end):
        pos = column.pos
        for read in column.pileups:
            if pos != start:
                continue
            if read.query_position is None:  # avoid reads without query_alignment_sequence
                continue
            current_read_base = read.query_sequence[read.query_position]
            counts['raw'][current_read_base] += 1
            if not read.is_duplicate and read.mapq >= min_mapping_quality:
                if read.is_reverse:
                    counts['hq-'][current_read_base] += 1
                else:
                    counts['hq+'][current_read_base] += 1
    return counts


def get_indel_from_cigar(cigar, ignore_softclip=False):
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


def chunks(num, array):
    """ split array for multiprocessing.
        ex.  > chunks(3, [1, 2, 3])
             [[1], [2], [3]]
             > chunks(3, [1, 2, 3, 4])
             [[1], [2], [3]]
             > chunks(2, [1, 2, 3, 4])
             [[1, 2], [3, 4]]
    """
    return [arri for arri in np.array_split(array, num) if len(arri) > 0]
