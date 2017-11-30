# -*- coding: utf-8 -*-

"""\
(c) 2015-2016 MGH Center for Integrated Diagnostics

"""

from __future__ import print_function
from __future__ import unicode_literals

from nose.tools import assert_dict_equal, assert_true, assert_equal, assert_raises
from vcf.model import _Record
import numpy as np
import os
import filecmp
import celltics.tools.vargroup as vg
from pkg_resources import resource_filename
from Bio import SeqIO


def test_get_indels_from_cigar():
    """Validates that we correctly parse indels from cigar strings"""
    cigar = [(4, 21), (0, 100), (2, 1), (0, 30)]
    results = {True: {99: ('D', 1)},
               False: {120: ('D', 1)}}

    for value in [True, False]:
        # pylint:disable=protected-access
        assert_dict_equal(results[value], vg.VariantGroup._get_indel_from_cigar(cigar, ignore_softclip=value))


def test_split_and_trim():
    """Validates that we correctly remove variants which don't meet a threshold"""
    records = [_Record(1, 0, '1', '', '', '', '', '', '', ''),
               _Record(2, 0, '2', '', '', '', '', '', '', ''),
               _Record(3, 0, '3', '', '', '', '', '', '', ''),
               _Record(4, 0, '4', '', '', '', '', '', '', '')]
    groups = {'1': vg.VariantGroup('1', records),
              '2': vg.VariantGroup('1', records), }
    arlen = len(records)
    groups['1'].coverage_array = np.zeros((arlen, arlen))
    groups['1'].existence_array = np.zeros((arlen, arlen))
    groups['1'].coverage_array[-1][0] = 100.
    groups['1'].existence_array[-1][0] = 50.
    groups['1'].coverage_array[-2][1] = 1000.
    groups['1'].existence_array[-2][1] = 1.
    groups['2'].coverage_array = groups['1'].coverage_array.copy()
    groups['2'].existence_array = groups['1'].existence_array.copy()

    groups['1'].set_filter_fq_pab(25)
    groups['2'].set_filter_fq_pab(0)

    correct = {'1': ['1', '4'],
               '2': ['1', '2', '3', '4']}
    removed = {'1': ['2', '3'],
               '2': []}
    keep, reject = {}, {}
    keep['1'], reject['1'] = groups['1'].split_and_trim()
    keep['2'], reject['2'] = groups['2'].split_and_trim()
    keep['1'] = keep['1'][0]
    keep['2'][0].unsplit(keep['2'][1:])
    obs_removed = {'1': [vlist.ID for vlist in reject['1']],
                   '2': [vlist.ID for vlist in reject['2']]}
    obs_correct = {'1': sorted([vlist.ID for vlist in keep['1'].variant_list]),
                   '2': sorted([vlist.ID for vlist in keep['2'][0].variant_list])}
    assert_dict_equal(correct, obs_correct)
    assert_dict_equal(removed, obs_removed)


def populate_array(size, number):
    """ Create array of size initialized with values = number """
    array = np.zeros((size, size))
    # pylint: disable=consider-using-enumerate
    for i in range(len(array)):
        for j in range(i):
            array[i][j] = number
    return array


def test_pagb():
    """ Test probability of a given b calculation """
    records = [_Record(1, 0, '1', '', '', '', '', '', '', ''),
               _Record(2, 0, '2', '', '', '', '', '', '', ''),
               _Record(3, 0, '3', '', '', '', '', '', '', ''),
               _Record(4, 0, '4', '', '', '', '', '', '', '')]
    group = vg.VariantGroup('1', records)
    group.coverage_array = populate_array(len(group.coverage_array), 400)
    group.existence_array = populate_array(len(group.coverage_array), 100)
    group.a_not_b_array = populate_array(len(group.coverage_array), 25)
    group.a_not_b_array[1][0] = 100
    group.a_not_b_array[2][0] = 25
    group.a_not_b_array[3][1] = 100
    group.b_not_a_array[3][0] = 25
    group.b_not_a_array[2][1] = 100
    group.b_not_a_array[3][1] = 100
    group.set_filter_fq_pagb(50, True)
    for i in range(len(group.filter)):
        for j in range(i, len(group.filter)):
            group.filter[i][j] = True
    assert_true(not group.filter.all())
    group.filter[3][1] = True
    assert_true(group.filter.all())


def test_add_min_read_filter():
    """ Test minimum read filter """
    records = [_Record(1, 0, '1', '', '', '', '', '', '', ''),
               _Record(2, 0, '2', '', '', '', '', '', '', ''),
               _Record(3, 0, '3', '', '', '', '', '', '', ''),
               _Record(4, 0, '4', '', '', '', '', '', '', '')]
    group = vg.VariantGroup('1', records)
    arlen = len(records)
    group.existence_array = np.zeros((arlen, arlen))
    correct = {'1': np.zeros((arlen, arlen)),
               '2': np.zeros((arlen, arlen))}
    group.existence_array[0][1] = 30
    group.existence_array[0][2] = 40
    group.existence_array[0][3] = 50
    group.existence_array[1][2] = 5
    group.existence_array[1][3] = 10
    group.existence_array[2][3] = 1
    correct['1'][0][2] = 1
    correct['1'][0][3] = 1
    correct['2'][0][1] = 1
    correct['2'][0][2] = 1
    correct['2'][0][3] = 1
    correct['2'][1][3] = 1
    correct = {key: correct[key] > 0 for key in correct}
    group.add_filter_min_reads(30)
    assert_true((correct['1'] == group.filter).all())
    group.reset_filter()
    group.add_filter_min_reads(5)
    assert_true((correct['2'] == group.filter).all())


def test_vargroup_pagb():
    """Validates the output of the variant grouper on VCF files using P(A|B) mode"""
    vcf_in = resource_filename('celltics.tests.data.files', 'vargroup_in.vcf')
    bam = resource_filename('celltics.tests.data.files', 'vargroup.bam')
    vcf_out = resource_filename('celltics.tests.data.files', 'vargroup_out.vcf')
    output_file = os.path.join(os.path.dirname(vcf_out), 'vargroup_test_output_pagb.vcf')
    vg.main(input_file=vcf_in, output_file=output_file, bam_file=bam, merge_distance=1000, fq_threshold=50,
            max_allele_fraction=0.01, min_alt_fq=0.07, min_alt_depth=5, no_filter=True, bam_filter_mode='max_pagb',
            write_mode='merged_only')
    assert_true(filecmp.cmp(vcf_out, output_file))


def test_vargroup_append():
    """Validates the output of the variant grouper appending to VCF files"""
    vcf_in = resource_filename('celltics.tests.data.files', 'vargroup_in.vcf')
    bam = resource_filename('celltics.tests.data.files', 'vargroup.bam')
    vcf_out = resource_filename('celltics.tests.data.files', 'vargroup_append_out.vcf')
    output_file = os.path.join(os.path.dirname(vcf_out), 'vargroup_test_append_output.vcf')
    vg.main(input_file=vcf_in, output_file=output_file, bam_file=bam, merge_distance=50, fq_threshold=50,
            min_reads=3, no_filter=True, write_mode='append',
            bam_filter_mode='pagb')
    assert_true(filecmp.cmp(vcf_out, output_file))


def test_vargroup():
    """Validates the output of the variant grouper on VCF files"""
    vcf_in = resource_filename('celltics.tests.data.files', 'vargroup_in.vcf')
    bam = resource_filename('celltics.tests.data.files', 'vargroup.bam')
    vcf_out = resource_filename('celltics.tests.data.files', 'vargroup_out.vcf')
    output_file = os.path.join(os.path.dirname(vcf_out), 'vargroup_test_output.vcf')
    vg.main(input_file=vcf_in, output_file=output_file, bam_file=bam, merge_distance=1000, fq_threshold=5,
            max_allele_fraction=0.01, min_alt_fq=0.07, min_alt_depth=5, no_filter=True, write_mode='merged_only',
            bam_filter_mode='pab')
    assert_true(filecmp.cmp(vcf_out, output_file))


def test_get_reference_seq():
    """Validates the reference sequence matches what is expected"""
    chrom = '1'
    start = 1
    end = 9
    expected_seq = 'gaaggaact'
    expected_ukn_seq = 'tgatcacagg'
    ref_seq = resource_filename('celltics.tests.data.files', 'mini_genome.fasta')
    mydict = SeqIO.to_dict(SeqIO.parse(ref_seq, 'fasta'))
    reference_seq = vg.get_reference_seq(chrom, start, end, seq_dict=mydict)
    assert_equal(expected_seq, reference_seq)
    reference_seq = vg.get_reference_seq('M', start, end, seq_dict=mydict)
    assert_equal(expected_ukn_seq, reference_seq)
    assert_raises(Exception, vg.get_reference_seq, chrom, 10, 11, seq_dict={'1': 'NNNNNN'})
