import os
import tempfile
import unittest

import pandas as pd
from Bio import SeqIO

from lib.constants import ALL_CONST_REGIONS_FILE_NAME
from lib.data import TemplateSeqFeature
from lib.file_io import (
    create_const_region_fasta,
    generate_index_to_sample_map,
    identify_used_fastq_files,
    parse_sample_file,
    parse_template_seq,
)
from test.test_constants import (
    TEST_PRIMER_INDEX_SEQ_FILE,
    TEST_RAW_FASTQ_FILE_DIR,
    TEST_SAMPLE_FILE,
    TEST_SEQUENCE_TEMPLATE,
)


class IOTests(unittest.TestCase):

    def _generate_template_seqs(self, is_paired_end=False, read_length=100):
        with open(TEST_SEQUENCE_TEMPLATE, 'r', encoding='utf-8') as f:
            sequence = f.readline().strip()
        args = {
            'is_paired_end': is_paired_end,
            'read_length': read_length,
        }
        return parse_template_seq(sequence, args)

    def test_parse_template_seq_single_end(self):
        parsed_template_seqs = self._generate_template_seqs(False, 100)
        self.assertTrue(len(parsed_template_seqs) == 1)
        parsed_template_seq = parsed_template_seqs[0]
        self.assertTrue(parsed_template_seq.get_expected_barcode_length() == 28)
        template_seq_array = parsed_template_seq.get_template_seq_array()
        expected_template_seq_array = [
            TemplateSeqFeature('U', 'UUUUUUUU', 0, 0, None),
            TemplateSeqFeature('D', 'DDDDDDDD', 0, 1, None),
            TemplateSeqFeature('CAAGCTTAGATCTGATATCGGTACCCAACAAACCACATTATGTCTTTAGCGATAACTTCGTATAGCATACATTATACGAAGTTAT', 'CAAGCTTAGATCTGATATCGGTACCCAACAAACCACATTATGTCTTTAGCGATAACTTCGTATAGCATACATTATACGAAGTTAT', 0, 2, 'const_region_0_2'),
            TemplateSeqFeature('X', 'XXXXXXXXXXXXXXXXXXXXXXXXXXXX', 0, 3, None),
            TemplateSeqFeature('GGTACCGATATCGGATCCGTCGACAAAAGCCTCCTTTAG', 'GGTACCGATATCGGATCCGTCGACAAAAGCCTCCTTTAG', 0, 4, 'const_region_0_4'),
            TemplateSeqFeature('D', 'DDDDDD', 0, 5, None),
            TemplateSeqFeature('U', 'UUUUUUUU', 0, 6, None),
        ]
        self.assertEqual(len(template_seq_array), len(expected_template_seq_array))
        for i, j in zip(template_seq_array, expected_template_seq_array):
            self.assertEqual(i, j)

    def test_parse_template_seq_paired_end(self):
        parsed_template_seqs = self._generate_template_seqs(True, 100)
        self.assertTrue(len(parsed_template_seqs) == 2)

    def test_parse_sample_file(self):
        parsed_template_seqs = self._generate_template_seqs(False, 100)
        sample_array = parse_sample_file(TEST_SAMPLE_FILE, parsed_template_seqs)
        index_to_sample_map = generate_index_to_sample_map(sample_array)
        self.assertEqual(len(sample_array), 4)
        sample_file_df = pd.read_csv(TEST_SAMPLE_FILE, sep='\t', header=None)
        sample_file_df.columns = ['sample', 'fastq', 'index1', 'index2']
        for i, sample in enumerate(sample_array):
            row = sample_file_df.iloc[i]
            self.assertEqual(sample.sample, row['sample'])
            self.assertEqual(sample.file_prefix, row['fastq'])
            self.assertEqual(sample.int_multi_bc_array, [row['index1'], row['index2']])
            self.assertTrue(sample.sample in [v for _, v in index_to_sample_map.items()])

    def test_identify_used_fastq_files(self):
        parsed_template_seqs = self._generate_template_seqs(False, 100)
        sample_array = parse_sample_file(TEST_SAMPLE_FILE, parsed_template_seqs)
        with tempfile.TemporaryDirectory() as temp_dir:
            args = {
                'output_dir': str(temp_dir),
                'fastq_dir': TEST_RAW_FASTQ_FILE_DIR,
                'is_paired_end': True,
                'resplit_fastq': False,
            }
            identify_used_fastq_files(sample_array, args)
            for sample in sample_array:
                for file in [sample.fwd_fastq, sample.rev_fastq, sample.bc_fastq, sample.r1_fastq, sample.r2_fastq, sample.umi_tab]:
                    self.assertIsNotNone(file)
                    self.assertTrue(os.path.exists(file))

    def test_create_const_region_fasta(self):
        parsed_template_seqs = self._generate_template_seqs(False, 100)
        with tempfile.TemporaryDirectory() as temp_dir:
            args = {
                'output_dir': str(temp_dir),
                'multiBC_fasta_file': TEST_PRIMER_INDEX_SEQ_FILE,
                'blast_path': '',
            }
            template_seq_lengths_dict = create_const_region_fasta(parsed_template_seqs, args)
            with open(args['multiBC_fasta_file'], 'r', encoding='utf-8') as infile:
                bc_fastq_lines = ''.join(infile.readlines())
            with open(args['output_dir'] + ALL_CONST_REGIONS_FILE_NAME, 'r', encoding='utf-8') as infile:
                const_regions_lines = ''.join(infile.readlines())
            for record in SeqIO.parse(args['output_dir'] + ALL_CONST_REGIONS_FILE_NAME, 'fasta'):
                self.assertEqual(template_seq_lengths_dict.get(record.id), len(record.seq))
            self.assertTrue(bc_fastq_lines in const_regions_lines)
            for template_seq in parsed_template_seqs:
                for seq in template_seq.get_template_seq_array():
                    if seq.fasta_feature_name is not None:
                        self.assertTrue(seq.fasta_feature_name in const_regions_lines)
                        self.assertTrue(seq.sequence in const_regions_lines)

    def test_generate_final_tables(self):
        # TODO: COMPLETE THIS METHOD!!
        self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
