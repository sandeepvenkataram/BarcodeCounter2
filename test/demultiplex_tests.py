from collections import defaultdict
from multiprocessing import Pool
import os
import tempfile
import unittest

import pandas as pd
from Bio import SeqIO

from lib.demultiplex import demultiplex_fastq, get_best_blast_match, get_subject_match_coordinates
from lib.file_io import create_const_region_fasta, generate_index_to_sample_map, identify_used_fastq_files, parse_sample_file, parse_template_seq
from test.test_constants import TEST_PRIMER_INDEX_SEQ_FILE, TEST_RAW_FASTQ_FILE_DIR, TEST_SAMPLE_FILE, TEST_SEQUENCE_TEMPLATE


class DemultiplexTests(unittest.TestCase):
    
    def _generate_template_seqs(self, is_paired_end=False, read_length=100):
        with open(TEST_SEQUENCE_TEMPLATE, 'r', encoding='utf-8') as f:
            sequence = f.readline().strip()
        args = {
            'is_paired_end': is_paired_end,
            'read_length': read_length,
        }
        return parse_template_seq(sequence, args)

    def test_get_best_blast_match(self):
        self.assertIsNone(get_best_blast_match(pd.DataFrame()))
        df = pd.DataFrame([["seq1", 1, 100, 5, 105, "1e-10"]], columns=["sseqid", "qstart", "qend", "sstart", "send", "evalue"])
        self.assertEqual(get_best_blast_match(df), ["seq1", 1, 100, 5, 105])
        df = pd.DataFrame([
            ["seq1", 1, 100, 5, 105, "1e-10"],
            ["seq2", 1, 100, 5, 105, "1e-5"]
        ], columns=["sseqid", "qstart", "qend", "sstart", "send", "evalue"])
        self.assertEqual(get_best_blast_match(df), ["seq1", 1, 100, 5, 105])
        df = pd.DataFrame([
            ["seq1", 1, 100, 5, 105, "1e-3"],
            ["seq2", 1, 100, 5, 105, "1e-5"]
        ], columns=["sseqid", "qstart", "qend", "sstart", "send", "evalue"])
        self.assertIsNone(get_best_blast_match(df))
        df = pd.DataFrame([
            ["seq1", 1, 100, 5, 105, "1e-3"],
            ["seq2", 1, 100, 5, 105, "1e-3"]
        ], columns=["sseqid", "qstart", "qend", "sstart", "send", "evalue"])
        self.assertIsNone(get_best_blast_match(df))

    def test_get_subject_match_coordinates(self):
        self.assertEqual(get_subject_match_coordinates(["seq1", "20", "80", "20", "80"], {"seq1": 150}), [0, 150, False])
        self.assertEqual(get_subject_match_coordinates(["seq1", "20", "80", "20", "80"], {"seq1": 95}), [0, 95, False])
        self.assertEqual(get_subject_match_coordinates(["seq1", "20", "80", "100", "5"], {"seq1": 150}), [-31, 84, True])
        self.assertEqual(get_subject_match_coordinates(["seq1", "30", "70", "10", "100"], {"seq1": 100}), [20, 70, False])
        with self.assertRaises(ValueError):
            get_subject_match_coordinates(["seq1", "10", "50", "5", "100"], {"seq1": 20})
        
    def test_identify_used_fastq_files(self):
        parsed_template_seqs = self._generate_template_seqs(False, 100)
        sample_array = parse_sample_file(TEST_SAMPLE_FILE, parsed_template_seqs)
        with tempfile.TemporaryDirectory() as temp_dir:
            args = {
                'output_dir': str(temp_dir),
                'fastq_dir': TEST_RAW_FASTQ_FILE_DIR,
                'multiBC_fasta_file': TEST_PRIMER_INDEX_SEQ_FILE,
                'is_paired_end': True,
                'resplit_fastq': True,
                'blast_path': '',
                'num_threads': 3,
            }
            template_seq_lengths_dict = create_const_region_fasta(parsed_template_seqs, args)
            identify_used_fastq_files(sample_array, args)
            sample_map = defaultdict(list)
            for sample in sample_array:
                if sample.resplit:
                    sample_map[sample.file_prefix].append(sample)
            sample_grouped_arrays = [v for _, v in sample_map.items()]
            used_fastq_files = [(x, args, parsed_template_seqs, template_seq_lengths_dict) for x in sample_grouped_arrays]  # need to add in args as a tuple as we are running this through pool.map
            self.assertTrue(len(used_fastq_files) > 0)
            with Pool(processes = int(args['num_threads'])) as pool:
                pool.map(demultiplex_fastq, used_fastq_files)
            for sample in sample_array:
                self.assertTrue(os.path.exists(args['output_dir']+sample_array[0].file_prefix+'_numReadsFoundPerSample.txt'))
                