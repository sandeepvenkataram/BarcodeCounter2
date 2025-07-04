from typing import Iterable
from Bio import SeqIO
import csv
import glob
import itertools as IT
import subprocess
import pandas as pd
from lib.constants import ALL_CONST_REGIONS_FILE_NAME
from lib.data import SampleFileMetadata, TemplateSequence


def parse_template_seq(sequence: str, args: dict) -> list[TemplateSequence]:
    template_array: list[TemplateSequence] = []
    assert all(y in args for y in ['is_paired_end', 'read_length']), 'Missing expected arguments when parsing the sequence template!'

    if (args['is_paired_end']):
        rlen = int(args['read_length'])
        assert rlen < len(sequence), 'You think you need paired-end analysis when the provided read length is longer than the length of the full sequence template??'
        template_array.append(TemplateSequence(sequence[0:rlen], 0))
        template_array.append(TemplateSequence(sequence[int(len(sequence)-rlen):], 1))
    else:
        template_array.append(TemplateSequence(sequence, 0))

    return template_array


def parse_sample_file(sample_file: str, template_array: list[TemplateSequence]) -> list[SampleFileMetadata]:
    sample_array = []
    num_inline_indices = 0
    for seq_array in template_array:
        for seq_feature in seq_array.get_template_seq_array():
            if seq_feature.is_index():
                num_inline_indices += 1
    sample_df: pd.DataFrame = pd.read_csv(sample_file, sep='\t', header=None)
    if sample_df.shape[1] < 3 or sample_df.shape[1] > 4:
        raise ValueError(
            f'The sample file {sample_file} doesnt have 3 or 4 columns, as expected')
    colnames = ['sample', 'fastq', 'index1']
    index_cols = ['index1']
    if len(sample_df.columns) == 4:
        colnames += ['index2']
        index_cols += ['index2']
    sample_df.columns = colnames
    sample_df = sample_df.drop_duplicates()
    for _, row in sample_df.iterrows():
        sample_array.append(SampleFileMetadata(
            row['sample'], row['fastq'], list(row[index_cols])))
    return sample_array


def generate_index_to_sample_map(inputs: Iterable[SampleFileMetadata]):
    return {input.file_prefix + '_' + '_'.join(input.int_multi_bc_array): input.sample for input in inputs}


def identify_used_fastq_files(sample_array: list[SampleFileMetadata], args: dict) -> None:
    for sample in sample_array:
        sample.identify_fastq_files(args)
        sample.touch_files(args)


def create_const_region_fasta(template_array: list[TemplateSequence], args: dict) ->dict:
    template_seq_lengths_dict = {}
    const_region_strings = []
        
    for template_seq in template_array:
        for seq_feature in template_seq.get_template_seq_array():
            if seq_feature.is_constant():
                const_region_strings.append(f'>{seq_feature.fasta_feature_name}\n{seq_feature.sequence}\n')
                
    with open(args['output_dir']+ALL_CONST_REGIONS_FILE_NAME, 'w', encoding='utf-8') as outfile:
        if (args['multiBC_fasta_file'] != None):
            with open(args['multiBC_fasta_file'], 'r', encoding='utf-8') as infile:
                outfile.write(''.join(infile.readlines())+'\n')
        outfile.write(''.join(const_region_strings))

    for record in SeqIO.parse(args['output_dir']+ALL_CONST_REGIONS_FILE_NAME, 'fasta'):
        template_seq_lengths_dict[record.id] = len(record.seq)
    blastCall = [args['blast_path']+'makeblastdb', '-in',
                 args['output_dir']+ALL_CONST_REGIONS_FILE_NAME, '-dbtype', 'nucl']
    subprocess.call(blastCall)
    return template_seq_lengths_dict


def generate_final_tables(args: dict):
    # print barcode counts as giant tab delimited table, with 1st column as barcode ID number and header file being the sample each column comes from
    patternString = args['output_dir']+'*_barcodeCounts*.tab'
    filenames = glob.glob(patternString)
    handles = [open(filename, 'r', encoding='utf-8') for filename in filenames]
    readers = [csv.reader(f, delimiter=',') for f in handles]
    filenames2 = ['BCID']
    filenames2.extend(filenames)

    BCNameToIdxDict = {}
    IdxToBCNameDict = {}
    with open(args['barcode_list_file'], 'r', encoding='utf-8') as infile:
        counter = 1
        for record in SeqIO.parse(infile, 'fasta'):
            # we have found a duplicate entry in the barcode list. quit!
            if (record.id in BCNameToIdxDict):
                raise ValueError('Duplicate entry '+record.id +
                                 ' found in input barcode list!')
            BCNameToIdxDict[record.id] = counter
            IdxToBCNameDict[counter] = record.id
            counter += 1

    with open(args['output_dir']+'allBarcodeCounts.tab', 'w', encoding='utf-8') as h:
        writer = csv.writer(h, delimiter='\t', lineterminator='\n', )
        writer.writerow(filenames2)
        i = 1
        for rows in IT.zip_longest(*readers, fillvalue=['']*2):
            combined_row = [IdxToBCNameDict[i]]
            for row in rows:
                row = row[:1]  # select the columns you want
                if len(row) == 1:
                    combined_row.extend(row)
                else:
                    combined_row.extend([''])
            writer.writerow(combined_row)
            i = i + 1
    for f in handles:
        f.close()
