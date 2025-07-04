
from typing import Union
import argparse
from collections import defaultdict

AVAILABLE_ARGUMENTS = [
    {'name_or_flags': '-fastqDir', 'dest': 'fastq_dir', 'help': 'directory location of fastq files'},
    {'name_or_flags': '-outputDir', 'dest': 'output_dir', 'help': 'location of output directory'},
    {'name_or_flags': '-templateSeq', 'dest': 'template_seq_file', 'help': 'Template sequence of amplicon locus. This file contains a single line with standard DNA bases. UMI (unique molecular identifier}, sequences are coded as U, multiplexing indices are coded as D and barcode loci coded as X. If these features have different lengths between samples, define the template using the longest possible length of each feature. Every feature annotated must be covered by the sequencing data, and no feature can span the exact middle of the template sequence when using paired end data.'},
    {'name_or_flags': '-sample', 'dest': 'sample_file', 'help': 'File defining the samples present in the sequencing data. This file is tab delimited, with no header line. The column values are: Sample Name\t File Prefix\t internal multiplexing barcode 1\t internal multiplexing barcode 2... The internal multiplexing barcode columns must correspond to the names of the sequences in the multiBCFasta file. Do not use spaces in any of the columns for file / directory names, as this tends to behave poorly.'},
    {'name_or_flags': '-multiBCFasta', 'dest': 'multiBC_fasta_file', 'help': 'A multi-line fasta file defining multiplexing tag sequences. Required if there are multiplexing tags within the amplicon sequence as defined by the templateSeq file.'},
    {'name_or_flags': '-barcodeList', 'dest': 'barcode_list_file', 'help': 'Optional fasta file specifying the barcodes present in the sample. If file is not supplied, unique barcodes will be identified de novo. The name for each sequence must be unique. If the file is being generated manually, the sequence must be a simple concatenation of all barcode regions as defined in the template sequence in 5'-3' order.'},

    {'name_or_flags': '-barcode5PrimeTrimLength', 'dest': 'barcode_5prime_trim_length', 'default': 0,  'help': "Number of bp to trim from the 5' end of each barcode. "},
    {'name_or_flags': '-barcode3PrimeTrimLength', 'dest': 'barcode_3prime_trim_length', 'default': 0,  'help': "Number of bp to trim from the 3' end of each barcode. "},
    {'name_or_flags': '-bcNGapLength', 'dest': 'bc_Ngap_length', 'default': 0,  'help': 'Number of bp of Ns to put between barcodes from forward and reverse reads. Use if the sequence is not covering the entire barcode and you you have provided the list of valid barcodes using the barcodeList argument.'},
    {'name_or_flags': '-blastPath', 'dest': 'blast_path', 'help': 'BLAST installation directory if it is not in the Path already', 'default': ''},
    {'name_or_flags': '-useBowtie2', 'dest': 'use_bowtie2', 'help': 'Flag to use Bowtie2 instead of the default BWA mem for barcode mapping. ', 'action': 'store_true', 'default': 'true'},
    {'name_or_flags': '-bowtie2Path', 'dest': 'bowtie2_path', 'help': 'Bowtie2 installation directory if it is not in the Path already', 'default': ''},
    {'name_or_flags': '-bwaPath', 'dest': 'bwa_path', 'help': 'BWA installation directory if it is not in the Path already', 'default': ''},
    {'name_or_flags': '-demultiplexOnly', 'dest': 'demultiplex_only', 'action': 'store_true',  'help': 'Use flag if you want to only split the raw fastq files and not continue with the rest of the barcode counting. This is useful when distributing demultiplexing across several machines, i.e. in a cluster.'},
    {'name_or_flags': '-numThreads', 'dest': 'num_threads', 'default': 1,  'help': 'Number of threads to be used for computation.'},
    {'name_or_flags': '-pairedEnd', 'dest': 'is_paired_end', 'action': 'store_true',  'help': 'Use if sequencing data is paired end'},
    {'name_or_flags': '-readLength', 'dest': 'read_length', 'type': int, 'default': 100,  'help': 'Expected length of each read from sequencing machine. Default = 100. Reduce this number from the true read length if necessary such that non-constant regions of the barcode locus are not shared between reads. This does not modify the input fastq files, but effectively truncates reads before processing'},
    {'name_or_flags': '-remapBarcodes', 'dest': 'remap_barcodes', 'action': 'store_true',  'help': 'Set to True if you want to remap barcodes even if the files already exist'},
    {'name_or_flags': '-reverseRead2', 'dest': 'reverse_read2', 'action': 'store_true',  'help': 'Set to True if you want to reverse-complement read2 when using paired-end data'},
    {'name_or_flags': '-resplitFastq', 'dest': 'resplit_fastq', 'action': 'store_true',  'help': 'Use flag if you want to resplitthe raw fastq files (i.e. if you have already done this and want to do it again},.'},
    {'name_or_flags': '-useUMI', 'dest': 'use_umi', 'action': 'store_true',  'help': 'Use flag if you want to remove PCR duplicate reads using UMI data'},
]

def generate_argparser(description: str, provided_arg_names: Union[list[str],  None]=None) -> defaultdict:
    '''Generate an argparser for barcode counter'''
    arg_parser = argparse.ArgumentParser(description=description)
    arg_names = [k for k in AVAILABLE_ARGUMENTS if k['name_or_flags'] in provided_arg_names] if provided_arg_names else [k for k in AVAILABLE_ARGUMENTS]
    if not arg_names:
        raise ValueError(f'None of the provided arg names {provided_arg_names}, are valid!')
    for arg in arg_names:
        name = arg['name_or_flags']
        del arg['name_or_flags']
        arg_parser.add_argument(name, **arg)
    args = arg_parser.parse_args()
    return defaultdict(lambda: None, vars(args))
