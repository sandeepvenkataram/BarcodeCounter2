
from collections import defaultdict
from multiprocessing import Pool
import subprocess
from pathlib import Path

from lib.argparse import generate_argparser
from lib.clustering import clusterBarcodesDNAClust
from lib.demultiplex import demultiplex_fastq
from lib.file_io import create_const_region_fasta, generate_final_tables, identify_used_fastq_files, parse_sample_file, parse_template_seq
from lib.mapping import map_barcodes


def main():
    """Main function"""
    ####
    # Read Command-line Arguments
    ####
    description = "Process Illumina PCR Amplicon Fastq files to count barcode frequencies. Must have blastn installed to run. Make sure there are no spaces in any file names or directory paths, the program will not work otherwise. Make sure that files are already split by their illumina indices (N700 / S500 index). \n\nGenerated Files:\n\nnumReadsFoundPerSample: One file per input fastq file, it details the assignment of each read to a combination of inline indices if they exist, and the subsequent filtering of reads. The file has 2 columns, the first being a unique identifier for the fastq file and inline index combination identified, and the second being an array of 6 numbers binning the reads into the following categories: \n\t1. correct reads used for subsequent mapping. \n\t2. There are too many Ns in the read.\n\t3. UMIs have been requested and not identified in the read/\n\t4. There is no barcode found in the read.\n\t5. The combination of indices do not match a sample specified in the sampleFile.\n\t6. A second check for finding a valid barcode in the read.\n\n For each sample in SampleFile with at least one barcode found a number of files are generated.\n R1(R2).fastq - raw (possibly truncated depending on the readLength parameter) reads associated with this sample\n barcode.fastq - the portion of the reads associated with all barcode sequences concatenated together, in the same order as the R1/R2.fastq files\n UMISeqs.tab - tab delmited UMI sequences if they exist and are being used, in the same order as the R1/R2.fastq files.\n readBarcodeID.txt - the ID number of the barcode in the clusteredBCs.fasta file that each read was mapped to, in the same order as the R1/R2.fastq files.\nbarcodeCalls.tab - Count of each barcode in the sample, removing UMI duplicates if asked for. Line 1 contains the counts for barcode 1 (defined in clusteredBCs.fasta), line 2 for barcode 2 etc.\n\nallBarcodeCalls.tab - tab delimited concatenation of all of the barcodeCalls.tab files, with a header row identifying which column comes from which sample."
    args = generate_argparser(description=description)
    args['bc_Ngap_length'] = int(args['bc_Ngap_length']) if args['bc_Ngap_length'] is not None else 0

    # create the output directory
    Path(args['output_dir']).mkdir(parents=True, exist_ok=True)


    with open(args['template_seq_file'], encoding='utf-8') as f:
        sequence = f.readline().strip()
    template_array = parse_template_seq(sequence, args)
    sample_array = parse_sample_file(args['sample_file'], template_array)
    template_seq_length_dict = create_const_region_fasta(template_array, args)
    

    ## Demultiplex data using multiprocessing if there are any to be demultiplexed
    identify_used_fastq_files(sample_array, args)
    sample_map = defaultdict(list)
    for sample in sample_array:
        if sample.resplit:
            sample_map[sample.file_prefix].append(sample)
    sample_grouped_arrays = [v for _, v in sample_map.items()]
    used_fastq_files = [(x, args, template_array, template_seq_length_dict) for x in sample_grouped_arrays]  # need to add in args as a tuple as we are running this through pool.map
    if(len(used_fastq_files)>0):
        with Pool(processes = int(args['num_threads'])) as pool:
           pool.map(demultiplex_fastq, used_fastq_files)
            
    ## If we are only trying to demultiplex, quit now
            
    if args['demultiplex_only']:
        # print("Terminating after splitting raw fastq files as requested.")
        return


    ## Cluster barcodes using DNAClust if necessary
        
    if args['barcode_list_file'] is None:
        clusterBarcodesDNAClust(args, sequence)


    ## Make database from barcode fasta file for mapping
    if(args['use_bowtie2']):
        subprocess.call([args['bowtie2_path']+"bowtie2-build",args['barcode_list_file'],args['barcode_list_file']])
    else:
        subprocess.call([args['bwa_path']+"bwa","index",args['barcode_list_file']])

    ## Map barcodes with multiprocessing

    mapping_args = [(sample, args) for sample in sample_array]
    with Pool(processes = int(args['num_threads'])) as pool:
        pool.map(map_barcodes, mapping_args)


    ## Generate final output table	

    generate_final_tables(args)

if __name__ == '__main__':
    main()