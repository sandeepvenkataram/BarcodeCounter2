
import glob
import pathlib
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def clusterBarcodesDNAClust(args: dict, template_sequence: str) -> dict:
    ##
    # concat all barcode fastq files by experiment into a single file for clustering, remove those sequences that appear less than 3 times
    ##
    all_files = glob.glob(args['output_dir']+'*_barcode.fastq')
    dedup_file_name = args['output_dir']+'allSamplesConcatDedup.fasta'
    read_count_file_name = args['output_dir']+'allSamplesConcatDedup.readCounts'
    dnaclust_output_filename = args['output_dir'] + 'allSamplesConcatDedup.dnaclustOut'
    clustered_bc_filename = args['output_dir']+'clusteredBCsDNAClust.fasta'
    expected_barcode_length = template_sequence.count('X')


    unique_bc_lines = {}
    for fname in all_files:
        with open(fname, encoding='utf-8') as infile_handle:
            for line in SeqIO.parse(infile_handle, 'fastq'):
                myseq = str(line.seq)
                if (myseq not in unique_bc_lines):
                    unique_bc_lines[myseq] = 0
                unique_bc_lines[myseq] += 1
    read_counter = 1
    
    with open(dedup_file_name, 'w', encoding='utf-8') as ofh, open(read_count_file_name, 'w', encoding='utf-8') as rcfh:
        seqs_to_write = []
        for line, bc in unique_bc_lines.items():
            # remove any reads with Ns in it (< .5% of reads) or sequences with too few reads or sequences with barcodes of 0 length (if they somehow got missed) or sequences that are too long (more than 50% longer than the expected sequence length)
            if 'N' not in line and bc > 3 and len(line) > 0 and len(line) <= int(1.5*expected_barcode_length):
                seqs_to_write.append(SeqRecord(Seq(line), id=str(
                    read_counter), name='', description=''))
                rcfh.write(str(bc)+'\n')
                read_counter += 1
        SeqIO.write(seqs_to_write, ofh, 'fasta')

    # use DNAclust to cluster reads
    cwd = pathlib.Path(__file__).parent.resolve()
    call_function = [str(cwd / '../bin/dnaclust'), '-s', '.95', '--approximate-filter','-k', '6', '-t', str(args['num_threads']), '-i', dedup_file_name, '>'+dnaclust_output_filename]
    subprocess.call(call_function)

    # create final barcode fasta file using the centers of the clusters found by DNAclust
    bcs_to_use = []
    with open(dnaclust_output_filename, 'r', encoding='utf-8') as infile:
        for line in infile:
            line = line.strip()
            line_split = line.split('\t')
            bcs_to_use.append(int(line_split[0]))

    bcs_to_use.sort()
    total_num_bcs = len(bcs_to_use)
    bcs_to_use_index = 0
    records_to_write = []
    with open(dedup_file_name, 'r', encoding='utf-8') as infile, open(clustered_bc_filename, 'w', encoding='utf-8') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            if (bcs_to_use_index < total_num_bcs and int(str(record.id)) == bcs_to_use[bcs_to_use_index]):
                bcs_to_use_index += 1
                record.id = str(bcs_to_use_index)
                records_to_write.append(record)
        SeqIO.write(records_to_write, outfile, 'fasta')

    args['barcode_list_file'] = clustered_bc_filename
    return args