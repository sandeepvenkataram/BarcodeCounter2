import os
import subprocess

from Bio import SeqIO

from lib.data import SampleFileMetadata


def map_barcodes(inputs: tuple[SampleFileMetadata, dict]):
    sample, args = inputs
    # only run on this sample if the output file doesn't exist or flag has been set
    index_string = sample.sample
    file_prefix = f'{args["output_dir"]}{index_string}'
    bc_fastq_file = f'{file_prefix}_barcode.fastq'
    bc_sam_file = f'{file_prefix}_barcode.sam'
    bc_id_file = f'{file_prefix}_readBarcodeID.txt'
    map_qual_file = f'{file_prefix}_readMappingQuality.txt'
    umi_seqs_file = f'{file_prefix}_UMISeqs.tab'
    barcode_counts_file = f'{file_prefix}__barcodeCounts.tab'
    unmapped_read_count_file = f'{file_prefix}_numUnmappedReads.txt'
    umi_dup_count_file = f'{file_prefix}_barcodeUMIDupCounts.tab'

    if os.path.isfile(bc_fastq_file) and (not os.path.isfile(barcode_counts_file) or args['remap_barcodes']):

        # make a dictionary to map barcode names in the barcode list fasta file to consecutive numbers for indexing in a vector.
        bc_name_to_idx_dict = {}
        total_num_bcs = 1
        with open(args['barcode_list_file'], 'r', encoding='utf-8') as infile:

            for record in SeqIO.parse(infile, 'fasta'):
                # we have found a duplicate entry in the barcode list. quit!
                if record.id in bc_name_to_idx_dict:
                    raise ValueError('Duplicate entry ' + record.id + ' found in input barcode list!')
                bc_name_to_idx_dict[record.id] = total_num_bcs
                total_num_bcs += 1

        if args['use_bowtie2']:  # bowtie2 call if flagged
            subprocess.call([args['bowtie2_path'] + 'bowtie2', '-L 10', '-q', '--very-sensitive-local', '-x ' + args['barcode_list_file'], '-U' + bc_fastq_file, '-S' + bc_sam_file])
        else:  # bwa mem call
            with open(bc_sam_file, 'w', encoding='utf-8') as outfile:
                subprocess.call([args['bwa_path'] + 'bwa', 'mem', '-k 10', '-y 12', args['barcode_list_file'], bc_fastq_file], stdout=outfile)

        # get the barcode match for each read and put into a single column output file. do the same for mapping quality
        with open(bc_id_file, 'w', encoding='utf-8') as outfile:
            subprocess.call("grep -v '^@' " + bc_sam_file + " | cut -f 3", stdout=outfile, shell=True)
        with open(map_qual_file, 'w', encoding='utf-8') as outfile:
            subprocess.call("grep -v '^@' " + bc_sam_file + " | cut -f 5", stdout=outfile, shell=True)

        bc_umi_map = {}
        bc_count_list = [0] * int(total_num_bcs - 1)
        bc_umi_dup_count_list = [0] * int(total_num_bcs - 1)
        total_unmapped_reads = 0

        # for each read bc / umi pair
        with open(bc_id_file, 'r', encoding='utf-8') as bch, open(map_qual_file, 'r', encoding='utf-8') as mqh, open(umi_seqs_file, 'r', encoding='utf-8') as umih:
            for bcid, mapq, umi_string in zip(bch, mqh, umih):
                bcid = bcid.strip()
                mapq = mapq.strip()
                if bcid == '*' or bcid == '' or mapq == '*' or int(mapq) < 20:
                    total_unmapped_reads = total_unmapped_reads + 1
                else:
                    # get the internal index corresponding to the matched barcode
                    bcid = bc_name_to_idx_dict[bcid]
                    mykey = str(bcid) + '\t' + umi_string
                    if mykey not in bc_umi_map or not args.UMI:
                        bc_umi_map[mykey] = 1
                        bc_count_list[int(bcid) - 1] = bc_count_list[int(bcid) - 1] + 1
                    else:
                        bc_umi_dup_count_list[int(bcid) - 1] = bc_umi_dup_count_list[int(bcid) - 1] + 1

        # write total count data to file
        with open(barcode_counts_file, 'w', encoding='utf-8') as outfile_handle:
            for count_val in bc_count_list:
                outfile_handle.write(str(count_val) + '\n')
        with open(unmapped_read_count_file, 'w', encoding='utf-8') as outfile_handle:
            outfile_handle.write(str(total_unmapped_reads))
        if args['use_umi']:
            with open(umi_dup_count_file, 'w', encoding='utf-8') as outfile_handle:
                for count_val in bc_umi_dup_count_list:
                    outfile_handle.write(str(count_val) + '\n')
