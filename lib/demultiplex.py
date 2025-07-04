from collections import defaultdict
from io import StringIO, TextIOWrapper
import os
from pathlib import Path
import tempfile
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip
import subprocess

import pandas as pd
from lib.constants import ALL_CONST_REGIONS_FILE_NAME, CONSTANT_REGION_BLAST_PARAMS, FILE_BUFFER_SIZE, MAX_N_IN_READS
from lib.data import IndexCounter, SampleFileMetadata, TemplateSeqFeature, TemplateSequence


def _generate_fastq_file_handle(filename: str):
    if not os.path.exists(filename):
        raise ValueError(f'File path {filename} not found!')
    if (filename.endswith('gz')):
        res = gzip.open(filename, 'rt')
    else:
        res = open(filename, encoding='utf-8')
    return res


def demultiplex_fastq(inputs: tuple[list[SampleFileMetadata], dict, list[TemplateSequence], list[list[str | None]], dict]):
    '''Given a list of SampleFileMetadata that all correspond to the same set of Fastq Files, demultiplex reads'''
    sample_metadata_list, args, template_array, template_seq_length_dict = inputs
    if len(sample_metadata_list) == 0:
        return
    first_sample = sample_metadata_list[0]
    outfile_name = args['output_dir']+first_sample.file_prefix+'_numReadsFoundPerSample.txt'
    # first check that all of the samples come from the same fastq files
    fwd_fastqs = {sample.fwd_fastq for sample in sample_metadata_list}
    rev_fastqs = {sample.rev_fastq for sample in sample_metadata_list}
    if len(fwd_fastqs) != 1 or len(rev_fastqs) != 1:
        raise ValueError(f'Samples: {inputs} do not all use the same set of fastq files')

    
    # outputs: barcode file fastq, total split fastq, umi sequences, unidentified reads
    index_counter = {}
    fwd_fastq_handle = _generate_fastq_file_handle(first_sample.fwd_fastq)
    rev_fastq_handle = _generate_fastq_file_handle(first_sample.rev_fastq) if first_sample.rev_fastq is not None else None
    fastq_generator = zip(SeqIO.parse(fwd_fastq_handle, 'fastq'), SeqIO.parse(rev_fastq_handle, 'fastq')) if rev_fastq_handle is not None else SeqIO.parse(fwd_fastq_handle, 'fastq')

    read_list = []

    bad_r1_filename, bad_r2_filename = args['output_dir']+first_sample.file_prefix+'_unmappedReads_R1.fastq', args['output_dir']+first_sample.file_prefix+'_unmappedReads_R2.fastq'
    bad_fwd_reads_handle = open(bad_r1_filename, 'w', encoding='utf-8')
    if rev_fastq_handle is not None:
        bad_rev_reads_handle = open(bad_r2_filename, 'w', encoding='utf-8')
    else:
        bad_rev_reads_handle = None
    

    for read_counter, reads in enumerate(fastq_generator):
        if rev_fastq_handle is not None:
            fwd_rec, rev_rec = reads
        else:
            fwd_rec, rev_rec = reads, None

        fwd_rec = fwd_rec[0:args['read_length']]

        if rev_rec is not None:
            rev_rec = rev_rec[0:args['read_length']]
            rev_rec = rev_rec.reverse_complement() if args['reverse_read2'] else rev_rec

        read_list.append([fwd_rec, rev_rec])
        if read_counter > 1 and read_counter % FILE_BUFFER_SIZE == 0:
            index_counter = demultiplex_fastq_helper(read_list, sample_metadata_list, index_counter, bad_fwd_reads_handle, bad_rev_reads_handle, args, template_array, template_seq_length_dict)
            read_list = []
    if len(read_list) > 0:
        index_counter = demultiplex_fastq_helper(read_list, sample_metadata_list, index_counter, bad_fwd_reads_handle, bad_rev_reads_handle, args, template_array, template_seq_length_dict)

    if rev_fastq_handle is not None:
        rev_fastq_handle.close()
        bad_rev_reads_handle.close()

    bad_fwd_reads_handle.close()
    fwd_fastq_handle.close()

    ##
    # finalizing for both single and paired end data, count total read statistics here
    ##

    rows = []
    for _, counter_obj in index_counter.items():
        rows.append(counter_obj.to_dict())
    pd.DataFrame(rows).to_csv(outfile_name, sep='\t', index=False)


# Helper function that processes a single read during demultiplexing
def demultiplex_fastq_helper(read_list: list[list[SeqRecord]], sample_metadata_list: list[SampleFileMetadata], index_counter: dict[str, IndexCounter], bad_fwd_reads_handle: TextIOWrapper, bad_rev_reads_handle: TextIOWrapper, args: dict,  template_array: list[TemplateSequence], template_seq_length_dict: dict):
    '''Go through each read pair, do QC and figure out which sample (i.e. a specific combination of index sequences) the read should be assigned to'''
    # these are all lists (or dictionaries of lists) of output so that we don't have to do too many seqio write calls since they are slow.
    bad_fwd_reads_list = []
    bad_rev_reads_list = []
    final_bc_record_list = defaultdict(list)
    fwd_record_list = defaultdict(list)
    rev_record_list = defaultdict(list)
    umi_list = defaultdict(list)

    index_to_sample_map = {sample.index_key(): sample.sample for sample in sample_metadata_list}

    expected_num_bcs = 0
    for template in template_array:
        for feature in template.get_template_seq_array():
            if feature.is_barcode():
                expected_num_bcs += 1

    extracted_regions = extract_regions_from_fastq(read_list, args, template_array, template_seq_length_dict)

    for reads, region_group in zip(read_list, extracted_regions):  # for each read we are processing
        fwd_record = reads[0]
        rev_record = reads[1]
        
        if region_group == []:
            # if we didn't find any matches for this read, dump it in the bad reads list!
            bad_fwd_reads_list.append(fwd_record)
            if args['is_paired_end']:
                bad_rev_reads_list.append(rev_record)
            continue

        identified_bc_seq_records, identified_umi_sequences, identified_index_bcs = region_group

        total_n_count = str(fwd_record.seq).count('N')
        if args['is_paired_end']:
            total_n_count = total_n_count + str(rev_record.seq).count('N')

        # this is a string that uniquely defines each sample that was multiplexed. This must correspond to the SampleFileMetadata.index_key() method
        sample_index = '_'.join([sample_metadata_list[0].file_prefix] + identified_index_bcs)
        if sample_index not in index_counter:
            index_counter[sample_index] = IndexCounter(sample_index)


        # See if we have problems with the reads
        
        
        if (sample_index not in index_to_sample_map):
            # if we couldn't find a sample associated with this multiplexing barcode combination
            flag = 4
        elif (len(identified_bc_seq_records) != expected_num_bcs):
            # if there is not the expected number of barcodes
            flag = 3
        elif (args['use_umi'] and len(identified_umi_sequences) == 0):
            # if there is no UMI and UMI is expected
            flag = 2
        elif (total_n_count > MAX_N_IN_READS):
            # if there are too many Ns in the read
            flag = 1
        else:
            flag=0

        index_counter[sample_index].update(flag)

        if flag == 0:  # the read matches a valid sample
            my_sample = index_to_sample_map[sample_index]

            # concat all BCs associated with this read to get final BC
            final_bc = None
            for bc_record in identified_bc_seq_records:
                if not bc_record:
                    continue
                if not final_bc:
                    final_bc = bc_record
                else:
                    curID = final_bc.id
                    if (args['bc_Ngap_length'] > 0):
                        NGap = SeqRecord(
                            Seq('N' * args['bc_Ngap_length']), id='NGap')
                        NGap.letter_annotations['phred_quality'] = [
                            40]*args['bc_Ngap_length']
                        final_bc = final_bc + NGap + bc_record
                    else:
                        final_bc = final_bc + bc_record
                    final_bc.id = curID
            umi_list[my_sample].append('\t'.join([str(x) for x in identified_umi_sequences]))
            final_bc_record_list[my_sample].append(final_bc)
            fwd_record_list[my_sample].append(fwd_record)
            if args['is_paired_end']:
                rev_record_list[my_sample].append(rev_record)
        else:
            bad_fwd_reads_list.append(fwd_record)
            if args['is_paired_end']:
                bad_rev_reads_list.append(rev_record)

    # write outputs to files
    SeqIO.write(bad_fwd_reads_list, bad_fwd_reads_handle, 'fastq')
    if args['is_paired_end']:
        SeqIO.write(bad_rev_reads_list, bad_rev_reads_handle, 'fastq')

    for my_sample, bc_record_list in final_bc_record_list.items():
        matching_sample_metadata = [sample for sample in sample_metadata_list if sample.sample == my_sample]
        if len(matching_sample_metadata) != 1:
            raise ValueError(f'Something is very wrong here with 0 / multiple matches! {matching_sample_metadata}')

        with open(args['output_dir']+my_sample+'_barcode.fastq', 'a+', encoding='utf-8') as f:
            # write BC portion of read to a fastq (concat among all BC features)
            SeqIO.write(bc_record_list, f, 'fastq')

        with open(args['output_dir']+my_sample+'_R1.fastq', 'a+', encoding='utf-8') as f:
            # write raw fwd read to sample read file
            SeqIO.write(fwd_record_list[my_sample],f, 'fastq')

        with open(args['output_dir']+my_sample+'_UMISeqs.tab', 'a+', encoding='utf-8') as f:
            f.write('\n'.join(umi_list[my_sample]))

        if args['is_paired_end']:
            with open(args['output_dir']+my_sample+'_R2.fastq', 'a+', encoding='utf-8') as f:
                # write raw rev read to sample read file
                SeqIO.write(rev_record_list[my_sample],f, 'fastq')
    return index_counter


def get_best_blast_match(blast_output_df: pd.DataFrame) -> pd.DataFrame | None:
    '''Return the blast hit with lowest e value, excluding the evalue. Requires a dataframe input in outfmt6 format, output is an array with the columns separated as strings.'''
    if blast_output_df.empty:  # if there are no hits, return empty array
        return None

    blast_output_df = blast_output_df[['sseqid', 'qstart', 'qend', 'sstart', 'send', 'evalue']]
    if blast_output_df.shape[0] == 1:  # if there is only one hit, return it
        return blast_output_df.iloc[0, :]
    # if the first hit is better than the second, return it
    evalues = blast_output_df['evalue'].tolist()[:1]
    if float(evalues[0]) < float(evalues[1]):
        return blast_output_df.iloc[1, :]
    return None  # there are multiple hits!


def get_subject_match_coordinates(top_blast_result, template_seq_length_dict) -> list:
    '''Extract coordinates of blast match dealing with reverse complemented sequences if necessary'''
    starting_coor = top_blast_result['qstart']
    ending_coor = top_blast_result['qend']
    if top_blast_result['sseqid'] not in template_seq_length_dict:
        raise ValueError('We dont know the length of this sequence!')
    if template_seq_length_dict[top_blast_result['sseqid']] < top_blast_result['qend']:
        raise ValueError('This sequence is somehow shorter than what was blast-ed!')

    rev = False
    if (top_blast_result['sstart'] < top_blast_result['send']):  # forward orientation
        starting_coor = starting_coor - top_blast_result['sstart']
        ending_coor = ending_coor + template_seq_length_dict[top_blast_result['sseqid']] - top_blast_result['send']
    else:
        rev = True
        starting_coor = starting_coor - template_seq_length_dict[top_blast_result['sseqid']] + top_blast_result['sstart'] - 1
        ending_coor = ending_coor + top_blast_result['send'] - 1
    return [starting_coor, ending_coor, rev]


def _generate_blast_results_df(read_seq_record_list: list[list[SeqIO.SeqRecord]], args: dict, template_array: list[TemplateSequence]):
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', prefix='fasta_file_for_blast_', delete_on_close=False) as file:
        seq_record_array = []
        for read_id, read_seq_record in enumerate(read_seq_record_list):
            for read_num in range(len(template_array)):
                seq_record_array.append(SeqRecord(read_seq_record[read_num].seq, id=f'read_{read_id}_{read_num}', name='', description=''))
        SeqIO.write(seq_record_array, file, 'fasta')
        file.close()
        # blast reads against constant and multiplexing index sequences
        blast_command = [args['blast_path']+'blastn', '-query', file.name, '-db', args['output_dir']+ALL_CONST_REGIONS_FILE_NAME]
        blast_command.extend(CONSTANT_REGION_BLAST_PARAMS)
        blast_result = subprocess.check_output(blast_command).decode('ascii').rstrip().split('\n')

    # put blast results into a dataframe with some metadata
    blast_result_df = pd.read_csv(StringIO('\n'.join(blast_result)), sep='\t', header=None)
    blast_result_df.columns = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')  # default columns for outfmt 6
    blast_result_df['read_id'] = blast_result_df['qseqid'].apply(lambda x: int(x.split('_')[1]))
    blast_result_df['read_num'] = blast_result_df['qseqid'].apply(lambda x: int(x.split('_')[2]))

    blast_result_df['is_const_region'] = blast_result_df['sseqid'].apply(lambda x: 'const_region' in x)
    return blast_result_df.sort_values(['read_id', 'evalue'])


def _calculate_feature_coordinates(template_seq_array: list[TemplateSeqFeature], template_seq_length_dict: dict, read_subdf: pd.DataFrame, read_length: int):
    prev_segment_starting_coord = -1
    reversed_read = False
    first_index = -1

    starting_coordinates = [None]*(len(template_seq_array)+1) # this tracks the boundaries within the read of all of the features
    good_coords = [0]*(len(template_seq_array)+1)
    identified_index_bcs = []
        
    # we need to figure out the locations of constant regions and index regions
    for i, template_seq_feature in enumerate(template_seq_array):  # for each feature in the template
        template_seq_is_index = template_seq_feature.is_index()
        if not template_seq_is_index and not template_seq_feature.is_constant():  # D is for index sequences and constant sequences are the only ones longer than a single character
            continue
        # if this is an indexing barcode location or a constant region
        
        # use an xor to get the appropriate rows
        feature_subdf = read_subdf.loc[read_subdf.apply(lambda x: x.is_const_region ^ template_seq_is_index, axis=1)]
        if template_seq_feature.fasta_feature_name is not None:
            feature_subdf = feature_subdf.loc[feature_subdf['sseqid'] == template_seq_feature.fasta_feature_name]
        top_blast_result = get_best_blast_match(feature_subdf)

        # if we don't have a good match for this region
        if not top_blast_result:
            continue

        if (first_index == -1):
            first_index = i
        # figure out where in the read this index region is hitting, set coordinates and track the first region we are mapping to account for possibly needing to reverse the match location
        match_coords = get_subject_match_coordinates(top_blast_result, template_seq_length_dict)

        if (prev_segment_starting_coord > match_coords[0] and not reversed_read):
            reversed_read = True

        # if this is the first segment we have found a hit for, figure out the right index where the segment is
        if (prev_segment_starting_coord == -1):
            prev_segment_starting_coord = match_coords[0]
        if (not reversed_read):  # if the read is not reversed, add in the coordinates properly
            starting_coordinates[i] = match_coords[0]
            starting_coordinates[i+1] = match_coords[1]
        else:
            if (first_index >= 0):  # if it is reversed and this is not the first segment we have hit, assume the previous segment we hit needs to have its indices swapped
                tmp_val = starting_coordinates[first_index]
                starting_coordinates[first_index] = starting_coordinates[first_index+1]
                starting_coordinates[first_index+1] = tmp_val
                first_index = -2
            starting_coordinates[i] = match_coords[1]
            starting_coordinates[i+1] = match_coords[0]
        good_coords[i] = 1
        good_coords[i+1] = 1

        if template_seq_is_index:
            identified_index_bcs.append(top_blast_result['sseqid'].values[0])

    # put the 0 in the right place depending on if the read is reversed or not relative to the template if we haven't found it already
    if (not reversed_read and starting_coordinates[0] is None):
        starting_coordinates[0] = 0
    if (reversed_read and starting_coordinates[len(starting_coordinates)-1] is None):
        starting_coordinates[len(starting_coordinates)-1] = 0

    # set the end of the read to be readLength in the right place depending on if the read is reversed or not relative to the template and if we haven't found it already
    if (reversed_read and starting_coordinates[0] is None):
        starting_coordinates[0] = read_length
    if (not reversed_read and starting_coordinates[len(starting_coordinates)-1] is None):
        starting_coordinates[len(starting_coordinates)-1] = read_length

    if 0 in good_coords:  # if we at least have a good coordinate for the 5' most position in the template
        # use expected length of features to fill in missing coordinates, templated off of mapped coordinates if possible
        for i, template_seq in enumerate(template_seq_array):
            if template_seq.is_constant():
                continue
            if starting_coordinates[i+1] is None and starting_coordinates[i] is not None and i in good_coords:
                # if we have a 5' non-inferred coordinate but not a 3' coordinate
                starting_coordinates[i+1] = starting_coordinates[i] + template_seq.seq_length
            elif starting_coordinates[i] is None and starting_coordinates[i+1] is not None and (i+1) in good_coords:
                # if we have a 3' non-inferred coordinate but not a 5' coordinate
                starting_coordinates[i] = starting_coordinates[i +1] - template_seq.seq_length
            elif starting_coordinates[i+1] is None and starting_coordinates[i] is not None:
                # if we have a 5' coordinate whether or not it is inferred, but not a 3' coordinate
                starting_coordinates[i+1] = starting_coordinates[i] + template_seq.seq_length
            else:
                continue
    return starting_coordinates, identified_index_bcs


# This gets UMI, multiplexing index and barcode regions from each read via blast

def extract_regions_from_fastq(read_seq_record_list: list[list[SeqIO.SeqRecord]], args: dict, template_array: list[TemplateSequence], template_seq_length_dict: dict):
    # make a fasta file from all reads we are processing and blast against database of all index and constant regions
    blast_result_df = _generate_blast_results_df(read_seq_record_list, args, template_array)
    #generate some more useful constants
    template_seq_arrays = [template.get_template_seq_array() for template in template_array]
    # now process each read pair
    final_return_val = [[]] * len(read_seq_record_list)
    
    for read_id, read_df in blast_result_df.groupby('read_id'):
        identified_index_bcs = []
        identified_umi_sequences = []
        identified_bc_seq_records = []
        fwd_record = read_seq_record_list[read_id ][0]
        rev_record = read_seq_record_list[read_id][1]

        # we now step through each template seq (fwd and rev) in order and see if we can find a match
        for read_num, template_seq_array in enumerate(template_seq_arrays):
            read_subdf = read_df.loc[read_df['read_num']==read_num]
            read = fwd_record if read_num==0 else rev_record

            if read_subdf.empty:
                continue

            # first make an empty array for coordinates of features

            starting_coordinates, new_index_bcs = _calculate_feature_coordinates(template_seq_array, template_seq_length_dict, read_subdf, int(args['read_length']))
            identified_index_bcs.extend(new_index_bcs)

            # now that we have all the coordinates, let us extract the sequences for each template feature
            for i, template_seq_feature in enumerate(template_seq_array):
                # set the end coordinate of this feature properly
                if (starting_coordinates[i] is None or starting_coordinates[i+1] is None):
                    continue
                start = min(starting_coordinates[i], starting_coordinates[i+1])
                end = min(max(starting_coordinates[i], starting_coordinates[i+1]), len(read.seq))

                if template_seq_feature.is_umi():  # extract UMI sequences if any
                    umi_seq = read.seq[start:end]

                    if (len(umi_seq) > 0 and args['use_umi']):
                        identified_umi_sequences.append(umi_seq)
                # extract coordinates of any BC region that exist
                if template_seq_feature.is_barcode():
                    if (read_id == 0):
                        start = start + int(args['barcode_5prime_trim_length'])
                    if ((read_id == 0 and not args['is_paired_end']) or (read_id == 1 and args['is_paired_end'])):
                        end = end - int(args['barcode_3prime_trim_length'])
                    mybc = read[start:end]
                    if (len(mybc.seq) > 0):
                        identified_bc_seq_records.append(mybc)

        return_val = [identified_bc_seq_records, identified_umi_sequences,
                     identified_index_bcs]
        final_return_val[read_id] = return_val
    return final_return_val
