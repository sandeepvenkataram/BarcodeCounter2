import glob
import os
import re
from dataclasses import dataclass

from lib.constants import VALID_FEATURE_TYPES

###########################################################################
## Struct Definitions
###########################################################################


@dataclass
class TemplateSeqFeature:
    code: str
    sequence: str
    read_num: int
    feature_position: int
    fasta_feature_name: str

    @property
    def seq_length(self) -> int:
        '''Length of the sequence'''
        return len(self.sequence)

    def is_index(self):
        return self.code == 'D'

    def is_umi(self):
        return self.code == 'U'

    def is_barcode(self):
        return self.code == 'X'

    def is_constant(self):
        return self.code == self.sequence and not (self.is_index() or self.is_barcode() or self.is_umi())


class TemplateSequence:
    '''Representation of a template sequence for read splitting'''

    sequence: str
    CONST_REGION_CODE = 'CONST_REGION'

    def __init__(self, sequence, read_num: int):
        self.valid_special_characters = VALID_FEATURE_TYPES.keys()
        self.sequence = sequence
        self.read_num = read_num
        if len(self.sequence) == 0:
            raise ValueError('Expected a non-empty sequence!')
        if len(re.sub('[ACGTacgtXUDNn]', '', self.sequence)) > 0:
            raise ValueError('Template sequence has illegal characters!\n')

    def get_expected_barcode_length(self) -> int:
        '''Get the total length of the concatenated barcode regions in this template'''
        return self.sequence.count('X')

    def get_template_seq_array(self) -> list[TemplateSeqFeature]:
        '''Get an array representation of the features in this template sequence'''
        seq_split = self.__split_sequence_by_group()
        res = []
        for i, seq in enumerate(seq_split):
            code = self.__get_char_class(seq[0])
            res.append(
                TemplateSeqFeature(
                    code=code if code != self.CONST_REGION_CODE else seq,
                    sequence=seq,
                    read_num=self.read_num,
                    feature_position=i,
                    fasta_feature_name=f'const_region_{self.read_num}_{i}' if code == self.CONST_REGION_CODE else None,
                )
            )
        return res

    def __get_char_class(self, my_seq: str):
        if my_seq[0] in self.valid_special_characters:
            return my_seq[0]
        return self.CONST_REGION_CODE

    def __split_sequence_by_group(self) -> list[str]:
        cur_char_idex = 0
        seq_split = []
        for i, my_char in enumerate(self.sequence):
            if self.__get_char_class(self.sequence[cur_char_idex]) == self.__get_char_class(my_char):
                continue
            seq_split.append(self.sequence[cur_char_idex:i])
            cur_char_idex = i
        seq_split.append(self.sequence[cur_char_idex:])
        return seq_split

    def __str__(self):
        return f'TemplateSequence Object for seq: {self.sequence}'


class SampleFileMetadata:
    sample: str
    file_prefix: str
    int_multi_bc_array: list[str]

    def __init__(self, sample: str, file_prefix: str, int_multi_bc_array: list[str]):
        self.sample = sample
        self.file_prefix = file_prefix
        self.int_multi_bc_array = int_multi_bc_array
        self.fwd_fastq = None
        self.rev_fastq = None
        self.bc_fastq = None
        self.r1_fastq = None
        self.r2_fastq = None
        self.umi_tab = None
        self.resplit = False

    def identify_fastq_files(self, args: dict):
        pattern_string = re.compile('.*' + self.file_prefix + '.*')
        read_files = glob.glob(args['fastq_dir'] + '*')
        read_files_2 = list(filter(pattern_string.match, read_files))
        read_files_2.sort()
        my_fwd = None
        my_rev = None
        if len(read_files_2) == 0:
            raise ValueError('No matching fastq files found for ' + self.file_prefix)
        if len(read_files_2) > 2:
            raise ValueError('More than two matching fastq files found for ' + self.file_prefix)
        my_fwd = read_files_2[0]
        if args['is_paired_end']:
            my_rev = read_files_2[1]

        self.fwd_fastq = my_fwd
        self.rev_fastq = my_rev

    def touch_files(self, args: dict):
        self.bc_fastq = args['output_dir'] + self.sample + '_barcode.fastq'
        self.r1_fastq = args['output_dir'] + self.sample + '_R1.fastq'
        self.r2_fastq = args['output_dir'] + self.sample + '_R2.fastq' if len(self.int_multi_bc_array) == 2 and args['is_paired_end'] else None
        self.umi_tab = args['output_dir'] + self.sample + '_UMISeqs.tab'
        for file in [self.bc_fastq, self.r1_fastq, self.r2_fastq, self.umi_tab]:
            if file is None:
                continue
            if not os.path.exists(file) or args['resplit_fastq']:
                self.resplit = True
                with open(file, 'w', encoding='utf-8'):
                    pass

    def index_key(self):
        return f'{self.file_prefix}_{'_'.join([str(x) for x in self.int_multi_bc_array])}'

    def __str__(self):
        return f'{self.sample}, {self.file_prefix}, {'_'.join(self.int_multi_bc_array)}'

    def __hash__(self):
        return str(self).__hash__()


class IndexCounter:

    def __init__(self, sample_key: str):
        self.sample_key = sample_key
        self.good_reads = 0
        self.n_reads = 0
        self.missing_umi = 0
        self.wrong_num_bcs = 0
        self.no_matching_sample = 0

    def update(self, flag: int):
        if flag == 0:
            self.good_reads += 1
        elif flag == 1:
            self.n_reads += 1
        elif flag == 2:
            self.missing_umi += 1
        elif flag == 3:
            self.wrong_num_bcs += 1
        elif flag == 4:
            self.no_matching_sample += 1
        else:
            raise ValueError(f'Invalid flag {flag} for index counting update!')

    def to_dict(self) -> dict:
        return {
            'sample': self.sample_key,
            'good_reads': self.good_reads,
            'n_reads': self.n_reads,
            'missing_umi': self.missing_umi,
            'wrong_num_bcs': self.wrong_num_bcs,
            'no_matching_sample': self.no_matching_sample,
        }
