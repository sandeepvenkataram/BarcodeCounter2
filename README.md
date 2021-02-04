# BarcodeCounter2

This software is designed to generate barcode frequency counts from trimmed / quality controlled Illumina sequencing data (fastq files). It can accept arbitrary barcode locus designs, allowing for multiple barcodes, UMI sequences and inline multiplexing indices. It can process both single-end and paired-end PCR amplicon data. It assumes that adapters were added by PCR and not blunt end ligation since it assumes a constant orientation of all reads. It uses python multiprocessing for speed.


# Dependencies:

Python 3

BioPython

DNAClust

BLAST suite

Bowtie2 or BWA

# Execution command for sample data

python3 barcodeCounter.py -fastqDir SampleData/rawFastqFiles/ -outputDir OutputDir/ -templateSeq SampleData/sequenceTemplate.txt -sample SampleData/sampleFile.txt -multiBCFasta SampleData/primerIndexSeq.fasta -pairedEnd -useUMI -numThreads 3

For detailed help information, run the program without any argument or use the -h flag.

# Required input arguments

Users need to supply:
1. A directory containing fastq/fastq.gz files containing their sequencing data
2. A sequence template file, defining the expected amplicon sequence in their data. The file contains a single line with the sequence and no header information. The sequence template uses non-standard letters to represent the UMI, multiplexing index and barcode regions within the sequence. See help information and sample files for details. This sequence should be oriented such that it is in the same 5'-3' orientation as the first read of the sequencing data.
3. A sample file defining the samples to be processed. This is a tab-delimited filewith columns Sample Name\t File Prefix\t internal multiplexing index 1 name\t internal multiplexing index 2 name\t internal multiplexing index 3 name.... The file prefix is a string that uniquely identifies a fastq file or pair of fastq files (for paired end data) that will be used to get the sequences for this sample. Multiple samples can be associated with a single fastq file prefix, as long as there is a unique combination of internal multiplexing indices for each one. 
4. If using multiplexing indices (as defined by the sequence template), users must supply a fasta file of sequences defining their expected sequences. If these indices are very short (<10bp), include some of the constant sequence 3' or 5' of the index to ensure that BLAST can uniquely find these sequences within your reads. The names of these sequences should correspond to the names given in the sample file.
5. Users can supply a fasta file defining the list of expected barcode sequences rather than having the software infer them de novo using a clustering algorithm. This is very useful if the same library is being sequenced repeatedly from different experiments.

# Definitions

UMI = Universal molecular identifier. Usually a randomized part of a primer that uniquely tags each amplified molecule, used to remove PCR duplicates from inflating amplicon counts.

Multiplexing index - a custom sequence within each primer to distinguish samples that are pooled together using the same Illumina indices. I typically do a two-step PCR reaction for my amplicons, the first one adds inline indices and UMIs, and the second uses Illumina Nextera primers to make the amplicon compatible with the sequencing machine. I have a panel of forward and reverse 1st step primers with distinct multiplexing indices so that I can increase the number of samples I submit in each lane.

Barcode - a randomized sequence within the genome / plasmid that you are attempting to genotype by amplicon sequencing, typically for the purpose of estimating the frequency of the genotype in a pooled library.


# Tips to make the software run smoothly

Do not run multiple instances of this program simultaneously. There are temporary files that are made in the execution directory, and multiple executions may overwrite them, causing serious errors. Instead, make a single large sample file and use multithreading to speed up processing time.

Ensure that all of your files are consistent with each other. The sequence template should be consistent with the reads (5' of template = 5' of read 1), the names of the multiplexing indices in the fasta file should match the sample sheet, and the sample sheet should define fasta file prefixes that are unique and match the fasta file names in the directory. I have tried to put in checks for all of these, but I can't guarantee that the code is foolproof.

Make sure you run the program with the "-pairedEnd" option if you have paired-end data! Also check the read length and set the parameter appropriately to make sure your reads are not getting truncated.

# Tips for experimental design

When designing multiplexing indices, do not reuse sequences between different indexing locations in the same read, or use sequences that are simply the reverse complement of each other. The software uses BLAST to match these sequences, and such reuse would either cause misassignment of reads to samples or incorrectly discard valid reads. 

In our experience, amplicon sequences appear to "recombine" on the Illumina sequencing machines (as of 2019). This can be mitigated by 1) submitting a number of samples on the same lane for whole-genome sequencing, rather than amplicon sequencing, to dilute the availability of identical DNA fragments for recombination (ideally >50% of the lane should be whole-genome) and 2) by associating each Illumina multiplexing index N7xx / N5xx with a specific multiplexing inline index on the other side of the amplicon such that recombination events can always be detected.

Minimize the number of PCR cycles used to generate the amplicons. More cycles = more chances for sequence errors and jack-potting events causing individual barcodes to be over- or under-represented. If you need more template, run multiple independent reactions and then pool them together to minimize this effect. You should also use high-fidelity polymerases to reduce errors.


# BarcodeExtractor

This software is a modification of barcodeCounter, whose purpose is to extract barcodes from FASTA formatted sequence files. This software is primarily used to extract barcodes from Sanger sequencing of individual barcoded clones.

# Required input arguments

Users need to supply:
1. An input FASTA file
2. A sequence template, as defined in barcodeCounter
3. An output directory
