# BarcodeCounter2

This software is designed to generate barcode frequency counts from trimmed / quality controlled Illumina sequencing data (fastq files). It can accept arbitrary barcode locus designs, allowing for multiple barcodes, UMI sequences and inline multiplexing indices. It can process both single-end and paired-end PCR amplicon data. It assumes that adapters were added by PCR and not blunt end ligation since it assumes a constant orientation of all reads. It uses python multiprocessing for speed.


# Dependencies:

Python 3

BioPython

Dada2 (R package, installed via bioconductor)

BLAST suite

# Execution command for sample data

python3 barcodeCounter.py -fastqDir SampleData/rawFastqFiles/ -outputDir OutputDir/ -templateSeq SampleData/sequenceTemplate.txt -sample SampleData/sampleFile.txt -multiBCFasta SampleData/primerIndexSeq.fasta -pairedEnd -useUMI -numThreads 3
