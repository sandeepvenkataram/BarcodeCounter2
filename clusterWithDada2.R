library(dada2)
library(ggplot2)
library(Biostrings)
packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments required: (input file).fastq, (errorSampleFile).fastq and outputPrefix", call.=FALSE)
} 

#predefined values
minReadFilter = 10

# get input file, file for making error model and output file prefix
fastqFileName = args[1]
fastqErrorSampleFileName = args[2]
outputPrefix = args[3]
filterFastqFilePath = file.path(paste0(outputPrefix,"BCs_filterAndTrimmed.fastq.gz"))
filterFastqErrorsFilePath = file.path(paste0(outputPrefix,"BCs_filterAndTrimmedForErrors.fastq.gz"))

#fitler and trim input!

filteredAndTrimmed = filterAndTrim(file.path(fastqFileName),filterFastqFilePath,maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE)
filteredAndTrimmedErrors = filterAndTrim(file.path(fastqErrorSampleFileName),filterFastqErrorsFilePath,maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE)

#learn errors, dereplicate, run dada and make output table
errF <- learnErrors(filterFastqErrorsFilePath, multithread=TRUE)
derepFs <- derepFastq(filterFastqFilePath, verbose=TRUE)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
seqTab <- makeSequenceTable(dadaFs)

#identify clustered barcodes and print to fasta output file
clusteredBCs = colnames(seqTab)
BCFasta = data.frame(x=paste0(">",seq(1:NROW(clusteredBCs))),seq=clusteredBCs)
		
write.table(file=paste0(outputPrefix,"clusteredBCs.fasta"),BCFasta,sep="\n",quote=F,col.names=F,row.names=F)
