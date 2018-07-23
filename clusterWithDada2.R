library(dada2)
library(ggplot2)
library(Biostrings)
packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two arguments required: (input file).fastq and outputPrefix", call.=FALSE)
} 

#predefined values
minReadFilter = 10

# get input file and output file prefix
fastqFileName = args[1]
outputPrefix = args[2]
filterFastqFilePath = file.path(paste0(outputPrefix,"BCs_filterAndTrimmed.fastq.gz"))

#fitler and trim input!

filteredAndTrimmed = filterAndTrim(file.path(fastqFileName),filterFastqFilePath,maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=FALSE)

#learn errors, dereplicate, run dada and make output table
errF <- learnErrors(filterFastqFilePath, multithread=FALSE)
derepFs <- derepFastq(filterFastqFilePath, verbose=TRUE)
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
seqTab <- makeSequenceTable(dadaFs)

#dada2 has a merged paired reads step that does counting and denoising. we are going to skip this since we effectively have single end data

#not going to remove chimeras as this is likely true chimeras anyways.

#where are all of the reads going?
clusteredBCs = colnames(seqTab)
BCFasta = data.frame(x=paste0(">",seq(1:NROW(clusteredBCs))),seq=clusteredBCs)
		
write.table(file=paste0(outputPrefix,"clusteredBCs.fasta"),BCFasta,sep="\n",quote=F,col.names=F,row.names=F)
#write.table(file=paste0(outputPrefix,"readTrack.tab"),track,sep="\t",quote=F)