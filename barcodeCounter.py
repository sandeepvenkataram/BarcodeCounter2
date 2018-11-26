###########################################################################
# September 7th, 2018
# Sandeep Venkataram, PhD.
# Postdoctoral Scholar, Kryazhimskiy Lab, 
# UCSD Division of Biological Sciences
# barcodeCounter.py 
# 
# This software is designed to generate barcode frequency counts from 
# trimmed / quality controlled Illumina sequencing data (fastq files).
# It can accept arbitrary barcode locus designs, allowing for multiple
# barcodes, UMI sequences and inline multiplexing indices.
# It can process both single-end and paired-end PCR amplicon data. It assumes
# that adapters were added by PCR and not blunt end ligation since it assumes
# a constant orientation of all reads. It uses python multiprocessing for speed
# 
# 
# There is legacy code to do clustering with Dada2. As it is too slow for typical 
# barcode calling purposes, it has been commented out.
# 
# Dependencies:
#  Python 3
#  BioPython
#  DNAClust
#  BLAST suite
#  bowtie2
#  linux environment for shell scripting (grep and cut)
# 
# Execution command for sample data
#	python3 barcodeCounter.py -fastqDir ../BCCounterTesting/rawFastqFiles/ -outputDir ../BCCounterTesting/testBCCounterOutputDir/ -templateSeq ../BCCounterTesting/sequenceTemplate.txt -sample ../BCCounterTesting/sampleFile.txt -multiBCFasta ../BCCounterTesting/primerIndexSeq.fasta -pairedEnd -useUMI -numThreads 3
#
###########################################################################


import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from collections import namedtuple
import csv
#from __future__ import division
import fileinput
import glob
import gzip
import itertools as IT
from multiprocessing import Pool, Lock
import os
import random
import re
import sys
import subprocess





###########################################################################
## Struct Definitions
###########################################################################
  
sampleStruct = namedtuple("sampleStruct","Sample FilePrefix intMultiBCArray")
pairedFastqFileStruct = namedtuple("pairedFastqFileStruct","FwdFastq RevFastq MatchPrefix")

###########################################################################
## Global Variables
###########################################################################

sampleArray = []
indexToSampleMap = {}
templateSeqLengthsDict = {}
templateSeqArray = [] #one to two inner arrays, describing the template corresponding to each read
validFeatureTypes = {"X" : "BC", "D" : "multiplexBC", "U" : "UMI"}
validFeatureCounts = {"X" : 0, "D" : 0, "U" : 0}
usedFastqFiles = []
constantRegionFastaFilenames = [] # same dimensions as templateSeqArray
constantRegionNames = [] # same dimensions as templateSeqArray
allConstRegionsFileName = "allCRs.fasta"
maxNInReads = 3
fileBufferSize = 10000 #make sure we aren't constantly accessing the disk but small enough to reasonably keep in memory
readsPerSampleForErrors = 10000

#Blast parameters
constantRegionsBlastParams = ["-word_size", "6","-outfmt","6","-evalue","1E0"]

###########################################################################
## Import arguments and define help output
###########################################################################


cmdLineArgParser = argparse.ArgumentParser(description="Process Illumina PCR Amplicon Fastq files to count barcode frequencies. Must have blastn and dada2 (R Bioconductor  package) installed to run. Make sure there are no spaces in any file names or directory paths, the program will not work otherwise. Make sure that files are already split by their illumina indices (N700 / S500 index). This code will not function if there are multiple barcodes or internal multiplexing barcodes within a single read. Please use the readLength flag to solve this issue. \n\nGenerated Files:\n\nnumReadsFoundPerSample: One file per input fastq file, it details the assignment of each read to a combination of inline indices if they exist, and the subsequent filtering of reads. The file has 2 columns, the first being a unique identifier for the fastq file and inline index combination identified, and the second being an array of 6 numbers binning the reads into the following categories: \n\t1. correct reads used for subsequent mapping. \n\t2. There are too many Ns in the read.\n\t3. UMIs have been requested and not identified in the read/\n\t4. There is no barcode found in the read.\n\t5. The combination of indices do not match a sample specified in the sampleFile.\n\t6. A second check for finding a valid barcode in the read.\n\n For each sample in SampleFile with at least one barcode found a number of files are generated.\n R1(R2).fastq - raw (possibly truncated depending on the readLength parameter) reads associated with this sample\n barcode.fastq - the portion of the reads associated with all barcode sequences concatenated together, in the same order as the R1/R2.fastq files\n UMISeqs.tab - tab delmited UMI sequences if they exist and are being used, in the same order as the R1/R2.fastq files.\n readBarcodeID.txt - the ID number of the barcode in the clusteredBCs.fasta file that each read was mapped to, in the same order as the R1/R2.fastq files.\nbarcodeCalls.tab - Count of each barcode in the sample, removing UMI duplicates if asked for. Line 1 contains the counts for barcode 1 (defined in clusteredBCs.fasta), line 2 for barcode 2 etc.\n\nallBarcodeCalls.tab - tab delimited concatenation of all of the barcodeCalls.tab files, with a header row identifying which column comes from which sample.")
cmdLineArgParser.add_argument("-fastqDir", dest="fastqDir", help="directory location of fastq files",required=True)
cmdLineArgParser.add_argument("-outputDir", dest="outputDir", help="location of output directory",required=True)
cmdLineArgParser.add_argument("-templateSeq", dest="templateSeqFile", help="Template sequence of amplicon locus. This file contains a single line with standard DNA bases. UMI (unique molecular identifier) sequences are coded as U, multiplexing indices are coded as D and barcode loci coded as X. If these features have different lengths between samples, define the template using the longest possible length of each feature. Every feature annotated must be covered by the sequencing data, and no feature can span the exact middle of the template sequence when using paired end data.",required=True)
cmdLineArgParser.add_argument("-sample", dest="sampleFile", help="File defining the samples present in the sequencing data. This file is tab delimited, with no header line. The column values are: Sample Name\t File Prefix\t internal multiplexing barcode 1\t internal multiplexing barcode 2... The internal multiplexing barcode columns must correspond to the names of the sequences in the multiBCFasta file. Do not use spaces in any of the columns for file / directory names, as this tends to behave poorly.",required=True)
cmdLineArgParser.add_argument("-readLength", dest="readLength", default=100,  help="Expected length of each read from sequencing machine. Default = 100. Reduce this number from the true read length if necessary such that non-constant regions of the barcode locus are not shared between reads. This does not modify the input fastq files, but effectively truncates reads before processing")
cmdLineArgParser.add_argument("-barcodeList", dest="barcodeListFile", help="Optional fasta file specifying the barcodes present in the sample. If file is not supplied, unique barcodes will be identified de novo. The name for each sequence should just be a number identifying the barcode ID.")
cmdLineArgParser.add_argument("-blastPATH", dest="blastPATH", help="BLAST installation directory if it is not in the PATH already", default="")
cmdLineArgParser.add_argument("-bowtie2PATH", dest="bowtie2PATH", help="Bowtie2 installation directory if it is not in the PATH already", default="")
cmdLineArgParser.add_argument("-DNAclustPATH", dest="DNAclustPATH", help="DNAclust installation directory if it is not in the PATH already", default="")
cmdLineArgParser.add_argument("-multiBCFasta", dest="multiBCFastaFile", help="A multi-line fasta file defining multiplexing tag sequences. Required if there are multiplexing tags within the amplicon sequence as defined by the templateSeq file.")
cmdLineArgParser.add_argument("-pairedEnd", dest="pairedEnd", action="store_true",  help="Use if sequencing data is paired end")
#cmdLineArgParser.add_argument("-RscriptExecPATH", dest="RscriptPATH", help="Rscript executable path if it is not in the PATH already. Must have Dada2 software package installed.", default="Rscript")
cmdLineArgParser.add_argument("-useUMI", dest="UMI", action="store_true",  help="Use flag if you want to remove PCR duplicate reads using UMI data")
cmdLineArgParser.add_argument("-numThreads", dest="numThreads", default=1,  help="Number of threads to be used for computation.")
cmdLineArgParser.add_argument("-skipSplitFastq", dest="skipSplitFastq", action="store_true",  help="Use flag if you want to skip the splitting of the raw fastq files (i.e. if you have already done this and do not want to redo it).")
cmdLineArgParser.add_argument("-demultiplexOnly", dest="demultiplexOnly", action="store_true",  help="Use flag if you want to only split the raw fastq files and not continue with the rest of the barcode counting. This is useful when distributing demultiplexing across several machines, i.e. in a cluster.")
cmdLineArgParser.add_argument("-remapBarcodes", dest="remapBarcodes", action="store_true",  help="Set to True if you want to remap barcodes even if the files already exist")


#cmdLineArgParser.add_argument("-rebarcoding", dest="rebarcodingFile", help="File defining timepoints within each experiment when new barcodes were introduced. This feature is not currently implemented")

args = cmdLineArgParser.parse_args()

###########################################################################
## Utility functions
###########################################################################


## Initializer for lock for multiprocessing
  
def init(l):
	global lock
	lock = l

	
	
## Calculate the number of lines in the file

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


	
## Print to std error

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

	
	
## Return the blast hit with lowest e value. Requires input in outfmt6 format, output is an array with the columns separated as strings. 

def getBestBlastMatch(blastOutputLocal):
	bestMatch = []
	currentBestEValue = 1
	blastOutputArray = blastOutputLocal
	if(len(blastOutputArray)==0): #if there are no hits, return empty array
		return bestMatch
	splitInfo = blastOutputArray[0].split("\t")
	if(len(splitInfo)<10): #if the line is poorly formatted, return empty array
		return bestMatch
	currentBestEValue = float(splitInfo[10])
	bestMatch = [splitInfo[1]]
	bestMatch.extend(list(map(int,splitInfo[6:10])))
	if(len(blastOutputArray)==1): #if there is only one hit, return it
		return bestMatch
	splitInfo = blastOutputArray[1].split("\t")
	if(currentBestEValue < float(splitInfo[10])): #if the first hit is better than the second, return it
		return bestMatch
	return "multiple" #there are multiple hits!


	
## Extract coordinates of blast match dealing with reverse complemented sequences if necessary

def getSubjectMatchCoordinates(topBlastResult):
	topBlastResultInts = [topBlastResult[0]]
	for x in topBlastResult[1:len(topBlastResult)]:
		topBlastResultInts.append(int(x))
	startingCoor = topBlastResultInts[1]
	endingCoor = topBlastResultInts[2]
	rev = False
	if(topBlastResultInts[3]<topBlastResultInts[4]): #forward orientation
		startingCoor = startingCoor - topBlastResultInts[3]
		endingCoor = endingCoor + templateSeqLengthsDict[topBlastResultInts[0]] - topBlastResultInts[4]
	else:
		rev = True
		startingCoor = startingCoor - templateSeqLengthsDict[topBlastResultInts[0]] + topBlastResultInts[3]-1
		endingCoor = endingCoor + topBlastResultInts[4]-1
	return [startingCoor, endingCoor, rev]

	
	
## Extract constant regions in template array, creates fasta file and blast database for each, then make a master database with all constant regions and multiplexing primer sequences

def createConstRegionFasta():
	readNumber = 0
	filenames = []
	for seqArray in templateSeqArray:
		readNumber = readNumber + 1
		seqNumber = 0
		constRegionFilenameArray = []
		constRegionNameArray = []
		for seq in seqArray:
			if(len(seq)>1):
				seqNumber = seqNumber+1
				filePrefix = "const_region_"+str(readNumber)+"_"+str(seqNumber)
				fileString = "."+filePrefix+".fasta"
				with open(fileString,"w") as outfile:
					outfile.write(">const_region_"+str(readNumber)+"_"+str(seqNumber)+"\n")
					outfile.write(seq+"\n")
				constRegionFilenameArray.append(fileString)
				constRegionNameArray.append(filePrefix)
				filenames.append(fileString)
				subprocess.call([args.blastPATH+"makeblastdb","-in",fileString,"-dbtype","nucl"])
			else:
				constRegionFilenameArray.append(None)
				constRegionNameArray.append(None)
		global constantRegionFastaFilenames
		constantRegionFastaFilenames.append(constRegionFilenameArray)
		global constantRegionNames
		constantRegionNames.append(constRegionNameArray)
	
	if(args.multiBCFastaFile != None):
		filenames.append(args.multiBCFastaFile)
	with open(allConstRegionsFileName, 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				seqName = ""
				for line in infile:
					if line[0:1] == ">":
						seqName = line[1:(len(line)+1)]
					else:
						templateSeqLengthsDict[seqName.strip()]=len(line.strip())
					outfile.write(line)
	blastCall = [args.blastPATH+"makeblastdb","-in",allConstRegionsFileName,"-dbtype","nucl"]
	subprocess.call(blastCall)
	

	
## This gets UMI, multiplexing index and barcode regions from each read via blast

def extractRegionsFromFastq(readSeqRecordList, prefixName):
	# make a fasta file from all reads we are processing and blast against database of all index and constant regions
	readSeqFileName = "."+prefixName+"_readSeq.fasta"
	with open(readSeqFileName,"w") as outfasta:
		for readID in range(0,len(readSeqRecordList)):
			for i in range(0,len(templateSeqArray)):
				outfasta.write(">read_"+str(readID)+"_"+str(i)+"\n")
				outfasta.write(str(readSeqRecordList[readID][i].seq)+"\n")
	
	# blast reads against constant and multiplexing index sequences
	blastCommand = [args.blastPATH+"blastn","-query",readSeqFileName,"-db",allConstRegionsFileName]
	blastCommand.extend(constantRegionsBlastParams)
	blastResult = subprocess.check_output(blastCommand).decode('ascii').rstrip().split("\n") 
	finalReturnVal = []
	notConstRegionRegex = re.compile("^((?!const_region).)*$")
	readDict = {}
	
	for line in blastResult: #store all blast hits in a dictionary for quick lookup
		mySplit = line.split("\t")
		curSample = mySplit[0]
		if curSample not in readDict:
			readDict[curSample]=[]
		readDict[curSample].append(line)
	
	for readID in range(0,len(readSeqRecordList)): #for each read pair we are supposed to process
		identifiedIndexBCs = []
		identifiedUMISequences = []
		identifiedBCSeqRecords = []
		flagged = 0
		fwdRecord = readSeqRecordList[readID][0]
		revRecord = None
		if args.pairedEnd:
			revRecord = readSeqRecordList[readID][1]
		
		for readNum in range(0,len(templateSeqArray)): #for each read in the read set we are supposed to be processing
			# first make an empty array for coordinates of features
			startingCoordinates = [None]*(len(templateSeqArray[readNum])+1)	
			myReadKey = "read_"+str(readID)+"_"+str(readNum)
			if (myReadKey not in readDict):
				continue
			readBlastResult=readDict[myReadKey] #get the blast result lines corresponding to this read
			
			blastLocationDict = {}
	
			for line in readBlastResult: #store all blast hits in a dictionary for quick lookup
				mySplit = line.split("\t")
				curSample = mySplit[1]
				if curSample not in blastLocationDict:
					blastLocationDict[curSample]=[]
				blastLocationDict[curSample].append(line)
			
			prevSegmentStartingCoord = -1
			reversedRead = False
			firstIndex = -1
			for i in range(0,len(templateSeqArray[readNum])): #for each feature in the template
				if(templateSeqArray[readNum][i]=="D" or (templateSeqArray[readNum][i]!="U" and templateSeqArray[readNum][i]!="X")): #if this is an indexing barcode location or a constant region
					myLocalHits = []
					topBlastResult = ""
					if(templateSeqArray[readNum][i]=="D"): # if this is a multiplexing index region
						myLocalHits = list(filter(notConstRegionRegex.match, readBlastResult))
						topBlastResult = getBestBlastMatch(myLocalHits) #get the best hitting indexing region
					else: #this is a constant region
						if constantRegionNames[readNum][i] not in blastLocationDict:
							continue
						myLocalHits = blastLocationDict[constantRegionNames[readNum][i]] #get blast hits for this specific constant region
						topBlastResult = getBestBlastMatch(myLocalHits)
					
					if(topBlastResult == [] or topBlastResult == "multiple"): #if we don't have a good match for this region
						continue
					
					if(firstIndex == -1):
						firstIndex = i
					#figure out where in the read this index region is hitting, set coordinates and track the first region we are mapping to account for possibly needing to reverse the match location
					matchCoords = getSubjectMatchCoordinates(topBlastResult)
					if(prevSegmentStartingCoord > matchCoords[0] and not reversedRead):
						reversedRead = True
					
					if(prevSegmentStartingCoord == -1): #if this is the first segment we have found a hit for, figure out the right index where the segment is
						prevSegmentStartingCoord = matchCoords[0]
					if (not reversedRead): #if the read is not reversed, add in the coordinates properly
						startingCoordinates[i]=matchCoords[0]
						startingCoordinates[i+1]=matchCoords[1]
					else:
						if(firstIndex >=0): #if it is reversed and this is not the first segment we have hit, assume the previous segment we hit needs to have its indices swapped
							tmpVal = startingCoordinates[firstIndex]
							startingCoordinates[firstIndex] = startingCoordinates[firstIndex+1]
							startingCoordinates[firstIndex+1] = tmpVal
							firstIndex = -2
						startingCoordinates[i]=matchCoords[1]
						startingCoordinates[i+1]=matchCoords[0]
					
					if(templateSeqArray[readNum][i]=="D"):
						identifiedIndexBCs.append(topBlastResult[0])
			if(not reversedRead and startingCoordinates[0] == None): #put the 0 in the right place depending on if the read is reversed or not relative to the template if we haven't found it already
				startingCoordinates[0] = 0
			if(reversedRead and startingCoordinates[len(startingCoordinates)-1] == None):
				startingCoordinates[len(startingCoordinates)-1] = 0
				
			if(reversedRead and startingCoordinates[0] == None): #set the end of the read to be readLength in the right place depending on if the read is reversed or not relative to the template and if we haven't found it already
				startingCoordinates[0] = int(args.readLength)
			if(not reversedRead and startingCoordinates[len(startingCoordinates)-1] == None):
				startingCoordinates[len(startingCoordinates)-1] = int(args.readLength)
			
			
			for i in range(0,len(templateSeqArray[readNum])): #now that we have all the coordinates, let us extract the sequences for each template feature
				#set the end coordinate of this feature properly
				if(startingCoordinates[i] == None or startingCoordinates[i+1] == None):
					continue
				start = min(startingCoordinates[i], startingCoordinates[i+1]) 
				end = max(startingCoordinates[i], startingCoordinates[i+1])
				maxLen = len(readSeqRecordList[readID][readNum].seq)
				if(end > maxLen):
					end = maxLen
					
				if(templateSeqArray[readNum][i]=="U"):#extract UMI sequences if any			
					UMIseq = readSeqRecordList[readID][readNum].seq[start:end]
					identifiedUMISequences.append(UMIseq)
					if(len(UMIseq)==0 and args.UMI):
						flagged = 2
				if(templateSeqArray[readNum][i]=="X"):#extract coordinates of any BC region that exist
					mybc = readSeqRecordList[readID][readNum][start:end]
					identifiedBCSeqRecords.append(mybc)
					if(len(mybc.seq)==0):
						flagged = 1
			
				if(flagged > 0 and 0==1): #print statements for bad reads, turned off for now
					lock.acquire()
					if(flagged == 1):
						print("BC has seq of 0 length!")
					if(flagged == 2):
						print("UMI has seq of 0 length!")
					print("Read num is: ", readNum)
					print(templateSeqArray[readNum])
					print(readSeqRecordList[readID][readNum].seq)
					print(start)
					print(end)
					print(startingCoordinates)
					print("\n".join(readBlastResult)+"\n")
					lock.release()
					sys.exit(0)
			if(1==0 and readID < 20):
				print("Processing a read from "+prefixName+"\t"+str(readID)+"\t"+str(readNum))
				print("\n".join(readBlastResult))
				print(templateSeqArray[readNum])
				print(readSeqRecordList[readID][readNum].seq+"\n"+str(startingCoordinates))
				print("UMIS:\n")
				print("\t".join(str(x) for x in identifiedUMISequences))
				print("BCs:\n")
				for rec in identifiedBCSeqRecords:
					print(rec)
				print(str(identifiedIndexBCs))
		returnVal = [identifiedBCSeqRecords,identifiedUMISequences,identifiedIndexBCs, fwdRecord, revRecord]
		finalReturnVal.append(returnVal)
	return(finalReturnVal)
	
	
	
	
###########################################################################	
## Functions for Parsing Input Files and Directories
###########################################################################


## read sample file provided by user and store in global array and map

def parseSampleFile():
	global sampleArray
	global indexToSampleMap
	numInlineIndices = 0
	for seqAr in templateSeqArray:
		for seq in seqAr:
			if(seq == "D"):
				numInlineIndices+=1
	with open(args.sampleFile) as infile: # go through every line in the sample file
		for line in infile:
			if(line.strip() == ""):
				continue
			lineSplit = line.strip().split("\t") #separate values in each line of the sample file
			if(len(lineSplit) != 2 + numInlineIndices): # if the line has not enough columns, error and quit
				eprint("Sample file line:\n"+line+"has incorrect number of columns!\n")
				eprint(str(len(lineSplit)))
				eprint(str(numInlineIndices))
				sys.exit(0)
			sampleIndexArray = "_".join(lineSplit[1:(len(lineSplit)+1)])
			sampleIdentityArray = lineSplit[0]
			
			if sampleIndexArray in indexToSampleMap: #if this particular file / index combination has been seen already, quit
				eprint("Sample file line:\n"+line+"has been assigned to a sample already!\n")
				sys.exit(0)
			
			#else catalog the sample
			
			sampleArray.append(sampleStruct(lineSplit[0],lineSplit[1],lineSplit[2:(len(lineSplit)+1)]))
			indexToSampleMap[sampleIndexArray] = sampleIdentityArray


## Read Template sequence provided by user and parse. 
#  This code works best if every UMI/Index/Barcode is only in a single read, so reads are effectively non-overlapping. 
#  Use -readLength to truncate the read if necessary			
			
def parseTemplateSeq():
	sequence = ""
	with open(args.templateSeqFile) as infile:
		sequence = infile.readline().strip()
	if(len(re.sub("[ACGTacgtXUDNn]","",sequence))>0): # if the sequence has disallowed characters, quit
		eprint("Template sequence has illegal characters!\n")
		sys.exit(0)
	if(args.multiBCFastaFile=="" and sequence.count("D")>0): #if there is no multiplexing file but we find multiplexing loci in the template, quit
		eprint("Multiplexing indices exist but no multiBCFastaFile provided!\n")
		sys.exit(0)
	#if it is paired end, split template using readSize, 
	templateArray = []
	if(args.pairedEnd):
		templateArray.append(sequence[0:int(args.readLength)])
		templateArray.append(sequence[int(len(sequence)-args.readLength):(len(sequence)+1)])
	else:
		templateArray.append(sequence)
	#collapse non-constant features into single character and add to global variable
	global templateSeqArray
	for seq in templateArray:
		seq = re.sub("D+","\tD\t",seq)
		seq = re.sub("U+","\tU\t",seq)
		seq = re.sub("X+","\tX\t",seq)
		seq = re.sub("\t\t","\t",seq)
		seq = re.sub("\t\t","\t",seq)
		seq = re.sub("\t\t","\t",seq)
		seq = re.sub("^\t","",seq)
		seq = re.sub("\t$","",seq)
		templateSeqArray.append(seq.split("\t"))

		
## Figure out which fastq files we are going to use, 
#  Pair forward and reverse reads as necessary, assume they are in alphabetical order (R1 before R2)

def identifyUsedFastqFiles():
	usedFastqFileDict = {}
	for sample in sampleArray:
		patternString = re.compile(".*"+sample.FilePrefix+".*.fastq(.gz)?")
		readFiles = glob.glob(args.fastqDir+"*")
		readFiles2 = list(filter(patternString.match,readFiles))
		readFiles2.sort()
		myfwd = None
		myrev = None
		if(len(readFiles2)==0 or (len(readFiles2)==1 and args.pairedEnd) or (len(readFiles2)==2 and not args.pairedEnd) or len(readFiles2)>2): #check that we found the expected number of fastq files (based on single or paired end data expected)
			eprint("Incorrect number of matching fastq files found for "+sample.FilePrefix)
			sys.exit(1)
		myfwd = readFiles2[0]
		if(args.pairedEnd):
			myrev = readFiles2[1]
		usedFastqFileDict[pairedFastqFileStruct(FwdFastq=myfwd,RevFastq=myrev,MatchPrefix=sample.FilePrefix)]=1
	global usedFastqFiles
	usedFastqFiles = usedFastqFileDict.keys()



###########################################################################
## Main Executable Functions
########################################################################### 


## Separate reads by inline and illumina indices, separate UMI and BC regions for each read as well. 
#  This is multiprocessed code, but is probably still the slowest part of the analysis
  

def demultiplexFastq(fastqPair):
	# outputs: barcode file fastq, total split fastq, umi sequences, unidentified reads
	indexCounter = {}
	badFwdReadsHandle = open(args.outputDir+fastqPair.MatchPrefix+"_unmappedReads_R1.fastq","w")
	
	fwdFastqHandle = open(fastqPair.FwdFastq)
	if(fastqPair.FwdFastq.endswith("gz")):
		fwdFastqHandle = gzip.open(fastqPair.FwdFastq, "rt")
	readCounter = 0
	readList = []
	if(args.pairedEnd):
		print("Processing new FASTQ file[s]!\nFWD file:\t"+fastqPair.FwdFastq+"\nREV file:\t"+fastqPair.RevFastq+"\n"+fastqPair.MatchPrefix+"\n")
		
		badRevReadsHandle = open(args.outputDir+fastqPair.MatchPrefix+"_unmappedReads_R2.fastq","w")
		revFastqHandle = open(fastqPair.RevFastq)
		if(fastqPair.RevFastq.endswith("gz")):
			revFastqHandle = gzip.open(fastqPair.RevFastq, "rt")
		
		for fwdRec, revRec in zip(SeqIO.parse(fwdFastqHandle,"fastq"), SeqIO.parse(revFastqHandle,"fastq")): #for each read, WE MAY WANT TO REPLACE THIS WITH GENERALFASTQITERATOR FOR PERFORMANCE! Zip should have OK performance, so that shouldn't be an issue
			readCounter = readCounter + 1
			fwdRec = fwdRec[0:args.readLength]
			revRec = revRec[0:args.readLength]
			readList.append([fwdRec,revRec])
			if(readCounter % fileBufferSize != 0): #if we are not at file buffer size, do not run parser since file io is super expensive and we want to minimize it
				continue
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, badFwdReadsHandle, badRevReadsHandle) #process buffered reads using helper method
			readList = []
				
		
		#If we are here, the file has finished parsing, now need to empty the buffer for the last time, redo what we just did.
		if(len(readList)>0):
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, badFwdReadsHandle, badRevReadsHandle)

		#close file handles!
		badRevReadsHandle.close()
		revFastqHandle.close()

	else: #same as above, but for single read sequencing. 
		#output file handles for unidentified reads
		print("Processing new FASTQ file[s]!\n"+fastqPair.FwdFastq+"\n"+fastqPair.MatchPrefix+"\n")
		
		for fwdRec in SeqIO.parse(fwdFastqHandle,"fastq"): 
			fwdRec = fwdRec[0:args.readLength]
			readCounter = readCounter + 1
			readList.append([fwdRec])#.reverse_complement(name=True,description=True,id=True,annotations=True)]) #reverse complement the reverse read to match with the template sequence
			if(readCounter % fileBufferSize != 0): #if we are not at file buffer size, do not run parser since file io is super expensive and we want to minimize it
				continue
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, badFwdReadsHandle, None) #process buffered reads using helper method
			readList = []
				
		
		#If we are here, the file has finished parsing, now need to empty the buffer for the last time, redo what we just did.
		if(len(readList)>0):
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, badFwdReadsHandle, None)
	
	##
	##finalizing for both single and paired end data, count total read statistics here
	##
	
	indexCounterSum = [0, 0, 0, 0, 0]	
	indexCounterHandle = open(args.outputDir+fastqPair.MatchPrefix+"_numReadsFoundPerSample.txt","w")
	for index in indexCounter.keys():
		indexCounterHandle.write(index+"\t"+str(indexCounter[index])+"\n")
		for i2 in range(0,5):
			indexCounterSum[i2] = indexCounterSum[i2] + indexCounter[index][i2]
	indexCounterHandle.write("Total\t"+str(indexCounterSum)+"\n")
	indexCounterHandle.close()
	badFwdReadsHandle.close()
	fwdFastqHandle.close()

	
## Helper function that processes a single read during demultiplexing
	
def demultiplexFastqHelper(readList, fastqPair, indexCounter, badFwdReadsHandle, badRevReadsHandle):
	
	extractedRegions = extractRegionsFromFastq(readList, fastqPair.MatchPrefix)
	#these are all lists (or dictionaries of lists) of output so that we don't have to do too many seqio write calls since they are slow.
	badFwdReadsList = []
	badRevReadsList = []
	finalBCRecordList = {}
	fwdRecordList = {}
	revRecordList = {}
	UMIList = {}
	
	expectedNumBCs = 0
	for subarray in templateSeqArray:
		for feature in subarray:
			if(feature == "X"):
				expectedNumBCs = expectedNumBCs+1
	for readID in range(0,len(extractedRegions)): #for each read we are processing
		fwdRecord = extractedRegions[readID][3]
		revRecord = extractedRegions[readID][4]
		totalNCount = str(fwdRecord.seq).count('N')
		if args.pairedEnd:
		    totalNCount = totalNCount  + str(revRecord.seq).count('N')
			
		#get the sample from the multiplexing barcodes
		myIndex = [fastqPair.MatchPrefix]
		myIndex.extend(extractedRegions[readID][2])
		myIndex = "_".join(myIndex) #this is a string that uniquely defines each sample that was multiplexed
		if myIndex not in indexCounter:
			indexCounter[myIndex] = [0, 0, 0, 0, 0]
		
		flag = 0
		if( totalNCount > maxNInReads): #if there are too many Ns in the read
			flag = 1
		if (args.UMI and len(extractedRegions[readID][1])==0): #if there is no UMI and UMI is expected
			flag = 2
		if (len(extractedRegions[readID][0]) != expectedNumBCs): #if there is not the expected number of barcodes
			flag = 3
		
		if(myIndex not in indexToSampleMap): #if we couldn't find a sample associated with this multiplexing barcode combination
			flag = 4
		if(flag == 0): #the read matches a valid sample
			mySample = indexToSampleMap[myIndex]
			if mySample not in finalBCRecordList:
				finalBCRecordList[mySample]=[]
				fwdRecordList[mySample]=[]
				revRecordList[mySample]=[]
				UMIList[mySample] = []
			
			#concat all BCs associated with this read to get final BC
			finalBC = None
			for bcRecord in extractedRegions[readID][0]:
				if(type(bcRecord).__name__ == "NoneType"):
					continue
				if(type(finalBC).__name__ == "NoneType"):
					finalBC = bcRecord
				else:
					finalBC = finalBC + bcRecord
			
			UMIList[mySample].append("\t".join(str(x) for x in extractedRegions[readID][1]))
			finalBCRecordList[mySample].append(finalBC)
			fwdRecordList[mySample].append(fwdRecord)
			if args.pairedEnd:
				revRecordList[mySample].append(revRecord)
		indexCounter[myIndex][flag] = indexCounter[myIndex][flag] + 1
		if(flag > 0):
			badFwdReadsList.append(fwdRecord)
			if args.pairedEnd:
				badRevReadsList.append(revRecord)
	#write outputs to files
	SeqIO.write(badFwdReadsList,badFwdReadsHandle,"fastq")
	if args.pairedEnd:
		SeqIO.write(badRevReadsList,badRevReadsHandle,"fastq")
	for mySample in finalBCRecordList.keys():
		barcodeFastqFileHandle = open(args.outputDir+mySample+"_barcode.fastq","a+")
		SeqIO.write(finalBCRecordList[mySample],barcodeFastqFileHandle,"fastq") #write BC portion of read to a fastq (concat among all BC features)
		barcodeFastqFileHandle.close()
		
		SampleSplitFastqFileHandleFWD = open(args.outputDir+mySample+"_R1.fastq","a+")
		SeqIO.write(fwdRecordList[mySample],SampleSplitFastqFileHandleFWD,"fastq") #write raw fwd read to sample read file
		SampleSplitFastqFileHandleFWD.close()
		
		UMIFileHandle = open(args.outputDir+mySample+"_UMISeqs.tab","a+")
		UMIFileHandle.write("\n".join(UMIList[mySample]))
		UMIFileHandle.close()
		
		if args.pairedEnd:
			SampleSplitFastqFileHandleREV = open(args.outputDir+mySample+"_R2.fastq","a+")
			SeqIO.write(revRecordList[mySample],SampleSplitFastqFileHandleREV,"fastq") #write raw rev read to sample read file
			SampleSplitFastqFileHandleREV.close()
	return indexCounter

	
## Do clustering across all samples using Dada2. 
#  WARNING: This is really slow, and probably won't work with large numbers of unique reads
#  This is currently not available without modifying the code
  
# def clusterBarcodesDada2():
	# ##
	# ## concat all barcode fastq files by experiment into a single file for clustering
	# ##
	# allFiles = glob.glob(args.outputDir+"*_barcode.fastq")
	# totalNumLines = 0
	# with open(args.outputDir+"allSamplesConcat.fastq","w") as outfile:
		# for fname in allFiles:
			# with open(fname) as infile:
				# for line in infile:
					# outfile.write(line)
					# totalNumLines += 1
	# totalNumReads = int(totalNumLines / 4)
	# #subsample 200k reads to make the dada2 error model, since this takes forever if you use the entire dataset
	# #code from https://pythonforbiologists.com/randomly-sampling-reads-from-a-fastq-file/
	# numberToSample = 200000
	# recordsToKeep = set(random.sample(range(totalNumReads + 1), numberToSample))
	# readNum = 0
	# with open(args.outputDir+"allSamplesConcat.fastq") as input, open(args.outputDir+"allSamplesConcatForErrors.fastq", "w") as output:
		# for line1 in input:
			# line2 = input.readline()
			# line3 = input.readline()
			# line4 = input.readline()
			# if readNum in recordsToKeep:
				# output.write(line1)
				# output.write(line2)
				# output.write(line3)
				# output.write(line4)
			# readNum += 1
	
	# #do the clustering with dada2 via r script, generates a fasta file of barcodes in output directory
	# callFunction = [args.RscriptPATH,os.path.dirname(sys.argv[0]) + "/clusterWithDada2.R",args.outputDir+"allSamplesConcat.fastq",args.outputDir+"allSamplesConcatForErrors.fastq", args.outputDir]
	# subprocess.call(callFunction)
	# args.barcodeListFile = args.outputDir+"clusteredBCsDada2.fasta"

	
## Do clustering across all samples using DNAClust. This is the default. 

	
def clusterBarcodesDNAClust():
	##
	## concat all barcode fastq files by experiment into a single file for clustering, remove those sequences that appear less than 3 times
	##
	allFiles = glob.glob(args.outputDir+"*_barcode.fastq")
	dedupFileName = args.outputDir+"allSamplesConcatDedup.fasta"
	readCountFileName = args.outputDir+"allSamplesConcatDedup.readCounts"
	dnaclustOutputFileName = args.outputDir+"allSamplesConcatDedup.dnaclustOut"
	clusteredBCFileName = args.outputDir+"clusteredBCsDNAClust.fasta"
	
	totalNumLines = 0
	uniqueBCLines = {}
	for fname in allFiles:
		with open(fname) as infileHandle:
			for line in SeqIO.parse(infileHandle,"fastq"):
				myseq = str(line.seq)
				if (myseq not in uniqueBCLines):
					uniqueBCLines[myseq] = 0
				uniqueBCLines[myseq] += 1
	readCounter = 1
	BCClusterListMap = {}
	readCounterMap = {}
	with open(dedupFileName,"w") as outfileHandle, open(readCountFileName,"w") as readCountFileHandle:
		for line in uniqueBCLines.keys():
			if 'N' not in line and uniqueBCLines[line] > 3: #remove any reads with Ns in it (< .5% of reads) or sequences with too few reads
				outfileHandle.write(">"+str(readCounter)+"\n")
				outfileHandle.write(line+"\n")
				BCClusterListMap[readCounter] = []
				readCounterMap[readCounter] = uniqueBCLines[line]
				readCountFileHandle.write(str(uniqueBCLines[line])+"\n")
				readCounter +=1
			
	# use DNAclust to cluster reads
	callFunction = [args.DNAclustPATH+"dnaclust", "-s",".95","--approximate-filter","-k","6","-t",str(args.numThreads),"-i",dedupFileName]
	with open(dnaclustOutputFileName,"w") as outfile:
		subprocess.call(callFunction,shell=True,stdout=outfile)
	
	# create final barcode fasta file using the centers of the clusters found by DNAclust
	bcsToUse = []
	with open(dnaclustOutputFileName,"r") as infile:
		for line in infile:
			line = line.strip()
			lineSplit = line.split("\t")
			bcsToUse.append(int(lineSplit[0]))

	bcsToUse.sort()
	totalNumBCs = len(bcsToUse)
	maxBCToUseIndex = bcsToUse[totalNumBCs-1]
	bcsToUseIndex = 0
	curIndex = 1
	with open(dedupFileName,"r") as infile, open(clusteredBCFileName,"w") as outfile:
		for headerline in infile:
			if(bcsToUseIndex > totalNumBCs - 1):
				break
			headerline = line.strip()
			seqline = infile.readline()
			seqline = seqline.strip()
			if(curIndex == bcsToUse[bcsToUseIndex]):
				outfile.write(">"+str(bcsToUseIndex+1)+"\n")
				outfile.write(seqline+"\n")
				bcsToUseIndex +=1
			curIndex +=1	
	
	args.barcodeListFile = clusteredBCFileName


## Map barcodes using bowtie2, multiprocessed code
  
def mapBarcodes(mySamp):
	#only run on this sample if the output file doesn't exist or flag has been set
	indexString = mySamp
	if (os.path.isfile(args.outputDir+indexString+"_barcode.fastq") and (not os.path.isfile(args.outputDir+indexString+"_readBarcodeID_bowtie2.txt") or args.remapBarcodes)):
	
		#make a dictionary to map barcode names in the barcode list fasta file to consecutive numbers for indexing in a vector.
		BCNameToIdxDict = {}
		IdxToBCNameDict = {}
		with open(args.barcodeListFile,"r") as infile:
			counter = 1
			for record in SeqIO.parse(infile,"fasta"):
				if(record.id in BCNameToIdxDict): #we have found a duplicate entry in the barcode list. quit!
					eprint("Duplicate entry "+record.id+" found in input barcode list!")
					sys.exit(1)
				BCNameToIdxDict[record.id]=counter
				counter +=1
		
		bcFastqFile = args.outputDir+indexString+"_barcode.fastq"
		bcSamFile = args.outputDir+indexString+"_barcode.sam"
		bcIDFile = args.outputDir+indexString+"_readBarcodeID_bowtie2.txt"
		mapQualFile = args.outputDir+indexString+"_readMappingQuality_bowtie2.txt"
		#bowtie2 call
		subprocess.call([args.bowtie2PATH+"bowtie2","-L 10","-q","-x "+args.barcodeListFile,"-U"+bcFastqFile,"-S"+bcSamFile])
		
		#get the barcode match for each read and put into a single column output file. do the same for mapping quality
		with open(bcIDFile,"w") as outfile:
			subprocess.call("grep -v '^@' "+bcSamFile+" | cut -f 3",stdout=outfile, shell=True)
		with open(mapQualFile,"w") as outfile:
			subprocess.call("grep -v '^@' "+bcSamFile+" | cut -f 5",stdout=outfile, shell=True)
		
		BCUMIMap = {}
		BCCountList = [0]*int(file_len(args.barcodeListFile)/2)
		BCUMIDupCountList = [0]*int(file_len(args.barcodeListFile)/2)
		totalUnmappedReads = 0
		bcIDFileHandle = open(bcIDFile,"r")
		mapQFileHandle = open(mapQualFile,"r")
		UMIFileHandle = open(args.outputDir+indexString+"_UMISeqs.tab","r")
		
		#for each read bc / umi pair
		for bcID, UMIstring, mapQ in zip(bcIDFileHandle, UMIFileHandle, mapQFileHandle):
			bcID = BCNameToIdxDict[bcID.strip()] #get the internal index corresponding to the matched barcode
			mapQ = mapQ.strip()
			if(bcID == "*" or bcID == "" or mapQ == "*" or int(mapQ) <= 20):
				totalUnmappedReads = totalUnmappedReads + 1
				continue
			else:
				mykey = bcID+"\t"+UMIstring
				if mykey not in BCUMIMap or not args.UMI:
					BCUMIMap[mykey]=1
					BCCountList[int(bcID)-1] = BCCountList[int(bcID)-1] + 1
				else:
					BCUMIDupCountList[int(bcID)-1] = BCUMIDupCountList[int(bcID)-1] + 1
		# the file is parsed, close filehandles
		bcIDFileHandle.close()
		UMIFileHandle.close()
		
		#write total count data to file
		with open(args.outputDir+indexString+"_barcodeCounts_bowtie2.tab","w") as outFileHandle:
			for countVal in BCCountList:
				outFileHandle.write(str(countVal)+"\n")
		with open(args.outputDir+indexString+"_numUnmappedReads_bowtie2.txt","w") as outFileHandle:
			outFileHandle.write(str(totalUnmappedReads))
		if(args.UMI):
			with open(args.outputDir+indexString+"_barcodeUMIDupCounts_bowtie2.tab","w") as outFileHandle:
				for countVal in BCUMIDupCountList:
					outFileHandle.write(str(countVal)+"\n")


## Create final concatenated table of read counts for every population/timepoint combination for every barcode

def generateFinalTables():
	#print barcode counts as giant tab delimited table, with 1st column as barcode ID number and header file being the sample each column comes from
	patternString = args.outputDir+"*_barcodeCounts*.tab"
	filenames = glob.glob(patternString)
	handles = [open(filename, 'r') for filename in filenames]    
	readers = [csv.reader(f, delimiter=',') for f in handles]
	filenames2 = ["BCID"]
	filenames2.extend(filenames)
	with  open(args.outputDir+"allBarcodeCounts.tab", 'w') as h:
		writer = csv.writer(h, delimiter='\t', lineterminator='\n', )
		writer.writerow(filenames2)
		i = 1
		for rows in IT.zip_longest(*readers, fillvalue=['']*2):
			combined_row = [i]
			for row in rows:
				row = row[:1] # select the columns you want
				if len(row) == 1:
					combined_row.extend(row)
				else:
					combined.extend([''])
			writer.writerow(combined_row)
			i = i + 1
	for f in handles:
		f.close()

		
		
		
###########################################################################
## Run Program
###########################################################################



## Call Initialization functions to parse user inputs

parseTemplateSeq()
parseSampleFile()
identifyUsedFastqFiles()
createConstRegionFasta()


## Demultiplex data using multiprocessing

l = Lock()
if(not args.skipSplitFastq):
	for mySamp in sampleArray: #make empty files for appending later on, so that we do not accidentally append reads into existing files.
		mySample = mySamp.Sample
		barcodeFastqFileHandle = open(args.outputDir+mySample+"_barcode.fastq","w")
		SampleSplitFastqFileHandleFWD = open(args.outputDir+mySample+"_R1.fastq","w")
		barcodeFastqFileHandle.close()
		SampleSplitFastqFileHandleFWD.close()
		UMIFileHandle = open(args.outputDir+mySample+"_UMISeqs.tab","w")
		UMIFileHandle.close()
		if(args.pairedEnd):
			SampleSplitFastqFileHandleREV = open(args.outputDir+mySample+"_R2.fastq","w")
			SampleSplitFastqFileHandleREV.close()
			
	with Pool(processes = int(args.numThreads), initializer = init, initargs = (l,)) as pool:
		pool.map(demultiplexFastq, usedFastqFiles)

		
## If we are only trying to demultiplex, quit now
		
if args.demultiplexOnly:
	print("Terminating after splitting raw fastq files as requested.")
	sys.exit(0)


## Cluster barcodes using DNAClust.
	
if args.barcodeListFile==None:
	clusterBarcodesDNAClust()


## Make database from barcode fasta file for mapping
  
subprocess.call([args.bowtie2PATH+"bowtie2-build",args.barcodeListFile,args.barcodeListFile])


## Map barcodes using bowtie2 with multiprocessing

uniqueSamples = {}
for sample in sampleArray:
	uniqueSamples[sample.Sample] = 1

with Pool(processes = int(args.numThreads)) as pool:
	pool.map(mapBarcodes, uniqueSamples.keys())


## Generate final output table	

generateFinalTables()
