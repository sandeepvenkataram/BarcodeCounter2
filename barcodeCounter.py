###########################################################################
# July 17th, 2018
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
# Dependencies:
#  Python 3
#  BioPython
#  Dada2 (R package, installed via bioconductor)
#  BLAST suite
#  bowtie2
#  linux environment for shell scripting (grep and cut)
# 
# Execution command for sample data
#	python3 barcodeCounter.py -fastqDir ../BCCounterTesting/rawFastqFiles/ -outputDir ../BCCounterTesting/testBCCounterOutputDir/ -templateSeq ../BCCounterTesting/sequenceTemplate.txt -sample ../BCCounterTesting/sampleFile.txt -multiBCFasta ../BCCounterTesting/primerIndexSeq.fasta -pairedEnd -useUMI -numThreads 3
###########################################################################



import argparse
from collections import namedtuple
import fileinput
import glob
import re
import sys
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import csv
import itertools as IT
import timeit
import time
import os
from multiprocessing import Pool, Lock



###########################################################################
## Struct Definitions
###########################################################################
  
sampleStruct = namedtuple("sampleStruct","Sample FilePrefix intMultiBCArray")
pairedFastqFileStruct = namedtuple("pairedFastqFileStruct","FwdFastq RevFastq MatchPrefix")

###########################################################################
## Global Variables
###########################################################################

sampleArray = []
sampleToIndexMap = {}
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
multiplexBCsBlastParams =  ["-word_size", "6","-outfmt","6","-evalue","1E-2"]
constantRegionsBlastParams = ["-word_size", "6","-outfmt","6","-evalue","1E0"]
BCRegionsBlastParams = ["-word_size", "6","-outfmt","6","-evalue","1E-4"]

###########################################################################
## Import arguments and define help output
###########################################################################


cmdLineArgParser = argparse.ArgumentParser(description="Process Illumina PCR Amplicon Fastq files to count barcode frequencies. Must have blastn and dada2 (R Bioconductor  package) installed to run. Make sure there are no spaces in any file names or directory paths, the program will not work otherwise.")
cmdLineArgParser.add_argument("-fastqDir", dest="fastqDir", help="directory location of fastq files",required=True)
cmdLineArgParser.add_argument("-outputDir", dest="outputDir", help="location of output directory",required=True)
cmdLineArgParser.add_argument("-templateSeq", dest="templateSeqFile", help="Template sequence of amplicon locus. This file contains a single line with standard DNA bases. UMI (unique molecular identifier) sequences are coded as U, multiplexing indices are coded as D and barcode loci coded as X. If these features have different lengths between samples, define the template using the longest possible length of each feature. Every feature annotated must be covered by the sequencing data, and no feature can span the exact middle of the template sequence when using paired end data.",required=True)
cmdLineArgParser.add_argument("-sample", dest="sampleFile", help="File defining the samples present in the sequencing data. This file is tab delimited, with no header line. The column values are: Sample Name\t File Prefix\t internal multiplexing barcode 1\t internal multiplexing barcode 2... The internal multiplexing barcode columns must correspond to the names of the sequences in the multiBCFasta file. Do not use spaces or underscores in any of the columns",required=True)
cmdLineArgParser.add_argument("-readLength", dest="readLength", default=100,  help="Expected length of each read from sequencing machine. Default = 100. Reduce this number from the true read length if necessary such that non-constant regions of the barcode locus are not shared between reads. This does not modify the actual reads")
cmdLineArgParser.add_argument("-barcodeList", dest="barcodeListFile", help="Optional fasta file specifying the barcodes present in the sample. If file is not supplied, unique barcodes will be identified de novo. The name for each sequence should just be a number identifying the barcode ID.")
cmdLineArgParser.add_argument("-blastPATH", dest="blastPATH", help="BLAST installation directory if it is not in the PATH already", default="")
cmdLineArgParser.add_argument("-bowtie2PATH", dest="bowtie2PATH", help="Bowtie2 installation directory if it is not in the PATH already", default="")
cmdLineArgParser.add_argument("-multiBCFasta", dest="multiBCFastaFile", help="A multi-line fasta file defining multiplexing tag sequences. Required if there are multiplexing tags within the amplicon sequence as defined by the templateSeq file.")
cmdLineArgParser.add_argument("-pairedEnd", dest="pairedEnd", action="store_true",  help="Use if sequencing data is paired end")
cmdLineArgParser.add_argument("-RscriptExecPATH", dest="RscriptPATH", help="Rscript executable path if it is not in the PATH already. Must have Dada2 software package installed.", default="Rscript")
cmdLineArgParser.add_argument("-useUMI", dest="UMI", action="store_true",  help="Use flag if you want to remove PCR duplicate reads using UMI data")
cmdLineArgParser.add_argument("-numThreads", dest="numThreads", default=1,  help="Number of threads to be used for computation.")
cmdLineArgParser.add_argument("-skipSplitFastq", dest="skipSplitFastq", action="store_true",  help="Use flag if you want to skip the splitting of the raw fastq files (i.e. if you have already done this and do not want to redo it).")
cmdLineArgParser.add_argument("-remapBarcodes", dest="remapBarcodes", default=False,  help="Set to True if you want to remap barcodes even if the files already exist")
cmdLineArgParser.add_argument("-useBowtieMapping", dest="useBowtieMapping", action="store_true",  help="Set flag if you want to map with bowtie instead of blast")

#cmdLineArgParser.add_argument("-rebarcoding", dest="rebarcodingFile", help="File defining timepoints within each experiment when new barcodes were introduced. This feature is not currently implemented")

args = cmdLineArgParser.parse_args()

###########################################################################
## Utility functions
###########################################################################

##calculate the number of lines in the file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## Print to std error
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

#return the blast hit with lowest e value. requries input in outfmt6 format, output is an array with the columns separated as strings. 
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

	#extract coordinates of blast match dealing with reverse complemented sequences if necessary
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

# extract constant regions in template array, creates fasta file and blast database for each
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
			
#this gets regions from many reads at once
def extractRegionsFromFastqSingleBlast(readSeqRecordList, prefixName):
	#first, make a fasta file from all reads we are processing and blast against database of all index and constant regions
	readSeqFileName = "."+prefixName+"_readSeq.fasta"
	with open(readSeqFileName,"w") as outfasta:
		for readID in range(0,len(readSeqRecordList)):
			for i in range(0,len(templateSeqArray)):
				outfasta.write(">read_"+str(readID)+"_"+str(i)+"\n")
				outfasta.write(str(readSeqRecordList[readID][i].seq)+"\n")
	
	blastCommand = [args.blastPATH+"blastn","-query",readSeqFileName,"-db",allConstRegionsFileName]
	blastCommand.extend(constantRegionsBlastParams)
	blastResult = subprocess.check_output(blastCommand).decode('ascii').rstrip().split("\n") 
	finalReturnVal = []
	notConstRegionRegex = re.compile("^((?!const_region).)*$")
	oneTime = 0
	twoTime = 0
	threeTime = 0
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
		
		for readNum in range(0,len(templateSeqArray)): #for each read in the read pair we are supposed to be processing
			# first make an empty array for coordinates of features and get the blast result lines for this read
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
			
			#timeB = time.time()
			prevSegmentStartingCoord = -1
			reversedRead = False
			firstIndex = -1
			for i in range(0,len(templateSeqArray[readNum])): #for each feature in the template
				if(templateSeqArray[readNum][i]=="D" or (templateSeqArray[readNum][i]!="U" and templateSeqArray[readNum][i]!="X")): #if this is an indexing barcode location
					myLocalHits = []
					topBlastResult = ""
					if(templateSeqArray[readNum][i]=="D"):
						myLocalHits = list(filter(notConstRegionRegex.match, readBlastResult))
						topBlastResult = getBestBlastMatch(myLocalHits) #get the best hitting indexing region
					else:
						if constantRegionNames[readNum][i] not in blastLocationDict:
							continue
						myLocalHits = blastLocationDict[constantRegionNames[readNum][i]] #get blast hits for this specific constant region
						topBlastResult = getBestBlastMatch(myLocalHits)
					
					if(topBlastResult == [] or topBlastResult == "multiple"):
						continue
					
					if(firstIndex == -1):
						firstIndex = i
					#figure out where in the read this index region is hitting, set coordinates and track what the index region name was
					matchCoords = getSubjectMatchCoordinates(topBlastResult)
					if(prevSegmentStartingCoord > matchCoords[0] and not reversedRead):
						reversedRead = True
					
					if(prevSegmentStartingCoord == -1): #if this is the first segment we have found a hit for, figure out the right index where the segment is
						prevSegmentStartingCoord = matchCoords[0]
					if (not reversedRead): #f the read is not reversed, add in the coordinates properly
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
			if(not reversedRead and startingCoordinates[0] == None): #put the 0 in the right place depending on if the read is reversed or not relative to the template
				startingCoordinates[0] = 0
			if(reversedRead and startingCoordinates[len(startingCoordinates)-1] == None):
				startingCoordinates[len(startingCoordinates)-1] = 0
				
			if(reversedRead and startingCoordinates[0] == None): #set the end of the read to be readLength in the right place depending on if the read is reversed or not relative to the template
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
					if(len(UMIseq)==0):
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

		returnVal = [identifiedBCSeqRecords,identifiedUMISequences,identifiedIndexBCs, fwdRecord, revRecord]
		finalReturnVal.append(returnVal)
	return(finalReturnVal)
	
		
	
###########################################################################	
## Functions for Parsing Input Files and Directories
###########################################################################

def parseSampleFile():
	global sampleArray
	global sampleToIndexMap
	global indexToSampleMap
	with open(args.sampleFile) as infile: # go through every line in the sample file
		for line in infile:
			if(line.strip() == ""):
				continue
			lineSplit = line.strip().split("\t") #separate values in each line of the sample file
			if(len(lineSplit)< 3): # if the line has not enough columns, error and quit
				
				eprint("Sample file line:\n"+line+"has incorrect number of columns!\n")
				sys.exit(0)
			sampleIndexArray = "_".join(lineSplit[1:(len(lineSplit)+1)])
			sampleIdentityArray = lineSplit[0]
			
			if sampleIndexArray in indexToSampleMap or sampleIdentityArray in sampleToIndexMap: #if this sample is indistinguishable from a previous sample, error and quit
				eprint("Sample file line:\n"+line+"has a duplicate sample identity or multiplexing index with a previous line!\n")
				sys.exit(0)
			
			#else catalog the sample
			
			sampleArray.append(sampleStruct(lineSplit[0],lineSplit[1],lineSplit[2:(len(lineSplit)+1)]))
			sampleToIndexMap[sampleIdentityArray] = sampleIndexArray
			indexToSampleMap[sampleIndexArray] = sampleIdentityArray

# THIS FUNCTION IS NOT USED AS REBARCODING IS NOT CURRENTLY SUPPORTED!!
#def parseRebarcodingFile():

def parseTemplateSeq():
	sequence = ""
	with open(args.templateSeqFile) as infile:
		sequence = infile.readline().strip()
	if(len(re.sub("[ACGTacgtXUDN]","",sequence))>0): # if the sequence has disallowed characters, quit
		eprint("Template sequence has illegal characters!\n")
		sys.exit(0)
	if(args.multiBCFastaFile=="" and sequence.count("D")>0): #if there is no multiplexing file but we find multiplexing loci in the template, quit
		eprint("Multiplexing indices exist but no multiBCFastaFile provided!\n")
		sys.exit(0)
	#if it is paired end, split template in half for separate annotation, assume no useful features overlap the middle of the template, may want to take a parameter that is read size and split that way so you get more information out of the read, but need to assume each feature is in a single read only
	templateArray = []
	if(args.pairedEnd):
		templateArray.append(sequence[0:int(args.readLength)])
		templateArray.append(sequence[int(len(sequence)-args.readLength):(len(sequence)+1)])
	else:
		templateArray.append(sequence)
	#collapse features into single character and add to global variable
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

# figure out which fastq files we are going to use, pair forward and reverse reads as necessary, assume they are in alphabetical order
def identifyUsedFastqFiles():
	usedFastqFileDict = {}
	for sample in sampleArray:
		patternString = re.compile(".*"+sample.FilePrefix+".*.fastq(.gz)?")
		readFiles = glob.glob(args.fastqDir+"*")
		readFiles2 = list(filter(patternString.match,readFiles))
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

#separate reads by inline and illumina indices, separate UMI and BC regions for each read as well. This is multiprocessed code. THIS IS THE SLOWEST PART OF THE CODE, TRYING TO OPTIMIZE!

def demultiplexFastqHelper(readList, fastqPair, indexCounter, UMISeqFileHandleMap, badFwdReadsHandle, badRevReadsHandle):
	
	extractedRegions = extractRegionsFromFastqSingleBlast(readList, fastqPair.MatchPrefix)
	#these are all lists (or dictionaries of lists) of output so that we don't have to do too many seqio write calls since they are slow.
	badFwdReadsList = []
	badRevReadsList = []
	finalBCRecordList = {}
	fwdRecordList = {}
	revRecordList = {}
	for readID in range(0,len(extractedRegions)): #for each read we are processing
		fwdRecord = extractedRegions[readID][3]
		revRecord = extractedRegions[readID][4]
		totalNCount = str(fwdRecord.seq).count('N')
		if args.pairedEnd:
		    totalNCount = totalNCount  + str(revRecord.seq).count('N')
		if( totalNCount > maxNInReads or (args.UMI and len(extractedRegions[readID][1])==0) or len(extractedRegions[readID][0]) == 0): #if there are too many Ns in the read, or there is no barcode or no UMI if UMI is expected, it is bad and dump it
			badFwdReadsList.append(fwdRecord)
			if args.pairedEnd:
				badRevReadsList.append(revRecord)
			continue
		
		#get the sample from the multiplexing barcodes
		myIndex = [fastqPair.MatchPrefix]
		myIndex.extend(extractedRegions[readID][2])
		myIndex = "_".join(myIndex) #this is a string that uniquely defines each sample that was multiplexed
		if myIndex not in indexCounter:
			indexCounter[myIndex] = 0
		indexCounter[myIndex] = indexCounter[myIndex] + 1
		if(myIndex not in indexToSampleMap): #if we couldn't find a sample associated with this multiplexing barcode combination, dump read into a dump file and continue on
			badFwdReadsList.append(fwdRecord)
			if args.pairedEnd:
			    badRevReadsList.append(revRecord)
			continue
		else: #the read matches a valid sample
			mySample = indexToSampleMap[myIndex]
			if mySample not in UMISeqFileHandleMap: #populate map of file handles if necessary
				UMISeqFileHandleMap[mySample] = open(args.outputDir+mySample+"_UMISeqs.tab","w")
			if mySample not in finalBCRecordList:
				finalBCRecordList[mySample]=[]
				fwdRecordList[mySample]=[]
				revRecordList[mySample]=[]
				
			
			#concat all BCs associated with this read to get final BC
			finalBC = None
			for bcRecord in extractedRegions[readID][0]:
				if(type(bcRecord).__name__ == "NoneType"):
					continue
				if(type(finalBC).__name__ == "NoneType"):
					finalBC = bcRecord
				else:
					finalBC = finalBC + bcRecord
			if(type(finalBC).__name__ == "NoneType" or len(finalBC.seq) == 0):
				badFwdReadsList.append(fwdRecord)
				if args.pairedEnd:
				    badRevReadsList.append(revRecord)
				continue
			UMISeqFileHandleMap[mySample].write("\t".join(str(x) for x in extractedRegions[readID][1])+"\n") #write UMI seqs tab delimited
			finalBCRecordList[mySample].append(finalBC)
			fwdRecordList[mySample].append(fwdRecord)
			if args.pairedEnd:
			    revRecordList[mySample].append(revRecord)
			
	
	SeqIO.write(badFwdReadsList,badFwdReadsHandle,"fastq")
	if args.pairedEnd:
		SeqIO.write(badRevReadsList,badRevReadsHandle,"fastq")
	for mySample in finalBCRecordList.keys():
		barcodeFastqFileHandle = open(args.outputDir+mySample+"_barcode.fastq","a+")
		SampleSplitFastqFileHandleFWD = open(args.outputDir+mySample+"_R1.fastq","a+")
		SeqIO.write(finalBCRecordList[mySample],barcodeFastqFileHandle,"fastq") #write BC portion of read to a fastq (concat among all BC features)
		SeqIO.write(fwdRecordList[mySample],SampleSplitFastqFileHandleFWD,"fastq") #write raw fwd read to sample read file
		barcodeFastqFileHandle.close()
		SampleSplitFastqFileHandleFWD.close()
		
		if args.pairedEnd:
			SampleSplitFastqFileHandleREV = open(args.outputDir+mySample+"_R2.fastq","a+")
			SeqIO.write(revRecordList[mySample],SampleSplitFastqFileHandleREV,"fastq") #write raw rev read to sample read file
			SampleSplitFastqFileHandleREV.close()
	return indexCounter

def demultiplexFastq(fastqPair):
	# outputs: barcode file fastq, total split fastq, umi sequences, unidentified reads
	indexCounter = {}
	if(args.pairedEnd):
		print("Processing new FASTQ file[s]!\nFWD file:\t"+fastqPair.FwdFastq+"\nREV file:\t"+fastqPair.RevFastq+"\n"+fastqPair.MatchPrefix+"\n")
		#output file handles for unidentified reads
		badFwdReadsHandle = open(args.outputDir+fastqPair.MatchPrefix+"_unmappedReads_R1.fastq","w")
		badRevReadsHandle = open(args.outputDir+fastqPair.MatchPrefix+"_unmappedReads_R2.fastq","w")
		
		#input handles of fastq files, doesn't matter if they are gz compressed or not
		fwdFastqHandle = open(fastqPair.FwdFastq)
		revFastqHandle = open(fastqPair.RevFastq)
		if(fastqPair.FwdFastq.endswith("gz")):
			fwdFastqHandle = gzip.open(fastqPair.FwdFastq, "rt")
		if(fastqPair.RevFastq.endswith("gz")):
			revFastqHandle = gzip.open(fastqPair.RevFastq, "rt")
		
		#file handle maps
		UMISeqFileHandleMap = {}
		readCounter = 0
		readList = []
		
		for fwdRec, revRec in zip(SeqIO.parse(fwdFastqHandle,"fastq"), SeqIO.parse(revFastqHandle,"fastq")): #for each read, WE MAY WANT TO REPLACE THIS WITH GENERALFASTQITERATOR FOR PERFORMANCE! Zip should have OK performance, so that shouldn't be an issue
			readCounter = readCounter + 1
			readList.append([fwdRec,revRec])#.reverse_complement(name=True,description=True,id=True,annotations=True)]) #reverse complement the reverse read to match with the template sequence
			if(readCounter % fileBufferSize != 0): #if we are not at file buffer size, do not run parser since file io is super expensive and we want to minimize it
				continue
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, UMISeqFileHandleMap, badFwdReadsHandle, badRevReadsHandle) #process buffered reads using helper method
			readList = []
				
		
		#If we are here, the file has finished parsing, now need to empty the buffer for the last time, redo what we just did. This should probably get refactored, but too lazy right now
		if(len(readList)>0):
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, UMISeqFileHandleMap, badFwdReadsHandle, badRevReadsHandle)
		#close file handles!
		badFwdReadsHandle.close()
		badRevReadsHandle.close()
		fwdFastqHandle.close()
		revFastqHandle.close()
		for indexString in UMISeqFileHandleMap.keys():
			UMISeqFileHandleMap[indexString].close()

	else: #same as above, but for single read sequencing. THIS HAS NOT BEEN TESTED!
		#output file handles for unidentified reads
		print("Processing new FASTQ file[s]!\n"+fastqPair.FwdFastq+"\n"+fastqPair.MatchPrefix+"\n")
		badFwdReadsHandle = open(args.outputDir+fastqPair.MatchPrefix+"_unmappedReads_R1.fastq","w")
		
		#input handles of fastq files, doesn't matter if they are compressed or not
		fwdFastqHandle = open(fastqPair.FwdFastq)
		if(fastqPair.FwdFastq.endswith("gz")):
			fwdFastqHandle = gzip.open(fastqPair.FwdFastq)
		
		#file handle maps
		UMISeqFileHandleMap = {}
		readCounter = 0
		readList = []
		
		for fwdRec in SeqIO.parse(fwdFastqHandle,"fastq"): 
			readCounter = readCounter + 1
			readList.append([fwdRec])#.reverse_complement(name=True,description=True,id=True,annotations=True)]) #reverse complement the reverse read to match with the template sequence
			if(readCounter % fileBufferSize != 0): #if we are not at file buffer size, do not run parser since file io is super expensive and we want to minimize it
				continue
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, UMISeqFileHandleMap, badFwdReadsHandle, None) #process buffered reads using helper method
			readList = []
				
		
		#If we are here, the file has finished parsing, now need to empty the buffer for the last time, redo what we just did. This should probably get refactored, but too lazy right now
		if(len(readList)>0):
			indexCounter = demultiplexFastqHelper(readList, fastqPair, indexCounter, UMISeqFileHandleMap, badFwdReadsHandle, None)
		#close file handles!
		badFwdReadsHandle.close()
		fwdFastqHandle.close()
		for indexString in UMISeqFileHandleMap.keys():
			UMISeqFileHandleMap[indexString].close()
	
	
#do clustering across all samples
def clusterBarcodes():
	#concat all barcode fastq files by experiment into a single file for clustering, make a separate file using a subsample of the reads from each sample for making the error model

	allFiles = glob.glob(args.outputDir+"*_barcode.fastq")
	with open(args.outputDir+"allSamplesConcat.fastq","w") as outfile, open(args.outputDir+"allSamplesConcatForErrors.fastq","w") as outfile2:
		for fname in allFiles:
			with open(fname) as infile:
				lineCounter = 0
				for line in infile:
					if lineCounter < readsPerSampleForErrors*4:
						outfile2.write(line)
					lineCounter = lineCounter + 1
					outfile.write(line)
	#do the clustering with dada2 via r script, generates a fasta file of barcodes in output directory
	callFunction = [args.RscriptPATH,os.path.dirname(sys.argv[0]) + "clusterWithDada2.R",args.outputDir+"allSamplesConcat.fastq",args.outputDir+"allSamplesConcatForErrors.fastq", args.outputDir]
	subprocess.call(callFunction)
	args.barcodeListFile = args.outputDir+"clusteredBCs.fasta"
	
#map barcodes, multiprocessed code

def mapBarcodesHelper(indexString, readNames, UMIStrings, BCUMIMap, BCCountList, outputIDFileHandle):
	blastCommand = [args.blastPATH+"blastn","-query",indexString+"readBC.fasta","-db",args.barcodeListFile]
	blastCommand.extend(BCRegionsBlastParams)
	blastResult = subprocess.check_output(blastCommand).decode('ascii').rstrip().split("\n")
	blastResultDict = {}
	for line in blastResult: #store all blast hits in a dictionary for quick lookup
		mySplit = line.split("\t")
		curSample = mySplit[0]
		if curSample not in blastResultDict:
			blastResultDict[curSample]=[]
		blastResultDict[curSample].append(line)
	
	for j in range(0,len(readNames)):
		readName = readNames[j]
		umi = UMIStrings[j]
		if readName not in blastResultDict:
			outputIDFileHandle.write("NA\n")
			continue
		readBlastResults = blastResultDict[readName]
		topBlastResult = getBestBlastMatch(readBlastResults)
		if(topBlastResult == [] or topBlastResult == "multiple"):
			outputIDFileHandle.write("NA\n")
		else:
			outputIDFileHandle.write(topBlastResult[0]+"\n")
			mykey = topBlastResult[0]+"\t"+umi
			if mykey not in BCUMIMap or not args.UMI:
				BCUMIMap[mykey]=1
				BCCountList[int(topBlastResult[0])-1] = BCCountList[int(topBlastResult[0])-1] + 1
	return [BCUMIMap, BCCountList]
	
def mapBarcodes(indexString):
	#get file handles, map read maps
	if (os.path.isfile(args.outputDir+indexString+"_barcode.fastq") and (not os.path.isfile(args.outputDir+indexString+"_readBarcodeID.txt") or args.remapBarcodes)):
		bcFastqFile = args.outputDir+indexString+"_barcode.fastq"
		bcFastqFileHandle = open(bcFastqFile,"r")
		outputIDFileHandle = 	open(args.outputDir+indexString+"_readBarcodeID.txt","w")
		UMIFileHandle = open(args.outputDir+indexString+"_UMISeqs.tab","r")
		BCUMIMap = {}
		BCCountList = [0]*int(file_len(args.barcodeListFile)/2)
		readFastaFileHandle=open(indexString+"readBC.fasta","w") #this file name is unique for each process, since each process called with a unique index string
		readCounter = 0
		readNames = []
		UMIStrings = []
		
		#for each read bc / umi pair
		for record, UMIstring in zip(SeqIO.parse(bcFastqFileHandle,"fastq"), UMIFileHandle):
			readCounter = readCounter + 1
			SeqIO.write(record,readFastaFileHandle,"fasta") #dump bc seqs to a fasta file
			readNames.append(record.id)
			UMIStrings.append(UMIstring)
			if(readCounter % fileBufferSize == 0): #if the buffer is full and we need to analayze the reads
				#close the bc fasta file, run blast
				readFastaFileHandle.close()
				[resOne,resTwo] = mapBarcodesHelper(indexString, readNames, UMIStrings, BCUMIMap, BCCountList, outputIDFileHandle)
				BCUMIMap = resOne
				BCCountList = resTwo
				readFastaFileHandle=open(indexString+"readBC.fasta","w")
				readNames = []
				UMIStrings = []
		# the file is parsed, clear the buffer doing the same thing as before
		readFastaFileHandle.close()
		[resOne,resTwo] = mapBarcodesHelper(indexString, readNames, UMIStrings, BCUMIMap, BCCountList, outputIDFileHandle)
		BCUMIMap = resOne
		BCCountList = resTwo
		outputIDFileHandle.close()
		bcFastqFileHandle.close()
		UMIFileHandle.close()
		
		with open(args.outputDir+indexString+"_barcodeCalls.tab","w") as outFileHandle:
			for countVal in BCCountList:
				outFileHandle.write(str(countVal)+"\n")
		
def mapBarcodesWithBowtie2(indexString):
	#get file handles, map read maps
	if (os.path.isfile(args.outputDir+indexString+"_barcode.fastq") and (not os.path.isfile(args.outputDir+indexString+"_readBarcodeID_bowtie2.txt") or args.remapBarcodes)):
		bcFastqFile = args.outputDir+indexString+"_barcode.fastq"
		bcSamFile = args.outputDir+indexString+"_barcode.sam"
		bcIDFile = args.outputDir+indexString+"_readBarcodeID_bowtie2.txt"
		
		subprocess.call([args.bowtie2PATH+"bowtie2","-L 10","-q","-x "+args.barcodeListFile,"-U"+bcFastqFile,"-S"+bcSamFile])
		with open(bcIDFile,"w") as outfile:
			subprocess.call("grep -v '^@' "+bcSamFile+" | cut -f 3",stdout=outfile, shell=True)
		
		BCUMIMap = {}
		BCCountList = [0]*int(file_len(args.barcodeListFile)/2)
		
		bcIDFileHandle = open(bcIDFile,"r")
		UMIFileHandle = open(args.outputDir+indexString+"_UMISeqs.tab","r")
		
		#for each read bc / umi pair
		for bcID, UMIstring in zip(bcIDFileHandle, UMIFileHandle):
			bcID = bcID.strip()
			if(bcID == "*" or bcID == ""):
				continue
			else:
				mykey = bcID+"\t"+UMIstring
				if mykey not in BCUMIMap or not args.UMI:
					BCUMIMap[mykey]=1
					BCCountList[int(bcID)-1] = BCCountList[int(bcID)-1] + 1
		# the file is parsed, clear the buffer doing the same thing as before
		bcIDFileHandle.close()
		UMIFileHandle.close()
		
		with open(args.outputDir+indexString+"_barcodeCalls_bowtie2.tab","w") as outFileHandle:
			for countVal in BCCountList:
				outFileHandle.write(str(countVal)+"\n")
		

def generateFinalTables():
	#print barcode counts as giant tab delimited table, with 1st column as barcode ID number and header file being the sample each column comes from
	patternString = args.outputDir+"*_barcodeCalls*.tab"
	filenames = glob.glob(patternString)
	handles = [open(filename, 'r') for filename in filenames]    
	readers = [csv.reader(f, delimiter=',') for f in handles]
	filenames2 = ["BCID"]
	filenames2.extend(filenames)
	with  open(args.outputDir+"allBarcodeCalls.tab", 'w') as h:
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

def init(l):
	global lock
	lock = l

parseSampleFile()
parseTemplateSeq()
identifyUsedFastqFiles()
createConstRegionFasta()
l = Lock()
if(not args.skipSplitFastq):
	with Pool(processes = int(args.numThreads), initializer = init, initargs = (l,)) as pool:
		pool.map(demultiplexFastq, usedFastqFiles)
if args.barcodeListFile==None:
	clusterBarcodes()

#map barcodes
if(not args.useBowtieMapping):
	subprocess.call(["makeblastdb","-in",args.barcodeListFile,"-dbtype","nucl"])
else:
	subprocess.call([args.bowtie2PATH+"bowtie2-build",args.barcodeListFile,args.barcodeListFile])
with Pool(processes = int(args.numThreads)) as pool:
	if(not args.useBowtieMapping):
		pool.map(mapBarcodes, sampleToIndexMap.keys())
	else:
		pool.map(mapBarcodesWithBowtie2, sampleToIndexMap.keys())

generateFinalTables()
