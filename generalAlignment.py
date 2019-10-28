"""
generalAlignment.py
v1: 9/19/2019 Christopher Rohde

***** THE ONLY THING THAT SHOULD BE MODIFIED IN THIS FILE IS THE INPUTDIR STATEMENT! *****

This routine is used to generally align a sequence to a reference.
This program assumes that the sequences are in a folder
	./Input/<INPUTDIR>/xxxx_seq (there can be multiple _seq directories, and you can kept the name that genewiz gives them)
and that there exists the files
	./Input/Sequences.txt and (optional) ./Input/References.fasta

References.txt is a fasta file containing the sequences to be aligned (CDS only, without the UTRs or tails)
There is also a master references file,
	./MasterReferences.fasta
which contains commonly used sequences. Note that if the same named entry exists in both References and MasterReferences,
the sequence in References will be used for alignment.

Sequences.txt contains whitespace separated values of the format
<name>	<reference sequence name>	<optional:tail type (TT_BsaI,TT_BbsI, None)> <optional2:UTRspecificationstring>
example:
mySeq_GeneA	GeneA #Defaults to TT_BsaI tails
mySeq_GeneB	GeneB TT_BsaI
mySeq_GeneC	GeneC TT_BsaI 5Hbb|3wpreA

The default tail configuration is TT_BsaI, which is (A9G)15-AAAAAA-GAGACC. Tail configurations can be added by modifying
the file ./TailReferences.fasta

The default UTR configuration is 5'-Hbb (specified as 5Hbb), 3'-Hbb (specified as 3Hbb). UTR configurations can be added
by modifying the file ./UTRReferences.fasta

Note that both UTRs must be specified, and separated by a pipe character (|)


***** THE ONLY THING THAT SHOULD BE MODIFIED IN THIS FILE IS THE INPUTDIR STATEMENT! *****
"""
INPUTDIR = '20191028 tgfB variants midis and YFP-TT1'

import os, shutil, glob, pathlib, datetime
from Bio import SeqIO
#The block below is required to import all of the nvLib elements
import pathlib,sys,os
nvLibPath = pathlib.Path.cwd().parent.joinpath('nvLib')
if not nvLibPath.exists():
	print("ERROR: Can't locate nvLib. Expected to find it at\n\t"+str(nvLibPath))
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, str(nvLibPath))
#imports below are nvLib elements
import novesliceElementTypeDetermination,seqFileProcessing,GAfileProcessing,inputFileProcessing

def main():
	
	ROOTINPUTDIR = 'Input'
	OUTPUTDIR = 'Aligned'
	REFERENCEDIR = 'References'
	TEMPDIR = 'TempFiles'
	MASTERREFERENCEFNAME = 'SequenceReferences.fasta'
	UTRREFERENCEFNAME = 'UTRreferences.fasta'
	TAILREFERENCEFNAME = 'TailReferences.fasta'
	pathForAlignmentFiles = pathlib.Path.cwd() / ROOTINPUTDIR / INPUTDIR
	pathForOutputFiles = pathlib.Path.cwd() / ROOTINPUTDIR / INPUTDIR / OUTPUTDIR
	pathForTempFiles = 	pathlib.Path.cwd() / TEMPDIR
	#Create the temp dir if it doesn't exist
	pathForTempFiles.mkdir(exist_ok = True)

	#Check that the input directories exist (Input/INPUTDIR/XXXX_seq)
	sequenceFileDirs = seqFileProcessing.checkForSequenceDirectories(pathForAlignmentFiles)
	if not sequenceFileDirs: return

	#If there are multiple sequence directories, combine all of the files into a single directory (the first one in the list)
	seqFileProcessing.combineSequenceFilesInOneDirectory(sequenceFileDirs)
	#Now we can use the first of the sequence directories since it will contain all of the files
	mainInputDir = sequenceFileDirs[0]

	#Master sequence reference file location (CDS, tail, and UTR sequences are in this directory)
	referenceDir = pathlib.Path.cwd() / REFERENCEDIR

	#Read in a file (Sequences.txt) containing the name, reference sequence name,
	#and (optionally) the tail and utr configuration
	sequencesNamesAndTypes = GAfileProcessing.parseGeneralAlignmentSequencesFile(pathForAlignmentFiles / 'Sequences.txt')
	if not sequencesNamesAndTypes: return
	#Get the local reference sequences, if any
	localReferenceFile = pathForAlignmentFiles / 'References.txt'
	localReferences = getReferencesDict(localReferenceFile) if localReferenceFile.exists() else {}

	#Make the output directory if it does not exist
	if not os.path.exists(pathForOutputFiles):
		print('Creating output directory\n\t'+str(pathForOutputFiles))
		os.mkdir(pathForOutputFiles)
		print('SUCCESS')
	else:
		print('Output will be stored in\n\t' + str(pathForOutputFiles))

	#Get all of the master sequence references, the possible UTRs, and possible tails
	masterReferences= getReferencesDict(referenceDir / MASTERREFERENCEFNAME)
	utrReferences 	= getReferencesDict(referenceDir / UTRREFERENCEFNAME)
	tailReferences 	= getReferencesDict(referenceDir / TAILREFERENCEFNAME)
	#Determine which of the local sequence references are present in the master reference, and raise a warning.
	sequenceReferences = GAfileProcessing.compareAndMergeReferenceDicts(masterReferences,localReferences)

	#TO_DO:
	#determine file groupings (which seq file elements we have)
	seqFileSet = seqFileProcessing.getRootNamesForSeqFiles(mainInputDir)
	#Determine which ones are missing from the Sequences.txt file and raise a warning (those items won't be aligned)
	seqFileSetWithInfoToAlign = GAfileProcessing.getFilesWithSequencesFileEntries(seqFileSet,sequencesNamesAndTypes)

	#Determine which ones have a reference not contained in the sequencesReferences dict
	seqFileSetWithReferenceSequences = fillInSequenceReferences(seqFileSetWithInfoToAlign,
										sequenceReferences, utrReferences, tailReferences)

	outputText = ''
	for pair in seqFileSetWithReferenceSequences:
		score = alignFilesToSequence(pair[0], mainInputDir, pathForTempFiles, pair[1],pathForOutputFiles)
		if score[0] == len(pair[1]):
			goodText = '\tPROBABLY GOOD'
		else:
			goodText = ''

		outStr = ('{:.<30}Longest match: {:5} / {:<5}\tOverall match: {:5} / {:<5}\t{}'.format(
			pair[0],score[0],len(pair[1]),score[1],len(pair[1]),goodText))
		print(outStr)
		outputText += str(outStr) + '\n'

	with open(pathForAlignmentFiles / 'AlignmentResults.txt','w') as f:
		scriptName = pathlib.Path(__file__).name
		f.write(scriptName + " alignment results performed at " + str(datetime.datetime.now()) + '\n')
		f.write(outputText)

	#Delete the temporary files
	for fi in pathForTempFiles.glob('*'):
		os.remove(fi)

def alignFilesToSequence(fileGroup, fileDir, tempDir, sequence, outputDir):
	fwdFiles,revFiles = seqFileProcessing.getSeqFileOrientations(fileDir,fileGroup)
	fragments = []
	for file in fwdFiles:
		fileName = fileGroup + '-' + file + '.seq'
		fragments.append(generateAlignmentFragment(fileDir / fileName , sequence, tempDir, reverseSequence = False))
	for file in revFiles:
		fileName = fileGroup + '-' + file + '.seq'
		fragments.append(generateAlignmentFragment(fileDir / fileName , sequence, tempDir, reverseSequence = True))
	#Pad the fragments to account for the fact that the M13F reference starts earlier than the reference sequence
	[paddedFragments,paddingLen] = padFragments(fragments)
	paddedReference = ' '*paddingLen + sequence
	
	outStr = processAlignmentFragmentsIntoSingleString(paddedReference,paddedFragments)
	#Generate alignment score string
	scoreStr = scoreFragments(paddedReference,paddedFragments)
	#calculate the score by finding the longest region between two 0s
	score = calculateLongerRunWithoutChar(scoreStr,'0')
	#Write some preamble text to the file
	#Write the alignment and save the file
	preambleText = ""

	allSequenceFiles = ", ".join(fwdFiles) + ", " + ", ".join(revFiles) 
	scriptName = pathlib.Path(__file__).name
	preambleText += scriptName + " alignment results performed at " + str(datetime.datetime.now()) + '\n'
	preambleText += 'Alignment of ' + fileGroup + ' using sequence files ' + allSequenceFiles + '\n'
	preambleText += 'Longest match: {} / {}\nOverall match: {} / {}'.format(score[0],len(sequence),score[1],len(sequence))
	preambleText += '\n\n'
	with open(outputDir.joinpath(fileGroup + "-aligned.txt"),'w') as f:
		f.write(preambleText)
		f.write(outStr)
		f.write('Score\t\t\t\t'+ scoreStr)
	return score

def padFragments(fragments):
	padding = [len(str(fragment[0].seq))-len(str(fragment[0].seq).lstrip('-')) for fragment in fragments]
	amountToPad = [max(padding) - el for el in padding]
	paddedFragments = []
	for i,fragmentGroup in enumerate(fragments):
		padStr = ' '*amountToPad[i]
		paddedFragments.append([padStr+fragmentGroup[0],
								padStr+fragmentGroup[1],
								padStr+fragmentGroup[2]])
	return([paddedFragments,max(padding)])

def scoreFragments(sequence,fragments):
	#Read all of the alignments as a list
	outStr = ''
	for i,referenceBase in enumerate(sequence):
		outVal = 0
		for fragmentGroup in fragments:
			if fragmentGroup[2][i] == '*':
				outVal += 1
		outStr += str(outVal)
	#Calculate longest fragment by locating all of the zeros in the score string
	return outStr

def calculateLongerRunWithoutChar(inStr, inCh):
	#note we add a zero to the end of the string to allow counting up to the end of the sequence
	locations = [i for i,ch in enumerate(inStr+'0') if ch==inCh]
	distanceToNext = [locations[i+1]-locations[i]-1 for i,pos in enumerate(locations[:-1])]
	return [max(distanceToNext),sum(distanceToNext)]

def processAlignmentFragmentsIntoSingleString(sequence, fragments):
	outStr = ''
	#Print each line in a file.  We need to pad the original reference sequences to account for
	#the M13F primer starting before the start of the reference
	outStr += "Master Reference\t"
	outStr += str(sequence) + '\n\n'
	for i,fragment in enumerate(fragments):
		outStr += 'Ref.' + str(i) + '\t\t\t\t'
		outStr += str(fragment[0].seq) + '\n'
		outStr += 'Seq.' + str(i) + '\t\t\t\t'
		outStr += str(fragment[1].seq) + '\n'
		outStr += 'Match.0' + str(i) + '\t\t\t'
		outStr += fragment[2] + '\n'
	return outStr

def generateAlignmentFragment(inFile,referenceSequence,tempDir,reverseSequence=False):
	with open(inFile,'r') as fi:
		mySeq = list(SeqIO.parse(fi, "fasta"))[0].seq
	
	outText = '>REFERENCE\n' + referenceSequence + '\n'
	if not reverseSequence:
		outText += '>' + inFile.name + '\n' + mySeq
	else:
		outText += '>' + inFile.name + '\n' + mySeq.reverse_complement()

	#Generate the file to be used as input to clustalw
	alignmentInputFileN = tempDir / inFile.stem
	with open(alignmentInputFileN,'w') as fo:
		fo.write(str(outText))
	#Get the command line for running clustalw and run it
	from Bio.Align.Applications import ClustalwCommandline
	#We greatly increased the gap opening penalty which seems to yield better local alignment
	clustalw_cline = ClustalwCommandline('clustalw2',infile=str(alignmentInputFileN),gapopen=50)
	x=clustalw_cline()
	return parseAlignmentFragment(inFile)

def parseAlignmentFragment(inFile):
	from Bio import AlignIO
	align = AlignIO.read('./TempFiles/'+ inFile.stem + '.aln','clustal')
	return [align[0],align[1],align._star_info]

def getReferencesDict(FILE):
	if not FILE.exists():
		print('ERROR: Cannot locate a references file.  Expected to find "' + 
			str(tailReferencesFi.relative_to(pathlib.Path.cwd())) + '"')
	referencesLi = inputFileProcessing.getAllFastaSequencesInFile(FILE)
	return inputFileProcessing.convertFastaListToDict(referencesLi)

def fillInSequenceReferences(seqFileSetWithInfoToAlign,	sequenceReferences, utrReferences, tailReferences):
	filePrefixAndSequenceToAlign = []
	for entry in seqFileSetWithInfoToAlign:
		[seq,tail,utr5,utr3] = [entry[1].upper(), entry[2].upper(), entry[3].upper(), entry[4].upper()]
		if seq not in sequenceReferences.keys():
			print('Error in entry ' + entry[0] + '. Sequence ' + entry[1] + ' not found in reference file')
			continue

		if tail not in tailReferences.keys() and tail!='NONE':
			print('Error in entry ' + entry[0] + '. Tail ' + entry[2] + ' not found in reference file')
			continue
		if utr5 not in utrReferences.keys() and utr5!='NONE':
			print('Error in entry ' + entry[0] + '. UTR ' + entry[3] + ' not found in reference file')
			continue
		if utr3 not in utrReferences.keys() and utr3!='NONE':
			print('Error in entry ' + entry[0] + '. UTR ' + entry[4] + ' not found in reference file')
			continue
		tailSeq = tailReferences[tail] if tail!='NONE' else ""
		utr5Seq = utrReferences[utr5] if utr5!='NONE' else ""
		utr3Seq = utrReferences[utr3] if utr3!='NONE' else ""

		sequence = utr5Seq + sequenceReferences[seq] + utr3Seq + tailSeq

		filePrefixAndSequenceToAlign.append([entry[0],sequence])
	return filePrefixAndSequenceToAlign

"""



def checkForSequenceDirectories(pathForAlignmentFiles):
	if not os.path.exists(pathForAlignmentFiles):
		print('Error: Input directory "' + pathForAlignmentFiles + '" does not exist.')
		return False
	sequenceFileDirs = [os.path.join(pathForAlignmentFiles,dir) for dir in os.listdir(pathForAlignmentFiles) if dir.endswith('_seq')]
	if not sequenceFileDirs:
		print('Error: No sequence file folders found in input directory "' + pathForAlignmentFiles + '"')
		return False
	return sequenceFileDirs

#Give a directory containing .seq files, determining the prefixes that contain either of the following groups of files
#	prefix-F1.seq, prefix-F2.seq, and prefix-R1.seq:				We can do a partial alignment.
#	all of the above files AND prefix-M13F.seq, prefix-M13R.seq:	We can do a full alignment
#Return a list of lists, each entry containing [prefix, TRUE/FALSE for Full/Partial alignment].
#Files that belong to a group not matching either of the above conditions are not returned.
def determineFileGroupings(sequenceFileDir):
	#Now that all of the sequence files are in a single directory, get the files
	files = glob.glob(os.path.join(sequenceFileDir,'*.seq'))
	#Strip out all of the path information
	fileNames = [os.path.basename(os.path.normpath(file)) for file in files]
	#Strip the primer identifier and extension from the file name (e.g. -F1.seq) to create a set of the potential alignments
	filePrefixesSet = set([file.rsplit('-')[0] for file in fileNames])
	#Convert back into a sorted list
	fileGroupings = list(filePrefixesSet)
	fileGroupings.sort()
	#For each list item, determine if we can do a full alignment (M13F, F1, F2, R1, M13R) + optional R2, 
	#or a partial alignment (F1,F2,R1), or we can't align
	returnedGroupings = []
	for grouping in fileGroupings:

		if all([grouping+'-F1.seq' in fileNames, grouping+'-F2.seq' in fileNames, grouping+'-R1.seq' in fileNames]):
			if all([grouping+'-M13F.seq' in fileNames , grouping+'-M13R.seq' in fileNames, grouping+'-R2.seq' in fileNames]):
				returnedGroupings.append([grouping, 2]) #Full alignment with M13F, M13R and R2
			if all([grouping+'-M13F.seq' in fileNames , grouping+'-M13R.seq' in fileNames]):
				returnedGroupings.append([grouping, 1]) #Full alignment with M13F and M13R only
			else:
				returnedGroupings.append([grouping, 0]) #Partial alignment
		else:
			print('{:.<30}Not enough files to do an alignment of file group'.format(grouping))
	return(returnedGroupings)

#Given the prefix for a file in the sequence file directory, read the M13R file and use that to determine the catalytic domain, if possible
def determineCatalyticDomains(sequenceFileDir, groupings):
	catDoms = []
	for entry in groupings:
		if not entry[1]:
			catDoms.append('FokI')
		else:
			fi = open(os.path.join(sequenceFileDir,entry[0] +'-M13R.seq'))
			mySeq = list(SeqIO.parse(fi,"fasta"))[0].seq  
			catDoms.append(novesliceElementTypeDetermination.determineCatalyticDomain(mySeq))
	return catDoms

def determineNovesliceConfigurations(sequenceFileDir, groupings):
	configs = []
	for entry in groupings:
		fi = open(os.path.join(sequenceFileDir,entry[0] +'-F1.seq'))
		mySeq = list(SeqIO.parse(fi,"fasta"))[0].seq
		configs.append(novesliceElementTypeDetermination.determineDnaBindingType(mySeq))
	return configs

def mergeAllInformationAboutGroupings(alignmentNamesAndTypes,configurations,catDoms):
	mergedList = []
	for i,entry in enumerate(alignmentNamesAndTypes):
		mergedList.append([entry[0],entry[1],configurations[i],catDoms[i]])
	return mergedList

def readInputFile(fileDir):
	inFile = os.path.join(fileDir,'Sequences.txt')
	if not os.path.exists(inFile):
		print("Error: Can't find file 'Sequences.txt' in directory '"+fileDir+"'\n(ensure the name of the file matches this EXACTLY)")
		return None
	sequencesNamesAndTypes = []
	f = open(inFile)
	for line in f:
		lineItems = line.split()
		if len(lineItems) < 2:
			print('Invalid line ' + line)
		else:
			if len(lineItems) == 3: #i.e. line species a TAL or noveslice configuration
				config = parseConfiguration(lineItems[2])
			if len(lineItems) == 2: #i.e. line doesn't specify a configuration, so see if you can infer it from its name
				config = parseGEPname(lineItems[0])
			if config:
#				print('{:<30}Alignment configuration determined'.format(lineItems[0]))
				sequencesNamesAndTypes.append([lineItems[0], lineItems[1],config[0],config[1]])
			else:
				print('{:<30}Error: Could not determine configuration (NoveSlice or TALEN)'.format(lineItems[0]))
	return sequencesNamesAndTypes

def parseConfiguration(inStr):
	if inStr[0].upper()=='T':
		return [False, 0]	#We have a TALEN, so riboslice_v2=False and the value of riboslice_50percent doesn't matter
	elif inStr[0].upper()=='N' or inStr[0].upper() == 'R':
		if inStr[-1].upper()=='A':
			return [True,1] #We have a NoveSlice A configuration, so riboslice_v2=True and riboslice_50percent=1
		elif inStr[-1].upper()=='B':
			return [True,2] #We have a NoveSlice B configuration, so riboslice_v2=True and riboslice_50percent=2
	return None

"""
def parseGEPname(inStr):
	"""Explanation of the regex (?i)_(T|ns|rs)\d{1,2}[L|F|R]_?([A|B]?)
	(?i)		The regex is case insensitive
	_ 			The first character matched must be a single underscore
	T|ns|rs		Match group 1: The underscore must be followed by either a 'T', a 'ns' or 'rs' (case insensitive)
					We will assume the type of sequence (TALEN or Noveslice) based on this unless something else is
					specified in the alignment file
	\d{1,2}  	The T/ns/rs must be followed by one or two digits.
					This assumes we won't generate more than 99 items with the same name (or we'd need more digits)
	[L|F|R]		The digits much be followed by an orientation (L or F or R).
	_?			The orientation is followed by zero or one underscores
	([A|B]?)	Match group 2: the orientation is followed by zero or one A or B indicating the noveslice orientation 
	"""
	import re
	pattern = '(?i)_(T|ns|rs)\d{1,2}[L|F|R]_?([A|B]?)'
	x=re.search(pattern,inStr)
	if x:
		if x[1].upper() == 'T':
			return [False, 0]	#We have a TALEN, so riboslice_v2=False and the value of riboslice_50percent doesn't matter
		elif x[1].upper() == 'NS' or x[1].upper() == 'RS':
			if x[2].upper() == 'A':
				return [True,1] #We have a NoveSlice A configuration, so riboslice_v2=True and riboslice_50percent=1
			elif x[2].upper() == 'B':
				return [True,2] #We have a NoveSlice B configuration, so riboslice_v2=True and riboslice_50percent=2
	print("{:.<30}No configuration specified and we can't infer the configuration from its name".format(inStr))
	return None #Otherwise there wasn't enough information in the gene editing protein name to determine its type

def determineAlignableSequences(alignmentNamesAndTypes,sequenceNamesAndTypes):
	alignmentNames = [x[0] for x in alignmentNamesAndTypes]
	sequenceNames = [x[0] for x in sequenceNamesAndTypes]

	#First determine which sequences have .seq files, but are not in Sequences.txt
	lackingRefSeq = set(alignmentNames) - set(sequenceNames)
	for el in sorted(lackingRefSeq):
		#if the item name ends with an underscore and a number (e.g. AAVS_T2R_1), then check if there's an item the
		#sequence file that matches the form without the underscore and number (e.g. AAVS_T2R)
		truncatedName = el.rsplit('_',1)[0]
		if truncatedName in sequenceNames:
			# print("{:.<30}Aligning based on sequence {}".format(el,truncatedName))
			#Note the .copy() here which is important to prevent modification of the original list entry
			longerEntry = sequenceNamesAndTypes[[x[0] for x in sequenceNamesAndTypes].index(truncatedName)].copy()
			longerEntry[0]=el
			sequenceNamesAndTypes.append(longerEntry)
		else:
			print("{:.<30}Files exist for alignment however no reference sequence found in Sequences.txt".format(el))
	#Next, determine which sequences have references but no sequence files (in Sequences.txt but no .seq files)
	lackingSeqFiles = set(sequenceNames) - set(alignmentNames)
	#We'll do the opposite of the above, and check if there are any alignments with an underscore, and if so, the
	#Sequences.txt file entry will be assumed to be for that
	sequenceNamesTruncated = [x[0].rsplit('_',1)[0] for x in sequenceNamesAndTypes]
	for el in sorted(lackingSeqFiles):
		if el not in sequenceNamesTruncated:
			print("{:.<30}Reference sequence found in Sequences.txt however not enough .seq files found".format(el))
	#Return all of the items that have both, as well as the parameters that will be used to align them
	#Note that we should re-determine sequenceNames since we altered it above
	sequenceNames = [x[0] for x in sequenceNamesAndTypes]
	alignableSequenceNames =  list(set(alignmentNames).intersection(set(sequenceNames)))
	alignableSequences = []
	#XXXXXXXXXXXXXXXXXX TO DO XXXXXXXXXXXXXXXXXXXX
	#Use the determined catalytic domain and 
	for el in sorted(alignableSequenceNames):
		alignParms = alignmentNamesAndTypes[[x[0] for x in alignmentNamesAndTypes].index(el)]
		seqParms = sequenceNamesAndTypes[[x[0] for x in sequenceNamesAndTypes].index(el)]
		#Check if the configuration specified in the Sequences.txt file is the same as we guess from looking at the translation of -F1.seq
		#If not, raise a warning
		if alignParms[2] != seqParms[3]:
			print('Warning for item ' + seqParms[0] + '\n\tDifferent configuration specified in Sequences.txt than was estimated based on examining the sequence')
			print('\tSpecified: ' + numberToNovesliceConfiguration(seqParms[3]) + '. Estimated: ' + numberToNovesliceConfiguration(alignParms[2])+ '. Defaulting to ' + numberToNovesliceConfiguration(seqParms[3]))
		#if we couldn't determine the catalytic domain, the routine returns None.  Default the catalytic domain to FokI in that case.
		if not alignParms[3]:
			catDom = 'FokI'
		else:
			catDom = alignParms[3]
		alignableSequences.append([seqParms[0],seqParms[1],alignParms[1],seqParms[2],seqParms[3],catDom])
	return alignableSequences

def numberToNovesliceConfiguration(num):
	if num == 0:
		return 'TALEN'
	elif num == 1:
		return 'NoveSlice A'
	elif num == 2:
		return 'NoveSlice B'
	else:
		return 'Unknown'


def alignIndividualSequence(recognitionSequence,
							filePrefix,
							inputDir,
							outputDir,
							nucleaseDomain='FokI',
							riboslice_v2=False,
							riboslice_50percent=1,
							alignmentLevel=0
							):

	#generateFullSequence returns a list [full sequence (seq), RVD location (string)]
	fullSeq = generateFullSequence(	recognitionSequence, nucleaseType = nucleaseDomain,
									riboslice_v2=riboslice_v2, riboslice_50percent=riboslice_50percent,
									include_pIDT=alignmentLevel)
	#generateAlignmentFragments returns 1-kb sequence fragments of the reference sequence
	#that follow the sequencing primers
	alignmentFragments = generateAlignmentFragments(fullSeq[0],alignmentLevel=alignmentLevel)
	#generateAlignments uses clustalw to align the files returned by genewiz to the fragments
	#created by generateAlignmentFragments.  It saved the alignment results to the temp dir
	generateAlignments(filePrefix,alignmentFragments,inputDir,alignmentLevel=alignmentLevel)
	#parseAlignments reads the temp files returned by clustalw and returns a list of lists
	#that contains the various alignments [[seqA_ref, seqA_genewiz, stars],[...]]
	alignments = parseAlignments(filePrefix,alignmentLevel=alignmentLevel)
	#Generate some preamble text
	preambleText = ""
	preambleText += filePrefix + "\n"
	preambleText += "Recognition sequence: " + recognitionSequence + "\n" 
	preambleText += "Nuclease Domain: " + nucleaseDomain + "\n"
	if riboslice_v2: preambleText += "DNA binding type: riboslice v2 (extra GHGG linkers)\n"
	else: preambleText += "DNA binding type: Sanjana2011 TALENs\n"
	preambleText += "Number of elements aligned = "
	if alignmentLevel == 0:
		preambleText += '3'
	elif alignmentLevel == 1:
		preambleText += '5'
	elif alignmentLevel == 2:
		preambleText += '6'
	else:
		preambleText += 'unknown alignment level specified'
	preambleText += '\n'
	import datetime
	preambleText += str(datetime.datetime.now()) + "\n"
	preambleText += "======================================================================\n"

	return printAlignments(filePrefix,outputDir,fullSeq[0],fullSeq[1],alignments,alignmentLevel=alignmentLevel,
					preambleText=preambleText)

def generateFullSequence(recognitionSequence, nucleaseType, riboslice_v2 = False, riboslice_50percent=1,include_pIDT = 0):
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna

	#Variables occasionally modified
	baseToRVDs = {'A':'NI','T':'NG','C':'HD','G':'NN'}

	#Constants
	#pIDTsmart
	pIDT_fivePrime = Seq('CCCGTGTAAAACGACGGCCAGTTTATCTAGTCAGCTTGATTCTAGCTGATCGTGGACCGGAAGGTGAGCCAGTGAGTTGATTGCAGTCCAGTTACGCTGGAGTCTGAGGCTCGTCCTGAATGATATGCGACCGCCGGAGGGTTGCGTTTGAGACGGGCGACAGATCCAGTCGCGCTGCTCTCGTCGATCCAAGCTT',generic_dna)
	pIDT_threePrime = Seq('GAATTCGGTGCGAGCGGATCGAGCAGTGTCGATCACTACTGGACCGCGAGCTGTGCTGCGACCCGTGATCTTACGGCATTATACGTATGATCGGTCCACGATCAGCTAGATTATCTAGTCAGCTTGATGTCATAGCTGTTTCCTGAGGCTCA',generic_dna)
	#UTRs
	five_prime_UTR = Seq('TAATACGACTCACTATAGGGACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGCTAGCCACC',generic_dna)
	three_prime_UTR = Seq('ACCGGTGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC',generic_dna)
	#Protein termini (N_terminus includes NLS / other upstream stuff)
	N_terminus = Seq('ATGGACTATAAGGACCACGACGGAGACTACAAGGATCATGATATTGATTACAAAGACGATGACGATAAGATGGCCCCAAAGAAGAAGCGGAAGGTCGGTATCCACGGAGTCCCAGCAGCCGTAGATTTGAGAACTTTGGGATATTCACAGCAGCAGCAGGAAAAGATCAAGCCCAAAGTGAGGTCGACAGTCGCGCAGCATCACGAAGCGCTGGTGGGTCATGGGTTTACACATGCCCACATCGTAGCCTTGTCGCAGCACCCTGCAGCCCTTGGCACGGTCGCCGTCAAGTACCAGGACATGATTGCGGCGTTGCCGGAAGCCACACATGAGGCGATCGTCGGTGTGGGGAAACAGTGGAGCGGAGCCCGAGCGCTTGAGGCCCTGTTGACGGTCGCGGGAGAGCTGAGAGGGCCTCCCCTTCAGCTGGACACGGGCCAGTTGCTGAAGATCGCGAAGCGGGGAGGAGTCACGGCGGTCGAGGCGGTGCACGCGTGGCGCAATGCGCTCACGGGAGCACCCCTCAA',generic_dna)
	C_terminus = Seq('TCAATCGTGGCCCAGCTTTCGAGGCCGGACCCCGCGCTGGCCGCACTCACTAATGATCATCTTGTAGCGCTGGCCTGCCTCGGCGGACGACCCGCCTTGGATGCGGTGAAGAAGGGGCTCCCGCACGCGCCTGCATTGATTAAGCGGACCAACAGAAGGATTCCCGAGAGGACATCACATCGAGTGGCA',generic_dna)
	#Previously built a construct with a different C-terminus before the StsI site.  Don't use this anymore.
	#C_terminus_StsI = Seq('TCAATCGTGGCCCAGCTTTCGAGGCCGGACCCCGCGCTGGCCGCACTCACTAATGATCATCTTGTAGCATTAGCGTGCCTCGGCGGACGACCCGCCTTGGATGCGGTGAAGAAGGGGCTCCCGCACGCGCCTGCATTGATTAAGCGGACCAACAGAAGGATTCCCGAGAGGACATCACATCGAGTGGCA',generic_dna)
	#Nuclease domains
	Linker = Seq('GGTTCC',generic_dna)
	FokI = Seq('CAACTCGTGAAGAGTGAACTTGAGGAGAAAAAGTCGGAGCTGCGGCACAAATTGAAATACGTACCGCATGAATACATCGAACTTATCGAAATTGCTAGGAACTCGACTCAAGACAGAATCCTTGAGATGAAGGTAATGGAGTTCTTTATGAAGGTTTATGGATACCGAGGGAAGCATCTCGGTGGATCACGAAAACCCGACGGAGCAATCTATACGGTGGGGAGCCCGATTGATTACGGAGTGATCGTCGACACGAAAGCCTACAGCGGTGGGTACAATCTTCCCATCGGGCAGGCAGATGAGATGCAACGTTATGTCGAAGAAAATCAGACCAGGAACAAACACATCAATCCAAATGAGTGGTGGAAAGTGTATCCTTCATCAGTGACCGAGTTTAAGTTTTTGTTTGTCTCTGGGCATTTCAAAGGCAACTATAAGGCCCAGCTCACACGGTTGAATCACATTACGAACTGCAATGGTGCGGTTTTGTCCGTAGAGGAACTGCTCATTGGTGGAGAAATGATCAAAGCGGGAACTCTGACACTGGAAGAAGTCAGACGCAAGTTTAACAATGGCGAGATCAATTTCCGCTCATAA',generic_dna)
	FokI_HE = Seq('CAACTCGTGAAGAGTGAACTTGAGGAGAAAAAGTCGGAGCTGCGGCACAAATTGAAATACGTACCGCATGAATACATCGAACTTATCGAAATTGCTAGGAACCCGACTCAAGACAGAATCCTTGAGATGAAGGTAATGGAGTTCTTTATGAAGGTTTATGGATACCGAGGGGAGCATCTCGGTGGATCACGAAAACCCGACGGAGCAATCTATACGGTGGGGAGCCCGATTGATTACGGAGTGATCGTCGACACGAAAGCCTACAGCGGTGGGTACAATCTTCCCATCGGGCAGGCAGATGAGATGCAACGTTATGTCGAAGAAAATCAGACCAGGAACAAACACATCAATCCAAATGAGTGGTGGAAAGTGTATCCTTCATCAGTGACCGAGTTTAAGTTTTTGTTTGTCTCTGGGCATTTCAAAGGCAACTATAAGGCCCAGCTCACACGGTTGAATCACATTACGAACTGCAATGGTGCGGTTTTGTCCGTAGAGGAACTGCTCATTGGTGGAGAAATGATCAAAGCGGGAACTCTGACACTGGAAGAAGTCAGACGCAAGTTTAACAATGGCGAGATCAATTTCCGCTCATAA',generic_dna)

	StsI = Seq('TGTGGCTCAGTTGTCAAGGCCTGATCCTGCATTGGCCGCTCTGACCAACGATCACCTGGTAGCTTTGGCTTGTTTAGGCGGACGGCCAGCTCTTGATGCTGTTAAGAAGGGATTGCCTCATGCTCCCGCCCTGATTAAGAGGACCAACCGACGTATTCCAGAGAGAACCTCCCACAGGGTGGCTGGATCTGTGTTGGAGAAGAGCGACATCGAGAAGTTCAAGAACCAGCTGAGGACCGAGCTCACCAACATCGACCACAGCTACCTTAAGGGCATTGACATCGCCTCCAAGAAGAAAACCAGCAACGTGGAGAACACCGAGTTCGAGGCCATCAGCACAAAGATCTTCACCGACGAGCTGGGCTTCAGTGGCAAACACTTAGGGGGCAGCAATAAGCCTGACGGACTTTTGTGGGACGATGACTGCGCCATCATTCTGGACAGCAAGGCTTACAGCGAGGGATTCCCTCTGACAGCCAGCCATACAGACGCTATGGGAAGGTACCTGAGGCAGTTCACAGAGCGCAAGGAGGAGATCAAGCCTACATGGTGGGACATTGCTCCTGAGCACCTGGACAACACCTACTTCGCCTACGTGTCTGGCAGCTTTAGCGGCAACTATAAGGAGCAGCTGCAGAAGTTTAGACAGGACACCAACCACCTTGGCGGAGCCTTGGAGTTCGTGAAGCTGTTACTTCTGGCTAACAACTACAAGACCCAGAAGATGAGCAAGAAGGAGGTCAAGAAGTCCATCCTGGACTACAACATCAGCTACGAGGAGTACGCCCCCTTGCTCGCTGAGATCGAGTGA')
	BbvI = Seq('AGTGGCTCAGTTGTCAAGGCCTGATCCAGCATTGGCCGCACTGACCAACGATCACCTGGTAGCTTTGGCTTGTCTTGGAGGTAGGCCCGCCTTGGATGCTGTTAAGAAGGGATTGCCTCATGCTCCTGCCTTGATTAAGAGGACCAACAGGCGGATTCCTGAGAGAACCTCTCATAGGGTGGCTGGCTCTGCTTTCTCTGGCTTGAACACCGAGTTTGTGGAGTCCTACCTGTACGAGGAGGATAACCTGAGGTTCGAGGACAAGACAGGCGAGGTCCTGAAGGCTATCGGCTTTGACGTGGAAATGAGGCCTAAGCCCGCTTCTATGGAGAGGACAGAGATCGAGATCATGGTGAAGTACGGCGACAGGCAGTGTGGCATCATCGACGCTAAGAACTACAGACAGAAGTTCGCCCTGAGCGCCAGCTTGACTTCTCATATGGCTTCCGAGTACATCCCCAACTATCAGGGCTACAAGGGCCTGAACGTGCAGTTCTTTGGATACGTGACAGCTGCCGATTTCAGCGGCGAGAAGAACCTGGAGAAGATCAGCAACAAGGTGCAGGAACACACCAGCAGCAGGGACATCAAGGGCTTGATGCTAAGCGCCAAGGTTTTGCTGGGCTTCCTGGACTACTGCCTCGAGAATGACATACCTGAGAACGAGAGGGTGAACCTGTTCATCAGAGCCGTGCAGAACAGGGGGTACAAGACCTTGGGAGAGATGTTGAAGGAGGCCAAGTACTGA')
#Previously had the StsI below.  Check why it's shorter?
#	StsI = Seq('GTGCTGGAGAAGAGCGACATCGAGAAGTTCAAGAACCAGCTGCGCACCGAGCTGACCAACATCGACCACAGCTACCTGAAGGGCATCGACATCGCCAGCAAGAAGAAGACCAGCAACGTGGAGAACACCGAGTTCGAGGCCATCAGCACCAAGATCTTCACCGACGAGCTGGGCTTCAGCGGCAAGCACCTGGGCGGCAGCAACAAGCCCGACGGCCTGCTGTGGGACGACGACTGCGCCATCATCCTGGACAGCAAGGCCTACAGCGAGGGCTTCCCCCTGACCGCCAGCCACACCGACGCCATGGGCCGCTACCTGCGCCAGTTCACCGAGCGCAAGGAGGAGATCAAGCCCACCTGGTGGGACATCGCCCCTGAACACCTGGACAACACCTACTTCGCCTACGTGAGCGGCAGCTTCAGCGGCAACTACAAGGAGCAGCTGCAGAAGTTCCGCCAGGACACCAACCACCTGGGCGGCGCCCTGGAGTTCGTGAAGCTGCTGCTGCTGGCCAACAACTACAAGACCCAGAAGATGAGCAAGAAGGAGGTGAAGAAGAGCATCCTGGACTACAACATCAGCTACTAA',generic_dna)
	#RVD sequences
	RVDtoSequence = {	'NI':Seq('AACATC',generic_dna),
						'NG':Seq('AACGGA',generic_dna),
						'HD':Seq('CATGAC',generic_dna),
						'NN':Seq('AACAAC',generic_dna) }
	RVD_common_S = Seq('GTCGTGGCAATTGCGAGC',generic_dna)
	
	RVD_common_E = Seq('GGGGGAAAGCAGGCACTCGAAACCGTCCAGAGGTTGCTGCCTGTGCTGTGCCAAGCGCACGG',generic_dna)

	#variable parts of repeat sequence, original Sanjana primers
	Posn3 = Seq('ACTTACGCCAGAGCAG',generic_dna) #Also used for position 9,15
	Posn4 = Seq('ACTAACCCCAGAGCAG',generic_dna) #Also used for position 10,16
	Posn5 = Seq('GTTGACCCCAGAGCAG',generic_dna) #Also used for position 11,17
	Posn6 = Seq('CCTGACCCCAGAGCAG',generic_dna) #Also used for position 12,18
	Posn7 = Seq('ACTGACACCAGAGCAG',generic_dna) #Also used for position 13,19
	Posn2 = Posn6 #Position 2 is the same in Sanjana
	Posn8 = Seq('ACTTACACCCGAACAA',generic_dna)
	Posn14 = Seq('CCTCACCCCAGAGCAG',generic_dna)
	#variable parts of repeat sequence, Riboslice v2 primers (these are overwriting the
	#Sanjana primers

	#I.e. riboslice RVDs are used
	if riboslice_v2:
#		print('in riboslice v_2')
		#riboslice_50percent != 2 means that we're have riboslice in the odd positions (==1 (rs in even positions) or ==0 (rs in all positions))
		if riboslice_50percent == 1:
			Posn3 = Seq('AGGACTTACGCCAGAGCAG',generic_dna) #Also used for position 9,15
			Posn5 = Seq('AGGGTTGACCCCAGAGCAG',generic_dna) #Also used for position 11,17
			Posn7 = Seq('AGGACTGACACCAGAGCAG',generic_dna) #Also used for position 13,19
#			print('changing positions 2,4,6,8,10,12,14,16,18')

		if riboslice_50percent == 2:
			Posn4 = Seq('AGGACTAACCCCAGAGCAG',generic_dna) #Also used for position 10,16
			Posn6 = Seq('AGGCCTGACCCCAGAGCAG',generic_dna) #Also used for position 2,12,18
			Posn2 = Seq('CCTGACCCCAGAGCAG',generic_dna)
			Posn8 = Seq('AGGACTTACACCCGAACAA',generic_dna)
			Posn14 = Seq('CGGCCTCACCCCAGAGCAG',generic_dna)
#			print('changing positions 3,5,7,9,11,13,15,17,19')


		RVD_common_E = Seq('GGGGGAAAGCAGGCACTCGAAACCGTCCAGAGGTTGCTGCCTGTGCTGTGCCAAGCGGGACACGG',generic_dna)
		RVD_common_Normal = Seq('GGGGGAAAGCAGGCACTCGAAACCGTCCAGAGGTTGCTGCCTGTGCTGTGCCAAGCGCACGG',generic_dna)

	#NB Code below is for the version with SG added
	# if riboslice_v2:
	# 	Posn3 = Seq('ATCAGGACTTACGCCAGAGCAG',generic_dna) #Also used for position 9,15		(AHGSG)
	# 	Posn4 = Seq('ATCAGGACTAACCCCAGAGCAG',generic_dna) #Also used for position 10,16 	(AHGSG)
	# 	Posn5 = Seq('GTCAGGGTTGACCCCAGAGCAG',generic_dna) #Also used for position 11,17 	(AHGSG)
	# 	Posn6 = Seq('GTCAGGCCTGACCCCAGAGCAG',generic_dna) #Also used for position 2,12,18 	(AHGSG)
	# 	Posn7 = Seq('ATCAGGACTGACACCAGAGCAG',generic_dna) #Also used for position 13,19 	(AHGSG)
	# 	Posn2 = Seq('CCTGACCCCAGAGCAG',generic_dna) #This is unchanged versus Sanjana 		(AHGSG)
	# 	Posn8 = Seq('ATCAGGACTTACACCCGAACAA',generic_dna) 									(AHGSG)
	# 	Posn8 = Seq('ATCAGGACTTACACCCGAACAA',generic_dna) 									(AHGSG)
	# 	Posn8 = Seq('ATCAGGACTTACACCCGAACAA',generic_dna) 									(AHGSG)
	# 	Posn14 = Seq('CTCCGGCCTCACCCCAGAGCAG',generic_dna) 									(AHGSG)

	#1/2 RVD sequences
	#Strictly speaking, the first A in halfRVD_S isn't part of the halfRVD (in terms of
	#codons), but it's easier to put it here.
	halfRVD_S = Seq('ACTCACGCCTGAGCAGGTAGTGGCTATTGCATCC',generic_dna)
	if riboslice_v2 and riboslice_50percent == 2:
		halfRVD_S = Seq('AGGACTCACGCCTGAGCAGGTAGTGGCTATTGCATCC',generic_dna) #GHGG version
		#halfRVD_S = Seq('ATCAGGACTCACGCCTGAGCAGGTAGTGGCTATTGCATCC',generic_dna) %SG version
	halfRVD_E = Seq('GGGGGCAGACCCGCACTGGAG',generic_dna)


	repeatRegion = Seq('',generic_dna)
	rvdString = ''
	for posn,base in enumerate(recognitionSequence):
		posn = posn+1 #Index from 1
		base = base.upper()
		#First verify that the sequence begins with a T
		if posn==1 and base!='T':
			print("Recognition sequence doesn't start with a T, quitting.")
			break
		#Determine the DNA representation of the RVD
		RVDSeq = RVDtoSequence[baseToRVDs[base]]
		#Determine what the variable part of the repeat sequence should be (this is required
		#due to the codon modification done to enable golden gate assembly and sequencing)
		if posn%6==0:
			firstRepeatSeqPart = Posn6
		elif posn%6==1:
			firstRepeatSeqPart = Posn7
		elif posn%6==2:
			if posn==2:
				firstRepeatSeqPart = Posn2
			elif posn==8:
				firstRepeatSeqPart = Posn8
			elif posn==14:
				firstRepeatSeqPart = Posn14
		elif posn%6==3:
			firstRepeatSeqPart = Posn3
		elif posn%6==4:
			firstRepeatSeqPart = Posn4
		elif posn%6==5:
			firstRepeatSeqPart = Posn5
		# if we're at the end of the sequence, only add the half RVD
		if posn == len(recognitionSequence):
			repeatRegion += halfRVD_S + RVDSeq + halfRVD_E
			rvdString+= "_" * len(halfRVD_S) + str(RVDSeq) + "_" * len(halfRVD_E)
		elif posn != 1:
			# otherwise add the full thing
			#XXYY00: Code for 50% substitutions
			if riboslice_v2==False:
				repeatRegion += firstRepeatSeqPart + RVD_common_S + RVDSeq + RVD_common_E
				rvdString+= "_" * len(firstRepeatSeqPart) + "_"*len(RVD_common_S) + \
						str(RVDSeq) + "_" * len(RVD_common_E)
			else:
				#These are the riboslice regions
				if (riboslice_50percent==1 and posn%2!=1) or (riboslice_50percent==2 and posn%2!=0) or riboslice_50percent==0:
						repeatRegion += firstRepeatSeqPart + RVD_common_S + RVDSeq + RVD_common_E
						rvdString+= "_" * len(firstRepeatSeqPart) + "_"*len(RVD_common_S) + \
								str(RVDSeq) + "_" * len(RVD_common_E)
#						print('Riboslice region ' + str(posn))
				#These are the talen-like regions
				else:
					repeatRegion += firstRepeatSeqPart + RVD_common_S + RVDSeq + RVD_common_Normal
					rvdString+= "_" * len(firstRepeatSeqPart) + "_"*len(RVD_common_S) + \
							str(RVDSeq) + "_" * len(RVD_common_Normal)
#					print('Normal region ' + str(posn))
	if nucleaseType.upper() == 'FOKI':
		NucleaseDomain = FokI
	elif nucleaseType.upper() == 'FOKI_HE':
		NucleaseDomain = FokI_HE
	elif nucleaseType.upper() == 'STSI':
		NucleaseDomain = StsI
	elif nucleaseType.upper() == 'BBVI':
		NucleaseDomain = BbvI

	#Generate the output sequence from the different fragments, in pIDT if that was desired
	outputSeq = five_prime_UTR + N_terminus + repeatRegion
	if nucleaseType.upper() == 'FOKI':
		outputSeq += C_terminus + Linker + NucleaseDomain + three_prime_UTR
	elif nucleaseType.upper() == 'FOKI_HE':
		outputSeq += C_terminus + Linker + NucleaseDomain + three_prime_UTR
	# elif nucleaseType.upper() == 'STSI':
	# 	outputSeq += C_terminus_StsI + Linker + NucleaseDomain + three_prime_UTR
	else:
		outputSeq += C_terminus + Linker + NucleaseDomain + three_prime_UTR
	if include_pIDT > 0:
		outputSeq = pIDT_fivePrime + outputSeq + pIDT_threePrime

	outputRVD = " "*(len(five_prime_UTR) + len(N_terminus)) + rvdString
	if include_pIDT > 0:
		outputRVD = " "*(len(pIDT_fivePrime)) + outputRVD

	return [outputSeq, outputRVD]

def generateAlignmentFragments(outputSeq,alignmentLevel=0):
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
	#Constants: Sequencing primers
	M13F_m21 = Seq('GTAAAACGACGGCCAGT',generic_dna)
	M13R = Seq('CAGGAAACAGCTATGAC',generic_dna)
	F1 = Seq('CCAGTTGCTGAAGATCGCGAAGC',generic_dna)
	F2 = Seq('ACTTACACCCGAACAAGTCG',generic_dna)
	R1 = Seq('TGCCACTCGATGTGATGTCCTC',generic_dna)
	R2 = Seq('CTGCCCGATGGGAAGATTGTAC',generic_dna)


	#Determine the position of the sequencing primers, and slice out a 1 kb region following
	#those primers to do the individual alignments to
	M13FLoc = outputSeq.find(M13F_m21)
	M13RLoc = outputSeq.find(M13R.reverse_complement())+len(M13R)
	F1Loc = outputSeq.find(F1)
	F2Loc = outputSeq.find(F2)
	R1Loc = outputSeq.find(R1.reverse_complement())+len(R1)
	R2Loc = outputSeq.find(R2.reverse_complement())+len(R2)

	S0 = outputSeq[M13FLoc:M13FLoc+1100]
	S1 = outputSeq[F1Loc:F1Loc+1100]
	S2 = outputSeq[F2Loc:F2Loc+1100]
	S3 = outputSeq[R1Loc-1100:R1Loc].reverse_complement()
	S4 = outputSeq[R2Loc-1100:R2Loc].reverse_complement()
	S5 = outputSeq[M13RLoc-1100:M13RLoc].reverse_complement()
	
	if alignmentLevel == 0:
		return [S1,S2,S3]
	elif alignmentLevel == 1:		
		return [S0,S1,S2,S3,S5]
	elif alignmentLevel == 2:
		return [S0,S1,S2,S3,S4,S5]

def generateAlignments(filePrefix,alignmentReferences,inputDir,alignmentLevel=0):
	from Bio import SeqIO
	if alignmentLevel == 2:
		fileSuffixes = ['M13F','F1','F2','R1','R2','M13R']
	elif alignmentLevel == 1:
		fileSuffixes = ['M13F','F1','F2','R1','M13R']
		#Genewiz changed the M13F primer from 'M13F(-21)'' to 'M13F' so we don't need to do the renaming anymore
#		fileSuffixes = ['M13F-21','F1','F2','R1','M13R']
		# #We need to rename (actually, we create a copy) of the 'M13F(-21)' file with the
		# #suffix M13F-21, otherwise clustalw won't process it.
		# import shutil
		# source = './inputFiles/' + filePrefix + '-M13F(-21).seq'
		# dest = './inputFiles/' + filePrefix + '-M13F-21.seq'
		# shutil.copyfile(source,dest)
	else:
		fileSuffixes = ['F1','F2','R1']
	for i,suf in enumerate(fileSuffixes):
		fileName = filePrefix + '-' + suf + '.seq'	
		fi = open(os.path.join(inputDir,fileName),'r')
		mySeq = list(SeqIO.parse(fi, "fasta"))[0].seq
		fi.close()
		outText = '>'+ filePrefix + '-' + suf + '_REFERENCE\n' + alignmentReferences[i] + '\n'
		outText += '>' + fileName + '\n' + mySeq
		#Generate the file to be used as input to clustalw
		alignmentInputFileN = './tempFiles/' + fileName
		fo = open(alignmentInputFileN,'w')
		fo.write(str(outText))
		fo.close()
		#Get the command line for running clustalw and run it
		from Bio.Align.Applications import ClustalwCommandline
		clustalw_cline = ClustalwCommandline('clustalw2',infile=alignmentInputFileN)
		x=clustalw_cline()

def parseAlignments(filePrefix,alignmentLevel=0):
	from Bio import AlignIO
	if alignmentLevel==2:
		fileSuffixes = ['M13F','F1','F2','R1','R2','M13R']
	elif alignmentLevel == 1:
#		fileSuffixes = ['M13F-21','F1','F2','R1','M13R']
		fileSuffixes = ['M13F','F1','F2','R1','M13R']
	else:
		fileSuffixes = ['F1','F2','R1']
	allAligns = []
	for i,suf in enumerate(fileSuffixes):
		fileName = filePrefix + '-' + suf + '.aln'
		align = AlignIO.read('./tempFiles/'+fileName,'clustal')
		allAligns.append([align[0],align[1],align._star_info])
	return allAligns

def printAlignments(filePrefix,outputDir,outputSeq,outputRVD,allAlignments,alignmentLevel=False,
					preambleText=""):
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
	#Constants: Sequencing primers
	M13F_m21 = Seq('GTAAAACGACGGCCAGT',generic_dna)
	M13R = Seq('CAGGAAACAGCTATGAC',generic_dna)
	F1 = Seq('CCAGTTGCTGAAGATCGCGAAGC',generic_dna)
	F2 = Seq('ACTTACACCCGAACAAGTCG',generic_dna)
	R1 = Seq('TGCCACTCGATGTGATGTCCTC',generic_dna)
	R2 = Seq('CTGCCCGATGGGAAGATTGTAC',generic_dna)

	#Determine the position of the sequencing primers
	M13FLoc = outputSeq.find(M13F_m21)
	M13RLoc = outputSeq.find(M13R.reverse_complement())+len(M13R)
	F1Loc = outputSeq.find(F1)
	F2Loc = outputSeq.find(F2)
	R1Loc = outputSeq.find(R1.reverse_complement())+len(R1)
	R2Loc = outputSeq.find(R2.reverse_complement())+len(R2)

	if alignmentLevel == 0:
		S1 = allAlignments[0][0].seq
		S1b = allAlignments[0][1].seq
		S1s = allAlignments[0][2]
		S2 = allAlignments[1][0].seq
		S2b = allAlignments[1][1].seq
		S2s = allAlignments[1][2]
		S3 = allAlignments[2][0].seq.reverse_complement()
		S3b = allAlignments[2][1].seq.reverse_complement()
		S3s = allAlignments[2][2][::-1] # reverse the string

		#Determine the longest consecutive sequence of matched alignments
		#First create equal length strings, properly aligned wrt the reference sequence
		St1 = ' ' * F1Loc + S1s + ' ' * (len(outputSeq)-len(S1s)-F1Loc)
		St2 = ' ' * F2Loc + S2s + ' ' * (len(outputSeq)-len(S2s)-F2Loc)
		St3 = ' ' * (R1Loc-len(S3s)) + S3s + ' ' * (len(outputSeq)-R1Loc)
		#Sum up the '*'s at each position
		import operator
		S=''.join([str(operator.countOf(tup,'*')) for tup in zip(St1,St2,St3)])
		#Find the longest stretch of non-zero elements
		maxLength = max([len(x) for x in S.split('0')])
	

	elif alignmentLevel == 1:
		S0 = allAlignments[0][0].seq
		S0b = allAlignments[0][1].seq
		S0s = allAlignments[0][2]
		S1 = allAlignments[1][0].seq
		S1b = allAlignments[1][1].seq
		S1s = allAlignments[1][2]
		S2 = allAlignments[2][0].seq
		S2b = allAlignments[2][1].seq
		S2s = allAlignments[2][2]
		S3 = allAlignments[3][0].seq.reverse_complement()
		S3b = allAlignments[3][1].seq.reverse_complement()
		S3s = allAlignments[3][2][::-1] # reverse the string
		S4 = allAlignments[4][0].seq.reverse_complement()
		S4b = allAlignments[4][1].seq.reverse_complement()
		S4s = allAlignments[4][2][::-1] # reverse the string

		St0 = ' ' * M13FLoc + allAlignments[0][2] + ' ' * (len(outputSeq)-len(allAlignments[0][2])-M13FLoc)
		St1 = ' ' * F1Loc + allAlignments[1][2] + ' ' * (len(outputSeq)-len(allAlignments[1][2])-F1Loc)
		St2 = ' ' * F2Loc + allAlignments[2][2] + ' ' * (len(outputSeq)-len(allAlignments[2][2])-F2Loc)
		St3 = ' ' * (R1Loc-len(allAlignments[3][2])) + allAlignments[3][2][::-1] + ' ' * (len(outputSeq)-R1Loc)
		St4 = ' ' * (M13RLoc-len(allAlignments[4][2])) + allAlignments[4][2][::-1] + ' ' * (len(outputSeq)-M13RLoc)

		import operator
		S=''.join([str(operator.countOf(tup,'*')) for tup in zip(St0,St1,St2,St3,St4)])	
		maxLength = max([len(x) for x in S.split('0')])

	elif alignmentLevel == 2:
		S0 = allAlignments[0][0].seq
		S0b = allAlignments[0][1].seq
		S0s = allAlignments[0][2]
		S1 = allAlignments[1][0].seq
		S1b = allAlignments[1][1].seq
		S1s = allAlignments[1][2]
		S2 = allAlignments[2][0].seq
		S2b = allAlignments[2][1].seq
		S2s = allAlignments[2][2]
		S3 = allAlignments[3][0].seq.reverse_complement()
		S3b = allAlignments[3][1].seq.reverse_complement()
		S3s = allAlignments[3][2][::-1] # reverse the string
		S4 = allAlignments[4][0].seq.reverse_complement()
		S4b = allAlignments[4][1].seq.reverse_complement()
		S4s = allAlignments[4][2][::-1] # reverse the string
		S5 = allAlignments[5][0].seq.reverse_complement()
		S5b = allAlignments[5][1].seq.reverse_complement()
		S5s = allAlignments[5][2][::-1] # reverse the string


		St0 = ' ' * M13FLoc + allAlignments[0][2] + ' ' * (len(outputSeq)-len(allAlignments[0][2])-M13FLoc)
		St1 = ' ' * F1Loc + allAlignments[1][2] + ' ' * (len(outputSeq)-len(allAlignments[1][2])-F1Loc)
		St2 = ' ' * F2Loc + allAlignments[2][2] + ' ' * (len(outputSeq)-len(allAlignments[2][2])-F2Loc)
		St3 = ' ' * (R1Loc-len(allAlignments[3][2])) + allAlignments[3][2][::-1] + ' ' * (len(outputSeq)-R1Loc)
		St4 = ' ' * (R2Loc-len(allAlignments[4][2])) + allAlignments[4][2][::-1] + ' ' * (len(outputSeq)-R2Loc)
		St5 = ' ' * (M13RLoc-len(allAlignments[5][2])) + allAlignments[5][2][::-1] + ' ' * (len(outputSeq)-M13RLoc)

		import operator
		S=''.join([str(operator.countOf(tup,'*')) for tup in zip(St0,St1,St2,St3,St4,St5)])	
		maxLength = max([len(x) for x in S.split('0')])
			
	textOut= 	preambleText + "Complete Sequence:\t" + outputSeq + "\n" + \
				"RVDs:\t\t\t\t" + outputRVD + "\n"
	if alignmentLevel > 0:
		textOut += 	"F0:\t\t\t\t\t" + str(S0).rjust(M13FLoc+len(S0)) + "\n" + \
					"S0:\t\t\t\t\t" + str(S0b).rjust(M13FLoc+len(S0)) + "\n" + \
					"A0:\t\t\t\t\t" + S0s.rjust(M13FLoc+len(S0)) + "\n"
	
	textOut +=	"F1:\t\t\t\t\t" + str(S1).rjust(F1Loc+len(S1)) + "\n" + \
				"S1:\t\t\t\t\t" + str(S1b).rjust(F1Loc+len(S1)) + "\n" + \
				"A1:\t\t\t\t\t" + S1s.rjust(F1Loc+len(S1)) + "\n" + \
				"F2:\t\t\t\t\t" + str(S2).rjust(F2Loc+len(S2)) + "\n" + \
				"S2:\t\t\t\t\t" + str(S2b).rjust(F2Loc+len(S2)) + "\n" + \
				"A2:\t\t\t\t\t" + S2s.rjust(F2Loc+len(S2)) + "\n" + \
				"R1:\t\t\t\t\t" + str(S3).rjust(R1Loc) + "\n" + \
				"S3:\t\t\t\t\t" + str(S3b).rjust(R1Loc) + "\n" + \
				"A3:\t\t\t\t\t" + S3s.rjust(R1Loc) + "\n"

	if alignmentLevel == 2:
		textOut +=	"R2:\t\t\t\t\t" + str(S4).rjust(R2Loc) + "\n" + \
					"S4:\t\t\t\t\t" + str(S4b).rjust(R2Loc) + "\n" + \
					"A4:\t\t\t\t\t" + S4s.rjust(R2Loc) + "\n"
		textOut +=	"M13R:\t\t\t\t" + str(S5).rjust(M13RLoc) + "\n" + \
					"S5:\t\t\t\t\t" + str(S5b).rjust(M13RLoc) + "\n" + \
					"A5:\t\t\t\t\t" + S5s.rjust(M13RLoc) + "\n"
				
	if alignmentLevel == 1:
		textOut +=	"M13R:\t\t\t\t" + str(S4).rjust(M13RLoc) + "\n" + \
					"S4:\t\t\t\t\t" + str(S4b).rjust(M13RLoc) + "\n" + \
					"A4:\t\t\t\t\t" + S4s.rjust(M13RLoc) + "\n"
	
	textOut += "Matches:\t\t\t" + S + "\n" + "Longest continuous match: " + str(maxLength)
	
	f = open(os.path.join(outputDir,filePrefix+'_alignment.txt'),'w')
	f.write(str(textOut))
	f.close()

	return(maxLength)

if __name__ == '__main__':
	main()