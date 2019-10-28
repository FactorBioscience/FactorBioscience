import pathlib

DEFAULT_TAIL = 'TT-BsaI'
DEFAULT_5UTR = '5Hbb'
DEFAULT_3UTR = '3Hbb'

"""
mySeq_GeneA	GeneA #Defaults to TT-BsaI tails
mySeq_GeneB	GeneB TT-BsaI
mySeq_GeneC	GeneC TT-BsaI 5Hbb|3wpreA
"""

#Process the input (Sequences.txt) file and return a tuple of the format
#(name, reference sequence name, tail type, 5'UTR type, 3'UTR type)
def parseGeneralAlignmentSequencesFile(infile):
	#check for the file
	if not infile.exists():
		print('Error: Input file "' + str(infile.relative_to(pathlib.Path.cwd())) + '" does not exist.')
		return None

	#load the file and parse out the items
	elementsToAlign = []
	with open(infile) as f:
		for i,line in enumerate(f):
			noComments = line.split('#',maxsplit=1)[0]
			args = noComments.split()

			if len(args) < 2:
				print('Warning: Line ' + str(i) + ' in file "' + str(infile.name) + '" does not contain a reference sequence:')
				print(str(i)+":\t"+line)
				continue
			if len(args) > 4:
				print('Warning: Line ' + str(i) + ' in file "' + str(infile.name) + '" contains too many items. It may be incorrectly processed:')
				print(str(i)+":\t"+line)
			UTRs = args[3].split('|') if len(args) > 3 else [DEFAULT_5UTR, DEFAULT_3UTR]
			if len(UTRs) != 2:
				print('Error processing UTRs in line ' + str(i) + ' in file "' + str(infile.name) + '"')
				print('Reverting to defaults (' + DEFAULT_5UTR + " | " + DEFAULT_3UTR + ')')
				UTRs = [DEFAULT_5UTR, DEFAULT_3UTR]
			tail = args[2] if len(args) > 2 else DEFAULT_TAIL
			elementsToAlign.append((args[0],args[1],tail,UTRs[0],UTRs[1]))
	return elementsToAlign

def compareAndMergeReferenceDicts(masterReferences,localReferences):
	sharedItems = [x for x in localReferences if x in masterReferences]
	if sharedItems:
		print('The following items were found in both the local and global reference files:')
		for item in sharedItems:
			print('\t'+item)
		print('Will use local references.')
	masterReferences.update(localReferences)
	return masterReferences

def getFilesWithSequencesFileEntries(seqFileSet,sequencesTxtEntries):
	#To handle cases like mySeq_2-F1.seq being an instance of mySeq in the Sequences.txt file, we will
	#first check to see if mySeq_2 is in seqFileSet, then check if mySeq is in it.
	#Note that this requires that an underscore be used to identify the clone number.
	sequencesTxtDict = {el[0]:[el[1],el[2],el[3],el[4]] for el in sequencesTxtEntries}

	filesWithSequencesFileEntries = []
	for fileGroup in sorted(seqFileSet):
		shorterName = fileGroup.rsplit('_',maxsplit=1)[0]
		if fileGroup in sequencesTxtDict.keys():
			el = sequencesTxtDict[fileGroup]
			filesWithSequencesFileEntries.append([fileGroup, el[0], el[1], el[2], el[3]])
		elif shorterName in sequencesTxtDict.keys():
			el = sequencesTxtDict[shorterName]
			filesWithSequencesFileEntries.append([fileGroup, el[0], el[1], el[2], el[3]])
		else:
			print('Warning. File group ' + fileGroup + ' lacks an entry in Sequences.txt. Cannot align.')

	#Determine the entries in Sequences.txt that are lacking a sequence file
	seqFileSetRootNames = {el.rsplit('_',maxsplit=1)[0] for el in seqFileSet}
	for entry in sequencesTxtEntries:
		if entry[0] not in seqFileSet and entry[0] not in seqFileSetRootNames:
			print('Warning. Sequences.txt entry ' + entry[0] + ' does not appear to have any corresponding .seq files.')

	return filesWithSequencesFileEntries
