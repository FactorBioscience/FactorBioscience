import unittest
import generalAlignment as ga
import GAfileProcessing as gaFP

class testgeneralAlignment(unittest.TestCase):

	def testCompareAndMergeReferenceDicts(self):
		print('X')
		#compareAndMergeReferenceDicts(masterReferences,localReferences)
		testMasterDict = {'itemA':1,'itemB':2,'itemC':'X','itemD':'Y'}
		testLocalDict = {'item1': 'A', 'item2':'B', 'itemA':100}
		expectedCombinedDict = {'itemA':100,'itemB':2,'itemC':'X','itemD':'Y', 'item1': 'A', 'item2':'B'}
		self.assertEqual(gaFP.compareAndMergeReferenceDicts(testMasterDict,{}),testMasterDict)
		self.assertEqual(gaFP.compareAndMergeReferenceDicts({},testLocalDict),testLocalDict)
		self.assertEqual(gaFP.compareAndMergeReferenceDicts(testMasterDict,testLocalDict),expectedCombinedDict)


	def testParseGeneralAlignmentSequencesFile(self):
		print('X')
		#parseGeneralAlignmentSequencesFile(infile)

	def testGetFilesWithSequencesFileEntries(self):
		print('X')
		#getFilesWithSequencesFileEntries(seqFileSet,sequencesTxtEntries)

	def test_scoreFragments(self):
		#First test the normal usage of scoreFragments
		reference = 'GGTTAACCGTGA'
		fragments = [
						['GGTTAACCGTGA',
						 'GGTTAACCGTCA',
						 '********** *'],
						['GGTTAACCGTGA',
						 'CGCTCAACCTCG',
						 ' * * * * *  '],
						['GGTTAA      ',
						 'GGTTAA      ',
						 '******      ']]
		scoreStr = '232323121201'
		self.assertEqual(ga.scoreFragments(reference,fragments),scoreStr)
		#Now ensure that the routine returns None if any of the input fragments differ in length
		fragmentsShort = [
						['GGTTAACCGTG',
						 'GGTTAACCGTCA',
						 '********** *']]
		self.assertIsNone(ga.scoreFragments(reference,fragmentsShort))

	def testCalculateLongestRunWithoutChar(self):
		testStrA = '1231232234212122132231231231231231231121212121'
		testStrB = '1231032204210212021302230123012301230123012301120121021201'
		testStrC = '1031232234212122132231231231231231231121212121'
		testStrD = ''
		self.assertEqual(ga.calculateLongestRunWithoutChar(testStrA,'0'),[46,46])
		self.assertEqual(ga.calculateLongestRunWithoutChar(testStrB,'0'),[4,44])
		self.assertEqual(ga.calculateLongestRunWithoutChar(testStrC,'0'),[44,45])
		self.assertEqual(ga.calculateLongestRunWithoutChar(testStrD,'0'),[0,0])





if __name__ == '__main__':
    unittest.main()