#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information:
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

'''Module containing utility classes and functions'''

import os.path, getopt, sys, traceback, subprocess
import Exceptions, Environment, Data

allResidues = ['ALA', 'ARG', 'ASP', 'ASN', 'CYS', 'GLY', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'PHE', 'SER', 'THR', 'TYR', 'TRP', 'VAL']
aminoAcidCodes = {'ALA':'A', 'ASP':'D', 'ARG':'R', 'ASN':'N', 'CYS':'C' , 'GLY':'G', 'GLN':'Q', 'GLU':'E', 'HIS':'H', 'ILE':'I', 
		'LYS':'K', 'LEU':'L', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

weights = {'ARG': 174, 'HIS': 155, 'LYS': 146, 'ASP': 133, 'GLU': 147, 'SER': 105,
	'THR': 119, 'ASN': 132, 'GLN': 146, 'CYS': 121, 'GLY': 75, 'PRO': 115,
	'ALA': 89, 'ILE': 131, 'LEU': 131, 'MET': 149, 'PHE': 165, 'TRP': 204, 'TYR': 181, 'VAL': 117}
	
volumes = {'ARG': 173, 'HIS': 153, 'LYS': 169, 'ASP': 111, 'GLU': 138, 'SER': 89,
	'THR': 116, 'ASN': 114, 'GLN': 144, 'CYS': 109, 'GLY': 60, 'PRO': 113,
	'ALA': 89, 'ILE': 167, 'LEU': 167, 'MET': 163, 'PHE': 190, 'TRP': 228, 'TYR': 194, 'VAL': 140}
	
surfaceAreas = {'ARG': 225, 'HIS': 195, 'LYS': 200, 'ASP': 150, 'GLU': 190, 'SER': 115,
	'THR': 140, 'ASN': 160, 'GLN': 180, 'CYS': 135, 'GLY': 75, 'PRO': 145,
	'ALA': 115, 'ILE': 175, 'LEU': 170, 'MET': 185, 'PHE': 210, 'TRP': 255, 'TYR': 230, 'VAL': 155}	

def WeightChangeForMutation(mutationCode):

	'''Calculates the change in Molecular Weight due to a mutation
	
	Params:
		mutationCode - must be in extended format (contain the wild-type residue code)
		
	Returns
		The difference in mw between new and old residues i.e. (New - Old)	
	'''	

	if IsReducedCode(mutationCode):
		mutationCode = ConvertMutationCodeFromReducedFormat(mutationCode, extended=True)

	components = ParseMutationCode(mutationCode)
	if(len(components) != 4):
		raise ValueError, 'Mutation code must be extended to determine weight change'
		
	return weights[components[3]] - weights[components[2]]
	
def VolumeChangeForMutation(mutationCode):

	'''Calculates the change in Molecular Weight due to a mutation
	
	Params:
		mutationCode - must be in extended format (contain the wild-type residue code)
		
	Returns
		The difference in mw between new and old residues i.e. (New - Old)	
	'''	

	if IsReducedCode(mutationCode):
		mutationCode = ConvertMutationCodeFromReducedFormat(mutationCode, extended=True)

	components = ParseMutationCode(mutationCode)
	if(len(components) != 4):
		raise ValueError, 'Mutation code must be extended to determine weight change'
		
	return volumes[components[3]] - volumes[components[2]]	


def CreateBarChartData(matrix, column, **options):

	'''Returns a binary string representing a bar chart image with the given options'''

	print options
	options["matrix"] = matrix
	options["column"] = column
	options["outputFile"] = 'picture.png'
	options["xaxis"] = None
	CreateBarChart(**options)
	print 'Reading chart'
	stream = open('picture.png')
	data = stream.read()
	stream.close()
	os.remove('picture.png')
	print 'Done'
	
	return data

def CreateBarChart(matrix, column, outputFile, xaxis=None, **options):

	'''Creates a bar chart out of the specified data
	
	Parameters:
		Matrix: A Matrix.Matrix instance containing the data to create the chart from
		column: The name of the column containing the value (yAxis) data
		outputFile: The file to write the chart to
		xaxis: The name of the column containing the xaxis values.
			If None the row index of each y value is used.
	
	Keyword Args (All currently required):
		title
		xticIncrement
		xlabel
		ylabel
		'''

	string = '''set datafile separator ','
set xlabel '%s'
#set xtics nomirror %d
set ylabel '%s'
set term png size 1280,800
set output '%s'
set boxwidth 0.5 relative
set style fill solid 0.8
plot '%s' u %d:%d with boxes t '%s' '''
		
	stream = open('.matrix.csv', 'w+')
	stream.write(matrix.csvRepresentation())	
	stream.close()

	#Have to add 1 since gnuplot indexes from 1 while pythong indexes from 0
	columnIndex = matrix.indexOfColumnWithHeader(column) + 1

	#If no xColumn is specified set it to 0 i.e. row index
	if xaxis == None:
		xaxisIndex = 0
	else:
		xaxisIndex = matrix.indexOfColumnWithHeader(xaxis) + 1
		
	string = string % (options['xlabel'], int(options['xticIncrement']), options['ylabel'], 
				outputFile, ".matrix.csv", xaxisIndex, int(columnIndex), options['title']) 		
	stream = open('.temp', 'w+')
	stream.write(string)
	stream.close()

	try:
		if options.has_key('gnuplotPath'):
			gnuplotPath = options['gnuplotPath']
		elif Environment.appConfiguration is not None:
			gnuplotPath = Environment.appConfiguration.get('PATHS', 'gnuplot')
		else:
			raise BaseException
	except:
		gnuplotPath = '/usr/local/bin'
			
	process = subprocess.Popen("%s/gnuplot .temp" % gnuplotPath, shell=True)
	process.wait()
	
	#print 'Done creating graph'
	
	os.remove('.temp')
	os.remove('.matrix.csv')
	
	return process.returncode

def GetTraceback():

	'''Returns the current stack traceback as a string'''
	
	#Get the traceback as a string = convoluted ....
	trace = traceback.format_list(traceback.extract_tb(sys.exc_info()[2]))
	trace = reduce(lambda x, y: x + y, trace)
	#Need to remove single quotes from the traceback	
	trace = trace.replace("'", " ")

	return trace

def InvertedCodeDictionary():

	singleLetterCodes = [aminoAcidCodes[code] for code in allResidues]
	invertedCodeDict = dict(zip(singleLetterCodes, allResidues))

	return invertedCodeDict
	
def IsExtendedMutationFormat(code):

	'''Returns True if code is in extended format.
	
	That is, it contains both the residue and the mutation name
	
	Parameters:
		code: A mutation code - can be normal or reduced'''
		
	isExtended = True	
	if IsReducedCode(code):
		#Check if the second last element is a string or a number
		try:		
			element = code[-2]
			int(element)
			isExtended = False
		except ValueError:
			pass
		except IndexError:
			raise Exceptions.MutationCodeFormatError, 'Code %s does not have enough elements' % code
	else:
		#Check if there are four components separated by colons
		components = code.split(':')
		if len(components) == 3:
			isExtended = False
	
	return isExtended
	
def IsExtendedResidueFormat(code):

	'''Returns True if code is in extended format.
	
	That is, it contains both the residue number and the residue name
	
	Parameters:
		code: A residue code - can be normal or reduced'''
		
	isExtended = True	
	if IsReducedCode(code):
		#Check if the last element is a string or a number
		element = code[-1]
		try:
			int(element)
			isExtended = False
		except ValueError:
			pass
	else:
		#Check if there are three components separated by colons
		components = code.split(':')
		if len(components) == 2:
			isExtended = False
	
	return isExtended	
	
def IsReducedCode(code):

	'''Returns True is code is in reduced format, false otherwise.
	
	A code is in reduced format simply if it contains no colons'''
	
	components = code.split(':')
	if len(components) == 1:
		return True
	else:
		return False	
	
def CreateResidueCode(chain, number, name=None):

	'''Creates a residue code from a chain id and residue number and optionally a name
	
	A normal residue code has the form ChainID:ResidueIndex where the residue index is a
	string of length 4 left padded with zeros.
	For example passing A, 4 return A:0004.
	
	If a name is passed an extended residue code is returned. This has the format
	ChainID:ResidueIndex(Padded):ResidueName'''
	
	if name is None:
		return  chain + ":" + "0"*(4 - len(str(number))) + str(number)
	else:
		return chain + ":" + "0"*(4 - len(str(number))) + str(number) + ":" + name
	
def ParseStandardResidueCode(residueCode, extended=False): 

	'''Parses a residue code returning the chain and residueIndex (as an int).
	
	Set extended to True if the code is in extended format,
	See CreateResidueCode for information on the residue code format.
	
	Raises an Exceptions.ResidueCodeError if the residue code is not in the correct format'''
	
	components = residueCode.split(':')
	try:
		if components[1] != '0000':
			data = (components[0], int(components[1].lstrip('0'))) + tuple(components[2:])
		else:
			data = (components[0], int(components[1])) + tuple(components[2:])	
	except BaseException, data:
		print data
		raise Exceptions.ResidueCodeFormatError, 'Residue Code is %s. Invalid for standard format' % residueCode
	
	if len(data) > 3:
		raise Exceptions.ResidueCodeFormatError, 'Residue Code is %s. Invalid for standard format - too many components' % residueCode
	
	return data
	
def ParseReducedResidueCode(code, extended=False):

	'''Returns the chain, the residue index and the residue name'''
	
	invertedMap = InvertedCodeDictionary()

	noChain = False
	if not code[0].isalpha():
		noChain = True
		chainId = ""
		indexStart = 0

	try:
		if not noChain:
			chainId = code[0]
			indexStart = 1

		if extended is True:
			index = code[indexStart:-1]
			return chainId, int(index), invertedMap[code[-1]]
		else:
			index = code[indexStart:]
			return chainId, int(index)
	except ValueError:
		error = 'Residue code %s is invalid for reduced format. %s is not a valid residue index' % (code, index)	
		raise Exceptions.ResidueCodeFormatError, error		
		
def ParseResidueCode(code):

	'''Parses residue code returning chain, residue index and the residue name if present
	
	code can be in any format. However if the code has an error it may be misidentifed'''
																	
	extended = IsExtendedResidueFormat(code)
	if IsReducedCode(code):
		return ParseReducedResidueCode(code, extended=extended)
	else:
		return ParseStandardResidueCode(code, extended=extended)	
	
			
def ParseStandardMutationCode(mutationCode, extended=False):

	'''Parses an mutation code returning the chain, residueIndex (as an int) and mutationName (non-extended format)
	
	If the mutationCode is in extended format this function additionally returns the residueName before the mutationName
	See CreateExtendedResidueCode for information on the extended residue code format.
	
	Raises an IndexError if the residue code is not in the correct format'''
	
	mutationCode = mutationCode.strip(" ")
	
	components = mutationCode.split(':')
	try:
		if components[1] != '0000':
			data = (components[0], int(components[1].lstrip('0'))) + tuple(components[2:])
		else:
			data = (components[0], int(components[1])) + tuple(components[2:])	
	except BaseException, data:
		raise Exceptions.MutationCodeFormatError, 'Mutation Code is %s. Invalid for standard format' % mutationCode

	if len(data) > 4 or len(data) < 3:
		raise Exceptions.MutationCodeFormatError, 'Mutation Code is %s. Invalid for standard format - incorrect number of components' % mutationCode
	
	return data
	

def ParseReducedMutationCode(code, extended=False):

	'''Parse the reduced mutation code.
	
	If the reduced code includes the original residue name set extendedFormat to be true.
	If the code is not in extendedFormat this function returns the chain, the residue index and the mutation name
	in that order. 
	Otherwise it returns he chain, the residue index, the residue name and the mutation name'''

	try:
		invertedMap = InvertedCodeDictionary()
		mutation = invertedMap[code[-1]]
		data = ParseReducedResidueCode(code[:-1], extended=extended)
		return data + (mutation,)
	except Exceptions.ResidueCodeFormatError, data:
		raise Exceptions.MutationCodeFormatError, data	
	except Exception, data:
		raise Exceptions.MutationCodeFormatError, 'Error with mutation code %s. %s' % (code, data)	
	
def ParseMutationCode(code):

	'''Parses mutation code returning chain, residue index, residue name if present, and mutation code
	
	code can be in any format'''
																	
	extended = IsExtendedMutationFormat(code)
	if IsReducedCode(code):
		return ParseReducedMutationCode(code, extended=extended)
	else:
		return ParseStandardMutationCode(code, extended=extended)	
	
def ConvertMutationCodeToReducedFormat(mutationCode, extended=False):
	
	'''Converts a mutation code to a reduced format
	
	 e.g. A:0025:ALA:GLY -> A25AG'''
	 
	if extended is True: 
		chain, residueIndex, residueName, mutation = ParseMutationCode(mutationCode)
		code = chain + str(residueIndex) + aminoAcidCodes[residueName] + aminoAcidCodes[mutation]
	else:
		chain, residueIndex, mutation = ParseMutationCode(mutationCode)
		code = chain + str(residueIndex) + aminoAcidCodes[mutation]
		
	return code	
	
def ConvertMutationCodeFromReducedFormat(mutationCode, extended=False):

	'''Converts a mutation code from reduced format
	
	 e.g. A25AG -> A:0025:ALA:GLY'''
	 
	components = ParseReducedMutationCode(mutationCode, extended)
	if len(components) == 3:
		return CreateResidueCode(chain=components[0], number=components[1]) + ":" + components[2]
	else:
		return CreateResidueCode(chain=components[0], number=components[1], name=components[2]) + ":" + components[3]	
	

def ConvertResidueCodeToReducedFormat(residueCode, extended=False):

	'''Converts an residue code to a reduced format
	
	 e.g. A:0025 -> A25 or A:0025:ALA -> A25A'''

	if extended is True:
		chain, residueIndex, residueName = ParseResidueCode(residueCode)
		code =  chain + str(residueIndex) + aminoAcidCodes[residueName]
	else:
		chain, residueIndex = ParseResidueCode(residueCode)
		code = chain + str(residueIndex)
		
	return code	
	
	
def CheckFile(path):
	
	'''Checks if path refers to an existing file.
	
	Parameters:
		path - A path to a file - Can be absolute or not.
		
	Returns:
		The absolute path to the file.
		
	Exceptions:
		Raises an IOError if the file does not exist or if its not a file'''
		
	#Make the filename absolute
	if not os.path.isabs(path):
		path = os.path.abspath(path)
	
	#Check the specified file exists.
	if not os.path.exists(path):
		raise IOError, 'File %s does not exist' % path
		
	#Check its a file
	if not os.path.isfile(path):
		raise IOError, 'Object at %s is not a file' % path	
	
	return path

def CheckDirectory(path):
	
	#Make the filename absolute
	if not os.path.isabs(path):
		path = os.path.abspath(path)
	
	#Checking write access is a pain in python - just let any write operations raise exceptions
	if not os.path.isdir(path):
		raise IOError, "%s is not a directory" % path
		
	return path
	

def WhatIfCommandBoilerPlate(inputFile, renumber=False):
	
	'''Returns the following string
	
	getmol PDBName
	Y
	1
	set
	%renum
	tot
	1
	
	The last three lines are omitted by default.
	PDBName is the last path component of inputFile'''
	
	command='getmol '+os.path.split(inputFile)[1]+' \n Y \n 1\n set \n'
	if renumber:
		command=command+' %renumb \n tot \n 1 \n '
		
	return command	
	
def WhatIfMakeMoleculeCommand(inputFile):

	'''Returns the following string
	
	%makmol 
	 PDBName
	 PDBName.new
	 
	 tot 0
	 '''

	outputFile=os.path.split(inputFile)[1]+'.new'
	command=' %makmol \n '+os.path.split(inputFile)[1]+'\n'+outputFile+' \n \n tot 0 \n'
	
	return command, outputFile

def CleanPDB2PQR(inputFile, outputFile, forceField="amber", removeWater=True, removeLigand=True, removeAltLoc=True, addHydrogens=True, correct=True):

	'''Cleans a PDB by using PDB2PQR
	
	See CleanPDB for argument details

	Note: With pdb2pqr you cannot remove-water or ligands.
		
	Errors: 
		Raises an exception if the inputFile is not a valid PDB file.'''

	import Protool

	try:
		command = 'pdb2pqr --chain --ff=%s' % (forceField)

		if removeWater == True:
			print >>sys.stderr, 'Warn: Currently PDB2PQR does not remove waters from pdb files'

		if removeLigand == True:
			print >>sys.stderr, 'Warn: Currently PDB2PQR can not be used to remove heterogens from PDB files'

		if addHydrogens == False:
			print >>sys.stderr, 'Warn: Turning of Hydrogen addition with PDB2PQR automatically turns of rotamer correction'
			command = command + " --clean "
		else:
			if correct == False:
				command = command + " --nodebump "

		#Protool ignores altlocs so we can use it to remove them 
		#Do this first as Protool as when protool reads then writes a pdb2pqr cleaned file 
		#it raises an error on reading it again 
		if removeAltLoc is True:
			pdb = Protool.structureIO()
			pdb.readpdb(inputFile)
			pdb.writepdb(outputFile)

			inputFile = outputFile

		command = command + ' %s %s' % (inputFile, outputFile)
		print 'Using: ', command

		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		stdout, stderr = process.communicate()

		if process.returncode != 0:
			raise ProteinDesignToolException, 'Error using pdb2pqr to clean pdb file %s' % inputFile
		
	except BaseException, data:
		print 'Encountered an exception cleaning pdb %s' % inputFile
		if stdout is not None:
			print 'PDB2PQR output follows:'
			print stdout
		
		raise 			

	print 'FINSIHED'

def CleanWHATIF(inputFile, outputFile, removeWater=True, removeLigand=True, removeAltLoc=True, addHydrogens=True, correct=True):

	'''Cleans a PDB using WHAT-IF		

	Errors: 
		Raises an exception if the inputFile is not a valid PDB file.'''

	import pKarun.WI_tools, Protool
	
	#First try to load the pdb using protool.
	#If its not valid 
	try:
		pdb = Protool.structureIO()
		pdb.readpdb(inputFile)
	except Exception, data:	
		raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (inputFile, data)

	try:
		logfile = None
		
		#Create the WHATIF command script
		command = WhatIfCommandBoilerPlate(inputFile)
		
		if removeWater is True:
			command = command + ' %DELWAT \n'
		if removeLigand is True:
			command = command + ' %DELLIG \n'
		if addHydrogens is True:
			command = command + ' %delhyd tot 0\n'
		if correct is True:
			command = command + ' %corall \n N\n'	
		if addHydrogens is True:
			command = command + ' setwif 339 1 \n %addhyd tot 0 \n tot 0\n'
			
		makeMoleculeCommand, whatIfOutputFile = WhatIfMakeMoleculeCommand(inputFile)
		command = command + makeMoleculeCommand
		
		#RunWhatIF - Different command needed if hydrogens are to be added
		if addHydrogens:
			logfile, newPDB = pKarun.WI_tools.RunWI(command,[whatIfOutputFile],None,[inputFile,'HTOP'])
		else:	
			logfile, newPDB = pKarun.WI_tools.RunWI(command,[whatIfOutputFile],None,[inputFile])

		#logfile is the output. newPDB is a dictionary with one entry.
		#The key is the new pdb name (relative and useless).
		#The value are the lines of the new pdb
		#Write the cleaned pdb to specified output location
		stream = open(outputFile, 'w+')		
		stream.writelines(newPDB.values()[0])
		stream.close()

		#Protool ignores altlocs so we can use it to remove them 
		if removeAltLoc is True:
			pdb = Protool.structureIO()
			pdb.readpdb(outputFile)
			pdb.writepdb(outputFile)
		
	except BaseException, data:
		print 'Encountered an exception cleaning pdb %s' % inputFile
		if logfile is not None:
			for line in logfile:
				print line
		
		raise 			

def CleanPDB(inputFile, outputFile, removeWater=True, removeLigand=True, removeAltLoc=True, addHydrogens=True, correct=True):

	'''Cleans a PDB by adding hydrogens, remove waters, ligands and altlocs and correct bad rotamers (default)
	
	The specific cleaning tasks done can be set using the function keywords arguments
	
	Parameters:
		inputFile - The location of the pdb to clean
		outputFile - The location where the cleaned pdb should be written.
		
	Errors: 
		Raises an exception if the inputFile is not a valid PDB file.'''

	cleanProgram = Environment.appConfiguration.pdbCleanProgram()
	if cleanProgram == "WHATIF":
		CleanWHATIF(inputFile, outputFile, 
			removeWater=removeWater, 
			removeLigand=removeLigand, 
			removeAltLoc=removeAltLoc, 
			addHydrogens=addHydrogens, 
			correct=correct)
	elif cleanProgram == "pdb2pqr": 
		if Environment.appConfiguration.has_option('PDBS', 'forceField'):
			forceField = Environment.appConfiguration.get('PDBS', 'forceField')
		else:
			forceField = 'amber'	

		CleanPDB2PQR(inputFile, outputFile, 
			forceField=forceField,
			removeWater=removeWater, 
			removeLigand=removeLigand, 
			removeAltLoc=removeAltLoc, 
			addHydrogens=addHydrogens, 
			correct=correct)
	else:
		print "Unknown pdb clean program %s specified in configuration - aborting" % cleanProgram

	
class CommandLineParser:

	'''Class for parsing the programs command line and providing access to the data'''

	def __fixPaths__(self):
	
		'''For each path supplied on the command line make it absolute'''
		
		#The options with paths are are configurationFile, ligand, mutationList
		if self.configurationFile() is not None:
			self.opt['--configurationFile'] = os.path.abspath(self.opt['--configurationFile'])

		if self.opt.has_key('--mutationList'):
			self.opt['--mutationList'] = os.path.abspath(self.opt['--mutationList'])
			
		if self.opt.has_key('--mutants'):
			self.opt['--mutants'] = os.path.abspath(self.opt['--mutants'])
			
		if self.binding() is True:
			self.opt['--ligand'] = os.path.abspath(self.opt['--ligand'])
			
		if self.opt.has_key('-p'):	
			self.opt['-p'] = os.path.abspath(self.opt['-p'])
			
		if self.opt.has_key('-o'):	
			self.opt['-o'] = os.path.abspath(self.opt['-o'])
			
		if self.opt.has_key('-w'):	
			self.opt['-w'] = os.path.abspath(self.opt['-w'])							

	def __init__(self):
	
		#Asssign an initial dictionary to opt.
		#This will causes all methods to behave in the same way if an option was not specified
		#or if parseCommandLine() was not called.
		self.opt = {}
	
		#getopt will append "--" to all the optional arguments
		#Therefore to check if they were present we have to add the -- here.
		#Ligand isn't included here because it requires a file
		self.optionalArguments = ["--pKa", "--scan", "--delta-pKa", "--mutation", 
					"--stability", "--mutationList", "--ligand",
					"--mutationQuality", "--modes", "--ionisableGroups=", 
					"--keepHeterogens"]
		self.environment = Environment.Environment()
	
	def parseCommandLine(self):

		'''Parses the command line
		
		Execeptions:
			Raises Exceptions.ArgumentError if there was a problem with the arguments supplied'''
	
		try:
			#getopt does not want the "-" or "--" when being called
			self.opt = getopt.getopt(sys.argv[1:], 
					"vw:p:o:hn:", 
					["pKa", "scan", "delta-pKa", "path", "modes", "stability", 
					"mutationList=", "mutation=", "mutants=", "ligand=", "mutationQuality="
					"ionisableGroups=", "configurationFile=", "create-configuration", "keepHeterogens"])
			self.opt = dict(self.opt[0])
		except getopt.GetoptError, data:
			raise Exceptions.ArgumentError, data

		if not self.opt.has_key("-w") and not self.createConfigurationRequested():
			raise Exceptions.ArgumentError, "Working directory (-w) must be specified"
		
		if not self.opt.has_key("-p") and not self.createConfigurationRequested():	
			raise Exceptions.ArgumentError, "PDB (-p) must be specified"
		
		#Check if any of the optional flags were specified		
		counts = [self.opt.keys().count(option) for option in self.optionalArguments]
		value = reduce(lambda x,y: x+y, counts)
		if value == 0 and  not self.createConfigurationRequested():
			self.environment.output('[PEAT-SA] Performing pka scan, stability analysis, ENM analysis and transition path finding')
			#Add the options arguments
			for option in self.optionalArguments:
				if option != '--mutationList':
					self.opt[option] = ""	

		#Add the output directory information if none was supplied	
		if not self.opt.has_key("-o"):
			self.opt["-o"] = os.getcwd()
			self.environment.output('[PEAT-SA] No output directory specified - Defaulting to %s' % self.opt["-o"])
		else:
			self.environment.output('[PEAT-SA] Output directory is %s' % self.opt["-o"])
	
		#Make any supplied paths absolute
		self.__fixPaths__()
		
	def helpRequested(self):
		
		'''Returns True if help was requested'''
		
		if self.opt.has_key("-h"):
			return True
		
		return False	
		
	def createConfigurationRequested(self):
	
		'''Returns True is the user asked for a configuration file to be made'''
		
		#Check if create-configuration was specfied
		if self.opt.has_key("--create-configuration"):
			return True
		
		return False		

	def outputDirectory(self):
	
		'''Returns the output directory to be used
		
		If none was specified returns the directory the program was started from.
		If the command line hasn't been parsed returns None'''
	
		if self.opt.has_key('-o'):
			return self.opt['-o']
		else:
			return None
		
	def workingDirectory(self):
	
		'''Returns the working directory to be used'''
	
		return self.opt["-w"]	
		
	def scan(self):	

		'''Returns True if an  scan was requested'''
	
		if self.opt.has_key('--scan') or self.opt.has_key('--delta-pKa'):
			return True
			
		return False	
		
	def scanMutation(self):
	
		'''Returns the value of the mutation argument is supplied.
		If not this always defaults to ALA.'''
		
		mutation = None
		if self.opt.has_key('--mutation'):
			mutation = self.opt['--mutation']
		else:
			mutation = 'ALA'
	
		return mutation
		
	def stability(self):	

		'''Returns True is a stability calculation was requested'''

		if self.opt.has_key('--stability'):
			return True
			
		return False	
		
	def binding(self):	

		'''Returns True is a binding calculation was requested'''

		if self.opt.has_key('--ligand'):
			return True
			
		return False	
		
	def modes(self):	

		if self.opt.has_key('--modes'):
			return True

		return False
		
	def transitionPath(self):
		
		return False
		
	def verbose(self):	

		'''Returns True if verbose output was requested'''

		if self.opt.has_key('-v'):
			return True
			
		return False

	def configurationFile(self):
	
		'''Returns the path to the specified configuration file or None if none was given'''
	
		#See was a specifc configuration file specified
		try:
			configurationFile = self.opt['--configurationFile']
		except KeyError:
			configurationFile = None
			
		return configurationFile	
		
	def pdbFile(self):
	
		'''Returns the pdb file the program is to be run on
		
		Raises a KeyError if not specified'''
	
		return self.opt["-p"]
		
	def mutationListFile(self):
	
		'''Returns the file containing the mutation list'''
		
		file = None
		if self.opt.has_key('--mutationList'):
			file = self.opt['--mutationList']
					
		return file	
					
	def mutants(self):
	
		'''Returns a mutant collection instance initialised with the contents of the mutants argument'''
	
		mutantCollection = None
		if self.opt.has_key('--mutants'):
			mutantDir = self.opt['--mutants']
			#remove the trailing slash if present so split will work
			if mutantDir[-1:] == os.sep:
				mutantDir = mutantDir[:-1]
				
			pathComponents = os.path.split(mutantDir)
			mutantCollection = Data.MutantCollection(name=pathComponents[1], location=pathComponents[0])
		
		return mutantCollection
			
	def pKa(self):
	
		if self.opt.has_key('--pKa'):
			return True
			
		return False
		
	def ligandFile(self):
	
		'''Returns the value of the ligand command line argument if it exists.'''
	
		if self.opt.has_key('--ligand'):
			return self.opt['--ligand']
			
		return None
		
	def removeHeterogens(self):
	
		'''If keepHeterogens command line argument was specified this returns False.
		
		Otherwise it returns True'''
		
		if self.opt.has_key('--keepHeterogens'):
			return False
			
		return True
		
	def outputName(self):
	
		'''Returns the values of the output name command line argument if it exists'''
	
		if self.opt.has_key('-n'):
			return 	self.opt['-n']
			
		return None
		
	def ionisableGroups(self):
	
		'''Returns the ionisable groups specified on the command line'''
		
		if self.opt.has_key('--ionisableGroups'):
			return self.opt['--ionisableGroups'].split(',')
			
		return None		
								
	def mutationQuality(self):

		'''Defines the threshold modelled mutations must be below to be used.
		
		Returns None if this was not specified
		
		Note: This is the limit on each subsitution'''
		
		if self.opt.has_key('--mutationQuality'):
			return float(self.opt['--mutationQuality'])
			
		return None		
																																			
																													
	
