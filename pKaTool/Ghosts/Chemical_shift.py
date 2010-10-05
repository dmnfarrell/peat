#!/usr/bin/env python
#
# pKaTool - scripts for analysing chemical shift perturbations
# Copyright (C) 2010 Predrag Kukic & Jens Erik Nielsen
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
# Email: Jens.Nielsen@gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#



cha2neuAcid = {'ASP':'ASH', 'GLU':'GLH',
	'CYM':'CYS', 'HIS':'HIP'} #HIS is HID in AMBER
	
neu2chaAcid = {'ASH':'ASP', 'GLH':'GLU',
	'CYS':'CYM', 'HIP':'HIS'}
	
cha2neuBase = {'LYS':'LYN'}
	
import sys, os
import string 
import subprocess
import optparse
from utilities_CSP import *


class CSP():
	'''
	Class CSP predicts CSPs (ghost titrations) caused by electric field effects
	using a PDB structure and charges assigned by PDB2PQR. It returns 
	dict1 = {'chainID:resNum':{titGroup1:span1, ..., titGroupN:spanN}, ...}
	dict2 = {'chainID:resNum':{pH1:CSP1, ..., pHN:CSPN}, ...}
	'''
	
	
	def __init__(self, pdb = None, nucleus = 'N', method = 'Coulomb', pKaDict = None,
					protDiel = 4.0, solvDiel = 80.0,):
		'''
		Initialize the class
		'''
		self.pdb = pdb
		self.nucleus = nucleus
		self.method = method
		self.pKaDict = pKaDict
		self.solvDiel = solvDiel
		self.protDiel = protDiel
	
	
	def parseTitGroups(self):
		'''
		Get the list of all titrating groups by comparing pKaDict and
		the groups parsed from PDB file
		'''
		import pKaTool.Ghost.get_dEF as get_dEF
		X = get_dEF.map_struct(self.pdb)
		titGroupsPDB = X.PI.get_titratable_groups()
		for titGroup in titGroupsPDB:
			if (not self.pKaDict.has_key(titGroup)) and (getResName(titGroup) in ['ASP', 'GLU', 'HIS']):
				pka = modelPKA[getResName(titGroup)]
				print '\n\tWarning: pKa value of %s has not been provided'%titGroup
				print '\t\t The model pka value of %.2f will be used'%pka
				self.pKaDict[titGroup] = pka
		
		#sorting pKaDict by pKa values in decending order 
		items = [(pKa, group) for group, pKa in self.pKaDict.items()]
		items.sort()
		pKaSort = [(group, pKa) for pKa, group in items]
		self.pKaSort = pKaSort
		return
	
	
	def whatIF(self):
		'''
		Run what if on the original pdb file
		'''
		try:
			import WHAT_IF
		except ImportError:
			import pKarun.WHAT_IF as WHAT_IF
		
		self.dirName = os.getcwd()
		os.chdir(self.tmpDirName)
		
		X = WHAT_IF.Whatif()
		X.RedirectOutput('whatif.log')
		command = 'GETMOL %s\n \n %%DELWAT\n %%DELLIG\n %%CORALL\n N\n %%MAKMOL\n'%self.pdb
		command += '\n %s\n Y\n TOT 0\n \n STOP\n Y\n'%self.pdbPath
		X.Run(command, 'WI.log')
		#
		#change all GLU and ASP to be in the neutral form
		#
		try:
			fd = open(self.pdbPath, 'r')
			lines = fd.readlines()
			fd.close()
		except IOError:
			print '\n\tFATAL: Cannot open %s file'%pdbCha 
			sys.exit(1)
		linesOutput = []
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				atomName = parseAtomName(line)
				#convert H names
				if conversionH.has_key(atomName):
					line = line[:12] + conversionH[atomName] + line[16:]
				group = getTitGroupID(line)
				#change resName for the provided titGroupSpan
				name = getResName(group)
				if cha2neuAcid.has_key(name):
					line = line[:17] + cha2neuAcid[name] + line[20:]
			linesOutput += [line]
		
		fd = open(self.pdbPath ,'w')
		fd.writelines(linesOutput)
		fd.close()
		return
	
	
	def getSpanDict(self):
		'''
		main driver that returns the dictionary of CSP spans
		'''
		self.parseTitGroups()
		self.workingDir = os.getcwd()
		pdb = os.path.split(self.pdb)[1]
		self.tmpDirName = '/tmp/%s'%pdb[:-4]
		if not os.path.isdir(self.tmpDirName):
			os.system('mkdir %s'%self.tmpDirName)
		os.system('cp %s %s'%(self.pdb,self.tmpDirName))
		self.pdbPath = os.path.join(self.tmpDirName, pdb)
		
		print '\n\tCalculating CSPs using %s method...'%self.method
		self.whatIF()
		pKaSort = self.pKaSort
		dictSpan = {}
		count = 0
		try:
			while 1:
				if pKaSort[count][1] > (pKaSort[count+1][1] - 0.1):
					print '\n\t\tWarning: Groups %s and %s have very close pKa values'%(pKaSort[count][0], pKaSort[count+1][0])
					print '\t\t	 and it is difficult to differentiate the CSPs they create'
					count += 1
				else:
					break
		except IndexError:
			pass
			
		dictCSPs = {}
		count = 0
		while (count < len(pKaSort)):
			titGroupSpan = [] # titGroup(s) which span will be calculated
			titGroupSpan.append(pKaSort[count])
			print '\n\t\tCalculating spans coming from %s...\n'%titGroupSpan[0][0]
			deltaCSP = self.calculateSpan(titGroupSpan)
			
			for chainID in deltaCSP.keys():
				for resNum in deltaCSP[chainID].keys():
					resID = '%s:%s'%(chainID.strip(), string.zfill(str(resNum), 4))
					if not dictCSPs.has_key(resID):
						dictCSPs[resID] = {}
					name = titGroupSpan[0][0]					
					dictCSPs[resID][name] = deltaCSP[chainID][resNum]
			count += 1
		try:
			os.chdir(self.dirName)
			os.system('rm -r %s'%self.tmpDirName)
		except:
			print '\n\tWarning: Cannot delete tmp directory %s'%self.tmpDirName
		return dictCSPs
	
	
	def calculateSpan(self, titGroupSpan):
		'''
		Calculates the CSP span coming from titGroups
		give in titGroupSpan.
		CSPs depend on the method used
		'''
		titGroupList = []
		for titGroup in titGroupSpan:
			titGroupList.append(titGroup[0])
		try:
			fd = open(self.pdbPath, 'r')
			lines = fd.readlines()
			fd.close()
		except IOError:
			print '\n\t\tFATAL: Cannot open %s file'%pdbCha 
			sys.exit(1)
		#
		#Print PDB files on low and high pH values
		#
		lineshighPH = []
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				group = getTitGroupID(line)
				#change resName for the provided titGroupSpan
				if group in titGroupList:
					name = getResName(group)
					coordinates = getCoordinates(line)
					if (cha2neuAcid.has_key(name)) or (neu2chaAcid.has_key(name)):
						linehighPH = [line[:17] + name + line[20:]]
					elif cha2neuBase.has_key(name):
						linehighPH = [line[:17] + cha2neuBase[name] + line[20:]]
					else:
						linehighPH = [line]
				else:
					linehighPH = [line]
				
				lineshighPH += linehighPH
				
		pdb = os.path.split(self.pdb)[1]
		firstTitGroup = getResNum(titGroupList[0])
		pdbPath = self.pdbPath
		pdbPathNew =  pdbPath[:-4] + '_%d.pdb'%firstTitGroup
		fd = open(pdbPathNew ,'w')
		fd.writelines(lineshighPH)
		fd.close()
		#
		# Run PDB2PQR
		#
		pqrPath = '%s.pqr'%pdbPath[:-4]
		pqrPathNew = '%s.pqr'%pdbPathNew[:-4]
		self.callpdb2pqr(pdbPath, pdbPathNew, pqrPath, pqrPathNew)
		#
		# calculate span(CSPs)
		#
		if self.method == 'PBE':
			from CSP_pbe import CSP_apbs
			Y = CSP_apbs(pqrNeu=pqrPath, pqrCha=pqrPathNew, nucleus=self.nucleus, center = coordinates, sdiel=self.solvDiel, pdiel=self.protDiel,)
			deltaCSP = Y.getCSP_apbs()
		else:
			from CSP_explicit import CSP_coulomb
			Y = CSP_coulomb(pqrNeu=pqrPath, pqrCha=pqrPathNew, nucleus=self.nucleus, diel=self.protDiel)
			deltaCSP = Y.getCSP_coulomb()
			
		self.pdbPath = pdbPathNew
		return deltaCSP
		
		
	def callpdb2pqr(self, pdbPath, pdbPathNew, pqrPath, pqrPathNew):
		'''
		Run pdb2pqr and convert pdb into pqr files
		'''
		try:
			from main import mainCommand
		except ImportError:
			from pdb2pqr.main import mainCommand
		#
		#using a system call
		#
		dir_name = os.getcwd()
		source = '%s/pqr.txt'%dir_name
		fd = open(source, 'w')
		fd.writelines([
			'python ~/python/pdb2pqr/pdb2pqr.py --ff=parse --chain %s %s\n'%(pdbPathNew, pqrPathNew),
			'python ~/python/pdb2pqr/pdb2pqr.py --ff=parse --chain %s %s\n'%(pdbPath, pqrPath),
			])
		fd.close()
		os.system('source %s > pqr.out'%source)
		#
		#without a system call
		#
		'''command1 = ['/home/predrag/python/pdb2pqr_svn/pdb2pqr.py', '--ff=amber', '%s'%pdbPathNew, '%s'%pqrPathNew]
		command2 = ['/home/predrag/python/pdb2pqr_svn/pdb2pqr.py', '--ff=amber', '%s'%pdbPath, '%s'%pqrPath]
		mainCommand(command1)
		mainCommand(command2)'''
		#
		#
		#
		return


	def getTitCurves(self, dictPh):
		'''
		#Return CSPs calculated at different pH values provided in dictPh
		'''
		resIDs = dictPh.keys()
		resIDs.sort()
		dictCSPs = self.dictCSPs
		dictTitCurves = {}
		for resID in resIDs:
			if not dictCSPs.has_key(resID):
				print "\n\tWarning: %s does not exist in the %s file"%(resID, self.pdb)
				continue
			if not dictTitCurves.has_key(resID):
				dictTitCurves[resID] = {}
			for pH in dictPh[resID]:
				CSP = 0
				for titGroup in dictCSPs[resID].keys():
					titName = string.split(titGroup, ':')[-1]
					if titName in ['HIS', 'LYS', 'ARG', 'NTERM', 'NTR']: 
						CSP += (dictCSPs[resID][titGroup])/(1+10**(pH - self.pKas[titGroup]))
					else:
						CSP += (dictCSPs[resID][titGroup])/(1+10**(self.pKas[titGroup] - pH))
						
				dictTitCurves[resID][pH] = CSP
		return dictTitCurves
	
	
	
def getTitGroupID(line):
	'''
	Return titGroupID A:resNum:resName
	'''
	chainID = parseChainID(line)
	resName = parseResName(line)
	if neu2chaAcid.has_key(resName):
		resName = neu2chaAcid[resName]
	resNum = line[22:26].strip()
	resNum = string.zfill(resNum, 4)
	return '%s:%s:%s'%(chainID, resNum, resName)
	
	
def parseCSP():
	'''
	Parse the options from the command line
	'''

	parser=optparse.OptionParser()
	
	#
	#set optparse options
	#
	parser.add_option(
		'--pdb',
		dest='pdb',
		default=None,
		type='string',
		help='<name of the pdb file>',
		)
	parser.add_option(
		'--nucleus',
		dest='nucleus',
		default='N',
		type='choice',
		choices=('N','H',),
		help='<nucleus type (N, H)>'
		)
	parser.add_option(
		'--method',
		dest='method',
		default='Coulomb',
		type='choice',
		choices=('Coulomb','PBE',),
		help='<method used for calculating electric field (Coulomb, PBE)>'
		)
	parser.add_option(
		'--pKaDict',
		dest='pKaDict',
		type='string',
		default="{}",
		help='<dictionary of pKa values {titGroup:pKa,...}>'
		)
	parser.add_option(
		'--protDiel',
		dest='protDiel',
		type='float',
		default=4.0,
		help='<protein dielectric constant>'
		)
	parser.add_option(
		'--solvDiel',
		dest='solvDiel',
		type='float',
		default=80.0,
		help='<solvent dielectric constant>'
		)
		
	(options, args,) = parser.parse_args()
	##
	##parse optparse options
	##
	pdb = options.pdb
	nucleus = options.nucleus
	method = options.method
	pKaDict = eval(options.pKaDict)
	protDiel = options.protDiel
	solvDiel = options.solvDiel#
	
	return pdb, nucleus, method, pKaDict, protDiel, solvDiel
	
	
def main():
	#main driver for testing the script
	pdb, nucleus, method, pKaDict, protDiel, solvDiel = parseCSP()
	Y = CSP(pdb, nucleus, method, pKaDict, protDiel, solvDiel)
	dictCSPs = Y.getSpanDict()
	
if __name__=='__main__':
	main()
		
