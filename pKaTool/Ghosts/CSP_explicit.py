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




import sys, os
import string 
import subprocess
import optparse
from utilities_CSP import *



class CSP_coulomb():
	'''
	class CS_coulomb predicts CSPs caused by electric field using simple Coulomb's law
	'''
	
	def __init__(self, pqrNeu, pqrCha, nucleus='N', diel=8.0):
		
		self.pqrNeu = pqrNeu
		self.pqrCha = pqrCha
		self.nucleus = nucleus
		self.diel = diel
		
		
	def parseCoordCharges(self, pqrFile, pdb=False):
		'''
		Parse coordinates and charges from pqr files
		'''
		try:
			fd = open(pqrFile, 'r');
			lines = fd.readlines()
			fd.close()
		except IOError:
			print '\n\t\tFATAL: Cannot open %s file'%pqrFile
			
		coordMol = {}
		chargeMol = {}
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				chainID = parseChainID(line)
				resNum = parseResNum(line)
				atomName = parseAtomName(line)
				resName = parseResName(line)
				if pdb == False:
				    charge = parsePQRCharge(line)
				if not coordMol.has_key(chainID):
					coordMol[chainID] = {}
					chargeMol[chainID] = {}
				if not coordMol[chainID].has_key(resNum):
					coordMol[chainID][resNum] = {}
					chargeMol[chainID][resNum] = {}
				coordinates = eval(getCoordinates(line))
				coordMol[chainID][resNum][atomName] = coordinates
				if pdb == False:
				    chargeMol[chainID][resNum][atomName] = charge
				
		return coordMol, chargeMol
	
	
	def getCSP_coulomb(self):
		'''
		Calculate potentials of all N, C (N and H) atoms
		'''
		coordCha, chargeCha = self.parseCoordCharges(self.pqrCha)
		coordNeu, chargeNeu = self.parseCoordCharges(self.pqrNeu)
		
		diel = self.diel
		deltaCSP = {}
		gradEF = {}
		chainIDs = coordCha.keys()
		chainIDs.sort()
		for chainIdCSP in chainIDs:
			deltaCSP[chainIdCSP] = {}
			gradEF[chainIdCSP] = {}
			resNumsCSP = coordCha[chainIdCSP].keys()
			resNumsCSP.sort()
			try:
				resNumsCSP.remove(1)
			except:
				pass
			#calculate CSP of each backbone nuclei
			for resNumCSP in resNumsCSP:
				E = 0.0
				for chainId in chainIDs:
					resNums = coordCha[chainId].keys()
					for resNum in resNums:
						for atomName in coordCha[chainId][resNum].keys():
							if atomName in backboneAtoms:
								continue
							if self.nucleus == 'N':
								try:
									coordCha[chainIdCSP][resNumCSP]['N']
									distChargeCha, cosAngleCha = getDistanceAngle(coordCha[chainIdCSP][resNumCSP]['N'],
									    coordCha[chainIdCSP][resNumCSP-1]['C'], coordCha[chainId][resNum][atomName])
								except:
									#print '\n\tWarning CSP_coulomb: Atom %d:%s does not exist in %s file'%(resNumCSP, self.nucleus, self.pqrCha)
									continue
								
								
							if self.nucleus == 'H':
								try:
									coordCha[chainIdCSP][resNumCSP]['H']
								except:
									#print '\n\tWarning CSP_coulomb: Atom %d:%s does not exist in %s file'%(resNumCSP, self.nucleus, self.pqrCha)#PROLINE
									continue
								distChargeCha, cosAngleCha = getDistanceAngle(coordCha[chainIdCSP][resNumCSP]['H'], coordCha[chainIdCSP][resNumCSP]['N'], coordCha[chainId][resNum][atomName])
								
							charge = chargeCha[chainId][resNum][atomName]
							E += fieldCoulomb(distChargeCha, cosAngleCha, charge)
							
						for atomName in coordNeu[chainId][resNum].keys():
							if atomName in backboneAtoms:
								continue
							if self.nucleus == 'N':
								try:
									coordNeu[chainIdCSP][resNumCSP]['N']
									distChargeNeu, cosAngleNeu = getDistanceAngle(coordNeu[chainIdCSP][resNumCSP]['N'], 
									    coordNeu[chainIdCSP][resNumCSP-1]['C'], coordNeu[chainId][resNum][atomName])
								except:
									#print '\n\tWarning CSP_coulomb: Atom %d:%s does not exist in %s file'%(resNumCSP, self.nucleus, self.pqrNeu)
									continue								
								
							if self.nucleus == 'H':
								try:
									coordNeu[chainIdCSP][resNumCSP]['H']
								except:
									#print '\n\tWarning CSP_coulomb: Atom %d:%s does not exist in %s file'%(resNumCSP, self.nucleus, self.pqrNeu)#PROLINE
									continue
								distChargeNeu, cosAngleNeu = getDistanceAngle(coordNeu[chainIdCSP][resNumCSP]['H'], coordNeu[chainIdCSP][resNumCSP]['N'], coordNeu[chainId][resNum][atomName])
								
							charge = chargeNeu[chainId][resNum][atomName]
							E -= fieldCoulomb(distChargeNeu, cosAngleNeu, charge)
							
				deltaCSP[chainIdCSP][resNumCSP] = E * nsp[self.nucleus] / diel
		return deltaCSP


def startGhost():
	'''
	Parsing parameters
	'''

	parser=optparse.OptionParser()
	
	#
	#set optparse options
	#
	parser.add_option(
		'--pqrNeu',
		dest='pqrNeu',
		default=None,
		type='string',
		help='<name of the pqr file with neutral moieties>',
		)
	parser.add_option(
		'--pqrCha',
		dest='pqrCha',
		default=None,
		type='string',
		help='<name of the pqr file with charged moieties>',
		)
	parser.add_option(
		'--nucleus',
		dest='nucleus',
		default='N',
		type='choice',
		choices=('N','H',),
		help='<nucleus type (N, H,)>'
		)
	parser.add_option(
		'--diel',
		dest='diel',
		type='float',
		default=8.0,
		help='<protein dielectric constant>'
		)
		
	(options, args,) = parser.parse_args()
	
	##
	##parse optparse options
	##
	pqrNeu=options.pqrNeu
	pqrCha=options.pqrCha
	nucleus=options.nucleus
	diel = options.diel
	
	return pqrNeu, pqrCha, nucleus, diel
	
	
def main():
	
	pqrNeu, pqrCha, nucleus, diel = startGhost()
	Y = CSP_coulomb(pqrNeu, pqrCha, nucleus, diel)
	deltaCSP = Y.getCSP_coulomb()
	
if __name__=='__main__':
	main()
		
