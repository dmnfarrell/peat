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


qe=1.602E-19
e0=8.85E-12
au_field = 5.14225E+11
kT_to_MV = 2.5692 #from kT/eA to MV/cm
MV_to_au = 5142.25 #from MV/cm to ppm/au
MV_to_au_efg = 971744.7 # from MV/cm2 to ppm/au efg

nsp = {'N':977.0, 'H':90}

conversionHv2v3 = {
	## v2 keys
	## v3 values
	'1HH1':'HH11', ## ARG
	'1HH2':'HH12', ## ARG
	'2HH1':'HH21', ## ARG
	'2HH2':'HH22', ## ARG
	'1HE2':'HE21', ## GLN
	'2HE2':'HE22', ## GLN
	'1HG1':'HG11', ## VAL
	'2HG1':'HG12', ## VAL
	'3HG1':'HG13', ## VAL
	'1HG2':'HG21', ## VAL
	'2HG2':'HG22', ## VAL
	'3HG2':'HG23', ## VAL
	'1HD1':'HD11', ## LEU
	'2HD1':'HD12', ## LEU
	'3HD1':'HD13', ## LEU
	'1HD2':'HD21', ## LEU,ASN
	'2HD2':'HD22', ## LEU,ASN
	'3HD2':'HD23', ## LEU
	}
	
conversionvHv3v2 = {
	## v3 keys
	## v2 values
	'HH11':'1HH1', ## ARG
	'HH12':'1HH2', ## ARG
	'HH21':'2HH1', ## ARG
	'HH22':'2HH2', ## ARG
	'HE21':'1HE2', ## GLN
	'HE22':'2HE2', ## GLN
	'HG11':'1HG1', ## VAL
	'HG12':'2HG1', ## VAL
	'HG13':'3HG1', ## VAL
	'HG21':'1HG2', ## VAL
	'HG22':'2HG2', ## VAL
	'HG23':'3HG2', ## VAL
	'HD11':'1HD1', ## LEU
	'HD12':'2HD1', ## LEU
	'HD13':'3HD1', ## LEU
	'HD21':'1HD2', ## LEU,ASN
	'HD22':'2HD2', ## LEU,ASN
	'HD23':'3HD2', ## LEU
	}
	
chargesSol = {
	'OW':-.834, ## TIP3
	'O':-.834,  ## TIP3
	'HW1':.417, ## TIP3
	'H1':.417,  ## TIP3
	'HW2':.417, ## TIP3
	'H2':.417,  ## TIP3
	'CL':-1,
	'Cl':-1,
	'Cl-':-1,
	}
chloride = [' Cl', ' CL', 'Cl-', 'CL-']
solvent = ['SOL', 'WAT']
modelPKA = {'ASP':4.0, 'GLU':4.4, 
	'HIS':6.3, 'CYS':8.7,'TYR':9.6, 
	'LYS':10.4, 'ARG':13.0
	}
backboneAtoms = [
    'N','H','H1','H2','H3','HXT', ## H3 if N-terminal (H,HXT in PDB)
    'CA','HA','HA2','HA3', ## HA2 if Gly (HA3 in PDB)
    'C','O','OC1','OC2','OXT', ## OC1,OC2 if C-terminal (O,OXT in PDB)
    ]
terminalAtoms = [
	'H1','H2','H3','H','HXT', ## N-terminal
	'C','O','OC1','OC2','OXT', ##C-terminal
	'C',"""O''""","""O'""", ##C-terminal
	]

distanceAtoms = {'ASP':'CG', 'GLU':'CD', 
	'HIS':'CE1', 'CYS':'SG','TYR':'CZ', 
	'LYS':'NZ', 'ARG':'CZ'
	}
	
import string


def getTitGroup(line):
	'''
	Return titGroupID A:resNum:resName
	'''
	chainID = line[21]
	resName = line[17:20].strip()
	resNum = line[22:26].strip()
	resNum = string.zfill(resNum, 4)
	return '%s:%s:%s'%(chainID, resNum, resName)
	
	
def getCoordinates(line):
	'''
	Return [x, y, z] coordinates of the atom
	'''
	x = line[30:38].strip()
	y = line[38:46].strip()
	z = line[46:54].strip()
	coordinates = '[%s, %s, %s]'%(x, y, z)
	return coordinates
	
	
def getChainID(titGroup):
	'''
	Return the chainID of the titratable group
	'''
	return string.split(titGroup, ':')[0]


def getResNum(titGroup):
	'''
	Return residue number of the titGrouop
	'''
	return int(string.split(titGroup, ':')[1])


def getResName(titGroup):
	'''
	Return residue name of the titGroup
	'''
	return string.split(titGroup, ':')[2]


def parseResName(line):
	'''
	Parse residue name of the line
	'''
	return line[17:20].strip()


def parseAtomName(line):
	'''
	Parse atom name of the line
	'''
	return line[12:16].strip()


def parseResNum(line):
	'''
	Parse residue number from the line
	'''
	return int(line[22:26].strip())


def parseChainID(line):
	'''
	Parse residue number from the line
	'''
	return line[21]


def parsePQRCharge(line):
	'''
	Parse charges from the line
	'''
	return (float(line[55:62].strip()))


def getDistance(coord1, coord2):
	'''
	Calculate distance between two dots in space
	'''
	import numpy
	import math
	
	xDiff = coord2[0]-coord1[0]
	yDiff = coord2[1]-coord1[1]
	zDiff = coord2[2]-coord1[2]
	Diff = [xDiff, yDiff, zDiff]
	dist = math.sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff)
	
	return Diff, dist


def getDistanceAngle(coord1, coord2, titCoord):
	'''
	Calculates bonds and angles between them
	'''
	import numpy
	bondVector, distBond = getDistance(coord1, coord2)
	chargeVector, distCharge = getDistance(coord1, titCoord)
	dp=numpy.dot(bondVector,chargeVector)
	cosAngle=dp/(distBond*distCharge)
	return distCharge, cosAngle


def fieldCoulomb(distCharge, cosAngle, charge):
	'''
	Calculate electric field using Coulomb's law
	'''
	import math
	E = 1.0/(4*math.pi*e0)*qe/(distCharge*1E-10)**2 *charge*cosAngle / au_field
	return E
	
