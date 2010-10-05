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
import optparse
from utilities_CSP import *


class CSP_md():
	'''
	class CSP_md predicts CSPs caused by electric field using MD simulations
	'''
	
	def __init__(self, pdb, top, start, end, step, nucleus, diel, resNums):
		
		self.pdb = pdb
		self.top = top
		self.start = start
		self.end = end
		self.step = step
		self.nucleus = nucleus
		self.diel = diel
		self.resNums = resNums
		
		
	def parseCharges(self):
		'''
		Parse charges from the topology file
		'''
		try:
			fd = open(self.top, 'r');
			lines = fd.readlines()
			fd.close()
		except IOError:
			print '\n\t\tFATAL: Cannot open %s file'%self.top
			sys.exit(1)
		chargeDict = {}
		for lineNum in range(len(lines)):
			if lines[lineNum].strip() == '[ atoms ]':
				for lineNum1 in range(lineNum+2, len(lines)):
					if lines[lineNum1].strip() == '':
						break
					line = lines[lineNum1].strip()
					resNum = int(line.split()[2])
					atomName = line.split()[4]
					if atomName in conversionvHv3v2:
						atomName = conversionvHv3v2[atomName]
					chargeValue = float(line.split()[6])
					if not resNum in chargeDict.keys():
						chargeDict[resNum] = {}
					chargeDict[resNum][atomName] = chargeValue
		
		return chargeDict
	
	
	def parseCoord(self, pdb):
		'''
		Parse coordinates from the pdb file
		'''
		try:
			fd = open(pdb, 'r');
			lines = fd.readlines()
			fd.close()
		except:
			print '\n\t\tFATAL: Cannot open %s file'%pdb
			sys.exit(1)
			
		coordMol = {}
		coordSol = {}
		coordCl = {}
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				resNum = parseResNum(line)
				resName = parseResName(line)
				atomName = parseAtomName(line)
				coordinates = eval(getCoordinates(line))
				if resName in solvent:
					if not resNum in coordSol.keys():
						coordSol[resNum] = {}
					coordSol[resNum][atomName] = coordinates
				elif resName in chloride:
					coordCl[resNum] = coordinates
				else:
					if not resNum in coordMol.keys():
						coordMol[resNum] = {}
					coordMol[resNum][atomName] = coordinates
		return coordMol, coordSol, coordCl
	
	
	def calculateCSP(self):
		'''
		calculates potentials of all N, C (N and H) atoms from num1 and num2 side chains, protein, water and Cl
		'''
		
		fd_progress = open("progress.log", 'a');
		fileFirst = self.pdb + '.pdb.' + str(self.start)
		coordMol, coordSol, coordCl = self.parseCoord(fileFirst)
		resNumsMol = coordMol.keys()
		resNumsMol.sort()
		resNumsSol = coordSol.keys()
		resNumsSol.sort()
		resNumsCl = coordCl.keys()
		resNumsCl.sort()
		chargeDict = self.parseCharges()
		diel = self.diel
		if len(self.resNums) > 0:
			resCalculate = self.resNums
		else:
			resCalculate = resNumsMol
		for resNumMol in resCalculate:
			if resNumMol == 1: continue 
			if 'H' not in coordMol[resNumMol].keys(): continue #PROLINE
			linesN = []
			linesH = []
			fileNameN = 'energyN%s.txt'%resNumMol
			fileNameH = 'energyH%s.txt'%resNumMol
			if (os.path.isfile(fileNameN)) and (os.path.isfile(fileNameH)): continue
			print >> fd_progress, "\nPredicting ghost titrations for residue %d ...\n"%(resNumMol)
			fdN = open(fileNameN, 'a')
			fdH = open(fileNameH, 'a')
			#N
			ENside = 0
			ENback = 0
			ENsol = 0
			ENcl = 0
			ENtotal =0#does not include Cl
			#H
			EHside = 0
			EHback = 0
			EHsol = 0
			EHcl = 0
			EHtotal =0#does not include Cl
			counter = 1
			structNum = self.start
			while structNum <= self.end:
				fileName = self.pdb + '.pdb.' + str(structNum)
				filePath = os.path.join(os.getcwd(), fileName)
				if not (os.path.isfile(fileName) or os.path.isfile(filePath)):
					print "WARNING: PDB file %s does not exist "%fileName
					structNum += self.step
					continue
				coordMol, coordSol, coordCl = self.parseCoord(fileName)
				print >> fd_progress, "\tParsing coordinates from structure %s..."%(fileName)
				for resNum in resNumsMol:
					for atomName in coordMol[resNum].keys():
						if (atomName in backboneAtoms) and (abs(resNum-resNumMol)<2): continue#too close to backbone atom
						#print resNum, atomName
						#if not atomName in chargeDict[resNum].keys():
						#	print 'Atom %s doesnt exist in topology file'%atomName
						#	continue #GLU charged doesn't have NE2 in top file
						try:
							coordN = coordMol[resNumMol]['N']
							coordC = coordMol[resNumMol-1]['C']
							coordH = coordMol[resNumMol]['H']
							distN, cosAngleN = getDistanceAngle(coordN,
								coordC, coordMol[resNum][atomName])
						except:
							continue
						try:
							distH, cosAngleH = getDistanceAngle(coordH,
								coordN, coordMol[resNum][atomName])
						except:
							continue
						charge = chargeDict[resNum][atomName]
						EN = fieldCoulomb(distN, cosAngleN, charge) * nsp['N'] / diel
						EH = fieldCoulomb(distH, cosAngleH, charge) * nsp['H'] / diel
						if atomName in backboneAtoms:
							ENback+=EN
							EHback+=EH
						else:
							ENside+=EN
							EHside+=EH
						ENtotal+=EN
						EHtotal+=EH
				
				for resNum in resNumsSol:
					for atomName in coordSol[resNum].keys():
						try:
							distN, cosAngleN = getDistanceAngle(coordN,
								coordC, coordSol[resNum][atomName])
						except:
							continue
						try:
							distH, cosAngleH = getDistanceAngle(coordH,
								coordN, coordSol[resNum][atomName])
						except:
							continue
						charge = chargesSol[atomName]
						EN = fieldCoulomb(distN, cosAngleN, charge) * nsp['N'] / diel
						EH = fieldCoulomb(distH, cosAngleH, charge) * nsp['H'] / diel
						ENsol+=EN
						EHsol+=EH
						ENtotal+=EN
						EHtotal+=EH
				
				for resNum in resNumsCl:
					try:
						distN, cosAngleN = getDistanceAngle(coordN,
							coordC, coordCl[resNum])
					except:
						continue
					try:
						distH, cosAngleH = getDistanceAngle(coordH,
							coordN, coordCl[resNum])
					except:
						continue
					charge = chargesSol['Cl-']
					ENcl+=fieldCoulomb(distN, cosAngleN, charge) * nsp['N'] / diel
					EHcl+=fieldCoulomb(distH, cosAngleH, charge) * nsp['H'] / diel
								
				lineN = '%5d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n'%(structNum, ENside/counter, ENback/counter, ENsol/counter, ENcl/counter, ENtotal/counter)
				fdN.write(lineN)
				lineH = '%5d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n'%(structNum, EHside/counter, EHback/counter, EHsol/counter, EHcl/counter, EHtotal/counter)
				fdH.write(lineH)
				structNum += self.step
				counter+=1
			fdN.close()
			fdH.close()
			print >> fd_progress, "\nDone...."
		fd_progress.close()
		return
	

	def plot(self):
		'''
		Plot the outcome of the MD simulations
		'''
		import matplotlib
		matplotlib.use('Agg')
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import numpy as np
		
		fileFirst = self.pdb + '.pdb.' + str(self.start)
		coordMol, coordSol, coordCl = self.parseCoord(fileFirst)
		resNumsMol = coordMol.keys()
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name = os.path.join(os.getcwd(),'figures/')
		name = ['side chain CSP', 'backbone CSP', 'solvent CSP', 'total CSP']
		color = ['r', 'b', 'c', 'k']
		for nucleus in ['N', 'H']:
			for resNumMol in resNumsMol:
				if resNumMol == 1: continue
				#PROLINE
				CSPside = []
				CSPback = []
				CSPsol = []
				CSPtotal = []
				if 'H' not in coordMol[resNumMol].keys(): continue #PROLINE
				fileName = 'energy%s%s.txt'%(nucleus, resNumMol)
				try:
					fd = open(fileName, "r")
					lines = fd.readlines()
					fd.close()
				except IOError:
					continue
				for line in lines:
					CSPside.append(float(line.split()[1]))
					CSPback.append(float(line.split()[2]))
					CSPsol.append(float(line.split()[3]))
					CSPtotal.append(float(line.split()[5]))
					
				numberFrames = len(lines)
				frames = range(1, numberFrames+1)
				plt.clf()
				plt.subplot(111)
				plt.grid(False)
				plt.hold(True)
				plt.title('CSP of %s%d'%(nucleus,resNumMol), fontsize='14')
				plt.xlabel('MD structure', verticalalignment = 'top', fontsize='12')
				plt.ylabel('predicted CSP[ppm/au]',fontsize='12')
				plt.plot(frames,CSPside,linewidth=2.0, label=name[0])
				plt.plot(frames,CSPback,linewidth=2.0, label=name[1])
				plt.plot(frames,CSPsol,linewidth=2.0, label=name[2])
				plt.plot(frames,CSPtotal,linewidth=2.0, label=name[3])
				plt.legend(shadow=True, fancybox=True, loc = 'best')
				leg = plt.gca().get_legend()
				ltext  = leg.get_texts()  # all the text.Text instance in the legend
				llines = leg.get_lines()  # all the lines.Line2D instance in the legend
				frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend
				frame.set_facecolor('0.95') # set the frame face color to light gray
				frame.set_alpha(0.9) # set the alpha value of the legend: it will be translucent
				plt.setp(ltext, fontsize='10')    # the legend text fontsize
				plt.setp(llines, linewidth=1.0)      # the legend linewidth
				figureName = dir_name+'nuclei%s%d.png'%(nucleus,resNumMol)
				plt.savefig(figureName, dpi = 600)
	
	
	
def parseGhost():
	'''
	Function for starting Ghost calculations
	'''
	
	print 
	print 'Ghost titration predictions using MD simulations'
	print

	parser=optparse.OptionParser()
	
	#
	#set optparse options
	#
	parser.add_option(
		'--pdb',
		dest='pdb',
		default=None,
		type='string',
		help='<name of the pdb file together with the path (this name is without number in the MD array and .pdb)>',
		)
	parser.add_option(
		'--top',
		dest='top',
		default=None,
		type='string',
		help='<name of the topology file from MD simulations>',
		)
	parser.add_option(
		'--start',
		dest='start',
		default='0',
		type='string',
		help='<number of the first structure>',
		)
	parser.add_option(
		'--end',
		dest='end',
		default='20000',
		type='string',
		help='<number of the last structure>',
		)
	parser.add_option(
		'--step',
		dest='step',
		default='1',
		type='string',
		help='<increment when using available structures>',
		)
	parser.add_option(
		'--nucleus',
		dest='nucleus',
		default='N',
		type='choice',
		choices=('N','H','C',),
		help='<nucleus type (N, H, C)>'
		)
	parser.add_option(
		'--exper_CS',
		dest='exper_CS',
		type='string',
		default="{29:-0.28, 30:-0.16, 32:-0.51, 33:-0.53, 34:0.46, 38:0.42, 42:0.2, 44:0.25, 45:-0.1, 46:0.23, 51:0.15, 52:0.3, 54:-0.15, 55:-0.13, 56:0.24, 58:-0.7, 59:-0.14, 60:-0.1, 91:0.2, 107:0.1, 108:0.76, 108:0.61, 111:0.81, 113:0.56, 114:0.31, 115:0.27}",
		help='<dictionary of experimental pka values {tit.group:pka,...}>'
		)
	parser.add_option(
		'--diel',
		dest='diel',
		type='float',
		default=2.0,
		help='<dielectric constant>'
		)
	parser.add_option(
		'--resNums',
		dest='resNums',
		type='string',
		default="[]",
		help='<nuclei for which we can calculate CSP "[9, 12]">'
		)
	
	(options, args,) = parser.parse_args()
	
	##
	##parse optparse options
	##
	pdb=options.pdb
	top = options.top
	start=int(options.start)
	end=int(options.end)
	step=int(options.step)
	nucleus=options.nucleus
	exper_CS=eval(options.exper_CS)
	diel=options.diel
	resNums = eval(options.resNums)
	
	
	return pdb, top, start, end, step, nucleus, diel, exper_CS, resNums
	
	
	
def main():
	
	pdb, top, start, end, step, nucleus, diel, exper_CS, resNums = parseGhost()
	Y = CSP_md(pdb, top, start, end, step, nucleus, diel, resNums)
	Y.calculateCSP()
	Y.plot()
	
if __name__=='__main__':
	main()
		
