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



class CSP_apbs():
	'''
	class CS_apbs predicts CSPs caused by electric field using APBS solver
	'''
	
	
	def __init__(self, pqrNeu, pqrCha, nucleus ='N', sdiel=80, pdiel=4, ratio=3.0, center="[0.0, 0.0, 0.0]",
		maps=0, radius="[0.0]", active_diel1=80.0, active_diel2=80.0, active_delta=1.0, active_coord="[0.0]",
		gauss=0, mu=977, sigma=97):
		'''
		Initialize the class
		'''
		
		self.pqrNeu = pqrNeu
		self.pqrCha = pqrCha
		self.nucleus = nucleus
		self.sdiel = sdiel
		self.pdiel = pdiel
		self.ratio = ratio
		self.center = eval(center)
		self.maps = maps
		self.radius = eval(radius)
		self.active_diel1 = active_diel1
		self.active_diel2 = active_diel2
		self.active_delta = active_delta
		self.active_coord = eval(active_coord)
		self.gauss = gauss
		self.mu = mu
		self.sigma = sigma
		
		if center == "[0.0, 0.0, 0.0]":
			print '\n\tWarning: Center of focussing has not been provided'
			print '\t         This can influence your results'
	
	
	def parseCoord(self, pdb):
		'''
		Parse coordinates from the file
		'''
		try:
			fd = open(pdb, 'r');
			lines = fd.readlines()
			fd.close()
		except IOError:
			print '\n\t\tFATAL: Cannot open %s file'%pdb
			
		count = 0
		coordMol = {}
		index = {}
		fd_out = open("%s_apbs.pqr"%pdb[:-4], 'w')
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				chainID = parseChainID(line)
				resNum = parseResNum(line)
				atomName = parseAtomName(line)
				resName = parseResName(line)
				if not coordMol.has_key(chainID):
					coordMol[chainID] = {}
					index[chainID] = {}
				if not coordMol[chainID].has_key(resNum):
					coordMol[chainID][resNum] = {}
					index[chainID][resNum] = {}
				coordinates = eval(getCoordinates(line))
				coordMol[chainID][resNum][atomName] = coordinates
				fd_out.write(line)
				index[chainID][resNum][atomName] = count
				count+=1
				if resNum == 1: continue
				if atomName == self.nucleus:
					if self.nucleus == 'N':
						try:
							diff, dist = getDistance(coordMol[chainID][resNum][atomName], coordMol[chainID][resNum-1]['C'])
						except:
							continue
						
					if self.nucleus == 'H':
						try:
							diff, dist = getDistance(coordMol[chainID][resNum][atomName], coordMol[chainID][resNum]['N'])
						except:
							continue
					#
					#coordinates of the dummy atom used for calculations of (gradient of) EF 
					#
					x = coordinates[0] + diff[0]/self.ratio
					y = coordinates[1] + diff[1]/self.ratio
					z = coordinates[2] + diff[2]/self.ratio
					xGrad = x+diff[0]/self.ratio		#used for the gradient of the field
					yGrad = y+diff[1]/self.ratio
					zGrad = z+diff[2]/self.ratio
					coordMol[chainID][resNum]['X'] = [x, y, z]
					line = "%s%s%s%8.3f%8.3f%8.3f  0.0000 0.0000%s"%(line[:13],'X',line[14:30],x,y,z,line[69:])
					fd_out.write(line)
					count+=1
					coordMol[chainID][resNum]['Y'] = [xGrad, yGrad, zGrad]
					line = "%s%s%s%8.3f%8.3f%8.3f  0.0000 0.0000%s"%(line[:13],'Y',line[14:30],xGrad,yGrad,zGrad,line[69:])
					fd_out.write(line)
					count+=1
		
		fd_out.close()
		return coordMol, index
	
	
	def calculatePotential(self, active_diel):
		'''
		Calculate potentials of all N, C (N and H) atoms
		'''
		coordCha, indexCha = self.parseCoord(self.pqrCha)
		coordNeu, indexNeu = self.parseCoord(self.pqrNeu)
		pdiel = self.pdiel
		sdiel = self.sdiel
		deltaCSP = {}
		gradEF = {}
		flag = 0 #0:neutral structure; 1:charged structure
		pqrs = ["%s_apbs.pqr"%self.pqrCha[:-4], "%s_apbs.pqr"%self.pqrNeu[:-4]]
		for pqr in pqrs:
			### writing down the dielectric maps of x, y and z ###
			if self.maps == 1:
				fp_maps = open("apbs_maps.in", "w")
				text = "read\n"
				text += "    mol pqr %s\n"%pqr
				text += "end\n"
				for i in range(4):
					dim = 4.0/pow(2,i)
					text += "elec\n"
					text += "    mg-manual\n"
					text += "    dime 65 65 65\n"
					text += "    grid %.2f %.2f %.2f\n"%(dim, dim, dim)
					text += "    gcent %.3f %.3f %.3f\n"%(self.center[0], self.center[1], self.center[2])
					text += "    mol 1\n"
					text += "    lpbe\n"
					text += "    bcfl sdh\n"
					text += "    ion charge 1 conc 0.150 radius 2.0\n"
					text += "    ion charge -1 conc 0.150 radius 2.0\n"
					text += "    pdie %5.2f\n"%pdiel
					text += "    sdie %5.2f\n"%sdiel
					text += "    srfm mol\n"
					text += "    chgm spl0\n"
					text += "    srad 1.4\n"
					text += "    swin 0.3\n"
					text += "    sdens 10.0\n"
					text += "    temp 298.15\n"
					text += "    write pot dx pot\n"
					text += "    write smol dx acc\n"
					text += "    write dielx dx xdiel%d\n"%i
					text += "    write diely dx ydiel%d\n"%i
					text += "    write dielz dx zdiel%d\n"%i
					text += "    write kappa dx kappa%d\n"%i
					text += "end\n"
				fp_maps.write(text)
				fp_maps.close()
				# run APBS #
				outputStream = open('/dev/null', 'w+')
				process = subprocess.Popen(['apbs apbs_maps.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
				process.wait()
				outputStream.close()
				#
				#changing a dielectric constant of the active site
				#
				(numberOfPoints, residum) = divmod(len(self.active_coord), 3)
				if numberOfPoints > 0:
					for i in range(4):
						for j in range(numberOfPoints):
							if len(self.radius) < j:
								self.radius.append(0.0)
							os.system('~/python/pKaTool/Ghost/3dmaps/diel_active --format=dx --input=xdiel%d.dx --output=xdiel%d.dx --diel=%.2f --xcoor=%.1f --ycoor=%.2f --zcoor=%.2f --radius=%.2f'%(i,i,active_diel, self.active_coord[j*3], self.active_coord[j*3+1], self.active_coord[j*3+2],self.radius[j]))
							os.system('~/python/pKaTool/Ghost/3dmaps/diel_active --format=dx --input=ydiel%d.dx --output=ydiel%d.dx --diel=%.2f --xcoor=%.1f --ycoor=%.2f --zcoor=%.2f --radius=%.2f'%(i,i,active_diel, self.active_coord[j*3], self.active_coord[j*3+1], self.active_coord[j*3+2],self.radius[j]))
							os.system('~/python/pKaTool/Ghost/3dmaps/diel_active --format=dx --input=zdiel%d.dx --output=zdiel%d.dx --diel=%.2f --xcoor=%.1f --ycoor=%.2f --zcoor=%.2f --radius=%.2f'%(i,i,active_diel, self.active_coord[j*3], self.active_coord[j*3+1], self.active_coord[j*3+2],self.radius[j]))
				
				
				bcfl = ["sdh", "map", "map", "map"]
				for i in range(4):
					fp = open("apbs.in", "w")
					text = "read\n"
					text += "    mol pqr %s\n"%pqr
					text += "    diel dx xdiel%d.dx ydiel%d.dx zdiel%d.dx\n"%(i, i, i)
					text += "    kappa dx kappa%d.dx\n"%i
					if i>0:
						text += "    pot dx pot%d.dx\n"%(i-1)
					text += "end\n"
					
					dim = 4.0/pow(2,i)
					text += "elec\n"
					text += "    mg-manual\n"
					text += "    dime 65 65 65\n"
					text += "    grid %.2f %.2f %.2f\n"%(dim, dim, dim)
					text += "    gcent %.3f %.3f %.3f\n"%(self.center[0], self.center[1], self.center[2])
					text += "    mol 1\n"
					text += "    lpbe\n"
					text += "    ion charge 1 conc 0.150 radius 2.0\n"
					text += "    ion charge -1 conc 0.150 radius 2.0\n"
					text += "    srfm mol\n"
					text += "    chgm spl0\n"
					text += "    srad 1.4\n"
					text += "    swin 0.3\n"
					text += "    sdens 10.0\n"
					text += "    temp 298.15\n"
					text += "    write pot dx pot%d\n"%i
					text += "    usemap diel 1-3\n"
					text += "    usemap kappa 1\n"
					text += "    pdie %5.2f\n"%pdiel
					text += "    sdie %5.2f\n"%sdiel
					if i>0:
						text += "    usemap pot 1\n"
					text += "    bcfl %s\n"%bcfl[i]
					text += "    calcenergy no\n"
					text += "    calcforce no\n"
					text += "end\n"
					
					fp.write(text)
					fp.close()
					# run APBS #
					outputStream = open('/dev/null', 'w+')
					process = subprocess.Popen(['apbs apbs.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
					process.wait()
					outputStream.close()
					os.system("rm io.mc")
					#os.system("apbs apbs.in > apbs.log")
					
					
			else:
				fp = open("apbs.in", "w")
				text = "read\n"
				text += "    mol pqr %s\n"%pqr
				text += "end\n"
				bcfl = ["sdh", "focus", "focus", "focus"]
				for i in range(4):
					dim = 4.0/pow(2,i)
					text += "elec\n"
					text += "    mg-manual\n"
					text += "    dime 65 65 65\n"
					text += "    grid %.2f %.2f %.2f\n"%(dim, dim, dim)
					text += "    gcent %.3f %.3f %.3f\n"%(self.center[0], self.center[1], self.center[2])
					text += "    mol 1\n"
					text += "    lpbe\n"
					text += "    bcfl %s\n"%bcfl[i]
					text += "    ion charge 1 conc 0.150 radius 2.0\n"
					text += "    ion charge -1 conc 0.150 radius 2.0\n"
					text += "    pdie %5.2f\n"%pdiel
					text += "    sdie %5.2f\n"%sdiel
					text += "    srfm mol\n"
					text += "    chgm spl0\n"
					text += "    srad 1.4\n"
					text += "    swin 0.3\n"
					text += "    sdens 10.0\n"
					text += "    temp 298.15\n"
					text += "    calcenergy no\n"
					text += "    calcforce no\n"
					text += "end\n"
				
				fp.write(text)
				fp.close()
				# run APBS #
				outputStream = open('/dev/null', 'w+')
				process = subprocess.Popen(['apbs apbs.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
				process.wait()
				outputStream.close()
				#os.system("apbs apbs.in > apbs.log")
				
			#
			#colect potentials
			#
			potential = []
			if self.maps == 1:
				try:
					fd = open("calc_0.txt",'r')
				except IOError:
					print '\n\tCannot find a list of potentials: calc_0.txt'
					sys.exit(1)
			else:
				try:
					fd = open("calc_3.txt",'r')
				except IOError:
					print '\n\tCannot find a list of potentials: calc_3.txt'
					sys.exit(1)
					
			lines = fd.readlines()
			fd.close()
			for line in lines:
				try:
					phi = float(line[:-1])
				except:
					phi = 0.0
				potential.append(phi)
				
			if flag == 0:
				index = indexCha
				coordMol = coordCha
			else:
				index = indexNeu
				coordMol = coordNeu
			
			chainIDs = index.keys()
			for chainID in chainIDs:
				if flag == 0:
					deltaCSP[chainID] = {}
					gradEF[chainID] = {}
				resNums = index[chainID].keys()
				resNums.sort()
				try:
					resNums.remove(1)
				except:
					pass
				for resNum in resNums:
					try:
						count = index[chainID][resNum][self.nucleus]
					except:
						print '\n\tWarning CSP_apbs: Atom %d:%s does not exist in %s file'%(resNum, self.nucleus, pqr)
						continue
						
					try:
						diff, dist = getDistance(coordMol[chainID][resNum][self.nucleus], coordMol[chainID][resNum]['X'])
					except:
						print '\n\tWarning CSP_apbs: Cannot find a distance between %s:%d:%s% and %s:%d:X'%(chainID,resNum,self.nucleus,chainID,resNum,)
						continue
					#
					#calculating CSP[ppm/au] and gradEF[MV/cm2]
					#
					if flag == 0:
						deltaCSP[chainID][resNum] = nsp[self.nucleus] * (potential[count+1] - potential[count]) * kT_to_MV / (dist * MV_to_au)
						gradEF[chainID][resNum] = ((potential[count+2] - potential[count+1])-(potential[count+1] - potential[count]))*kT_to_MV*100/(pow(dist,2)*MV_to_au_efg)
					else:
						deltaCSP[chainID][resNum] -= nsp[self.nucleus]*(potential[count+1] - potential[count])*kT_to_MV/(dist*MV_to_au)
						gradEF[chainID][resNum] -= ((potential[count+2] - potential[count+1])-(potential[count+1] - potential[count]))*kT_to_MV*100/(pow(dist,2)*MV_to_au_efg)
			flag += 1
			
		for chainID in chainIDs:
			for count in range(1,max(index[chainID].keys())):
				if not deltaCSP[chainID].has_key(count):
					deltaCSP[chainID][count] = 0.0
		return deltaCSP
		
		
	def getCSP_apbs(self):
		'''
		main driver for CSP calculations
		'''
			
		if self.maps == 1:
			deltaCSP = {}
			active_diel = self.active_diel1
			while (active_diel <= self.active_diel2):
				deltaActiveCSP = self.calculatePotential(active_diel)
				deltaCSP[active_diel] = deltaActiveCSP
				active_diel += self.active_delta	
			
		else:
			deltaCSP = self.calculatePotential(80.0)
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
		'--sdiel',
		dest='sdiel',
		type='float',
		default=80.0,
		help='<solvent dielectric constant>'
		)
	parser.add_option(
		'--pdiel',
		dest='pdiel',
		type='float',
		default=4.0,
		help='<protein dielectric constant>'
		)
	parser.add_option(
		'--ratio',
		dest='ratio',
		type='float',
		default=3.0,
		help='<ratio to get dx>'
		)
	parser.add_option(
		'--center',
		dest='center',
		type='string',
		default="[5.17, 15.73, 26.25]",
		help='<center of the pbe>'
		)
	parser.add_option(
		'--maps',
		dest='maps',
		type='int',
		default=0,
		help='<0:do not use 3d maps, 1:using 3d maps in apbs>'
		)
	parser.add_option(
		'--radius',
		dest='radius',
		type='string',
		default="[]",
		help='<radii for defining the active site>'
		)
	parser.add_option(
		'--active_diel1',
		dest='active_diel1',
		type='float',
		default=80.0,
		help='<starting dielectric constant of the active site>'
		)
	parser.add_option(
		'--active_diel2',
		dest='active_diel2',
		type='float',
		default=80.0,
		help='<ending dielectric constant of the active site>'
		)
	parser.add_option(
		'--active_delta',
		dest='active_delta',
		type='float',
		default=1.0,
		help='<delta dielectric constant of the active site>'
		)
	parser.add_option(
		'--active_coord',
		dest='active_coord',
		type='string',
		default="[0.0]",
		help='<points for defining the active site>'
		)
	parser.add_option(
		'--gauss',
		dest='gauss',
		type='int',
		default = 0,
		help='<if you want to change NSP using normal distribution>'
		)
	parser.add_option(
		'--mu',
		dest='mu',
		type='float',
		default = 977.0,
		help='<mean value of Gaussian distribution>'
		)
	parser.add_option(
		'--sigma',
		dest='sigma',
		type='float',
		default = 97.7,
		help='<standard deviation of Gaussian distribution>'
		)
		
	(options, args,) = parser.parse_args()
	
	##
	##parse optparse options
	##
	pqrNeu=options.pqrNeu
	pqrCha=options.pqrCha
	nucleus=options.nucleus
	sdiel = options.sdiel
	pdiel=options.pdiel
	ratio = options.ratio
	center = options.center
	maps = options.maps
	radius = options.radius
	active_diel1 = options.active_diel1
	active_diel2 = options.active_diel2
	active_delta = options.active_delta
	active_coord = options.active_coord
	gauss = options.gauss
	mu = options.mu
	sigma = options.sigma
	
	return pqrNeu, pqrCha, nucleus, sdiel, pdiel, ratio, center, maps, radius, active_diel1, active_diel2, active_delta, active_coord, gauss, mu, sigma
	
	
def main():
	
	pqrNeu, pqrCha, nucleus, sdiel, pdiel, ratio, center, maps, radius, active_diel1, active_diel2, active_delta, active_coord, gauss, mu, sigma = startGhost()
	Y = CSP_apbs(pqrNeu, pqrCha, nucleus, sdiel, pdiel, ratio, center, maps, radius, active_diel1, active_diel2, active_delta, active_coord, gauss, mu, sigma)
	deltaCSP = Y.getCSP_apbs()
	
if __name__=='__main__':
	main()
		
