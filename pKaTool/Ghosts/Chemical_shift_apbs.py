#!/usr/bin/python
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



MV_to_au = 5142.25
MV_to_au_efg = 971744.7
conversionH = {
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

chargesSol = {
	'OW':-.834, ## TIP3
	'HW1':.417, ## TIP3
	'HW2':.417, ## TIP3
	'CL':-1,
	'Cl':-1,
	}

backboneAtom = [
    'N','H','H1','H2','H3', ## H3 if N-terminal (H,HXT in PDB)
    'CA','HA','HA2', ## HA2 if Gly (HA3 in PDB)
    'C','O','OC1','OC2', ## OC1,OC2 if C-terminal (O,OXT in PDB)
    ]



import sys, os, subprocess
import string
import optparse



class CS_apbs_analysis():
	'''
	class CS_explicit_analysis predicts CSPs caused by electric field using MD simulations
	'''
	
	def __init__(self, pqrNeu, pqrCha, atom, exper_CS, threshold, plot, diel_low, diel_high, diel_increment, sdiel, ratio, center, error_CS, maps, radius,  active_diel1, active_diel2, active_delta, active_coord):
		
		self.pqrNeu = pqrNeu
		self.pqrCha = pqrCha
		self.nucleus = atom
		self.diel_low = diel_low
		self.diel_high = diel_high
		self.diel_increment = diel_increment
		self.sdiel = sdiel
		self.ratio = ratio
		self.center = center
		self.CSP_exp = exper_CS
		self.error_CS = error_CS
		self.maps = maps
		self.radius = radius
		self.active_diel1 = active_diel1
		self.active_diel2 = active_diel2
		self.active_delta = active_delta
		self.active_coord = active_coord
		
	
	def parseCoord(self, pdb):
		'''
		Parse coordinates from the file
		'''
		count = 0
		coordMol = {}
		index = {}
		fd1 = open("%s_apbs.pqr"%pdb[:-4], 'w')
		try:
			fd = open(pdb, 'r');
			lines = fd.readlines()
			fd.close()
		except:
			print 'Cannot open the structure %s'%(pdb)
		for line in lines:
			record = line[:6].strip()
			if record == 'ATOM':
				resNum = int(line[22:26])
				atomName = line[12:16].strip()
				if atomName in conversionH.keys():
					atomName = conversionH[atomName]
					
				resName = line[17:20]
				if not resNum in coordMol.keys():
						coordMol[resNum] = {}
						index[resNum] = {}
				x = float(line[30:38])
				y = float(line[38:46])
				z = float(line[46:54])
				coordMol[resNum][atomName] = [x, y, z]
				
				
				fd1.write(line)
				index[resNum][atomName] = count
				count+=1
				if resNum == 1: continue
				if atomName == self.nucleus:
					if (self.nucleus == 'N'):
						try:
							diff, dist = self.distance(coordMol[resNum][atomName], coordMol[resNum-1]['C'])
						except:
							continue
						
					if (self.nucleus == 'H'):
						try:
							diff, dist = self.distance(coordMol[resNum][atomName], coordMol[resNum]['N'])
						except:
							continue
						
					x += diff[0]/self.ratio
					y += diff[1]/self.ratio
					z += diff[2]/self.ratio
					xGrad = x+diff[0]/self.ratio		#used for the gradient of the field
					yGrad = y+diff[1]/self.ratio
					zGrad = z+diff[2]/self.ratio
					coordMol[resNum]['X'] = [x, y, z]
					line = "%s%s%s%8.3f%8.3f%8.3f  0.0000 0.0000%s"%(line[:13],'X',line[14:30],x,y,z,line[69:])
					fd1.write(line)
					count+=1
					coordMol[resNum]['Y'] = [xGrad, yGrad, zGrad]
					line = "%s%s%s%8.3f%8.3f%8.3f  0.0000 0.0000%s"%(line[:13],'Y',line[14:30],xGrad,yGrad,zGrad,line[69:])
					fd1.write(line)
					count+=1
					
		
		fd1.close()
					
		return coordMol, index
	
	
	def distance(self, coord1, coord2):
		
		import numpy
		import math
		
		xDiff = coord2[0]-coord1[0]
		yDiff = coord2[1]-coord1[1]
		zDiff = coord2[2]-coord1[2]
		Diff = [xDiff, yDiff, zDiff]
		dist = math.sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff)
		
		return Diff, dist
	
	
	
	def calculatePotential(self, active_diel):
		'''
		calculates potentials of all N, C (N and H) atoms
		'''
		
		diel = self.diel_low
		sdiel = self.sdiel
		coordCha, indexCha = self.parseCoord(self.pqrCha)
		coordNeu, indexNeu = self.parseCoord(self.pqrNeu)
		deltaEF = {}
		gradxEF = {}
		difCSP={}
		difCSP_total={}
		difCSP_total_scale={}
		while (diel<=self.diel_high):
			deltaEF[diel] = {}
			gradxEF[diel] = {}
			difCSP[diel]={}
			difCSP_total[diel]=0.0
			difCSP_total_scale[diel]=0.0
			flag = 0 #0 means neutral structure and 1 means charged structure
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
						text += "    ion charge 1 conc 0.050 radius 2.0\n"
						text += "    ion charge -1 conc 0.050 radius 2.0\n"
						text += "    pdie %5.2f\n"%diel
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
					outputStream = open('/dev/null', 'w+')
					process = subprocess.Popen(['/home/predrag/bin/apbs/apbs apbs_maps.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
					process.wait()
					outputStream.close()
					#os.system("/home/predrag/bin/apbs/apbs apbs_maps.in > apbs.log")
					
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
							text += "    pot dx pot.dx\n"
						text += "end\n"
						
						dim = 4.0/pow(2,i)
						text += "elec\n"
						text += "    mg-manual\n"
						text += "    dime 65 65 65\n"
						text += "    grid %.2f %.2f %.2f\n"%(dim, dim, dim)
						text += "    gcent %.3f %.3f %.3f\n"%(self.center[0], self.center[1], self.center[2])
						text += "    mol 1\n"
						text += "    lpbe\n"
						text += "    ion charge 1 conc 0.050 radius 2.0\n"
						text += "    ion charge -1 conc 0.050 radius 2.0\n"
						text += "    srfm mol\n"
						text += "    chgm spl0\n"
						text += "    srad 1.4\n"
						text += "    swin 0.3\n"
						text += "    sdens 10.0\n"
						text += "    temp 298.15\n"
						text += "    write pot dx pot\n"
						text += "    usemap diel 1-3\n"
						text += "    usemap kappa 1\n"
						text += "    pdie %5.2f\n"%diel
						text += "    sdie %5.2f\n"%sdiel
						if i>0:
							text += "    usemap pot 1\n"
						text += "    bcfl %s\n"%bcfl[i]
						text += "    calcenergy no\n"
						text += "    calcforce no\n"
						text += "end\n"
						
						fp.write(text)
						fp.close()
						outputStream = open('/dev/null', 'w+')
						process = subprocess.Popen(['/home/predrag/bin/apbs/apbs apbs.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
						process.wait()
						outputStream.close()
						os.system("rm io.mc")
						
						
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
						text += "    ion charge 1 conc 0.050 radius 2.0\n"
						text += "    ion charge -1 conc 0.050 radius 2.0\n"
						text += "    pdie %5.2f\n"%diel
						if self.maps == 2:
							text += "    sdie %5.2f\n"%diel
						else:
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
					outputStream = open('/dev/null', 'w+')
					process = subprocess.Popen(['/home/predrag/bin/apbs/apbs apbs.in'], stdout=outputStream, stderr=subprocess.STDOUT, shell=True)
					process.wait()
					outputStream.close()
					#os.system("/home/predrag/bin/apbs/apbs apbs.in > apbs.log")
					
				
				#
				#colect results
				#
				potential = []
				if self.maps == 1:
					fd = open("calc_0.txt",'r')
				else:
					fd = open("calc_3.txt",'r')
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
					
				resNums = index.keys()
				try:
					resNums.remove(1)
				except:
					pass
					
				
				for resNum in resNums:
					try:
						count = index[resNum][self.nucleus]
					except:
						#resNums.remove(resNum)
						if resNum in self.CSP_exp.keys():
							del self.CSP_exp[resNum]
						continue
					try:
						diff, dist = self.distance(coordMol[resNum][self.nucleus], coordMol[resNum]['X'])
					except:
						#resNums.remove(resNum)
						if resNum in self.CSP_exp.keys():
							del self.CSP_exp[resNum] 
						continue
					#
					#EF in MV/cm
					#
					if flag == 0:
						deltaEF[diel][resNum] = (potential[count+1] - potential[count])*2.5692/(dist*MV_to_au)
						gradxEF[diel][resNum] = ((potential[count+2] - potential[count+1])-(potential[count+1] - potential[count]))*2.5692*100/(pow(dist,2)*MV_to_au_efg)	#MV/cm2
					else:
						deltaEF[diel][resNum] -= (potential[count+1] - potential[count])*2.5692/(dist*MV_to_au)
						gradxEF[diel][resNum] -= ((potential[count+2] - potential[count+1])-(potential[count+1] - potential[count]))*2.5692*100/(pow(dist,2)*MV_to_au_efg)	#MV/cm	
					
				flag += 1
				
			for resNum in resNums:
				if resNum in self.CSP_exp.keys():
					if self.nucleus == 'N':
						difCSP[diel][resNum] = abs(977*deltaEF[diel][resNum] - self.CSP_exp[resNum])
					if self.nucleus == 'H':
						difCSP[diel][resNum] = abs(188*deltaEF[diel][resNum] - self.CSP_exp[resNum])

					difCSP_total[diel] += difCSP[diel][resNum]
					difCSP_total_scale[diel] += difCSP[diel][resNum]/difCSP[self.diel_low][resNum]

			diel+=self.diel_increment
		
		if self.maps == 2:
			fd_dif = open("difference_Coulomb%.1f.txt"%active_diel, "w")
			fd_cs = open("absoluteCS_Coulomb%.1f.txt"%active_diel, "w")
		else:
			fd_dif = open("difference%.1f.txt"%active_diel, "w")
			fd_cs = open("absoluteCS%.1f.txt"%active_diel, "w")
			
		for resNum in resNums:
			
			try:
				count = index[resNum][self.nucleus]
			except:
				continue
				
			diel = self.diel_low
			line_dif = ''
			line_cs = ''
			flag_dif =0
			if resNum in self.CSP_exp.keys():
				flag_dif =1
				line_dif = '%5d'%resNum
				line_cs = '%5d'%resNum
			while (diel<=self.diel_high):
				if flag_dif ==1:
					line_dif += '  %10.6f'%difCSP[diel][resNum]
					if self.nucleus == 'N':
						cs_absolute = 977*deltaEF[diel][resNum]
					if self.nucleus == 'H':
						cs_absolute = 188*deltaEF[diel][resNum]
					line_cs += '  %10.6f'%cs_absolute
				diel+=self.diel_increment
			
			if line_dif != '':
				line_dif += '\n'
				fd_dif.write(line_dif)
			if line_cs != '':
				line_cs += '\n'
				fd_cs.write(line_cs)
		fd_cs.close()
		
		#write total difference in electric field.
		line_dif = 'total'
		line_dif_scale = 'scale'
		diel = self.diel_low
		while (diel<=self.diel_high):
			line_dif += '  %10.6f'%difCSP_total[diel]
			line_dif_scale += '  %10.6f'%difCSP_total_scale[diel]
			diel+=self.diel_increment
		line_dif += '\n'
		fd_dif.write(line_dif)
		fd_dif.write(line_dif_scale)
		fd_dif.close()
		return 
		
		
		
	def run(self):
		'''
		main driver for potential calcuation
		'''
		if self.maps == 1:
			active_diel = self.active_diel1
			while (active_diel <= self.active_diel2):
				self.calculatePotential(active_diel)
				active_diel += self.active_delta
			
		else:
			self.calculatePotential(active_diel=0.0)
		
		self.maps=2
		self.calculatePotential(active_diel=0.0)
		
		return
		
		
		
	
def startGhost():
	'''
	Function for starting Ghost calculations
	'''
	
	print 
	print 'Ghost titration predictions using explicit model'
	print

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
		'--atom',
		dest='atom',
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
		'--error_CS',
		dest='error_CS',
		type='string',
		default="{29:0.03, 30:0.02,32:0.03,33:0.03,34:0.03,38:0.03,42:0.02,44:0.02,45:0.02,46:0.05,51:0.03,52:0.06,54:0.03,55:0.03,56:0.04,58:0.08,59:0.04,60:0.04,91:0.04,107:0.02,108:0.04,109:0.03,111:0.03,113:0.02,114:0.02,115:0.02}",
		help='<dictionary of experimental pka values {tit.group:pka,...}>'
		)
	parser.add_option(
		'--threshold',
		dest='threshold',
		default=0.15,
		type='float',
		help='<NRMSD between calculated and experimental curves - default 0.15>'
		)
	parser.add_option(
		'--plot',
		dest='plot',
		default=1,
		type='choice',
		choices=('0','1'),
		help='<plot figures (0, 1)>'
		)
	parser.add_option(
		'--sdiel',
		dest='sdiel',
		type='float',
		default=80.0,
		help='<the surface dielectric constant>'
		)
	parser.add_option(
		'--diel1',
		dest='diel1',
		type='float',
		default=1.0,
		help='<down boundary of the dielectric constant>'
		)
	parser.add_option(
		'--diel2',
		dest='diel2',
		type='float',
		default=2.0,
		help='<upper boundary of the dielectric constant>'
		)
	parser.add_option(
		'--diel_delta',
		dest='diel_delta',
		type='float',
		default=1.0,
		help='<increment of dielectric values>'
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
		help='<0:PBE without 3d maps, 1:PBE using 3d maps in apbs, 2:Coulombs law using PBE pdiel=sdiel>'
		)
	parser.add_option(
		'--radius',
		dest='radius',
		type='string',
		default="[6.0]",
		help='<radii for defining the active site>'
		)
	parser.add_option(
		'--active_diel1',
		dest='active_diel1',
		type='float',
		default=4.0,
		help='<starting dielectric constant of the active site>'
		)
	parser.add_option(
		'--active_diel2',
		dest='active_diel2',
		type='float',
		default=4.0,
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
	(options, args,) = parser.parse_args()
	
	##
	##parse optparse options
	##
	pqrNeu=options.pqrNeu
	pqrCha=options.pqrCha
	atom=options.atom
	exper_CS=eval(options.exper_CS)
	threshold=options.threshold
	plot=options.plot
	sdiel = options.sdiel
	diel_low=options.diel1
	diel_high=options.diel2
	diel_increment=options.diel_delta
	ratio = options.ratio
	center = eval(options.center)
	error_CS=eval(options.error_CS)
	maps = options.maps
	radius = eval(options.radius)
	active_diel1 = options.active_diel1
	active_diel2 = options.active_diel2
	active_delta = options.active_delta
	active_coord = eval(options.active_coord)
	
	
		
	
	return pqrNeu, pqrCha, atom, exper_CS, threshold, plot, diel_low, diel_high, diel_increment, sdiel, ratio, center, error_CS, maps, radius, active_diel1, active_diel2, active_delta, active_coord
	
	
	
	
def main():
	
	pqrNeu, pqrCha, atom, exper_CS, threshold, plot, diel_low, diel_high, diel_increment, sdiel, ratio, center, error_CS, maps, radius, active_diel1, active_diel2, active_delta, active_coord = startGhost()
	Y = CS_apbs_analysis(pqrNeu, pqrCha, atom, exper_CS, threshold, plot, diel_low, diel_high, diel_increment, sdiel, ratio, center, error_CS, maps, radius, active_diel1, active_diel2, active_delta, active_coord)
	Y.run()
	
	
	
if __name__=='__main__':
	main()
		
