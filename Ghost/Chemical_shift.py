#!/usr/bin/python

#Chemical shift analysys class


import sys, os
import string
from pylab import *

class CS_analysis():
	'''Class CS_analysis predicts chemical shifts caused by "through space" Electric field effect.
	Predictions are made using PB calculations (APBS and Chresten's solver) and Coulomb's law.
	As an output dielectric map and map of moving residues are produced'''
	
	def __init__(self, pdb, ph_values, exper_pKa, method='Coulomb', atom='H'):
		'''pdb is pdb file of the protein whose chemical shifts are predicted
		pH_values is a list of pH values for which the predicted CS is calculated
		exper_pKa are pKa values of titratable groups extracted from titration curves
		method="PBE"
		method="APBS"
		method="Coulomb"
		atom is the name of atom which chemical shift is predicted
		atom='N'
		atom='H'
		atom='C'
		'''
		
		self.pdb=pdb
		self.ph_values=ph_values
		self.exper_pKa=exper_pKa
		self.method=method
		self.atom=atom
		self.diel_low=4.0
		self.diel_high=30.0
		self.diel_increment=0.5
		
		
		
	def format(value):
		if value:
			return '%6.3f' %value
		else:
			return 'NA'
			
			
	def run_Ghost(self,exper_CS=None):
		#
		# Main driver for predicting ghost titrations
		#
		self.get_dCS()
		self.calculate_shift()
		self.plot_curves(exper_CS)
		if exper_CS!=None:
			self.distance(exper_CS)
			self.set_newBfactor(threshold=0.15)
			
		return
			
			
	def get_dCS(self):
		"""Calculate the change in chemical shift due to a full charge on each titratable group
		and pickle the dictionary"""
		import cPickle
		import get_dEF
		X=get_dEF.map_struct(self.pdb)
		X.build_Hs()
		residues=X.PI.residues.keys()
		residues.sort()
		#
		# Get the titratable groups
		#
		titgroups=X.PI.get_titratable_groups()
		titgroups.sort()
		print '\nTitratable groups: ',titgroups
		self.residues=residues
		self.titgroups=titgroups
		
		
		directory=self.method+'_potmap/'
		dir_name=os.path.join(os.getcwd(),directory)
		dict_name=dir_name+self.atom+'chemical_shift_'+self.pdb[:-4]+'_'+self.method+'.dat'
		#
		#if dictionary havent been made yet the calculation is run
		#
		if not os.path.isfile(dict_name):
			import pKa.pKD_tools
			#
			#initialize the dictionary which collects all data in one place 
			#
			dict_CS={}
			diel=self.diel_low
			while diel<=self.diel_high:
				dict_CS[diel]={}
				for residue in self.residues:
					dict_CS[diel][residue]={}
				for titgroup in self.titgroups:
					#SHOULD INCLUDE PARALLELISM HERE
					titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
					charge=X.PI.titgroups[titgroup_type]['charge']
					if self.method=='Coulomb':
						for residue in self.residues:
							dCS=X.get_dCS(residue+':'+self.atom,titgroup, charge=charge, method=self.method, diel=float(diel))
							dict_CS[diel][residue][titgroup]=dCS
					else:
						residue=self.residues[1]
						dCS=X.get_dCS(residue+':'+self.atom,titgroup, charge=charge, method=self.method, diel=float(diel))
						for residue in self.residues:
							dict_CS[diel][residue][titgroup]=dCS[residue]
							
				diel+=self.diel_increment
			fd=open(dict_name,'w')
			cPickle.dump(dict_CS,fd)
			
		else:
			fd=open(dict_name,'r')
			print '\nLoading existing dictionary dict_CS:'
			print dict_name
			dict_CS=cPickle.load(fd)
			
		fd.close()
		self.dict_CS = dict_CS
		print "\nCalculating ghost titrations for dielectric constants in range (%4.2f,%4.2f) with step %4.2f\n" %(self.diel_low,self.diel_high,self.diel_increment)
		
	def calculate_shift(self):
		'''This function predicts chemical shifts of all residues according to pKa values of all titratable groups'''
		
		delta_CS={}
		for residue in self.residues:
			delta_CS[residue]={}
			diel=self.diel_low
			while diel <= self.diel_high:
				delta_CS[residue][diel]={}
				for ph in range(len(self.ph_values)):
					delta=0.00
					for titgroup in self.titgroups:
						if self.dict_CS[diel][residue][titgroup]!=(0.00 or None):
							if (titgroup[-3:]=='HIS') or (titgroup[-5:]=='NTERM') or (titgroup[-3:]=='LYS') or (titgroup[-3:]=='ARG'):
								delta+=(self.dict_CS[diel][residue][titgroup])/(1+10**(self.ph_values[ph]-self.exper_pKa[titgroup]))
							else:
								delta+=(self.dict_CS[diel][residue][titgroup])/(1+10**(self.exper_pKa[titgroup]-self.ph_values[ph]))
						else:
							pass
					delta_CS[residue][diel][ph]=delta
				diel+=self.diel_increment
		
		self.delta_CS=delta_CS
		
		
	def preprocess(self,dict):
		'''Finds the minimum of the curve and translate the curve downward for that value'''
		
		minimum=min(dict.values())
		for ph in range(len(self.ph_values)):
			dict[ph]-=minimum
		return dict.values()
		
		
		
	def plot_curves(self,exper_CS):
		'''This function plots for every residue both experimental CS curves and all calculated CS curves one for every dielectric constant
		Only the curves for residues which NMR data exist are plotted'''
		print "Ploting ghost titration curves in figures/"
		if exper_CS != None:
			residues=exper_CS.keys()
		else:
			residues=self.residues
		residues.sort()
		for residue in residues:
			clf()
			figure(num=1,dpi=400,figsize=(10,5),frameon=False)
			grid(True)
			hold=(True)
			if exper_CS != None:
				y=self.preprocess(exper_CS[residue])
				y=array(y)
				plot(self.ph_values,y,'r',linewidth=1.5)
				
			xlabel('pH values')
			diel=self.diel_low
			while diel <= self.diel_high:
				z=self.preprocess(self.delta_CS[residue][diel])
				z=array(z)
				plot(self.ph_values,z,'b',linewidth=0.8)
				diel+=self.diel_increment
			
			if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
				os.system('mkdir figures')
			dir_name=os.path.join(os.getcwd(),'figures/')
			file_name=dir_name+'%sres%s.png'%(self.pdb,residue)
			#axis([2, 8, 0, 2.5])
			savefig(file_name)
			
			
	def distance(self,exper_CS):
		'''This function finds the NRMSD between calculated and experimental curves'''
		
		
		distance={}
		ultimate_distance={}
		ultimate_diel={}
		residues=exper_CS.keys()
		residues.sort()
		for residue in residues:
			exper=self.preprocess(exper_CS[residue])
			maximum=max(exper)
			distance[residue]={}
			diel=self.diel_low
			while diel <= self.diel_high:
				delta=self.preprocess(self.delta_CS[residue][diel])
				dist=0.00
				for ph in range(len(self.ph_values)):
					dist+=(delta[ph]-exper[ph])**2
				distance[residue][diel]=sqrt(dist/len(self.ph_values))/maximum
				diel+=self.diel_increment
			
			ultimate_distance[residue]=min(distance[residue].values())
			diel=diel=self.diel_low
			while diel <= self.diel_high:
			#for diel in range(self.diel_low, self.diel_high+self.diel_increment, self.diel_increment):#
				if ultimate_distance[residue]==distance[residue][diel]:
					ultimate_diel[residue]=diel
					break
				else:
					diel+=self.diel_increment
		self.distance=distance
		self.ultimate_distance=ultimate_distance
		self.ultimate_diel=ultimate_diel
		
		
	def zero_Bfactor(self):
		'''function that zeros the value of Bfactor coloumn in pdb file'''
		
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name=os.path.join(os.getcwd(),'figures/')
		
		if os.path.isfile(self.pdb):
			fd=open(self.pdb)
			lines=fd.readlines()
			fd.close()
		else:
			raise FileNotFoundError, self.pdb
		
		file_name=dir_name+self.pdb[:-4]+'_zero'+'.pdb'
		fd=open(file_name,'w')
		
		for line in lines:
			type=string.strip(string.upper(line[:6]))
			if type=='ATOM':
				line=line[:60]+'  0.00'+line[66:]
			print >> fd, line
		fd.close()
		
		
	def set_newBfactor(self,threshold=0.15):
		'''reads distance values from ultimate_distance and put them as Bfactor value in the pdb file
		Threshold is the threshold nrmsd distance between calculated and experimental CS curves
		If nrmsd is bigger than threshold than the residue will be printed in output file'''
		
		dir_name=os.path.join(os.getcwd(),'figures/')
		
		self.zero_Bfactor()
		pdbfile=dir_name+self.pdb[:-4]+'_zero'+'.pdb'
		if os.path.isfile(pdbfile):
			fd=open(pdbfile)
			lines=fd.readlines()
			fd.close()
		
		residues_below=[]
		file_name=dir_name+self.pdb[:-4]+'_replace'+'.pdb'
		file_name1=dir_name+'residues_below_nrmsd_%9s.txt'%(self.pdb)
		fd=open(file_name,'w')
		fd1=open(file_name1,'w')
		
		for line in lines:
			type=string.strip(string.upper(line[:6]))
			if type=='ATOM':
				num=string.strip(line[22:27])
				residue='A:'+string.zfill(num,4)
				if self.ultimate_distance.has_key(residue):
					new_Bfactor=str(self.ultimate_distance[residue])[:4]
					line=line[:62]+new_Bfactor+line[66:]
					if self.ultimate_distance[residue]<=threshold:
						if not residue in residues_below:
							print >>fd1, residue,'%4.2f  %3.1f'%( self.ultimate_distance[residue], self.ultimate_diel[residue])
							residues_below.append(residue)
			print >> fd, line
		self.residues_below=residues_below
		fd.close()
		fd1.close()
		
		
	def accessibility(self):
		'''Calls whatif and save the accessibility of residues in a file'''
		
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name=os.path.join(os.getcwd(),'figures/')
		
		source='%swhatif_%s.txt'%(dir_name,self.pdb)
		fd=open(source,'w')
		fd.writelines(['/local/predrag/bin/whatif/DO_WHATIF.COM <<EOF\n',
			'GETMOL %s\n'%self.pdb,
			'\n',
			'%DELWAT\n',
			'%DELLIG\n',
			'%SETACC\n',
			'NOWAT 0\n',
			'\n',
			'%SHOACC\n',
			'TOT 0\n',
			'STOP\n',
			'Y\n'
			])
		fd.close()
		os.system('source %s > %saccessibility_%s.out' %(source,dir_name,self.pdb))
		fd=open('%saccessibility_%s.out' %(dir_name,self.pdb),'r')
		lines=fd.readlines()
		fd.close()
		for row in range(len(lines)):
			line=string.strip(lines[row])
			if line=='Res# Res    PDB#        Tot. Acc.     Back.     Side.':
				break
		fd=open('%saccessibility_%s.out' %(dir_name,self.pdb),'w')
		for i in range(row+1,len(lines)):
			line=string.strip(lines[i])
			if not line=='':
				print >> fd, line
			else:
				break
		fd.close()
		
	
	def print_accessibility(self,threshold=0.15,access='total'):
		'''Plot the influence of the accessibility on dielectric constant
		It does it only for residues which nrmsd is below the threshold'''
		if access=='total':
			colomn=-3
		elif access=='backbone':
			colomn=-2
		elif access=='side chain':
			colomn=-1
		else:
			print
			print 'Usage: accessibility could take only following values:'
			print 'total'
			print 'background'
			print 'side chain'
			os._exit(0)
			
		dir_name=os.path.join(os.getcwd(),'figures/')
		file_name='%saccessibility_%s.out' %(dir_name,self.pdb)
		if not os.path.isfile(file_name):
			self.accessibility()
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name=os.path.join(os.getcwd(),'figures/')
		fd=open('%saccessibility_%s.out' %(dir_name,self.pdb),'r')
		lines=fd.readlines()
		fd.close()
		acc_total={}
		acc_side={}
		acc_back={}
		plot2D={}
		plot2D_1={}
		for line in lines:
			line=string.strip(line)
			num=string.split(line)[3]
			residue='A:'+string.zfill(num,4)
			acc_total[residue]=float(string.split(line)[-3])
			acc_back[residue]=float(string.split(line)[-2])
			acc_side[residue]=float(string.split(line)[-1])
			
		for residue in self.residues:
			if (residue in self.ultimate_distance.keys()) and (self.ultimate_distance[residue]<=threshold):
				if access=='total':
					acc=acc_total[residue]
				elif access=='backbone':
					acc=acc_back[residue]
				elif access=='side chain':
					acc=acc_side[residue]
				plot2D[acc]=self.ultimate_diel[residue]
				
				if self.ultimate_distance[residue]<=threshold/2:
					plot2D_1[acc]=self.ultimate_diel[residue]
					
		clf()
		figure(num=1,dpi=400,figsize=(10,5),frameon=False)
		grid(True)
		hold(True)
		xlabel('%s access.[A^2]'%access)
		ylabel('effective dielectric constant')
		plot(plot2D.keys(),plot2D.values(),'bo',label='%4s<NRMSD<=%4s'%(str(threshold/2),str(threshold)))
		plot(plot2D_1.keys(),plot2D_1.values(),'ro',label='NRMSD<=%4s'%str(threshold/2))
		
		file_name=dir_name+'%s_accessibility_vs_diel.png'%access
		legend()
		savefig(file_name)
		
		
		#3D plot
		X=[]
		Y=[]
		Z=[]
		X1=[]
		Y1=[]
		Z1=[]
		for residue in self.residues:
			if (residue in self.ultimate_distance.keys()) and (self.ultimate_distance[residue]<=threshold):
				if self.ultimate_distance[residue]<=threshold/2:
					X1.append(acc_back[residue])
					Z1.append(self.ultimate_diel[residue])
				else:
					X.append(acc_back[residue])
					Z.append(self.ultimate_diel[residue])
				interactions=[]
				for titgroup in self.titgroups:
					if (titgroup[:6] in self.ultimate_distance.keys()):
						if (titgroup[-3:]=='HIS') or (titgroup[-3:]=='ASP') or (titgroup[-3:]=='GLU'):
							diel=self.ultimate_diel[residue]
							interactions.append(self.dict_CS[diel][residue][titgroup])
				interactions.sort()
				for titgroup in self.titgroups:
					if interactions[0]==self.dict_CS[diel][residue][titgroup]:
						if self.ultimate_distance[residue]<=threshold/2:
							Y1.append(acc_side[titgroup[:6]])
						else:
							Y.append(acc_side[titgroup[:6]])
						break

		clf()
		import matplotlib.axes3d as p3
		fig=figure(num=1,dpi=400,figsize=(10,5),frameon=False)
		ax = p3.Axes3D(fig)
		ax.scatter3D(array(X),array(Y),array(Z),color='b')
		ax.scatter3D(array(X1),array(Y1),array(Z1),color='r')
		ax.set_xlabel('back access residue')
		ax.set_ylabel('side chain access titgroup')
		ax.set_zlabel('effective diel')
		file_name=dir_name+'3Daccessibility_vs_diel.png'
		savefig(file_name)
		
		
	def plot_individual_shift(self):
		'''Plot chemical shifts due to a full charge on each titratable group individually
		Also plots the unexpected curve shapes - curves that do not decline while increasing the pH value'''
		
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name=os.path.join(os.getcwd(),'figures/')
		
		for residue in self.residues:
			clf()
			figure(num=1,dpi=400,figsize=(10,5),frameon=False)
			grid(True)
			hold(True)
			xlabel('effective dielectric constant')
			for titgroup in self.titgroups:
				diel=self.diel_low
				y=[]
				m=0
				while diel <= self.diel_high:
					if self.dict_CS[diel][residue][titgroup]==None:
						m=1
						break
					else:
						y.append(self.dict_CS[diel][residue][titgroup])
						diel+=self.diel_increment
				if m==0:
					plot(y,'b',linewidth=0.8)
			file_name=dir_name+'individual_shift_%s.png'%residue
			savefig(file_name)
	
	def important_titgroup(self,level=0.5):
		'''For all nuclei finds the most influential titratable groups on a given pH value'''
		
		if not os.path.isdir(os.path.join(os.getcwd(),'figures/')):
			os.system('mkdir figures')
		dir_name=os.path.join(os.getcwd(),'figures/')
		
		CS={}
		file_name=dir_name+'important_titgroup_%9s.txt'%(self.pdb)
		fd=open(file_name,'w')
		for residue in self.residues_below:
			print >> fd, '\nResidue:\n',residue
			diel=self.ultimate_diel[residue]
			print >> fd, 'diel_const: ', diel
			print >> fd, '\nTit group    value\n'
			CS[residue]={}
			for titgroup in self.titgroups:
				if self.dict_CS[diel][residue][titgroup]!=(0.00 or None):
					if (titgroup[-3:]=='GLU') or (titgroup[-5:]=='NTERM') or (titgroup[-5:]=='CTERM') or (titgroup[-3:]=='ASP'):
						CS[residue][titgroup]=abs(self.dict_CS[diel][residue][titgroup])
						if CS[residue][titgroup]>level:
							print >> fd, '%12s %4.2f'%(titgroup,CS[residue][titgroup]) 
				else:
					pass
			
		fd.close()
		




def main():
	#need pdb file and experimental chemical shift data extracted from Ekin
	#pdbs=['PKA.1BEB_G64D.pdb','1BEB_dimer_G64D.pdb','1DV9-model0.pdb','PKA.2BLG_dimer.pdb']
	pdbs=['1BEB_dimer_G64D.pdb']
	for pdb in pdbs:
		ph_values=[2.38,2.88,3.38,3.88,4.38,4.88,5.38,5.88,6.38,6.88,7.38,7.88]
		#pKa 1BEB dimer#
		#exper_pKa={ 'A:0002:TYR':15.00,'A:0001:NTERM': 8.24,'A:0005:GLN:NTERM':15.00,'A:0005:NTERM':15.00,'A:0162:CTERM': 3.75,'A:0160:CTERM':0.00,'A:0160:CYS:CTERM':0.00,'A:0141:LYS': 11.80, 'A:0040:ARG': 15.00, 'A:0112:GLU': 5.38, 'A:0070:LYS': 11.70, 'A:0028:ASP': 4.69, 'A:0051:GLU': 3.85, 'A:0074:GLU': 5.16, 'A:0085:ASP': 3.62, 'A:0137:ASP':2.67, 'A:0042:TYR': 9.80, 'A:0020:TYR': 9.80, 'A:0075:LYS': 12.50, 'A:0157:GLU': 6.05, 'A:0161:HIS': 7.16, 'A:0114:GLU': 3.36, 'A:0045:GLU': 3.99, 'A:0062:GLU': 1.41, 'A:0101:LYS': 10.70, 'A:0135:LYS': 13.00, 'A:0064:ASP': 4.25, 'A:0065:GLU': 5.24, 'A:0158:GLU': 5.03, 'A:0102:TYR': 9.80, 'A:0129:ASP': 3.85, 'A:0098:ASP': 2.53, 'A:0130:ASP': 4.75, 'A:0127:GLU': 4.66, 'A:0014:LYS': 10.80, 'A:0033:ASP': 1.70, 'A:0096:ASP': 2.22, 'A:0108:GLU': 2.75, 'A:0069:LYS': 10.80, 'A:0138:LYS': 11.20, 'A:0055:GLU': 2.06, 'A:0047:LYS': 13.00, 'A:0008:LYS': 10.80, 'A:0099:TYR': 9.80, 'A:0060:LYS': 11.10, 'A:0100:LYS': 10.30, 'A:0148:ARG': 17.50, 'A:0124:ARG': 16.40, 'A:0134:GLU': 3.10, 'A:0053:ASP': 3.99, 'A:0091:LYS': 11.40, 'A:0089:GLU': 10.29, 'A:0077:LYS': 10.50, 'A:0146:HIS': 7.36, 'A:0044:GLU': 6.69, 'A:0131:GLU': 5.42, 'A:0083:LYS': 10.20, 'A:0011:ASP': 2.89, 'B:0001:NTERM': 8.24,'B:0005:GLN:NTERM':15.00,'B:0005:NTERM':15.00,'B:0162:CTERM': 3.75,'B:0160:CTERM':0.00,'B:0160:CYS:CTERM':0.00,'B:0141:LYS': 11.80, 'B:0040:ARG': 15.00, 'B:0112:GLU': 5.38, 'B:0070:LYS': 11.70, 'B:0028:ASP': 4.69, 'B:0051:GLU': 3.85, 'B:0074:GLU': 5.16, 'B:0085:ASP': 3.62, 'B:0137:ASP':2.67, 'B:0042:TYR': 9.80, 'B:0020:TYR': 9.80, 'B:0075:LYS': 12.50, 'B:0157:GLU': 6.05, 'B:0161:HIS': 7.16, 'B:0114:GLU': 3.36, 'B:0045:GLU': 3.99, 'B:0062:GLU': 1.41, 'B:0101:LYS': 10.70, 'B:0135:LYS': 13.00, 'B:0064:ASP': 4.25, 'B:0065:GLU': 5.24, 'B:0158:GLU': 5.03, 'B:0102:TYR': 9.80, 'B:0129:ASP': 3.85, 'B:0098:ASP': 2.53, 'B:0130:ASP': 4.75, 'B:0127:GLU': 4.66, 'B:0014:LYS': 10.80, 'B:0033:ASP': 1.70, 'B:0096:ASP': 2.22, 'B:0108:GLU': 2.75, 'B:0069:LYS': 10.80, 'B:0138:LYS': 11.20, 'B:0055:GLU': 3.01, 'B:0047:LYS': 13.00, 'B:0008:LYS': 10.80, 'B:0099:TYR': 9.80, 'B:0060:LYS': 11.10, 'B:0100:LYS': 10.30, 'B:0148:ARG': 17.50, 'B:0124:ARG': 16.40, 'B:0134:GLU': 3.10, 'B:0053:ASP': 3.99, 'B:0091:LYS': 11.40, 'B:0089:GLU': 10.29, 'B:0077:LYS': 10.50, 'B:0146:HIS': 7.36, 'B:0044:GLU': 6.69, 'B:0131:GLU': 5.42, 'B:0083:LYS': 10.20, 'B:0011:ASP': 2.89}
		#pKa 1BEB dimer#only from calculated
		exper_pKa={ 'A:0002:TYR':15.00,'A:0001:NTERM': 8.24,'A:0005:GLN:NTERM':15.00,'A:0005:NTERM':15.00,'A:0162:CTERM': 3.75,'A:0160:CTERM':0.00,'A:0160:CYS:CTERM':0.00,'A:0141:LYS': 11.80, 'A:0040:ARG': 15.00, 'A:0112:GLU': 4.91, 'A:0070:LYS': 11.70, 'A:0028:ASP': 4.69, 'A:0051:GLU': 4.01, 'A:0074:GLU': 4.74, 'A:0085:ASP': 3.62, 'A:0137:ASP':0.44, 'A:0042:TYR': 9.80, 'A:0020:TYR': 9.80, 'A:0075:LYS': 12.50, 'A:0157:GLU': 6.05, 'A:0161:HIS': 7.16, 'A:0114:GLU': 4.21, 'A:0045:GLU': 4.91, 'A:0062:GLU': 1.41, 'A:0101:LYS': 10.70, 'A:0135:LYS': 13.00, 'A:0064:ASP': 4.03, 'A:0065:GLU': 4.54, 'A:0158:GLU': 5.03, 'A:0102:TYR': 9.80, 'A:0129:ASP': 3.85, 'A:0098:ASP': 1.65, 'A:0130:ASP': 3.41, 'A:0127:GLU': 4.38, 'A:0014:LYS': 10.80, 'A:0033:ASP': 1.70, 'A:0096:ASP': 2.56, 'A:0108:GLU': 1.47, 'A:0069:LYS': 10.80, 'A:0138:LYS': 11.20, 'A:0055:GLU': 2.06, 'A:0047:LYS': 13.00, 'A:0008:LYS': 10.80, 'A:0099:TYR': 9.80, 'A:0060:LYS': 11.10, 'A:0100:LYS': 10.30, 'A:0148:ARG': 17.50, 'A:0124:ARG': 16.40, 'A:0134:GLU': 3.79, 'A:0053:ASP': 4.28, 'A:0091:LYS': 11.40, 'A:0089:GLU': 10.29, 'A:0077:LYS': 10.50, 'A:0146:HIS': 8.23, 'A:0044:GLU': 5.49, 'A:0131:GLU': 4.57, 'A:0083:LYS': 10.20, 'A:0011:ASP': 3.05, 'B:0001:NTERM': 8.24,'B:0005:GLN:NTERM':15.00,'B:0005:NTERM':15.00,'B:0162:CTERM': 3.75,'B:0160:CTERM':0.00,'B:0160:CYS:CTERM':0.00,'B:0141:LYS': 11.80, 'B:0040:ARG': 15.00, 'B:0112:GLU': 5.38, 'B:0070:LYS': 11.70, 'B:0028:ASP': 4.69, 'B:0051:GLU': 3.85, 'B:0074:GLU': 5.16, 'B:0085:ASP': 3.62, 'B:0137:ASP':2.67, 'B:0042:TYR': 9.80, 'B:0020:TYR': 9.80, 'B:0075:LYS': 12.50, 'B:0157:GLU': 6.05, 'B:0161:HIS': 7.16, 'B:0114:GLU': 3.36, 'B:0045:GLU': 3.99, 'B:0062:GLU': 1.41, 'B:0101:LYS': 10.70, 'B:0135:LYS': 13.00, 'B:0064:ASP': 4.25, 'B:0065:GLU': 5.24, 'B:0158:GLU': 5.03, 'B:0102:TYR': 9.80, 'B:0129:ASP': 3.85, 'B:0098:ASP': 2.53, 'B:0130:ASP': 4.75, 'B:0127:GLU': 4.66, 'B:0014:LYS': 10.80, 'B:0033:ASP': 1.70, 'B:0096:ASP': 2.22, 'B:0108:GLU': 2.75, 'B:0069:LYS': 10.80, 'B:0138:LYS': 11.20, 'B:0055:GLU': 3.01, 'B:0047:LYS': 13.00, 'B:0008:LYS': 10.80, 'B:0099:TYR': 9.80, 'B:0060:LYS': 11.10, 'B:0100:LYS': 10.30, 'B:0148:ARG': 17.50, 'B:0124:ARG': 16.40, 'B:0134:GLU': 3.10, 'B:0053:ASP': 3.99, 'B:0091:LYS': 11.40, 'B:0089:GLU': 10.29, 'B:0077:LYS': 10.50, 'B:0146:HIS': 7.36, 'B:0044:GLU': 6.69, 'B:0131:GLU': 5.42, 'B:0083:LYS': 10.20, 'B:0011:ASP': 2.89}
		#pKa 2BLG
		#exper_pKa={ 'A:0002:TYR':15.00,'A:0001:NTERM': 8.24,'A:0005:GLN:NTERM':15.00,'A:0005:NTERM':15.00,'A:0162:CTERM': 3.75,'A:0160:CTERM':0.00,'A:0160:CYS:CTERM':0.00,'A:0141:LYS': 11.80, 'A:0040:ARG': 15.00, 'A:0112:GLU': 5.38, 'A:0070:LYS': 11.70, 'A:0028:ASP': 2.14, 'A:0051:GLU': 3.85, 'A:0074:GLU': 5.16, 'A:0085:ASP': 7.52, 'A:0137:ASP':2.67, 'A:0042:TYR': 9.80, 'A:0020:TYR': 9.80, 'A:0075:LYS': 12.50, 'A:0157:GLU': 4.94, 'A:0161:HIS': 7.16, 'A:0114:GLU': 3.36, 'A:0045:GLU': 3.99, 'A:0062:GLU': 3.34, 'A:0101:LYS': 10.70, 'A:0135:LYS': 13.00, 'A:0064:ASP': 4.25, 'A:0065:GLU': 5.24, 'A:0158:GLU': 6.65, 'A:0102:TYR': 9.80, 'A:0129:ASP': 2.21, 'A:0098:ASP': 2.53, 'A:0130:ASP': 4.75, 'A:0127:GLU': 4.66, 'A:0014:LYS': 10.80, 'A:0033:ASP': 6.10, 'A:0096:ASP': 2.22, 'A:0108:GLU': 2.75, 'A:0069:LYS': 10.80, 'A:0138:LYS': 11.20, 'A:0055:GLU': 3.52, 'A:0047:LYS': 13.00, 'A:0008:LYS': 10.80, 'A:0099:TYR': 9.80, 'A:0060:LYS': 11.10, 'A:0100:LYS': 10.30, 'A:0148:ARG': 17.50, 'A:0124:ARG': 16.40, 'A:0134:GLU': 3.10, 'A:0053:ASP': 3.99, 'A:0091:LYS': 11.40, 'A:0089:GLU': 3.74, 'A:0077:LYS': 10.50, 'A:0146:HIS': 7.36, 'A:0044:GLU': 6.69, 'A:0131:GLU': 5.42, 'A:0083:LYS': 10.20, 'A:0011:ASP': 2.89, 'B:0001:NTERM': 8.24,'B:0005:GLN:NTERM':15.00,'B:0005:NTERM':15.00,'B:0162:CTERM': 3.75,'B:0160:CTERM':0.00,'B:0160:CYS:CTERM':0.00,'B:0141:LYS': 11.80, 'B:0040:ARG': 15.00, 'B:0112:GLU': 5.38, 'B:0070:LYS': 11.70, 'B:0028:ASP': 2.14, 'B:0051:GLU': 3.85, 'B:0074:GLU': 5.16, 'B:0085:ASP': 7.52, 'B:0137:ASP':2.67, 'B:0042:TYR': 9.80, 'B:0020:TYR': 9.80, 'B:0075:LYS': 12.50, 'B:0157:GLU': 4.94, 'B:0161:HIS': 7.16, 'B:0114:GLU': 3.36, 'B:0045:GLU': 3.99, 'B:0062:GLU': 3.34, 'B:0101:LYS': 10.70, 'B:0135:LYS': 13.00, 'B:0064:ASP': 4.25, 'B:0065:GLU': 5.24, 'B:0158:GLU': 6.65, 'B:0102:TYR': 9.80, 'B:0129:ASP': 2.21, 'B:0098:ASP': 2.53, 'B:0130:ASP': 4.75, 'B:0127:GLU': 4.66, 'B:0014:LYS': 10.80, 'B:0033:ASP': 6.10, 'B:0096:ASP': 2.22, 'B:0108:GLU': 2.75, 'B:0069:LYS': 10.80, 'B:0138:LYS': 11.20, 'B:0055:GLU': 3.01, 'B:0047:LYS': 13.00, 'B:0008:LYS': 10.80, 'B:0099:TYR': 9.80, 'B:0060:LYS': 11.10, 'B:0100:LYS': 10.30, 'B:0148:ARG': 17.50, 'B:0124:ARG': 16.40, 'B:0134:GLU': 3.10, 'B:0053:ASP': 3.99, 'B:0091:LYS': 11.40, 'B:0089:GLU': 3.47, 'B:0077:LYS': 10.50, 'B:0146:HIS': 7.36, 'B:0044:GLU': 6.69, 'B:0131:GLU': 5.42, 'B:0083:LYS': 10.20, 'B:0011:ASP': 2.89}
		#pKa 2BLG#only from calculated
		#exper_pKa={ 'A:0002:TYR':15.00,'A:0001:NTERM': 8.24,'A:0005:GLN:NTERM':15.00,'A:0005:NTERM':15.00,'A:0162:CTERM': 3.75,'A:0160:CTERM':0.00,'A:0160:CYS:CTERM':0.00,'A:0141:LYS': 11.80, 'A:0040:ARG': 15.00, 'A:0112:GLU': 4.32, 'A:0070:LYS': 11.70, 'A:0028:ASP': 2.14, 'A:0051:GLU': 4.35, 'A:0074:GLU': 4.57, 'A:0085:ASP': 7.52, 'A:0137:ASP':1.15, 'A:0042:TYR': 9.80, 'A:0020:TYR': 9.80, 'A:0075:LYS': 12.50, 'A:0157:GLU': 4.94, 'A:0161:HIS': 7.16, 'A:0114:GLU': 4.22, 'A:0045:GLU': 4.29, 'A:0062:GLU': 3.34, 'A:0101:LYS': 10.70, 'A:0135:LYS': 13.00, 'A:0064:ASP': 4.29, 'A:0065:GLU': 4.29, 'A:0158:GLU': 6.65, 'A:0102:TYR': 9.80, 'A:0129:ASP': 2.21, 'A:0098:ASP': 1.73, 'A:0130:ASP': 4.12, 'A:0127:GLU': 4.38, 'A:0014:LYS': 10.80, 'A:0033:ASP': 6.10, 'A:0096:ASP': 2.92, 'A:0108:GLU': 1.50, 'A:0069:LYS': 10.80, 'A:0138:LYS': 11.20, 'A:0055:GLU': 3.52, 'A:0047:LYS': 13.00, 'A:0008:LYS': 10.80, 'A:0099:TYR': 9.80, 'A:0060:LYS': 11.10, 'A:0100:LYS': 10.30, 'A:0148:ARG': 17.50, 'A:0124:ARG': 16.40, 'A:0134:GLU': 3.59, 'A:0053:ASP': 4.48, 'A:0091:LYS': 11.40, 'A:0089:GLU': 3.74, 'A:0077:LYS': 10.50, 'A:0146:HIS': 7.55, 'A:0044:GLU': 5.92, 'A:0131:GLU': 3.55, 'A:0083:LYS': 10.20, 'A:0011:ASP': 3.41, 'B:0001:NTERM': 8.24,'B:0005:GLN:NTERM':15.00,'B:0005:NTERM':15.00,'B:0162:CTERM': 3.75,'B:0160:CTERM':0.00,'B:0160:CYS:CTERM':0.00,'B:0141:LYS': 11.80, 'B:0040:ARG': 15.00, 'B:0112:GLU': 5.38, 'B:0070:LYS': 11.70, 'B:0028:ASP': 2.14, 'B:0051:GLU': 3.85, 'B:0074:GLU': 5.16, 'B:0085:ASP': 7.52, 'B:0137:ASP':2.67, 'B:0042:TYR': 9.80, 'B:0020:TYR': 9.80, 'B:0075:LYS': 12.50, 'B:0157:GLU': 4.94, 'B:0161:HIS': 7.16, 'B:0114:GLU': 3.36, 'B:0045:GLU': 3.99, 'B:0062:GLU': 3.34, 'B:0101:LYS': 10.70, 'B:0135:LYS': 13.00, 'B:0064:ASP': 4.25, 'B:0065:GLU': 5.24, 'B:0158:GLU': 6.65, 'B:0102:TYR': 9.80, 'B:0129:ASP': 2.21, 'B:0098:ASP': 2.53, 'B:0130:ASP': 4.75, 'B:0127:GLU': 4.66, 'B:0014:LYS': 10.80, 'B:0033:ASP': 6.10, 'B:0096:ASP': 2.22, 'B:0108:GLU': 2.75, 'B:0069:LYS': 10.80, 'B:0138:LYS': 11.20, 'B:0055:GLU': 3.01, 'B:0047:LYS': 13.00, 'B:0008:LYS': 10.80, 'B:0099:TYR': 9.80, 'B:0060:LYS': 11.10, 'B:0100:LYS': 10.30, 'B:0148:ARG': 17.50, 'B:0124:ARG': 16.40, 'B:0134:GLU': 3.10, 'B:0053:ASP': 3.99, 'B:0091:LYS': 11.40, 'B:0089:GLU': 3.47, 'B:0077:LYS': 10.50, 'B:0146:HIS': 7.36, 'B:0044:GLU': 6.69, 'B:0131:GLU': 5.42, 'B:0083:LYS': 10.20, 'B:0011:ASP': 2.89}
		import cPickle
		fd=open('/home/people/predrag/python/pKaTool/Ghost/experim_chem_shift_N.dat')
		exper_CS=cPickle.load(fd)
		fd.close()
		Y=CS_analysis(pdb,ph_values,exper_pKa,method='APBS', atom='N')
		Y.run_Ghost()
		
		#Y.print_accessibility(threshold=0.15,access='total')
		#Y.important_titgroup(level=0.1)
		#Y.plot_individual_shift()
if __name__=='__main__':
	main()
