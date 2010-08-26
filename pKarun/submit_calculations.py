#!/usr/bin/python

'''
Script that automatically runs WHATIF pKa calcualtions on the cluster

			Predrag Kukic
			Jens Erik Nielsen
'''




import optparse
import sys, os
import string


def callWhatif(pdb):
	'''
	Runs Whatif and cleans the pdb file
	'''
	dir_name = os.getcwd()
	source = '%s/whatif.txt'%dir_name
	fd = open(source, 'w')
	fd.writelines(['/home/predrag/whatif/DO_WHATIF.COM <<EOF\n',
		'GETMOL %s\n'%pdb,
		'\n',
		'%DELWAT\n',
		'%DELLIG\n',
		'%CORALL\n',
		'N\n',
		'%MAKMOL\n',
		'\n',
		'pka.%s\n'%pdb,
		'TOT 0\n',
		'\n',
		'STOP\n',
		'Y\n'
		])
	fd.close()
	os.system('source %s > whatif.out'%source)
	new_pdb = 'pka.%s'%pdb,
	return new_pdb




def callPkarun(pdb):
	'''
	Runs pKa calculations
	'''
	dir_name = os.getcwd()
	source = '%s/pKarun.txt'%dir_name
	fd = open(source, 'w')
	fd.writelines(['ln -s /home/predrag/whatif/dbdata/TOPOLOGY.H .\n',
	'ln -s /home/predrag/whatif/dbdata/DELRAD.DAT .\n',
	'ln -s /home/predrag/whatif/dbdata/DELCRG.DAT .\n',
	'ln -s /home/predrag/whatif/dbdata/OPLS.RAD .\n',
	'ln -s /home/predrag/whatif/dbdata/OPLS.CRG .\n',
	'python /home/predrag/python/pKarun/pKarun_main.py %s -dbcrit 1000 -indi 8 -auto 2\n'%pdb,
	])
	fd.close()
	os.system('source %s' %source)
	





def startCalculations():
	'''
	parses the pdb file specified
	'''
	parser = optparse.OptionParser()
	
	#
	#set optparse options
	#
	parser.add_option(
	'--pdb',
	dest = 'pdb',
	default = None,
	type = 'string',
	help = '<pdb file for pKa calculations>',
	)
	
	(options, args, ) = parser.parse_args()
	
	#
	#parse optparse options
	#
	pdb = options.pdb
	
	return pdb
	
	
	
	
def main():
	pdb = startCalculations()
	pdb = callWhatif(pdb)
	callPkarun(pdb)
	
if __name__=='__main__':
	main()


