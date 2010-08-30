#! /bin/env python
'''Tool for filtering a csv file where each row corresponds to data for a protein mutant.
There are four filter types possible - 
Positions: Match entries based on residue indexes;
Mutations: Match entries based on mutations e.g. 46L.;
Substitutions: Match entries based on substitutions e.g. ALA,GLY.;
GroupSubs: Match entries based on group exchanges e.g. charged-polar,nonpolar-polar.
Exacltly how a record is determined to have matched the criteria is determined by the match-mode flag.
'''

import optparse, sys
import PEAT_SA.Core as Core

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-f", "--file", dest="file",
                  help="A csv file containing a column labled mutations", metavar="FILE")
parser.add_option("-s", "--selectionMode",
                  dest="filterMode",
	          default='Positions',
                  help="One of 'Positions', 'Mutations', 'Substitutions' or 'GroupSubs'",
                  metavar='SELECTION_MODE')		  		  		  		  		  		  
parser.add_option("-c", "--criteria",
                  dest="criteria",
                  help="A comman separated list of the criteria for the given filter mode. "
		  "This will be residue codes for 'Mutations' e.g. ALA,GLY, indexes for Positions e.g. 13,14 etc",
                  metavar='CRITERIA')		  		  		  		  		  		  
parser.add_option("-n", "--negate", action="store_true",
		  dest="negate", default=False,
		  help="Replace the specified criteria with its opposite. E.g. specifying '-s Mutations -c ALA,GLY -n' "
		  "replaces 'ALA,GLY' with all amino-acids codes EXCEPT ALA and GLY.")
parser.add_option("-m", "--matchMode",
                  dest="matchMode", default="Exclusive",
                  help="The match mode defines exactly how a record is determined to have met the criteria "
		  "Values: 'Any', 'All', 'Exclusive'\nDefault %default",
                  metavar='MATCH_MODE')		  		  		  		  		  		  
parser.add_option("", "--chains",
                  dest="chains",
                  help="Restrict filtering to a specified chain (i.e. for homodimers etc.)",
                  metavar='CHAINS')		  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if not sys.stdin.isatty():
	csv = sys.stdin.read()
	matrix = Core.Matrix.matrixFromCSVRepresentation(csv)
elif options.file is None:
        print 'CSV file must be provided'
        sys.exit(1)
else:
	matrix = Core.Matrix.matrixFromCSVFile(options.file)
	
posmode = False
submode = False
mutmode = False
exmode = False
if options.filterMode == 'Positions':
	posmode = True
elif options.filterMode ==  'Substitutions':
	submode = True
elif options.filterMode == 'Mutations':
	mutmode = True
elif options.filterMode == 'GroupSubs':
	exmode = True
else:
	print 'Error - Unknown filter mode, %s, specified' % options.filterMode
	sys.exit(1)

if posmode and options.chains is None:
	#This is to be backward compatible with some HIV protease scripts
	options.chains = 'A,B'

if posmode is True:
	indexes = options.criteria.split(',')
	chains = options.chains.split(',')
	positions = []
	for chain in chains:
		positions.extend([chain+index for index in indexes]) 

	if options.negate is True:
		residues = matrix.mutatedResidues()
		for position in positions:
			try:
				residues.remove(position)
			except ValueError:
				sys.stderr.write('Position %s not mutated in any entry - skipping\n' % position)

		positions = residues

	submatrix = matrix.entriesWithMutatedResidues(positions, matchMode=options.matchMode)

elif submode is True:
	subs = options.criteria.split(',')
	if options.negate is True:
		allSubs = matrix.substitutions()
		for sub in subs:
			allSubs.remove(sub)
		subs = allSubs 

	submatrix = matrix.entriesWithSubstitutions(subs, matchMode=options.matchMode)

elif mutmode is True:
	mutations = options.criteria.split(',')
	if options.negate is True:
		allMutations = matrix.mutationCodes()
		for mutation in mutations:
			allMutations.remove(mutation)
		mutations = allMutations 

	submatrix = matrix.entriesWithMutations(mutations, matchMode=options.matchMode)

elif exmode is True:
	exchanges = [el.strip() for el in options.criteria.split(',')]
	exchanges = [el.split('-') for el in exchanges]
	if options.negate is True:
		allExchanges = matrix.typeSubstitutions()
		allExchanges = [el.split('-') for el in allExchanges] 
		for exchange in exchanges:
			allExchanges.remove(exchange)
		exchanges = allExchanges 

	submatrix = matrix.entriesWithTypeSubstitutions(exchanges, matchMode=options.matchMode)

if submatrix is not None:
	print submatrix.csvRepresentation(),
	sys.stderr.write('Percentage extracted %5.3lf\n' % (100.0*submatrix.numberOfRows()/matrix.numberOfRows()))
else:
	sys.stderr.write('Error: Nothing passed filter\n')
