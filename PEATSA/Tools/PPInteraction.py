#! /usr/bin/env python
'''Calculates Delta Delta G_bind(WT->M) for a set of mutants of a protein-protein complex using PEATSA data

Result is a csv file containing the Delta Delta G_bind(WT->M) values for each mutant

NOTE: The chain IDs and residue numbering of the proteins in the complex and individual
results files must be identical.'''
import PEATSA.Core as Core
import optparse, sys, operator
import Protool

def MutationsInPartner(code, partner):

	retval = True

	try:
		data = partner.dataForMutation(code)
	except IndexError:
		retval = False

	return retval

def GetCrossParternData(set, pdb, partner1, partner2):

	codes = set.chainCodes(reduced=True, pdb=pdb)
	print 'Per chain codes:',
	print codes
	for code in codes:
		chain = code[0]
 
	values = []
	for code in codes:
		try: 
			data = partner1.dataForMutation(code)
			values.append(data[-1])
			#print '\tFound data for %s in partner 1: ' % code,
			#print data
		except IndexError:
			try:
				data = partner2.dataForMutation(code)
				values.append(data[-1])
				#print '\tFound data for %s in partner 2: ' % code, 
				#print data
			except IndexError:
				pass

	return values


if __name__ == "__main__":

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-p", "--protein", dest="protein",
			  help="File containing pdb structure of complex.", metavar="PROTEIN")
	parser.add_option("-c", "--complex-stability", dest="complex",
			  help="Csv file containing complex stability results.", metavar="COMPLEX")
	parser.add_option("-x", "--partner1-stability",
			  dest="partner1",
			  help="CSV file containing the stability result for the (free) first partner of the complex",
			  metavar='PARTNER_1')
	parser.add_option("-y", "--partner2-stability",
			  dest="partner2",
			  help="CSV file containing the stability result for the (free) second partner of the complex",
			  metavar='PARTNER_2')
			  
	(options, args) = parser.parse_args()

	run = True
	if options.protein is None:
		print 'Structure of the complex must be provided'
		run = False

	if options.complex is None:
		print 'Stability results for the complex must be provided'
		run = False
			  
	if options.partner1 is None:
		print 'Stability results for the first partner must be provided'
		run = False
		
	if options.partner2 is None:
		print 'Stability results for the second partner must be provided'
		run = False

	if not run:
		
		sys.exit(1)

	complex = Core.Matrix.matrixFromCSVFile(options.complex)
	partner1 = Core.Matrix.matrixFromCSVFile(options.partner1)
	partner2 = Core.Matrix.matrixFromCSVFile(options.partner2)
	pdb = Protool.structureIO()
	pdb.readpdb(options.protein)

	combined = zip(complex.mutations, complex.total)
	print 'Data found for %d mutants of the complex' % len(combined)

	interactionEnergies = []
	for element in combined:
		set = Core.Data.MutationSet(code=element[0])
		code = "+".join(set.reducedMutationCodes(pdb))

		values = []
		
		if MutationsInPartner(code, partner1):
			data = partner1.dataForMutation(code)
			values.append(data[-1])
		elif MutationsInPartner(code, partner2): 
			data = partner2.dataForMutation(code)
			values.append(data[-1])
		else:
			values = GetCrossParternData(set, pdb, partner1, partner2)	

		if len(values) == 0:
			print 'No free-partner data available for mutation %s' % set
			print 'The chain ids and residue numbering of the proteins in the complex and free structures must be identicial'
			raise BaseException

		interactionEnergy = element[1] -  reduce(operator.add, values)
		interactionEnergies.append(interactionEnergy)

	interactionEnergies = zip(complex.mutations, interactionEnergies)
	m = Core.Matrix.Matrix(rows=interactionEnergies, headers=['Mutations', 'Total'])
	f = open('results.csv', 'w+')
	f.write(m.csvRepresentation())
	f.close()

	print '\nHotspots (Delta Delta G > 1Kcal/mol):\n'
	for row in m:
		if row[1] > 4:
			print '%s: %s kj/mol' % tuple(row)

	print '\nNeutral (Delta Delta G < 1Kcal/mol):\n'
	for row in m:
		if row[1] < 4:
			print '%s: %s kj/mol' % tuple(row)
