#! /bin/env python
import os, math, operator
import optparse, sys
import PEAT_SA.Core as Core

#In KJ per mol kelvin
GasConstant = 8.314472E-3

def VitalitiesForDrug(drugData, substrateData, combination='Average', verbose=True):	
		
	print '\nCalculating vitalities for %s' % drugData
	print '\nUsing %s of substrate values' % combination
		
	#Get Drug and substrate binding - The drug data is a list of (mutationCode, bindingValue) pairs
	#Make MutationCode:Value dictionaries for the substrates
	drugBindingData = zip(drugData.deltaBindingResults.mutations, drugData.deltaBindingResults.total)
	substrateBindingData = [dict(zip(substrate.deltaBindingResults.mutations, substrate.deltaBindingResults.total)) for substrate in substrateData]

	vitalities = []
	missingData = 0
	substrateUsed = []
	for mutationCode, drugValue in drugBindingData:
		#Find substrate with maximum delta delta binding

		try:
			substrateValues = [bindingData[mutationCode] for bindingData in substrateBindingData]
		except:
			if verbose:
				print 'Unable to get substrate value for %s' % mutationCode
			missingData = missingData + 1
			continue
		#if verbose:
		#	print '\t      Values %s' % (substrateValues)

		if combination == 'Average':
			substrateValue = reduce(operator.add, substrateValues)/len(substrateValues)
			maxSub = 'None'
		elif combination == 'Max':
			substrateValue = max(substrateValues)
			maxSub = substrateData[substrateValues.index(substrateValue)].name
		elif combination == 'Min':
			substrateValue = min(substrateValues)
			maxSub = substrateData[substrateValues.index(substrateValue)].name
		elif combination == 'Boltzmann':
			#print '%s' % substrateValues
			maxSub = 'None'
			weights = [math.exp(-1*value/(GasConstant*300)) for value in substrateValues]
			partition = reduce(operator.add, weights)
			weights = [value/partition for value in weights]
			#print 'Wieghts %s' % newValues 
			newValues = map(lambda x,y: x*y, weights, substrateValues)
			#print 'Wieghted %s' % newValues 
			substrateValue = reduce(operator.add, newValues)
			#print 'Final value %lf' % substrateValue

		#Calculate the vitality
		vitality = (drugValue - substrateValue)
		#vitality = math.exp(vitality)
		#Hack since mutations aren't in correct format - using _ instead of +
		mutationCode = mutationCode.replace('_', '+')
			
		vitalities.append((mutationCode, vitality))
		if verbose:
			print '\t%s\t%10.3lf\t%10.3lf (%s) \t%10.3lf' % (mutationCode, drugValue, substrateValue, maxSub, vitality)
	
	print 'Unable to retrieve substrate data for %d mutations' % missingData
	
	if verbose:
		sortedVitalities = []
		sortedVitalities.extend(vitalities)
		sortedVitalities.sort(cmp, lambda x: x[1], True)
		for data in sortedVitalities:
			print '%s %lf' % data
	
	vitalities = Core.Matrix.PEATSAMatrix(rows=vitalities, headers=['Mutations', 'Vitality'])	
		
	return vitalities

def parseCommandLine():

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-d", "--drug", dest="drugResults",
			  help="The drug results directory.", metavar="DRUG")
	parser.add_option("-s", "--substrate",
			   dest="substrateResults",
			  help="A comman delimited list of the substrates to use",
			  metavar='SUBSTRATES')
	parser.add_option("-v", "--verbose", action="store_true",
			   dest="verbose", default=False,
			  help="Verbosity")
	parser.add_option("-o", "--outputFile", 
			   dest="outputFile", default='Vitalities.csv',
			  help="Where to put the results")
									  
	(options, args) = parser.parse_args()

	return options
	
if __name__ == "__main__":

	options = parseCommandLine()

	if options.drugResults is None:
		print 'Bound configuration must be provided'
		sys.exit(1)
	else:
		split = os.path.split(options.drugResults)
		drugData = Core.Data.DataSet(name=split[1], location=split[0])
			  
	if options.substrateResults is None:
		print 'Substrate results must be provided'
		sys.exit(1)
	else:
		directories = options.substrateResults.split(',')
		directories = [dir.strip() for dir in directories]
		data = [os.path.split(dir) for dir in directories]
		substrateData = [Core.Data.DataSet(name=dir[1], location=dir[0]) for dir in data]

	#The calculation
	result = VitalitiesForDrug(drugData, substrateData, combination='Average', verbose=options.verbose)

	#Now get Max and Min vitalties
	vitalities = VitalitiesForDrug(drugData, substrateData, combination='Max', verbose=options.verbose)
	result.addColumn(vitalities.column(1))
	result.setColumnHeader(2, 'Max')

	vitalities = VitalitiesForDrug(drugData, substrateData, combination='Min', verbose=options.verbose)
	result.addColumn(vitalities.column(1))
	result.setColumnHeader(3, 'Min')

	vitalities = VitalitiesForDrug(drugData, substrateData, combination='Boltzmann', verbose=options.verbose)
	result.addColumn(vitalities.column(1))
	result.setColumnHeader(4, 'Boltzmann')
				
	stream = open(options.outputFile, 'w+')
	stream.write(result.csvRepresentation())
	stream.close()
	print 'Wrote results to %s' % options.outputFile
	
	#Get the stabilities for this drug
	stabilityOut = drugData.stabilityResults.filterColumnByValue(columnHeader='Total', value=10)
	print "%d mutations destabilise protein by more then 10 KJ/mol" % len(stabilityOut)	
	destabilisingMutations = [element[0] for element in stabilityOut]
	for mutation in destabilisingMutations:
		print mutation
		
		
