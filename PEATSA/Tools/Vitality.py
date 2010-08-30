import os, math
import PEATSA.Core as Core

#In KJ per mol kelvin
GasConstant = 8.314472E-3

def IsVitalityHighResistance(vitality):
	
	if vitality < 2.3:
		return False
	else:
		return True
		
def IsHIVScoreHighResistance(score):

	if score >= 20:
		return True
	else:
		return False
		
def VitalityToResistance(vitality):

	if vitality > 2.3:
		return 2
	elif vitality > 1.1:
		return 1
	#elif vitality < 0.3:
	#	return -1
	else:
		return 0
		
def HIVScoreToResistance(score):

	if score >= 20:
		return 2
	elif score >= 10:
		return 1
	elif score < 10:
		return 0

def VitalitiesForDrug(drugDataName, drugDataLocation, substrateDataNames, substrateDataLocation, expectedValues=None, separateChains=False, verbose=False):	
		
	print '\nCalculating vitalities for %s' % drugDataName	
		
	#Load Drug data
	drugData = Core.Data.DataSet(name=drugDataName, location=drugDataLocation)

	#Load Substrate Data
	substrateData = [Core.Data.DataSet(name=value, location=substrateDataLocation) for value in substrateDataNames]

	#Get Drug and substrate binding - The drug data is a list of (mutationCode, bindingValue) pairs
	#Make MutationCode:Value dictionaries for the substrates
	drugBindingData = zip(drugData.deltaBindingResults.mutations, drugData.deltaBindingResults.total)
	substrateBindingData = [dict(zip(substrate.deltaBindingResults.mutations, substrate.deltaBindingResults.total)) for substrate in substrateData]

	vitalities = []
	missingData = 0
	for mutationCode, drugValue in drugBindingData:
		#Find substrate with maximum delta delta binding
		try:
			substrateValues = [bindingData[mutationCode] for bindingData in substrateBindingData]
		except:
			if verbose:
				print 'Unable to get substrate value for %s' % mutationCode
			missingData = missingData + 1
			continue

		maxValue = min(substrateValues)
		maxSub = substrateDataNames[substrateValues.index(maxValue)]

		#Calculate the vitality
		vitality = (1/(GasConstant*300))*(drugValue - maxValue)
		vitality = math.exp(vitality)
		#Hack since mutations aren't in correct format - using _ instead of +
		mutationCode = mutationCode.replace('_', '+')
			
		vitalities.append((mutationCode, vitality))
		#print '%s\t%10.3lf\t%10.3lf (%s) \t%10.3lf' % (mutationCode, drugValue, maxValue, maxSub, vitality)
	
	print 'Unable to retrieve substrate data for %d mutations' % missingData
	
	#Combine the values for Chains A & B
	if separateChains:
		chainAValues = [element for element in vitalities if element[0][0] == 'A']
		chainBValues = [element for element in vitalities if element[0][0] == 'B']
		combinedVitalities = []
		for pair in zip(chainAValues, chainBValues):
			combinedValue = pair[0][1] + pair[1][1]
			#print '%s\tA: %-10.3lf\tB: %-10.3lf\tTotal: %-10.3lf' % (pair[0][0][1:], pair[0][1], pair[1][1], combinedValue)
			combinedVitalities.append((pair[0][0][1:], combinedValue))

		vitalities = combinedVitalities
	
	if verbose:
		sortedVitalities = vitalities
		sortedVitalities.sort(cmp, lambda x: x[1], True)
		print "\nSorted Vitalites:"
		for data in sortedVitalities:
			print '%s %lf' % data
	
	vitalities = Core.Matrix.PEATSAMatrix(rows=vitalities, headers=['Mutations', 'Vitality'])	
		
	return vitalities
	
if __name__ == "__main__":
	
	substrateData = ['1KJ7_N25D.peatsa', '1KJH_N25D.peatsa', '1KJF_N25D.peatsa']
	drugData = ['1HSGClean.peatsa', '1MUIClean.peatsa', '1OHRClean.peatsa', '2AQUClean.peatsa']
	pdbDrugMap = {"1MUIClean":"AB1" , "2AQUClean":"DR7", "1OHRClean":"1UN", "1HSGClean":"MK1"}
	drmData = Core.Matrix.matrixFromCSVFile('DRM.csv')
	
	results = {}
	for data in drugData:
		vitalities = VitalitiesForDrug(data, '.', substrateData, '.')
		stream = open(os.path.splitext(data)[0], "w+")
		stream.write(vitalities.csvRepresentation())
		stream.close()
		results[os.path.splitext(data)[0]] = vitalities
	
	rows = 2*[[0,0]]
	totals = Core.Matrix.Matrix(rows=rows, headers=['Non', 'High Resistance'])	
				
	for data in results.keys():
		vitalities = results[data]
		print '\n', data
		print vitalities
		print "Total mutations %d" % (vitalities.numberOfRows())
		
		resistanceLevels = vitalities.mapColumn(columnHeader="Vitality", mapFunction=VitalityToResistance)
		
		resistance = [element[1] for element in resistanceLevels]
		highResistance = filter(lambda x: x==2, resistance)
		lowResistance = filter(lambda x: x==1, resistance)
		noResistance = filter(lambda x: x==0, resistance)
		hyperSusceptible = filter(lambda x: x==-1, resistance)
		
		print "High Resistance %d (%5.3lf)" % (len(highResistance), 100*len(highResistance)/float(vitalities.numberOfRows()))
		print "Low Resistance %d (%5.3lf)" % (len(lowResistance), 100*len(lowResistance)/float(vitalities.numberOfRows()))
		print "No Resistance %d (%5.3lf)" % (len(noResistance), 100*len(noResistance)/float(vitalities.numberOfRows()))
		print "Hyper Susceptible %d (%5.3lf)" % (len(hyperSusceptible), 100*len(hyperSusceptible)/float(vitalities.numberOfRows()))
		
		#Get the stabilities for this drug
		drugData = Core.Data.DataSet(name="%s.peatsa" % data, location='.')
		stabilityOut = drugData.stabilityResults.filterColumnByValue(columnHeader='Total', 
						value=10)
		print "%d mutations destabilise protein by more then 10 KJ/mol" % len(stabilityOut)	
		
		#Hack
		stabilityOut = [[element[0].replace('_', '+'), element[1]] for element in stabilityOut]
		
		destabilisingMutations = [element[0] for element in stabilityOut]
		#print destabilisingMutations
		#	 print [element[0] for element in resistanceLevels]
		resistanceLevels = [element for element in resistanceLevels if element[0] not in destabilisingMutations]
		print "Removed %d mutations from vitality due to stability considerations" % (vitalities.numberOfRows() - len(resistanceLevels))
		resistance = [element[1] for element in resistanceLevels]
		highResistance = filter(lambda x: x==2, resistance)
		lowResistance = filter(lambda x: x==1, resistance)
	
		print "High Resistance Filtered %d (%5.3lf)" % (len(highResistance), 100*len(highResistance)/float(vitalities.numberOfRows()))
		print "Low Resistance Filtered %d (%5.3lf)" % (len(lowResistance), 100*len(lowResistance)/float(vitalities.numberOfRows()))
		
		#Get the DRMs for this drug
		#name = pdbDrugMap[data]
		#data = drmData.mapColumn(columnHeader=name, mapFunction=HIVScoreToResistance)
		
		#resistanceLevels = dict(resistanceLevels)
		#rows = 2*[[0,0]]
		#matrix = Core.Matrix.Matrix(rows=rows, headers=['Non', 'High Resistance'])
		
		#for element in data:
		#	actual = element[1]
		#	try:
		#		predicted =  resistanceLevels[element[0]]
		#	except:
		#		print 'Missing data for %s' % element[0]
		#		continue
			
		#	if actual == 2 and predicted == 2:
		#		matrix.matrix[1][1] = matrix.matrix[1][1] + 1
		#	elif actual == 2 and (predicted == 0 or predicted == 1):	
		#		matrix.matrix[0][1] = matrix.matrix[0][1] + 1
		#	elif (actual == 0 or actual == 1) and predicted == 2:	
		#		matrix.matrix[1][0] = matrix.matrix[1][0] + 1
		#	else:
		#		matrix.matrix[0][0] = matrix.matrix[0][0] + 1
		#
		#print	matrix.csvRepresentation()
		
	#	for i in range(2):
	#		for j in range(2):
	#			totals.matrix[i][j] = totals.matrix[i][j] + matrix.matrix[i][j]
		
	#print totals.csvRepresentation()	
		
		
		
		
