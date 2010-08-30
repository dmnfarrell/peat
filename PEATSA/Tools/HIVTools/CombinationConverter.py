#! /usr/bin/env python
import sys
import PEAT_SA.Core as Core
import Protool
import itertools

def getPathSequence(combinations):

	path = []
	currentSet = set(combinations[0].split(','))
	path.append(combinations[0])
	for i in range(1, len(combinations)):
		newSet = set(combinations[i].split(','))
		newElement = newSet.difference(currentSet)
		path.append(list(newElement)[0])
		currentSet = newSet
	
	return path

def getTypePath(combinations, typeMap):

	path = getPathSequence(combinations)
	print 'Mutation accumulation patern ', path
	types = []
	for el in path:
		try:
			types.append(typeMap[el])
		except KeyError:
			types.append('N/A')
		
	print types

	return types
	

def codesToCombinations(mutationCodes):

	'''Converts a mutation code, which involves both chains, to a list of non-chain specific codes
	
	e.g. A84V+B84V+A32F+B32F => 84V, 32F'''

	holder = []
	for code in mutationCodes:
		mutationSet = Core.Data.MutationSet(code)
		holder.append(list(set([code[1:] for code in mutationSet.reducedMutationCodes()])))

	return holder

def observedPathsForCombination(combination, observedCombinations):

	print '\nSearching for observed paths to combination %s' % combination
	numberOfMutations = len(combination)
	pathways = itertools.permutations(combination, numberOfMutations)

	checked = 0
	found = 0
	observedPathways = []
	for pathway in pathways:
		#print 'Putative pathway %s' % list(pathway)
		parts = []
		for j in range(1,numberOfMutations + 1):
			observed = False
			sub = pathway[:j]
			#print '\tChecking subpath %s is observed' % list(sub)
			subPerms = itertools.permutations(sub)
			#Check if this sub combination 
			for subPerm in subPerms:
				subPerm = ','.join(subPerm)
				if observedCombinations.count(subPerm) == 1:
					#print '\tObserved Sub %s!' % subPerm
					parts.append(subPerm)
					observed = True
					break
			
			if observed is False:
				break

		if observed:
			found = found + 1
			observedPathways.append(parts)

		checked = checked + 1

	print '%d putative pathways. %d found\n' % (checked, found)

	return observedPathways

def vitalityProfileForPath(path, vitalities, fold):

	print 'Vitalities :',
	values = []
	for combination in path:
		print vitalities[combination],
		values.append(vitalities[combination])
	
	print '\n',
	print 'Foldn :',
	folds = []
	for combination in path:
		print fold[combination],
		folds.append(fold[combination])

	print '\n'
	
	return values, folds

		
#Read in types
typeData = Core.Matrix.matrixFromCSVFile(sys.argv[2])
typeIndex = typeData.indexOfColumnWithHeader('Type')

#Get all entries for specified drug
drugName = sys.argv[4]
trimMatrix = Core.Matrix.PEATSAMatrix(rows=[[0]*9], headers=typeData.columnHeaders())
drugNameIndex = typeData.indexOfColumnWithHeader('Drug Name')
for row in typeData:
	if row[drugNameIndex] == drugName:
		trimMatrix.addRow(row)

#Read in combinations 
combinationData = Core.Matrix.matrixFromCSVFile(sys.argv[1])
mutationCodes = combinationData.column(0)
combinations = codesToCombinations(mutationCodes)
print combinations
vitalities = combinationData.columnWithHeader(drugName+'Vitality')
fold = combinationData.columnWithHeader(drugName+'Fold')

pdb = Protool.structureIO()
pdb.readpdb(sys.argv[3])
types = []

combinationStrings = [','.join(combo) for combo in combinations]
#Skip WT
mutations = trimMatrix.columnWithHeader('Mutations')[1:]
mutations = [(el[:-2] + el[-1]) for el in mutations]
typeMap = dict(zip(mutations, trimMatrix.column(typeIndex)))
filteredPaths = []
for  combination in combinations:
	paths = observedPathsForCombination(combination, combinationStrings)
	for path in paths:
		accumulationPattern = getPathSequence(path) 
		if accumulationPattern[-1][:2] == '46' and len(accumulationPattern) > 1:
			print 'Found paths ending with mutation to 46'
			filteredPaths.append(path)

results = []
for path in filteredPaths:
	#typePath = getTypePath(path, typeMap)
	profile, foldres = vitalityProfileForPath(path, dict(zip(combinationStrings, vitalities)), dict(zip(combinationStrings, fold)))
	if profile[-2:].count('') == 0:
		mutation = path[-2]
		entry = [mutation, profile[-1] - profile[-2], foldres[-1]/foldres[-2]]
		if results.count(entry) == 0:
			results.append(entry)
	else:
		print 'Skipping - Missing data\n'

for m in results:
	print m



