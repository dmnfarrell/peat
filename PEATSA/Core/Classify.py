#! /bin/env python
import numpy, operator, os, sys
import mvpa.datasets
import mvpa.measures.noiseperturbation
import mvpa.clfs.lars
import PEATSA.Core as Core

class StructureClassifier(object):
    """Classification algorithm for analysis of structural mutations"""
    def __init__(self):
        return
        
    
    def createDataSet(self, matrix, features, recordRetriever, responseVariable):
    
        '''Creates a mvpva.dataset.Dataset instance from the provided data
    
        Parameters:
            matrix - A matrix where each row corresponds to a mutant.
                The mutant is identified by the column labeled 'Mutations'
                One of the other columns of the matrix must contain the data on the response variable to be examined
            features -
                A list containing a set of features. A feature is some property of each mutant.
                The point of the classifications is to see which features influence the response variable
            recordRetreiver -
                A function that take a feature and returns all mutants that have that feature.
            responseVariable -
                A string which identifies the column in matrix containing the response variable to tbe examined'''
    
        #Create the dataset matrix
        rows = []
        mutations = list(matrix.mutations)
        for i in range(len(mutations)):
            rows.append(len(features)*[0])
    
        mutationColumnIndex = matrix.indexOfColumnWithHeader('Mutations')
        for i in range(len(features)):
            #Get all the mutants that have feature
            feature = features[i]
            if type(feature) != list:
                feature = [feature]
    
            matching = recordRetriever(feature)
            for row in matching:
                mutation = row[mutationColumnIndex]
                index = mutations.index(mutation)
                rows[index][i] = 1
    
        responses = matrix.__getattr__(responseVariable)
        responses = [abs(response) for response in responses]
    
        #Create the classifier data set
        matrix = Core.Matrix.Matrix(rows=rows, headers=features)
        stream = open('.temp.csv', 'w+')
        stream.write(matrix.csvRepresentation())
        stream.close()
    
        matrix = numpy.loadtxt('.temp.csv', delimiter=',', skiprows=1)
        data = mvpa.datasets.Dataset(samples=matrix, labels=responses)
    
        os.remove('.temp.csv')
    
        return data
        

    def compareResidues(self, x, y):    
        retval = cmp(x[1], y[1])
        if retval == 0:
            retval = cmp(x[2], y[2])    
        return -1*retval
        
 
    def runAnalysis(self, data, features):
    
        '''Runs the LARS algorithm on data and returns features strongly associated with the response variable
    
        Parameters
            data - A mvpa.datasets.Dataset instance.
            features - A list of features being used for the classification process.
                This is not necessary for the classification process since data already
                contains the information.
                However it is necessary to map back the result of the classification process to the features'''
    
        classifier = mvpa.clfs.lars.LARS(model_type='lar', trace=False)
        classifier.train(data)
    
        sen = classifier.getSensitivityAnalyzer(force_training=True)
        perturber = mvpa.measures.noiseperturbation.NoisePerturbationSensitivity(sen)
        pertuberResults = perturber(data)
    
        result = []
        errorCutoff = 1
        setCutoff = 3
        found = []
        systematicGroups = []
        outliers = []
        for i in range(len(features)):
            vector = data.samples[:,i]
            sum = reduce(operator.add, vector)
            weight = classifier.weights[i]
            entry = [features[i], classifier.weights[i], sum]
            if weight != 0:
                print entry, " ", pertuberResults[i]
                found.append(entry)
            if sum >= setCutoff and weight >= errorCutoff:
                systematicGroups.append(entry)
            elif sum <= setCutoff and weight >= errorCutoff:
                outliers.append(entry)
    
        systematicGroups.sort(self.compareResidues)
        outliers.sort(self.compareResidues)
    
        if len(systematicGroups) is not 0:
            print 'Features possible associated with systematic error'
            for res in systematicGroups:
                print res
        else:
            print 'No systematic groups identified'
    
        if len(outliers) is not 0:
            print 'Features associated with possible outliers'
            for res in outliers:
                print res
        else:
            print 'No outlier groups identified'
    
        return found, systematicGroups, outliers
        

    def intraGroupDependance(self, matrix, group, featureFinder, allFeatures, responseVariable):
    
        '''Compares the error of entries in matrix that have features in groups
        
        Params:
            matrix - The samples - Each sample is a mutation with a value for the responseVariabe
            group - A list of features e.g. residue positions, substitutions. This is a subset of allFeatures
            featureFinder - A method that when past a feature will return a list of entries in matrix with that feature
            allFeatures - The complete list of features (which group is a subset of)
        '''
            
    
        features = group[0]
        print 'Group - %s' % features
        
        if len(features) == 0:
            print 'Only one feature in group - cannot perform intra-group dependence'
            return
        
        print '%10s\t%10s\t%10s\t%10s' % ('Feature One', 'Feature Two', 'Samples', 'Mean')
        #Iterate over the group features
        for feature in features:
            #Get all data that has that feature
            data = getattr(matrix, featureFinder)([feature])
            for featureTwo in features:
                if featureTwo == feature:
                    continue
                    
                #Get the data that has both features
                submatrix = getattr(data, featureFinder)([featureTwo])
                if submatrix is not None:
                    error = [abs(el) for el in submatrix.columnWithHeader(responseVariable)]
                    meanError = reduce(operator.add, error)/float(len(error))
                    print '%10s\t%10s\t%10d\t%10.3lf' % (feature, featureTwo, len(error), meanError)
    
            #Now find mean of mutants with this mutation but NO other systematic mutation
            check = []
            check.extend(features)
            check.pop(check.index(feature))
            subset = [code for code in allFeatures if code not in check]
    
            submatrix = getattr(matrix, featureFinder)([feature], matchMode='All')
            submatrix = getattr(submatrix, featureFinder)(subset, matchMode='Exclusive')
            if submatrix is not None:
                error = [abs(el) for el in submatrix.columnWithHeader(responseVariable)]
                meanError = reduce(operator.add, error)/float(len(error))
                print '%10s\t%10s\t%10d\t%10.3lf' % (feature, 'Extra-Group', len(error), meanError)
                
    

    def createGroups(self, matrix, allFeatures, systematicFeatures, outlierFeatures, featureFinder, multi=False):
    
        for feature in systematicFeatures:
            allFeatures.remove(feature)
    
        for feature in outlierFeatures:
            allFeatures.remove(feature)
    
        groups = {}
        #If the indivifaul features do not involve multiple positions or multiple subs its straight forward
        if not multi:
            groups['good'] = [allFeatures, getattr(matrix, featureFinder)(allFeatures, matchMode="Exclusive")]
            groups['systematic'] = [systematicFeatures, getattr(matrix, featureFinder)(systematicFeatures, matchMode="Any")]
            groups['outlier'] = [outlierFeatures, getattr(matrix, featureFinder)(outlierFeatures, matchMode="Any")]
            return groups
    
        #If the individual features are e.g. double mutants its more difficult
        #We then have to search for entries matching each feature one at a time
    
        rows = []
        if len(systematicFeatures) is not 0:
            for feature in systematicFeatures:
                if type(feature) is not list:
                    feature = [feature]
                data = getattr(matrix, featureFinder)(feature, matchMode="All")
                rows.extend(data.mutations)
    
            #Different features may match the same row - have to make a non-redundant set
            rows = list(set(rows))
            rows = [matrix.dataForMutation(mutation) for mutation in rows]
    
            groups['systematic'] = [systematicFeatures, Core.Matrix.PEATSAMatrix(rows=rows, headers=matrix.columnHeaders())]
            systematicEntries = groups['systematic'][1].mutations
        else:
            groups['systematic'] = [systematicFeatures, None]
            systematicEntries = []
    
        rows[:] = []
        if len(outlierFeatures) is not 0:
            for feature in outlierFeatures:
                if type(feature) is not list:
                    feature = [feature]
                data = getattr(matrix, featureFinder)(feature, matchMode="All")
                rows.extend(data.mutations)
    
            rows = list(set(rows))
            rows = [matrix.dataForMutation(mutation) for mutation in rows]
            groups['outlier'] = [outlierFeatures, Core.Matrix.PEATSAMatrix(rows=rows, headers=matrix.columnHeaders())]
            outlierEntries = groups['outlier'][1].mutations
        else:
            groups['outlier'] = [outlierFeatures, None]
            outlierEntries = []
    
        rows[:] = []
        mutations = matrix.mutations
        for i in range(len(mutations)):
            if mutations[i] not in systematicEntries and mutations[i] not in outlierEntries:
                rows.append(matrix.row(i))
    
        groups['good'] = [allFeatures, Core.Matrix.PEATSAMatrix(rows=rows, headers=matrix.columnHeaders())]
    
        return groups
        

    def commonElements(self, matrixOne, matrixTwo):
    
        mutantsOne = matrixOne.mutations
        mutantsTwo = matrixTwo.mutations
    
        rows = []
        for i in range(len(mutantsTwo)):
            if mutantsTwo[i] in mutantsOne:
                rows.append(matrixTwo.row(i))
    
        result = None
        if len(rows) is not 0:
            result = Core.Matrix.PEATSAMatrix(rows, headers=matrixTwo.columnHeaders())
    
        return result
    

    def compareGroups(self, groupOne, groupTwo, idOne, idTwo, featureFinder, responseVariable):
    
            if groupOne['systematic'][1]==None:
                    return
            print '%10s\t%10s\t%10s\t%10s' % (idOne, idTwo, 'Samples', 'Mean')
            for keyOne in groupOne.keys():
                    data = groupOne[keyOne][1]
                    if data is None:
                            continue
    
                    for keyTwo in groupTwo.keys():
                            if groupTwo[keyTwo][1] is not None:
                                    features = groupTwo[keyTwo][0]
                                    #subdata = getattr(data, featureFinder)(features, matchMode='Any')
                                    subdata = self.commonElements(data, groupTwo[keyTwo][1])
    
                                    if subdata is not None:
                                            elements = subdata.numberOfRows()
                                            mean = reduce(operator.add, [abs(el) for el in subdata.columnWithHeader(responseVariable)])/elements
                                    else:
                                            mean = 0
                                            elements = 0
    
                                    if elements != 0:
                                            print '%10s\t%10s\t%10d\t%10.3lf' % (keyOne, keyTwo, elements, mean)
                                    else:
                                            print '%10s\t%10s\tNo Samples ...' % (keyOne, keyTwo)


    def doRun(self, matrix, responseVariable):
        """Run standard features - position, mutations, substitutions"""
        print matrix.columnHeaders()
        errorIndex = matrix.indexOfColumnWithHeader(responseVariable)
       
        print 'Positions \n'
        features = matrix.mutatedResidues()
        data = self.createDataSet(matrix, features, matrix.entriesWithMutatedResidues, responseVariable)
        allPosition, systematicPosition, outliersPosition = self.runAnalysis(data, features)
    
        '''print '\nSubstitutions\n'
        features = matrix.substitutions()
        data = self.createDataSet(matrix, features, matrix.entriesWithSubstitutions, responseVariable)
        systematicSubstitution, outliersSubstitution = self.runAnalysis(data, features)
    
        print '\nPaired Positions\n'
        features = matrix.pairedMutatedResidues()
        data = self.createDataSet(matrix, features, matrix.entriesWithMutatedResidues, responseVariable)
        systematicDualPos, outliersDualPos = self.runAnalysis(data, features)
    
        print '\nMutations\n'
        features = matrix.mutationCodes()
        data = self.createDataSet(matrix, features, matrix.entriesWithMutations, responseVariable)
        systematicMuts, outliersMuts = self.runAnalysis(data, features)'''
    
        #subsitutionPositionDependance(systematicPosition, systematicSubstitution, matrix)
        #subsitutionPositionDependance(outliersPosition, systematicSubstitution, matrix)
    
        posGroups = self.createGroups(matrix,
                matrix.mutatedResidues(),
                [pos[0] for pos in systematicPosition],
                [pos[0] for pos in outliersPosition],
                'entriesWithMutatedResidues')
    
        '''subGroups = self.createGroups(matrix,
                matrix.substitutions(),
                [pos[0] for pos in systematicSubstitution],
                [pos[0] for pos in outliersSubstitution],
                'entriesWithSubstitutions')
    
        dualGroups = self.createGroups(matrix,
                matrix.pairedMutatedResidues(),
                [pos[0] for pos in systematicDualPos],
                [pos[0] for pos in outliersDualPos],
                'entriesWithMutatedResidues',
                multi=True)
    
        mutGroups = self.createGroups(matrix, 
                matrix.mutationCodes(), 
                [pos[0] for pos in systematicMuts], 
                [pos[0] for pos in outliersMuts], 
                'entriesWithMutations',
                multi=True)'''
    
        print '\nIntra-Group\n'
    
        self.intraGroupDependance(matrix, posGroups['systematic'], 'entriesWithMutatedResidues',
                                     matrix.mutatedResidues(), responseVariable)
        print '\n'
        '''self.intraGroupDependance(matrix, subGroups['systematic'], 'entriesWithSubstitutions',
                                    matrix.substitutions(), responseVariable)
        print '\n'
        self.intraGroupDependance(matrix, mutGroups['systematic'], 'entriesWithMutations', matrix.mutationCodes(), responseVariable)
        
        print '\nCo-dependance\n'
        self.compareGroups(posGroups, subGroups, 'Position', 'Substitution', 'entriesWithSubstitutions', responseVariable)
        print '\n'
        self.compareGroups(posGroups, mutGroups, 'Position', 'Mutations', 'entriesWithMutations', responseVariable)
        print '\n'
        self.compareGroups(posGroups, posGroups, 'Position', 'Position', 'entriesWithMutatedResidues', responseVariable)
        print '\n'
        self.compareGroups(subGroups, subGroups, 'Substitution', 'Substitution', 'entriesWithSubstitutions', responseVariable)
        print '\n'
        self.compareGroups(dualGroups, posGroups, 'Double Pos', 'Single Pos', 'entriesWithSubstitutions', responseVariable)
        print '\n'
        self.compareGroups(dualGroups, subGroups, 'Double Pos', 'Subs', 'entriesWithSubstitutions', responseVariable)
        '''        
        return allPosition, systematicPosition, outliersPosition
        
