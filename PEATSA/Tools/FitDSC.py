#! /usr/bin/env python
'''Fits data from a dsc experiment. The data must be in a csv file with the temperatures 
in one column and the differential heat-capacity in the other. The heat-capacity is expected 
to be in mols/K'''

import optparse, operator
import PEATSA.Core as Core
from PEATDB.Ekin.Fitting import Fitting
import numpy, math
import scipy.interpolate
import scipy.integrate

class DSCFitter:

	def _chopData(self, start, end):

		'''Returns a dict containing the temperature and value data for the given region'''

		aDict = {}
		aDict['temperatures'] =	self.temperatures[start:end]
		aDict['data'] =	self.data[start:end]

		return aDict

	def _indexesForTemperatures(self, data, temperatures):

		'''Finds the indexes of the elements in temperatures which are closes to the values in data

		Parameters:
			data: A tuple of temperatures
			temperatures: A list of temperatures

		Returns:
			A tuple containing the indexes of the elements in temperatures that are closest to the va'''
		

		indexes = []
		for el in data:
			for i in range(len(temperatures)):
				if temperatures[i] > el:
					#print el, temperatures[i], temperatures[i-1]
					indexes.append(i - 1)
					break

			if i == len(temperatures) - 1:
				indexes.append(i)
		
		if len(indexes) == 3:
			indexes.pop()
		
		return indexes

	def _normaliseData(self):

		'''Performs some transformations on the dsc data

		Makes the lowest point on the curve 0 - this is for the integration
		Converts Heat Capacity (cal/K) to specific heat-capacity (kcal.mol/K)
		Convertis celsius to Kelvin'''

		if self.verbose is True:
			print 'Converting cal/K to kcal/(K.mol)'
			print 'Converting temperature to Kelvin'

		minValue = min(self.data)
		self.data = [el - minValue for el in self.data]
		self.data = [el/(1000*self.mol) for el in self.data]
		self.temperatures = [temp + 273.15 for temp in self.temperatures]

	def __init__(self, matrix, temperatureColumn, dataColumn, foldedRange, unfoldedRange, mol, verbose=False):

		'''Returns a DSCFitter instance initialised to fit the data in matrix

		Note the DSC measurements are assumed to be in calories/kelvin and the temperature to be in degrees

		Parameters:
			matrix - A Core.Matrix.Matrix instance containing the DSC data
			temperatureColumn - The column in matrix containing the temperature data
			dataColumn - The column in matrix containing the heat-capacity measurements.
			foldedRange - The range of temperatures to use to fit the folded region
			unfoldedRange - The range of temperatures to use to fit the unfolded region
			mol - The amount of the substance being scanned in mols.
		'''

		self.verbose = verbose
		self.mol = mol
		self.temperatures = matrix.columnWithHeader(temperatureColumn)
		self.temperatures = [float(el) for el in self.temperatures]
		self.data = matrix.columnWithHeader(dataColumn)
		self.data = [float(el) for el in self.data]

		self._normaliseData()
		self.setRanges(foldedRange, unfoldedRange)
		self.useFolded = False

	def __str__(self):

		string = 'DSC Fitter.\nCurve Regions:\n'
		string = string + ("\tFolded region %5.2lf - %5.2lf Kelvin (Rows %d to %d). Width: %-5.2lf K\n" % 
				(self.foldedRange[0], self.foldedRange[1], self.foldedIndexes[0], 
				self.foldedIndexes[1], self.foldedRange[1] - self.foldedRange[0]))
		string = string + ("\tTranition region %5.2lf - %5.2lf Kelvin. (Rows %d to %d). Width: %-5.2lf K\n" % 
				(self.foldedRange[1], self.unfoldedRange[0], self.foldedIndexes[1], self.unfoldedIndexes[0], 
				self.unfoldedRange[0] - self.foldedRange[1]))
		string = string + ("\tUnfolded region %5.2lf - %5.2lf Kelvin. (Rows %d to %d). Width: %-5.2lf K\n" % 
				(self.unfoldedRange[0], self.unfoldedRange[1], self.unfoldedIndexes[0], self.unfoldedIndexes[1],
				self.unfoldedRange[1] - self.unfoldedRange[0]))

		if self.baselineFitted:
			string = string + '\nFitt Data:\n'	
			string = string + '\tFolded Region: Slope: %4.3E Intercept: %4.3E\n'	% (self.foldedSlope, self.foldedIntercept)
			string = string + '\tUnfolded Region: Slope: %4.3E Intercept: %4.3E\n' % (self.unfoldedSlope, self.unfoldedIntercept)
			string = string + '\tFolded Offset: %4.3E kcal/(mol.K)\n' % self.foldedOffset
			string = string + '\tDelta Specific Heat: %4.3E kcal/(mol.K)\n' % self.deltaHC

		return string

	def setRanges(self, foldedRange, unfoldedRange):

		'''Sets the ranges in the curve corresponding to folded, unfolded and transition regions

		Parameters:
			foldedRange: A list or tuple containin the indexes of the row in the matrix
			corresponding to the start and end of the folded segment
			unfoldedRange: As above except for the unfolded segment'''

		foldedRange = [temp + 273.15 for temp in foldedRange]
		unfoldedRange = [temp + 273.15 for temp in unfoldedRange]

		foldedIndexes = self._indexesForTemperatures(foldedRange, self.temperatures)
		unfoldedIndexes = self._indexesForTemperatures(unfoldedRange, self.temperatures)
		transitionIndexes = (foldedIndexes[1], unfoldedIndexes[0])

		if(unfoldedIndexes[0] == unfoldedIndexes[1]):
			raise Exception, 'Not enough data in unfolded region to fit'

		self.foldedRange = foldedRange
		self.unfoldedRange = unfoldedRange
		self.foldedIndexes = foldedIndexes
		self.unfoldedIndexes = unfoldedIndexes
		self.transitionIndexes = transitionIndexes

		self.foldedRegion = self._chopData(self.foldedIndexes[0], self.foldedIndexes[1])
		self.unfoldedRegion = self._chopData(self.unfoldedIndexes[0], self.unfoldedIndexes[1])
		self.transitionRegion = self._chopData(self.foldedIndexes[1], self.unfoldedIndexes[0])

		self.baselineFitted = False
		self.linearFitted = False

	def progressBaseline(self):

		'''Calculates the progress baseline for the dsc peak defined by the matrix data

		Parameters:
			deltaCp - Estimate of the heat-capacity different between the folded and unfolded states
			offset - The value of the heat-capacity just before the transition'''

		#The linear segements must be fitted to get the progress baseline
		if self.linearFitted == False:
			self._fitLinear()

		temperatures = self.transitionRegion['temperatures']
		data = self.transitionRegion['data']
			
		integrals = []
		spline = scipy.interpolate.UnivariateSpline(temperatures, data) 	
		initialTemp = temperatures[0]
		for temp in temperatures[1:]:
			integrals.append(spline.integral(initialTemp, temp))

		integrals.insert(0,0)
		totalIntegral = integrals[-1]
		fractions = [el/totalIntegral for el in integrals]
		baseline = []
		for i in range(len(fractions)):
			fraction = fractions[i]
			temperature = temperatures[i]
			#Cp(t) = Cn(T) + fraction*deltaCp - where Cn(T) is the extrapoloated HC of the folded state
			#print 'Fraction %lf D-SHC %lf Slope %lf Temp %lf Offset %lf' % (fraction, self.deltaHC, self.foldedSlope, temperature, self.foldedOffset)
			value = fraction*self.deltaHC + self.foldedSlope*temperature + self.foldedIntercept
			baseline.append(value)

		lnk = [(math.log(fraction) - math.log(1-fraction)) for fraction in fractions[1:-1]]
		lnk.insert(0,'NA')
		lnk.append('NA')

		#The progress baseline we have the fraction unfolded protein
		#Kunfold = f/(1-f) 
		rows = zip(temperatures, integrals, fractions, baseline, lnk)
		return Core.Matrix.Matrix(rows=rows, headers=['Temperature', 'Integral', 'AreaFraction', 'Baseline', 'lnku'])

	def _fitLinear(self):

		'''Private method: Fits the linear segments of the data'''

		if self.linearFitted is False:
			expdata = zip(self.foldedRegion['temperatures'], self.foldedRegion['data'])
			d,f = Fitting.doFit(model='Linear', conv=1E-11, expdata=expdata, silent=True)
			self.foldedSlope = d[0]
			self.foldedIntercept = d[1]

			self.foldedOffset = self.foldedSlope*self.temperatures[self.foldedIndexes[1]] + self.foldedIntercept

			if self.useFolded is True or (len(self.unfoldedRegion['temperatures']) < 2):
				self.unfoldedSlope = self.foldedSlope
				unfoldedHC = self.data[self.unfoldedIndexes[0]]
				foldedHC = self.foldedSlope*self.temperatures[self.unfoldedIndexes[0]] + self.foldedIntercept
				self.unfoldedIntercept = self.foldedIntercept + unfoldedHC - foldedHC
				#print 'Unfolded intercept derived from Folded ', self.unfoldedIntercept
			else:
				expdata = zip(self.unfoldedRegion['temperatures'], self.unfoldedRegion['data'])
				d,f = Fitting.doFit(model='Linear', conv=1E-11, expdata=expdata, silent=True)
				self.unfoldedSlope = d[0]
				self.unfoldedIntercept = d[1]
				#print 'Unfolded intercept derived from fit ',self.unfoldedIntercept

			#Calculate the heat capacity difference at the starting of the unfolded region
			unfoldedHC = self.data[self.unfoldedIndexes[0]]
			foldedHC = self.foldedSlope*self.temperatures[self.unfoldedIndexes[0]] + self.foldedIntercept
			self.deltaHC = unfoldedHC - foldedHC

			self.linearFitted = True

	def curve(self):

		rows = zip(self.temperatures, self.data)
		return Core.Matrix.Matrix(rows=rows, headers=['Temperature', 'Value'])

	def baseline(self):

		'''Returns a matrix representing the baseline calculated for the data.'''

		if self.baselineFitted is False:

			#The linear fit may already be done if progressBaseline was called previously
			if self.linearFitted is False:
				self._fitLinear()

			baselineValues = [(self.foldedSlope*temp + self.foldedIntercept) for temp in self.foldedRegion['temperatures']]
			totalBaselineRows = zip(self.foldedRegion['temperatures'], baselineValues)

			progressBaseline = self.progressBaseline()
			baselineValues = progressBaseline.columnWithHeader('Baseline')
			totalBaselineRows.extend(zip(self.transitionRegion['temperatures'], baselineValues))

			baselineValues = [(self.unfoldedSlope*temp + self.unfoldedIntercept) for temp in self.unfoldedRegion['temperatures']]
			totalBaselineRows.extend(zip(self.unfoldedRegion['temperatures'], baselineValues))

			self.baselineMatrix = Core.Matrix.Matrix(rows=totalBaselineRows, headers=['Temperature', 'Value'])
			self.baselineFitted = True

		return self.baselineMatrix

	def normalisedCurve(self):

		'''Normalises the dsc data by subtracting the baseline.

		This essentially makes the curve independant of specific heat capacities of the states

		Returns:
			A matrix containing the normalised curve'''
		
		#Subtract the baseline segments
		baseline = self.baseline()
		baselineValues = baseline.columnWithHeader('Value')

		#The baseline may not start at the first point in the graph
		#Have to take this into account
		offset = self.foldedIndexes[0]
		start, end = self.foldedIndexes
		start = start - offset
		end = end - offset
		foldedSegment = map(operator.sub, self.foldedRegion['data'], baselineValues[start:end]) 
		dscCurveRows = zip(self.foldedRegion['temperatures'], foldedSegment)

		#Subtract the progress baseline
		start, end = self.transitionIndexes
		start = start - offset
		end = end - offset
		transitionSegment = map(operator.sub, self.transitionRegion['data'], baselineValues[start:end]) 
		dscCurveRows.extend(zip(self.transitionRegion['temperatures'], transitionSegment))

		#Subtract the unfolded baseline
		start, end = self.unfoldedIndexes
		start = start - offset
		end = end - offset
		unfoldedSegment = map(operator.sub, self.unfoldedRegion['data'], baselineValues[start:end]) 
		dscCurveRows.extend(zip(self.unfoldedRegion['temperatures'], unfoldedSegment))

		return Core.Matrix.Matrix(rows=dscCurveRows, headers=['Temperature', 'Value'])

	def fit(self):

		'''Fits the normalised cureve to a two-state and non-two-state dsc model'''

		if len(self.unfoldedRegion['temperatures']) < 2 or self.useFolded is True:
			method = 'FoldedOnly'			
		else:
			method = 'Both'

		dscCurve = self.normalisedCurve()
		value = [el*4.184 for el in dscCurve.value]
		
		fittingParameters,fittingInstance=Fitting.doFit(expdata=zip(dscCurve.temperature,value),model='DSC2state', silent=True)

		vantHoff = fittingParameters[1]/4.184
		meltingTemp = fittingParameters[0]
		res1 = ['TwoState', method, self.foldedRange[1], self.unfoldedRange[0], 
			meltingTemp, vantHoff, 'N/A', fittingParameters['error']/(4.184*4.184), vantHoff/meltingTemp, '1.0']

		fittingParameters,fittingInstance=Fitting.doFit(expdata=zip(dscCurve.temperature,value), model='DSCindependent', silent=True)

		vantHoff = fittingParameters[2]/4.184
		calorimetric = fittingParameters[0]*vantHoff
		meltingTemp = fittingParameters[1]
		res2 = ['NonTwoState', method, self.foldedRange[1], self.unfoldedRange[0], 
				meltingTemp, vantHoff, calorimetric, fittingParameters['error']/(4.184*4.184), vantHoff/meltingTemp, calorimetric/vantHoff]

		return [res1, res2]

	def fitHeaders(self):
	
		return ['Model', 'LinearMethod', 'Tstart', 'Tend', 'MeltingTemp', 'VantHoff', 'Calormetric', 'Error', 'Entropy', 'Cooperativity']

if __name__ == "__main__":

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-i", "--input", dest="file",
			  help="A csv file containing the dsc data", metavar="FILE")
	parser.add_option("-t", "--temperatureColumn", dest="tempCol", default='Temperature',
			  help="The column containing the temperature data ", metavar="TEMP")
	parser.add_option("-d", "--dataColumn", dest="dataCol", default='Value',
			  help="The column containing the dsc data", metavar="DATA")
	parser.add_option("-f", "--foldedRange", dest="foldedRange",
			  help="The range of temperatures defining the folded region", metavar="FOLD")
	parser.add_option("-u", "--unfoldedRange", dest="unfoldedRange",
			  help="The range of temperatures defining the unfolded region", metavar="END")
	parser.add_option("-c", "--check-interval", dest="checkInterval", type='float', default=2,
			  help="The fitter will modify the end folded and start unfolded temperatures within this range to find better fits", 
			  metavar="INTERVAL")
	parser.add_option("-s", "--check-step", dest="checkStep", type='float', default=1,
			  help="The fitter step-size used to generate new end folded and start unfolded temperatures values", metavar="STEP")
	parser.add_option("-m", "--molar", dest="mol", type='float',
			  help="The amount of protein in the dsc cell in mols", metavar="MOL")
	parser.add_option("", "--twiddle", action="store_true",
			  dest="twiddle", default=False,
			  help="Perturb the definitions of the folded and unfolded regions and refit.\n"
			  "The number of refits is determines by the check paramters")
	parser.add_option("", "--folded-only", action="store_true",
			  dest="foldedOnly", default=False,
			  help="Fits the unfolded region using the slope obtained from the folded region")

	(options, args) = parser.parse_args()

	m = Core.Matrix.matrixFromCSVFile(options.file)
	foldedRange = [float(el) for el in options.foldedRange.split(",")]
	unfoldedRange = [float(el) for el in options.unfoldedRange.split(",")]

	fitter = DSCFitter(m, options.tempCol, options.dataCol, foldedRange, unfoldedRange, options.mol)
	print fitter

	twoState = []
	nonTwoState = []
	if options.twiddle is True:
		end = foldedRange[1]
		foldedRange = numpy.arange(end - options.checkInterval, end + options.checkInterval, options.checkStep)

		start = unfoldedRange[1]
		unfoldedRange = numpy.arange(start - options.checkInterval, start + options.checkInterval, options.checkStep)

		for end in foldedRange:
			for start in unfoldedRange: 
				try:
					fitter.setRanges([foldedRange[0], end], [start, unfoldedRange[1]])
					print '\n'
					print fitter
					results = fitter.fit()
					twoState.append(results[0])
					nonTwoState.append(results[1])
				except Exception ,data:
					print data
	else:
		results = fitter.fit()
		twoState.append(results[0])
		nonTwoState.append(results[1])
		print fitter

	results = Core.Matrix.Matrix(rows=twoState, headers=fitter.fitHeaders())
	print results.csvRepresentation()

	results = Core.Matrix.Matrix(rows=nonTwoState, headers=fitter.fitHeaders())
	print results.csvRepresentation()

	stream = open('IntegrationData.csv' ,'w+')
	stream.write(fitter.progressBaseline().csvRepresentation())
	stream.close()

	stream = open('Baseline.csv' ,'w+')
	stream.write(fitter.baseline().csvRepresentation())
	stream.close()

	stream = open('NormalisedCurve.csv' ,'w+')
	stream.write(fitter.normalisedCurve().csvRepresentation())
	stream.close()

	stream = open('Curve.csv' ,'w+')
	stream.write(fitter.curve().csvRepresentation())
	stream.close()

