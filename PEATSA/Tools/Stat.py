#! /bin/env python
#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
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
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#
'''Program for compiling statistics on the fit between predicted and observed values for a set of samples.
Note the fit is assumed to be 1:1 i.e. no fitting is performed.
The program calculates the correlation and R (reduced chisquared) value of the fit, along with mean and std. dev of the error.
The percentage of outliers in the data set can be specified, or a range of outlier percentages can be tested.
For each set the program also tests if the residuals are normally distributed (required for chisquared statistic).
Under assumption that the residuals should be normally distributed this can identify where the actual outliers are.
'''

import PEATSA.Core as Core
import operator, optparse
import numpy, math, sys
import scipy.stats.stats as stats
import scipy.stats.distributions as distributions
import scipy.stats

def sampleStatistics(array):

	array = numpy.array(array)
	mean = numpy.mean(array)
	stdev = numpy.std(array)
	squared =  array*array
	sse = numpy.sum(squared)
	rmse = math.sqrt(sse/len(array))
	
	return (mean, stdev, rmse)
	
def normalityTests(array):

	#Histogram
	#print scipy.stats.describe(array)

	histogram, xedges = numpy.histogram(array, 100, new=True)
	histogram = Core.Matrix.Matrix(rows=zip(list(xedges), list(histogram)))

	cumulative = [histogram.element(0,1)]
	for i in range(1, histogram.numberOfRows()):
		cumulative.append(cumulative[i-1] + histogram.element(i,1))

	length = float(cumulative[-1])
	cdf = [el/length for el in cumulative]

	histogram.addColumn(cumulative, -1)
	histogram.setColumnHeader(-1, 'Cumulative')

	histogram.addColumn(cdf, -1)
	histogram.setColumnHeader(-1, 'CDF')

	mean = numpy.mean(array)
	stdev= numpy.std(array)
	normalised = (array - mean)/stdev	

	#Shapiro-Wilk
	#Compares the expected slope of a Q-Q plot to a least squared fit of the slope
	#If the errors are normal the slope of the Q-Q will be the standard-deviation of the errors
	#In this case we have normalised the data-set to set the deviation to one
	shapiro = scipy.stats.shapiro(normalised)
	dangostino = scipy.stats.normaltest(normalised)

	#Kurtosis and Skew Test
	#print 'Kurtosis-Test %6.2lf %E' % scipy.stats.kurtosistest(normalised)

	#Kolomogrov-Smirnov
	kolomogrov = scipy.stats.kstest(normalised, 'norm')

	#Q-Q
	orderStat = range(1, len(array) + 1)
	length = float(len(array))
	for i in range(len(array)):
		orderStat[i] = orderStat[i]/length
		
	invnormStat = [distributions.norm.ppf(x) for x in orderStat]
	normalised.sort()
	data = zip(invnormStat, normalised)
	qq = Core.Matrix.Matrix(rows=data, headers=['Theoretical', 'Experimental'])

	return shapiro + dangostino + kolomogrov + (histogram, qq)


def correlationStatistics(predictions, observations, obsError, predError):

	predictions = numpy.array(predictions)
	observations = numpy.array(observations)

	correlation = numpy.corrcoef(predictions, observations, rowvar=0)
	correlation = correlation[1][0]

	if predError != 0:
		obsError = math.sqrt(obsError*obsError + predError*predError)

	#Chi-Squared
	#Assuming the observations have an error that is normaly distributed with deviation obsError
	#Checks if the errors are actually distributed around the fit-line with this deviation
	#Technically returns the probability that the observed distribution of errors comes from the supposed distribution
	errors = observations - predictions
	squaredErrors = numpy.square(errors)
	mse = numpy.mean(squaredErrors)
	chisquared = mse*len(squaredErrors)/(math.pow(obsError, 2))
	chisquaredProb = stats.chisqprob(chisquared, len(squaredErrors) - 1)
	degreesOfFreedom = len(squaredErrors) - 1
	reducedChisquared = chisquared/degreesOfFreedom

	return (correlation, chisquared, reducedChisquared, chisquaredProb)

if __name__ == "__main__":

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-f", "--file", dest="file",
			  help="A csv file containing observations and predictions for a set of samples", metavar="FILE")
	parser.add_option("-p", "--predicted",
			  dest="predicted",
			  help="The column containing the predictions",
			  metavar='PREDICTED')		  		  		  		  		  		  
	parser.add_option("-e", "--experimental",
			  dest="experimental",
			  help="The column containing the experimental values",
			  metavar='EXPERIMENTAL')		  		  		  		  		  		  
	parser.add_option("-s", "--start",
			  dest="start",
			  default=0,
			  help="The program will calculate statistics starting with this percentage of outliers removed",
			  metavar='start')		  		  		  		  		  		  
	parser.add_option("-l", "--end",
			  dest="end",
			  default=20,
			  help="The program will calculate statistics ending with this percentage of outliers removed",
			  metavar='end')		  		  		  		  		  		  
	parser.add_option("-t", "--step",
			  dest="step",
			  default=5,
			  help="The program will test sets start + n*step, while the result is less than the end percentage ",
			  metavar='end')		  		  		  		  		  		  
	parser.add_option("", "--scale-exp",
			  dest="scaleExperimental",
			  default=1,
			  help="A scale factor that will be applied to the experimental column",
			  metavar='EXPSCALE')		  		  		  		  		  		  
	parser.add_option("", "--scale-pred",
			  dest="scalePredicted",
			  default=1,
			  help="A scale factor that will be applied to the predicted column",
			  metavar='PREDSCALE')		  		  		  		  		  		  
	parser.add_option("", "--obs-error",
			  dest="observedError",
			  default=1,
			  help="The error in the observations",
			  metavar='OBSERROR')		  		  		  		  		  		  
	parser.add_option("", "--pred-error",
			  dest="predictedError",
			  default=0,
			  help="The error in the predictions",
			  metavar='OBSERROR')		  		  		  		  		  		  
	parser.add_option("", "--qq", action="store_true",
			  dest="outputQQ", default=False,
			  help="Output QQ matrix")
	parser.add_option("", "--hist", action="store_true",
			  dest="outputHist", default=False,
			  help="Output error histogram")
	parser.add_option("", "--sets", action="store_true",
			  dest="outputSets", default=False,
			  help="Output csv files for each set tested")
	parser.add_option("", "--outliers", action="store_true",
			  dest="outputOutliers", default=False,
			  help="Output csv files containing the outliers removed from each set tested")
	parser.add_option("", "--suppressHeaders", action="store_true",
			  dest="suppressHeaders", default=False,
			  help="Don't print the headers of the results matrix to stdout")

	(options, args) = parser.parse_args()

	if not sys.stdin.isatty():
		csv = sys.stdin.read()
		m = Core.Matrix.matrixFromCSVRepresentation(csv)
	elif options.file is None:
		print 'CSV file must be provided'
		sys.exit(1)
	else:
		m = Core.Matrix.matrixFromCSVFile(options.file)

	if options.experimental is None:
		print 'Experimental column must be specified'

	if options.predicted is None:
		print 'Predicted column must be specified'

	exp = m.columnWithHeader(options.experimental)
	exp = [el*float(options.scaleExperimental) for el in exp]
	predicted = m.columnWithHeader(options.predicted)
	predicted = [el*float(options.scalePredicted) for el in predicted]

	#Calculate the errors andd sort by them
	error = map(operator.sub, predicted, exp)
	m.addColumn(error)
	m.setColumnHeader(m.numberOfColumns() - 1, 'InternalError')

	absError = map(lambda x: abs(x), error)
	m.addColumn(absError)
	m.setColumnHeader(m.numberOfColumns() - 1, 'AbsError')

	squaredError = map(lambda x: x*x, error)
	m.addColumn(squaredError)
	m.setColumnHeader(m.numberOfColumns() - 1, 'SortableSquaredError')

	m.sort(columnHeader='SortableSquaredError', descending=False)

	#Get the sorted data
	exp = m.columnWithHeader(options.experimental)
	predicted = m.columnWithHeader(options.predicted)
	error = m.internalError

	#Divisions to be tested
	percentages = range(100 - int(options.end), 101 - int(options.start), int(options.step))
	percentages = [el/100.0 for el in percentages]
	percentages.reverse()
	divisions = [int(math.ceil(m.numberOfRows()*fraction)) for fraction in percentages]

	rows = []
	for el in zip(percentages, divisions):

		percentage = el[0]
		division = el[1]

		correlation, chisquared, r, pvalue = correlationStatistics(exp[:division], predicted[:division], 
								float(options.observedError), float(options.predictedError))
		mean, stdev, rmse = sampleStatistics(error[:division])
		data = normalityTests(error[:division])
		qq = data[-1]
		histogram = data[-2]

		rows.append([percentage, division, correlation, mean, stdev, rmse, chisquared, r, 
				pvalue, data[1], data[3], data[5]])

		if options.outputQQ:
			f = open('QQ%s.csv' % percentage, 'w+')
			f.write(qq.csvRepresentation())
			f.close()

		if options.outputHist:
			f = open('ErrorHist%s.csv' % percentage, 'w+')
			f.write(histogram.csvRepresentation())
			f.close()

		if options.outputSets:
			f = open('DataSet%s.csv' % percentage, 'w+')
			f.write(m[:division].csvRepresentation())
			f.close()
		
		if options.outputOutliers and (division != m.numberOfRows()):
			f = open('Outliers%s.csv' % percentage, 'w+')
			f.write(m[division:].csvRepresentation())
			f.close()
		
	headers = ['Percentage', 'Samples', 'Correl', 'MeanError', 'StdevError', 'RMSE', 'ChiSquared', 	
			'ReducedChi', 'ChiProb', 'ShaprioProb', 'DAngostinoProb', 'KSProb']
	matrix = Core.Matrix.Matrix(rows=rows, headers=headers)
	includeHeaders = not options.suppressHeaders
	print matrix.csvRepresentation(includeHeaders=includeHeaders),

