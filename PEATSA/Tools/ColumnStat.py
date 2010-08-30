#!/bin/env python
#
#  ColumnStat.py
#  ProteinDesignTool
#
#  Created by Michael Johnston on 18/06/2010.
#  Copyright (c) 2010 UCD. All rights reserved.
#
import PEAT_SA.Core as Core
import optparse
import numpy, math, sys
import scipy.stats.stats as stats
import scipy.stats.distributions as distributions
import scipy.stats

def normalityTests(array):

	#Histogram
	#print scipy.stats.describe(array)

	histogram, xedges = numpy.histogram(array, 100, new=True)
	histogram = Core.Matrix.Matrix(rows=zip(list(xedges), list(histogram)))

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

if __name__ == "__main__":
	
	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-f", "--file", dest="file",
			  help="A csv file", metavar="FILE")
	parser.add_option("-c", "--columns",
			  dest="columns",
			  help="A comma separated list of the indexes of the columns to calculate statistics for",
			  metavar='start')
	parser.add_option("", "--qq", action="store_true",
			  dest="outputQQ", default=False,
			  help="Output QQ matrix for each column")
	parser.add_option("", "--hist", action="store_true",
			  dest="outputHist", default=False,
			  help="Output histogram for each column")		  

	(options, args) = parser.parse_args()

	if not sys.stdin.isatty():
		csv = sys.stdin.read()
		m = Core.Matrix.matrixFromCSVRepresentation(csv)
	elif options.file is None:
		print 'CSV file must be provided'
		sys.exit(1)
	else:
		m = Core.Matrix.matrixFromCSVFile(options.file)

	if options.columns is None:
		print >>sys.stderr, 'No columns specified - defaulting to all'
		indexes = range(m.numberOfColumns())
	else:			
		indexes = [int(el.strip()) for el in options.columns.split(',')]	

	rows = []
	for columnIndex in indexes:
		
		#Returns length (min, max), mean, variance, skewness, kurtosis
		#col = [float(el) for el in m.column(columnIndex)]
		col = m.column(columnIndex)
		stats = scipy.stats.describe(col)
		row = list(stats[1])
		row.extend(stats[2:4])
		row[3] = math.sqrt(row[3])
		
		#Returns ShapiroValue, Shapiro Prob, Dangostino Value, Dangostino Prob
		#Kolomogrov Value, Kolomogrov Prob, histogram, qq
		stats = normalityTests(col)
		qq = stats[-1]
		histogram = stats[-2]
		row.extend(stats[:-2])

		name = m.headerForColumn(columnIndex)
		row.insert(0, name)
		rows.append(row)
		
		if options.outputQQ:
			f = open('QQ%s.csv' % name, 'w+')
			f.write(qq.csvRepresentation())
			f.close()

		if options.outputHist:
			f = open('ErrorHist%s.csv' % name, 'w+')
			f.write(histogram.csvRepresentation())
			f.close()

	headers = ['Name', 'Min', 'Max', 'Mean', 'Stdev', 'Shapiro', 'ShaprioProb', 
			'DAngostino', 'DAngostino', 'KS', 'KSProb']
	matrix = Core.Matrix.Matrix(rows=rows, headers=headers)
	print matrix.csvRepresentation(),


