#! /bin/env python
# Written by Michael Johnston, 2010
'''Performs PCA on an input matrix. Rows of the matrix correspond to samples, while
columns correspond to the variables. The range of columns containing the variales to be used can be specified.
Outputs eigenvectors to Eigenvectors.csv and the original data transformed to a new basis is output to NewBasis.csv'''
import numpy
import sys, optparse
from PEAT_SA import Core

def writeToFile(matrix, filename):
		
	stream = open(filename, 'w')
	stream.write(matrix.csvRepresentation())
	stream.close()

def pca(m):

	'''Performs pca on the Core.Matrix.Matrix instance m.

	Returns: 
		eigenvalues - A numpy 1-D array.
		eigenvectors - A numpy 2-D array
		transformedData - A numpy 2-D array'''

	print >>sys.stderr, 'Calculating mean vector'
	data = numpy.array(m.matrix)
	average = numpy.zeros(m.numberOfColumns())
	for row in data:
		average = average + row

	average /= m.numberOfRows()
	temp = zip(m.columnHeaders(), average)
	print >>sys.stderr, 'Result: '
	for el in temp:
		print >>sys.stderr, '\t%s: %lf' % tuple(el)

	print >>sys.stderr, '\nMean-Centering'
	data = data - numpy.tile(average, [m.numberOfRows(),1])

	print >>sys.stderr, 'Calculating covariance matrix'
	cov = numpy.cov(data, rowvar=0)

	print >>sys.stderr, 'Performing eigenvalue decomposition'
	eigenvalues, eigenvectors = numpy.linalg.linalg.eig(cov)
	eigenvectors = eigenvectors.astype(numpy.float32)

	print >>sys.stderr, 'Sorting'
	x = range(len(eigenvalues))
	x.sort(lambda x,y: cmp(eigenvalues[x], eigenvalues[y]), reverse=True)
	eigenvalues = eigenvalues[x]
	eigenvectors = eigenvectors[:,x]

	print >>sys.stderr, 'Complete'

	z = numpy.dot(data, eigenvectors)
	
	return eigenvalues, eigenvectors, z

if __name__ == "__main__":

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-f", "--file", dest="file",
			  help="A csv files", metavar="FILE")
	parser.add_option("-s", "--startColumn", dest="start", default=0, type="int",
			  help="The index of the column contains the first variable. Defaults to %default", metavar="START")
	parser.add_option("-l", "--length", dest="length", type="int",
			  help="The number of columns containing variables. Defaults to the number of columns - START", metavar="LENGTH")
	parser.add_option("", "--useComponents", dest="components",
			  help="Options: A comma separated list of component indexes (0 for 0th component)\n" 
			"The input matrix will be reduced to only include contributions along specified components", 
			metavar="COMPONENTS")

	(options, args) = parser.parse_args()

	if options.file is None:
		print 'CSV file must be provided'
		sys.exit(1)
	else:
		r = Core.Matrix.matrixFromCSVFile(options.file)

	if options.length is None:
		options.length = r.numberOfColumns() - options.start

	project = False
	if options.components is not None:
		project = True
		components = [int(el.strip()) for el in options.components.split(',')]

	#Get sub matrix containing just samples
	end = options.length + options.start
	m = r[:, options.start:end]

	print >>sys.stderr, 'There are %d samples and %d variables (dof)' % (m.numberOfRows(), m.numberOfColumns())

	eigenvalues, eigenvectors, z = pca(m)

	print 'Eigenvalues: ', eigenvalues

	#Write out vectors
	ev = Core.Data.Matrix.Matrix(rows=list(eigenvectors))
	ev.addColumn(m.columnHeaders(), 0)
	headers = ['Variables']
	for i in range(m.numberOfColumns()):
		headers.append('Ev%d' % i)

	ev.setColumnHeaders(headers)
	writeToFile(ev, 'Eigenvectors.csv')

	#Write out new basis
	basis = Core.Matrix.Matrix(rows=list(z))
	basis.addColumn(r.column(0), 0)
	basis.addColumn(r.column(1), 1)
	headers.pop(0)
	headers = r.columnHeaders()[:2] + tuple(headers)
	basis.setColumnHeaders(headers)
	writeToFile(basis, 'NewBasis.csv')

	if project:
		data = numpy.array(m.matrix)
		rows = []
		for row in data:
			result = numpy.zeros(m.numberOfColumns())
			#For each chosen component find the projection
			#Get the corresponding vector in the original space
			for component in components:
				projection = numpy.dot(row, eigenvectors[:, component])
				result += projection*eigenvectors[:, component]
			rows.append(list(result))

		projections = Core.Matrix.Matrix(rows=rows)
		projections.addColumn(r.column(0), 0)
		projections.addColumn(r.column(1), 1)
		projections.setColumnHeaders(r.columnHeaders())
		writeToFile(projections, 'Projections.csv')

