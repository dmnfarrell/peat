#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
import math, numpy, sys, os, copy
import csv
import matplotlib  
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from PEATDB.Ekin.Fitting import Fitting
from PEATSA import Core


class PCAPlugin(Plugin):
    """PCA plugin"""
    #capabilities = ['gui']
    menuentry = 'PCA Plugin'
    gui_methods = {} 
    about = 'This plugin allows you to do PCA'    
    
    def main(self, parent=None, DB=None):       
        if parent==None:
            self.DB = DB       
        else:
            return

    def writeToFile(self, matrix, filename):                    
        stream = open(filename, 'w')
        stream.write(matrix.csvRepresentation())
        stream.close()
        return
	
    def doPCA(self,m,standardize=True):
        '''Performs pca on the Core.Matrix.Matrix instance m.

        Returns: 
                eigenvalues - A numpy 1-D array.
                eigenvectors - A numpy 2-D array
                transformedData - A numpy 2-D array'''

        print >>sys.stderr, 'Calculating mean vector'
        data = numpy.array(m.matrix)
        if standardize==True:
            data = self.standardize(data)   
        average = numpy.zeros(m.numberOfColumns())
        for row in data:                
            row = numpy.array(row, dtype=numpy.float)               
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

    def standardize(self, data):
        """standardize data"""
        import scipy.stats as st
        newdata = copy.deepcopy(data)
        i=0
        for col in zip(*data):
            newdata[:,i] = st.zs(col)            
            i+=1
        print newdata    
        return newdata
    
    def plotResults(self,evals,evecs,b,m):
        """Plot results to help visualize components"""
        data = numpy.array(m.matrix)
        labels = m.columnHeaders()
        plt.rc('font', family='monospace')

        #plot first 2 PCs in 3d score plot        
        '''x,y,z = b[:,0], b[:,1], b[:,2]
        f=plt.figure()
        ax = Axes3D(f)
        ax.scatter(x,y,zs=z,marker='o',lw=2,alpha=0.5,c='b',s=30)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')'''
        
        #f.subplots_adjust(hspace=1,wspace=1)
        
        f=plt.figure()
        i=1
        length = len(data[0])
        if length>6: length=6
        for i in range(0,length):
            ax=f.add_subplot(3,3,i+1)
            c=0
            lines=[]
            for v in zip(*data):     
                c+=1
                if c>10: break
                #v = [float(j)/(max(v)-min(v)) for j in v]               
                l=ax.plot(b[:,i],v,'x',mew=1,alpha=0.2)
                lines.append(l)
                ax.set_title('Ev%s' %str(i+1))
                
            i+=1    
        f.legend(lines,labels,loc='lower right')        
        ax=f.add_subplot(337)
        ind=numpy.array(range(len(evals)))
        ax.plot(ind,evals,'-o',lw=2)
        ax.set_xlabel('Eigenvalues')
        f.savefig('PCAresults.png')
        f.subplots_adjust(hspace=0.4,wspace=0.4)
        print 'Eigenvalues: ', evals
        print 'Eigenvectors: ', evecs
        plt.show()
        return
    
    def test(self):
        
        features=['stab','act','solv','res','asa']        
        x = numpy.random.normal(2, 6, 500)
        y = numpy.random.normal(4, 1, 500)
        #y = [i+numpy.random.normal(2,0.3) for i in x]        
        z = [i+numpy.random.normal(2,0.24) for i in y]
        #z = numpy.random.normal(4, 1, 500)
        s = numpy.random.gamma(4, 1, 500)
        t = numpy.random.gamma(4, 1, 500)
        
        filename = 'testdata.csv'
        f=open(filename,'w')
        cw=csv.writer(f)
        cw.writerow(features)
        for i in zip(x,y,z,s,t):
            cw.writerow(i)
        f.close()
        
        '''A,X = Fitting.doFit(expdata=zip(x,y), model='Linear',silent=True)
        fitx = numpy.arange(min(x),max(x),1)
        fity = X.getFitLine(fitx)
        A,X1 = Fitting.doFit(expdata=zip(x,z), model='Linear',silent=True)
        fitz = X.getFitLine(fitx)'''
        
        f=plt.figure()
        ax = Axes3D(f)
        ax.scatter(x,y,zs=z,marker='o',lw=2,alpha=0.5,c='b',s=30)
        #ax.plot(fitx,fity,zs=fitz,alpha=0.6,lw=2,c='r',label='fit xyz')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.legend()
        f.subplots_adjust(hspace=1,wspace=1) 
        
        m = Core.Matrix.matrixFromCSVFile(filename)
        evals, evecs, z = self.doPCA(m)
        self.plotResults(evals, evecs, z, m)
        return
    
def main():
    import os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    parser.add_option("-t", "--test", dest="test", action='store_true',
                       help="test func", default=False) 
    parser.add_option("-s", "--start", dest="start", default=0, type="int",
                       help="start")
    parser.add_option("-e", "--end", dest="end", default=0, type="int",
                       help="end")
    parser.add_option("-z", "--standardize", dest="standardize", action='store_true',
                       help="end", default=False)    
    opts, remainder = parser.parse_args()
    
    P = PCAPlugin()
    if opts.file != None and os.path.exists(opts.file):
        r = Core.Matrix.matrixFromCSVFile(opts.file)
        if opts.start != None:
            m = r[:, opts.start:]
        print 'There are %d samples and %d variables (dof)' % (m.numberOfRows(), m.numberOfColumns())       
        evals, eigenvectors, z = P.doPCA(m, opts.standardize)
        P.plotResults(evals, eigenvectors, z, m)
        
	#Write out vectors
	ev = Core.Data.Matrix.Matrix(rows=list(eigenvectors))
	ev.addColumn(m.columnHeaders(), 0)
	headers = ['Variables']
	for i in range(m.numberOfColumns()):
            headers.append('Ev%d' % i)
	ev.setColumnHeaders(headers)
	P.writeToFile(ev, 'Eigenvectors.csv')

	#Write out new basis
	basis = Core.Matrix.Matrix(rows=list(z))
	basis.addColumn(r.column(0), 0)
	#basis.addColumn(r.column(1), 1)
	headers.pop(0)
	headers = r.columnHeaders()[:1] + tuple(headers)
	basis.setColumnHeaders(headers)
	P.writeToFile(basis, 'NewBasis.csv')
	
    if opts.test == True:
        P.test()   
         
if __name__ == '__main__':
    main()
