#! /bin/env python

import PEATSA.Core as Core
import PEATSA.Core.Classify as Classify
import sys

if __name__ == '__main__':   
        
        responseVariable = sys.argv[2]
        matrix = Core.Matrix.matrixFromCSVFile(sys.argv[1])
        C = Classify.StructureClassifier()
        C.doRun(matrix, responseVariable)
        
