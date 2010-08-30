#! /bin/env python

import PEAT_SA.Core as Core
import PEAT_SA.Core.Classify as Classify
import sys

if __name__ == '__main__':   
        
        responseVariable = sys.argv[2]
        matrix = Core.Matrix.matrixFromCSVFile(sys.argv[1])
        C = Classify.StructureClassifier()
        C.doRun(matrix, responseVariable)
        
