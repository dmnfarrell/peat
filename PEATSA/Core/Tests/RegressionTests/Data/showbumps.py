#!/bin/env python

import Protool
X=Protool.structureIO()
X.readpdb('1crn.pdb')
import Protool.mutate
M=Protool.mutate.Mutate(X)
M.show_all_bumps()
