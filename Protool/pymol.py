#!/usr/bin/env python

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] 
import sys, time, os
from pymol import cmd
#pymol.finish_launching()
pymol.cmd.load('2LZT_H.pdb')
