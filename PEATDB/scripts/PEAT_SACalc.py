#!/usr/bin/env python
#
# This file is part of Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#

import sys,os
import string
import PEAT_DB.Database as DB
import PEAT_SA.Core as Core
from PEATDB.Base import PDatabase

class PEATrunner(object):
    """Gets the mutants in PEAT_DB and performs PEAT_SA calculations on them.
       This data can be saved back into the DB"""

    def __init__(self,datadir=None):
        """Initialise PT as a DB object if a datadir is given"""

        self.datadir=datadir
        if datadir == None:
            self.PT = None
            return
        else:
            self.PT = DB.Database(self.datadir, Tk=False)
            self.listRecs()
            return

    def setRefProtein(self):
        """Get the reference protein"""

        self.refprot = self.PT.get_refprotein()
        if refprot is None:
            print 'There is no ref protein'
            return
        else:
            print 'ref protein:', self.refprot

        '''self.DBI is an instance of Database
        You can set the refprotein with this command:
        self.DBi.set_refprotein('wt')'''


    def getMutants(self):
        """Get mutant structures and pass them to PEAT_SA"""
        PT=self.PT
        bold = "\033[1m"; reset = "\033[0;0m"
        #logfile for failed proteins
        failed=open("failed.txt", 'w')
        #pdb directory
        pdbdir = os.path.join(os.path.split(self.datadir)[0])
        name = os.path.split(self.datadir)[1]

        #get wt structure and put in protool instance ?
        wtpdblines = PT.DB['wt']['Structure']

        #save it to a pdbfile
        fd=open(os.path.join(pdbdir,'wt.pdb'),'wb')
        for line in wtpdblines:
            fd.write(line)
        fd.close()

        #set the ref protein?
        PT.set_refprotein('wt')
        self.refprot = PT.get_refprotein()
        print 'Got ref protein', self.refprot

        mutantCollection = Core.Data.MutantCollection(pdbFile='wt.pdb', name=name,
                                                    location=pdbdir, clean='no')
        print 'created mutant collection instance'

        i=0
        for protein in PT.proteins:
            #get operations using ref protein
            print bold + protein + reset
            is_parent, operations = PT.isparent(protein,self.refprot)
            if not operations:
                print 'no alignment for', protein
                failed.write(protein+'\n')
                continue
            print 'got operations', operations

            sys.stdout = open("log.txt", 'a')
            i+=1
            try:
                pdblines, X = PT.model_on_the_fly(PT.DB[protein]['Structure'])
            except:
               print 'failed to get model'
            sys.stdout = sys.__stdout__

            #get mutant info from operations and pass to mutantcollection
            mutationSet = Core.Data.MutationSet(name=protein)
            for op in operations:

                ch, ind, rn, mt = self.getmutantInfo(op)
                if len(ch) != 1:
                    continue
                #does one for each point mutation in operations
                mutationSet.addMutation(chain=ch, residueIndex=ind, mutation=mt)

                #uncomment this to only test with ALA mutants
                #if mt != 'ALA':
                #    print 'non ALA mutation, skipping'
                #    mutationSet=None
                #    break
            if mutationSet!=None:
                mutantCollection.addMutant(mutant=X, mutationSet=mutationSet)
                print 'adding mutation'

        print 'Got mutant collection:'
        print mutantCollection.mutantFiles()
        failed.close()
        return mutantCollection

    def getmutantInfo(self, ops):
        """Get mutant info from protool operations eg. ['A:0166:THR:PHE']"""
        print ops.split(':')
        chain, resindex, resname, mutation = ops.split(':')
        return chain, resindex, resname, mutation

    def performCalculations(self, mCollection, stab=True, pka=True):
        """Perform the PEAT_SA calculations"""

        pdTool = Core.ProteinDesignTool.ProteinDesignTool(
                                configurationFile=None,
                                workingDirectory='/local/farrell/temp',
                                pdbFile='/local/farrell/wt.pdb',
                                outputDirectory='/local/farrell',
                                dataName='10R')

        pdTool.configuration.set('PKA SCAN PARAMETERS',  'recalc_intpka_dist', '10')
        #pdTool.configuration.set('PKA SCAN PARAMETERS',  'pKMCsteps', '20000')
        #pdTool.configuration.writeToFile('pdtool_config.txt')
        names = [mutation.name for mutation in mCollection.mutations()]
        if stab == True:
            #dstability calculations
            pdTool.runStabilityCalculation(mCollection.mutantFiles(), verbose=True)
            pdTool.dataDirectory.stabilityResults.addColumn(names, index=0)
            pdTool.dataDirectory.stabilityResults.setColumnHeader(0, 'Name')
        if pka == True:
            #dpKa calculations
            pdTool.runScan(mCollection, verbose=True)
            pdTool.dataDirectory.scanResults.addColumn(names, index=0)
            pdTool.dataDirectory.scanResults.setColumnHeader(0, 'Name')
        pdTool.dataDirectory.synchronize()
        return


def main():
    "Run the application"
    import sys
    from optparse import OptionParser
    parser = OptionParser()

    parser.add_option("-d", "--dir", dest="dir",
                        help="open a PEAT project dir", metavar="FILE")
    parser.add_option("-m", "--mut", dest="mut",
                        help="use mutant collection", metavar="FILE")
    parser.add_option("-s", "--stab", action='store_true', dest="stab",
                        help="do stability calculations", default=False)
    parser.add_option("-p", "--pka", action='store_true', dest="pka",
                        help="do dpka calculations", default=False)

    opts, remainder = parser.parse_args()
    print opts

    if opts.dir != None:
        P = PEATrunner(datadir=opts.dir)
        #get the mutants from PEATDB and do mutantcollection
        m = P.getMutants()
    elif opts.mut != None:
        #no dir given, we don't need the peat project if mutantcollection available
        P = PEATrunner()
        m = Core.Data.MutantCollection(name=opts.mut, location='/local/farrell')
        print m
    if opts.stab == False and opts.pka == False:
        print 'Not doing any calculations, use -s or -p'
    else:
        P.performCalculations(m, stab=opts.stab, pka=opts.pka)

    return

if __name__ == '__main__':
    main()
