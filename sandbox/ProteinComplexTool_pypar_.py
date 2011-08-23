#! /usr/bin/env python

import os, sys, optparse
import PEATSA.Core as Core
import PEATSA.Core.Matrix
import MySQLdb
import matplotlib.pyplot as plt
import numpy as np


"""

Requires ProteinComplexTool_execute.py to do calculations!


"""


class ProteinComplexTool:
    def __init__(self):
        return
    
    def DeltaStability(self,inputFile, mutationList, configurationFile, workingDirectory, outputDirectory):

            '''Calculates the stability difference between a protein and set of mutants

            Parameters:
                    inputFile: A PDB file of the protein
                    mutationList: A list of Data.MutationSet instances. Each represents a mutant of the protein.
                    configurationFile: The location of a proteinDesignTool.conf file - defaults to home directory.
                    workingDirectory: Where the calculation will be run. 
                    outputDirectory: Where the results will be written.

            Returns
                    A Data.DataSet instance containing one matrix, stabilityResults.
                    Each row of this matrix corresponds to a mutant defined in the mutationList argument.'''

            #Create the ProteinDesignTool instance
            tool = Core.ProteinDesignTool.ProteinDesignTool(configurationFile, 
                            workingDirectory=workingDirectory,
                            pdbFile=inputFile, 
                            outputDirectory=outputDirectory,
                            removeHeterogens=True)

            #The above cleans the pdb file and copies it to the working directory.
            #Use this pdb from now on.
            inputFile = tool.pdbFile

            #Create the mutants
            mutantCollection = Core.Data.MutantCollection(pdbFile=inputFile,mutationList=mutationList,location=outputDirectory,temporary=True)

            #Run stability calculation
            #The results are added to the ProteinDesignTool instance's dataDirectory attribute
            #This is an instance of Data.DataSet class
            tool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles())

            #Clean up - Deletes files copied to the working directory for Uffbaps
            tool.cleanUp()

            return tool.dataDirectory

    def remALT(self,pdb): 
        
        '''Removes alternative residues from the working pdb. Replaces the working pdb.'''
        
        import Protool

        x = Protool.structureIO()
        x.readpdb('%s.pdb' % (pdb))
        x.RemoveALT()
        x.writepdb('%s.pdb' % (pdb), dont_write_HETATMS=1)
        print "[ProteinComplexTool] Alternative Residues removed."

    def splitter(self,pdbDir,pdb,reactions_list,cur,db):

        ''' Takes the reaction list provided by the user, ensures there is no duplication
            in the list and converts the list into the command required by MDAnalysis to split
            the pdb. If the user does not provide a list, it was automatically split the pdb
            into constituent chains.

            Returns a list of the chains split from the pdb.'''
        
        import string

        if reactions_list == ['']:
            # query the database
            environment.output(cur.execute("SELECT DISTINCT Chain_ID from pdb where PDB_ID = '%s';" % (pdb)))
            a = cur.fetchall() # fetch results
            print 'a', a
            expr=[]
            chains = [i[0] for i in a]
            for i in chains:
                s=["segid "+i]
                expr.append(s) 
            b = str(a) # convert from tuple to string
            exclude = set(string.punctuation) # set of punctutation characters
            b = ''.join(ch for ch in b if ch not in exclude) # remove punctuation from b
            e = ''.join(b.split(' '))
            self._do_split(pdbDir, pdb, expr, e)
            return e
        else:
            expr1=[]
            for c in reactions_list:
                if len(c)>1:
                    s = ["segid "+i for i in c]
                    expr1.append(s)
                else:
                    expr1.append(["segid "+c])
            self._do_split(pdbDir,pdb, expr1, reactions_list)
            return reactions_list

    def _do_split(self,pdbDir,pdb, expr, e):
        import MDAnalysis

        u = MDAnalysis.Universe(pdbDir, permissive=False)
        for i in range(len(expr)):
            Z = u.selectAtoms(*expr[i])
            Z.write('%s_%s.pdb' % (pdb,e[i]))
            print '[ProteinComplexTool] Extracted chain(s) %s from %s' % (e[i], pdb)

    
    def displayResults(self,pdb,split_list,comp_list,mutList,reaction,cur,db):
        
        ''' Reads the results from the MySQL db based on the split_list. The ddG binding is calculated here
            by subtracting the dGs of all the extracted chains (reactants) from the complex (products). The ddGs
            are then stored in a new table named based on the reaction list and scan type. 
            
            Automatic output is generated by ddG(multichain complex, line graph) or dG (single chain complex, 
            bar chart)'''

        width=0.5
        cur.execute("SELECT * FROM results_%s_%s;" % (split_list[0], mutList))
        complexResults = cur.fetchall()
        mutations = [i[0] for i in complexResults] # Mutation list
        complexScores = [i[1] for i in complexResults] # dG scores of pdb complex
        count = len(mutations) # Number of calcs
        ind = np.arange(count)
        """
        nums = []
        for y in muts:
            nums.append(''.join(i for i in y if i.isdigit()))
            """
        if len(split_list)>1:    # For binding calcs, no matter in what order chains were split
            chainResults = []
            for i in split_list[1:]:
                cur.execute("select * from results_%s_%s;" % (i,mutList))
                chainResults.append(cur.fetchall())
            
            chainScores = [i[1] for y in chainResults for i in y] # dG scores of chains split from pdb 
            ddG = []
            cur.execute("create table if not exists ddG_%s_%s(mutation VARCHAR(10), ddG FLOAT);" % (pdb, comp_list))
            for i in range(len(complexScores)):
                ddG.append(complexScores[i] - chainScores[i]) # ddG = dG complex - sum of dG of split chains
            for i in range(len(mutations)):
                print "ddG %s %s" % (mutations[i], ddG[i])
                cur.execute("insert into ddG_%s_%s (mutation, ddG) VALUES (%s%s%s, %s%s%s);" % (pdb, comp_list, '"', mutations[i], '"', '"',ddG[i],'"'))
            #plt.bar(ind+(width/2), ddG, width,color='b')
            #plt.axhline(linewidth=2, color='r')
            plt.figure(num=None, figsize=(16, 14), dpi=80)
            for i in range(len(mutations)):
                if mutations[i].startswith('A'):
                    plt.scatter(chainScores[i], complexScores[i], s=25,c='b',marker='o', linewidths=None)
                    #plt.text(chainScores[i], complexScores[i], mutations[i])
                if mutations[i].startswith('B'):
                    plt.scatter(chainScores[i], complexScores[i], s=25,c='r',marker='o', linewidths=None)
                    #plt.text(chainScores[i], complexScores[i], mutations[i])
                if mutations[i].startswith('C'):
                    plt.scatter(chainScores[i], complexScores[i], s=25,c='g',marker='o', linewidths=None)
                    #plt.text(chainScores[i], complexScores[i], mutations[i])
                if mutations[i].startswith('D'):
                    plt.scatter(chainScores[i], complexScores[i], s=25,c='y',marker='o', linewidths=None)
                    #plt.text(chainScores[i], complexScores[i], mutations[i])
            plt.plot([min(complexScores)*3,max(complexScores)*1.1],[min(complexScores)*3,max(complexScores)*1.1], 'r-')
            plt.xlabel("dG kJ/mol of the chains")
            plt.ylabel("dG kJ/mol of the complex")
            #plt.xlabel("Position of Mutation")
            #plt.ylabel("ddG kJ/mol")
            #plt.title("ddG Binding calculations for ALA scan of %s %s" % (split_list[0], reaction))
            plt.title("dG stability for ALA scan of %s %s complex vs chains" % (split_list[0], reaction))
        else:
            for i in range(len(mutations)):
                print " dG %s, %s" % (mutations[i], complexScores[i])
            plt.bar(ind,complexScores,width,color='r')
            plt.title("dG Stability calculations for ALA scan of %s %s" % (split_list[0], reaction))

        #plt.xticks(ind+(width/2), mutations, rotation=90, fontsize=8)
        #plt.figure(num=None, figsize=(8, 6), dpi=80)
        plt.show()
        sys.exit()

    def run_calcs(self, split_list, mutList, numProc):
        for i in split_list:
            w_pdb = os.path.join(os.getcwd(),'%s.pdb' % (i))
            os.system("mpirun -np %d python ProteinComplexTool_execute.py -p %s -d %s -m %s" %(opts.numProc, i, w_pdb, opts.mutList))
            print "[ProteinComplexTool] Calculations completed."
        
    
def main():
    import pypar

    run = ProteinComplexTool()
    
    # Option to select pdb, config file, working dir etc.. 
    parser = optparse.OptionParser()
    # PDB option
    parser.add_option("-p", "--pdb", 
                      help="Choose all or a pdb id", 
                      dest="pdb", default ="all")
    # Configuration File 
    parser.add_option("-c", "--configurationFile", 
                      help="Location of configuration file", 
                      dest="configFile", default="/home/satnam/proteinDesignTool.conf")
    # Output Directory
    parser.add_option("-o", "--outputDirectory", 
                      help="Location of output directory", 
                      dest="outputDir", default=os.getcwd())
    # Working Directory 
    parser.add_option("-w", "--workingDirectory", 
                      help="Location of working directory", 
                      dest="workingDir", default=os.getcwd())
    # Mutation List or ALA scan option
    parser.add_option("-m", "--mutationList", 
                      help="Location of mutation list file", 
                      dest="mutList", default="ALA")
    # Choose option for user-defined calculations
    parser.add_option("-u", "--userCalcs", 
                      help="Choose True or False if you would like to specifiy the calculations, otherwise each chain will be split", 
                      dest="userCalcOpt", default='True')
    # Show Results
    parser.add_option("-s", "--showResults", 
                      help="Shows previous results? True or False. If they don't exist, they will be calculated.", 
                      dest="showResults", default=True)
    # Delete results from database
    parser.add_option("-d", "--deleteResults", 
                      help="Deletes all results for the specified pdb from the database. Default False.", 
                      dest="deleteResults", default='False')
    # Number of CPUs for execute script
    parser.add_option("-n", "--numProc",
                      help="Number of Processors allocated for execute script.",
                      dest="numProc", default=8)
    

    (opts, args) = parser.parse_args()
    
    db = MySQLdb.connect(host="localhost", 
                             user = "root", 
                             passwd = "samsung", 
                             db = "sat")
    cur = db.cursor()
    
    # Connect ot the database and show PDBs
    cur.execute("SELECT VERSION()")
    ver = cur.fetchone()
    print "[ProteinComplexTool] MySQLdb connection successful!"
    print "[ProteinComplexTool] MySQL server version: %s" % ver[0]
    cur.execute("SELECT distinct PDB_ID from pdb;")
    a = cur.fetchall()
    b = ','.join([i[0] for i in a])
    print "[ProteinComplexTool] PDBs in local db"
    print b

    # PDB filename and directory handling
    pdb = opts.pdb
    pdbfile = ''.join((pdb,'.pdb'))
    pdbDir = os.path.join(opts.outputDir,pdbfile)

    # Delete tables already in database
    if opts.deleteResults != 'False':
        cur.execute("SHOW tables like 'results_%s%s';" % (pdb, '%'))
        drop_tables=cur.fetchall()
        for i in drop_tables:
            for y in i:
                cur.execute("DROP TABLE %s;" % (y))
                print "[ProteinComplexTool] Results for %s deleted" % (pdb)
    else:
        pass
    
    # Query database and displays entity and chain information of pdb
    cur.execute("SELECT distinct Entity_ID, Chain_ID, Chain_name, type from pdb where PDB_ID = %s;", (pdb))
    entity = [] # entities in the pdbfile
    chains = []

    for i in cur.fetchall():
        print "Entity: %s, Chain Name: %s, Type: %s, Chain ID: %s" % (i[0], i[2], i[3], i[1])  
        entity.append(i[0])
        chains.append(i[1])               
    entity.sort()

    # Remove alternative residues
    run.remALT(pdb)

    # User defined interactions. If left blank will extract each chain from pdb.
    if opts.userCalcOpt == 'True':
        print "[ProteinComplexTool] What components are consumed (enter chain IDs in the form AB+C+D):"
        reactants = sys.stdin.readline()
        reactants = reactants.rstrip("\n")
        print "[ProteinComplexTool] What products are produced (enter chain IDs in the form ABC+D):"
        products = sys.stdin.readline()
        products = products.rstrip("\n")
        if reactants == '':
            reactants = '+'.join(chains)
        else:
            pass
    
    reaction = reactants+'-->'+products
    print reaction

    if products == '':
        products = ''.join(chains)
    else:
        pass

    reactants_list = reactants.split('+')
    products_list = products.split('+')

    # Split the pdb into chains, returns chains that have been split (A,B etc)
    split_reactants = run.splitter(pdbDir,pdb,reactants_list,cur,db)        
    split_products = run.splitter(pdbDir,pdb,products_list,cur,db)
    comp_list = split_products + split_reactants
    comp_list = '_'.join(comp_list)

    split_list = []
    split_list_products = []
    split_list_reactants = []

    for i in split_reactants:
        s = pdb+'_'+i
        split_list_reactants.append(s)

    for i in split_products:
        s = pdb+'_'+i
        split_list_products.append(s)

    # Comp_list, used for naming of ddG table - pdb_prods_reacts_mut-type
    comp_list = comp_list +'_'+os.path.split(opts.mutList)[1]
    
    # Ensuring no duplicate reactions are run, espicially when user leaves input blank.
    splitlist = split_list_products + split_list_reactants
    for i in splitlist:
        if i not in split_list:
            split_list.append(i)  

    
    # Displays results results
    if opts.showResults == 'True':
        run.run_calcs(split_list, opts.mutList, opts.numProc)
        run.displayResults(pdb, split_list, comp_list, opts.mutList, reaction, cur, db)
    else:
        run.run_calcs(split_list, opts.mutList, opts.numProc)
        """
        # Runs calculations and send the names of the split fes to the execute script
        for i in split_list:
            w_pdb = os.path.join(os.getcwd(),'%s.pdb' % (i))
            os.system("mpirun -np %d python ProteinComplexTool_execute.py -p %s -d %s -m %s" %(opts.numProc, i, w_pdb, opts.mutList))
            print "[ProteinComplexTool] Calculations completed."
        #run.displayResults(pdb, split_list, comp_list, reaction, cur, db)
        """
    

if __name__ == '__main__':
    main()
