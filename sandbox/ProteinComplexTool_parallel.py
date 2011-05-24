#! usr/bin/env python

import optparse, os, csv, glob, sys
import MySQLdb
import PEATSA.Core as Core
import PEATSA.Core.Matrix
import matplotlib.pyplot as plt
import numpy as np
import Environment


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

    def remALT(self,pdbfile,environment): 
        import Protool

        x = Protool.structureIO()
        x.readpdb('%s.pdb' % (pdbfile))
        x.RemoveALT()
        x.writepdb('%s.pdb' % (pdbfile), dont_write_HETATMS=1)
        environment.output('Removed alternate residues')

    def splitter(self,pdbDir,pdb,reactions_list,cur,db,environment):
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
            self.do_split(pdbDir, pdb, expr, e)
            return e
        else:
            expr1=[]
            for c in reactions_list:
                if len(c)>1:
                    s = ["segid "+i for i in c]
                    expr1.append(s)
                else:
                    expr1.append(["segid "+c])
            environment.wait()
            self.do_split(pdbDir,pdb, expr1, reactions_list, environment)
            return reactions_list

    def do_split(self,pdbDir,pdb, expr, e,environment):
        import MDAnalysis

        u = MDAnalysis.Universe(pdbDir, permissive=False)
        for i in range(len(expr)):
            environment.output(expr[i])
            Z = u.selectAtoms(*expr[i])
            Z.write('%s_%s.pdb' % (pdb,e[i]))
            environment.output('Extracted chain(s) %s from %s' % (e[i], pdb))

    def createMutlist(self,pdb):
        mutList = Core.Data.CreateScanList(pdb, mutation='ALA', skipResidueTypes=['ALA', 'GLY'])

        return mutList

    def displayResults(self,pdb,split_list,comp_list,cur,db, environment):
        
        width=0.5
        environment.output(cur.execute("SELECT * FROM results_%s;" % (split_list[0])))
        complexResults = cur.fetchall()
        mutations = [i[0] for i in complexResults] # Mutation list
        complexScores = [i[1] for i in complexResults] # dG scores of pdb complex
        count = len(mutations) # Number of calcs
        ind = np.arange(count)
        
        if len(split_list)>1:    # For binding calcs, no matter in what order chains were split
            chainResults = []
            for i in split_list[1:]:
                cur.execute("select * from results_%s;" % (i))
                chainResults.append(cur.fetchall())
            chainScores = [i[1] for y in chainResults for i in y] # dG scores of chains split from pdb 
            ddG = []
            cur.execute("create table if not exists ddG_%s_%s(mutation VARCHAR(10), ddG FLOAT);" % (pdb, comp_list))
            for i in range(len(complexScores)):
                ddG.append(complexScores[i] - chainScores[i])
            for i in range(len(mutations)):
                 environment.output("ddG %s %s" % (mutations[i], ddG[i]))
                 cur.execute("insert into ddG_%s_%s (mutation, ddG) VALUES (%s%s%s, %s%s%s);" % (pdb,comp_list, '"', mutations[i], '"', '"',ddG[i],'"'))
            plt.plot(ind+(width/2), ddG, 'o-')
            plt.axhline(linewidth=2, color='r')
            plt.title("ddG Binding calculations for ALA scan of %s" % (split_list[0]))
        else:
            for i in range(len(mutations)):
                environment.output("%s, %s" % (mutations[i], complexScores[i]))
            plt.bar(ind,complexScores,width,color='r')
            plt.title("dG Stability calculations for ALA scan of %s" % (split_list[0]))

        plt.xticks(ind+(width/2), mutations, rotation=90, fontsize=8)
        plt.show()
        sys.exit()

def main():
    
    # Run program
    environment = Environment.Environment()
    # Connect to local database containing info about BMP pdbs
    db = MySQLdb.connect(host="localhost", user = "root", passwd = "samsung", db = "sat")
    cur = db.cursor()
    cur.execute("SELECT VERSION()")
    ver = cur.fetchone()
    environment.output( "MySQLdb connection successful")
    environment.output("MySQL server version: %s" % ver[0])

    # Show pdbs in database
    cur.execute("SELECT distinct PDB_ID from pdb;")
    a = cur.fetchall()
    b = ','.join([i[0] for i in a])
    environment.output("PDBs in database:")
    environment.output(b) 


    # Option to select pdb, config file, working dir etc.. 
    parser = optparse.OptionParser()
    # PDB option
    parser.add_option("-p", "--pdb", help="Choose all or a pdb id", dest="pdb", default ="all")
    # Mutation List or ALA scan option
    parser.add_option("-m", "--mutationList", help="Location of mutation list file", dest="mutList", default="ALA")
    # Configuration File 
    parser.add_option("-c", "--configurationFile", help="Location of configuration file", dest="configFile", default="/home/satnam/proteinDesignTool.conf")
    # Output Directory
    parser.add_option("-o", "--outputDirectory", help="Location of output directory", dest="outputDir", default=os.getcwd())
    # Working Directory 
    parser.add_option("-w", "--workingDirectory", help="Location of working directory", dest="workingDir", default=os.getcwd())
    # Choose option for user-defined calculations
    parser.add_option("-u", "--userCalcs", help="Choose True or False if you would like to specifiy the calculations, otherwise each chain will be split", dest="userCalcOpt", default=False)
    # Show Results Option
    parser.add_option("-s", "--showResults", help="Shows previous results? True or False. If they don't exist, they will be calculated.", dest="showResults", default=True)
    # Delete results from database
    parser.add_option("-d", "--deleteResults", help="Deletes all results for the specified pdb from the database. Default False.", dest="deleteResults", default=False)

    (opts, args) = parser.parse_args()

    # Instantiate the class
    run = ProteinComplexTool()
    
    # pdb name/file handling
    pdb = opts.pdb
    pdbFile = ''.join((pdb,'.pdb'))
    pdbDir = os.path.join(opts.outputDir,pdbFile)
    environment.output(pdbDir) 

    # Checking if user selected PDB is in the database
    """
    if opts.pdb != None:
        if opts.pdb not in b:
            raise sys.exit('PDB not in Database, choose one from list')
            """
    if opts.pdb in b:
        environment.output('PDB in Database')
        environment.output('Checking what calculations can be performed')


    # Check what calcs can be done with user defined PDB
    cur.execute("SELECT distinct Entity_ID, Chain_ID, Chain_name, type from pdb where PDB_ID = %s%s%s;" % ('"',pdb,'"'))
    entity = [] # entities in the pdbfile
    chains = []


    for i in cur.fetchall():
        environment.output("Entity: %s, Chain Name: %s, Type: %s, Chain ID: %s" % (i[0], i[2], i[3], i[1]))    
        entity.append(i[0])
        chains.append(i[1])               

    entity.sort()
    
    # Delete results
    if opts.deleteResults == 'True':
        cur.execute("SHOW tables like 'results_%s%s';" % (pdb, '%'))
        drop_tables=cur.fetchall()
        environment.output(drop_tables)
        for i in drop_tables:
            cur.execute("DROP TABLE '%s';" % (i))
            environment.output("Results for %s deleted" % i)
    else:
        pass

    # Remove Alternate Residues from pdb, will overwrite the file
    run.remALT(pdb,environment)    

    # User defined splitting of chains from PDB, can be left 
    # blank and the PDB will be split to individual chains
    reactions_list = ['']
    
    if opts.userCalcOpt != 'False':
        environment.output("What components are consumed (enter chain IDs in the form AB+C+D):")
        reactants = sys.stdin.readline()
        reactants = reactants.rstrip("\n")
        environment.output("What products are produced (enter chain IDs in the form ABC+D):")
        products = sys.stdin.readline()
        products = products.rstrip("\n")
        #reactants = raw_input("What components are consumed (enter chain IDs in the form AB+C+D):")
        #products = raw_input("What products are produced (enter chain IDs in the form ABC+D):")
    else:
        pass
    
    # If user leaves input blank, then the default is to calculate 
    # every chain individually vs complex  
    if reactants == '':
        reactants = '+'.join(chains)
    else:
        pass
    
    if products == '':
        products = ''.join(chains)
    else:
        pass
    
    reactants_list = reactants.split('+')
    products_list = products.split('+')
                 
    # Split the pdb into chains, returns chains that have been split (A,B etc)
    split_reactants = run.splitter(pdbDir,pdb,reactants_list,cur,db,environment)
    split_products = run.splitter(pdbDir,pdb,products_list,cur,db,environment)
    comp_list =  split_products + split_reactants
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
    
    splitlist = split_list_products + split_list_reactants
    for i in splitlist:
        if i not in split_list:
            split_list.append(i)

    # Split_list is a list of the pdb and the individual pdbs 
    # that have been split
    environment.output(split_list)
        
    # Show results
    if opts.showResults == 'True':
        count = 0
        environment.output(cur.execute("show tables;"))
        tables = cur.fetchall()
        resTable = "".join(("results_",pdb))
        for i in tables:
            for y in i:
                if y.startswith(resTable):
                   count +=1
        if count != 0:
            run.displayResults(pdb,split_list,comp_list,cur,db,environment)
    else:
        pass

    comp_list = comp_list +'_'+os.path.split(opts.mutList)[1]
    # Run the calculations
    # Load and check mutant list given by user, else do ALA scan
    
    if opts.mutList != "ALA":
        mfile = Core.Data.MutationListFile(filename=opts.mutList,create=True)
        mfile.removeDuplicates(autoUpdate=False)
        mutList = mfile.mutantList()    
    else:
        for i in split_list:
            w_pdb = os.path.join(opts.outputDir,'%s.pdb' % (i))
            mutList = Core.Data.CreateScanList(pdbFile=w_pdb, mutation='ALA', skipResidueTypes=['ALA', 'GLY'])
        
    muts = []
    score = []
    for i in split_list:
        w_pdb = os.path.join(opts.outputDir,'%s.pdb' % (i))
        #mutList = Core.Data.CreateScanList(pdbFile=w_pdb, mutation='ALA', skipResidueTypes=['ALA', 'GLY'])
        results = run.DeltaStability(inputFile=w_pdb, 
                  mutationList=mutList, 
                  configurationFile=opts.configFile, 
                  workingDirectory=opts.workingDir, 
                  outputDirectory=opts.outputDir) 
        for mutant in range(results.stabilityResults.numberOfRows()):
            #j = results.stabilityResults[mutant][0]
            #k = results.stabilityResults[mutant][-1]
            #cur.execute("insert into results (PDB_ID, mutation, score) VALUES (%s,%s,%s);", (i,j,k))")
            environment.output(cur.execute("insert into results (PDB_ID, mutation, score) VALUES (%s,%s,%s);", (i,results.stabilityResults[mutant][0],results.stabilityResults[mutant][-1])))
            print "keep going"
        #environment.wait()
        environment.output("Calculated %s stability and results added to database" % (i))
    # Display results
    #run.displayResults(pdb,split_list,comp_list,cur,db,environment)
    
    

if __name__=='__main__':
    main()


