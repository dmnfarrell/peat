import PEATSA.Core as Core
import PEATSA.Core.Matrix
import optparse
import pypar, os
import MySQLdb


"""

Requires ProteinComplexTool_pypar.py to do calculations!


"""

          
proc = pypar.size()
myid = pypar.rank()
node = pypar.get_processor_name()

def isRoot(myid):
    if myid == 0:
        return True
    else:
        return False

def DeltaStability(inputFile, mutationList, configurationFile, workingDirectory, outputDirectory):

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


def do_run(pdb, i, cur, db, mutationList):
    
    if mutationList != "ALA":
        mfile = Core.Data.MutationListFile(filename=mutationList,create=True)
        mfile.removeDuplicates(autoUpdate=False)
        mutList = mfile.mutantList()
        if isRoot(myid):
            print mfile.numberOfMutants()
    else:
        mutList = Core.Data.CreateScanList(pdbFile=i, mutation='ALA', skipResidueTypes=['ALA', 'GLY'])
    
    results = DeltaStability(inputFile = i, mutationList = mutList, configurationFile='/home/satnam/proteinDesignTool.conf', workingDirectory = os.getcwd(), outputDirectory = os.getcwd())
    
    # Results are submitted to results_pdb+chain and only by one processor
    if isRoot(myid):
        cur.execute("create table if not exists results_%s_%s(mutation VARCHAR(20), score FLOAT);" % (pdb,os.path.split(mutationList)[1]))
        for mutant in range(results.stabilityResults.numberOfRows()):
            cur.execute("insert into results_%s_%s (mutation, score) VALUES (%s%s%s, %s%s%s);" % (pdb,os.path.split(mutationList)[1], '"', results.stabilityResults[mutant][0], '"', '"', results.stabilityResults[mutant][-1],'"'))
        print "Calculated %s stability and results added to database" % (pdb)
            
    pypar.finalize()


def main():
    
    # Ensure all Processors are ready
    pypar.barrier()
    print "Processor %d is ready" % (myid)
    
    # Connect to MySQL db
    db = MySQLdb.connect(host="localhost", 
                         user = "root", 
                         passwd = "samsung", 
                         db = "sat")
    cur = db.cursor()


    # Option parser from wrapper script
    parser = optparse.OptionParser()
    # PDB
    parser.add_option("-p", "--pdb", 
                      help="Choose all or a pdb id", 
                      dest="pdb", default ="all")
    # PDB directory
    parser.add_option("-d", "--dir", 
                      help="i", 
                      dest="i", default ="all")

    parser.add_option("-m", "--mutationList", 
                      help="Location of mutation list file", 
                      dest="m", default="ALA")
    
    (opts, args) = parser.parse_args()
    
    # Run calculations
    do_run(opts.pdb, opts.i, cur, db, opts.m)

    # Finalize and exit
    pypar.finalize()

if __name__ == '__main__':
    main()
