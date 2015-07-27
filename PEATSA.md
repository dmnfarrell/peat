## Introduction ##

PEATSA is a part of the PEAT suite of programs for storing, analysing and predicting biophysical data on proteins It allows you to analyse the effect of each amino acid residue on the biophysical characteristics of a protein. Currently PEAT\_SA computes changes in,

  * pKa values
  * Protein stability
  * Protein-ligand binding affinity
for a given set of point mutations.

PEATSA unifies the UFFBAPS and pKD programs and performs most of the necessary preparation steps e.g. adding hydrogens, automatically. This makes it easier for users to run the calculations they want as little to no technical knowledge about these programs is required to use PEATSA.

## Using PEATSA ##

There are currently three ways of using PEATSA. PEATSA is accessible via a web-server on enzyme.ucd.ie
We have also created a module that allows clients to interact with remote installations of PEATSA

## Installing PEATSA on a local machine ##
The rest of this document describes how to install PEATSA locally. In the following the three CalculationTypes are referred to as - Scan, Stability and Binding.

## Obtaining PEATSA ##

PEATSA is licensed under the GPL and the source code is available in our GoogleCode repository.
Currently we do not have a convenient installation procedure for PEATSA and installation will therefore involve a great deal of manual tweaking. Potential users should contact Jens.Nielsen@ucd.ie to arrange installation. However we have also made it possible for client to interact with remote a PEATSA Server which means a single PEATSA installation can be created that can be accessed by multiple users in different locations.
Once you have obtained PEATSA, look at the InstallationNotes to familiarize yourself with the package dir layout.

## Installing PEATSA ##

PEATSA has a number of dependencies on open-source third party packages which must be installed before beginning:

  * numpy
  * python 2.5
  * Boost C++ libs
  * pdb2pqr
  * APBS with python wrappers enabled

In addition pdb2pqr can be replaced with WHAT IF and APBS with DelPhi if desired. However note that WHAT IF and DelPhi are not open-source and licenses for them must be obtained separately.
PEATSA also depends on a number of other programs in the PEAT suite:
UFFBAPS, FFF, pKaTool

These packages must be manually compiled for PEATSA to work. All are part of the PEAT suite and are hence their source code is present when you check out our repository.
PEATSAParallel describes how to get PEATSA running in a parallel environment.

### Do it yourself ###
These are the instructions for installing PEATSA from scratch yourself.
First Steps
First make sure you have all the required third-party dependencies listed in InstallationNotes. This is probably the most difficult step.

Getting the software
Run the following command:
```
svn co http://peat.googlecode.com/svn/trunk/ PEAT
```
Building
```
cd to PEAT/UFFBAPS and run make
cd to PEAT/pKaTool/ and run make -f Makefile.Linux
cd to PEAT/FFF and
make depends
make
```

Finally
Two options -
Add PEATSA/Core/ to your $PATH
Create an alias to PEATSA/Core/ProteinDesignTool.py in a directory in your path.
Now your ready to move on to PreparingARun.

## The Configuration File ##

Before running PEATSA you need to setup your PEATSA ConfigurationFile This usually only has to be done once ever. Amongst other things the configuration file tells PEATSA where it can find its component programs - UFFBAPS, pKD - and which programs to use.
Whenever you run PEATSA it looks for a proteinDesignTool.conf file first in the directory you ran the program from, then in your $HOME directory. Otherwise it falls back on default values. So to get the program to work you need to create a personal instance of this file.
This can be done by running PEATSA as follows -
ProteinDesignTool --create-configuration
This will create a file called proteinDesignTool.conf in the directory where you ran the above command . The file is divided into sections each containing a number of //option=value// pairs. You need to set two of these.
Here's an example of the file:

```
[PKA SCAN METHOD]
MCsteps = 0 
dpKa_method = MC
tabulated = 1 
PBEsolver = APBS

[PKA SCAN PARAMETERS]
recalc_intpka_dist = 10
generate_mutations = False
save_temp_files = False
pHstep = 0.01
use_titration_curves = 1 
calc_dpka = 1 
pKMCsteps = 200000
pHstop = 12.0
recalc_intpka = 1 
pHstart = 0.1 
mutation_quality=0.5

[PKA SCAN OUTPUT]
save_solutions = None
verbose = 1 

[PATHS]
pKa Tools = /home/johnston/PEAT/
uffbaps = /home/johnston/PEAT/UFFBAPS/

[WORKING DIRECTORY]

copyIfExists = 1 
copyLocation = /home/johnston/PEATSACopies/
overwriteExistingCopies = 0 
useUniqueID = 1 

[PARALLEL]
executable = /home/johnston/PEAT/PEAT_SA/Core/ProteinDesignTool.py 
queue = medium
numberOfNodes = 4 
processorsPerNode = 4 
runDirectory = . 
logDirectory = . 
```

The only options that need to be set are those in the PATHS section. Here you specify the locations of PEAT\_SA's component programs - pKD (the pKaTool option) and UFFBAPS. The above example contains the paths as they would appear if you have installed PEAT\_SA as described in the InstallationNotes.
The other three sections, PKA SCAN METHOD, PKA SCAN PARAMETERS and PKA SCAN OUTPUT, are options that will be passed to the pKD program - see its documentation for more about these options.

## Setting up WHAT IF ##

If you are using WHATIF you need to create a small file detailing where it is.
Here are the steps:
  * Create a file called .WHAT\_IF in your home directory
  * Edit it and add two lines: The location of the WHAT\_IF executable, The location of the WHAT\_IF config file - this is called WHATIF.fig

Heres an example
```
/software/bin/DO_WHATIF.COM
/software/bin/
```

## Running a calculation ##

The PEATSA interface makes things easy for you because it handles the details of running the UFFBAPS and pKD programs for many protein mutants, gathers the results, and outputs them in an easy to use format.
NOTE: Before you can run PEATSA on a protein you need to created a valid configuration file.
Let's say you have a pdb file called 1crn.pdb you want to run the program on. The first thing to note is that PEATSA uses two directories for its output - the WorkingDirectory and the OutputDirectory. In this example I will assume the working directory is called $HOME/1crn and the output directory is the one where the command is being run. The working directory must already exist.
You can run PEATSA on this data by issuing the following command:
```
ProteinDesignTool -p 1crn.pdb -w $HOME/1crn --mutation=ALA --stability
```
This will run a stability alanine scan calculation for the specified PDB file. Note you don't have to do any copying of files or other setup tasks - PEATSA takes care of all this for you. See ResultsFiles for what you get back in each case.
PEATSA has a number of other command line options which you can use to refine its behaviour. To see these options do
```
ProteinDesignTool --help
```
Check the various CalculationTypes for more specific examples.
If PEATSA doesn't work for some reason file a Getting Help/Bugs or Getting Help/Bugs and we'll help you overcome the problem.

## PEATSA Output ##

### Log Files ###
PEATSA outputs some log files which contain details on the calculations it has performed. The exact files output will depend on the calculations run.
Normally most users do not have to be concerned with the log files. However they may interest advanced users and can provide important information when it comes to identifying what went wrong with a calculation.

Scan.log

This file is output for each Scan calculation in the directory where the command was run. Note that subsequent runs in the same directory will overwrite any previously existing Scan.log. This file contains information on the progress of the pKD algorithms.
Results Files
For each type of calculation you can run - scan, stability, binding - you get one results file. Each results file is in csv (comma separated value) format which can be read by all spreadsheet programs and is best way to view these files.

The files are -

  * ScanResults.csv
  * StabilityResults.csv
  * BindingResults.csv

They are output to a directory by default called '$(PDBId).peatsa/' where $(PDBId) is the name of the file containing the protein structure (without extension) passed to the program. A different name can be specified using the '-n' command line option.
You also get another directory called '$(PDBId).mutants' which contains the mutant models PEATSA created. The top level of this directory contains two files 'Mutate.log' and 'Scores.csv'. The first of these contains information on the modelling process and the second contains the bump-scores for each mutant PEATSA attempted to model. Also in this directory are two sub-dirs - Mutants and Ligands. The first contains the mutant pdb files, the second any ligand files that were also passed to the program.

### Stability results ###
See Stability for information on the calculation that outputs this file.
The stability results file is a table in csv format and is always called StabilityResults.csv. Each table row corresponds to the mutation of one residue to alanine and tells you how this mutation changes the stability (the free-energy difference between the wild-type and mutant) of the protein.
The free-energy difference is the sum of the contributions from a number of different factors. These are -

  * Van der Waals
  * Electrostatic
  * Hydrogen-Bond - abbreviated H-Bond
  * Desolvation
  * Backbone

The energies are in KJ.mol-1. More information on these terms can be found in the UFFBAPS description.
An example output file - taken from a run on 1crn.pdb - is shown below
```
#  
# Generated  by PEAT_SA on Thursday 19/06/2008 12:25:44
# Values are in KJ.mol-1
# This file is best viewed in a spreadsheet program e.g. Excel
#
Chain Id, Residue Number, Residue Name, Van der Waals, Electrostatic, H-Bond,  Desolvation, Backbone Entropy, Total
T, 0001, THR, 1.57, 0.0, 0.0, -0.94, 0.0, 0.631
T, 0002, THR, 1.84, 0.0, 0.0, 0.04, 1.26, 3.133
T, 0003, CYS, -3.57, 0.0, 0.0, 1.95, 1.52, -0.1
T, 0004, CYS, -4.97, 0.0, 0.0, 3.99, 2.79, 1.812
T, 0005, PRO, 6.19, 0.0, 0.0, -0.81, -1.95, 3.429
```

The first three row entries define the residue that has been mutated to alanine. For example the first row corresponds to THR1 in chain T. The next columns give the contribution of each factor to the free-energy difference between the original structure and the alanine mutant (THR 1 -> ALA). The file is best viewed in a spreadsheet program as the relationship between these values and the column headers is much more obvious.

### Binding results ###

Notes

**Problems with mol2 files**

Currently there are a few formatting issues when using mol2 files for ligand binding calculations with PEAT-SA. Chimera creates a mol2 file that begins like this:
```
@<TRIPOS>MOLECULE
lysmgm_ts_2b.pdb
2067 2089 130 0 0
PROTEIN
AMBER ff03.r1
To run properly you need to remove the 2nd line so it looks like:
@<TRIPOS>MOLECULE
2067 2089 130 0 0
PROTEIN
AMBER ff03.r1
Yasara on the other hand (or at least some versions of it) doesn't add a Substructure section at the end of the file which should look something like this:
@<TRIPOS>SUBSTRUCTURE
    1 NDG501      1 RESIDUE           4 A     NDG     1 ROOT
    2 NAG502     16 RESIDUE           4 A     NAG     2
    3 NAG503     30 RESIDUE           4 A     NAG     2
    4 NAG504     44 RESIDUE           4 A     NAG     1
```

Here is an example of a correct file. For more info on mol2 files generally see http://tripos.com/mol2/mol2_format3.html