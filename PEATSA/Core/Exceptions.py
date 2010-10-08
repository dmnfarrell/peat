#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
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
import Utilities

'''Defines PEAT-SA specific exceptions'''

class ProteinDesignToolException(BaseException):
	
	'''Generic exception representing a problem with the tool'''
	
	errorString = "Detected a generic error from the program - this is probably a bug"
	
class EnvironmentError(ProteinDesignToolException):

	'''Exception raised when something is wrong with the environment for running a calculation'''

	errorString = "Detected an error in the tools environment"
	
class ConfigurationError(EnvironmentError):

	'''Raised when there is a problem with the configuration of the run'''
	
	errorString = "Detected an error with the tools configuration" 
	
class WorkingDirectoryError(EnvironmentError):

	'''Raises when there is a problem with the working directory specified for a run'''
	
	errorString = "Detected an error with the tools working directory"
	
class ScanError(ProteinDesignToolException):

	'''Raised when something goes wrong with the scan process'''

	errorString = "Encountered an error while attempting to perform an alanine scan"
	
	pass
	
class MissingPKADataError(WorkingDirectoryError):

	'''Raised when a scan setup failed because the initial pKa data is missing'''
	
	errorString = "Could not setup up scan run because required pKa data was missing"
	
class StabilityError(ProteinDesignToolException):

	'''Raised when something is goes wrong with a stability calculation'''
	
	errorString = "Encountered an error while attempting to perform stability calculation"
	
	pass
	
class ArgumentError(ProteinDesignToolException):

	'''Raised when something is wrong with the command line args'''
	
	errorString = "Encountered a problem with the command line arguments provided"
	
class MutantCollectionError(ProteinDesignToolException):

	'''Raised when there is a problem creating a MutantCollection instance'''
	
	errorString = "Encountered an error when creating a MutantCollection"
	
class MutantModellingError(ProteinDesignToolException):

	'''Raised when there is a problem modelling a mutant'''
	
	errorString = "Encountered an error when modelling a mutant"	
	
class DataAlreadyExistsError(ProteinDesignToolException):

	'''Raised when an attempt is made to add already existing data to a PEAT-SA output directory'''
	
	errorString = "Encountered a problem when adding data to the output directory"
	
class UnableToRemoveDataError(ProteinDesignToolException):

	'''Raised when an attempt to delete an existing PEATSA output data directory fails'''
	
	errorString = "Encountered a problem when attempting to remove the data directory"
	
class FileFormatError(ProteinDesignToolException):

	'''Raised when something is wrong with the format of one of the files supplied'''
	
	errorString = "Encountered a problem with the format of an input file"	
	
class MutationListError(ProteinDesignToolException):

	'''Raised when there is a generic error with a mutation list file'''
	
	errorString = "Encountered a problem with a mutation list"
	
class MutationListEntryError(MutationListError):

	'''Raised when an element cannont be added to the mutation list'''
	
	errorString = "Could not add supplied mutant to the mutation list"
	
class MutationListFileFormatError(MutationListError):

	'''Raised when there is the format of a supplied mutation list file is invalid'''
	
	errorString = "Detected problem with the format of the mutation list file"
	
class MutationListDuplicatEntryError(MutationListError):

	'''Raised when two instances of the same mutant are in a list'''
	
	errorString = "Detected a dupicate entry in the list"
	
class MutationListValidationError(MutationListError):

	'''Indicated some of the mutants in a mutation list cannot be applied to a given pdb'''
	
	errorString = "Failed to validate mutation list against supplied pdb"

class ResidueCodeFormatError(ProteinDesignToolException):

	'''Raises when an supplied residue code could not be parsed into its components'''
	
	errorString = "Encountered error when attempting to parse a residue code"
	
class MutationCodeFormatError(ProteinDesignToolException):

	'''Raises when an supplied mutation code could not be parsed into its components'''
	
	errorString = "Encountered error when attempting to parse a mutation code"
	
class MutationCodeApplicationError(ProteinDesignToolException):

	'''Raised when the mutation defined by a mutation code cannot be applied to a pdb 
	
	This is because the chain or residue index defined by the code does not 
	exist in the pdb the code is being applied to'''
	
	def __init__(self, mutationCode, message=None):
	
		components = Utilities.ParseMutationCode(mutationCode)
		if len(components) == 4:
			chain, residueIndex, wt, mutation = components 
		else:
			chain, residueIndex, mutation = components
				
		self.residueCode = Utilities.CreateResidueCode(chain, residueIndex)
		self.mutationCode = mutationCode
		self.errorString = "An error was encountered when searching for residue %s in the pdb structure -\n" % self.residueCode
		if message==None:
			self.errorString = self.errorString +  "Specified residue does not exist"
		else:
			self.errorString = self.errorString + message	

		print self.errorString
	
	
	
	
	
