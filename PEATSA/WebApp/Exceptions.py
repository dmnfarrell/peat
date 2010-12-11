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
'''Contains WebApp specific exceptions'''
import PEATSA.Core as Core

class DatabaseRetrievalError(Core.Exceptions.ProteinDesignToolException):

	'''Raised when something happens when trying to retrive data from the web data base'''
	
	errorString = "Encountered a problem when attempting to retrieve database data"
	
class SubmissionException(Core.Exceptions.ProteinDesignToolException):

	'''Generic error with the submitted data'''
	
	errorString = "Error occured while submitting job"

class FileFormatError(Core.Exceptions.ProteinDesignToolException):

	'''Something wrong with a file submitted to the web-app'''
	
	errorString = "Error with submitted file"
