#!/usr/bin/env python
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
'''Cgi Script for dynamically creating a mutant pdb file used in a PEAT-SA job'''
import cgi, cgitb, tempfile, os, sys
import PEATSA.WebApp as WebApp
import PEATSA.Core as Core

#URL'S have query string jobId=$jobID&mutationCode=$mutationCode	
if __name__ == "__main__":

	#Turn on exception logging
	cgitb.enable()
	formData = cgi.FieldStorage()

	#Redirect stdout
	oldout = sys.stdout
	newout = tempfile.NamedTemporaryFile()
	sys.stdout = newout
	
	#Get Data
	try:
		#Set environment
		configuration = WebApp.UtilityFunctions.DefaultWebAppConfiguration()
		os.environ['HOME'] = configuration.get('PATHS', 'whatif')
		
		jobId = formData.getvalue('jobId')
		mutationCode = formData.getvalue('mutationCode')
			
		#Create Mutation set
		mutationSet = Core.Data.MutationSet(code=mutationCode)
				
		#Get the pdb
		connection = WebApp.UtilityFunctions.ConnectionFromDefaultConfiguration()
		job = WebApp.Data.Job(jobID=jobId, connection=connection)
		structure = job.structure()
		
		#Write it out
		temp = tempfile.NamedTemporaryFile()
		name = temp.name
		temp.write(structure)
		temp.flush()

			
		#Create Mutant
		collection = Core.Data.MutantCollection(pdbFile=name, 
					location='/tmp/', 
					name=name, 
					maxOverlap=0.5,
					temporary=False)
		collection.createMutants(mutationList=[mutationSet]) 

		if len(collection.mutantFiles()) == 0:
			raise Exception, "Modelling Failure - %s" % temp.name
	
		#Get PDB
		filename = collection.fileForMutationSet(mutationSet)
		stream = open(filename)
		string = stream.read()
		stream.close()

		
		#Clean up
		temp.close()
		sys.stdout = oldout
		newout.close()
		
		#Download
		print "Content-length: %d", len(string)
		print "Content-type: application/pdb"
		print "Content-Disposition: attachment; filename=%s\n\n" % mutationSet.filename(reduced=True)
		print string
	except Exception, data:
		temp.close()
		sys.stdout = oldout
		newout.close()
		url = WebApp.JobSubmission.ConstructErrorURL(domain="MutantCreationErrorDomain",
			description="Encounter a %s exception while attempting to create mutant" % Exception,
			detailedDescription="Exception data - %s. Mutant %s" % (data, mutationCode),
			recoverySuggestion="This may be bug - please report to site administrators")
		#Redirect to the webpage indicated by resultURL
		print "Location: %s\n\n" % url
	
		
				
