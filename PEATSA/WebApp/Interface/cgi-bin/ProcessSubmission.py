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
'''Cgi Script executed when a user submits a job request from the WebApp main page'''
import cgi, cgitb
import PEATSA.WebApp as WebApp

#This script processes the input data submitted from FrontPage.php
#The result in every case is an url.
#This will either be to Results.php if the job was launched or Error.php otherwise	
if __name__ == "__main__":

	#print 'ContentType: plain/text\n'
	#Turn on exception logging
	cgitb.enable()

	#The FormData class process and checks the form data.
	formData = WebApp.JobSubmission.FormData(cgi.FieldStorage())
	
	#If the data submitted by the user is ok we proceed to run the job
	if formData.error() is not None:
		resultURL = formData.errorURL()
	else:	
		#The JobConstructor class sets up everything required to run a PEATSA calc 
		#based on the information provided by the FormData instance.
		#This involves
		#	- adding an entry for the job to the WebApp db 
		#	- writing out the uploaded files to the correct locations
		#	- constructing the command line for the backend executable.
		jobConstructor = WebApp.JobSubmission.JobConstructor(formData, construct=True)
		
		#If nothing went wrong with setting up the job then we can run the backend
		#The backend will be either WebApp.WebScript.py or WebApp.RemoteLaunch.py
		if jobConstructor.error() is None:	
			resultURL = jobConstructor.runBackend()
		else:
			resultURL = jobConstructor.errorURL()
			
	#Redirect to the webpage indicated by resultURL
	print "Location: %s\n\n" % resultURL
	
