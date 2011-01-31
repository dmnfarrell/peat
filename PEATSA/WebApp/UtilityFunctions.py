#! /usr/bin/env python
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
'''Utility functions used by other WebApp classes'''
import os, MySQLdb, datetime, md5, time, random
import PEATSA.Core as Core

def ResourceDirectory():

	'''Returns the path to the resource directory for the Backend'''

	#Assign the location of the module to a variable
	moduleLocation = os.path.dirname(__file__)
	moduleLocation = os.path.abspath(moduleLocation)

	#Create the resource path
	resources = os.path.join(moduleLocation, 'Resources')
		
	return resources
	
def ExecutablePath():

	'''Returns the absolute path to PDTWebScript.py.
	
	This is the backend executable for the WebApp which actually runs the job'''
	
	moduleLocation = os.path.dirname(__file__)
	moduleLocation = os.path.abspath(moduleLocation)	
	
	executablePath = os.path.join(moduleLocation, 'WebScript.py')
	
	return executablePath

def CreateDateString():

	'''Creates a date string in the format yyyy-mm-dd hh:mm:ss
	
	This format is acceptable for insertion as a DATETIME value in SQL'''

	date = datetime.datetime.today()
	string = "%s-%s-%s %s:%s:%s" % (date.year, date.month, date.day, date.hour, date.minute, date.second)
	return string

def GetNewSID(tag):

	'''Returns a unique string for use as a Session ID'''
	
	t1 = time.time()
	t2 = t1 + random.random()
	base = md5.new( tag + str(t1 +t2) )
	sid = tag + '_' + base.hexdigest()
	return sid
	
def ConnectionFromConfiguration(configuration):

	'''Creates and returns a connection based on the information in the configuration object.
	
	The configuration object should be an instance of the Core.Environment.Configuration class.
	The underlying configuration file should contain a section 'DATABASE' with the value
	host, user, password and database'''

	connection = MySQLdb.connect(host=configuration.get("DATABASE", "host"),
				 user=configuration.get("DATABASE","user"),
				 passwd=configuration.get("DATABASE", "password"), 
				 db=configuration.get("DATABASE", "database"))	
				 
	return connection			 
				 
def ConnectionFromDefaultConfiguration():

	'''Create and returns a connection based on the information in the default webapp configuration file.
	
	This is located in Backend/Resources/webApplication.conf'''

	return ConnectionFromConfiguration(DefaultWebAppConfiguration())	
	
def DefaultWebAppConfiguration():
	
	'''Returns the default configuration object based on the default configuration file.	
	
	This is located in Backend/Resources/webApplication.conf'''
	
	configurationFile = os.path.join(ResourceDirectory(), 'webApplication.conf')	
	return Core.Environment.Configuration(filename=configurationFile)

def SendNotificationEmail(job):

	'''Sends an email to the address associated with job is done

	Note: A mail is only sent if no previous mail has been sent.

	Parameters:
	job - A WebApp.Data.Job instance representing a PEATSA job

	'''

	if job.isEmailSent():
		return

	import smtplib
	from email.MIMEText import MIMEText
	from email.MIMEMultipart import MIMEMultipart

	message = MIMEMultipart()
	message['From'] = 'peatadmin@ucd.ie'
	message['Reply-To'] = 'peatadmin@ucd.ie'
	message['To'] = job.email()

	if job.error() is None:
		message['Subject'] = 'PEATSA Job Finished'
		text = "Your job has finished.\n\n"
		text = text + "See http://peat.ucd.ie/PEATSA/Pages/Results.php?jobId=%s" % job.identification
	else:
		message['Subject'] = 'PEATSA Job Error'
		text = "Unfortunately we encountered an error when running your job.\n\n"
		text = text + "Please see http://peat.ucd.ie/PEATSA/Pages/Results.php?jobId=%s for more information." % job.identification
			
	text = text + "\n\nRegards,\n\tThe PEAT-SA developers."

	bodyText = MIMEText(text)
	message.attach(bodyText)

	# Send the email
	try:
		server=smtplib.SMTP('mail.ucd.ie')
		server.sendmail('peatadmin@ucd.ie',job.email(),message.as_string())
		server.quit()  
	except Exception, data:
		print 'Exception sending emaili to', job.email()
		print data

	job.setEmailSent()
	
