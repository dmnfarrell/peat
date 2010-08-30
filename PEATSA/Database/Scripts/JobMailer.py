#! /usr/bin/python
import sys, time, optparse, os
import PEATSA.Core as Core
import PEATSA.WebApp as WebApp
import ConstructJob
import MySQLdb

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-c", "--configurationFile", dest="configurationFile",
                  help="A PEATSA configuration file.", metavar="CONF")
		  
(options, args) = parser.parse_args()

if options.configurationFile is None:
        print 'Configuration file must be provided'
        sys.exit(1)

configuration = Core.Environment.Configuration(filename=options.configurationFile)
jobTable = configuration.get('DATABASE', 'jobTable')

while(1):
	connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
	jobManager = WebApp.Data.JobManager(connection=connection, jobTable=jobTable)
	selectQuery = """SELECT JobId FROM %s WHERE SentMail='0' AND State='Finished' AND NOT Email='Unknown'""" % (jobTable)
	cursor = connection.cursor()
	cursor.execute(selectQuery)		
	ids = [el[0] for el in cursor.fetchall()]
	for id in ids:
		job = WebApp.Data.Job(id, connection)
		print 'Sending mail for job %s to %s' % (job.identification, job.email())
		WebApp.UtilityFunctions.SendNotificationEmail(job)
	connection.close()
	time.sleep(30)
