#! /usr/bin/python
import sys, time, optparse, os
import PEATSA.Core as Core
import PEATSA.WebApp as WebApp
import ConstructJob

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-c", "--configurationFile", dest="configurationFile",
                  help="A PEATSA configuration file.", metavar="CONF")
parser.add_option("-d", "--directory", dest="dir",
                  help="Where to write to job data.", metavar="DIR")
		  
(options, args) = parser.parse_args()

if options.configurationFile is None:
        print 'Job id must be provided'
        sys.exit(1)
else:
	options.configurationFile = os.path.abspath(options.configurationFile)

if options.dir is None:
        print 'No output directory specified - default to current directory'
	options.dir = os.getcwd()
else:
	options.dir = os.path.abspath(options.dir)

if not os.path.exists(options.dir):
	os.mkdir(options.dir)

configuration = Core.Environment.Configuration(filename=options.configurationFile)
jobTable = configuration.get('DATABASE', 'jobTable')
stream = open(os.path.join(options.dir, 'Run.log'), 'w+')
sys.stdout = stream

print 'Starting run daemon - run directory %s' % options.dir
stream.flush()

while(1):
	#Connection might not be valid when we get back here each loop so have to recreate it
	connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
	jobManager = WebApp.Data.JobManager(connection=connection, jobTable=jobTable)
	waitingJobs = jobManager.jobsInState('Ready')
	for id in waitingJobs:
		outputDirectory = os.path.join(options.dir, id)
		print 'Constructing job ', id
		#Set it to use the configuration we're passed here rather than the default config
		pbsFile = ConstructJob.ConstructJob(id, outputDirectory, configuration, connection, jobConfigurationFile=options.configurationFile)
		print 'Submitting %s' %  pbsFile
		job = WebApp.Data.Job(id, connection)	
		job.setState('Launched')
		os.system('qsub %s' % pbsFile)
		print 'Done'
	stream.flush()	
	connection.close()
	time.sleep(30)
			
