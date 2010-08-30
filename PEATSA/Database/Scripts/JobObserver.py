#! /usr/bin/python
'''JobObserver monitors the queue position of jobs that have been queued but are not running.'''
import sys, time, optparse, os
import PEAT_SA.Core as Core
import PEAT_SA.WebApp as WebApp
import subprocess

def getQueueStatistics(id):
	
	'''Returns a dictionary of queue statistics for job with id.

	Note id must be a pbs id (not a PDT server id'''

	command = 'qstat -f1 %d' % id
	task = subprocess.Popen(args=command, stdout=subprocess.PIPE, shell=True)
	task.wait()
	components = task.stdout.readlines()
	components.pop(0)
	components.pop()
	components = [el.split("=") for el in components]
	stats = {}
	for data in components:
		data = [el.strip() for el in data]
		stats[data[0]] = data[1]

	return stats

def getJobPositions(queue):

	print queue
	noneString = 'no idle jobs in queue\n'
	command = 'diagnose -p %s' % queue
	task = subprocess.Popen(args=command, stdout=subprocess.PIPE, shell=True)
	task.wait()
	components = task.stdout.readlines()

	queuedJobs = []
	if components[0] == command:
		pass
	else:
		#Strip out non-job lines
		components = components[5:-5]
		ids = []
		for component in components:
			data = component.split()
			ids.append(int(data[0].strip()))
				
		for id in ids:
			stat = getQueueStatistics(id)	
			queuedJobs.append(os.path.splitext(stat['Job_Name'])[0])

	return queuedJobs


if __name__ == "__main__":

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-c", "--configurationFile", dest="configurationFile",
			  help="A PEAT-SA configuration file." 
			  "This must contain the queue attribute in the the PARALLEL section", 
			  metavar="CONF")
			  
	(options, args) = parser.parse_args()

	if options.configurationFile is None:
		print 'Configuration file must be specified'
		sys.exit(1)
	else:
		options.configurationFile = os.path.abspath(options.configurationFile)

	configuration = Core.Environment.Configuration(filename=options.configurationFile)
	jobTable = configuration.get('DATABASE', 'jobTable')
	queue = configuration.get('PARALLEL', 'queue')

	print 'Starting observer daemon'

	while(1):
		#Connection might not be valid when we get back here each loop so have to recreate it
		connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
		jobManager = WebApp.Data.JobManager(connection=connection, jobTable=jobTable)
		queuedJobs = jobManager.jobsInState('Ready')
		print 'Queued jobs'
		print queuedJobs 
		positions = getJobPositions('servers')
		print 'Queued jobs positions', positions
		#Set a status message for jobs in the queue
		#Jobs that leave the queue will have their status message updated
		#by WebApp.py
		for id in queuedJobs:
			try:
				position = positions.index(id)
			except ValueError:	
				position = -1

			#There are two cases where the job may not be present
			#The job was started since we retrieved its status
			#The job hasn't been enqueued yet
			if position == -1:
				if runningJobs.count(id) != -1:
					job.setQueueStatus('Your job is running')
				else:
					job.setQueueStatus('Your job is about to be enqueued')
			if position == 0:
				job.setQueueStatus('Your job is next to be run')
			elif position == 1:
				job.setQueueStatus('There is one job ahead of you in the queue')
			else:
				job.setQueueStatus('There are %d jobs ahead of you in the queue', position)

		connection.close()
		time.sleep(30)
			
