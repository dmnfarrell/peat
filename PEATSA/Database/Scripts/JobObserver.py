#! /usr/bin/python
'''JobObserver monitors the queue position of jobs that have been queued but are not running.'''
import sys, time, optparse, os
import PEAT_SA.Core as Core
import PEAT_SA.WebApp as WebApp
import subprocess

def getQueueStatisitcs(queue):

	'''Returns a array of dictionaries. 

	Each dictionary corresponds to a job in queue.
	The keys in each dictionary are 
		- JobId
		- Name
		- User
		- Time 
		- State 
	
	Due to restrictions of qstat output the job name is not available.
	User getJobStatistics - passing the JobId to get detailed stats.
	Note: The order of jobs in the returned dictionary correpsonds to
	the order of jobs in the queue'''

	command = 'qstat %s' % queue
	task = subprocess.Popen(args=command, stdout=subprocess.PIPE, shell=True)
	task.wait()
	components = task.stdout.readlines()
	#Discard first two lines - headers and separator
	components = components[2:]
	components = [el.split() for el in components]
	#The components are JobId, Name, User, Time. State and Queue
	#However only ~10 chars are output for name which isn't enough
	#Thereore we leave it out
	keys = ['JobId', 'User', 'Time', 'State']
	stats = []
	for data in components:
		data = [el.strip() for el in data]
		data.pop(1)
		data.pop()
		stats.append(dict(zip(keys, data)))

	return stats


def getJobStatistics(id):
	
	'''Returns a dictionary of PBS statistics for job with id.

	Note id must be a pbs job id e.g. 11021.enzyme'''

	command = 'qstat -f1 %s' % id
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
	
	'''Returns an FIFO ordered array of the names of jobs waiting in the given queue

	Note: The name is usually equivalent to the name of the jobs pbs file
	without the path extension'''

	stats = getQueueStatisitcs(queue)
	queuedJobsIds = [el['JobId'] for el in stats if el['State'] == 'Q']

	queuedJobNames = []
	for id in queuedJobsIds:
		stat = getJobStatistics(id)	
		queuedJobNames.append(os.path.splitext(stat['Job_Name'])[0])

	return queuedJobNames

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

	while(1):
		#Connection might not be valid when we get back here each loop so have to recreate it
		connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
		jobManager = WebApp.Data.JobManager(connection=connection, jobTable=jobTable)

		#This function returns an array of the ids of jobs in the queue.
		#The elements in the array are in FIFO order
		positions = getJobPositions(queue)

		#Set a status message for jobs in the queue
		#Jobs that leave the queue will have their status message updated by WebApp.py
		for position in range(len(positions)):
			id = positions[position]
			job = WebApp.Data.Job(id, connection)

			if position == 0:
				job.setQueueStatusMessage('Your job is next to be run')
			elif position == 1:
				job.setQueueStatusMessage('There is one job ahead of you in the queue')
			else:
				job.setQueueStatusMessage('There are %d jobs ahead of you in the queue' % position)

		connection.close()
		time.sleep(10)
			
