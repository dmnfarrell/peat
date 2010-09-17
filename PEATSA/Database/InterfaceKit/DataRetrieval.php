<?php
	/**
	  Protein Engineering Analysis Tool Structure Analysis (PEATSA)
	  Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
	 
	  Author: Michael Johnston
	 
	  This program is free software: you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published by
	  the Free Software Foundation, either version 3 of the License, or
	  (at your option) any later version.

	  This program is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	  GNU General Public License for more details.
	 
	  You should have received a copy of the GNU General Public License
	  along with this program.  If not, see <http://www.gnu.org/licenses/>.
	 
	  Contact information:
	  Email: Jens.Nielsen_at_gmail.com
	  Normal mail:
	  Jens Nielsen
	  SBBS, Conway Institute
	  University College Dublin
	  Dublin 4, Ireland
	**/

	//Contains functions for retrieving job data from a PEATSA database
	//An open connection for the database to access can be provided to each function
	//In addition parameters for a default connection (session connection) can be set.
	//If this is present, and no conection is passed to a function, it will use these parameters to
	//create a connection.
	//Every function returns an errorArray containing information on what happened
	include_once 'Session.php';
	include_once 'Utilities.php';
	
	function get_pdb_code_for_job($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}

		$databaseFields = get_database_info();
		$query = "SELECT PDBID FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
		
		return $row['PDBID'];
	}
	
	function is_email_sent($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		$databaseFields = get_database_info();
		$query = "SELECT SentMail FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$errorArray = create_error_array($jobId,
							 "PDT.DatabaseErrorDomain",
							 "Database query failed",
							 mysql_error(),	
							 "This could be a bug");
			$error = 1;
		}
		
		if($error == 1)
			return $errorArray;
		
		$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
		
		return $row['SentMail'];
	}
	
	//Returns True if setting succesful.
	//Otherwise returns an errorArray
	function set_email_sent($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		$databaseFields = get_database_info();
		$query = "UPDATE Jobs SET SentMail='1' WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		//queryId is True on success on False on failure
		if (!$queryId) 
		{ 
			$errorArray = create_error_array($jobId,
							 "PDT.DatabaseErrorDomain",
							 "Database query failed",
							 mysql_error(),	
							 "This could be a bug");
			$error = 1;
		}
		
		if($error == 1)
			return $errorArray;
		
		return $queryId;		
	}
	
	function get_email_for_job($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		$databaseFields = get_database_info();
		$query = "SELECT Email FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$errorArray = create_error_array($jobId,
						       "PDT.DatabaseErrorDomain",
						       "Database query failed",
						       mysql_error(),	
						       "This could be a bug");
			$error = 1;
		}
		
		if($error == 1)
			return $errorArray;
		
		$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
		
		return $row['Email'];
	}
	
	function get_queue_status_for_job($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		$databaseFields = get_database_info();
		$query = "SELECT QueueStatusMessage FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$errorArray = create_error_array($jobId,
							 "PDT.DatabaseErrorDomain",
							 "Database query failed",
							 mysql_error(),	
							 "This could be a bug");
			$error = 1;
		}
		
		if($error == 1)
			return $errorArray;
		
		$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
		
		return $row['QueueStatusMessage'];
	}
	
	function get_calculation_state_for_job($jobId, &$error, $connection = NULL)
	{
		//Connect to the database
		$error = 0;
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		//Perform the query
		$databaseFields = get_database_info();
		$query = "SELECT Stability, Scan, Binding FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);

		if (!$queryId) 
		{ 
			$errorArray = create_error_array($jobId,
						     "PDT.DatabaseErrorDomain",
						     "Database query failed",
						     mysql_error(),	
						     "This could be a bug");
			$error = 1;
		}
		
		if($error == 1)
			return $errorArray;
		
		//If the status query was successful, retreive and process the data.
		//Check if anything went wrong
		$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
		
		if($row['Scan'] == "")
		{
			//If the value of Scan isn't Selected or NotSelected something is wrong!
			//No entry for the job exists in the db - 
			//Create an error info array
			$error = 1;
			$errorArray = create_error_array($jobId,
						     "PDT.DatabaseErrorDomain",
						     'Unable to retreive job data',
						     "Job $jobId not present in database",	
						     "This could be a bug - Please contact us with the above information");
		}
		else
		{
			//If we are here an entry for the job was presented in the database - now get any error data
			//Retrieve ErrorDescription and DetailedDescription even though they may be empty
			//To save having to get them later if an error is found to have occurred
			$query = "SELECT Error, ErrorDescription, DetailedDescription FROM Jobs WHERE JobID = '$jobId'";
			mysql_select_db($databaseFields["Database"], $connection);
			$queryId = mysql_query($query, $connection);

			$errorInfo = mysql_fetch_array($queryId, MYSQL_ASSOC);
			
			//Check if there was an error
			if($errorInfo["Error"] == 1)
			{
				$error = 1;
				//Add linebreaks for any newlines in the error description strings
				$errorArray = create_error_array($jobId,
							     "PDT.DatabaseErrorDomain",
							     nl2br($errorInfo["ErrorDescription"]),
							     nl2br($errorInfo["DetailedDescription"]),	
							     "Possible error in data submitted to job");
			}
		}
		
		//Check the status of the last check
		if($error == 1)
		{
			return $errorArray;		
		}
		else
		{			
			return $row;
		}
	}	
	
	//Returns True if there is data for the specified calculation False if not
	//If there is a problem $error is set to 1 and this function returns an errorArray

	function data_exists_for_calculation($jobId, $calculation, &$error, $connection = NULL)
	{
		$error = 0;
		
		//Get calculation information
		if($calculation != "")
		{
			$name = $calculation."Results";
		}
		else
		{
			$error = 1;
			$errorData = create_error_array($jobId, 
							'PDT.NavigationErrorDomain',
							'Required data is missing.',
							'No calculation information supplied.',
							'');
			
			return $errorData;
		}	
		
		//Connect to the database
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
		
		if($error == 1)
			return $connection;
		
		//Get the data for the calculation
		$databaseFields = get_database_info();
		$query = "SELECT Size, Content FROM Data WHERE JobID = '$jobId' AND MatrixName = '$name'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if(!$queryId) 
		{ 
			$error = 1;
			$errorData = create_error_array($jobId, 
							'PDT.DatabaseErrorDomain',
							'Error when accessing the database',
							mysql_error(),
							'Possible incorrect database parameters');
			$retval = $errorData;
		}
		else if(mysql_num_rows($queryId) == 0)
		{
			$retval = 0;
		}
		else
		{
			$retval = 1;
		}
		
		//Close the connection to the database
		mysql_close($connection);
		return $retval;		
	}
		
	//Returns an array containing a csv formatted string for the data calculation and string size.
	//The array keys are 'content' and 'size'
	//If there is a problem $error is set to 1 and this function returns an errorArray
	function get_data_for_calculation($jobId, $calculation, &$error, $connection = NULL)
	{
		$error = 0;
		
		//Get calculation information
		if($calculation != "")
		{
			$name = $calculation."Results";
		}
		else
		{
			$error = 1;
			$errorData = create_error_array($jobId, 
					'PDT.NavigationErrorDomain',
					'Required data is missing.',
					'No calculation information supplied.',
					'');
			
			return $errorData;
		}	
		
		//Connect to the database
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
			
		if($error == 1)
			return $connection;
			
		//Get the data for the calculation
		$databaseFields = get_database_info();
		$query = "SELECT Size, Content FROM Data WHERE JobID = '$jobId' AND MatrixName = '$name'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$error = 1;
			$errorData = create_error_array($jobId, 
							'PDT.DatabaseErrorDomain',
							'Error when accessing the database',
							mysql_error(),
							'Possible incorrect database parameters');
			$retval = $errorData;
		}
		else
		{
			//Check the query returned something
			if(mysql_num_rows($queryId) == 0)
			{
				//No entry exists for this job 
				//Someone probably navigated directly to the page and either
				//A) The job wasn't finished.
				//B) The job id is invalid.
				$error = 1;
				$errorData = create_error_array($jobId, 
								'PDT.DatabaseErrorDomain',
								"No $calculation data exists for this job ID ($jobId).<br>"."Either the job isn't finished, a $calculation calculation was not run, or the job ID is invalid.",
								'Check job ID');
				
				$retval = $errorData;
			}
			else
			{
				$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
				$content = $row["Content"];
				$size = $row["Size"];
				$retval = array('content' => $content, 'size' => $size);		
			}
		}	

		//Close the connection to the database
		mysql_close($connection);
		return $retval;
	}
	
	//Returns an array containing svg image data for the calculation and image size.
	//The array keys are 'content' and 'size'
	//If there is a problem $error is set to 1 and this function returns an errorArray
	function get_image_for_calculation($jobId, $calculation, &$error, $connection = NULL)
	{
		$error = 0;
		
		//Get calculation information
		if($calculation != "")
		{
			$name = $calculation."Image";
		}
		else
		{
			$error = 1;
			
			$errorData = create_error_array($jobId, 
							'PDT.InvalidArgumentErrorDomain',
							"No calculation specified ($jobId).<br>"."A calculation to fetch the data for must be provided",
							"");
			return $errorData;				
		}
		
		//Connect to the database
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
				
		if($error == 1)
			return $connection;

		//Get the data for the calculation
		$databaseFields = get_database_info();
		$query = "SELECT Size, Content FROM Images WHERE JobID = '$jobId' AND Name = '$name'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$error = 1;
			$errorData = create_error_array($jobId, 
							'PDT.DatabaseErrorDomain',
							'Error when accessing the database',
							mysql_error(),
							'Possible incorrect database parameters');
			$retval = $errorData;
		}
		else
		{
			//Check the query returned something
			if(mysql_num_rows($queryId) == 0)
			{
				//No entry exists for this job 
				//Either
				//A) The job wasn't finished.
				//B) The job id is invalid.
				$error = 1;
				$errorData = create_error_array($jobId, 
								'PDT.DatabaseErrorDomain',
								'Error when accessing the database.',
								"No $calculation data exists for this job ID ($jobId).<br>"."Either the job isn't finished, a $calculation calculation was not run, or the job ID is invalid.",
								'Check job ID');
				$retval = $errorData;				
				
			}
			else
			{
				$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
				$content = $row["Content"];
				$size = $row["Size"];
				$retval = array('content' => $content, 'size' => $size);		
			}
		}	

		//Close the connection to the database
		mysql_close($connection);
		return $retval;
	}	
	
	//Returns an string containing PDB strucutre data for the job 
	//If there is a problem $error is set to 1 and this function returns an errorArray
	function get_structure_for_job($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		
		//Connect to the database
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}		
		
		if($error == 1)
			return $connection;
		
		//Get the data for the calculation
		$databaseFields = get_database_info();
		$query = "SELECT Structure FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$error = 1;
			//Can't connect to the db
			$errorData = create_error_array($jobId, 
							'PDT.DatabaseErrorDomain',
							'Error when accessing the database',
							mysql_error(),
							'Possible incorrect database parameters');
			$retval = $errorData;
		}
		else
		{
			//Check the query returned something
			if(mysql_num_rows($queryId) == 0)
			{
				//No entry exists for this job 
				//Someone probably navigated directly to the page and either
				//A) The job wasn't finished.
				//B) The job id is invalid.
				$error = 1;
				$errorData = create_error_array($jobId, 
								'PDT.DatabaseErrorDomain',
								'Error when accessing the database.',
								"No structure data exists for this job ID ($jobId).<br>"."Either the job isn't finished, or the job ID is invalid.",
								'Check job ID');
				
				$retval = $errorData;
			}
			else
			{
				$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
				$retval = $row["Structure"];
			}
		}	
		
		//Close the connection to the database
		mysql_close($connection);
		return $retval;
	}
	
	//Returns an string containing ligand strucutre data for the job 
	//If there is a problem $error is set to 1 and this function returns an errorArray
	function get_ligand_for_job($jobId, &$error, $connection = NULL)
	{
		$error = 0;
		
		//Connect to the database
		if($connection == NULL)
		{
			$connection = session_db_connection($error);
		}
				
		if($error == 1)
			return $connection;
		
		//Get the data for the calculation
		$databaseFields = get_database_info();
		$query = "SELECT Ligand FROM Jobs WHERE JobID = '$jobId'";
		mysql_select_db($databaseFields["Database"], $connection);
		$queryId = mysql_query($query, $connection);
		
		if (!$queryId) 
		{ 
			$error = 1;
			$errorData = array();
			//Can't connect to the db
			//Can't connect to the db
			$errorData = create_error_array($jobId, 
							'PDT.DatabaseErrorDomain',
							'Error when accessing the database',
							mysql_error(),
							'Possible incorrect database parameters');
			$retval = $errorData;
		}
		else
		{
			//Check the query returned something
			if(mysql_num_rows($queryId) == 0)
			{
				//No entry exists for this job 
				//Someone probably navigated directly to the page and either
				//A) The job wasn't finished.
				//B) The job id is invalid.
				$error = 1;
				$errorData = create_error_array($jobId, 
								'PDT.DatabaseErrorDomain',
								'Error when accessing the database.',
								"No ligand data exists for this job ID ($jobId).<br>"."Either the job isn't finished, a binding calculation was not run, or the job ID is invalid.",
								'Check job ID');
				$retval = $errorData;				
				
			}
			else
			{
				$row = mysql_fetch_array($queryId, MYSQL_ASSOC);
				$retval = $row["Ligand"];
			}
		}	
		
		//Close the connection to the database
		mysql_close($connection);
		return $retval;
	}
	
?>
