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

	//Contains functions for enabling Database/InterfaceKit to be used by different
	//applications on the same server (which can then access different databases).
	//
	//Also allows same server to access multiple different databases depending on user preferences.
	//The only two functions for use outside are *_peatsa_database_*
	//
	//Also contains the funtion create_error_array since its needed here
	//
		
	include_once('Config.php');
	
	function create_error_array($jobId = 'None', $domain, $description, $detailedDescription, $recoverySuggestion)
	{
		$data = array("jobId" => $jobId,
			      "domain" => $domain,
			      "description" =>$description,
			      "detailedDescription" => $detailedDescription,
			      "recoverySuggestion" => $recoverySuggestion);
		
		return $data;
	}
	
	//Starts a PEAT_SADB session for visitors to $domain
	//All vistors to $domain will use the same db session info
	function start_peatsa_database_session($domain, $time = NULL)
	{
		//See if the user already has a session running
		$params = session_get_cookie_params();
		
		if(empty($params['domain']))
		{
			//No session data present - send a cookie to the user for $domain.
			//Anytime they return to $domain on this server this cookie will be loaded
			session_set_cookie_params(14*24*60*60, $domain);
		}
		
		//Looks for a session stored with name PEAT_SADB on the users host that
		//is also from $domain.
		if(session_id() == "")
		{
			session_name('PEAT_SADB');
			session_start();
		}
	}

	//Sets the session configuration file to $filename
	//This will be the configuration returned by get_info if no_other arguments are passed
	function set_peatsa_database_session_configuration_file($filename, $override = 0)
	{
		if(empty($_SESSION["configurationFile"]) || $override == 1)
		{
			//If no session is started this creates a session id
			//which will be either
			//1. Stored in Cookie on the users computer
			//2. Stored in the URL
			if(session_id() == "")
			{
				session_name('PEAT_SADB');
				session_start();
			}
			$_SESSION["configurationFile"] = $filename;
			session_commit();
		}
	}
	
	function session_configuration_file()
	{
		//This retreives the session id from the URL or a cookie.
		//and recreates the session.
		if(session_id() == "")
		{
			session_name('PEAT_SADB');
			session_start();
		}
			
		return $_SESSION["configurationFile"];
	}
	
	//Returns a configuraion object for $filename (which must be an IniFile format)
	//If $filename is NULL reads the configuration for the session configuration file
	function get_configuration($filename = NULL)
	{
		if($filename == NULL)
		{
			$filename = session_configuration_file();
		}
		
		$parser = new Config();
		$root = &$parser->parseConfig($filename, "IniFile");
		return $root;
	}
	
	//Extracts datatbase information contained in a configuration file to an array
	//The configuration file must have a specific format - see README
	//If no filename is provided, the default session configuration file is used (if defined)
	function get_database_info($filename = NULL)
	{
		$configuration = get_configuration($filename);		
		$databaseSection = $configuration->getItem("section", "DATABASE");
		
		$host = $databaseSection->getItem("directive", "host")->getContent();
		$database = $databaseSection->getItem("directive", "database")->getContent();
		$user =  $databaseSection->getItem("directive", "user")->getContent();
		$password = $databaseSection->getItem("directive", "password")->getContent();	
		
		return array("Host" => $host, 
			     "Database" => $database,
			     "User" => $user,
			     "Password" => $password);
	}

	//Creates a dsn array from the information in the configuration file
	//The configuration file must have a specific format - see README
	//If no filename is provided, the default session configuration file is used (if defined)
	function get_dsn_info($filename = NULL)
	{
		$configuration = get_configuration($filename);		
		$databaseSection = $configuration->getItem("section", "DATABASE");
		
		$host = $databaseSection->getItem("directive", "host")->getContent();
		$database = $databaseSection->getItem("directive", "database")->getContent();
		$user =  $databaseSection->getItem("directive", "user")->getContent();
		$password = $databaseSection->getItem("directive", "password")->getContent();	
		
		return array("hostspec" => $host, 
			     "database" => $database,
			     "username" => $user,
			     "password" => $password);
	}
	
	//Returns the session db paramters
	//This is an array containining the entries Host, User and Password
	//This function is essentially a wrapper arround get_database_info to increase readability
	function session_db_parameters()
	{
		return get_database_info();
	}

	function session_dsn_parameters()
	{
		$dsn = get_dsn_info();
		return $dsn;
	}
	
	//Returns a connection to the session db
	//Note: Each call to this function returns a unique connection
	function session_db_connection(&$error)
	{
		//Returns the db defined for this session
		$databaseFields = session_db_parameters();
		$connection = mysql_connect($databaseFields["Host"], $databaseFields["User"], $databaseFields["Password"]);
		
		if(!$connection)
		{
			$error = 1;
			//Can't connect to the db
			return create_error_array('PDT.DatabaseErrorDomain',
					'Error when accessing the database.',
					mysql_error(),
					"The database is probably down for maintainence. Please try again later.<br>`");
		}
		else
		{
			return $connection;
		}
	}
?>
