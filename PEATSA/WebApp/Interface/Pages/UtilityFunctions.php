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

	//Contains various UtilityFunctions used by the PEAT-SA web-app
	include_once('Config.php');
	include_once 'PEATSA/Database/InterfaceKit/DataRetrieval.php';

	function send_job_completion_mail($jobID)
	{
		include('Mail.php');
		
		$error = 0;
		$email = get_email_for_job($jobID, $error);
		//Don't do anything if we can't get an email
		if($error == 1)
		{
			return;
		}
		
		//Check if we've already sent it
		$sent = is_email_sent($jobID, $error);
		if($error == 1)
		{
			//We failed to find out if this was already sent or not
			//Ignore for now
			return;
		}		
		
		//It will be unknown if the user didn't put a email
		if(($email != "Unknown") && !$sent)
		{
			$smtp = Mail::factory('smtp', array ('host' => 'mail.ucd.ie', 'auth' => false));
			
			//admin
			$headers = array("From"=>"peatadmin@ucd.ie", "Subject"=>"PEATSA Job Finished");
			$body = "Your job has finished\n\nSee http://enzyme.ucd.ie/PEATSA/Pages/Results.php?jobId=$jobID\nRegards,\n\nThe PEAT-SA developers.";
			$mail = $smtp->send($email, $headers, $body);
			set_email_sent($jobID);
		}
	}
	
	function send_job_error_mail($jobID)
	{
		//notify admin that someone has registered and send user a mail to verify address.
		include('Mail.php');
		
		$error = 0;
		$email = get_email_for_job($jobID, $error);
		
		//Don't do anything if we can't get an email
		if($error == 1)
		{
			return;
		}
		
		//Check if we've already sent it
		$sent = is_email_sent($jobID, $error);
		if($error == 1)
		{
			//We failed to find out if this was already sent or not
			//Ignore for now
			return;
		}
		
		//It will be unknown if the user didn't put a email
		if(($email != "Unknown") && !$sent)
		{
			$smtp = Mail::factory('smtp', array ('host' => 'mail.ucd.ie', 'auth' => false));
			
			//admin
			$headers = array("From"=>"peatadmin@ucd.ie", "Subject"=>"PEATSA Job Error");
			$body = "Unfortunately we encountered an error when running your job\n\nPlease see http://enzyme.ucd.ie/PEATSA/Pages/Results.php?jobId=$jobID for more information,\nRegards,\n\nThe PEAT-SA developers.";
			$mail = $smtp->send($email, $headers, $body);
			set_email_sent($jobID);
		}
	}

	function create_error_url_from_array($data)
	{
		$query = http_build_query($data, '&amp');
		//$parts = array("scheme" => "http",
		//				"host" => "localhost",
		//				"port" => "8888",
		//				"path" => "Pages/Error.php",
		//				"query" => $query);
		$server = $_SERVER['SERVER_NAME'];
		$port = $_SERVER['SERVER_PORT'];
		$errorURL = "http://$server:$port/PEATSA/Pages/Error.php?$query";
		
		return $errorURL;
	}
	
	function create_error_url($jobId, $domain, $description, $detailedDescription, $recoverySuggestion)
	{
		$data = array("jobId" => $jobId,
			      "domain" => $domain,
			      "description" =>$description,
			      "detailedDescription" => $detailedDescription,
			      "recoverySuggestion" => $recoverySuggestion);
		
		$errorURL = create_error_url_from_array($data);
		return $errorURL;
	}
	
	
	function amino_acids()
	{
		return array("ALA", "ASP", "ARG", "ASN", "CYS", "GLY", "GLU", "GLN",
			     "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", 
			     "SER", "TYR", "TRP", "THR", "VAL");
	}
?>
