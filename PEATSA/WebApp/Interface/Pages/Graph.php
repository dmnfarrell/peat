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

	/*
 	 * This page retrieves an image for a specified calculation of a job
	 * and creates a page for it.
	 * The required URL query data is the job id and the name of the calculation
	 * whose results are to be retrieved.
	 *
	 * The error handling is done in much the same way as Results.php - 
	 * Check there for more details.
	 */

	//Check required data is present before doing anything
	//This has to be the very first thing in the file for it to work.
	include_once 'UtilityFunctions.php';
	include_once 'PEATSA/Database/InterfaceKit/DataRetrieval.php';

	//First check that a job id was passed
	if(!array_key_exists('jobId', $_GET))
	{
			$recoverySuggestion = 'Did you try to navigate directly to Graph.php?<br>';
			$recoverySuggestion = $recoverySuggestion."If you haven't submitted a job then go to the main page first.<br><br>";

			//Redirect to the error page
			$errorURL = create_error_url($_GET['jobId'],
							'PDT.NavigationErrorDomain',
							'Required data is missing.',
							'No job identification present in request to Graph.php.',
							$recoverySuggestion);
			header("Location: $errorURL \n\n");
			exit();
	}
	
	//Check what type of data was selected
	$calculation = $_GET["calculation"];
	$outputName = $calculation."Results.png";
	
	/*
	 * Next retrieve the data
	 */

	$error = 0;
	$jobId = $_GET['jobId'];
	//This returns an errorURLif it fails
	//and an array containing two key, content and size, if it succeeds
	$data = get_image_for_calculation($jobId, $calculation, &$error);

	//If no error has been detected output the data
	if($error == 0)
	{
		//Prepare and output content for download
		$size = $data["size"];
		header("Content-length: $size");
		header("Content-type: image/png");
		echo $data["content"];
		exit();
	}
	else
	{
		$errorURL = create_error_url_from_array($data);
		//Redirect to the error page
		header("Location: $errorURL \n\n");
		exit();
	}
?>

