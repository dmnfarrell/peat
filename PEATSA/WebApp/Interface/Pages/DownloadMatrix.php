<?php
	/*
 	 * This page retrieves a results file for a specified calculation of a job
	 * and downloads it to the users computer.
	 * The required URL query data is the job id and the name of the calculation
	 * whose results are to be retrieved.
	 *
	 * The error handling is done in much the same way as Results.php - 
	 * Check there for more details.
	 */

	//Check required data is present before doing anything
	//This has to be the very first thing in the file for it to work.
	include 'PEATSA/Database/InterfaceKit/DataRetrieval.php';

	//First check that a job id was passed
	if(!array_key_exists('jobId', $_GET))
	{
			$recoverySuggestion = 'Did you try to navigate directly to Download.php?<br>';
			$recoverySuggestion = $recoverySuggestion."If you haven't submitted a job then go to the main page first.<br><br>";

			//Redirect to the error page
			$errorURL = create_error_url($_GET["jobId"],
							'PDT.NavigationErrorDomain',
							'Required data is missing.',
							'No job identification present in request to Download.php.',
							$recoverySuggestion);
			header("Location: $errorURL \n\n");
			exit();
	}
	
	//Check what type of data was selected
	$calculation = $_GET["calculation"];
	$outputName = $calculation."Results.csv";
	
	/*
	 * Next retrieve the data
	 */

	$error = 0;
	$jobId = $_GET['jobId'];
	//This returns an errorURLif it fails
	//and an array containing two key, content and size, if it succeeds
	$data = get_data_for_calculation($jobId, $calculation, &$error);

	//If no error has been detected output the data
	if($error == 0)
	{
		//Prepare and output content for download
		$size = $data["size"];
		header("Content-length: $size");
		header("Content-type: application/csv");
		header("Content-Disposition: attachment; filename=$outputName");
		echo $data["content"];
		exit();
	}
	else
	{
		$errorURL = $data;
		//Redirect to the error page
		header("Location: $errorURL \n\n");
		exit();
	}
?>
