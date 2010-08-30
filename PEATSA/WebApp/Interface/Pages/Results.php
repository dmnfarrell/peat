<?php
	/*
	 * This page displays the status and results of a submitted job.
	 * The job is identified by its jobId which is passed in the URL to this page.
	 * The status is tracked by reading data from the PEAT-SA mysql database.
	 * The database table Jobs contains a row detailing the status of each calculation requested.
	 * While the job is detected as still running the page automatically refreshes itself.
	 *
	 * On detecting any error the flag $error is set to one and the information variables 
	 * $errorDomain, $errorDescription, $detailedDescription & $recoverySuggestion are set also.
	 * At the end of this php script if $error is one an URL is created for the error display page,
 	 * using the information on the problem, and the brower redirected there.
	 */

	//Check the required data is present, and the database is accessible, before doing anything.
	//This has to be the very first thing in the file for it to work.
	//The required data is a jobID and an entry in the database for the job.
	include 'PageComponents.php';
	include_once 'PEATSA/Database/InterfaceKit/Utilities.php';

	if(!array_key_exists('jobId', $_GET))
	{
		$recoverySuggestion = 'Did you try to navigate directly to Results.php?<br>';
		$recoverySuggestion = $recoverySuggestion."If you haven't submitted a job then go to the main page first.<br><br>";
		$recoverySuggestion = $recoverySuggestion."If you were trying to access your results you need to know your job id.<br>";
		$recoverySuggestion = $recoverySuggestion."If you don't know this, then we're afraid you'll have to resubmit your job ...";

		//Redirect to the error page
		$errorURL = create_error_url($_GET["jobId"],
						'PDT.NavigationErrorDomain',
						'Required data is missing',
						'No job identification present in request to Results.php',
						$recoverySuggestion);
		header("Location: $errorURL \n\n");
		exit();
	}
	else
	{
		$error = 0;		
		//Get the state of each requested calculation
		$states = get_calculation_state_for_job( $_GET['jobId'], $error);
		if($error == 0)
		{	
			//If the job is still running set the page to refresh		
			$running = number_running_calculations($states);
			if($running != 0)
			{
				header('refresh: 10;');
			}
			else
			{
				//Send an email saying its finished
				//If the user did give an email this won't do anything
				send_job_completion_mail($_GET['jobId']);
			}
		}
		else
		{
			send_job_error_mail($_GET['jobId']);
			$errorURL = create_error_url_from_array($states);
			header("Location: $errorURL \n\n");
			exit();
		}
	}
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<!-- This is the results page for a PDT job 
	 It displays a table containing the files generated for the job so far -->
<html>

	<head>
	<?php
		if($running != 0)
		{
			echo "<title>PEAT-SA (Beta) - Processing Job</title>";
		}
		else
		{
			echo "<title>PEAT-SA (Beta) - Job Finished</title>";
		}
		
	?>	
		
		<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
		<link href="ResultsPageStyle.css" rel="stylesheet" type="text/css">
		<script type="text/javascript" src="../../jmol/Jmol.js"></script>
		<script type="text/javascript" src="../../js/mootools-1.2.1-core-nc.js"></script>
		<script type="text/javascript" src="../../js/mootools-1.2-more.js"></script>
		<script type="text/javascript" src="http://www.google.com/jsapi"></script>
		<script type="text/javascript" src="ResultsGraphs.js"></script>
	</head>
	<body onload=InitialiseGraph()>
	<div id="container">
	<?php
		add_navigation_pane();
		if($running == 0)
		{
			echo '<div class="jmol-controls">';
			$jobId = $_GET['jobId'];

			echo '<form align=left enctype="multipart/form-data" name="ColourPDB" title="ColourPDB">';
			echo '<h3>JMol Controls</h3>';
			echo '<p class="jmol-controls">';
			
			//Create the dynamic javascript
			//First create an URL to pass to JSGetCalculationData.php
			$temp = array('jobId' => $jobId);
			$array = array_merge($temp, $states);
			$query = http_build_query($array, '&amp');
			$server = $_SERVER['SERVER_NAME'];
			$port = $_SERVER['SERVER_PORT'];
			$javascriptURL =  "http://$server:$port/PEAT_SA/Pages/JSGetCalculationData.php?$query";
			
			//Next call JSCalculationData - This will dynamically create the requried JS functions
			//It also contains the color_calculation function
			echo "<script  type=\"text/javascript\" src='$javascriptURL'></script>";
			echo "<script  type=\"text/javascript\" src='ColourStructure.js'></script>";
			
			//Create the buttons
			echo "<input type='radio' name='color' value='structure' onclick=colorStructure(this.value)>";
			echo "<label class='jmol-controls'>Color Structure</label><br><br>";

			$i = 1;
			foreach($states as $calculation => $value)
			{
				if(($value == "Finished") && ($calculation != "Scan"))
				{
					echo "<input type='radio' name='color' value=$calculation onclick=colorStructure(this.value)>";
					echo "<label class='jmol-controls'>Color $calculation</label><br><br>";
					$i++;
				}
			}	
			
			echo '<b style="line-height:200%;">Rotate</b>:<br> click + move pointer<br>';
			echo '<b style="line-height:200%;">Zoom</b>:<br> shift + click + move pointer<br>';
			echo '<b style="line-height:200%;">Translate</b>:<br> shift + double click + move pointer<br>';
			echo '</p></form>';
			echo'</div>';
		}
		
	?>
	<div id="content" class="resultsPage">
	<?php
		/*
		 * All the information on the job was retrieved above. 
		 * Now create the html page to display it.
		 * The page consists of three parts
		 * A - A Jmol Applet (If Job is finished) Or a Spinny wheel (If its not)
		 * B - A Table showing the status of each job and providing links to results
		 * C - A final message
		 */
		
		//Display nice spinny wheel if the job is running
		$jobId = $_GET['jobId'];
		if($running != 0)
		{
			echo '<img src="Resources/loadinfo.net.gif" height="48" width="48" class="loading"><br>';
		}
		else
		{
			echo '<div id="jmol-display">';
				echo '<div id="Applet">';
				//Insert the jmol applet
				add_jmol_view($jobId);
				echo '</div>';
				//Initial inclusion of scale
				echo '<div id="Scale" class="Scale" style="display: none">';
				echo '<img class="Scale" src="Resources/Scale.png"><br>';
				echo '<div id="maxvalue" style="position: absolute; display: block;top:0px;left:42px">Value One</div>';
				echo '<div id="midvalue" style="position: absolute; display: block;top:142px;left:42px">Value One</div>';
				echo '<div id="minvalue" style="position: absolute; display: block;top:282px;left:42px">Value One</div>';
				echo '<div style="display: block; margin-top:15px;">kj.mol<sup>-1</sup></div>';
				echo '</div>';
			echo '</div>';
		}
		
		//Next create a table displaying the calculation state
		//1 - Check if any of the calculations are finished 
		//2 - Create links to the download page for any calculations that are
		//3 - Populate a table with the status of each job
		
		//Table layout
		echo '<div class="table">';
		echo '<br><table class="results" border=1><tbody><tr class="results"><th class="results">Task</th><th class="results">Status</th>';
		echo '<th class="results">Results</th></tr>';
		
		$value = get_mutant_creation_state($states);
		
		echo "<tr class='results'>";
		echo "<td class='results'>Mutant Creation</td>";
		echo "<td class='results' id='status'`>$value</td>";
		echo "<td class='results'`>None</td>";
		echo "</tr>";
		
		foreach($states as $calculation => $value)
		{
			if(($value != 'Finished') && ($value != 'NotSelected'))
			{
				//FIXME: Hack - Replace 'Scan' with deltapKa
				if($calculation == 'Scan')
				{
					$calculation = '&Delta;pKa';
				}
			
				echo "<tr class='results'>";
				echo "<td class='results'>$calculation</td>";
				echo "<td class='results' id='status'`>$value</td>";
				echo "<td class='results'`>None</td>";
				echo "</tr>";
			}
			else if($value == 'Finished')
			{
				//Construct download url for csv results file - just have to pass the data along	
				$array = array('jobId' => $jobId, 'calculation' => $calculation);
				$query = http_build_query($array, '&amp');
				$server = $_SERVER['SERVER_NAME'];
				$port = $_SERVER['SERVER_PORT'];
				
				$downloadURL = "http://$server:$port/PEAT_SA/Pages/DownloadMatrix.php?$query";
				$graphLink = "(<a title='Click for interactive results graph' href=javascript:ShowResultsGraph('$calculation')>Graph</a>)";

				//FIXME: Hack - Replace 'Scan' with deltapKa
				if($calculation == 'Scan')
				{
					$calculation = '&Delta;pKa';
				}
				
				$name = $calculation."Results.csv";
				$downloadLink = "<a href='$downloadURL' type='text/plain'>$name</a>";
		
				echo "<tr class='results'>";
				echo "<td class='results'>$calculation</td>";
				echo "<td class='results' id='status'>Finished</td>";
				echo "<td class='results'>$downloadLink $graphLink</td>";
				echo "</tr>";
			}
		}

		echo '</tbody>';
		echo '</table>';
		echo '</div>';
		
		if($running != 0)
		{
			//Get the status of the job
			//The following function returns a message based on the jobs current state
			//Job Queued - A string giving the position of the job in the queue
			//Job Running - 'Job is running'
			//Job not queued - 'Unknown'
			$error = 0;
			$message = get_queue_status_for_job($jobId, $error);
			if($message == "Unknown")
			{
				$message = "Your job is being processed";
			}
		
			echo '<div id="graphContainer">';
			
			echo "<p class='center'>$message.<br>";
			echo 'The table will be updated automatically with links to the results as the selected calculations finish.</p>';	
			echo "<p class='center'>";
			echo "Note: If you navigate away from this page you will need to save a link to it in order to retrieve your results.</p>";
		}
		else
		{
			echo '<div id="graphContainer">';
			echo '<p class="center">Your job is finished.<br>';
			echo "This page, and your results, will be accessible for two weeks.</p>";
		}
		
		echo '<div style="display: none;" id="resultsGraph" class="resultsGraph" name="someLink"></div>';

		echo "</div>";	
		echo "</div>";	

		add_colophon();
	echo "</div>";	
		
    ?>
    </body>
</html>
