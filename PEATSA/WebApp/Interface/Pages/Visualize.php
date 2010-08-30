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

	$error = 0;
	$jobId = $_GET['jobId'];
	$states = get_calculation_state_for_job($jobId, $error);
	if($error == 1)
	{
		$errorURL = create_error_url_from_array($states);
		header("Location: $errorURL \n\n");
		exit();
	}
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<!-- This is the results page for a PDT job 
	 It displays a table containing the files generated for the job so far -->
<html>

	<head>
	<?php
		echo "<title>PEAT-SA (Beta) - Visualization Test</title>";
	?>	
		<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
		<script type="text/javascript" src="ResultsGraphs.js"></script>
		<script type="text/javascript" src="http://www.google.com/jsapi"></script>
		<script type="text/javascript">
    
		// Load the Visualization API and the piechart package.
		google.load('visualization', '1', {'packages':['piechart','areachart', 'columnchart', 'imagebarchart']});

		// Set a callback to run when the Google Visualization API is loaded.
		google.setOnLoadCallback(drawChart);

		var data = 0;
		
		function drawChart() 
		{
			// Create our data table.
			var calculationData = dataForCalculation('Stability');
			data = new google.visualization.DataTable();

			data.addColumn('string', 'Mutation');
			data.addColumn('number', 'Total');
			for(var i=1; i < calculationData.length; i++)
			{
				calculationData[i][1] = parseFloat(calculationData[i][1]);
				data.addRow(calculationData[i]);
			}

			// Instantiate and draw our chart, passing in some options.
			var chart = new google.visualization.ImageBarChart(document.getElementById('chart_div'));
			chart.draw(data, {width: 400, height: 240, is3D: true, title: 'Results'});
		}

		function setGraph(graph)
		{
			console.log(data);
			console.log(graph);
			if(graph == "imageBarChart")
			{
				var chart = new google.visualization.ImageBarChart(document.getElementById('chart_div'));
			}
			else if(graph == "piechart")
			{
				var chart = new google.visualization.PieChart(document.getElementById('chart_div'));
			}
			else if(graph == "areaChart")
			{
				 var chart = new google.visualization.AreaChart(document.getElementById('chart_div'));
			}
			else if(graph == "columnChart")
			{
				 var chart = new google.visualization.ColumnChart(document.getElementById('chart_div'));
			}
			chart.draw(data, {width: 740, height: 500, legend: 'none', fontSize: 8, title: 'Stability'});

		}
		</script>

	</head>
	<body>
	<div id="container">
	<?php
		add_navigation_pane();

		//Create the dynamic javascript by passing an URL to JSGetCalculationData.php
		//First create the URL
		$temp = array('jobId' => $jobId);
		$array = array_merge($temp, $states);
		$query = http_build_query($array, '&amp');
		$server = $_SERVER['SERVER_NAME'];
		$port = $_SERVER['SERVER_PORT'];
		$javascriptURL =  "http://$server:$port/PEAT_SA/Pages/JSGetCalculationData.php?$query";
			
		//Next call JSCalculationData - 
		echo "<script  type=\"text/javascript\" src='$javascriptURL'></script>";
		echo "<script  type=\"text/javascript\" src='ColourStructure.js'></script>";

		echo '<form align=left enctype="multipart/form-data" name="Grapher" title="Grapher">';
		echo '<h3>Visualise Controls</h3>';
		echo '<p class="jmol-controls">';

		echo "<select name='residue' onclick=setGraph(document.Grapher.residue.value)>";
		echo "<option value='piechart'>Pie Chart</option>";
		echo "<option value='imageBarChart'>Bar Chart</option>";
		echo "<option value='areaChart'>Area Chart</option>";
		echo "<option value='columnChart'>Column Chart</option>";
		echo "</select>";

		echo "</p>";
		echo "</form>";
	?>
	<div id="content">
	<div id="chart_div"></div>
	</div>
	</div>
</body>
</html>		
