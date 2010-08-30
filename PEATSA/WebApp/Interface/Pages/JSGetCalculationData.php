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

	//PHP 
	//Contains two functions that dynamically creates javascript functions
	// that returns the data for each calculation as an array
	//If the calculation is not run the created function returns nil
	Header("content-type: application/x-javascript");
	include_once "PEAT_SA/Database/InterfaceKit/DataRetrieval.php";
	include_once "PEAT_SA/Database/InterfaceKit/Utilities.php";
	
	//A php function that creates a javascript function! - for finished jobs
	//1. Get the data
	//2. Get the totals
	//3. Create a string with the totals.
	//4. Write javascript to assigne the string to a var
	function output_complete_javascript_function($jobId, $calculation)
	{
		echo "\nfunction create_"."$calculation"."_data_array(){\n";
		
		$error = 0;
		$retval = get_data_for_calculation($jobId, $calculation, $error);
		if($error == 0)
		{
			$data = matrix_from_csv_representation($retval['content']);
			
			if($calculation != "Scan")
			{
				$totalColumn = count($data[0]) - 1; 
				echo "var data = new Array();\n";
				foreach($data as $index => $value)
				{
					echo "var row = new Array();\n";
					echo "row[0] = '$value[0]';\n";
					echo "row[1] = '$value[$totalColumn]';\n";
					echo "data[$index] = row;\n";
				}
				echo "data.splice(0,1);";
			}
			else
			{
				echo "var data = new Array();\n";
				foreach($data as $index => $value)
				{
					echo "var row = new Array();\n";
					foreach($value as $counter => $result)
					{
						echo "row[$counter] = '$result';\n";
					}
					echo "data[$index] = row;\n";
				}
			}
			
			echo "return data;\n";
		}
		else
		{
			//Construct javascript redirect
			$errorURL = $retval;
			echo "document.location.href = '$errorURL';\n";
		}
		
		echo "}\n\n";		
	}
	
	//Returns an empty array
	//Used for calculations that were not requested but we need a wrapper for
	//to avoid undefined errors
	function output_empty_javascript_function($calculation)
	{
		echo "\nfunction create_"."$calculation"."_data_array(){\n";
		
		echo "var data = new Array();\n";
		echo "return data;\n";
		
		echo "}\n\n";		
	}
	
	//The Stability function
	if($_GET["Stability"] == "Finished")
	{
		output_complete_javascript_function($_GET["jobId"], "Stability");
	}	
	else
	{
		output_empty_javascript_function("Stability");
	}
	
	//The Scan function
	if($_GET["Scan"] == "Finished")
	{
		output_complete_javascript_function($_GET["jobId"], "Scan");
	}	
	else
	{
		output_empty_javascript_function("Scan");
	}
	
	//The Scan function
	if($_GET["Binding"] == "Finished")
	{
		output_complete_javascript_function($_GET["jobId"], "Binding");
	}	
	else
	{
		output_empty_javascript_function("Binding");
	}
	
?>
