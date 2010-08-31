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

	//Contains various utility functions used by other php code
	//Including processing data e.g. csv files, retrieved from the database etc.
	include_once('Config.php');
	include_once 'Session.php';
		
	function matrix_from_csv_representation($string)
	{
		//Steps 
		//Write string to a temp file
		//Read back in each row using fgetcsv
		
		$handle = tmpfile();
		fwrite($handle, $string);
		fseek($handle, 0);
		
		$matrix = array();
		$count = 0;
		while(($row = fgetcsv($handle)) !== FALSE)
		{
			$matrix[$count] = $row;
			$count++;
		}
		
		fclose($handle);
		return $matrix;
	}
	
	//Escapes a string so php does not interpret newlines
	function escape_newlines($string)
	{
		//Remove double quoted newlines and replace with single-quoted newlines
		$array = explode("\n", $string);
		$string = implode('\n', $array);
		$string = $string.'\n';
		
		return $string;
	}
	
	//Extracts the first line from a csv file and returns an array of the elements
	//along with a string containing the body of the csv data
	function extract_csv_headers($string)
	{
		$array = explode("\n", $string);
		$headers = explode(",", $array[0]);
		foreach($headers as $index => $value)
		{
			$value = trim($value);
			$headers[$index] = $value;
		}
		
		unset($array[0]);
		$string = implode("\n", $array);
		return array("headers"=>$headers, "body"=>$string);
	}
	
	function number_running_calculations($data)
	{
		$running = 0;
		foreach($data as $row => $value)
		{
			if(($value != "Finished") && ($value != "NotSelected"))
			{
				$running = $running + 1;
			}
		}
		
		return $running;
	}
	
	function get_mutant_creation_state($states)
	{
		//If no calculations are running then mutant creation is finished
		//Note: A calculation is only considered not running when it is finished
		if(number_running_calculations($states) == 0)
		{
			$retval = 'Finished';
		}
		else
		{
			$retval = 'Running';
			foreach($states as $calculation => $value)
			{
				//If any calculation is Queued the job isn't
				//running and the mutant creation step is also queued
				if($value == 'Queued')
				{
					$retval = 'Queued';
					break;
				}
				else if($value == 'Running')
				{
					//If any calculation is running
					//then mutant creation is finished
					$retval = 'Finished';
					break;
				}
			}
		}
		
		return $retval;
	}
?>
