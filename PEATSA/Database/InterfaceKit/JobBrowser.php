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

	include_once 'Structures/DataGrid.php';
	include_once 'HTML/Table.php';
	include_once 'Formatting.php';

	function format_job_browser_calculation_column($params, $args = array())
	{
		extract($params);
		extract($args);
		$state = $record[$calculation];
		if($state == 'Finished')
		{
			$server = $_SERVER['SERVER_NAME'];
			$port = $_SERVER['SERVER_PORT'];
			$array = array("jobId"=>$record['jobid'], "calculation"=>$calculation);
			$query = http_build_query($array, '&amp');
			//FIXME: Need to dynamically set URL
			$link = "http://$server:$port/Test/Pages/DisplayTable.php?$query";
			return "<a href=$link>View Data</a>";
		}
		else
		{
			return "N/A";
		}	
		
	}
	
	function format_job_browser_structure_column($params, $args = array())
	{
		extract($params);
		$server = $_SERVER['SERVER_NAME'];
		$port = $_SERVER['SERVER_PORT'];
		$array = array("jobId"=>$record['jobid']);
		$query = http_build_query($array, '&amp');
		//FIXME: Need to dynamically set URL
		$link = "http://$server:$port/Test/Pages/DisplayStructure.php?$query";
		$code = $record['pdbid'];
		return "<a href=$link>$code</a>";
	}
	
	function format_job_browser_columns($browser)
	{
		//Format the mutationCode column so you can click on the code to generate the mutant
		$calculation = "stability";
		$column = $browser->getColumnByName($calculation);
		$column->setFormatter('format_job_browser_calculation_column', array("calculation"=>$calculation));
		
		$calculation = "binding";
		$column = $browser->getColumnByName($calculation);
		$column->setFormatter('format_job_browser_calculation_column', array("calculation"=>$calculation));
		
		$calculation = "scan";
		$column = $browser->getColumnByName($calculation);
		$column->setFormatter('format_job_browser_calculation_column', array("calculation"=>$calculation));
		
		$column = $browser->getColumnByName('pdbid');
		$column->setFormatter('format_job_browser_structure_column');
	}

	function job_browser($table, $params, $limit = 12, $statement = 'None')
	{
		//Test implementation - 
		// Setup your database connection
		$dsn = 'mysql://'.$params['username'].':'.$params['password'].'@'.$params['hostspec'].'/'.$params['database'];
		$options = array('dsn' => $dsn);

		// Bind a basic SQL statement as datasource
		$browser = new Structures_DataGrid($limit);
		if($statement == 'None')
		{
			$statement = "SELECT JobID, Date, PDBID, Stability, Binding, Scan FROM $table";
		}
		$error = $browser->bind($statement, $options);	
		//FIXME: Check error
		
		format_job_browser_columns($browser);
		set_datagrid_attributes($browser->getRenderer());

		return $browser;
	}

?>

