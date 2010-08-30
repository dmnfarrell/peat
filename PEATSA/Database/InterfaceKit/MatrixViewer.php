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

	include_once 'DataRetrieval.php';
	include_once 'Formatting.php';
	include_once 'Structures/DataGrid.php';
	include_once 'Structures/DataGrid/DataSource/CSV.php';
	include_once 'HTML/Table.php';

	function addMutationGenerationLink($params, $args = array())
	{
		extract($params);
		extract($args);
		$server = $_SERVER['SERVER_NAME'];
		$port = $_SERVER['SERVER_PORT'];
		$array = array("jobId"=>$jobId, "mutationCode"=>$record[0]);
		$query = http_build_query($array, '&amp');
		$link = "http://$server:$port/PEAT_SA/cgi-bin/GetMutant.py?$query";
		$code = $record[0];
		return "<a href=$link>$code</a>";
	}
	
	//Returns a Structures_DataGrid object initialised to display $matrix
	//$matrix is a csv file in a string
	function matrix_viewer($matrix)
	{
		//Parse the csv string into headers and content
		$array = extract_csv_headers($matrix);
		$headers = $array["headers"];
		$matrix = $array["body"];			
		
		//Set up the dataSource and dataGrid
		$dataSource = new Structures_DataGrid_DataSource_CSV();
		$dataSource->bind($array["body"]);	
		
		$matrixViewer = new Structures_DataGrid(50);
		$matrixViewer->bindDataSource($dataSource);
		
		//Set the column headers
		$columns = $matrixViewer->getColumns();
		foreach($columns as $index => $column)
		{
			$column->setLabel($headers[$index]);
		}
		
		//Format the mutationCode column so you can click on the code to generate the mutant
		$mutationColumn = $columns[0];
		$mutationColumn->setField('mutationCode');
		$mutationColumn->setFormatter('addMutationGenerationLink', array('jobId'=>$jobId));
		
		set_datagrid_attributes($matrixViewer->getRenderer());	
		
		return $matrixViewer;
	}
?>


