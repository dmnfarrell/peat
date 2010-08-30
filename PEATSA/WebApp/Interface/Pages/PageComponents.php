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

	//Contains functions that include various components on the web-app pages
	include_once 'UtilityFunctions.php';
	include_once 'PEATSA/Database/InterfaceKit/DataRetrieval.php';

	//Adds the navigation pane to a page
	function add_navigation_pane()
	{
		$server = $_SERVER['SERVER_NAME'];
		$port = $_SERVER['SERVER_PORT'];
		echo '<div id="navigation">';
		echo '<p class="navigation">';
		echo "<a class='navigation' href='http://$server:$port/PEATSA/Pages/FrontPage.php'>Home</a><br><br>";
		echo "<a class='navigation' href='http://enzyme.ucd.ie'>Group Site</a><br><br>";
		echo "<a class='navigation' href='http://$server:$port/PEATSA/Pages/MutationList.php'>Mutation Lists</a><br><br>";
		echo "<a class='navigation' href='http://$server:$port/PEATSA/Pages/Methods.php'>Methods</a><br><br>";
		echo "<a class='navigation' href='http://$server:$port/PEATSA/Pages/Technology.php'>Technology</a><br><br>";
		echo "<a class='navigation' href='http://$server:$port/PEATSA/Pages/FAQ.php'>FAQ</a><br>";
		echo '</p>';
		echo '<div id=sfiLogoDiv>';
		echo '<img id=sfiLogo alt="SFI Logo" src="Resources/sfi_cmyk.png">';
		echo '</div>';
		echo '</div>';	
	}

	//Adds the colophon to a page
	function add_colophon()
	{
		echo '<div id="colophon">';
		echo '<p class="colophon">';
		echo 'Copyright 2009-2010 @ Michael Johnston & Jens Erik Nielsen';	
		echo '</p>';
		echo '</div>';
	}
	
	//Adds a JMol applet showing the structure used in job with identification $jobID to a page
	function add_jmol_view($jobID)
	{
		/*
		 * This php function writes out the javascript code necessary for 
 		 * displaying the structure used in the job in JMol
		 *
  		 * However we can't use the $structure string directly as the newlines are not escaped.
		 * This causes a problem when we try to do
		 *	echo "var string = \"".$structure."\";\n";
		 * as php will interpret the "\n" characters when it executes the echo,
		 * leading to a syntax error in the Javascript as $structure will appear as a multiline string.
		 * Instead we have to first remove all "\n" characters.
		 * We then add the newlines back using single quoated string so php wont interpret them
		 */	
		$error = 0;
		$structure = get_structure_for_job($jobID, $error);
		$structure = escape_newlines($structure);
		$structureCode = get_pdb_code_for_job($jobID, $error);
		echo '<p class="center" style="font-size: 16px; margin-left:auto; margin-right:auto">';
		echo "$structureCode";
		echo "</p>";
		$script = "\"set perspectiveDepth on;ribbons;cpk off;wireframe off;\"";
		echo '<script language="JavaScript" type="text/javascript">';
		echo "var string = \"".$structure."\";\n";
		echo "jmolInitialize('../../jmol');\n";
		echo "jmolAppletInline([400, 400], string, $script);\n"; 
		
		//If there is a ligand insert it
		//Note: This doesn't work because jmolLoadInlineScript doesn't
		//work until the applet is initialised
		/*$error = 0;
		$states = get_calculation_state_for_job($jobID, $error);
		if(($error == 0) && ($states['Binding'] != 'NotSelected'))
		{
			$ligand = get_ligand_for_job($jobID);
			$ligand = escape_newlines($ligand);
			$script = "\"select 2.0; wireframe on;\"";
			echo "var string = \"".$ligand."\";\n";
			echo "jmolLoadInline(string, '0');\n"; 
			//echo "jmolLoadInlineScript(string, $script, '0');\n"; 
			echo '</script>';
		}*/				

		echo '</script>';  	
	}	
?>
