<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>

	<head>
		<title>PEAT-SA (Beta) - Home</title>
		<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
		<link href="FrontPageStyle.css" rel="stylesheet" type="text/css">
		<script type="text/javascript" src="../../js/mootools-1.2.1-core-nc.js"></script>
		<script type="text/javascript" src="../../js/mootools-1.2-more.js"></script>
		<script type="text/javascript" src="../../js/parseuri.js"></script>
		<script type="text/javascript" src="ResultsGraphs.js"></script>
		<script>
			//Hide the mutation list box on load
			function HideBoxes()
			{
				//This wraps these elements in sliding boxes
				//and retracts those boxes
				HideElement("mutationListBox");
				HideElement("calculationHelpBox");
				HideElement("mutationHelpBox");
				//This sets the elements display to 'block'
				//Initially they are 'none' so they don't 
				//appear at all when the page loads
				ShowContent("mutationListBox");
				ShowContent("calculationHelpBox");
				ShowContent("mutationHelpBox"); 
			}
		</script>
		
		<script type="text/javascript">

		  //Analytics
		  var _gaq = _gaq || [];
		  _gaq.push(['_setAccount', 'UA-21087105-2']);
		  _gaq.push(['_trackPageview']);

		  (function() {
		    var ga = document.createElement('script'); 
		    ga.type = 'text/javascript'; 
		    ga.async = true;
		    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
		    
		    var s = document.getElementsByTagName('script')[0];
		    s.parentNode.insertBefore(ga, s);
		  })();

		</script>	

	</head>
	<body onload=HideBoxes()>
	<div id="container">
	<?php 
		include 'PageComponents.php';
		include_once 'PEATSA/Database/InterfaceKit/Session.php';
		//Start a session for this visitor to this PEATSA server (the domain has to be set on a per-installation basis)
		//This is required to enable two "web-applications" on the same server use
		//PEATSA/Database/InterfaceKit with different default databases.
		//Note - The argument is here is the path on the server the session cookie will be associated with
		//It assumed everything is installed under /PEATSA/ on the server
		start_peatsa_database_session('/PEATSA/');
		
		//Set the database configuration file to use with this server
		//For all visitors to this site the data retrieval functions will use this file
		//to determine the db to connect to.
		//Also must be set on a per installation basis
		set_peatsa_database_session_configuration_file('../../Resources/webApplication.conf');
		
		add_navigation_pane();
	?>	
	<div id="content">
		<h1><span>	
		PEAT-SA
		</span></h1>
		<h3>Discover the effect of mutations</h3>
		<p>
		The PEAT-SA webserver predicts the effect of mutations on a number of properties of a protein. <br>
	         It is developed by the <a href="http://enzyme.ucd.ie/">Nielsen Group</a> in University College Dublin.
		</p>
		<!-- enctype required for file uploads - not well documented anywhere! -->
		<form enctype="multipart/form-data" action="/PEATSA/cgi-bin/ProcessSubmission.py" method="post" name="PDT">
		
		<!-- Calculation Section -->
		<p style="margin-bottom: 1em">
                <b style="font-size: 14px; color: #FFFFFF">1. Choose Calculations</b>
		<a name='calculationHelpBox' id='calculationHelpBoxLink' onclick=ToggleElement("calculationHelpBox")>
		<img title="Show Calculation Help" src="Resources/HelpIcon.png" class="helpImage"/></a>
		</p>
		<div class="helpBox" id="calculationHelpBox">
			<b>&Delta;pKa</b> calculates how a mutation affects the pKa values of the specified ionisable residues in the protein.<br>
			<b>Stability</b> calculates the change in stability of the protein due to each mutation. <br>
			If <b>Binding</b> is selected, and a ligand is supplied, &Delta;G<sub>bind</sub> and &Delta;&Delta;G<sub>bind</sub>
			is calculated for each mutant.	<br>
		</div>
		<p style="margin-top: 1em">	
			<input type="checkbox" name="calculation" value="scan"> - &Delta;pKa		
			<label style="padding-left: 5em; width:10em">Ionisable Groups:</label> 
			<input type="text" value="All" size=15 id="ionisableGroups" name="ionisableGroups">  
			<br>
			<input type="checkbox" name="calculation" value="stability" checked="checked"> - Stability
			<br>
			<input type="checkbox" name="calculation" value="binding"> - Binding
		</p>
		
		<!-- Mutation Section -->
		<p style="margin-bottom: 1.0em"><b style="font-size: 14px; color: #FFFFFF;">2. Specify Mutations</b>
		<a name='mutationHelpBox' id='mutationHelpBoxLink' onclick=ToggleElement("mutationHelpBox")>
		<img title="Show Mutation Help" src="Resources/HelpIcon.png" class="helpImage"></a>
		</p>
		<div class="helpBox" id='mutationHelpBox'>
			There are two methods for specifying the mutations to be tested:<br><br>
			<b>Residue Scan</b> will calculate the effect of a single point mutation of each residue to the amino-acid<br>
			 chosen from the drop-down list.<br>
			<b>Mutation List</b> allows you to supply a file defining an arbitrary set of mutations to be tested.<br>
			The <a href="/PEATSA/Pages/MutationList.php">Mutation List Format</a> page explains the files format.
		</div>
		<p style="margin-bottom: 0.7em; margin-top: 1.0em">
			<input type="radio" name="mutation" value="mutationFile"><label>Mutation List</label>
			<input type="file" name="mutationFile" size=10>
			<br>
			<a name='mutationListBox'id='mutationListBoxLink' onclick=ToggleElement("mutationListBox")>  Or create a mutation list now ... [+]</a>
		</p>
		<div id="mutationListBox">
			<textarea name="mutationListArea" id='mutationListArea' cols="75" rows="7">Enter mutation list here. Please ensure format is correct.
The Mutation List Format page (link in navigation panel) explains the format in detail.
If a file is uploaded it will override a list entered here. (Note delete these lines before writing).		
			</textarea>
		</div>
		<p style="margin-top: 0em">
		<input type="radio" name="mutation" value="scan" checked="checked"><label>Residue Scan</label>
		<?php
			//Get the amino acids
			$data = amino_acids();
			
			//Sort the data directories and add them to a pop up list
			sort($data);
			echo "<select name='residue'>";
			foreach($data as $index => $value)
			{
				echo "<option value='$value'>$value</option>";
			}
			echo "</select>";
		?>
		<br><br><br>

		<!-- File Upload Section -->
		<b style="font-size: 14px; color: #FFFFFF">3. Upload Files</b><br><br>
		<strong><label>PDB file </label></strong>
		<input style="margin-bottom: 1.5em;" type="file" name="pdbFile" size="10">
                or provide a pdb/pKD ID (
		<?php
			$server = $_SERVER['SERVER_NAME'];
			$port = $_SERVER['SERVER_PORT'];
			echo "<a href='http://$server:$port/PEATSA/Pages/FAQ.php#pKaCode'>What's this?</a>";
		?>	
		)
		<input type="text" size=10 id="pdbCode" name="pdbCode">  
		<br>
		<strong><label>Ligand file (mol2 format) </label></strong>
		<input type="file" name="ligandFile" size="10">
		<br>
		<strong><label style='width:30em'>Email address (for notification of job completion):</label></strong>
		<input type="text" size=25 id="email" name="email" style="margin-top: 1.5em; margin-left: 0em">
		<br>
		<input type="submit" name="submitCalculation" value="Submit Job">	
		<br>
		</form>
	</div>
		<?php	
			add_colophon();
		?>
	</div>	
	</body>
</html>
