<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>

	<head>
		<title>PEAT-SA (Beta) s- Mutation Lists</title>
		<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
	</head>
	<body>
	<div id="container">
	<?php 
	include 'PageComponents.php';
	add_navigation_pane();
	?>	
	<div id="content">
		<h1><span>	
		PEAT-SA (Beta)
		</span></h1>
		<h3>Mutation Lists</h3>
		<p>
		PEAT-SA creates mutants of a supplied protein and calculates the effect of the mutations on the proteins properties.
		A <i>mutation list</i> allows you to specify the exact mutants you wish to examine.
		You can supply a mutation list to the program via the text-box on the submission page or by uploading a file.
		<br><br>
		There are two different formats for writing mutation lists, <em>Single Point Mutation</em> and <em>Standard</em>, and the following sections
		explain each.
		Which you use depends on your requirements but you <emph>cannot</emph> mix the two formats.
		</p>
		<p>
		<b style="text-decoration: underline;font-size: 14px;">Single Point Mutants</b><br><br>	
		Specifying a set of sinlge point mutants (mutants that contain only one mutated residue) is the simplest case.
		Each line in the mutation list has the following format:
		<br><br>
		<i>Chain:ResidueNumber, Mutation, Mutation ...</i>
		<br><br>
		For example,
		<br><br>
		<i>A:0001, ASP, HIS, ILE</i>
		<br><br>
		This concise format tells PEAT-SA to create three mutants of the uploaded protein by mutating the first residue of chain A.
		In the first mutant the residue will be changed to Aspartate, in the second to Histidine and so on.
		All 19 possible single point mutations can be specified for a given residue. 
		<b>Note: </b>The mutations must be specified using three-letter amino-acid codes. The codes do not
		have to be in upper case.
		<br>
		</p>
		<p>
		<b style="text-decoration: underline;font-size: 14px;">Standard Format</b><br><br>	
		The standard format is more verbose than the SPM format but allows specification of mulit-mutation mutants.
		Each line in the mutation list defines one mutant and has the following format:
		<br><br>
		<i>Chain:ResidueNumber:Mutation+Chain:ResidueNumber:Mutation+...</i>
		<br><br>
		For example,
		<br><br>
		<i>A:0001:ASP+B:0001:ASP<br>A:0012:VAL+A:0124:ILE</i>
		<br><br>
		The above lines tell PEAT-SA to create two mutants. The first has the first residue in both chains mutanted to ASP.
		The second has two mutations in chain A - residue 12 is mutated to VAL and residue 124 is mutated to ILE. <br>	
		You can also use the more condensed format,
		<br><br>
		<i>A1D+B1D<br>A12V+A124I</i>
		<br><br>
		<b>Note:</b> Although any number of mutations may be specified per mutant, the more there are the less reliable the results will be.
		</p>
	</div>
		<?php	
			add_colophon();
		?>
	</div>	
	</body>
</html>
