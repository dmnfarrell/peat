<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>

	<head>
		<title>PEAT-SA (Beta) - FAQ</title>
		<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
		<link href="FAQStyle.css" rel="stylesheet" type="text/css">
		<script type="text/javascript" src="../../js/mootools-1.2.1-core-nc.js"></script>
		<script type="text/javascript" src="../../js/mootools-1.2-more.js"></script>
		<script type="text/javascript" src="../../js/parseuri.js"></script>
		<script type="text/javascript" src="ResultsGraphs.js"></script>
		<script>
			//Hide all the faq answers on load
			function HideFAQDivs()
			{
				var divs = new Array("acronym", "calculationTime", "state", "pKaCode", "calculationDifference", 
						"calculationFormat", "weird", "chainbreak-error", "mainchain-error", "binding", "coloring",
						"mutationList", "graph", "mutationColumn",
						"resultsColumns", "colorScale", "jobTime", "units", "pka-error");
			
				//If the div's name is in the anchor don't hide it
				//This allows urls linking to the question to arrive here with the question open
				data = parseUri(location.href);
			
				for(var i=0; i<divs.length; i++)
				{
					if(data.anchor != divs[i])
					{
						HideElement(divs[i]);
					}
					
					document.getElementById(divs[i]).style.display = "block";
				}
			}
		</script>	
	</head>
	<body onload=HideFAQDivs()>
	<div id="container">
	<?php 
	include 'PageComponents.php';
	add_navigation_pane();
	?>	
	<div id="content">
		<h1><span>	
		PEAT-SA (Beta)
		</span></h1>
		<h3 style="font-size: 16px;">Frequently Asked Questions</h3>
		<p>Click on a question for an answer. <br></p>
		<b style="font-size: 14px; margin-left:3em;">General</b><br><br>
		
		<a name='acronym' class='faq-question' onclick=ToggleElement("acronym")>What does PEAT-SA stand for?</a><br>
		<div id="acronym" class='faq-answer'>
		PEAT-SA stands for Protein Engineering and Analysis Tool - Structure Analysis. It is part of the PEAT suite
		of products developed by the Nielsen group which also includes <a href="http://enzyme.ucd.ie/main/index.php/PEAT_DB">PEAT-DB<a>.
		</div>
		
		
		<a name='state' class='faq-question' onclick=ToggleElement("state")>What is PEAT-SA's development state?</a><br>
		<div id="state" class='faq-answer'>
		Currently we're in a 'beta' state. This means that the web-server works in all the cases we've tested
		but we expect there to be some bugs and rough-edges feature wise. Furthermore, the site may occasionally be down, 
		or act strangely as it is under active development.<br>
		If you encounter any problems or have any suggestions on extra features/improvements you would 
		like to see, don't hesitate to contact us.
		</div>
		
		<a name='calculationTime' class='faq-question' onclick=ToggleElement("calculationTime")>How long does a calculation take?</a><br>
		<div id="calculationTime" class='faq-answer'>
		The calculation time depends on the protein but rough estimates would be (per mutant submitted)
		<ul>
		<li> Mutant Generation - 0.7 secs (200 residues) -> 1.5 secs (500 residues) </li>
		<li> Stability & Binding -  0.7 secs (200 residues) -> 1.9 secs (500 residues) </li>
		<li> pKa - One ionisable residue - 22.3 secs (64 residues) -> 45 secs (200 residues) </li>
		<li> pKa - Two ionisable residues - 25.2 secs (64 residues) -> 64 secs (200 residues) </li>
		</ul>
		So calculating the effect of mutating each residue to alanine in a 200 residue protein on stability and binding
		would take ~ 200*(0.7+0.7+0.7) = 7 minutes. &Delta;pKa calculations are more expensive. 
		However the time per ionsiable group decreases as you specify more groups.
		If you specify 1 ionsiable group the &Delta;pKa calculation would take 200*45 = 2.5 hours. 
		Specfying two ionisable groups would take 3.6 hours.
		</div>
		
		<a name='jobTime' class='faq-question' onclick=ToggleElement("jobTime")>Is there a time limit on jobs?</a><br>
		<div id="jobTime" class='faq-answer'>
		Yes - the time limit for a job is 12 hours. If your calculation exceeds this time your results will be lost.
		Therefore it is better to divide your top into smaller tasks if you think it may take around 12 hours.
		</div>
		
		<a name='pKaCode' class='faq-question' onclick=ToggleElement("pKaCode")>What is a pKD ID</a><br>
		<div id="pKaCode" class='faq-answer'>
		When you run a pKa calculation for a protein using <a href="http://enzyme.ucd.ie/cgi-bin/pKD/server_start.cgi">pKD Server</a> 
		the server assigns the job a unique identification. We call this a pKD ID. To perform a &Delta;pKa calculation
		with PEAT-SA you must supply the pKD ID of a pKD server job instead of a PDB ID or file.
		</div>
		
		<a name='mutationList' class='faq-question' onclick=ToggleElement("mutationList")>What is a Mutation List?</a><br>
		<div id="mutationList" class='faq-answer'>
		A mutation list is a list of mutation codes, each of which defines a mutant of a protein structure.
		The <a href="/ProteinDesignTool/Pages/MutationList.php">Mutation List Format</a> page explains the format in detail.
		</div>
		
		<a name='calculationDifference' class='faq-question' onclick=ToggleElement("calculationDifference")>
		Is there a difference between doing calculations together to doing them separately?</a><br>
		<div id="calculationDifference" class='faq-answer'>
		Yes. For a binding calculation the mutants are modelled taking into account the ligand.
		If a stability and/or a &Delta;pKa calculation are also requested they use these mutants.
		</div>
		
		<a name='calculationFormat' class='faq-question' onclick=ToggleElement("calculationFormat")>
		What format to I get the calculation results in?</a><br>
		<div id="calculationFormat" class='faq-answer'>
		The results for each calculation are returned in csv (comma separated value) format, which
		can be directly loaded by any spreadsheet application.
		</div>
		
		<a name='weird' class='faq-question' onclick=ToggleElement("weird")>
		Why do the pages look strange?</a><br>
		<div id="weird" class='faq-answer'>
		If the website layout looks strange its is most likely due to the browser you are using.
		Older browsers do not correcly support CSS (Cascading Style Sheets) which are used to define the site layout.
		The website has been tested, and renders correctly, on Firefox 3+, Internet Explorer 7 & 8, Safari 3+ and Google Chrome.
		</div>

		<br><b style="font-size: 14px; margin-left:3em;">Errors</b><br><br>

		<a name='mainchain-error' class='faq-question' onclick=ToggleElement("mainchain-error")>
		I get the error "The supplied structure is missing main chain heavy atoms" but they all seem to be there?</a><br>
		<div id="mainchain-error" class='faq-answer'>
		The webserver expects every residue to have a main-chain oxygen called 'O'. 
		However some programs produce PDB files where the oxygens in the C-Terminal are called 'OT1' and 'OT2' (or some variation).
		This causes PEATSA to think an atom is missing. The easiest way around this problem is to use the standard PDB naming for
		these oxygens 'O' and 'OXT'.
		</div>

		<a name='chainbreak-error' class='faq-question' onclick=ToggleElement("chainbreak-error")>
		I get the error "The supplied structure contains at least one chain break" but there is none?</a><br>
		<div id="chainbreak-error" class='faq-answer'>
		In the situation where there is no actual chain break this often means that two chains have been given the same chain id in the PDB file. 
		PEATSA checks for chain-breaks by seeing if more than one C-Terminal is present for a particular chain id after it has cleaned the structure,
		and this will happen in the above situation.
		</div>
		
		<a name='pka-error' class='faq-question' onclick=ToggleElement("pka-error")>
		I get the error "No pKa data corresponding to supplied code [some-code] is present" when trying to run a &Delta;pKa calculation?</a><br>
		<div id="pka-error" class='faq-answer'>
		In order to run a &Delta;pKa you must previously have prepared your protein using <a href="http://enzyme.ucd.ie/cgi-bin/pKD/server_start.cgi">pKD Server</a>.
		This essentially precalculates information that can be used to accelerate any subsequent perturbation calculation
		on that structure.
		When you submit your protein to pKDServer your job gets assigned a pKD ID (see question above). It is this code you
		must specify when submitting a &Delta;pKa calculation, in order that the correct information can be found.
		</div>

		<br><b style="font-size: 14px; margin-left:3em;">Results</b><br><br>

		<a name='units' class='faq-question' onclick=ToggleElement("units")>
		What units are the results given in?</a><br>
		<div id="units" class='faq-answer'>
		The results for binding and stability calculations are in kj/mol while the &Delta;pKa results are in pKa units.
		1 pKa unit is equivalent to a free-energy difference of 5.7 kj/mol (at 300K).
		</div>
		
		<a name='binding' class='faq-question' onclick=ToggleElement("binding")>
		I get strange binding results, why?</a><br>
		<div id="binding" class='faq-answer'>
		If the total &Delta;&Delta;G<sub>bind</sub>'s seem physically strange, you should check that:
		<ul>
		<li> That the position of the ligand defined in the mol2 file is correct i.e. is in the binding site</li>
		<li> The charges of the ligand atoms are specified. </li>
		</ul>
		We find that when using some molecular visualisation programs to extract a ligand as a mol2 file the coordinates
		are not saved relative to the original protein position by default. 
		This leads the ligands position to be outside of the protein. 
		</div>

		<a name='coloring' class='faq-question' onclick=ToggleElement("coloring")>
		Coloring the structure doesn't work/is very slow?</a><br>
		<div id="coloring" class='faq-answer'>
		We have observed this problem with Firefox 3 on the Mac and it may affect other browsers.
		Due to the speed with which the browser executes javascript it can take up to 15secs for the 
		residues in the jmol structure to become coloured.
		</div>
		
		<a name='graph' class='faq-question'  onclick=ToggleElement("graph")>
		Why can't I colour/see a graph of the &Delta;pKa results?</a><br>
		<div id="graph" class='faq-answer'>
		Answer Coming Soon.
		</div>
		
		<a name='mutationColumn' class='faq-question'  onclick=ToggleElement("mutationColumn")>
		How do I interpret the string in the <em>Mutations column</em> of a results file?</a><br>
		<div id="mutationColumn" class='faq-answer'>
		The string in the mutations column indicates the mutant the results correspond to.
		It consists of a code for each mutation in the mutant concatenated by '+'.
		The code for each mutant has the form:<br>
		(ChainID)(ResidueIndex)(OriginalResidue)(MutatedResidue) e.g. A54AV
		</div>
				
		<a name='resultsColumns' class='faq-question'onclick=ToggleElement("resultsColumns")>
		Why aren't there separate columns for chain/residue number etc.?</a><br>
		
		<div id="resultsColumns" class='faq-answer'>
		Although such columns make sense if a mutant only contains a single mutation, they
		can't be used for multiple-mutation mutants. Therefore for consistency we use one format that
		fits all cases.
		</div>
		
		<a name='colorScale' class='faq-question' onclick=ToggleElement("colorScale")>
		Nothing is coloured when I click Colour Structure/Binding?</a> <br>
		<div id="colorScale" class='faq-answer'>
		Unfortunately the colouring only works for single point mutations, since there is no easy way to represent
		the values of double/triple mutants by colouring residues. A similar problem arises 
		if the same residue mutated twice in one job - the value corresponding to only one of these mutants can
		be shown.
		</div>
		<br><br>
	</div>
		<?php	
			add_colophon();
		?>
	</div>	
	</body>
</html>
