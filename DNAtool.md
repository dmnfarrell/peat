<img src='http://peat.googlecode.com/svn/wiki/images/DNAtool_logo.png' align='right' width='100'>
<b>DNAtool</b> is a stand-alone subcomponent of PEATDB that contains standard functions for DNA sequence analysis, and is designed for the researcher who works with smallish pieces of DNA containing a single gene.<br>
<br>
<h2>Introduction</h2>

This is a tool designed for biologists to speed up the process of primer design, restriction digest analysis, and sequencing results, but also allow for easy interpretation of cumbersome nucleotide sequences.<br>
<br>
DNA tool is accessed from the main PEAT window, using the toolbar at the top.  Once selected, two windows will automatically open, a DNA sequence manipulation window, and a primer database window. Provided you have uploaded a nucleotide sequence for your protein, you should see both the nucleotide sequence in green, and the corresponding amino acid sequence in black, below.<br>
<br>
<h2>Viewing a Nucleotide Sequence</h2>

There are a number of small changes that can be made to the nucleotide sequence, to make it easier to analyse.  The individual bases can be coloured according to nucleotide type, A (Blue), T (Red), C (Orange), G (Black).  Otherwise, nucleotides can be coloured as triplets, or viewed as a single colour. 	The spacing of the sequence can be altered to allow for easier reading, and the text size, font type, and font style, can all be altered.<br>
By default, the nucleotide sequence is numbered at every tenth residue, and highlighted with a red marker making it easy to determine you location in the sequence.  The corresponding amino acid sequence is also numbered, at every fifth residue. The entire length of the sequence can be viewed by using the scroll bar at the bottom of the window.<br>
<br>
Restriction sites in the nucleotide sequence are highlighted above the sequence with the name of the corresponding restriction enzyme.  The exact restriction digest site can be seen by left clicking on any one of the yellow restriction enzyme names.  In addition, selecting “digest details” in the upper right hand corner of the window provides a list of all restriction enzymes that cut, and the nucleotide positions of those cuts.<br>
<br>
<h2>Primer Design</h2>

To design a new mutagenic primer for your sequence, select “primer design” form the toolbar, and “mutagenic primer” from the drop down menu.  A new window called “mutagenic primer design”, will now open.  Inputting the number of the amino acid you want to mutate will bring you to that location in the sequence, and highlight the amino acid name in red.  Select the amino acid you wish to mutate too from the drop-down menu in the labelled, “New AA” and the desired melting temperature for this primer.  The default setting is 65oC.  Once you are happy with your selections click the “design primer” button and DNA Tool will generate a list of suitable primers.  All primers generated will be output in the window below.  Left clicking on a primer will align the nucleotide sequence of the primer with the original sequence.  Mis-matches in the sequence are highlighted in red, in addition to the codon change required to generate the mutation, the addition or removal of a restriction site was performed automatically.  The details of the primer, including the melting temperature, and the restriction site changes are displayed in a window beside the list of generated primers.<br>
<br>
<h2>DNA Sequencing Results/Analysis</h2>

DNATool allows for fast analysis of sequenced DNA results.  A newly sequenced forward and reverse DNA strand can be overlaid on the wild-type or mutant sequence in addition to the primer used to create that mutation.<br>
Select the “sequence analysis” tab from the toolbar, and “load multiple DNA sequences”.  Browse to the folder containing your sequenced DNA files, and select the files of interest.  These files will now appear in the sequences window.  Select both the forward and reverse sequences by pressing the Ctrl key.  Both sequences will now upload into the “DNA Sequence Manipulation” window.  To upload the mutagenic primer that generated the mutation, used the primer database window.<br>
<br>
Again, mis-matches are automatically highlighted in red, and changes in the amino acid sequence are shown beneath the nucleotide sequence.  The example shows a Aspartic acid to Alanine mutation at position 18 in the amino acid sequence.  The primer is coloured blue, the forward and reverse strand of sequenced DNA in purple, the WT nucleotide sequence is green.<br>
<br>
<h2>Screenshots</h2>

<img src='http://peat.googlecode.com/svn/wiki/images/DNAtool_mainwindow.png' width='600'>

Main window.<br>
<br>
<img src='http://peat.googlecode.com/svn/wiki/images/Mutagenic_primer_design.png' width='600'>

Mutagenic primer design.<br>
