<h1>DataPipeline</h1>
<img src='http://peat.googlecode.com/svn/wiki/images/PipeLine_logo.png' align='right' width='300'>
<br>
<h2>Introduction</h2>

DataPipeline is a python desktop and command line application that uses the fitting and plotting libraries from PEAT to automate the import of raw data in a variety of formats. It is also used for transforming the imported data through various stages of fitting to achieve a final measured parameter. This is useful for handling large amounts of data (e.g from csv files) in a consistent way without having to store everything in spreadsheets.<br>
<br>
<h2>Rationale</h2>
Raw data from biology, chemistry and just about any laboratory experiments comes in a large variety of formats. Often the data format is proprietary but can often be converted to plain text or csv files for processing by the user. The problem is that the user then has to put the data into a spreadsheet and make sense of it. For one file this might be trivial, but for a lot of files that need processing it becomes laborious. Even when automated in a spreadsheet the workflow can be very confusing.<br>
<br>
<h2>Features</h2>

<ul><li>configuration file provides flexible import of raw data in text format<br>
</li><li>file names can be parsed and grouped by their labels (see examples)<br>
</li><li>possible to fit raw data and then chain these results to a further round of fitting<br>
</li><li>add your own non-linear fitting models using a special module<br>
</li><li>results output as csv files and plots if desired<br>
</li><li>works from command line as well as desktop application<br>
</li><li>programmers can add their own importers<br>
</li><li>apply pre-defined or custom filters/functions to the data set and chain them together to form a pipeline</li></ul>

<img src='http://peat.googlecode.com/svn/wiki/images/datapipeline_workflow.png' align='right' width='500'>

<h2>Workflow</h2>

Shown at left is the general workflow for the software. Boxes in dashed outline are optional steps.<br>
<br>
<h2>Formats</h2>

The majority of text readable experimental data from, say, biochemical assays have formatting that can be categorised into pre-defined foramts. For example a table of data can either be presented with independent datasets arranged in '''rows or columns'''. Other variations are then subsets of these cases. (For our purposes 'dataset' refers to a single set of x-y points. The whole purpose of the exercise is to extract these x-y points.) The format is specified in the configuration file under the 'format' keyword, see Configuration. The typical formats you would generally expect to find are built into DataPipeline, that should cover a large percentage of possible cases that will arise. See DataPipelineFormats.<br>
<br>
The chart below illustrates the concept graphically. You will be quite quickly able to see which one of these cases corresponds to your own data file. More complex combinations are possible of course, but this will be rarely the case for most users. For very specific and usual formats a custom importer class can be written to handle it and integrated into the application. Knowledge of python is required for this, see the API section below.<br>
<br>
<img src='http://peat.googlecode.com/svn/wiki/images/Dataformats_overview.png' />

<h3>Grouped Data</h3>
A special case that is sometimes encountered in biological assays are files with data grouped by some experimental condition for a range of samples. In this case the data is multidimensional.<br>
<br>
<h2>Configuration</h2>

A default configuration file is written by to the .pipeline folder in the users home directory and custom files can be adapted from this one and placed wherever the user wishes. A typical configuration file is as follows:<br>
<br>
<pre><code>[base]<br>
format = databyrow<br>
rowstart = 0<br>
colstart = 0<br>
rowend = 0<br>
colend = 0<br>
colheaderstart = 0<br>
rowheaderstart = 0<br>
rowheader = <br>
colheader = <br>
colheaderlabels = <br>
rowheaderlabels = <br>
rowrepeat = 0<br>
colrepeat = 0<br>
delimeter = ,<br>
workingdir = /home/user/.pipeline/workingdir<br>
functionsconf = /home/user/.pipeline/functions.conf<br>
ignorecomments = 1<br>
decimalsymbol = .<br>
preprocess = <br>
xformat = <br>
yformat = <br>
groupbyfields = 0<br>
<br>
[files]<br>
groupbyname = 0<br>
parsenamesindex = 0<br>
parsemethod =<br>
replicates = 0<br>
extension = .txt<br>
<br>
[fitting]<br>
xerror = 0<br>
yerror = 0<br>
iterations = 50<br>
modelsfile = <br>
<br>
[models]<br>
model1 = <br>
<br>
[variables]<br>
variable1 = <br>
<br>
[functions]<br>
function1 = <br>
<br>
[excel]<br>
sheet = 0<br>
numsheets = 1<br>
<br>
[plotting]<br>
saveplots = 0<br>
fontsize = 9<br>
normalise = 0<br>
grayscale = 0<br>
alpha = 0.8<br>
font = sans-serif<br>
markersize = 25<br>
linewidth = 1<br>
showerrorbars = 0<br>
dpi = 100<br>
marker = o<br>
markers = -,o<br>
legend = 0<br>
<br>
[custom]<br>
<br>
</code></pre>

<h3>Explanation of selected options</h3>

Some options are self-explanatory, such as the plotting options. The functions section is explained in DataPipelineFunctions. The remainder are detailed below.<br>
<br>
<table><thead><th> <b>option</b></th><th> <b>possible values</b> </th><th> <b>explanation</b> </th></thead><tbody>
<tr><td>workingdir    </td><td>any valid path          </td><td>the folder where all results are placed</td></tr>
<tr><td>format        </td><td>see formats diagram above</td><td>general structure of the data fall into predefined categories</td></tr>
<tr><td>decimalsymbol </td><td>. or ,                  </td><td>symbol used to indicate decimal point</td></tr>
<tr><td>delimiter     </td><td>any symbol except a numerical value</td><td>separator between data</td></tr>
<tr><td>rowstart      </td><td>any integer value       </td><td>row where the data starts, including x labels</td></tr>
<tr><td>colstart      </td><td>any integer value       </td><td>column where the data starts, including y labels</td></tr>
<tr><td>rowend        </td><td>any integer value       </td><td>row where the data ends, optional</td></tr>
<tr><td>colend        </td><td>any integer value       </td><td>column where the data ends, optional</td></tr>
<tr><td>rowheaderlabel</td><td>comma separated list    </td><td>provide row names manually</td></tr>
<tr><td>colheaderlabel</td><td>comma separated list    </td><td>provide column names manually</td></tr>
<tr><td>colrepeat     </td><td>0 or any value >1       </td><td>indicates that sets of data are grouped in evenly spaced columns</td></tr>
<tr><td>rowrepeat     </td><td>0 or any value >1       </td><td>indicates that sets of data are grouped in evenly spaced rows</td></tr>
<tr><td>ignorecomments</td><td>0 or 1                  </td><td>ignore lines starting with #</td></tr>
<tr><td>xformat       </td><td>see 'specifying time formats' below</td><td>specify if x values in time formats</td></tr>
<tr><td>yformat       </td><td>as above                </td><td>specify if y values in time formats</td></tr>
<tr><td>xerror        </td><td>any decimal             </td><td>experimental error on all x values</td></tr>
<tr><td>yerror        </td><td>any decimal             </td><td>experimental error on all x values</td></tr>
<tr><td>groupbyfield  </td><td>0 or 1                  </td><td>tells the program to group input data by field names, see below</td></tr>
<tr><td>iterations    </td><td>any value greater than 1</td><td>number of rounds of fitting</td></tr>
<tr><td>modelsfile    </td><td>any valid filename      </td><td>file to load you models from</td></tr>
<tr><td>models        </td><td>any valid model name    </td><td>model to use for fitting see ModelDesign</td></tr>
<tr><td>groupbyname   </td><td>0 or 1                  </td><td>activate the grouping by file name functionality, explained below</td></tr>
<tr><td>extension     </td><td>any file extension      </td><td>filter input files by extension, default is .txt</td></tr>
<tr><td>parsenamesindex</td><td>any valid index or empty</td><td>index(es) to parse filenames for data labels</td></tr>
<tr><td>parsemethod   </td><td>numeric,text or both    </td><td>method to parse filenames, will detect only numbers, text or both</td></tr></tbody></table>

<h3>Specifying time formats</h3>
The xformat and yformat keywords allow you to specify if either dimension of your data are in units of time. These can then be converted to integer values of seconds for numerical analysis. For example if you x values increasing in amounts of 10 seconds in the folllowing format: 0:10, 0:20, 0:30. You enter <b>xformat = %M:%S</b> in the configuration file.<br>
<br>
<h3>Grouping by field names</h3>

The <b>groupbyfield</b> parameter tells the program to try to re-arrange the imported data so that mutliple files are grouped into one or more field names. The field names are the row or column header labels depending on you input format. In this way you can combine many individual input files into a one or more output csv files (and an ekin project). An example is given in DataPipelineCaseStudy2.<br>
<br>
<h3>Grouping by file labels</h3>

This functionality is enabled by setting the <b>groupbyname</b> keyword to 1. It allows data to be re-grouped according to a text or numerical label that can be extracted from the file names. The index of this label is specified by the <b>parsenamesindex</b> keyword and <b>parsemethod</b> tells the parser whether to extract numbers or text only or alphanumeric words.<br>
<br>
Example:<br>
<pre><code>Parsing the filename 'AA_2-12-12_ph8.0_rep2.txt'<br>
<br>
The options parsenamesindex=3, parsemethod=numeric will detect the label '8.0'<br>
The options parsenamesindex=0, parsemethod=text detects the label 'AA'.<br>
</code></pre>

Note that currently only a single file label can be parsed but in future multiple indexes will allow data to be grouped by more than one label.<br>
<br>
This functionality may be used in conjunction with some model fitting. For example if we had a set of files named according to a temperature the data was taken at and each file contained some other assay that could be fit to some model. We would parse the file names, import, fit each file data and regroup to achieve a set of temperatures vs. the fitted parameters. See also DataPipelineCaseStudy3 for a worked example.<br>
<br>
<h2>Model Fitting</h2>

Part of the application involves fitting imported data and grouping the fitted parameters into new sets of data which can in turn be fit to another model. Custom fitting models can be created using the ModelDesign module, where further information is provided.<br>
<br>
<h2>Filtering and pre-processing</h2>

See DataPipelineFunctions for information on specifying custom filters for the data. See also DataPipelineCaseStudy1.<br>
<br>
<h2>Desktop Application</h2>

The desktop application has a simple layout that shows a preview of the currently loaded file, a log window, file queue and a plot preview area that shows how the data will be formatted after import for the currently loaded file. Below is a screenshot of the main window. The application includes access to the ModelDesign program, a simple text editor for editing the configuration file and a batch file renaming dialog. The currently loaded files and configuration can be saved together as a project for convenience. The program can also be executed from the command line if a desktop environment is not available.<br>
<br>
<img src='http://peat.googlecode.com/svn/wiki/images/datapipeline_scr1.png' />

<h2>API</h2>

DataPipeline is written in python and uses modules from the PEATDB project for plotting and fitting. For those interested in the programming interface and adding new importers, see DataPipelineAPI.<br>
<br>
<h2>Other links</h2>

DataPipelineInstallation<br>
<br>
<h2>Citing</h2>

Please reference the software using this citation:<br>
<br>
<code>DataPipeline: Automated importing and fitting of large amounts of biophysical data. D Farrell and JE Nielsen. J.Comput.Chem. 2012 DOI: 10.1002/jcc.23066</code>

Link to this paper: <a href='http://onlinelibrary.wiley.com/doi/10.1002/jcc.23066/abstract'>http://onlinelibrary.wiley.com/doi/10.1002/jcc.23066/abstract</a>

A pdf of this paper is available <a href='https://www.researchgate.net/publication/229160677_DataPipeline_Automated_importing_and_fitting_of_large_amounts_of_biophysical_data/file/9fcfd5016bd33cfb20.pdf?ev=pub_int_doc_dl'>here</a>