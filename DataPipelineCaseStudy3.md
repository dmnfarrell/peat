## DataPipeline Case Study 3: Custom handling of multi-dimensional Kinetics data ##

### Introduction ###
The advent of high throughput protein and enzyme characterization technology has created the possibility of performing accurate biophysical measurements for a very large number of mutant proteins. The following example formed part of such a study: we studied the pH-dependence of kcat and KM for an enzyme using measurements of initial rates at 8 substrate concentrations at each of 10 pH values in triplicate. The ultimate goal was to track the change in kinetic parameters kcat/KM with pH in order to detect changes in the pKa value at the active site. We therefore had to fit three successive models to the data.

<img src='http://peat.googlecode.com/svn/wiki/images/datapipelinecasestudy3.png' align='right' width='500'>

<h3>Format</h3>
Since each initial rate measurement consisted of 90 absorbance/concentration measurements each enzyme clone resulted in 90x8x10x3 = 21,600 data points to be analyzed per mutant. The entire data set consisted of 100 mutant enzymes thus giving rise to approximately 2 million data points, which all must be analyzed consistently to give accurate and comparable values for the parameters of interest. The study was previously done using a set of spreadsheets with macros for plotting and fitting, resulting in speed issues and problems when any changes needed to be made by the actual users.<br>
<br>
<h3>Workflow</h3>
The workflow was as follows:<br>
<ol><li>Import the time vs. absorbance data at each substrate concentration for every variant, in triplicate.<br>
</li><li>Fit every raw dataset to a linear model to obtain a velocity at every concentration<br>
</li><li>Grouping these velocities versus their corresponding substrate concentrations gives us a dataset that we can fit to the Michaelis-Menten model and a yield KM value. Combined with a known enzyme concentration we can find kcat/KM<br>
</li><li>The final step was to group the kcat/KM results by pH and fit to the Henderson-Hasselbalch equation and extract a pKa value.</li></ol>

To handle the data a custom importer was written since it did not match any of our standard formats. (A more general version of this importer will be added as one of the standard formats). Each data file contained a set of columns corresponding to 12 variants with absorbance values for each time point on the y-axis. However the time points were also grouped for the 8 substrate concentrations, giving multi-dimensional data per file which needed to be extracted. The figure right shows how the data looked and also gives an overview of how the data was organized for fitting.<br>
We used the fit propagation functionality described elsewhere. Replicates were treated as independent and therefore fit separately. Note also that in our worked example we fit to KM values since finding the kcat requires an intermediate calculation step which we performed with a small amount of extra custom code. We are currently working on a more flexible way to specify automatic intermediate processing steps between successive iterations of fitting, in addition to the pre-processing step mentioned already.<br>
<br>
<h3>Configuration</h3>
The relevant configuration settings are shown here. Fit propagation is done by enabling the <b>groupbyname</b> keyword along wit specifying the required models and parameters. Also not the <b>xformat</b> option specifying that these are in a time format.<br>
<br>
<pre><code>[base]<br>
format = kineticsdata<br>
rowstart = 3<br>
colstart = 0<br>
rowend = 0<br>
colend = 12<br>
rowheader = 3.2,1.6,0.8,0.4,0.2,0.1,0.05,0.025<br>
colheader = wt 5,wt 3,wt 2,68 5,68 3,68 2,138 5,138 3,138 2,248 5,248 3,248 2<br>
rowrepeat = 9<br>
delimeter = tab<br>
decimalsymbol = ,<br>
xformat = %M:%S<br>
<br>
[files]<br>
groupbyname = 1<br>
parsenamesindex = 2<br>
replicates = 0<br>
<br>
[fitting]<br>
yerror = 0.01<br>
iterations = 50<br>
<br>
[models]<br>
model1 = linear<br>
model2 = Michaelis-Menten<br>
model3 = sigmoid<br>
<br>
[variables]<br>
variable1 = a<br>
variable2 = Km<br>
variable3 = tm<br>
</code></pre>

<h3>Sample data</h3>
Sample data for a set of 12 variants over three replicates is available for download as a zip file <a href='http://peat.googlecode.com/files/kinetics_data.zip'>here</a> and this case study can be tested in the application by loading one of the replicate folders and the configuration file included.