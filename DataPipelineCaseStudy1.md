## DataPipeline Case Study 1: Spectral data processing ##

To illustrate the filtering capabilities of DataPipeline we generated simulated spectral data such as is found in mass spectrometry type experiments.

<img src='http://peat.googlecode.com/svn/wiki/images/datapipelinecasestudy1.png' align='right' width='400'>

The purpose was to demonstrate that basic peak detection could be carried out using a number of the pre-defined functions that are available by default. The simulated data were generated as follows:<br>
<ul><li>A selection of Gaussian peaks at random intervals and widths is created<br>
</li><li>A changing baseline, given by a power-law as a function of x is added<br>
</li><li>Normally distributed noise is added to the baseline (in this case with a standard deviation of 8% of the baseline value)</li></ul>

The number of peaks, widths and range of points is not important for our purposes, the success of the algorithm mainly being a function of the noise present. The result for one set of points is shown below, with the data at each step illustrated. The functions were applied in the following order:<br>
<ul><li>Smoothing: a fixed window smoothing filter with window length adjusted to the noise characteristics.<br>
</li><li>Baseline removal: the function baselinecorrection removes the baseline by dividing the data into segments of equal size, using the minimum of the y values for that range, then interpolating over the input x data again to yield the baseline.<br>
</li><li>Peak detection: the function detectpeaks finds peaks by searching for values which are surrounded by lower or larger values for maxima and minima respectively. The lookahead parameter determines the distance to look ahead from a peak candidate to then next one. The delta parameter specifies a minimum difference between a peak and the following points, before a peak may be considered a peak. This is useful to hinder the algorithm from picking up false peaks towards the end of the signal.</li></ul>

The following settings were used in the functions.conf file for the function parameters:<br>
<pre><code>[detectpeaks]<br>
lookahead = 3<br>
delta = 20<br>
[smooth]<br>
window = hanning<br>
windowlen = 7<br>
[baselinecorrection]<br>
segmentsize = 0.3<br>
</code></pre>

Note that this is a very basic example with test data and real experimental datasets would need to be treated perhaps with more specific pre-preprocessing filters. But the general concept remains the same. Also note that the peak maxima are not conserved with the fixed-window smoothing filter we have used. In practice a continuous wavelet transform can be applied to the data that accounts for the varying peak widths and also helps to avoid false-positives due to noise peaks. Such an algorithm could be easily added as a filter to DataPipeline.<br>
<br>
See also <a href='DataPipelineCaseStudy2.md'>DataPipelineCaseStudy2</a>