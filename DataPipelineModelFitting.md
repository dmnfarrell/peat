# DataPipeline model fitting #

Part of the functionality of DataPipeline is to fit imported data and group the fitted parameters into new sets of data which can in turn be fit to another model and so on. This page discusses the method by which data is fit and errors propagated, if necessary.

## Fitting ##

Fitting of data is done using a non-linear fitting algorithm of the kind found in many software packages. Note that the speed of the fitter will vary depending on the number of points per dataset and the specific model being fit. In addition convergence of many non-linear fitting algorithms are highly dependant on good initial guess for the model. New models and initial guesses can be tested/created in the ModelDesign utility. These are saved to a file which then specified in the configuration.

## Estimating experimental uncertainty ##

To get an estimate of the effects of error on the data points used in fitting, we adopt a monte carlo approach that samples the errors from multiple random fits. This is done as follows:
  * The data points are randomly adjusted within the range of the error on x and/or y values
  * Th model is refit and parameter estimates stored
  * The above is repeated multiple times
  * We then use the standard deviation on the mean of the parameters from these runs as our uncertainty. This assumes a normal distribution of the errors.

Note that this process will extend the execution time substantially as it requires at least 10-20 iterations per dataset to get a reasonable estimate.

## Propagation of fitted parameters ##

We use a recursive fitting procedure to propagate the fitted parameters from raw datasets and fit these together in another model. The figure below shows an example. Three models are successively fit to the data. The raw data is fit to a linear model and the slopes are grouped by their labels into curves that are fit to a michaelis-menten model. These in turn are grouped to be fit to a sigmoidal function.

![http://peat.googlecode.com/svn/wiki/images/fit_propagation.png](http://peat.googlecode.com/svn/wiki/images/fit_propagation.png)

For those interested, the programming scheme for this is shown below. The data is passed to the fitting procedure as a nested dictionary, the bottom level of which contains the sets of raw data as x-y lists.

<img src='http://peat.googlecode.com/svn/wiki/images/Recursive_fitting.png' width='500'>

<h2>Propagation of parameter errors</h2>

With the above procedure, the user may often require an estimate of how the uncertainties on the original data points affect the final fitted parameters when one or more rounds of fitting are done. The errors are propagated as follows:<br>
<br>
<ol><li>use the original errors on the raw data points to estimate the experimental uncertainty on the required fitted parameter as described already<br>
</li><li>if these parameters are passed to another model, we use their uncertanties as the y errors for the next estimate<br>
</li><li>this is repeated for as many models as required</li></ol>
