# Model design in DataPipeline #

ModelDesign is a module of the PEAT suite of applications that is used to create non-linear fitting models for use in the Ekin application or elsewhere. The fitting equation, parameters and start values are entered in the dialogs as text and the model may be tested interactively on the users sample data.

The application is useful for:
  * quickly transferring and evaluating equations from the literature
  * testing the usefulness of a particular model for specific data or differing amounts of noise
  * rapidly comparing different versions of a model
  * finding the appropriate start values or guesses often required in non-linear fitting

This application is currently available in Ekin and DataPipeline or can be run standalone.

See DataPipelineModelFitting for more information on how the models are used on that application.

## Usage ##
<img src='http://peat.googlecode.com/svn/wiki/images/Modeldesignapp.png' align='right' width='600'>

The user interface is simple and mostly self explanatory. On the left side of the window are entry fields for equation, parameters, description and parameter start values. The right side shows a plot of the current test data. When the test model button is pressed a fitter is created from the current settings and used to to fit the data. The plot is updated at each iteration to show the progress. Current fitted parameters values are updated in the plot also. The concept being that if a fit cannot converge the user can interactively alter the start values to improve performance.<br>
<br>
One or more models are stored in a python dictionary which can be saved when necessary. This dictionary can be read by Ekin or DataPipeline, applications that perform model fitting. Or users can simply make note of the details and use them in their own fitting application. To create a new model, just press the 'New Model' button, name the model and some default settings are added to the entry fields.<br>
<br>
The application opens with some default test datasets. Your own test data can be added by importing a csv file or loading a previously saved Ekin project. Multiple sets of data can be imported and moved through with the navigation buttons above the plot.<br>
<br>
<h2>Writing Equations</h2>

Models are written using the python syntax, so anyone familiar with that language will find it straightforward. See <a href='http://docs.python.org/tutorial/introduction.html#using-python-as-a-calculator'>this link</a> for an introduction for using python as a calculator. The equation must include the 'x' symbol at least once since this represents the independent variable in your data. The model is then fit so as to minimize the difference between your y values and the equation for each corresponding x value. This is called least squares fitting.<br>
<br>
<h2>Equation Examples</h2>
Several basic examples are embedded into the program by default. These are meant to provide an illustration of how models are constructed. They are:<br>
<br>
<table><thead><th>name</th><th>equation</th><th>guess values</th></thead><tbody>
<tr><td>Linear</td><td>a*x+b   </td><td>a=min(y); b=max(y)-min(y)/max(x)-min(x)</td></tr>
<tr><td>Gaussian</td><td>a*exp(-(pow((x-b),2)/(pow(2*c,2))))</td><td>None        </td></tr>
<tr><td>Sigmoid</td><td>bottom+(top-bottom)/(1+exp((tm-x)/slope))</td><td>tm=(max(x)-min(x))/2+min(x); bottom=min(y); top=max(y)</td></tr></tbody></table>

<h2>Model Selection</h2>

Before deciding on curve fitting approaches for their data users may wish to determine the most appropriate models to utilize. The ability to compare two or more models is provided in the application for this purpose. This function is evoked by pressing the 'compare models' button and is currently quite basic. It performs the test with the currently loaded data and simply reports the best model with the p-value. Note that you should generally perform this test a number of times with different data points to be sure that the result is consistent and not just a result of the distribution of points in one set of data.<br>
<br>
<h3>Method used</h3>
The approach used is sometimes referred to as the <b>extra sum-of-squares F-test</b> and is detailed in the <a href='http://www.graphpad.com/index.cfm?cmd=library.page&pageID=12&categoryID=6'>Graphpad Practical Guide to Curve Fitting</a>. (An F-test is any statistical test in which the test statistic has an F-distribution under some null hypothesis). The basic idea of the test is as follows: Given two related models A and B, B is a nested version of A with one more parameter. B will always be able to fit the data at least as well as A. Typically B will give a better (i.e. lower sum of squares error) fit to the data than A. To choose between models, one needs to know if B gives a significantly better fit to the data. The method requires that models are nested, that is, related in some way and with increasing number of parameters and that the same data is used to test all models.