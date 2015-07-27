# DataPipeline functions and filters #

These are sets of predefined functions that can be applied to the data one after another in the order they are specified in the functions section of the configuration file.

## Usage ##

Predefined functions are used by specfying them by name in the _functions_ section in the order in which they are to be carried out. In this way the functions can be chained together. For example, to differentiate and then smooth the raw data, the user would specify `function1=differentiate` and `function2=smooth` in the functions section. Fitting could then be done on the resulting data as normal, if required.

## Specifying function parameters ##

Most of the filters have parameters (or 'arguments') that would normally be passed to the python method programmatically. For example a smoothing filter has a window width parameter and so on. To adjust these options a separate configuration file is provided with sections for each function and the default list of parameters. This file is called `functions.conf` and is written to the .pipeline folder the first time the program is run. Usually it is easier to edit this default file, though it can be copied to another location if required.

## Current default functions ##

| **name** | **explanation** | **parameters** |
|:---------|:----------------|:---------------|
|differentiate|find numerical derivative of y data|none            |
|smooth    |smooth using convolution of various windows with signal|window(hamming,hanning,blackman,barlett)|
|gaussiansmooth|smooth with gaussian window|degree(any integer)|
|fouriernoisefilter|low or high pass filter|                |
|savitzkygolayfilter|smoothing using polynomial regression|windowsize(odd integer),order(polynomial fit),deriv(take derivative)|
|baselinecorrection|find y data baseline using interpolation and subtract|segmentsize(integer>0)|
|detectpeaks|detect peaks in noisy data|lookahead(integer>0),delta()|
|removeoutliers|remove outliers from normally distributed data|percentile(fraction between 0 and 1)|
|normalise |normalise y data to input range|value(any integer)|
|omitrange |omit x or y range of data points|xmin,xmax,ymin,ymax - ranges to omit|

## Adding custom functions ##

For now, functions can simply be added by inserting new methods into the Processor class found in Processing.py. Also add the name of the function to the classes `predefined` attribute.
New functions should be specified in the following form:
```
def somefunction(self, x, y):
    #process values here, if x or y is unchanged, just return it
    return newx,newy
```

## Function options ##
The default options are held in functions.conf, which is written on first use of the program. The following are default parameters:
```
[differentiate]

[detectpeaks]
lookahead = 5
delta = 20

[gaussiansmooth]
degree = 6

[smooth]
window = hanning
windowlen = 11

[fouriernoisefilter]

[savitzkygolayfilter]
windowsize = 11
order = 1
deriv = 0

[baselinecorrection]
segmentsize = 0.2

[removeoutliers]
percentile = 0.95

[normalise]
value = 1

[omitrange]
xmin = 200
xmax = 400
ymin = 
ymax = 
```