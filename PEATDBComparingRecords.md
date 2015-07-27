# Comparing Records in PEATDB #

Numerical tabulated data is stored in PEATDB using our [[Ekin](Ekin.md)] format, these data can be plotted and compared in any combination using the compare records function, called from the Tools menu. This is useful for example when you wish to compare a plot for some residue for multiple mutants with the wild type. Note that a 'dataset' here refers to one set of X-Y data points, i.e. one plot.

## General Usage ##

  * Simply select the record(s) and field to compare and then use the rightmost dialog buttons to select the particular 'datasets' that you want to filter out. This is done by entering a wildcard value in the data labels box, or the exact name of the label you want to plot. Pressing return or 'preview' shows which labels are found, if any.
  * A ''dataset'' constitutes an individual label for each x-y data inside an ekin project. So multiple records may have the same field holding the same dataset. We use this functionality when comparing wt vs mutant for properties tracked across the same residue.
  * Press display to view the plots, which can be overlayed in one plot (for each field) or in many sub plots.
  * If you select more than one record/field, the plots can be grouped by field OR record, depending on how you want to slice the data.

## Comparing across a number of records and fields at once ##

One of the purposes of this module is to allow users to view a possibly large dataset in ways they might otherwise not be able to look at it. Therefore you might want to compare between several fields on 1 record or you might want to compare between several records. In later case you would want to group the plots by record. In the former you would group by field. Use the '''group by''' radiobutton to choose which way to slice the data.

The figure shows comparisons of several datasets for 2 a wt and mutant taken for 3 fields. This example shows plots of 2 kinds of NMR titration data for several residues and a field with kinetics data. Note that all of these datasets are fitted to models in ekin, the red curves.

## Annotate Plots ##

You may add basic annotations to the plots from the top toolbar menu.

## Publication Quality Plots ##

Since we use matplotlib for plotting, it is possible to produce plots of high quality.
  * Appearance and rendering may be improved if you tick the '''use latex''' option inside the plot options dialog. This will use LaTex to produce properly rendered greek symbols and math expressions, however drawing of plots will take somewhat longer, so be patient.
  * Saving the plots at 300 dpi is usually required for publications. See example in the gallery below.

## Using LaTex ##

This requires that some supporting bits of software are present on your system. In linux, usually these dependencies will be installed with the matplotlib packages. On windows you will need dvipng, a working LaTex install and ghostscript. See http://matplotlib.sourceforge.net/users/usetex.html
We recommend you install MikTex if you wish to have access to latex on windows. Download at http://miktex.org/2.8/setup

## Examples ##