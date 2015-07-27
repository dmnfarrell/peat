

## Introduction ##
PEAT\_SA allows you to analyse the effect of each amino acid residue on the biophysical characteristics of a protein. This is a PEAT plugin that allows users to submit to and retrieve results from the PEATSA server. Jobs are given a unique ID and their status can be checked periodically. When ready results can be retrieved and displayed in a table. Thes values can be exported or sent to a new labbook sheet or merged with an existing one, alongside experimental values with matching mutation field.

## Usage ##
### Submit jobs ###
The submit job button is used to create new jobs and allows you to provide the mutations, pdb file and ligand file for binding calculations.
You need to provide a structure and some mutations at minimum. If your DB has a reference protein (often a wild-type), you can choose to use this. Otherwise you can load one from disk or type the 4 letter PDB code and it will be retrieved from the pdb site.
Mutations can be typed in the dialog box and/or loaded from a text file. Records that have correctly filled mutations field can also be directly placed into the dialog.

The mutation format required is that used in PEAT, of the form:
A10A or A:0010:ALA and A10A+A11N for multiple mutations
Job name should be something that you will recognise later on, so you can retrieve the correct job results. Once jobs are submitted the results will be kept on the server indefinitely, until you remove them. All that PEATDB stores is a reference to the ID on the server.
You may optionally tell PEAT that one of your columns contains the reference experimental data you want to correlate to the predictions. This will be saved and used to make plots when the results are ready.

### Getting results ###

To check a job is finished you can periodically return and view details to show the status. When the job is done, select it from the list of jobs and click 'View Results'. If you have previously provided the column containing the reference experimental data when submitting, PEAT will remember this and automatically plot the correlation for you. Otherwise you will be asked for the column. If you don't yet have this data in the DB, you can enter/import it into PEAT and always retrieve the results later.
The 'Manage Results' command brings up another dialog with a table of results for each mutant is displayed. You can save this as a csv file, send it to a new labbook sheet in your current db.
Most usefully, you can also merge the results with table (labbook) of corresponding mutations that may have experimental values, for quick comparison.

### Viewing a Correlation ###
When the plugin has a list of exp vs predicted values, it can make a correlation plot. (This function can also be accessed separately from the plugins menu for general use). The correlation will always be show with predicted values on the x-axis vs your chosen experimental data column on the y-axis. The points are clickable. When selected the mutation code and name of the corresponding mutant is displayed in the plot.
The function of the plot is to allow the experimentalist to see not only how well the predictions work, but for any point what are the main effects causing the (de)stabilising effects, for example.

## Configuration ##

PEATSA uses a configuration file usually kept in the users home directory. This file contains paths to dependencies/python classes used for running a full PEATSA installation. The first time you open the plugin this file is created if it doesn't exist. For using the plugin the only section of the file to worry about is the DATABASE section. You can edit the file from the plugin dialog by clicking the 'create/edit conf file' button. Make sure the database sections reads as follows. Then save the file:
```
[DATABASE]
database =  DBSAInterface
user = peatdb
password = 123
host = enzyme.ucd.ie 
```

## Requirements ##

You should install the PEATSA python classes alongside PEATDB if you don't already have it. The PEATSA package is available on svn using:
```
svn co http://peat.googlecode.com/svn/trunk/PEATSA peat-read-only
svn co http://peat.googlecode.com/svn/trunk/UFFBAPS peat-read-only
```

See Also

If you do not wish to use this plugin you may also use the [PEAT\_SA webinterface](http://peat.ucd.ie/ProteinDesignTool/Pages/FrontPage.php)