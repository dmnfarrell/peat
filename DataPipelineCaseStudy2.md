## DataPipeline Case Study 2: Grouping multiple files by their field names ##

This example shows how to import the data from many files into one by using the group keyword. The data in this case were pH NMR titration data files. Each file represents the data for one protein residue and each contains 3 columns: the x-data (pH points) and a column for 1H and 15N measurements.

Example contents of one file:
```
7.913 127.098 8.836  
7.494 127.122 8.839  
7.113 127.146 8.843  
7.067 127.156 8.839  
6.594 127.167 8.841  
6.379 127.171 8.845  
6.351 127.176 8.852  
6.143 127.177 8.844  
5.849 127.192 8.850  
5.803 127.186 8.854  
5.653 127.205 8.856  
5.400 127.222 8.860  
5.192 127.218 8.866  
5.139 127.233 8.866  
4.903 127.254 8.871  
4.730 127.260 8.875  
4.539 127.270 8.879  
4.293 127.269 8.879  
4.063 127.263 8.882  
```

### Configuration ###

This kind of data is straightforward but it needs a few special configuration settings:

  * The format is databycolumn, but there are no column headings. So we set these manually using the _colheaderlabels_ keyword. Note that adding only the first column name means the second is ignored.
  * The other objective is to group the results into files for the 1H and 15N columns respectively. This is achieved using the _groupbyfields_ keyword. This causes the program to attempt to put the data from each file into a (usually smaller) number of files based on their field labels (column or row). Without this step we would end up with a file for each input, which might not be convenient.
  * The files are named with a .inp  extension and we need to put this in the configuration so the program picks them up.

The required configuration file has these settings:

```
[base]
colheaderlabels = 15N,1H
groupbyfields = 1
[files]
extension = .inp
```

The result is a csv file and ekin project file for each column.

### Sample data ###
Sample data files for this case study and the appropriate configuration file are available [here](http://peat.googlecode.com/files/nmr_data.zip).