**Labbook** is a simple spreadsheet like application. The Labbook allows the user to import his organisation lab data from an excel spreadsheet and retains some spreadsheet-like functions. However it does not attempt to replicate the sophisticated functions of a spreadsheet.
Most importantly, the data can be edited and shared with other users by uploading it to the central repository when used inside PEATDB. Labbook can also be used as a standalone application.

## Main features ##

  * Multiple sheets can be held in one file
  * Multiple nested tables can be held per cell
  * Table can be viewed as single sheet or over multiple pages
  * Columns can be sorted in ascending/descending order
  * Cell background or text can be coloured
  * Fill down function
  * Move columns around
  * Rename columns
  * Sort columns
  * Supports rudimentary formulas
  * Cells can hold binary files (PEATDB only)

## Field types ##

  * text
  * number
  * file : images, pdfs, docs etc. (only useable when using labbok within PEAT)
  * table : a nested table, editing this will open another instance of labbook, which can in turn store sub tables and so on
  * ekin : for ekin data. An instance of ekin will open, you will be asked for the correct mode initially.

## Open a Labbook file (.labbook) ##

  1. Go to the "Project"  menu at the top.
  1. Select "Open".
  1. A file dialog pops up, select the file (pickle file)

## Save Labbook when in PEAT ##

To save the changes permanently with a PEAT project the Labbook data can be saved to the local database. This is simply accomplished by pressing the 'return to database' button when closing labbook. If you do not want to return the changes, simply close the window.

## Import a csv file ##

Select a file with the .csv extension this file is in a comma separated list format which can be obtained by saving an excel spreadsheet in the CSV format. The first line in the excel spreadsheet or the .csv file is treated as the column names.

## Export a csv file ##

Will export the whole visible cell data as a csv file
Import an external fileset
This function can only be carried out when using labbook inside a PEAT DB, since imported files are added to svn control and sent to the repository separately. The labbook retains references to the files in each cell and when the user clicks the cell, they can open the file. Most file types can be added and generally are treated as binary files by svn.

## Images in Cells ##

Images in labbook are treated specially, giving a thumbnail preview of the images when the mouse moves over the cell. This functionality requires the PIL imaging module.

## Screenshot ##

<img src='http://peat.googlecode.com/svn/wiki/images/Labbook_Scr1.png' width='600'>