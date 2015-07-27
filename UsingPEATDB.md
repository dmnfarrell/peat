# Working with PEATDB #

This article gives a general introduction to working with PEATDB (or just PEAT). The client is designed with protein engineering in mind, but in principle can be used for other purposes that require data storage and visualization. The plugin system allows users to write their own functionality using the API.

## Introduction ##

The PEAT window is arranged rather like a spreadsheet in a table. The application allows rudimentary functions like sorting, re-ordering of columns, resizing of columns, recoloring cells and so on. For large numbers of records (i.e. >1000), table rendering will get slower and it is recommended to use a paging view. This can be selected from the menu.
Table preferences are accessible from the right click menu or from the main menu, settings->table prefs. This allows users to control the size of fonts, row heights etc.

## Database Functionality ##

PEAT uses an object database backend called ZODB. This allows us to store objects without using a traditional RDBMS like MySQL. Users can choose to create a local database file or connect remotely to a server storing the DB. Before connecting to the server, the database has to be created by someone with the correct privileges. By default PEAT uses the relstorage backend. See more details here on setting up a server.

## Creating Databases ##

For a local database, you can simply choose 'new DB' and save the file as you would any document. All changes are made locally and cannot be shared with other users at the same time. This is single user mode.
For a remote/shared database. (This refers primarily to the relstorage backend as that is default). The project can be either setup on the mysql server using an interface like phpmyadmin or created directly in the application. This requires that the user access to a mysql account (but NOT root) that allows them to create new databases. See Creating a new DB in PEAT.

### Database size ###
You can manage the database size by periodically doing 'packing'. The database grows larger partly because it retains multiple versions of older transactions, which is what allows undo. However you will find after some time that old transactions can't be undone anyway, so you may want to remove versions older than a certain time. Those changes will be removed from the change log. This is called packing the database. It can be done from the application, called from the settings menu.
For the default relstorage backend, packing can also be done using the zodbpack script called outside peat. The script requires a configuration file to specify the DB. Here is a typical config file:

```
<relstorage>
 pack-gc true
 pack-duty-cycle 0.9
 <mysql>
   db zodb
   user username
   passwd password
 </mysql>
</relstorage>
```

The script is called using:
```
zodbpack -d <days> <conffile>
```

### Size Limits ###
There is no limit in principle - except the size of your filesystem and RAM (for keeping the object index). Single filestorages with 20GB and more are common practice. But of course: never forget to pack the storages on a regular basis. In terms of practicality, the PEAT client may have some specific issues: such as searching a database of 100,000 records; there will be some delay as the field data is retrieved from storage to the client and parsed.

## View database changes ##

You can view the records that have been changed since your last save by choosing View->Show Changed Items. Currently this simply gives a list of the records you have changed.

## Handling Conflicts ##

PEAT is not a relational database system (RDBMS). Unlike a database, the ZODB mode of operation is nonlocking (i.e. it uses optimistic concurrency control). That means that if two users are accessing the same record item, nothing prohibits both of them from making changes to that file. In fact, this method generally works quite well, assuming a reasonable degree of organisation between and within groups collaborating on the same data concurrently. In the rare cases when a conflict occurs, users are notified and forced to manually intercede and resolve the conflict in the context of the particular data. PEAT provides the user with the choice of reverting their changes and then accepting or overwriting the conflicting changes.

## The Change log ##

This shows a record of all the changes made to the DB in reverse order. It shows who made the transaction and when and a comment if available. Users should make comments when saving to the DB, if they are collaborating with other people who might need to see what has been changed later. Some changes can be undone, if no subsequent commits have altered the same records. The undo is counted as a new revision and an entry for this transaction also appears in the change log.

## Storing Binary Files ##

PEAT allows you to store binary files to the Database (called blobs in ZODB terminology). These could be images, zip, pdf files, word documents or even videos that are too large to insert in the database the normal way. Most any file type is supported.
The files are attached to specific record/field, but stored in the file system directly. PEAT includes a column type 'file' specifically designed for this. You can also add files inside the Labbook. Each cell corresponds to one file.
For local projects the files are kept in a folder called projectname.fs.blob where projectname is the name of your DB. If you move the project files elsewhere, obviously you need to move this folder too.
For remote Databases the files are stored remotely but a copy of the files you download when you view the files is kept locally too. Any files you add will therefore be available to other users connecting to the DB. (See PEATDB Settings and PEAT Server Setup)
You can also import multiple files from a folder at once. Choose import external/binary fileset. The dialog asks you to select a directory and filter by the type of file to use if needed. Simply click 'Import' when ready. Note that the DB might take some moments to save when you have imported many large files.

## Searching ##

Simple searching of the database can currently be carried out by filtering on the table. Multiple filters can be applied joined together with logical operators. The user interface is intuitive and self explanatory. Filtering currently works on simple column types like text, number. For more advanced field types such as ekin, see PEAT Advanced Search.

## Comparing Records ##

See article on Comparing Records in PEAT

Importing Data

You can currently import CSV (comma seperated values) files into PEAT. A dialog will be shown allowing you to tell the program what the formatting on your data is, allowing some flexibility on what PEAT can handle. More complex importing of other data could be handled by writing your own plugin.

## Memory Usage ##

Loading many objects into memory will eventually cause swapping if the DB is large. Therefore each database has a cachesize setting to limit the memory usage. It should be set appropriate to your needs, the application will not be able to easily guess this value, so it is left to the user. Usually the person managing or creating the database should set this and other users won't have to worry about it. If you have many small-medium sized records you can make the value fairly high (~500) for efficiency. If you have very large records, set it low (~50). You can work out a value based on your own RAM size. Each record/field item roughly corresponds to one object.
To change the cache size, go to Settings->DB Specific Settings and enter a new value and press Apply. You will need to save the database for the change to be remembered.

## Sample Projects ##

|name|purpose|address:project:port|
|:---|:------|:-------------------|
|hewlsample|contains NMR and sequence data for Hen Egg White Lysozyme|enzyme.ucd.ie:hewlsample:8080|
|filestest|contains arbitrary binary file data|	enzyme.ucd.ie:filestest:8080|
|Titration\_DB|NMR titration curves for many proteins|enzyme.ucd.ie:titration\_db:8080|
|Large|Large test set of ~5000 randomly generated records with some dummy fields|enzyme.ucd.ie:large:8080|

## Plugins ##

See PEATDBPlugins.