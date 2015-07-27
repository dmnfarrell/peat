# Using DataPipeline from the command line #

If you have installed the program from the source distribution on linux the command `PipelineCommand` will be available at the console. This launches and executes the pipeline with the given options without having to use the desktop application at all. You might want to do this if you already have the correct configuration and raw data files ready but don't have, or need, access to a graphical desktop.

## Usage ##

```
PipelineCommand -h
```

will give you all options. Current options are:

```
Options:
  -h, --help            show this help message and exit
  -c FILE, --conf=FILE  Provide a conf file
  -f FILE, --file=FILE  Raw file
  -d DIRECTORY, --dir=DIRECTORY
                        Folder of raw files
  -p FILE, --project=FILE
                        Project file
```

## Windows ##
Windows users can also use the command line app even if they have installed the application from the MSI installer. It can be found in the installation folder (usually C:\Program files\DataPipeline). We recommend creating a text file with the command and full path and whatever options you want, as above. Full paths to files should be used.