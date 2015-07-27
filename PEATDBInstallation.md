# Installation instructions for PEATDB #

## Windows ##

Just download the MSI installer from the downloads section. This will place shortcuts on your desktop. (You can install the source package instead, as described below, but this is a little trickier in windows).

## Linux ##

The easiest way is to use the package manager on your distribution to install python and setuptools (or pip). You can then simply type the following at the command line:
```
easy_install PEATDB
```
or
```
pip install PEATDB
```

**Or** install the python source package from the downloads section. Unzip the file and type:
```
python setup.py install
```
Any of these commands will also install all the other required python packages like numpy and matplotlib. If the command fails and reports problems installing ZODB3, use the package manager to install python-dev (or python-devel) and gcc. Then run the command again.

## Mac ##

We do not currently support a standalone installer for OS X. But the easy\_install method above should work.