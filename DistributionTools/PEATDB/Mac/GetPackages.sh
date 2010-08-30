#! /bin/bash
#This script gets relstorage, ZODB and all zope packages, and MySQLdb
#The following two packages must be built manually
#   matplotlib 0.99 - not on easy install
#   numpy - easy install goes crazy on this
#Following two lines required or MySQL-python will fail
export PYTHONPATH=$PWD/ExternalPackages/:$PYTHONPATH
export PATH=/usr/local/mysql/bin/:$PATH
#Ensure no zope packages are installed in standard locations
#If there are there will be import problems - its possible this can be resolved by changing
#py2app build flags but haven't figures this out yet
/usr/bin/easy_install-2.5  --install-dir=$PWD/ExternalPackages relstorage
/usr/bin/easy_install-2.5  --install-dir=$PWD/ExternalPackages MySQL-python
