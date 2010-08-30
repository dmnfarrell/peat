#! /bin/bash
PEATSOURCE=$1
PACKAGEDIR=$1/DistributionTools/PEATDB/Mac/
EXTERNALPACKAGES=$PACKAGEDIR/ExternalPackages/
VERSION=$2
export PYTHONPATH=$EXTERNALPACKAGES:$PEATSOURCE:$PYTHONPATH

echo -e "Building package PEATDB${VERSION} from $PEATSOURCE in $PACKAGEDIR \n"

echo Updating source ...
#svn update $1
cd $PACKAGEDIR
rm -rf build/ dist/ PEATDB/ 
echo -e "\nBuilding application ..."
python2.5 setup.py py2app --dist-dir PEATDB -P PEATDB-Info.plist -i PEATDB,Protool,Pmw,MySQLdb,ZODB,relstorage,matplotlib --use-pythonpath &> buildoutput
echo -e "\nCreating disk-image..."
hdiutil create PEATDB${VERSION}.dmg -srcfolder PEATDB -volname PEATDB$VERSION -ov
echo -e "\nInternet-enabling disk-image..."
hdiutil internet-enable -yes PEATDB${VERSION}.dmg 
echo -e "\nDone"

