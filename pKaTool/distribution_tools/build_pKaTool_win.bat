cd c:\python2.5\Lib\site-packages\
del /Q build
del /Q dist
c:\python25_SVN\python.exe c:\python25_SVN\Lib\site-packages\pKaTool_setup.py bdist_wininst --install-script install_win_pKaTool.py
