#!python
#
# $Id: install_win_pKaTool.py 208 2007-01-31 12:45:26Z nielsen $
#
# Windows install script for pKaTool
#

import os, sys, string

def create_shortcut_safe(target,description,link,args,workdir,icon):
    if os.path.isfile(link):
        os.unlink(link)
        file_created(link)
    create_shortcut(target,description,link,args,workdir,icon)
    return

def install():
    #
    # Find the install dir
    #
    installdir=None
    import sys, os
    for dir in sys.path:
        rdir=os.path.join(dir,'pKaTool')
        if os.path.isdir(rdir):
             installdir=dir
             break
    if not installdir:
        return
    #
    # Perform post installation tasks
    #
    ostyp=None
    var_s=['OSTYPE','OS']
    for var in var_s:
        if os.environ.has_key(var):
            ostyp=string.lower(os.environ[var])
    if not ostyp:
        return
    #
    # Windows?
    #
    if string.find(ostyp,'windows')!=-1:
        #
        # Yup
        #
        #
        # Create a shortcut on all desktops
        #
        pwd=os.getcwd()
        desktop=get_special_folder_path("CSIDL_COMMON_DESKTOPDIRECTORY")
        #os.chdir(desktop)
        pKa_dir=os.path.join(installdir,'pKaTool')
        mainscript=os.path.join(pKa_dir,'pKaTool.py')
        print 'pKaTool has been installed in %s.' %(pKa_dir)
        #
        # Find the location of the python executable
        #
        print
        print 'Creating shortcuts...',
        sys.stdout.flush()
        #python_exec=sys.executable
        choices=['python.exe','python']
        for choice in choices:
            exec_try=os.path.join(sys.prefix,choice)
            if os.path.isfile(exec_try):
                python_exec=exec_try
                break
        
        #
        # Desktop shortcut
        #
        target=os.path.join(desktop,'pKaTool.lnk')
        icon=os.path.join(pKa_dir,'pKaTool.ico')
        create_shortcut_safe(python_exec,'pKaTool',target,mainscript,desktop,icon)
        file_created(target)
        #
        # Put it in the Windows start menu
        #
        startmenu=get_special_folder_path("CSIDL_STARTMENU")
        pKa_startmenu_dir=os.path.join(startmenu,'programs','pKaTool')
        if not os.path.isdir(pKa_startmenu_dir):
            os.mkdir(pKa_startmenu_dir)
            directory_created(pKa_startmenu_dir)
        #
        # Shortcuts
        #
        pKa_link=os.path.join(pKa_startmenu_dir,'pKaTool.lnk')
        icon=os.path.join(pKa_dir,'pKaTool.ico')
        create_shortcut_safe(python_exec,'pKaTool',pKa_link,mainscript,desktop,icon)
        #
        # pKa system shortcut
        #
        pKasys_link=os.path.join(pKa_startmenu_dir,'pKa_system.lnk')
        icon=os.path.join(pKa_dir,'pKaTool.ico')
        create_shortcut_safe(python_exec,'pKa_system',pKasys_link,os.path.join(pKa_dir,'pKa_system.py'),desktop,icon)
        #
        # EM_effect
        #
        print 'Creating EM_effect shortcut'
        EM_effect_link=os.path.join(pKa_startmenu_dir,'EM_effect.lnk')
        icon=os.path.join(pKa_dir,'pKaTool.ico')
        create_shortcut_safe(python_exec,'EM_effect',EM_effect_link,os.path.join(pKa_dir,'EM_effect.py'),desktop,icon)
        print 'Done'
    else:
        #
        # Nothing to do
        #
        pass
    return


def remove():
    return


# main()
if len(sys.argv) > 1:
    if sys.argv[1] == '-install':
        install()
    elif sys.argv[1] == '-remove':
        remove()
    else:
        print "Script was called with option %s" % sys.argv[1]

