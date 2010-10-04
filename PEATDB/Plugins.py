#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Written by Damien Farrell, Feb 2010

#see http://lucumr.pocoo.org/2006/7/3/python-plugin-system

import sys
import os

class Plugin(object):
    """Base Plugin class, should be inherited from by any plugin"""
    capabilities = []
    requires = []
    
    def __init__(self, DB=None, parent=None):
        self.DB = DB
        self.parent = parent        
        return

    def _getmethods(self):
        """Get a list of all available public methods"""
        import inspect     
        mems = inspect.getmembers(self, inspect.ismethod)
        methods = [m for m in mems if not m[0].startswith('_')]
        return methods
    
    def __repr__(self):
        return '<%s %r>' % (
            self.__class__.__name__,
            self.capabilities)

    def loadDB(self, local=None):
        from PEATDB.Base import PDatabase
        if local != None:
            self.DB = PDatabase(local=local)
        return
    
def load_plugins(plugins):
    failed = []
    for plugin in plugins:        
        try:
            __import__(plugin, None, None, [''])           
        except:
            #print 'failed to load %s plugin' %plugin
            failed.append(plugin)
    return failed

def init_plugin_system(folders):
    for folder in folders:        
        if not os.path.exists(folder):
            continue
        if not folder in sys.path:
             sys.path.insert(0, folder)
        plugins = parsefolder(folder)
        failed = load_plugins(plugins)
    return failed

def find_plugins():   
    return Plugin.__subclasses__()

def parsefolder(folder):
    """Parse for all .py files in plugins folder or zip archive"""
    filenms=[]
    homedir = os.path.expanduser("~")
    if os.path.isfile(folder):        
        #if in zip file, we need to handle that (installer distr)   
        import zipfile
        zf = zipfile.ZipFile(folder,'r')
        dirlist = zf.namelist()       
        for x in dirlist:          
            if 'plugins' in x and x.endswith('.py'):                
                zf.extract(x)
        zf.close()
        #copy plugins to home dir where they will be found
        shutil.copytree('plugins', os.path.join(homedir,'plugins'))
   
    elif os.path.isdir(folder):
        dirlist=os.listdir(folder)        
        filenm=""
        for x in dirlist:
             filenm=x
             if(filenm.endswith("py")):
                 filenms.append(filenm.strip('.py'))
        filenms.sort()        
        filenameslist=[]
        filenameslist=[os.path.basename(y) for y in filenms]        
        return filenameslist     
    
_instances = {}

def get_plugins_by_capability(capability):
    result = []
    for plugin in Plugin.__subclasses__():
        if capability in plugin.capabilities:
            if not plugin in _instances:
                _instances[plugin] = plugin()
            result.append(_instances[plugin])
    return result

def describe_class(obj):
    """ Describe the class object passed as argument,
       including its methods """

    import inspect 
    methods = []
    cl = obj.__class__
    print 'Class: %s' % cl.__name__
    count = 0   
    for name in cl.__dict__:
       item = getattr(cl, name)       
       if inspect.ismethod(item):
           count += 1
           #describe_func(item, True)
           methods.append(item)

    if count==0:
      print '(No members)'
    return methods   


def describe_func(obj, method=False):
   """ Describe the function object passed as argument.
   If this is a method object, the second argument will
   be passed as True """
   
   try:
       arginfo = inspect.getargspec(obj)
   except TypeError:
      print 
      return
   
   args = arginfo[0]
   argsvar = arginfo[1]

   if args:
       if args[0] == 'self':
           wi('\t%s is an instance method' % obj.__name__)
           args.pop(0)

       wi('\t-Method Arguments:', args)

       if arginfo[3]:
           dl = len(arginfo[3])
           al = len(args)
           defargs = args[al-dl:al]
           wi('\t--Default arguments:',zip(defargs, arginfo[3]))

   if arginfo[1]:
       wi('\t-Positional Args Param: %s' % arginfo[1])
   if arginfo[2]:
       wi('\t-Keyword Args Param: %s' % arginfo[2])
    
