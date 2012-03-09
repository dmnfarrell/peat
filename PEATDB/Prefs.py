# Manage personal preferences
#
# (C) Copyright Jens Erik Nielsen, University College Dublin 2005
#

import os
from Tkinter import *
import Pmw

class Preferences:
    """Manage personal preferences"""
    def __init__(self,program,defaults):
        """Find and load the preferences file"""
        
        filename = '.'+program+'_preferences'
        dirs = self.get_dirs()
        self.noprefs = False
        try:
            for ldir in dirs:
                self.peatpath = os.path.join(ldir, '.peatdb')
                fn=os.path.join(self.peatpath, filename)

                if os.path.isfile(fn):
                    self.load_prefs(fn)
                    self.save_prefs()
                    return
                else:
                    self.noprefs = True
            if self.noprefs == True:
                raise
        except:
            # If we didn't find a file then set to default and save
            print 'Did not find preferences!'
            self.prefs = defaults.copy()
            print dirs
            self.peatpath = os.path.join(dirs[0], '.peatdb')
            self.pref_file = os.path.join(self.peatpath,filename)
            self.prefs['_prefdir'] = self.peatpath
            self.prefs['_preffile'] = self.pref_file
            self.save_prefs()

            # Can we set more variables?
            # Defaults savedir?
            if os.environ.has_key('HOMEPATH'):
                self.prefs['datadir']=os.environ['HOMEPATH']
            if os.environ.has_key('HOME'):
                self.prefs['datadir']=os.environ['HOME']

            # Use 'my documents' if available
            if hasattr(self.prefs,'datadir'):
                mydocs=os.path.join(self.prefs['datadir'],'My Documents')
                if os.path.isdir(mydocs):
                    self.prefs['datadir']=mydocs
            # Always save
            self.save_prefs()
        return

    def __del__(self):
        # Make sure we save the file when killed
        self.save_prefs()
        return

    def set(self,key,value):
        # Set a key
        self.prefs[key]=value
        self.save_prefs()
        return

    def get(self,key):
        # Get a value
        if self.prefs.has_key(key):
            return self.prefs[key]
        else:
            raise NameError,'No such key'
        return

    def has_key(self,key):
        """No we have this key"""
        return self.prefs.has_key(key)

    def delete(self,key):
        if self.prefs.has_key(key):
            del self.prefs[key]
        else:
            raise 'No such key',key
        self.save_prefs()
        return

    def get_dirs(self):
        """Compile a prioritised list of all dirs"""
        dirs=[]
        keys=['HOME','HOMEPATH','HOMEDRIVE']
        import os, sys
        for key in keys:
            if os.environ.has_key(key):
                dirs.append(os.environ[key])
        #
        if os.environ.has_key('HOMEPATH'):
            # windows
            dirs.append(os.environ['HOMEPATH'])

        # Drives
        possible_dirs=["C:\\","D:\\","/"]
        for pdir in possible_dirs:
            if os.path.isdir(pdir):
                dirs.append(pdir)

        # Check that all dirs are real
        rdirs=[]
        for dirname in dirs:
            if os.path.isdir(dirname):
                rdirs.append(dirname)
        return rdirs

    def load_prefs(self,filename):        
        # Load prefs        
        self.pref_file=filename
        import pickle
        try:
            fd=open(filename)
            self.prefs=pickle.load(fd)
            fd.close()
        except:
            fd.close()
            fd=open(filename,'rb')
            self.prefs=pickle.load(fd)
            fd.close()
        return

    def save_prefs(self):
        # Save prefs
        if not os.path.exists(self.peatpath):
            os.mkdir(self.peatpath)            
        import pickle
        fd=open(self.pref_file,'w')
        pickle.dump(self.prefs,fd)
        fd.close()
        return

class preferences_dialog:

    def __init__(self,parent,parentframe=None,subset='PEAT',callback=None):
        """Open the settings dialog"""
        self.parent=parent
        if parentframe!=None:
            self.settings = Frame(master=parentframe,relief=RAISED)
            self.settings.pack(fill=BOTH)
        else:
            self.settings=Toplevel()
            self.settings.title('PEAT settings')
        self.balloon=Pmw.Balloon(self.settings)
        import os
        home = os.path.expanduser("~")
        blobdir = os.path.join(home, '.peatblob')
        if subset=='PEAT':
            variables=[['username','','textbox',''],
                       ['password','','password',''],
                       ['blobdir',blobdir,'textbox','blob directory for remote DBs, if using relstorage(mysql) this should be shared fs'],
                       ['promptforlogcomments',True,'boolean','Prompt for log comments'],
                       ['showDialogsinSidePane',True,'boolean','Show certain dialogs in sidepane by default'],
                       ['thumbsize','200','textbox','Thumbnail size for external files'],
                       ['molgraphApplication','pymol',['yasara','vmd','pymol','rasmol','other'],'Molecular graphics app'],
                       ['molgraphAppPath',True,'textbox','Path to your molecular graphics application']]

        # Put lots of choices up
        row=0
        vars={}
        big_choice={}
        self.balloon = Pmw.Balloon(self.settings)
        for varname,default,choices,helptxt in variables:
            if not self.parent.preferences.prefs.has_key(varname):
                self.parent.preferences.set(varname,default)

            # Find out which type of preference we have
            if type(choices)==type([]):
                # List of choices
                var_value=self.parent.preferences.get(varname)
                vars[varname]=StringVar()
                vars[varname].set(var_value)
                big_choice[varname]={'type':'options','choices':[]}
                #for choice in choices:
                #    big_choice[varname]['choices']=choice)

                optmenu = Pmw.OptionMenu (self.settings,
                            labelpos = 'w',
                            label_text = varname,
                            menubutton_textvariable = vars[varname],
                            items = choices,
                            menubutton_width = 10 )
                optmenu.grid(row=row,column=0,columnspan=2)
                if helptxt!='':
                    self.balloon.bind(optmenu, helptxt)
            elif choices=='boolean':
                var_value=self.parent.preferences.get(varname)
                vars[varname]=BooleanVar()
                lbl = Label(self.settings,text=varname)
                lbl.grid(row=row,column=0)
                col=1
                vars[varname].set(var_value)
                Checkbutton(self.settings,variable=vars[varname]).grid(row=row,column=col)
                col=col+1
                big_choice[varname]={'type':'boolean'}

            elif choices=='textbox' or choices=='password':
                # Free text with a default value
                var_value=self.parent.preferences.get(varname)
                vars[varname]=StringVar()
                vars[varname].set(var_value)
                lbl = Label(self.settings,text=varname)
                lbl.grid(row=row,column=0)
                if choices == 'password':  s='*'
                else: s=None
                Entry(self.settings,textvariable=vars[varname],
                      bg='white',width=15,show=s).grid(row=row,column=1,columnspan=4)
                big_choice[varname]={'type':'textbox'}

                # Make a dropdown list of previous choices
                try:
                    self.parent.preferences.get(varname+'_previous')
                except:
                    self.parent.preferences.set(varname+'_previous',[default])
                prev_choices=self.parent.preferences.get(varname+'_previous')
                self.mb  =  Menubutton (self.settings,text="->",relief=RAISED )
                self.mb.grid(row=row,column=5)
                self.status_menu=Menu(self.mb,tearoff=0)

                # Print the project names
                for setting in prev_choices:
                    self.status_menu.add_radiobutton(label=setting,
                                            variable=vars[varname],
                                            value=setting,
                                            indicatoron=1)
                self.mb['menu']=self.status_menu
                self.balloon.bind(self.mb,'Previous values')
            row=row+1
            if helptxt!='':
                self.balloon.bind(lbl, helptxt)

        # Functions for saving settings
        def save_settings():
            for varname in big_choice.keys():
                if big_choice[varname]['type']=='options':
                    value=vars[varname].get()
                    self.parent.preferences.set(varname,value)
                elif big_choice[varname]['type']=='boolean':
                    value=vars[varname].get()
                    self.parent.preferences.set(varname,value)
                elif big_choice[varname]['type']=='textbox':
                    value=vars[varname].get()
                    self.parent.preferences.set(varname,value)
                    # Save the previous value
                    prev_vals=self.parent.preferences.get(varname+'_previous')
                    if not value in prev_vals:
                        prev_vals.append(value)
                    self.parent.preferences.set(varname+'_previous',prev_vals)

                else:
                    raise Exception('Unknown preference type')
            self.settings.destroy()
            self.parent.preferences.save_prefs()
            if callback != None:
                callback()
            return

        def cancel():
            self.settings.destroy()
            return

        # Buttons for saving or cancelling
        bf = Frame(self.settings); bf.grid(row=row+1,column=0,
                                           columnspan=4,padx=2,pady=2)
        Button(bf,text='Save settings',command=save_settings).pack(side=RIGHT,fill=BOTH)
        Button(bf,text='Close', command=cancel).pack(side=RIGHT,fill=BOTH)
        return

