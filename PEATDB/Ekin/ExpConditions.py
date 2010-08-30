from Tkinter import *
import pdb
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-2006
# All rights reserved
#

class ExpConditions:

    """Stores and shows dialogs for the experimental conditions in Ekin"""
    def __init__(parent=None):
        self.parent=parent
        self.Meta_Data = []
        self.Operator_Meta = {'operator':'','date':'','status':''}
        self.exp_type = ''
        self.text_frame=None
        self.data_definition()
        self.protocol_textdump=[]
        print "init finished !"
        #use for recording menu radio button choices where needed
        self.choiceVars= {}

    #
    # Commits current exp conditions to data structure
    #
    def commit_exp_cond(self):
        import copy,os

        i =0
        for J in self.Meta_Data:

            if self.Meta_Data[i].has_key('type'):
                if self.Meta_Data[i]['type']=='file':
                    file_name = copy.copy(self.entries[i].get())
                    if os.path.isfile(file_name):
                        try:
                            fd=open(file_name,'r')
                            lines =fd.readlines()
                            self.Meta_Data[i]['value'] = lines
                            fd.close()
                        except:
                            print "could not read file ",file_name
            #if radiomenubutton is type get contents of corresponding stringvars
            #elif self.Meta_Data[i].has_key('choices'):
            #    self.Meta_Data[i]['value'] = copy.copy(self.tv[i].get())
            elif self.Meta_Data[i].has_key('Wavelengths') or self.Meta_Data[i].has_key('choices'):
                self.Meta_Data[i]['value'] = copy.copy(self.wavelengths.get())
            else:
                self.Meta_Data[i]['value'] = copy.copy(self.entries[i].get())
            i = i + 1

        # last 3 entries are from Operator_Meta
        k = len(self.entries) - len(self.Operator_Meta.keys())+1

        self.Operator_Meta['operator']=self.entries[k].get()
        self.Operator_Meta['date']=self.entries[k+1].get()
        self.set_status()

        self.meta_win.destroy()

    #
    #Show dynamic experimental conditions dialog, depends on mode
    #
    def show_exp_dialog(self):

        import time
        #print "showing dialog Meta data ",self.Meta_Data
        self.tv = []
        self.labels = []
        self.entries = []
        self.menus = []
        self.units = []
        textVar = StringVar()
        self.status_type = StringVar()
        self.wavelengths=StringVar()
        self.status = ['Assigned','In progress','Completed']

        row_c = 0
        column = 0

        for I in self.Meta_Data:

            if I['value']!=None and I['value']!='':
                textVar.set(I['value'])
                self.tv.append(textVar)
            else:
                textVar.set(I['default'])
                self.tv.append(textVar)

            #self.tv[row_c].set(textVar.get())
            v = self.tv[row_c].get()
            #print 'v',v

            self.labels.append(Label(self.meta_win,text= I['Name']))
            self.labels[row_c].grid(row=row_c,column=column,sticky='ne',padx=2,pady=2)
            if I['Name']=='Wavelengths' and I.has_key('choices'):
                values=I['choices']

                #print values

                self.wavelengths.set(v)
                m_button=Menubutton(self.meta_win,textvariable=self.wavelengths,relief=RAISED,width=20)
                m_menu=Menu(m_button,tearoff=0)
                m_button['menu']=m_menu

                for c in values:
                    m_menu.add_radiobutton(label=c,variable=self.wavelengths,
                                                    value=c,
                                                    indicatoron=1)
                self.entries.append(m_button)
            else:
                self.entries.append(Entry(self.meta_win,width=20))
                self.entries[row_c].select_clear()
            if I.has_key('choices'):
                #set the radiobutton widget to the right value
                print  'self.tv[row_c]',self.tv[row_c].get()
            else:
                if I['Name']=='Experimental protocol':
                    self.protocol_textdump=I['value']
                    self.entries[row_c].insert(END,v)
                else:
                    self.entries[row_c].insert(0,v)

            self.entries[row_c].grid(row=row_c,column=2,sticky='nw',padx=2,pady=2)

            if I.has_key('units'):
                unit = I['units']
                self.units.append(Label(self.meta_win,text= unit).grid(row=row_c,column=3,sticky='nw',padx=2,pady=2))
            elif I['Name']=='Experimental protocol':
                loadbutton=Button(self.meta_win,text='Load file',command=self.Load_File)
                loadbutton.grid(row=row_c+1,column=0,sticky='news',padx=2,pady=2)
                self.fileShow = self.entries[row_c]
                #show file button
                showF= Button(self.meta_win,text='Display',command=self.Show_File)
                showF.grid(row=row_c+1,column=2,padx=2,pady=2,sticky='news')

            else:
                unit = ''
                self.units.append(Label(self.meta_win,text= unit).grid(row=row_c,column=3,sticky='ne'))

            row_c = row_c + 1

        # operator
        self.labels.append(Label(self.meta_win,text= 'operator'))
        self.labels[row_c].grid(row=row_c+1,column=column,sticky='ne',padx=2,pady=2)

        self.entries.append(Entry(self.meta_win,width=10))
        self.entries[row_c].select_clear()
        if self.Operator_Meta['operator'] == '':
            self.entries[row_c].insert(0,'')
        else:
            self.entries[row_c].insert(0,self.Operator_Meta['operator'])
        self.entries[row_c].grid(row=row_c+1,column=2,sticky='news',padx=2,pady=2)

        # date
        row_c = row_c + 1
        self.labels.append(Label(self.meta_win,text= 'date'))
        self.labels[row_c].grid(row=row_c+1,column=column,sticky='ne',padx=2,pady=2)
        if self.Operator_Meta['date'] =='':
            date = time.asctime()
        else:
            date = self.Operator_Meta['date']
        print "date ",date
        self.entries.append(Entry(self.meta_win,width=20))
        self.entries[row_c].select_clear()
        if self.Operator_Meta['date'] == '':
            self.entries[row_c].insert(0,str(date))
        else:
            self.entries[row_c].insert(0,self.Operator_Meta['date'])
        self.entries[row_c].grid(row=row_c+1,column=2,sticky='ne',padx=2,pady=2)

        # status
        row_c = row_c + 1
        self.labels.append(Label(self.meta_win,text= 'status'))
        self.labels[row_c].grid(row=row_c+1,column=column,sticky='ne')

        if self.Operator_Meta['status'] != '':
            self.status_type.set(self.Operator_Meta['status'])

        self.status_button=Menubutton(self.meta_win,textvariable=self.status_type,relief=RAISED,width=20)
        self.status_menu=Menu(self.status_button,tearoff=0)
        self.status_button['menu']=self.status_menu

        for status in self.status:
            self.status_menu.add_radiobutton(label=status,
                                            variable=self.status_type,
                                            value=status,
                                            indicatoron=1,
                                            command=self.set_status)

        self.status_button.grid(row=row_c+1,column=2)


        OKbutton = Button(self.meta_win,command=self.commit_exp_cond,text='OK')
        OKbutton.grid(row=row_c+2,column=0,columnspan=2,sticky='news',padx=4,pady=4)
        Closebutton = Button(self.meta_win,command=self.commit_exp_cond,text='Close')
        Closebutton.grid(row=row_c+2,column=2,columnspan=1,sticky='news',padx=4,pady=4)
        return

    def update_textvar(self):
		i=0
		for t in self.tv:
			v=t[i].get()
			print v
			i=i+1
		return

    def set_status(self):

        self.Operator_Meta['status']= self.status_type.get()
        print "status ",self.status_type.get()
        return

    def Load_File(self):
        """Load the exp protocol from a file and dump the text in a list to be stored"""
        import tkFileDialog, os
        self.protocol = tkFileDialog.askopenfile(mode='r',filetypes= [("all files", "*")], title="Choose A File")

        if self.protocol == None:
            return
        import fileinput
        self.protocol_textdump=[]
        self.fileShow.delete(0, END)
        for line in fileinput.input(self.protocol.name):
            S=line.rstrip('\n')
            self.protocol_textdump.append(line)
            self.fileShow.insert(END,S+'.')
            #print 'appending',line

        if self.protocol != None:
            for I in self.Meta_Data:
                if I['Name'] == 'Experimental protocol':
                    I['value'] = self.protocol_textdump

        self.meta_win.focus_set()

        return

    def Show_File(self):
        '''Show contents of the protocol data in text window'''
        import textFrame
        self.text_frame=textFrame.textFrame(self)
        try:
            #self.text_frame.load_text(self.fileShow.get())
            self.text_frame.load_text(self.protocol_textdump)
        except:
            import tkMessageBox
            tkMessageBox.showwarning('No file selected',
                                     'Please load a file first',
                                     parent=self.meta_win)

        return

    def set_metadata(self,event=None):
        """Set up exp meta data dialog window"""
        self.fileShow = ''
        self.meta_win=Toplevel()
        self.meta_win.title("Experimental conditions")
        self.xsize = 400
        if self.currentmode == 'Protein Stability' or self.currentmode == 'pH-stability':
            self.ysize = 460
        else:
            self.ysize = 380
        #assign from the Ekin_main variable
        self.exp_type = self.currentmode
        self.meta_win.geometry("%dx%d" %(self.xsize+5,self.ysize+5))
        #print "parent ",self.meta_win.winfo_toplevel()
        if self.Meta_Data == []:
            self.data_definition(self.exp_type)
        self.show_exp_dialog()
        self.protocol = ''

        return

    def get_metadata(self,MetaDir):
        """Gets the meta data list from the project data"""
        keys = MetaDir.keys()

        for I in keys:
            if I == 'Operator':
                self.Operator_Meta = MetaDir[I]
                #print "Operator ",self.Operator_Meta

            elif I == 'Meta_Data':
                self.Meta_Data = MetaDir[I]

        #for I in self.Meta_Data:
        #     print I['Name'],':',I['value']

        return


