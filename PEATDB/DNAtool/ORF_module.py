#
# $Id: ORF_module.py 272 2005-11-10 02:04:11Z nielsen $
#
# ORF_module - routines for selecting and manipulating an ORF
#
# (c) Copyright 2005 Jens Erik Nielsen, University College Dublin
#

from Tkinter import *

class ORF_handler(Frame):

    def __init__(self,parent_class):
        #
        # Set data
        #
        self.parent_class=parent_class
        self.aa_seq=None
        #
        # Windows
        #
        self.ORF_win=None
        #
        # Select ORF
        #
        self.select_ORF()
        self.parent_class.master.wait_window(self.ORF_win) # Wait for the ORF window to close
        return

    #
    # ----------
    #

    def find_ORF(self):
        #
        # Find an ORF
        #
        self.parent_class.data['ORFS']={}
        #
        # Loop over all aa data
        #
        frame_num=0
        for frame in self.parent_class.data['AA_seqs1']:
            frame_num=frame_num+1
            #
            # Should we look at all frames?
            #
            frame_sel=self.frame.get()
            if frame_num!=frame_sel and frame_sel!=0:
                continue
            #
            # yes, this one is ok
            #
            if not self.parent_class.data['ORFS'].has_key(frame_num):
                self.parent_class.data['ORFS'][frame_num]={}
            #
            # Find the ORFs for this frame
            #
            ORF_num=0
            position=1
            ORF_start=1
            ORF_stop=None
            in_ORF=None
            stop_only_on_stopcodon=self.stop_condition.get()
            for AA in frame:
                #
                # Are we within the boundaries set by the Scales?
                #
                scale_start=self.ORF_start.get()
                scale_stop=self.ORF_stop.get()
                #print 'fr_num: %2d, pos: %3d, start: %3d, stop: %3d' %(frame_num,position,scale_start,scale_stop)
                if position>=scale_stop or position<scale_start:
                    #
                    # We are outside the ranges set by the scales
                    #
                    # Update the position
                    #
                    position=position+1
                    if position>scale_stop:
                        #
                        # If we are stopping early then make sure to add the ORF we were in
                        #
                        if in_ORF:
                            #
                            # Can we stop on any AA?
                            #
                            stop_possible=1
                            if stop_only_on_stopcodon==1:
                                if AA!='*':
                                    stop_possible=None
                            #
                            # Add the ORF if allowed
                            #
                            if stop_possible:
                                ORF_stop=position-1
                                ORF_num=ORF_num+1
                                self.parent_class.data['ORFS'][frame_num][ORF_num]={'start':ORF_start,'stop':ORF_stop,'length':ORF_stop-ORF_start+1}
                                in_ORF=None
                        #
                        # If we're too far along then we break
                        #
                        break
                    else:
                        #
                        # We haven't reached the start yet so we continue
                        #
                        # Skip to next AA
                        #
                        continue
                #
                # Could this be a start position?
                #
                possible_start=None
                #
                # self.start_contition==1 means we require Met to start an ORF
                #
                if self.start_condition.get()==1:
                    if AA=='M':
                        possible_start=1
                else:
                    if AA!='*':
                        possible_start=1
                #
                # Is this a stop position?
                #
                must_stop=None
                if AA=='*':
                    must_stop=1
                #
                # Should we extend the ORF, start or stop?
                #
                if in_ORF:
                    #
                    # If we're in an ORF
                    #
                    if must_stop:
                        #
                        # Here we only stop if on a stop codon
                        #
                        ORF_stop=position
                        ORF_num=ORF_num+1
                        self.parent_class.data['ORFS'][frame_num][ORF_num]={'start':ORF_start,'stop':ORF_stop,'length':ORF_stop-ORF_start+1}
                        in_ORF=None
                        if self.display_ORF.get()==0:
                            break
                    else:
                        pass
                else:
                    if possible_start:
                        ORF_start=position
                        in_ORF=1
                    else:
                        pass
                #
                # Update counter 
                #
                position=position+1
            #
            # Make sure we capture the last ORF
            #
            if in_ORF:
                #
                # If we're in an ORF
                #
                #
                # Can we stop on any AA?
                #
                stop_possible=1
                if stop_only_on_stopcodon==1:
                    if AA!='*':
                        stop_possible=None
                #
                # Add the ORF if allowed
                #
                if stop_possible:
                    ORF_stop=position
                    ORF_num=ORF_num+1
                    print 'Adding final ORF',{'start':ORF_start,'stop':ORF_stop,'length':ORF_stop-ORF_start+1}
                    self.parent_class.data['ORFS'][frame_num][ORF_num]={'start':ORF_start,
                                                                        'stop':ORF_stop,
                                                                        'length':ORF_stop-ORF_start+1}
                    in_ORF=None
            #
            # Next frame
            #
        #
        # Find the longest ORF
        #
        max_ORF=0
        frame_sel=None
        ORF_sel=None
        for frame in self.parent_class.data['ORFS'].keys():
            #print 'Checking frame %3d' %(frame)
            for ORF_num in self.parent_class.data['ORFS'][frame].keys():
                if self.parent_class.data['ORFS'][frame][ORF_num]['length']>max_ORF:
                    max_ORF=self.parent_class.data['ORFS'][frame][ORF_num]['length']
                    frame_sel=frame
                    ORF_sel=ORF_num
        #
        # Store the data in the data array
        #
        # parent_class is the DNAtool class..
        #
        self.parent_class.data['ORF_selected']={'frame':frame_sel,'number':ORF_sel}
        self.parent_class.data['ORF_selected']['start']=self.parent_class.data['ORFS'][frame_sel][ORF_sel]['start']
        start=self.parent_class.data['ORF_selected']['start']
        self.parent_class.data['ORF_selected']['stop']=self.parent_class.data['ORFS'][frame_sel][ORF_sel]['stop']
        stop=self.parent_class.data['ORF_selected']['stop']
        self.parent_class.data['ORF_selected']['length']=self.parent_class.data['ORFS'][frame_sel][ORF_sel]['length']
        self.parent_class.data['ORF_selected']['aaseq']=self.parent_class.data['AA_seqs1'][frame_sel-1][start-1:stop]
        self.parent_class.data['ORF_selected']['aaseq3']=self.parent_class.data['AA_seqs3'][frame_sel-1][start-1:stop]
        #self.parent_class.data['frame_sel']=frame_sel
        #self.parent_class.data['ORF_sel']=ORF_sel
        return

    
    #
    # -------------
    #

    def print_cur_ORF(self):
        #
        # Print the current ORF
        #
        #ORF_num=self.parent_class.data['ORF_sel']
        #frame=self.parent_class.data['frame_sel']
        #found_ORF=1
        #if not self.parent_class.data['ORFS'].has_key(frame):
        #    found_ORF=None
        #elif len(self.parent_class.data['ORFS'][frame].keys())==0:
        #    found_ORF=None
        #
        # we found an ORF so we can set the parameters
        #
 ##        #if found_ORF:
##             #start=self.parent_class.data['ORFS'][frame][ORF_num]['start']
##             #stop=self.parent_class.data['ORFS'][frame][ORF_num]['stop']
##             #length=self.parent_class.data['ORFS'][frame][ORF_num]['length']
##             self.aa_seq=self.parent_class.data['AA_seqs1'][frame-1][start-1:stop]
##             #self.DNA_seq=self.parent_class.data['AA_seqs1'][frame-1][start-1:stop]
##             #self.parent_class.data['ORF_selected']={'start':start,'stop':stop,'frame':frame,'ORF_num':ORF_num}
##             #frame_s=str(frame)
##         else:
##             #
##             # No ORF here...
##             #
##             start=0
##             stop=0
##             length=0
##             frame_s='-'
##             self.aa_seq='No ORF found with current parameters'
##             self.parent_class.data['ORF_selected']=None
        #
        # Set the paremeters for printing
        #
        data=self.parent_class.data['ORF_selected']
        frame_s=str(data['frame'])
        #
        row=3
        l1=Label(self.ORF_win,text='Frame %s' %str(data['frame']))
        l1.grid(row=row,column=0)
        #
        l4=Label(self.ORF_win,text='Length: %4d' %(data['length']))
        l4.grid(row=row,column=1)
        #
        row=row+1
        l2=Label(self.ORF_win,text='Start: %4d' %data['start'])
        l2.grid(row=row,column=0)
        #
        row=row+0
        l3=Label(self.ORF_win,text='END: %4d' %data['stop'])
        l3.grid(row=row,column=1)
        #
        # Set the amino acid sequence
        #
        self.aaseq=self.parent_class.data['ORF_selected']['aaseq']
        #
        #
        #
        # Print the AA sequence in a textbox
        #
        scrollbar=Scrollbar(self.ORF_win,orient='vertical')
        scrollbar.grid(row=0,column=3,rowspan=6,sticky='NEWS')
        #
        # Textbox
        #
        self.textbox=Text(self.ORF_win,background='white',foreground='black',width=50,state=NORMAL,exportselection=1,yscrollcommand=scrollbar.set)
        self.textbox.grid(row=0,column=2,rowspan=6,sticky='NEWS')
        self.textbox.config(state=NORMAL)
        self.textbox.insert(END,self.aaseq)
        #
        # Set the scale
        #
        if data['length']>0:
            self.ORF_stop.set(data['stop'])
            self.ORF_start.set(data['start'])
        return

    #
    # -------------
    #

    def select_ORF(self):
        #
        # Open window
        #
        self.ORF_win=Toplevel()
        self.ORF_win.title('Select Open Reading Frame')
        self.ORF_win.geometry('+100+400')
        #
        # Pulldown menu
        #
        menu=Menu(self.ORF_win)
        #
        # File menu
        #
        ORF_menu=Menu(menu,tearoff=0)
        self.frame=IntVar()
        ORF_menu.add_radiobutton(label='Any',variable=self.frame,value=0)
        ORF_menu.add_radiobutton(label='Frame 1',variable=self.frame,value=1,command=self.update_ORF)
        ORF_menu.add_radiobutton(label='Frame 2',variable=self.frame,value=2,command=self.update_ORF)
        ORF_menu.add_radiobutton(label='Frame 3',variable=self.frame,value=3,command=self.update_ORF)
        menu.add_cascade(label='Frame',menu=ORF_menu)
        # ---------------------------------------------------------
        #
        # ORF start menu
        #
        start_menu=Menu(menu,tearoff=0)
        self.start_condition=IntVar()
        start_menu.add_radiobutton(label='Met',variable=self.start_condition,value=1,command=self.update_ORF_reset_count)
        start_menu.add_radiobutton(label='Any',variable=self.start_condition,value=0,command=self.update_ORF_reset_count)
        #
        # Add the start_menu to the main menu
        #
        menu.add_cascade(label='ORF start',menu=start_menu)
        # ---------------------------------------------------------
        #
        # ORF stop menu
        #
        stop_menu=Menu(menu,tearoff=0)
        self.stop_condition=IntVar()
        stop_menu.add_radiobutton(label='STOP Codon',variable=self.stop_condition,value=1,command=self.update_ORF_reset_count)
        stop_menu.add_radiobutton(label='Any',variable=self.stop_condition,value=0,command=self.update_ORF_reset_count)
        #
        # Add the start_menu to the main menu
        #
        menu.add_cascade(label='ORF stop',menu=stop_menu)
        # ---------------------------------------------------------
        #
        # ORF display menu
        #
        stop_menu=Menu(menu,tearoff=0)
        self.display_ORF=IntVar()
        stop_menu.add_radiobutton(label='First ORF',variable=self.display_ORF,value=0,command=self.update_ORF)
        stop_menu.add_radiobutton(label='Longest ORF',variable=self.display_ORF,value=1,command=self.update_ORF)
        #
        # Add the start_menu to the main menu
        #
        menu.add_cascade(label='Display',menu=stop_menu)
        # ---------------------------------------------------------
        self.ORF_win.config(menu=menu)
        #
        # Scales for selecting start and stop
        #
        self.ORF_start=IntVar()
        self.ORF_stop=IntVar()
        row=0
        l=Label(self.ORF_win,text='ORF start:',relief='ridge')
        l.grid(row=row,column=0)
        frame=self.frame.get()
        scl=Scale(self.ORF_win,from_=1,to=len(self.parent_class.data['AA_seqs1'][frame]),resolution=1,
                  orient='horizontal',relief='ridge',variable=self.ORF_start,label='frame number')
        scl.grid(row=row,column=1, sticky='wens')
        #
        #
        row=1
        l2=Label(self.ORF_win,text='ORF end:',relief='ridge')
        l2.grid(row=row,column=0)
        scl2=Scale(self.ORF_win,from_=1,to=len(self.parent_class.data['AA_seqs1'][1]),resolution=1,
                  orient='horizontal',relief='ridge',variable=self.ORF_stop,label='frame number')
        scl2.grid(row=row,column=1, sticky='wens')
        #
        # Button for Updating selection
        #
        btn=Button(self.ORF_win,text='Update',command=self.update_ORF)
        btn.grid(row=2,column=0)
        btn2=Button(self.ORF_win,text='I found my ORF!',command=self.pick_ORF,fg='green')
        btn2.grid(row=2,column=1)
        #
        # Set defaults
        #
        self.ORF_stop.set(len(self.parent_class.data['AA_seqs1'][1]))
        self.ORF_start.set(1)
        self.start_condition.set(0)
        self.stop_condition.set(0)
        self.display_ORF.set(1)
        self.frame.set(0)
        #
        # Update
        #
        self.update_ORF()
        #
        # Done
        #
        return

    #
    # -----------------
    #

    def update_ORF_reset_count(self,junk=None):
        #
        # Reset the values for the Scales so we search the entire sequence again
        #
        frame=self.frame.get()
        self.ORF_start.set(1)
        self.ORF_stop.set(len(self.parent_class.data['AA_seqs1'][frame]))
        self.update_ORF()
        return

    #
    # -----------------
    #

    def update_ORF(self,junk=None):
        #
        # Make sure that the ranges are ok
        #
        if self.ORF_stop.get()<self.ORF_start.get():
            self.ORF_stop.set(self.ORF_start.get())
        #
        # Update the ORF selection
        #
        self.find_ORF()
        self.print_cur_ORF()
        return

    #
    # -----------------
    #

    def pick_ORF(self,junk=None):
        #
        # ORF picked - Get number of first AA
        #
        first_aanum=None
        import tkSimpleDialog
        first_aanum=tkSimpleDialog.askinteger('Number AA seq',
                                              'Give the number you want to use for the first amino acid\n (negative values ok)',
                                              parent=self.ORF_win,initialvalue=1)
        if not first_aanum:
            if self.parent_class.data.has_key('ORF_selected'):
                del self.parent_class.data['ORF_selected']
        else:
            if not self.parent_class.data.has_key('ORF_selected'):
                self.parent_class['ORF_selected']={}
            self.parent_class.data['ORF_selected']['aastart_number']=first_aanum
        #
        # Return
        #
        self.ORF_win.destroy()
        return 
#
# ----------
#

if __name__=="__main__":
    print
    print 'Do not run this module on its own'
    print
