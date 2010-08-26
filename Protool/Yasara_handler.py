#
#  Yasara.py
#  Python_Xcode
#
#  Created by Jens on 31/05/2010.
#  Copyright (c) 2010 University College Dublin. All rights reserved.
#
class yasara_handler:

    def __init__(self,yasaradir='/local/nielsen/bin/yasara/',console=True):
        #
        # Default on amylase, mac below
        #
        import os
        if not os.path.isdir(yasaradir):
            yasaradir='/Users/nielsen/desktop/yasara.app'
        dirname=os.path.split(yasaradir)[1]
        if dirname.lower()=='yasara.app':
            yasaradir=os.path.join(yasaradir,'yasara')
        #
        import sys,os
        sys.path.append(os.path.join(yasaradir,'pym'))
        sys.path.append(os.path.join(yasaradir,'plg'))
        import yasaramodule as yasara
        self.yasara=yasara
        if not console:
            self.yasara.info.mode='txt'
            self.yasara.info.licenseshown=0
            self.yasara.Console('Off')
        return

    def load_mol(self,pdbfile,center=None):
        #self.yasara.run('DelAll')
        obj=self.yasara.LoadPDB(pdbfile,center=center)
        self.yasara.run('Style Stick')
        self.yasara.run('HideRes Hoh')
        res=self.yasara.ColorObj(obj,'Grey')
        self.yasara.run('HUD Off') 
        return obj

    def AlignObj(self,obj1,obj2):
        """Align all objects with MOTIF. Return rmsd"""
        resultlist=self.yasara.AlignObj(obj1,obj2,method='MOTIF')
        return resultlist[0]
        
    def buildobj(self,lines):
        """Build a yasara object from a list of pdblines"""
        return self.yasara.BuildPDB(lines)
        
    def col_res(self,obj,residues,color):
        """Color a number of residues a color"""
        if type(residues) is type([]):
            for res in residues:
                self.yasara.ColorRes("%d Obj %d" %(res,obj),color)
        else:
            for res in residues.keys():
                #print "%d Obj %d" %(res,obj),int(180*residues[res])
                self.yasara.ColorRes("%d Obj %d" %(res,obj),int(120*residues[res]))
        return
