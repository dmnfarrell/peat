#!/usr/bin/python
#
# Pymol classes for structural mapping
# UCD @2010

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import sys, time, os
import pickle
import pymol
__main__.pymol = pymol
import pymol.cmd as cmd
from pymol import stored

pymol.finish_launching()

   
class structMapper(object):
    """Map features on structures using pymol"""
    colors = ['bluewhite', 'red',
        'br0','orange','lightteal','yellow',
        'purple','darksalmon','deepolive',
        'limegreen','pink']
    
    titatoms = {'ASP':'CG', 'GLU':'CG', 
                'HIS':'CE1', 'CYS':'SG','TYR':'CZ', 
                'LYS':'NZ', 'ARG':'CZ'} 
                
    def __init__(self, **kwargs):  
     
        return

    def settings(self):
        cmd.hide('lines')
        cmd.show('cartoon')
        cmd.hide('everything','resn hoh')
        cmd.color('green','all')
        cmd.bg_color('white')
        cmd.set('ray_trace_fog',0)
        cmd.set('stick_transparency',0.2)
        cmd.set('transparency',0.8)
        cmd.set('sphere_scale', 0.5)
        cmd.set('stick_radius', 0.5)
        cmd.set('label_position', (1,2,1))
        cmd.set('label_color', 'black')
        cmd.set('label_size', 10)
        cmd.set('label_font_id', 7)
        cmd.turn('z',90)        
        cmd.turn('x',60)
        return
        
    def mapMutations(self):
        """map frequency of mutations for 10R"""
        from Correlation import CorrelationAnalyser
        corr = CorrelationAnalyser()
        from PEATDB.Base import PDatabase
        corr.DB = PDatabase(server='localhost', username='farrell',
                             password='123', project='novo')
        freq = corr.analyseMutations(sheet='stability7', filterby=('fit',['ok','?']),
                                xcol='Total')
        print freq
        mx=float(max(freq.values()))
        self.settings()
        for f in freq:
            resnum = f
            val=freq[f]
            x=[1.0,val/mx,0.1]
            print val,x
            cmd.set_color(str(val)+'c',x)
            cmd.select('r','resi %s' %(resnum))
            cmd.color(str(val)+'c','r')
            #cmd.show('sticks','r')
            cmd.label('r and name ca', str(val))
          
        cmd.select('asite', 'resi 35, resi 64, resi 143')        
        cmd.color('blue','asite')
        cmd.show('sticks','asite')
        cmd.label('asite and name n', 'resn+resi')
        
        '''from SpectrumBar import spectrumbar
        spectrumbar([1.0,1.0,1.0],radius=1.0,name=spectrumbar,
            head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0)'''  
        self.save()
        return
        
    def load(self,pdbfile):
        self.pdbfile=pdbfile
        cmd.load(pdbfile, 'pdb')
        return
         
    def save(self):
        name = os.path.splitext(self.pdbfile)[0]
        cmd.set('ray_trace_mode',0)
        cmd.png(name+'.png',800,800)
        cmd.save(name+'.pse')        
        return


def main():
    from optparse import OptionParser
    parser = OptionParser()
    S = structMapper()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a pdb")
    opts, remainder = parser.parse_args()
    S.load(opts.file)
    S.mapMutations()
    #S.save()
    cmd.quit()
    return
    
if __name__ == '__main__':
    main()

