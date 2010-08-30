#!/usr/bin/python
#
# Test pymol API
# UCD @2010

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pickle
import pymol
import pymol.cmd as cmd
from pymol import stored

pymol.finish_launching()
pdbfile = os.path.abspath(sys.argv[1])
ghostsfile = os.path.abspath(sys.argv[2])
try:
    ghostsfile1 = os.path.abspath(sys.argv[3])
except:
    ghostsfile1 = None
name = os.path.splitext(ghostsfile)[0] 

colors = ['bluewhite', 'red',
 	'br0','orange','lightteal','yellow',
 	'purple','darksalmon','deepolive',
 	'limegreen','pink']

titatoms = {'ASP':'CG', 'GLU':'CG', 
            'HIS':'CE1', 'CYS':'SG','TYR':'CZ', 
            'LYS':'NZ', 'ARG':'CZ'}

pkatoms = {'pKa':'n','pKa1':'c','pKa2':'ca','pKa3':'cb'}
        
sourceclrs = {}
cmd.load(pdbfile, 'pdb')
#cmd.hide('lines')
#cmd.show('sticks')
#cmd.show('ribbon')
#cmd.show('surface')
cmd.hide('everything','resn hoh')
cmd.color('green','all')
cmd.bg_color('white')
cmd.set('ray_trace_fog',0)
cmd.set('stick_transparency',0.2)
#cmd.set('transparency',0.8)
cmd.set('sphere_scale', 0.4)
cmd.set('stick_radius', 0.3)
cmd.set('label_position', (1,2,1))
cmd.set('label_color', 'black')
cmd.set('label_size', 10)
cmd.set('label_font_id', 7)

cmd.turn('y',90)

#load ghosts and display
G = pickle.load(open(ghostsfile,'r'))
#create color index for source residues
i=0
for r in G:
    for p in G[r].keys():
        if 'pK' in p:  
            s = G[r][p]           
            if not s in sourceclrs:
                sourceclrs[s] = colors[i]
                i+=1
            if i>=len(colors): i=0
    
def draw(G):    
    for r in G:
        for p in G[r].keys():
            if 'pK' in p:        
                resnum = G[r]['resnum']
                source = G[r][p]
                ch,tresnum,tres = source.split(':')
                tresnum=int(tresnum)
                srcatom=titatoms[tres]
                clr = sourceclrs[source]
                print resnum,tresnum,clr
                labelatom=pkatoms[p]
                cmd.select('bb','resi %s and name %s' %(resnum,labelatom))
                cmd.show('spheres','bb')
                cmd.color(clr,'bb')
                atom = titatoms[tres] 
                cmd.select('tit','resi %s' %(tresnum))
                cmd.color(clr,'tit')
                cmd.show('sticks','tit')
                cmd.label('tit and name %s' %atom,'resn+resi')
            
draw(G)
if ghostsfile1 != None:
    G = pickle.load(open(ghostsfile1,'r'))
    draw(G)
    
cmd.set('ray_trace_mode',0)
cmd.png(name+'.png',800,800)
cmd.save(name+'.pse')
cmd.quit()

 
def function1(userSelection):
    return
cmd.extend( "yourFunction", function1 );


