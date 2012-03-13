#!/usr/bin/env python

import pKaTool.pKaIO, sys

IO=pKaTool.pKaIO.pKaIO(sys.argv[1])
matrix=IO.read_matrix()

exp={':0035:GLU':{':0052:ASP':1.2,':0119:ASP':0.3,':0007:GLU':0.3}}
for T1 in exp.keys():
    for T2 in exp[T1].keys():
        if not exp.has_key(T2):
            exp[T2]={}
        exp[T2][T1]=exp[T1][T2]


matrix2={}
tTGs=matrix.keys()
TGs=[]
for tg in tTGs:
    if tg.find('TYR')==-1:
        TGs.append(tg)
for TG1 in TGs:
    matrix2[TG1]={}
    for TG2 in TGs:
        value=0.0
        if exp.has_key(TG1):
            if exp[TG1].has_key(TG2):
                import math
                value=math.log(10)*exp[TG1][TG2]
        #value=abs((matrix[TG1][TG2][0]+matrix[TG2][TG1][0])/2.0)
        #if value<0.7:
        #    value=value
        matrix2[TG1][TG2]=value

import TwoDplots

TwoDplots.heatmap(matrix2,title='Interaction energy matrix',firstkey='TG1',secondkey='TG2',zlabel='Interaction energy (kT)',firstticks=sorted(TGs),secondticks=sorted(TGs))