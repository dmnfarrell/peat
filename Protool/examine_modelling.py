#!/usr/bin/env python
import os
dirname=os.path.join(os.getcwd(),'scores')
files=os.listdir(dirname)
data=[]
x=[]
y=[]
for filename in files:
    fn=os.path.join(dirname,filename)
    if fn[-4:]!='_old':
        fd=open(fn)
        import pickle
        score=pickle.load(fd)
        fd.close()
        if not os.path.isfile(fn+'_old'):
            continue
        fd=open(fn+'_old')
        oss=pickle.load(fd)
        fd.close()
        if oss<10.0:
            x.append(oss)
            y.append(score)
        if score>1.0 and oss<20.0:
            data.append([score,oss,filename])
data.sort()
for dp in data:
    print dp

import pylab
pylab.plot(x,y,'bo')
pylab.show()
        
