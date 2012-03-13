#!/usr/bin/env python

N_Ghosts={':0035:GLU':{29:-0.28, 30:-0.16, 32:-0.51, 33:-0.53, 34:0.46, 38:0.42, 42:0.2, 44:0.25, 45:-0.1, 46:0.23, 54:-0.15, 55:-0.13, 56:0.24, 58:-0.7, 91:0.2, 107:0.1, 108:0.76, 113:0.56, 114:0.31, 115:0.27},
    ':0052:ASP':{33:-0.21, 35:0.45, 36:-0.2, 38:0.11, 48:0.32, 51:0.83, 54:-0.28, 56:0.28,71:-0.13, 109:0.1},
    ':0015:HIS':{10:-0.4, 11:0.44, 16:0.09, 18:-0.2, 19:0.22, 87:-0.19, 94:-0.19},
    ':0119:ASP':{77:-0.15, 97:0.1, 99:-0.42, 100:0.71, 102:-0.7, 104:-0.24, 107:0.42}}
    
N_error={':0035:GLU':{29:0.03, 30:0.02,32:0.03,33:0.03,34:0.03,38:0.03,42:0.02,44:0.02,45:0.02,46:0.05,51:0.03,52:0.06,54:0.03,55:0.03,56:0.04,58:0.08,59:0.04,60:0.04,91:0.04,107:0.02,108:0.04,109:0.03,113:0.02,114:0.02,115:0.02},
':0052:ASP':{33:0.03, 35:0.04, 36:0.04, 38:0.02, 48:0.05, 51:0.05, 54:0.02, 56:0.02, 71:0.01, 109:0.01}}
    
H_Ghosts={}
H_error={}
#
# Reformat into restraint file for epsmap determination
#
residues=[]
titgroups=[]
for DC in [N_Ghosts,H_Ghosts]:
    for tg in DC.keys():
        if not tg in titgroups:
            titgroups.append(tg)
        for res in DC[tg].keys():
            if not res in residues:
                residues.append(res)
#
# ------
#            
            
print 'Titratable groups',titgroups
print 'Residues with ghost information',sorted(residues)

import string
ghosts=[]
errors=[]
print '# Titgroup Residue N_Ghost  N_error  H_Ghost  H_Error'
for tg in titgroups:
    for residue in sorted(residues):
        line='%8s %8s ' %(tg,':'+string.zfill(residue,4))
        found=False
        #print tg, residue
        for atom,DC,EDC in [['N',N_Ghosts,N_error],['H',H_Ghosts,H_error]]:
            for name,DDD in [['ghost',DC],['error',EDC]]:
                if DDD.has_key(tg):
                    if DDD[tg].has_key(residue):
                        line=line+'%7.4f  ' %float(DDD[tg][residue])
                        found=True
                        if name=='ghost':
                            ghosts.append(abs(float(DDD[tg][residue])))
                        else:
                            errors.append(float(DDD[tg][residue]))
                    else:
                        line=line+'%7s  ' %('None')
                else:
                    line=line+'%7s  ' %('None')
        if found:
            print line

print 'Number of ghosts: %d, average span: %5.3f' %(len(ghosts),sum(ghosts)/len(ghosts))