#!/usr/bin/env python
import os
#
# Reformat experimentl restraints into a single python pickle file
#
dirnames=['Fergal','Helen','HA'] # Directories that hold the files
res=['no_ghost','questionable','rock_solid']
atoms=['H','N','HA']

exp_data={}
for dirname in dirnames:
    for restraint in res:
        for atom in atoms:
            if dirname=='HA' and atom!=dirname:
                continue
            if dirname!='HA' and atom=='HA':
                continue
            #
            # Loop over all atoms that we are looking at
            #
            filename=os.path.join(os.getcwd(),dirname,'%s_%s.csv' %(restraint,atom))
            if not os.path.isfile(filename):
                print 'Could not find: %s' %filename
                continue
            #
            # read the file and store the ghosts
            #
            print 'Reading',filename
            fd=open(filename)
            lines=fd.readlines()
            fd.close()
            #
            header=lines[0].split(',')
            for line in lines[1:]:
                if line.strip()=='*':
                    continue
                sp=line.split(',')
                #print sp
                import string
                nucleus=':'+string.zfill(sp[0],4)
                count=0
                for titgroup in header[1:]:
                    titgroup=titgroup.strip().upper().replace('"','')
                    T={'D':'ASP','H':'HIS','E':'GLU','CTERM':':0129:CTERM'}
                    if T.has_key(titgroup):
                        titgroup=T[titgroup]
                    else:
                        #print T.keys(),titgroup[0]
                        if T.has_key(titgroup[0]):
                            titgroup=':'+string.zfill(titgroup[1:],4)+':'+T[titgroup[0]]
                    count=count+1
                    #print titgroup,header[count],count
                    #
                    # Make sure we hae the right keys in the dictionary
                    #
                    if not exp_data.has_key(titgroup):
                        exp_data[titgroup]={}
                    if not exp_data[titgroup].has_key(nucleus):
                        exp_data[titgroup][nucleus]={}
                    #
                    # Deal with the individual files
                    #
                    if restraint=='no_ghost':
                        #
                        # if we have a zero in the no_ghost file, then it means that there is no ghost titration 
                        #
                        if sp[count].strip()=='0':
                            if exp_data[titgroup][nucleus].has_key('atom'):
                                print 'Duplicate restraint %s for titgroup %s, nucleus: %s, atom %s' %(restraint,titgroup,nucleus,atom)
                                print exp_data[titgroup][nucleus][atom]
                                raise Exception()
                            exp_data[titgroup][nucleus][atom]='absent'
                    elif restraint=='questionable':
                        if sp[count].strip()!='':
                            if exp_data[titgroup][nucleus].has_key('atom'):
                                print 'Duplicate restraint for titgroup %s, nucleus: %s, atom %s' %(restraint,titgroup,nucleus,atom)
                                print exp_data[titgroup][nucleus][atom]
                                raise Exception()
                            exp_data[titgroup][nucleus][atom]=('q'+sp[count].strip()).replace('"','')
                    elif restraint=='rock_solid':
                        if sp[count].strip()!='':
                            if exp_data[titgroup][nucleus].has_key('atom'):
                                print 'Duplicate restraint %s for titgroup %s, nucleus: %s, atom %s' %(restraint,titgroup,nucleus,atom)
                                print exp_data[titgroup][nucleus][atom]
                                raise Exception()
                            exp_data[titgroup][nucleus][atom]=sp[count].strip().replace('"','')

#
# Write the pickle file
#                    
import pickle
fd=open('new_restraints.pickle','w')
pickle.dump(exp_data,fd)
fd.close()

#
# Count the number of restraints
#
print exp_data.keys()
count=0
countab=0
atom_count={'H':0,'N':0,'HA':0}
for tg in exp_data.keys():
    for nucleus in exp_data[tg].keys():
        for atom in exp_data[tg][nucleus].keys():
            if exp_data[tg][nucleus][atom].lower()!='absent':
                count=count+1
                atom_count[atom]=atom_count[atom]+1
            else:
                countab=countab+1
print count,'positive restraints'
print countab,'absent restraints'
print count+countab
print atom_count

        