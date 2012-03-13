#!/usr/bin/env python

def main():
    import sys, csv
    csvfile=sys.argv[1]
    f=open(csvfile,'r')
    cr = csv.reader(f)
    tmp_names = cr.next()
    names=[]
    for name in tmp_names:
        names.append(name.lower().strip())
    data={}
    #
    for n in names:
        data[n.lower()]=[]
    for row in cr:
        for i in range(len(row)):
            x=row[i]
            n=names[i]
            data[n].append(x)
    #
    # Reformat the data
    #
    #
    cpmg={}
    for index in range(100):
        for ph in ['54','65','77','85']:
            if not cpmg.has_key(ph):
                cpmg[ph]={}
            thisrecord=[]
            for tp in ['res','chg','sdev']:
                dname='%s%s' %(tp,ph)
                if len(data[dname])>index:
                    thisrecord.append(data[dname][index])
            if len(thisrecord)>2:
                res=thisrecord[0]
                if res not in ['','pA','kEX']:
                    cpmg[ph][res]=thisrecord[1:][:]
    #
    # Now find maximal changes between pH values
    #
    resdata={}
    for pH in cpmg.keys():
        for res in cpmg[pH].keys():
            ires=int(res)
            if not resdata.has_key(ires):
                resdata[ires]=[]
            resdata[ires].append(float(cpmg[pH][res][0]))
    #
    print '# Residue, max change due to CPMG'
    residues=sorted(resdata.keys())
    for residue in residues:
        data=resdata[residue]
        print '%3d, %7.3f' %(residue,0.05*abs(max(data)-min(data)))
    
    #print cpmg
            
                
    return
    
    
if __name__=='__main__':
    main()
        
    
    
    