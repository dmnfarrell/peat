#!/usr/bin/env python

class muts:

    def __init__(self,mutfile):
        """Decompose mutation effects"""
        data=self.read_muts(mutfile)
        print data.keys()
        backup=data.copy()
        #
        # find single muts
        #
        unknown=len(data.keys())
        last_unknown=10000
        results={}
        round=1
        while unknown<last_unknown:
            print unknown
            decomp=self.single_muts(data)
            for mutation,dtm in decomp:
                if not results.has_key(mutation):
                    results[mutation]=[]
                results[mutation].append([dtm,round])
            data=self.substitute(data,decomp)
            last_unknown=unknown
            unknown=len(data.keys())
            #
            round=round+1
        #
        # Print the results
        #
        names=[]
        xs=[]
        for mutation in sorted(results.keys()):
            dtms=[]
            for dtm,round in results[mutation]:
                dtms.append(dtm)
            dtms.sort()
            maxdiff=abs(dtms[0]-dtms[-1])
            avg=sum(dtms)/float(len(dtms))
            txt='%10s: avg: %7.2f, maxdiff; %7.2f, ' %(mutation,avg,maxdiff)
            names.append(mutation)
            xs.append(avg)
            for dtm,round in results[mutation]:
                txt=txt+'%7.2f (%2d), ' %(dtm,round)
            print txt
        #
        # Print the clones we could not reduce
        #
        print 'Mutant proteins that cannot be further decomposed'
        for mutant in sorted(data.keys()):
            print '%10s, %7.2f %25s org. mut: %40s' %(mutant,data[mutant]['dtm'],data[mutant]['longform'],backup[mutant]['longform'])
        #
        # Plot weights against what we found here
        #
        weights=self.read_weights('/Volumes/10R_data/weights.csv')
        ws={}
        count=0
        for residue in weights['res']:
            ws[residue[1:]]=float(weights['weight'][count])
            count=count+1
        print ws
        ys=[]
        for name in names:
            if ws.has_key(name):
                ys.append(ws[name])
            else:
                ys.append(0.0)
        import pylab
        pylab.plot(xs,ys,'ro',label='data')
        #
        # Do linear fit
        #
        import numpy
        x=numpy.array(xs)
        y=numpy.array(ys)
        A=numpy.vstack([x,numpy.ones(len(x))]).T
        a,b=numpy.linalg.lstsq(A,y)[0]
        r=numpy.corrcoef(x,y)[0][1]
        # y= ax+b
        xs.sort()
        pylab.plot([xs[0],xs[-1]],[a*xs[0]+b,a*xs[-1]+b],'g-',label='Fit. R=%5.2f' %r)
        pylab.xlabel('Deduced average dT50 (K)')
        pylab.ylabel('Weight from LARS analysis')
        pylab.legend(loc=4)
        pylab.show()
        return
        
    #
    # -----
    #
        
    def read_weights(self,weightfile):
        import csv
        f=open(weightfile,'r')
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
                n=names[i]
                data[n].append(row[i])
        return data
        
    #
    # -----
    #
        
    def substitute(self,data,decomp):
        """Substitute the new data in"""
        for mutation,dtm in decomp:
            for protein in data.keys():
                #print data[protein]['longform']
                muts=data[protein]['longform'].split('+')
                if mutation in muts:
                    muts.remove(mutation)
                    data[protein]['dtm']=data[protein]['dtm']-dtm
                import string
                data[protein]['longform']=string.join(muts,'+')
                if len(muts)==0:
                    del data[protein]
        return data
        
    #
    # ------
    #
    
    def single_muts(self,data):
        """Search for mutant proteins that just contain a single mutation"""
        enes=[]
        for protein in data.keys():
            #print data[protein]['longform']
            num=len(data[protein]['longform'].split('+'))
            if num==1:
                enes.append([data[protein]['longform'],data[protein]['dtm']])
        return enes
        
    #
    # ------
    #

    def read_muts(self,csvfile):    
        """Load the experimental data"""
        import csv
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
                n=names[i]
                data[n].append(row[i])
        D={}
        #print data.keys()
        count=0
        for protein in data['name']:
            D[protein]={}
            for prop in ['longform','dtm','predicted']:
                if prop=='longform':
                    nm=[]
                    for mut in data[prop][count].split('+'):
                        mut=mut[1:-2]
                        nm.append(mut)
                    import string
                    mut=string.join(nm,'+')
                    D[protein][prop]=mut
                else:
                    D[protein][prop]=float(data[prop][count])
            count=count+1
        return D


if __name__=='__main__':
    X=muts('/Volumes/10R_data/10R_stabs.csv')