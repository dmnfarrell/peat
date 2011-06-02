#!/usr/bin/env python
import numpy as np

class analyze_explore:

    def __init__(self,options):
        #
        # Read the dir
        #
        import large_dataset_manager
        self.DM=large_dataset_manager.data_manager(options.dir)
        #
        # Read all datasets
        #
        matrix={}
        count=0
        data=self.DM(count)
        while data:
            count=count+1
            data=self.DM(count)
            if count>1000:
                data=None
            print '\b.',
            import sys
            sys.stdout.flush()
            #
            if data:
                sum_intpkadiff,sum_intene,tit_fit,pHact_fit=self.analyze_point(data)
                if not matrix.has_key(sum_intpkadiff):
                    matrix[sum_intpkadiff]={}
                if not matrix[sum_intpkadiff].has_key(sum_intene):
                    matrix[sum_intpkadiff][sum_intene]=[0,0]
                matrix[sum_intpkadiff][sum_intene][0]=matrix[sum_intpkadiff][sum_intene][0]+tit_fit
                matrix[sum_intpkadiff][sum_intene][1]=matrix[sum_intpkadiff][sum_intene][1]+pHact_fit
        #
        # Calculate ratio
        #
        minS=9999.9
        maxS=-9999.9
        for ip in matrix.keys():
            for ie in matrix[ip].keys():
                minS=min(minS,ie)
                maxS=max(maxS,ie)
                #
                vals=matrix[ip][ie]
                matrix[ip][ie]=vals[0]
        # Fill in empty values
        minF=min(matrix.keys())
        maxF=max(matrix.keys())
        print minF,maxF
        print minS,maxS
        #stop
        for ip in np.arange(minF,maxF+1,1):
            if not matrix.has_key(ip):
                matrix[ip]={}
            for ie in range(minS,maxS+1,1):
                if not matrix[ip].has_key(ie):
                    matrix[ip][ie]=0.0
        import TwoDplots
        TwoDplots.heatmap(matrix,firstkey='sum intpka diff',secondkey='sum intene (kT)')
                
        print 'I found %7.2e datasets' %count
        return
        
    #
    # ----
    #
    
    def analyze_point(self,data):
        """Analyze this datapoint"""
        numgroups=data['numgroups']
        intpkas=[]
        intenes=[]
        grptypes=[]
        config=data['configuration']
        count=0
        for group in range(numgroups):
            intpkas.append(config[count+1])
            grptypes.append(config[count])
            count=count+2
        #
        for tc in range(count,len(config)):
            intenes.append(config[tc])
        #
        IPdiff=0.0
        for IP in intpkas[1:]:
            IPdiff=IPdiff+abs(intpkas[0]-IP)/(float(len(intpkas[1:])))
        sum_IE=sum(intenes)
        # Get the fit to HH for the titration curves
        titcurvs=data['titcurvs']
        tc_slope=self.fit_curves(titcurvs)
        return float('%5.0f' %IPdiff),int(sum_IE),tc_slope,0.0
        
    #
    # -----
    #
        
    def fit_curves(self,curves):
        """Fit the titration curves to the HH equation"""
        error=[]
        import string, pKarun, pKarun.pKa_general
        PKana=pKarun.pKa_general.pKanalyse()
        for group in curves.keys():
            solution,sq=PKana.fit_to_henderson(curves[group])
            sq=solution[0]
            try:
                sq=float(sq)
            except:
                sq=0.0
            error.append(abs(sq))
        return sum(error)/len(error)
        



if __name__=="__main__":
    print
    print 'Do analyses of system exploration'
    print
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    #parser.add_option('-n','--groupnum',dest='numgroups',type='int',action='store',default=3,help='Number of titratable groups. Default: %default')
    #parser.add_option('-s','--step',dest='step',action='store',type='float',default=1.0,help='Stepsize for each parameter in pKa units. Default: %default')
    parser.add_option('-d','--dir',dest='dir',action='store',type='string',default='resultdir',help='Directory where results will be stored. Default: %default')
    (options, args) = parser.parse_args()
    #
    # Call the class
    #
    X=analyze_explore(options)