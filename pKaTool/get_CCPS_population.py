#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
# Copyright (C) 2010 Jens Erik Nielsen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: 
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland

class values:
    def __init__(self):
        return
        
def get_CCPS_population(pdbfile,ccps,MCsteps=200000,pHstep=0.1,pHstart=0.1,pHstop=12.0):
    V=values()
    V.MCsteps=MCsteps
    V.pHstep=pHstep
    V.pHstop=pHstop
    V.pHstart=pHstart
    V.ccps=ccps
    V.show=False
    V.list_groups=False
    V.show_plot=False
    return _get_CCPS_population(V,[pdbfile])
    

def _get_CCPS_population(options,args):
    """recalculate titration curves and get the pH-dependent population of select protonation states"""
    import pKaIO
    pdbfile=args[0]
    MCsteps=options.MCsteps
    pHstep=options.pHstep
    pHstart=options.pHstart
    pHend=options.pHstop
    IO=pKaIO.pKaIO(pdbfile)
    import pKa.pKa_MC 
    print 'Recalculating wild type titration curves'
    pKa_params={'pHstart':pHstart,'pHstop':pHend,'pHstep':pHstep,'pKMCsteps':MCsteps,'verbose':1}
    pKaCALC=pKa.pKa_MC.pKa_calculation_class(pdbfile,pKa_info=None,params=pKa_params)
    pKaCALC.set_MC_CPP()
    #
    #
    #
    groups=pKaCALC.MC.desolv.keys()
    groups.sort()
    if options.list_groups:
        print 'The PDB file contains these titratable groups:\n',groups
        return
    #
    # Set the protonation state
    #
    CCPS={}
    for spec in options.ccps:
        split=spec.split('=')
        group=split[0]
        charge=abs(int(split[1]))
        if CCPS.has_key(group):
            raise Exception('You cannot specify two protonation states for a single group: %s' %group)
        if not group in groups:
            raise Exception('Unknown group: %s' %group)
        CCPS[group]=charge
    #
    # Construct monitor_states        
    #
    monitor_states=[]
    for group in groups:
        if CCPS.has_key(group):
            monitor_states.append(CCPS[group])
        else:
            monitor_states.append(-99)
    #
    # Calculate the pKa values
    #
    print 'pHrange: step: %.2f start: %.2f, end: %.2f. MCsteps: %e ' %(pHstep,pHstart,pHend,MCsteps)
    wtpkas,titration_curves=pKaCALC.MC.calc_pKas(MCsteps,pHstep,pHstart,pHend,0,
        monitor_states=monitor_states)
    pHs=pKaCALC.MC.pH[:]
    CCPS=pKaCALC.MC.CCPS_pop[:]
    if options.show_plot:
        import matplotlib.pyplot as plt
        plt.plot(pHs,CCPS,'ro-',linewidth=2.0)
        plt.show()
    return pHs,CCPS

    
if __name__=='__main__':
    print
    print 'Calculate the pH-dependent population of the CCPS'
    print 'Jens Erik Nielsen, 2009'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <pdbfile>',version='%prog 1.0')
    parser.add_option('-p',"--protonation_state",type='string',dest='ccps',action='append',
                      help='Protonation state that should be monitored, e.g. A:0035:GLU=-1. -p must be specified multiple times for protonation states involving more than one residue.')
    parser.add_option('-n',"--pHstep",type='float',dest='pHstep',action='store',
                      help='pHstep when calculating titration curves. Deafult: %default',default=0.1)
    parser.add_option('-s',"--pHstart",type='float',dest='pHstart',action='store',
                      help='pHstart when calculating titration curves. Deafult: %default',default=0.1)
    parser.add_option('-e',"--pHstop",type='float',dest='pHstop',action='store',
                      help='pHstop when calculating titration curves. Deafult: %default',default=12.0)
    parser.add_option('-m',"--MCsteps",type='int',dest='MCsteps',default=200000,action='store',
                    help='Number of Monte Carlo steps per pH value. Default: %default')
    parser.add_option('-l',"--list_groups",dest='list_groups',default=False,action='store_true',
                    help='List titratable groups in file. Default: %default')
    parser.add_option('-x','--plot',dest='show_plot',default=False,action='store_true',
                    help='Plot the pH-dependence of the CCPS population using matplotlib. Default: %default')
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify a PDB file')
    #
    # Call main
    #
    _get_CCPS_population(options,args)
