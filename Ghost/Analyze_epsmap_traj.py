#!/usr/bin/env python

class analyse_epsmap:

    def __init__(self,mindir,options):
        """Set dir and other values"""
        self.topdir=mindir
        self.options=options
        self.epsmap=None
        for snapshot in range(options.begin,options.end):
            self.get_map(snapshot)
        return
        
        
    def get_map(self,snapshot_number):
        """Save a dx map for this epsmap"""
        print 'Snapshot',snapshot_number
        import os
        ss_file=os.path.join(self.topdir,'trajectory','map_step_%d.pickle' %snapshot_number)
        fd=open(ss_file)
        import pickle
        ss=pickle.load(fd)
        fd.close()
        #
        # See if we have a selv.epsmap instance
        #
        if not self.epsmap:
            import Epsmap_class
            startdir=os.path.join(self.topdir,'start_values')
            names={}
            for name in ['x','y','z']:
                names[name]=os.path.join(startdir,'%sdiel.dx' %name)
            self.epsmap=Epsmap_class.Epsmap_class(names['x'],names['y'],names['z'])
            self.cubes=self.epsmap.get_cubes(sidelength=ss['sidelength'],eps=80.0)
        #
        # Set the eps values
        #
        for cube in ss['cubeeps'].keys():
            self.epsmap.set_cubeeps(cube,ss['cubeeps'][cube])
        #
        # Write the shifted epsmaps
        #
        resultdir=os.path.join(self.topdir,'dxmaps')
        if not os.path.isdir(resultdir):
            os.mkdir(resultdir)
        print 'Writing snapshot number %d to %s' %(snapshot_number,resultdir)
        self.epsmap.set_and_write_epsmap(resultdir,'map_%d' %snapshot_number)
        return
        
        
if __name__=='__main__':
    import optparse
    parser=optparse.OptionParser()
    parser.add_option('-b','--begin',dest='begin',type='int',help='Beginning snapshot. Default: %default',default=1)
    parser.add_option('-e','--end',dest='end',type='int',help='Ending snapshot. Default: %default',default=2)
    parser.add_option('-n','--num',dest='num',type='int',help='Number of snapshots to analyze. Overrides --end. Default: %default',default=1)
    #
    # Parse the input
    #
    (options, args) = parser.parse_args()
    
    if options.begin!=1 and options.end==2:
        options.end=options.begin+options.num
    elif options.num!=1:
        options.end=options.start+options.num
    #
    # Initialize the class
    #
    import os
    X=analyse_epsmap(os.getcwd(),options)
    
