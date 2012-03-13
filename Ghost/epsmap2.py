
class epsmap2:

    def fill_actsite(self,depth,eps,mask):
        """Fill the active site to a specific depth with a dielectric value.
        Change only epsvalues that are different than mask"""
        bottom='A:0058:CD1'
        top=['A:0112:NH1','A:0062:CH2','A:0047:OG1'] # midpoint is max dist from bottom
        crds=[]
        for c in top:
            crds.append(P.getcoords(c))
        midpoint=sum(crds)/len(crds)
        print 'Midpoint',midpoint
        maxdist=midpoint-P.getcoords(bottom)
        print 'Max dist',maxdist
        #
        # Loop over all cubes and see which are in between bottom and midpoint
        #
        for cube in self.epsmap.cubes.keys():
            print self.epsmap.cubes[cube]
            stop
            
        return
        