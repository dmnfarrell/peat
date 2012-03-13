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
import sys
sys.path.append('/Users/nielsen/bin/APBS/tools/python/vgrid')
from vgrid import *
NMAX=5

class Epsmap_class:

    def __init__(self,xmap,ymap,zmap):
        #print
        #print 'Reading initial Epsmaps'
        #print
        self.maps={}
        for name,emap in [['x',xmap],['y',ymap],['z',zmap]]:
            print 'Reading',emap
            data=[]
            value=0.0
            
            startVio()
            import sys
            from sys import stdout, stderr
            grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data)
            Vgrid_readDX(grid, "FILE", "ASC", "", emap)

            nx = grid.nx
            ny = grid.ny
            nz = grid.nz
            hx = grid.hx
            hy = grid.hy
            hzed = grid.hzed
            xmin = grid.xmin
            ymin = grid.ymin
            zmin = grid.zmin
            # Get the grid data
            grid_data=[]
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):                    
                        inval=0.0
                        pt = [xmin + i*hx, ymin + j*hy, zmin + k*hzed]
                        ret, value = Vgrid_value(grid, pt, inval)
                        grid_data.append(value)
                        if not ret:
                            raise Exception('Grid point not found')
            #
            # Create our own epsmaps
            #
            self.maps[name]=grid_class(nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin,grid_data)
            delete_vgrid(grid)
        #
        # Create a coarser grid for finding the volume that has the biggest effect on the results
        #
        self.basemap=self.maps['x'] # Use x-shifted map
        #print 'Initial epsmaps have been read'
        return

    #
    # ------
    #
        
    def get_cubes(self,sidelength=5,eps=80.0,clearmaps=True):
        """Initialize the cube arrays"""
        self.cubes={}
        self.cubelength={}
        self.sidelength=sidelength
        map=self.basemap
        import numpy
        self.cubegrid={} # Cube grid holds the cube number as a function of x, y and z counts (not coordinates!)
        cubecount=0
        countx=0
        #print numpy.arange(map.xmin+0.5*sidelength,map.xmin+map.nx*map.hx,sidelength)
        #print map.xmin,map.xmin+map.nx*map.hx,sidelength
        for x in numpy.arange(map.xmin+0.5*sidelength,map.xmin+map.nx*map.hx,sidelength):
            county=0
            self.cubegrid[countx]={}
            for y in numpy.arange(map.ymin+0.5*sidelength,map.ymin+map.ny*map.hy,sidelength):
                self.cubegrid[countx][county]={}
                countz=0
                for z in numpy.arange(map.zmin+0.5*sidelength,map.zmin+map.nz*map.hz,sidelength):
                    self.cubes[cubecount]={'coord':numpy.array([x,y,z]),'eps':eps,'cubegrid':[countx,county,countz]}
                    self.cubegrid[countx][county][countz]=cubecount
                    # Update counters
                    cubecount=cubecount+1
                    countz=countz+1
                #
                county=county+1
            #
            countx=countx+1
        #print '%d cubes initialized with sidelength %5.2f and eps: %5.2f' %(cubecount,sidelength,eps)
        #print 'nx: %d, ny: %d, nz: %d' %(countx,county,countz)
        #
        # Clear all the epsmaps
        #
        if clearmaps:
            self.clearmaps(value=eps)
        return sorted(self.cubes.keys())        
        
    #
    # -------
    #
    
    def set_cubeeps(self,cubenum,value):
        """Set the eps of the cube to value. This only sets the value in self.cubes"""
        #print 'Setting cubenum %d to %5.2f' %(cubenum,value)
        self.cubes[cubenum]['eps']=value
        return
        
    def get_cubeeps(self,cubenum):
        """Return the eps of this cube"""
        return self.cubes[cubenum]['eps']
        
    #
    # -------
    #
    
    def set_and_write_epsmap(self,dirname,name):
        """Set the values in all three maps according to self.cubes, and write the maps"""
        cubes=self.cubes.keys()
        for mapname in self.maps.keys():
            self.maps[mapname].clearmap(80.0)
            for cube in cubes:
                if int(self.cubes[cube]['eps'])!=80:
                    cubecenter=self.cubes[cube]['coord']
                    #print 'Setting box with center', cubecenter
                    #print 'to EPS = %5.2f' %(self.cubes[cube]['eps'])
                    #print
                    self.maps[mapname].set_box(cubecenter[0],cubecenter[1],cubecenter[2],self.sidelength,self.cubes[cube]['eps'])
                    #self.set_cubeeps(cube,self.cubes[cube]['eps'])
        # Write
        #self.write_maps('/Users/nielsen/lib/development/python/pKa/Ghost_full/testmap','testmap')
            
        return self.write_maps(dirname,name)
        
    #
    # ------
    #
        
    def get_neighbours(self,cubenum):
        """Get the neighbours to this cube using the cubegrid indeces"""
        index=self.cubes[cubecount]['cubegrid']
        neighbours=[]
        for x in range(index[0]-1,index[0]+2):
            for y in range(index[1]-1,index[1]+2):
                for z in range(index[2]-1,index[2]+2):
                    try:
                        neighbours.append(self.cubegrid[x][y][z])
                    except:
                        print 'Cube outside grid',x,y,z
        return neighbours
        
    #
    # ------
    #
        
    def write_maps(self,dirname,name):
        filenames=[]
        for dimension in ['x','y','z']:
            map=self.maps[dimension]
            filename='%s/%s%s.dx' %(dirname,dimension,name)
            self.writemap(map,filename)
            filenames.append(filename)
        return filenames
    
    #
    # ------
    #
        
    def writemap(self,map,filename='tester.dx'):
        """Write this map in dx format"""
        mydata=[]
        #count=0
        for k in range(map.nz):
            for j in range(map.ny):
                for i in range(map.nx):
                    #index = map.get_index(i,j,k)
                    #print index,count,i,j,k
                    #if index!=count:
                    #    raise Exception('Error in reformatting map')
                    mydata.append(map.get_value(i,j,k))
                    #count=count+1
        #
        # Create the grid
        #
        #print mydata
        mygrid = Vgrid_ctor(map.nx, map.ny, map.nz,
                            map.hx, map.hy, map.hz,
                            map.xmin, map.ymin, map.zmin, mydata)
        Vgrid_writeDX(mygrid, "FILE", "ASC", "", filename,'What is this', null_array())
        delete_vgrid(mygrid)
        #print 'Writemap wrote %s' %filename
        return
        
    #
    # -----
    #
    
    def clearmaps(self,value=80):
        """Clear the potential maps"""
        for mapname in self.maps.keys():
            self.maps[mapname].clearmap(value=value)
        return
    
#
# -------
#    
        
class grid_class:

    """Class for handling 3D grids"""

    def __init__(self,nx,ny,nz,hx,hy,hz,xmin,ymin,zmin,data):
        """Initialize the map"""
        self.data=[]
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.hx=hx
        self.hy=hy
        self.hz=hz
        self.xmin=xmin
        self.ymin=ymin
        self.zmin=zmin
        self.data=data[:]
        return
        
    def get_index(self,i,j,k):
        return (k*self.nx*self.ny + j*self.nx + i)
        
    def get_coord(self,i,j,k):
        """Get the x,y,z coordinate of i,j,k"""
        import numpy
        return numpy.array(x*self.hx+self.xmin,y*self.hy+self.ymin,z*self.hz+self.zmin)
        
    def get_value(self,i,j,k):
        """Get the value at i,j,k (index values not coordinates"""
        index=self.get_index(i,j,k)
        return float(self.data[index])

    def get_coord_value(self,x,y,z):
        """Get the value at coordinates x, y, z
        No interpolation is done - so use only for approximate calcs..."""
        i=int((x-self.xmin)/self.hx)
        j=int((y-self.ymin)/self.hy)
        k=int((z-self.zmin)/self.hz)
        return self.get_value(i,j,k)

    def clearmap(self,value):
        """Fill the map with this value"""
        numvals=len(self.data)
        self.data=[value]*numvals
        return
        
    def set_value(self,i,j,k,value):
        index=self.get_index(i,j,k)
        self.data[index]=value
        return
        
    def set_sphere(self,x,y,z,radius,value):
        """Set all grid points within radius of (x,y,z) to value"""
        xf=int((x-self.xmin)/self.hx)
        yf=int((y-self.ymin)/self.hy)
        zf=int((z-self.zmin)/self.hz)
        
        dx=radius/self.hx+2
        dy=radius/self.hy+2
        dz=radius/self.hz+2
        for tx in range(xf-dx,xf+dx):
            for ty in range(yf-dy,yf+dy):
                for tz in range(zf-dz,zf+dz):
                    if self.dist(self.get_coord(tx,ty,tz),self.get_coord(x,y,z))<=radius:
                        if self.is_valid(tx,ty,tz):
                            self.set_value(tx,ty,tz,value)
        #
        # Done
        #
        return
    
    #
    # -----
    #        
    
    def set_box(self,x,y,z,sidelength,value):
        """Set all grid points within the box to value"""
        xf=int((x-self.xmin)/self.hx)
        yf=int((y-self.ymin)/self.hy)
        zf=int((z-self.zmin)/self.hz)
        #
        #print 'index in 65 cubed map',xf, yf, zf
        dx=abs(int(0.5*sidelength/self.hx))+1
        dy=abs(int(0.5*sidelength/self.hy))+1
        dz=abs(int(0.5*sidelength/self.hz))+1
        for tx in range(xf-dx,xf+dx):
            if tx<0 or tx>=self.nx:
                continue
            for ty in range(yf-dy,yf+dy):
                if ty<0 or ty>=self.ny:
                    continue
                for tz in range(zf-dz,zf+dz):
                    if tz<0 or tz>=self.nz:
                        continue
                    #if self.is_valid(tx,ty,tz):
                    self.set_value(tx,ty,tz,value)
                    #else:
                    #    print 'Not a valid grid point',tx,ty,tz
                    #    #raise Exception()
        return 

    #
    # ------
    #
        
    def is_valid(self,x,y,z):
        """Is this set of coordinates within the map/grid?"""
        if x>=0 and x<self.nx:
            if y>=0 and y<self.nx:
                if z>=0 and z<self.nz:
                    return True
        return False
   
    #
    # ------
    #
                            
    def dist(self,coord1,coord2):
        """Get the distance between two coordinates"""
        vector=coord1-coord2
        s=0.0
        for element in vector:
            s=s+element*element
        import math
        return math.sqrt(s)
        
                        
        
        
        
    
        
