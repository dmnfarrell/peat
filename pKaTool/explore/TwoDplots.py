#!/usr/bin/env python
# Copied and modified from matplotlib reference
#
#  2Dplots.py
#  
#
#  Created by Jens on 30/05/2010.
#  Copyright (c) 2010 University College Dublin. All rights reserved.
#

import numpy 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

"""

Example: suppose you want red to increase from 0 to 1 over the bottom
half, green to do the same over the middle half, and blue over the top
half.  Then you would use:

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0))}

If, as in this example, there are no discontinuities in the r, g, and b
components, then it is quite simple: the second and third element of
each tuple, above, is the same--call it "y".  The first element ("x")
defines interpolation intervals over the full range of 0 to 1, and it
must span that whole range.  In other words, the values of x divide the
0-to-1 range into a set of segments, and y gives the end-point color
values for each segment.

Now consider the green. cdict['green'] is saying that for
0 <= x <= 0.25, y is zero; no green.
0.25 < x <= 0.75, y varies linearly from 0 to 1.
x > 0.75, y remains at 1, full green.

If there are discontinuities, then it is a little more complicated.
Label the 3 elements in each row in the cdict entry for a given color as
(x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
value is interpolated between y1[i] and y0[i+1].

Going back to the cookbook example, look at cdict['red']; because y0 !=
y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
but then it jumps down, so that for x from 0.5 to 1, red increases from
0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
back to 0, and ramps back to 1 as x goes from 0.5 to 1.

row i:   x  y0  y1
                /
               /
row i+1: x  y0  y1

Above is an attempt to show that for x in the range x[i] to x[i+1], the
interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
never used.

"""



cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }

cdict2 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 1.0),
                   (1.0, 0.1, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.1),
                   (0.5, 1.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

cdict3 = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

def heatmap(values,title='Effect of mutations',firstkey='X',secondkey='Y',colorbar=True,zlabel='Z',firstticks=False,secondticks=False):
    # Now we will use this example to illustrate 3 ways of
    # handling custom colormaps.
    # First, the most direct and explicit:

    blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

    # Second, create the map explicitly and register it.
    # Like the first method, this method works with any kind
    # of Colormap, not just
    # a LinearSegmentedColormap:

    blue_red2 = LinearSegmentedColormap('BlueRed2', cdict2)
    #plt.register_cmap(cmap=blue_red2)

    # Third, for LinearSegmentedColormap only,
    # leave everything to register_cmap:

    #plt.register_cmap(name='BlueRed3', data=cdict3) # optional lut kwarg
    #
    # Reformat the data
    #
    X=[]
    Y=[]
    Z=[]
    xvals=values.keys()
    xvals.sort()
    for xval in xvals:
        xs=[]
        ys=[]
        zs=[]
        yvals=values[xval].keys()
        yvals.sort()
        for yval in yvals:
            xs.append(xval)
            ys.append(yval)
            zs.append(values[xval][yval])
        X.append(xs)
        Y.append(ys)
        Z.append(zs)

    #plt.figure(figsize=(10,10))
    plt.imshow(Z,interpolation='nearest') #,vmin=0.0,vmax=2.8)#, interpolation='nearest', cmap=blue_red1)
    plt.xlabel(secondkey)
    plt.ylabel(firstkey)
    #
    # Plot ticks if so asked
    #
    if firstticks:
        if firstticks is True:
            plt.yticks(range(len(xvals)),xvals)
        else:
            plt.yticks(firstticks[0],firstticks[1])
    if secondticks:
        if secondticks is True:
            plt.xticks(range(len(yvals)),yvals,rotation='vertical')
        else:
            plt.xticks(secondticks[0],secondticks[1],rotation='vertical')
    if colorbar:
        cbar=plt.colorbar()
        cbar.set_label(zlabel)
    #plt.xlim([-4,4])

    # Add the title
    plt.title(title)
    import numpy
    #plt.xticks(numpy.arange(0,121,10))
    #plt.yticks(numpy.arange(0,165,10))
    plt.savefig('heatmap.png',dpi=300)
    plt.show()
    return
    
if __name__=='__main__':
    values={}
    for x in range(100):
        values[x]={}
        for y in range(100):
            values[x][y]=x*y
    heatplot(values)
