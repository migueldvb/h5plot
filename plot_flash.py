#!/usr/bin/python2.5
"""
Plot FLASH HDF5 data using matplotlib
"""

import sys
import os
from scipy import mgrid, ndimage
from optparse import OptionParser
from tables import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-p", "--polar", action="store_true", dest="polar", default=False,
        help="plot in polar coordinates")
parser.add_option("-c", "--cartesian", action="store_false", dest="polar",
        help="plot in Cartesian coordinates")
parser.add_option("-l", "--log", action="store_true", dest="log", default=True, help="plot in log scale")
parser.add_option("--linear", action="store_false", dest="log", help="plot in linear scale")
parser.add_option("-n", "--ns", dest="ns", type="int", default=128, help="number of colors")
parser.add_option("-f", "--file", dest="filename", help="Input file to read data from")
parser.add_option("-s", "--save", action="store_true", default=False, dest="save", help="Save to output image")
parser.add_option("-o", "--output", action="store", default=".", dest="outdir", help="Save to directory")
parser.add_option("-d", "--dist", dest="dist", default=1, type="float", help="distance physical units")
parser.add_option("-b", "--bar", action="store_true", default=False, dest="bar", help="Print color bar")
parser.add_option("-v", "--vect", action="store_true", default=False, dest="vect", help="Print velocity arrows")
parser.add_option("-e", "--ext", default="png", dest="ext", help="Extension of output file")
parser.add_option("-x", "--noaxis", action="store_false", default="True", dest="axis", help="Do not print axis")
(options, args) = parser.parse_args()

# Read HDF5 data
h5file = openFile(options.filename, "r")
print "opening file", options.filename
coord = h5file.getNode('/coordinates').read()
size = h5file.getNode('/block size').read()
bnd_box = h5file.getNode( '/bounding box').read()
dens_data = h5file.getNode('/dens').read()
velx_data = h5file.getNode('/velx').read()
vely_data = h5file.getNode('/vely').read()
node_type = h5file.getNode('/node type').read()
lrefine = h5file.getNode('/refine level').read()
index_good = np.where(node_type == 1)
lwant = max(lrefine[index_good])

nxb, nyb = 8, 8
xrange = [min(bnd_box[:,0,0]), max(bnd_box[:,0,1])]
yrange = [min(bnd_box[:,1,0]), max(bnd_box[:,1,1])]
top_blocks = np.where(lrefine == 1)
ntopx = np.where(bnd_box[:,1,0].flat[top_blocks] == yrange[0])
ntopx = ntopx[0].size
ntopy = np.where(bnd_box[:,0,0].flat[top_blocks] == xrange[0])
ntopy = ntopy[0].size
dx_fine = (xrange[1]-xrange[0])/(ntopx*nxb*2**(lwant-1))
dy_fine = (yrange[1]-yrange[0])/(ntopy*nyb*2**(lwant-1))
nx = long (ntopx*nxb*2**(lwant-1))
ny = long (ntopy*nyb*2**(lwant-1))
plot_var = np.zeros((nx,ny))
velx = np.zeros((nx,ny))
vely = np.zeros((nx,ny))
# plot_var = np.arange(nx*ny).reshape(nx,ny)

x = (np.arange(nx)+.5)*dx_fine + xrange[0]
y = (np.arange(ny)+.5)*dy_fine + yrange[0]

# loop through good blocks
for cur_blk in index_good[0]:
    scaling = 2**(lwant - lrefine[cur_blk])
    # find out where the master array should live
    xind = np.where(x > bnd_box[cur_blk,0,0])
    xind = xind[0][0]
    yind = np.where(y > bnd_box[cur_blk,1,0])
    yind = yind[0][0]
    xspan = scaling*nxb
    xend = xind + xspan
    yspan = scaling*nyb
    yend = yind + yspan
    if scaling > 1:
        xgrid,ygrid=mgrid[0:xspan,0:yspan]
        plot_var[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( dens_data[cur_blk,0,:,:], np.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
        velx[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( velx_data[cur_blk,0,:,:], np.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
        vely[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( vely_data[cur_blk,0,:,:], np.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
    else:
        plot_var[xind:xend,yind:yend] = dens_data[cur_blk,0,:,:].transpose()
        velx[xind:xend,yind:yend] = velx_data[cur_blk,0,:,:].transpose()
        vely[xind:xend,yind:yend] = vely_data[cur_blk,0,:,:].transpose()
h5file.close()

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# plot contours using pyplot
if options.log: plot_var = np.log10(plot_var) # log scale
if options.polar == False: # cartesian
    r,t = np.meshgrid(options.dist*x,np.append(y,y[-1]+dy_fine))
    xarr = r*np.cos(t)
    yarr = r*np.sin(t)
    levmin,levmax = np.floor(np.log10(plot_var.min())), np.ceil(np.log10(plot_var.min())+1)
    lev_exp = np.arange(levmin,levmax,(levmax-levmin)/options.ns)
    levs = np.power(10, lev_exp)
    xlim = options.dist*xrange[1]
    plt.contourf(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose(),options.ns)
#   plt.pcolor(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose())
#   plt.contourf(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose(),levs,
#           norm=matplotlib.colors.LogNorm())
#           locator=matplotlib.ticker.LogLocator())
#   plt.pcolormesh(xarr,yarr,plot_var.transpose())
#           norm=colors.LogNorm())
#   plt.axis('scaled')
    plt.axes().set_aspect('equal')
    if options.axis:
        plt.axis([-xlim,xlim,-xlim,xlim])
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
    else:
        plt.axis('off')
else: # polar
#   plt.contourf(x*options.dist,y,plot_var.transpose(),options.ns)
#   plt.axis([options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]])
    plt.imshow(plot_var.transpose(), aspect="auto", interpolation="hanning",
             extent=(options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]))
    plt.xlabel("r [AU]")
    plt.ylabel("azimuth")
if options.bar: plt.colorbar()
# if options.vect: plt.quiver(x*options.dist,y,velx,vely) # plot arrows
if options.outdir != "." or options.save:
    plt.savefig(options.outdir+"/"+os.path.basename(options.filename)+"."+options.ext)
#             bbox_inches="tight")
else:
    plt.show()
