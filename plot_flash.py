#!/usr/bin/python

from tables import *
from numpy import *
import sys
from plplot import *
from scipy import ogrid, mgrid, ndimage


def cmap1_init():
    i = [0.0,    # left boundary
           1.0]   # right boundary
    h = [240,    # blue -> green -> yellow ->
           0]     # -> red
    l = [0.6, 0.6]
    s = [0.8, 0.8]
    plscmap1n (256)
    plscmap1l(0, i, h, l, s)
    plscol0 (0,255,255,255)
    plscol0 (1,0,0,0)

filename = sys.argv[1]
h5file = openFile(filename, "r")
coord = h5file.getNode('/coordinates').read()
size = h5file.getNode('/block size').read()
bnd_box = h5file.getNode( '/bounding box').read()
amr_data = h5file.getNode('/dens').read()
node_type = h5file.getNode('/node type').read()
lrefine = h5file.getNode('/refine level').read()
index_good = numpy.where(node_type == 1)
lwant = max(lrefine[index_good])

nxb = 8
nyb = 8
xrange = [min(bnd_box[:,0,0]), max(bnd_box[:,0,1])]
yrange = [min(bnd_box[:,1,0]), max(bnd_box[:,1,1])]
top_blocks = numpy.where(lrefine == 1)
ntopx = numpy.where(bnd_box[:,1,0].flat[top_blocks] == yrange[0])
ntopx = ntopx[0].size
ntopy = numpy.where(bnd_box[:,0,0].flat[top_blocks] == xrange[0])
ntopy = ntopy[0].size
dx_fine = (xrange[1]-xrange[0])/(ntopx*nxb*2**(lwant-1))
dy_fine = (yrange[1]-yrange[0])/(ntopy*nyb*2**(lwant-1))
nx = long (ntopx*nxb*2**(lwant-1))
ny = long (ntopy*nyb*2**(lwant-1))
plot_var = numpy.zeros((nx,ny))
# plot_var = numpy.arange(nx*ny).reshape(nx,ny)

x = (numpy.arange(nx)+.5)*dx_fine + xrange[0]
y = (numpy.arange(ny)+.5)*dy_fine + yrange[0]

for cur_blk in index_good[0]:
    scaling = 2**(lwant - lrefine[cur_blk])
    # find out where the master array should live
    xind = numpy.where(x > bnd_box[cur_blk,0,0])
    xind = xind[0][0]
    yind = numpy.where(y > bnd_box[cur_blk,1,0])
    yind = yind[0][0]
    xspan = scaling*nxb
    xend = xind + xspan
    yspan = scaling*nyb
    yend = yind + yspan
    if scaling > 0:
        xgrid,ygrid=mgrid[0:xspan,0:yspan]
#         coeffs = ndimage.spline_filter(amr_data[cur_blk,0,:,:])
        plot_var[xind:xend,yind:yend] = ndimage.map_coordinates \
        (amr_data[cur_blk,0,:,:], array([ygrid/scaling, xgrid/scaling]), \
        prefilter=False)
    else:
        plot_var[xind:xend,yind:yend] = amr_data[cur_blk,0,:,:]
 
h5file.close()

zmin = min(log(plot_var.flat))
zmax = max(log(plot_var.flat))
ns = 16
shedge = zmin + (zmax - zmin) * (arange(ns))/(ns-1.)
fill_width = 2
cont_color = 0
cont_width = 0

plsdev("xwin")
cmap1_init()
plinit()
plenv(float(xrange[0]), float(xrange[1]), float(yrange[0]), float(yrange[1]),0,0)
plshades (log(plot_var), float(xrange[0]), float(xrange[1]), float(yrange[0]), float(yrange[1]), shedge, fill_width, 1)
# plshades (log(plot_var), .2,2,-3,3, shedge, fill_width, 1)
plend()
