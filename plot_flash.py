#!/usr/bin/python

import numpy
import sys, os
from pylab import *
from plplot import *
from scipy import ogrid, mgrid, ndimage
from tables import *
# from cmap1_init import *

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
    if scaling > 1:
        xgrid,ygrid=mgrid[0:xspan,0:yspan]
        plot_var[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( amr_data[cur_blk,0,:,:], array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
    else:
        plot_var[xind:xend,yind:yend] = amr_data[cur_blk,0,:,:].transpose()
 
h5file.close()

dist = 70
r,t = meshgrid(dist*x,y)
xarr = r*cos(t)
yarr = r*sin(t)
ns = 128
# im = imshow(log(plot_var.transpose()), aspect="auto",\
#         extent=(xrange[0],xrange[1],yrange[0],yrange[1]))
# im = contourf(xarr,yarr,log(plot_var.transpose()),ns)
dist = 0.042
im = contourf(x*dist,y*dist,log(plot_var.transpose()),ns)
colorbar(im)
xlabel("x [AU]")
ylabel("y [AU]")
savefig (os.path.basename(filename)+".png")

# zmin = min(log(plot_var.flat))
# zmax = max(log(plot_var.flat))
# ns = 64
# shedge = zmin + (zmax - zmin) * (arange(ns))/(ns-1.)
# x.shape = (-1,1)
# xg = x*cos(y)
# yg = x*sin(y)
# fill_width = 2
# cont_color = 0
# cont_width = 0
# plsdev("png")
# plsfnam(os.path.basename(filename)+".png")
# cmap1_init()
# plinit()
# cart=True
# xmin = float(xrange[0]*70.)
# if cart:
#     plenv(-xmin,xmin,-xmin,xmin,1,0)
#     plshades (log(plot_var), -xmin,xmin,-xmin,xmin, shedge, fill_width, 1, \
#         pltr2, xg, yg, 2)
#     pllab("x","y","")
# else:
#     plenv(xrange[0], float(xrange[1]), float(yrange[0]), float(yrange[1]),0,0)
#     plshades (log(plot_var), float(xrange[0]), float(xrange[1]), \
#         float(yrange[0]), float(yrange[1]), shedge, fill_width, 1)
# # plshades (log(plot_var), .2,2,-3,3, shedge, fill_width, 1)
# plend()
