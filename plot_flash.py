#!/usr/bin/python2.5
"""
Plot FLASH HDF5 data using matplotlib
"""

import numpy, sys, os
from scipy import mgrid, ndimage
from tables import *
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as P

parser = OptionParser()

parser.add_option("-p", "--polar", action="store_true", dest="polar", default=False,
		help="plot in polar coordinates")
parser.add_option("-c", "--cartesian", action="store_false", dest="polar",
		help="plot in Cartesian coordinates")
parser.add_option("-l", "--log", action="store_true", dest="log", default=True, help="plot in log scale")
parser.add_option("--linear", action="store_false", dest="log", help="plot in linear scale")
parser.add_option("-n", "--ns", dest="ns", type="int", default=128, help="number of colors")
parser.add_option("-f", "--file", dest="filename", help="Input file to read data from")
parser.add_option("-s", "--save", action="store_true", default=False, dest="output", help="Save to output image")
parser.add_option("-d", "--dist", dest="dist", default=1, type="float", help="distance physical units")
parser.add_option("-b", "--bar", action="store_true", default=True, dest="bar", help="Print color bar")
parser.add_option("-v", "--vect", action="store_true", default=False, dest="vect", help="Print velocity arrows")

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
velx = numpy.zeros((nx,ny))
vely = numpy.zeros((nx,ny))
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
        ( dens_data[cur_blk,0,:,:], numpy.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
        velx[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( velx_data[cur_blk,0,:,:], numpy.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
        vely[xind:xend,yind:yend] = ndimage.map_coordinates \
        ( vely_data[cur_blk,0,:,:], numpy.array([xgrid/scaling, ygrid/scaling] ), \
        prefilter=False).transpose()
    else:
        plot_var[xind:xend,yind:yend] = dens_data[cur_blk,0,:,:].transpose()
        velx[xind:xend,yind:yend] = velx_data[cur_blk,0,:,:].transpose()
        vely[xind:xend,yind:yend] = vely_data[cur_blk,0,:,:].transpose()
h5file.close()

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

if options.log: plot_var = numpy.log10(plot_var) # log scale
if options.polar == False: # polar
	r,t = numpy.meshgrid(options.dist*x,numpy.append(y,y[-1]+dy_fine))
	xarr = r*numpy.cos(t)
	yarr = r*numpy.sin(t)
	levmin,levmax = numpy.floor(numpy.log10(plot_var.min())), numpy.ceil(numpy.log10(plot_var.min())+1)
	lev_exp = numpy.arange(levmin,levmax,(levmax-levmin)/options.ns)
	levs = numpy.power(10, lev_exp)
	P.contourf(xarr,yarr,numpy.append(plot_var,plot_var[:,0:1],axis=1).transpose(),options.ns)
# 	P.contourf(xarr,yarr,numpy.append(plot_var,plot_var[:,0:1],axis=1).transpose(),levs,
# 			norm=matplotlib.colors.LogNorm())
# 			locator=matplotlib.ticker.LogLocator())
# 	P.pcolormesh(xarr,yarr,plot_var.transpose())
# 			norm=colors.LogNorm())
	xlim = options.dist*xrange[1]
	P.axis([-xlim,xlim,-xlim,xlim])
	P.axis('scaled')
	P.xlabel("x [AU]")
	P.ylabel("y [AU]")
else: # cartesian
# 	P.contourf(x*options.dist,y,plot_var.transpose(),options.ns)
# 	P.axis([options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]])
	P.imshow(plot_var.transpose(), aspect="auto", interpolation="hanning",
	         extent=(options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]))
	P.xlabel("r [AU]")
	P.ylabel("azimuth")
if options.bar: P.colorbar()
# if options.vect: P.quiver(x*options.dist,y,velx,vely) # plot arrows
if options.output:
	P.savefig (os.path.basename(options.filename)+".png")
else:
	P.show()
