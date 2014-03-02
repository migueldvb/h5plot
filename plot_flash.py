#!/usr/bin/python
"""Plot FLASH HDF5 data using matplotlib"""

import sys
import os
import math
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
parser.add_option("-q", "--equal", action="store_true", default=False, dest="equal", help="Make axis equal")
parser.add_option("-x", "--noaxis", action="store_false", default=True, dest="axis", help="Do not print axis")
parser.add_option("-t", "--thumb", action="store_true", default=False, dest="thumb", help="create thumbnail")
parser.add_option("--stream", action="store_true", default=False, dest="stream", help="plot streamlines")
parser.add_option("--block", action="store_true", default=False, dest="block",
        help="Print grid structure")
parser.add_option("--slice", dest="slice", type="float", default=1, help="zoom in")
(options, args) = parser.parse_args()

def stream(pxinit,pyinit,dirt):
# Read initial position and direction
    segment_nr = 4e2
    iter_nr = 1e2
    dt = 1e-4
    drint = 0.01

    rx = xrange[0] +(pxinit+0.5)*dx_fine
    ry = yrange[0] +(pyinit+0.5)*dy_fine

    xpos = np.array((rx,rx))
    ypos = np.array((ry,ry))

    for i in np.arange(1,segment_nr):
        for j in np.arange(1,iter_nr):
            vel = get_velocity(rx,ry)
            if (geometry == "polar"):
                ry = ry + dirt*dt*vel[1]/rx
            elif (geometry == "Cartesian"):
                ry = ry + dirt*dt*vel[1]
            rx = rx + dirt*dt*vel[0]
#             if (( drint < math.sqrt((xpos[1]- rx)**2+(ypos[1]- ry)**2) and (geometry == "Cartesian"))): break
#             break if (( drint < sqrt((xpos[1]- rx)**2+(xpos[1]*ypos[1]- rx*ry)**2) && (geometry == "polar")))

        xpos[0] = xpos[1]
        xpos[1] = rx
        ypos[0] = ypos[1]
        ypos[1] = ry

        if ((rx < xrange[0]) or (rx > xrange[1])): break
#         if (((ry < yrange[0]) or(ry > yrange[1]-0.5*dy_fine)) and (geometry == "Cartesian")): break
#       last if (((ry < yrange[0]) || (ry > yrange[1]-0.5*dy_fine)) and (geometry == "polar"))
#       last if ((0.5*drint) > sqrt((xpos[0]- rx)**2+(ypos[0]- ry)**2))
        if (((0.5*drint) > np.sqrt((xpos[0]- rx)**2+(xpos[0]*ypos[0]- rx*ry)**2)) and (geometry == "polar")): break
        if ((ry < yrange[0]) and (geometry == "polar")): ry = 2.*pi+ry
        if ((ry > yrange[1]) and (geometry == "polar")): ry = ry-2.*pi

#         if geometry == "Cartesian": ypos = np.array((ry,ry))

#         if geometry == "Cartesian":
#             xcart = np.array((xpos[0]*math.cos(ypos[0]), xpos[1]*math.cos(ypos[1])))
#             ycart = np.array((xpos[0]*math.sin(ypos[0]), xpos[1]*math.sin(ypos[1])))
#             plt.plot(xcart,ycart,'black')
#         else:
        plt.plot(options.dist*xpos,options.dist*ypos,'black')

def get_velocity(rxg,ryg):
   ix = int((rxg - xrange[0] - 0.5*dx_fine)/dx_fine)# + 1
   qx = (rxg - xrange[0] - (ix+0.5)*dx_fine)/dx_fine

   iy = int(((ryg - yrange[0] - 0.5*dy_fine)/dy_fine))# + 1
   qy = (ryg - yrange[0] - (iy+0.5)*dy_fine)/dy_fine

   # Interpolate to calculate the velocity
   ix1 = ix+1
   iy1 = iy+1
#    if ((iy >= ny) and (geometry == "polar")): iy = iy%ny 
#    if ((iy < 0) and (geometry == "polar")): iy = iy+ny 
#    if ((iy1 >= ny) and (geometry == "polar")): iy1 = iy1%ny 
#    if ((iy1 < 0) and (geometry == "polar")): iy1 = iy1+ny 
#    if ((ix >= nx) and (geometry == "polar")): ix = ix%nx 
#    if ((ix < 0) and (geometry == "polar")): ix = ix+nx 
#    if ((ix1 >= nx) and (geometry == "polar")): ix1 = ix1%nx 
#    if ((ix1 < 0) and (geometry == "polar")): ix1 = ix1+nx 

   vel = np.array((velx[ix,iy] * (1-qx)*(1-qy) + velx[ix1,iy] * qx*(1-qy) + velx[ix1,iy1] * qx*qy + velx[ix,iy1] * (1-qx)*qy,
        vely[ix,iy] * (1-qx)*(1-qy) + vely[ix1,iy] * qx*(1-qy) + vely[ix1,iy1] * qx*qy + vely[ix,iy1] * (1-qx)*qy))
   return vel

# Read HDF5 data
h5file = openFile(options.filename, "r")
print("opening file", options.filename)
coord = h5file.getNode('/coordinates').read()
size = h5file.getNode('/block size').read()
bnd_box = h5file.getNode( '/bounding box').read()
dens_data = h5file.getNode('/dens').read()
velx_data = h5file.getNode('/velx').read()
vely_data = h5file.getNode('/vely').read()
node_type = h5file.getNode('/node type').read()
lrefine = h5file.getNode('/refine level').read()
gid = h5file.getNode('/gid').read()
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
nx = int(ntopx*nxb*2**(lwant-1))
ny = int(ntopy*nyb*2**(lwant-1))
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

geometry = "Cartesian"
# plot contours using matplotlib.pyplot
if options.log: plot_var = np.log10(plot_var) # log scale
if options.thumb: plt.figure(figsize=(5,4))
if options.polar: # polar to cartesian
    r,t = np.meshgrid(options.dist*x,np.append(y,y[-1]+dy_fine))
    xarr = r*np.cos(t)
    yarr = r*np.sin(t)
    levmin,levmax = np.floor(np.log10(plot_var.min())), np.ceil(np.log10(plot_var.min())+1)
    lev_exp = np.arange(levmin,levmax,(levmax-levmin)/options.ns)
    levs = np.power(10, lev_exp)
    xlim = options.dist*xrange[1]
#     plt.contourf(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose(),options.ns)
#   plt.pcolor(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose())
#   plt.contourf(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1).transpose(),levs,
#           norm=matplotlib.colors.LogNorm())
#           locator=matplotlib.ticker.LogLocator())
    plt.pcolormesh(xarr,yarr,plot_var.transpose())
#           norm=colors.LogNorm())
#   plt.axis('scaled')
    plt.axis([-xlim,xlim,-xlim,xlim])
else: # do not transform coordinates 
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
#     plt.contourf(x*options.dist,y,plot_var.transpose(),options.ns)
#   plt.axis([options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]])
#     range = [options.slice*i for i in [options.dist*xrange[0],options.dist*xrange[1],yrange[0],yrange[1]]]
    range = [options.slice*options.dist*i for i in [xrange[0],xrange[1],yrange[0],yrange[1]]]
    plt.imshow(np.flipud(plot_var[int(nx*(1.-options.slice)/2.):int(nx*(1+options.slice)/2.),
        int(ny*(1.-options.slice)/2.):int(ny*(1+options.slice)/2.)].transpose()),
        aspect="auto", interpolation="hanning", extent=range)
if options.block:
    for cur_blk in index_good[0]: # find boundaries of current block
        if gid[cur_blk,6] > 0 or coord[cur_blk,0]<range[0] or coord[cur_blk,0]>range[1] or \
            range[2]>coord[cur_blk,1] or coord[cur_blk,1]>range[3]: continue # cur_block has children nodes or is out of range
        x0 = bnd_box[cur_blk,0,0]
        x1 = bnd_box[cur_blk,0,1]
        y0 = bnd_box[cur_blk,1,0]
        y1 = bnd_box[cur_blk,1,1]
        plt.plot(options.dist*np.array([x0,x0,x1]), options.dist*np.array([y0,y1,y1]),'k')
    plt.axis(range)
plt.xlabel("x [R$_{\odot}$]")
plt.ylabel("y [R$_{\odot}$]")
if options.bar: plt.colorbar()
if options.equal: plt.axes().set_aspect('equal')
if not options.axis: plt.axis('off')
if options.vect: # plot arrows
    ygrid, xgrid = np.meshgrid(x, y)
    plt.quiver(xgrid[40::80,40::80], ygrid[40::80,40::80],
            velx[40::80,40::80], vely[40::80,40::80],
            pivot='mid')
if options.stream:
    for i in np.arange(ny*.25,ny*0.75,100):
        stream(nx*.5,i,1)
        stream(nx*.5,i,-1)
        stream(i,ny*.5,1)
        stream(i,ny*.5,-1)
    plt.axis(range)
if options.outdir != "." or options.save:
    plt.savefig(options.outdir) # bbox_inches="tight")
#     if options.outdir.find("png") > 0:
#         plt.savefig(options.outdir) # bbox_inches="tight")
#     else:
#         plt.savefig(options.outdir+"/"+os.path.basename(options.filename)+"."+options.ext) # bbox_inches="tight")
else:
    plt.show()

