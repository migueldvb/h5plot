#!/usr/bin/python
"""Plot FLASH HDF5 data using matplotlib

Examples
--------
flash.py -f mira-wind-r20AU_hdf5_plt_cnt_0090
"""

import os
from argparse import ArgumentParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import flashhdf5

parser = ArgumentParser()
parser.add_argument("-p", "--polar", action="store_true", dest="polar", default=False,
        help="plot in polar coordinates")
parser.add_argument("-c", "--cartesian", action="store_false", dest="polar",
        help="plot in Cartesian coordinates")
parser.add_argument("-l", "--log", action="store_true", dest="log", default=True, 
        help="plot in log scale")
parser.add_argument("--linear", action="store_false", dest="log", help="plot in linear scale")
parser.add_argument("-n", "--ns", dest="ns", type=int, default=128, help="number of colors")
parser.add_argument("-f", "--file", dest="filename", help="Input file to read data from")
parser.add_argument("-s", "--save", action="store_true", default=False, dest="save", 
        help="Save to output image")
parser.add_argument("-o", "--output", action="store", default="", dest="out", help="Save to directory")
parser.add_argument("-d", "--dist", dest="dist", default=1, type=float, help="distance physical units")
parser.add_argument("-b", "--bar", action="store_true", default=False, dest="bar",
        help="Print color bar")
parser.add_argument("-v", "--vect", action="store_true", default=False, dest="vect", 
        help="Print velocity arrows")
parser.add_argument("-e", "--ext", default="png", dest="ext", help="Extension of output file")
parser.add_argument("-x", "--noaxis", action="store_false", default=True, dest="axis", 
        help="Do not print axis")
parser.add_argument("-t", "--thumb", action="store_true", default=False, dest="thumb", 
        help="create thumbnail")
parser.add_argument("--block", action="store_true", default=False, dest="block",
        help="Print grid structure")
parser.add_argument("--slice", dest="slice", type=float, default=1., help="zoom in")
parser.add_argument("--fliplr", action="store_true", dest="fliplr", default=False, 
        help="flip left right")
args = parser.parse_args()

hdf5 = flashhdf5.FlashHDF5(args.filename)
plot_var = hdf5.get_var('dens')
nx, ny = plot_var.shape
xrange = hdf5.xrange
yrange = hdf5.yrange

# plot contour maps using matplotlib.pyplot
if args.log: plot_var = np.log10(plot_var) # log scale
if args.thumb: plt.figure(figsize=(5,4))

if args.polar: 
    # convert polar to Cartesian
    r,t = np.meshgrid(args.dist*hdf5.x, np.append(hdf5.y, hdf5.y[-1]+hdf5.dy_fine))
    xarr = r*np.cos(t)
    yarr = r*np.sin(t)
    levmin, levmax = np.floor(np.log10(plot_var.min())), np.ceil(np.log10(plot_var.min())+1)
    lev_exp = np.arange(levmin, levmax, (levmax-levmin)/args.ns)
    levs = np.power(10, lev_exp)
    xlim = args.dist*xrange[1]
#     plt.contourf(xarr, yarr, np.append(plot_var,plot_var[:,0:1],axis=1),args.ns)
#   plt.pcolor(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1))
#   plt.contourf(xarr,yarr,np.append(plot_var,plot_var[:,0:1],axis=1),levs,
#           norm=matplotlib.colors.LogNorm())
#           locator=matplotlib.ticker.LogLocator())
    plt.pcolormesh(xarr, yarr, plot_var)
#           norm=colors.LogNorm())
#   plt.axis('scaled')
    plt.axis([-xlim,xlim,-xlim,xlim])
    plt.xlabel("$x$ [R$_{\odot}$]")
    plt.ylabel("$y$ [R$_{\odot}$]")
    plt.axes().set_aspect('equal')
else: # do not transform coordinates 
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
#   plt.contourf(x*args.dist,y,plot_var,args.ns)
#   plt.axis([args.dist*xrange[0],args.dist*xrange[1],yrange[0],yrange[1]])
    range = [args.slice*args.dist*i for i in xrange+yrange]
    if args.fliplr: plot_var = np.fliplr(plot_var)
    plt.imshow(plot_var[int(nx*(1.-args.slice)/2.):int(nx*(1+args.slice)/2.),
        int(ny*(1.-args.slice)/2.):int(ny*(1+args.slice)/2.)],
        aspect="auto", interpolation="hanning", extent=range)
    plt.xlabel("r [R$_{\odot}$]")
    plt.ylabel("$\phi$ [${\pi}$]")

if args.block:
    for cur_blk in index_good[0]: # find boundaries of current block
        if gid[cur_blk,6] > 0 or coord[cur_blk,0]<range[0] or coord[cur_blk,0]>range[1] or \
            range[2]>coord[cur_blk,1] or coord[cur_blk,1]>range[3]: continue # cur_block has children nodes or is out of range
        x0 = bnd_box[cur_blk,0,0]
        x1 = bnd_box[cur_blk,0,1]
        y0 = bnd_box[cur_blk,1,0]
        y1 = bnd_box[cur_blk,1,1]
        plt.plot([x0,x0,x1],[y0,y1,y1],'k')
    plt.axis(range)

if not args.axis: plt.axis('off')
if args.bar: plt.colorbar()
# if args.vect: plt.quiver(x*args.dist,y,velx,vely) # plot arrows
if args.out != "":
    plt.savefig(args.out)
elif args.save:
    plt.savefig(os.path.basename(args.filename)+"."+args.ext) # bbox_inches="tight")
else:
    plt.show()
