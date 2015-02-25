#!/usr/bin/python
"""
Plot FLASH HDF5 data using matplotlib

Examples
--------
flash.py -f mira-wind-r20AU_hdf5_plt_cnt_0090
"""

import os
from argparse import ArgumentParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import flashhdf5

# Parse command line strings
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
parser.add_argument("-o", "--output", action="store", default="", dest="out",
        help="Save to directory")
parser.add_argument("-d", "--dist", dest="dist", default=1, type=float,
        help="distance in physical units")
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
parser.add_argument("-q", "--equal", action="store_true", default=False,
        dest="equal", help="Make axes equal")
parser.add_argument("--xz", action="store_true", default=False, dest="xz",
        help="Plot xz plane")
args = parser.parse_args()

# Read HDF5 data file
hdf5 = flashhdf5.FlashHDF52D(args.filename)
plot_var = hdf5.get_var('dens')
nx, ny = plot_var.shape
xrange = hdf5.xrange
yrange = hdf5.yrange
# convert to real coordinates
_range = [i for i in xrange+yrange]
_drange = [args.slice*args.dist*i for i in xrange+yrange]

# plot contour maps using matplotlib.pyplot
fig = plt.figure()
ax = fig.add_subplot(111)

if args.log: plot_var = np.log10(plot_var) # log scale
if args.thumb: plt.figure(figsize=(5,4)) # thumbnail output

if args.polar:
    # convert polar to Cartesian coordinates
    r, t = np.meshgrid(args.dist*hdf5.x, np.append(hdf5.y, hdf5.y[-1]+hdf5.dy_fine))
    xarr = r*np.cos(t)
    yarr = r*np.sin(t)#*np.cos(80*np.pi/180)
    levmin, levmax = np.floor(np.log10(plot_var.min())), np.ceil(np.log10(plot_var.min())+1)
    lev_exp = np.arange(levmin, levmax, (levmax-levmin)/args.ns)
    levs = np.power(10, lev_exp)
    xlim = args.dist*xrange[1]
    plt.pcolormesh(xarr, yarr, plot_var.transpose())
#             norm=matplotlib.colors.LogNorm())
#         plt.axis('scaled')
    if args.xz:
        ax.add_patch(matplotlib.patches.Wedge((0,0), args.dist*2, 45, 135,
            ec='none', fc='white'))
        ax.add_patch(matplotlib.patches.Wedge((0,0), args.dist*2, 225, 315,
            ec='none', fc='white'))
        ylim = np.sin(np.pi/4.)*xlim
        plt.axis([-xlim, xlim, -ylim, ylim])
        plt.xlabel("$x$ [AU]")
        plt.ylabel("$z$ [AU]")
    else:
        plt.axis([-xlim,xlim,-xlim,xlim])
        plt.xlabel("$x$ [AU]")
        plt.ylabel("$y$ [AU]")
#         plt.xlabel("$x$ [R$_{\odot}$]")
#         plt.ylabel("$y$ [R$_{\odot}$]")
    plt.axes().set_aspect('equal')
else: # do not transform coordinates 
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
#   plt.axis([args.dist*xrange[0],args.dist*xrange[1],yrange[0],yrange[1]])
    if args.fliplr: plot_var = np.fliplr(plot_var)
    plt.imshow(np.flipud(plot_var[int(nx*(1.-args.slice)/2.):int(nx*(1+args.slice)/2.),
        int(ny*(1.-args.slice)/2.):int(ny*(1+args.slice)/2.)].transpose()),
        aspect="auto", interpolation="hanning", extent=_drange)
    plt.xlabel("r [R$_{\odot}$]")
    plt.ylabel("$\phi$ [${\pi}$]")

if args.block: # Print grid structure
    # find boundaries of current block
    for cur_blk in hdf5.index_good[0]:
        # cur_block has children nodes or is out of range
        if hdf5.gid[cur_blk,6] > 0 or \
                hdf5.coord[cur_blk,0] < _range[0] or \
                hdf5.coord[cur_blk,0] > _range[1] or \
                hdf5.coord[cur_blk,1] < _range[2] or \
                hdf5.coord[cur_blk,1] > _range[3]: continue
        x0 = hdf5.bnd_box[cur_blk,0,0]*args.dist
        x1 = hdf5.bnd_box[cur_blk,0,1]*args.dist
        y0 = hdf5.bnd_box[cur_blk,1,0]
        y1 = hdf5.bnd_box[cur_blk,1,1]
        if args.polar:
            # plot radial spokes
            plt.plot([x0*np.cos(y0), x1*np.cos(y0)], 
                    [x0*np.sin(y0), x1*np.sin(y0)], 'k')
            # plot azimuthal arcs
            arc = matplotlib.patches.Arc((0., 0.), 2*x0, 2*x0, 0., np.degrees(y0),
                    np.degrees(y1))
            ax.add_patch(arc)
        else:
            plt.plot([x0,x0,x1], [y0,y1,y1], 'k')
    if args.polar:
        ax.add_patch(matplotlib.patches.Circle((0,0), 4*args.dist, fc='none'))
    else:
        plt.axis(_range)

if args.equal: plt.axes().set_aspect('equal')
if not args.axis: plt.axis('off')
if args.bar: plt.colorbar()
# if args.vect: plt.quiver(x*args.dist,y,velx,vely) # plot arrows
if args.out:
    plt.savefig(args.out)
elif args.save:
    plt.savefig(os.path.basename(args.filename)+"."+args.ext,
            transparent=True) # bbox_inches="tight")
else:
    plt.show()
hdf5.close()
