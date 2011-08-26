#!/usr/bin/python
"""
Plot FLASH HDF5 data using matplotlib and mayavi

Examples
--------
>>> flash.py -f mira-wind-r20AU_hdf5_plt_cnt_0090
"""

from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import flashhdf5
from enthought.mayavi import mlab

# Parse command line strings
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename", help="Input file to read data from")
args = parser.parse_args()

# Read 2D data file
mirahdf5 = flashhdf5.FlashHDF5('/mnt/cdrom/mako/data6/mira-wind-70AU/t-3e2/mira-wind-70AU_hdf5_plt_cnt_0412')
dens = np.log10(mirahdf5.get_var('dens'))
r, t = np.meshgrid(mirahdf5.x, np.append(mirahdf5.y, mirahdf5.y[-1]+mirahdf5.dy_fine))
xarr = r*np.cos(t)
yarr = r*np.sin(t)

# Read 3D data file
hdf5 = flashhdf5.FlashHDF53D(args.filename)
plot_var = np.log10(hdf5.get_var('dens'))
print plot_var.shape

# plt.contourf(plot_var[:,8,:])
# plt.show()

# mlab.contour3d(plot_var)
# mlab.pipeline.volume(mlab.pipeline.scalar_field(plot_var), vmin=0.2, vmax=0.8)
mlab.surf(dens, warp_scale="auto")
# mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(dens))
#                                     plane_orientation='x_axes',
#                                     slice_index=40)
# mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(plot_var),
#                                     plane_orientation='y_axes',
#                                     slice_index=20)
mlab.outline()
# mlab.savefig('Mira.png', size=(400,300))
mlab.show()
