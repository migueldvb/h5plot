#!/usr/bin/python
"""Plot FLASH HDF5 data using matplotlib

Examples
--------
flash.py -f mira-wind-r20AU_hdf5_plt_cnt_0090
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

hdf5 = flashhdf5.FlashHDF53D(args.filename)
plot_var = np.log10(hdf5.get_var('dens'))
print plot_var.shape
# plt.contourf(plot_var[:,8,:])
# plt.show()
# mlab.contour3d(plot_var)
# mlab.pipeline.volume(mlab.pipeline.scalar_field(plot_var), vmin=0.2, vmax=0.8)
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(plot_var),
                                    plane_orientation='x_axes',
                                    slice_index=40)
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(plot_var),
                                    plane_orientation='y_axes',
                                    slice_index=20)
mlab.outline()
mlab.show()
