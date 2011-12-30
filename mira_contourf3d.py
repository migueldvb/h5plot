"""
.. versionadded:: 1.1.0
   This demo depends on new features added to contourf3d.
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import flashhdf5
import os

filename = '/home/miguel/project/mira/public_html/figures/r30AU/mira-wind-r30AU_hdf5_plt_cnt_0100'
hdf5 = flashhdf5.FlashHDF5(filename)
Z = hdf5.get_var('dens')
_xrange = hdf5.xrange
_yrange = hdf5.yrange
# convert to real coordinates
_range = [70*i for i in _xrange+_yrange]
X, Y = np.meshgrid(70*hdf5.x, hdf5.y)#, hdf5.y[-1]+hdf5.dy_fine))
# X = r*np.cos(t)
# Y = r*np.sin(t)
Z = np.log10(Z).transpose()

fig = plt.figure()
ax = fig.gca(projection='3d')
# X, Y, Z = axes3d.get_test_data(0.05)
cset = ax.contourf(X, Y, Z, 128, zdir='x', offset=10)
cset = ax.contourf(X, Y, Z, 128, zdir='y', offset=4)
cset = ax.contourf(X, Y, Z, 128, zdir='z', offset=-7)
ax.plot_surface(X, Y, Z, rstride=32, cstride=32, alpha=0.3)

ax.set_xlabel('X')
ax.set_xlim(10, 140)
ax.set_ylabel('Y')
ax.set_ylim(-4, 4)
ax.set_zlabel('Z')
ax.set_zlim(-7, -1.5)

# plt.show()
plt.savefig(os.path.basename(filename), dpi=300)
