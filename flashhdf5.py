#!/usr/bin/python
"""
FLASH HDF5 data class

This class is ported from the IDL routines in FLASH2.5
source/fidlr3/merge_amr.pro
source/fidlr3/read_amr.pro
"""

from scipy import mgrid, ndimage
import tables
import h5py
import numpy as np
from io import StringIO

class FlashHDF53D(object):
    """Read FLASH HDF5 data"""
    def __init__(self, filename):
        self.filename = filename
        self.h5file = tables.openFile(self.filename, "r")
        self.h5f = h5py.File(self.filename, 'r')

    def close(self):
        self.h5file.close()
        self.h5f.close()

    def n_dim(self):
        """Number of dimensions"""
        self.ndim = self.h5f['coordinates'].shape[1]

    def get_var(self, var, axis=1, zslice=0):
        """
        Interpolate data to a uniform grid

        Parameters
        ----------
        axis : int
            axis to get a slice from (0, 1, 2)
        zslice : int
            cell number in axis
        """
        print("opening file", self.filename)
        # Read node data
        self.coord = self.h5f['/coordinates']
        self.ndim = self.coord[1]
        self.bnd_box = self.h5file.getNode( '/bounding box').read()
        plot_data = self.h5file.getNode('/{0}'.format(var)).read()
        self.node_type = self.h5file.getNode('/node type').read()
        self.lrefine = self.h5file.getNode('/refine level').read()
        # number of cells in each direction
        try:
            nxb = self.h5file.getNode('/integer scalars')[0][1]
            nyb = self.h5file.getNode('/integer scalars')[1][1]
            nzb = self.h5file.getNode('/integer scalars')[2][1]
        except tables.exceptions.NoSuchNodeError:
            # FLASH2.5 simulation parameters record
            nxb = self.h5file.getNode('/simulation parameters')[0]['nxb']
            nyb = self.h5file.getNode('/simulation parameters')[0]['nyb']
            nzb = self.h5file.getNode('/simulation parameters')[0]['nzb']
#         nxb, nyb, nzb = 8, 8, 8

        self.index_good = np.where(self.node_type == 1)
        lmax = max(self.lrefine[self.index_good])
        top_blocks = np.where(self.lrefine == 1)
        self.xrange = [min(self.bnd_box[:,0,0]), max(self.bnd_box[:,0,1])]
        self.yrange = [min(self.bnd_box[:,1,0]), max(self.bnd_box[:,1,1])]
        self.zrange = [min(self.bnd_box[:,2,0]), max(self.bnd_box[:,2,1])]
        # find the number of level 1 blocks whose lower y coord is the minimum
        # y value and whose lower z coord is the minimum z value
        ntopx = np.where((
                self.bnd_box[:,1,0].flat[top_blocks] == self.yrange[0]) \
                & (self.bnd_box[:,2,0].flat[top_blocks] == self.zrange[0]))
        ntopx = ntopx[0].size
        ntopy = np.where((
                self.bnd_box[:,0,0].flat[top_blocks] == self.xrange[0]) \
                & (self.bnd_box[:,2,0].flat[top_blocks] == self.zrange[0]))
        ntopy = ntopy[0].size
        ntopz = np.where((
                self.bnd_box[:,0,0].flat[top_blocks] == self.xrange[0]) \
                & (self.bnd_box[:,1,0].flat[top_blocks] == self.yrange[0]))
        ntopz = ntopz[0].size
        self.dx_fine = (self.xrange[1]-self.xrange[0])/(ntopx*nxb*2**(lmax-1))
        self.dy_fine = (self.yrange[1]-self.yrange[0])/(ntopy*nyb*2**(lmax-1))
        self.dz_fine = (self.zrange[1]-self.zrange[0])/(ntopz*nzb*2**(lmax-1))
        nx = int(ntopx*nxb*2**(lmax-1))
        ny = int(ntopy*nyb*2**(lmax-1))
        nz = int(ntopy*nzb*2**(lmax-1))
        plot_var = np.zeros((nx, ny, nz))
        # face-centered coordinates
        self.x = (np.arange(nx)+.5)*self.dx_fine + self.xrange[0]
        self.y = (np.arange(ny)+.5)*self.dy_fine + self.yrange[0]
        self.z = (np.arange(nz)+.5)*self.dz_fine + self.zrange[0]

        # loop through good blocks
        for cur_blk in self.index_good[0]:
            scaling = 2**(lmax - self.lrefine[cur_blk])
            # find out where the master array should live
            xind = np.where(self.x > self.bnd_box[cur_blk,0,0])
            xind = xind[0][0]
            yind = np.where(self.y > self.bnd_box[cur_blk,1,0])
            yind = yind[0][0]
            zind = np.where(self.z > self.bnd_box[cur_blk,2,0])
            zind = zind[0][0]
            xspan = scaling*nxb
            xend = xind + xspan
            yspan = scaling*nyb
            yend = yind + yspan
            zspan = scaling*nzb
            zend = zind + zspan
            if scaling > 1:
                # Map array data by interpolation
                xgrid, ygrid, zgrid = mgrid[0:xspan, 0:yspan, 0:zspan]
                plot_var[xind:xend, yind:yend, zind:zend] = \
                    ndimage.map_coordinates(plot_data[cur_blk,:,:,:],
                    np.array([xgrid/scaling, ygrid/scaling, zgrid/scaling]),
                    prefilter=False).transpose()
            else:
                plot_var[xind:xend, yind:yend, zind:zend] = \
                    plot_data[cur_blk,:,:,:].transpose()
        return plot_var

class FlashHDF52D(FlashHDF53D):
    pass

    def get_var(self, var):
        """
        Interpolate data to a uniform grid
        """
        print("opening file", self.filename)
        # Read node data
        self.coord = self.h5f['/coordinates']
        size = self.h5f['/block size']
        self.bnd_box = self.h5f['/bounding box']
        plot_data = self.h5file.getNode('/{0}'.format(var)).read()
        node_type = self.h5file.getNode('/node type').read()
        self.lrefine = self.h5file.getNode('/refine level').read()
        self.gid = self.h5file.getNode('/gid').read()
        # number of cells in each direction
        try:
            nxb = self.h5file.getNode('/integer scalars')[0][1]
            nyb = self.h5file.getNode('/integer scalars')[1][1]
        except tables.exceptions.NoSuchNodeError:
            # FLASH2.5 simulation parameters record
            nxb = self.h5file.getNode('/simulation parameters')[0]['nxb']
            nyb = self.h5file.getNode('/simulation parameters')[0]['nyb']

        self.index_good = np.where(node_type == 1)
        lmax = max(self.lrefine[self.index_good])
        top_blocks = np.where(self.lrefine == 1)
        self.xrange = [min(self.bnd_box[:,0,0]), max(self.bnd_box[:,0,1])]
        self.yrange = [min(self.bnd_box[:,1,0]), max(self.bnd_box[:,1,1])]
        # find the number of level 1 blocks whose lower y coord is the minimum
        # y value and whose lower z coord is the minimum z value
        ntopx = np.where((self.bnd_box[:,1,0].flat[top_blocks] ==
            self.yrange[0]))
        ntopx = ntopx[0].size
        ntopy = np.where((self.bnd_box[:,0,0].flat[top_blocks] ==
            self.xrange[0]))
        ntopy = ntopy[0].size
        self.dx_fine = (self.xrange[1]-self.xrange[0])/(ntopx*nxb*2**(lmax-1))
        self.dy_fine = (self.yrange[1]-self.yrange[0])/(ntopy*nyb*2**(lmax-1))
        nx = int(ntopx*nxb*2**(lmax-1))
        ny = int(ntopy*nyb*2**(lmax-1))
        plot_var = np.zeros((nx, ny))
        # face-centered coordinates
        self.x = (np.arange(nx)+.5)*self.dx_fine + self.xrange[0]
        self.y = (np.arange(ny)+.5)*self.dy_fine + self.yrange[0]

        # loop through good blocks
        for cur_blk in self.index_good[0]:
            scaling = 2**(lmax - self.lrefine[cur_blk])
            # find out where the master array should live
            xind = np.where(self.x > self.bnd_box[cur_blk,0,0])
            xind = xind[0][0]
            yind = np.where(self.y > self.bnd_box[cur_blk,1,0])
            yind = yind[0][0]
            xspan = scaling*nxb
            xend = xind + xspan
            yspan = scaling*nyb
            yend = yind + yspan
            if scaling > 1:
                # Map array data by interpolation
                xgrid, ygrid = mgrid[0:xspan, 0:yspan]
                plot_var[xind:xend, yind:yend] = \
                ndimage.map_coordinates(plot_data[cur_blk,0,:,:],
                        np.array([xgrid/scaling, ygrid/scaling] ),
                        mode='nearest', prefilter=False).transpose()
            else:
                plot_var[xind:xend,yind:yend] = plot_data[cur_blk,0,:,:].transpose()
        return plot_var

class FlashNBody:
    """Read N Body data"""
    def __init__(self, filename):
        self.filename = filename
        f = open(filename)
        f.readline()
        f.readline()
        nbo_curr = np.genfromtxt(StringIO(f.readline()))
        self.primary = nbo_curr[:3:2]
        self.secondary= nbo_curr[1:4:2]
        f.close()
