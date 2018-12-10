# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

#####################################
# NOTE:
# This reader is under development,
# and presently not fully functional
#####################################

import logging

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator

from basereader import BaseReader, pyproj

from mesh import Mesh

class Reader(BaseReader):

    def _read_mesh(self, filename):
        ds = Dataset(filename)
        print filename
        logging.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in ds.variables:
            #print var_name
            var = ds.variables[var_name]
            #if var.ndim > 1:
            #    continue  # Coordinates must be 1D-array
            attributes = var.ncattrs()
            if 'standard_name' in attributes:
                standard_name = var.__dict__['standard_name']
            if 'long_name' in attributes:
                long_name = var.__dict__['long_name']
#            if 'axis' in attributes:
#                axis = var.__dict__['axis']
            if 'units' in attributes:
                units = var.__dict__['units']
#            if '_CoordinateAxisType' in attributes:
#                CoordinateAxisType = var.__dict__['_CoordinateAxisType']
            if var_name == 'x' or \
                    long_name == 'x-coordinates':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                x = var[:]*unitfactor
                self.unitfactor = unitfactor
                self.numx = var.shape[0]
            if var_name == 'y' or \
                    long_name == 'y-coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                y = var[:]*unitfactor
                self.numy = var.shape[0]
            if var_name == 'depth' or long_name == 'Bathymetry':
                self.zname = var_name
                if 'positive' not in var.ncattrs() or \
                        var.__dict__['positive'] == 'up':
                    z = var[:]
                else:
                    z = -var[:]
            if var_name == 'ele':
                if 'start_index' in attributes:
                    start_index = ds.variables[var_name].__dict__['start_index']
                else:
                    start_index = 0
                #print start_index
                ele = var[:] - start_index
        
        ds.close()
        
        if 'x' not in locals():
            raise ValueError('Did not find x-coordinate variable')
        if 'y' not in locals():
            raise ValueError('Did not find y-coordinate variable')
        if 'ele' not in locals():
            raise ValueError('Did not find ele variable')
        if 'z' not in locals():
            raise ValueError('Did not find depth variable')
        
        msh = Mesh()
        msh.setxye(x,y,z,ele)
        msh.project_HTRS_WGS84()
        self.msh = msh
        self.xmin = msh.x.min()
        self.xmax = msh.x.max()
        self.ymin = msh.y.min()
        self.ymax = msh.y.max()
        #self.zmin = msh.z.min()
        #self.zmax = msh.z.max()
            
        
    def __init__(self, filename=None, name=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        # Stavit cu false i sam raditi prostornu interpolaciju
        self.return_block = False

        variable_aliases = {'dahv_u': 'x_sea_water_velocity',
                            'dahv_v': 'y_sea_water_velocity'}
        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                from glob import glob
                self._read_mesh(glob(filestr)[0])
                
                logging.info('Opening files with MFDataset')
                self.Dataset = MFDataset(filename)
            else:
                logging.info('Opening file with Dataset')
                self._read_mesh(filename)
                self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        # We are reading and using lon/lat arrays,
        # and not any projected coordinates
        self.proj4 =  '+proj=latlong'
        #self.proj4 = '+init=epsg:3765'

        self.variable_mapping = {}
        #self.variable_missing_value = {}
        
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, self.zname]:
                continue
            var = self.Dataset.variables[var_name]
            #attributes = var.ncattrs()
            if var_name == 'time':
                # Read and store time coverage (of this particular file)
                time = var[:]
                time_units = var.__dict__['units']
                self.times = num2date(time, time_units)
                self.start_time = self.times[0]
                self.end_time = self.times[-1]
                if len(self.times) > 1:
                    self.time_step = self.times[1] - self.times[0]
                else:
                    self.time_step = None
            if var_name in variable_aliases:
                self.variable_mapping[variable_aliases[var_name]] = str(var_name)


        self.variables = self.variable_mapping.keys()

        #self.xmin = self.lon.min()
        #self.xmax = self.lon.max()
        #self.ymin = self.lat.min()
        #self.ymax = self.lat.max()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()


    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        """Method which must be invoked by any reader (subclass).

        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.
            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).
            x, y: float or ndarrays; coordinates of requested points in the
                Spatial Reference System (SRS) of the reader (NB!!)
            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)
            block: bool, see return below

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neichbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """

#        print "get_variables: ", requested_variables, time, x, y, z, block
#        self.params = (requested_variables, time, x, y, z, block)        
#        return
#        (requested_variables, time, x, y, z, block) = self.params
    
    
        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        variables={}
        for par in requested_variables:
            var_name = self.variable_mapping[par]
            var = self.Dataset.variables[var_name]
            if var.ndim == 1:
                data = var[:]
            elif var.ndim == 2:
                data = var[indxTime,:]
            elif var.ndim == 3:
                data = var[indxTime,0,:]
            else:
                raise ValueError('Wrong dimension of %s: %i' %
                                 (var_name, var.ndim))
                
            res = self.msh.interp_scalar(data, x, y)
            variables[par] = res
            
        return variables
    
        

    def get_variables_block_true(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        """Method which must be invoked by any reader (subclass).

        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.
            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).
            x, y: float or ndarrays; coordinates of requested points in the
                Spatial Reference System (SRS) of the reader (NB!!)
            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)
            block: bool, see return below

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neichbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """

#        print "get_variables: ", requested_variables, time, x, y, z, block
#        self.params = (requested_variables, time, x, y, z, block)        
#        return
#        (requested_variables, time, x, y, z, block) = self.params
    
    
        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)


        # Finding a subset around the particles, so that
        # we do not interpolate more points than is needed.
        # Performance is quite dependent on the given buffer,
        # but it should not be made too small to make sure
        # particles are inside box
        buffer = .1  # degrees around given positions

        lonmin = x.min() - buffer
        lonmax = x.max() + buffer
        latmin = y.min() - buffer
        latmax = y.max() + buffer
        c = np.where((self.msh.x > lonmin) &
                     (self.msh.x < lonmax) &
                     (self.msh.y > latmin) &
                     (self.msh.y < latmax))[0]

        # Making a lon-lat grid onto which data is interpolated
        lonstep = .01  # hardcoded for now
        latstep = .01  # hardcoded for now
        lons = np.arange(lonmin, lonmax, lonstep)
        lats = np.arange(latmin, latmax, latstep)
        lonsm, latsm = np.meshgrid(lons, lats)

        # Initialising dictionary to contain data
        variables = {'x': lons, 'y': lats, 'z': z,
                     'time': nearestTime}

        # Reader coordinates of subset
        for par in requested_variables:
            var_name = self.variable_mapping[par]
            var = self.Dataset.variables[var_name]
            if var.ndim == 1:
                data = var[c]
            elif var.ndim == 2:
                data = var[indxTime,c]
            elif var.ndim == 3:
                data = var[indxTime,0,c]
            else:
                raise ValueError('Wrong dimension of %s: %i' %
                                 (var_name, var.ndim))

            if 'interpolator' not in locals():
                logging.debug('Making interpolator...')
                interpolator = LinearNDInterpolator((self.msh.y[c],
                                                     self.msh.x[c]),
                                                    data)
            else:
                # Re-use interpolator for other variables
                interpolator.values[:,0] = data

            variables[par] = interpolator(latsm, lonsm)

        #print variables
        return variables
