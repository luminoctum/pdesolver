#! /usr/bin/env python2.7
from pylab import *
from netCDF4 import *
from pycli.ode.gridop import *
import json
from collections import OrderedDict
json.encoder.FLOAT_REPR = lambda o: format(o, '.3E')

class ncfile:
    def __init__(self, fname, mode):
        self.file = Dataset(fname, mode, format = 'NETCDF3_64BIT')
    def add_dim(self, name, value, long_name, units): 
        dim = self.file.createDimension(name, len(value))
        dim = self.file.createVariable(name, 'f4', (name,))
        if len(value) > 0: dim[:] = value
        dim.long_name = long_name
        dim.units = units
        if name == 'time': self.time = dim
    def add_var(self, name, dim, value, long_name, units): 
        var = self.file.createVariable(name, 'f4', dim)
        if dim[0] == 'time':
            if len(dim) == 2: var[0,:] = value
            if len(dim) == 3: var[0,:,:] = value
            if len(dim) == 4: var[0,:,:,:] = value
        else:
            if len(dim) == 1: var[:] = value
            if len(dim) == 2: var[:,:] = value
            if len(dim) == 3: var[:,:,:] = value
        var.long_name = long_name
        var.units = units
    def add_atts(self, dicatts):
        for key, value in dicatts.iteritems():
            setattr(self.file, key, value)
    def __getitem__(self, var):
        return self.file.variables[var]
    def __del__(self):
        self.file.close()

class InitBase:
    var = {}
    def __init__(self, setups, varlist, output = 'dynamics.nc'):
        self.nx         = setups['num_grids_in_x']
        self.ny         = setups['num_grids_in_y']
        self.xlen       = setups['length_x']
        self.ylen       = setups['length_y']
        self.start      = setups['time_start']
        self.step       = setups['time_step']
        self.end        = setups['time_end']
        self.frame      = setups['steps_per_frame']
        setups["model_output"] = output;
        self.setups     = setups
        self.output     = output
        self.varlist    = varlist
    def initialize(self):
        self.set_variables()
        self.write_ncfile()
        self.write_control_file()
    def set_variables(self): pass
    def write_ncfile(self):
        file = ncfile(self.output, 'w')
        dims = {}
        for varlist in self.varlist:
            if varlist[1][0] == 'd':
                if varlist[0] == 'time':
                    file.add_dim(varlist[0], [], varlist[2], varlist[3])
                else:
                    file.add_dim(varlist[0], self.var[varlist[0]], varlist[2], varlist[3])
                dims[varlist[1][1]] = varlist[0]
            else:
                dim_name, dim_shape = (), ()
                for key in varlist[1]: 
                    dim_name = dim_name + (dims[key],)
                    if dims[key] != 'time':
                        dim_shape = dim_shape + (len(self.var[dims[key]]),)
                if varlist[0] in self.var.keys():
                    file.add_var(varlist[0], dim_name, self.var[varlist[0]].T, varlist[2], varlist[3])
                else:
                    file.add_var(varlist[0], dim_name, zeros(dim_shape), varlist[2], varlist[3])
        file.time[:] = 0
        file.add_atts(self.setups)
    def write_control_file(self):
        file = open('control.json', 'w')
        file.write(str(self))
        file.close()
    def __str__(self):
        return json.dumps(self.setups, indent = 4, separators = (', ', ' : '))
