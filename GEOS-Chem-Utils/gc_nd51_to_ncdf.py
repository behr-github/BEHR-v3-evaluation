#!/usr/bin/env python2
from __future__ import print_function

import argparse
from bpch import bpch
from datetime import datetime
import netCDF4 as ncdf
import numpy
import os
from warnings import warn

import pdb


def check_version(version, max_version, msg=None):
    major, minor = [int(x) for x in version.split('.')[:2]]
    max_major, max_minor = [int(x) for x in max_version.split('.')[:2]]
    if major > max_major or minor > max_minor:
        if msg is None:
            msg = 'Incompatible version: {} > {}'.format(version, max_version)
        warn(msg)


def copy_nc_dims(bpch_obj, nc_obj):
    # The dimensions in the BPCH pseudo-netCDF structure don't have much information in themselves, but they seem to
    # correspond to variables in the rest of the file. All we need for dimension information is their size. A size of 0
    # indicates unlimited to netCDF4.
    for dim_name, dim in bpch_obj.dimensions.iteritems():
        if dim.isunlimited():
            dim_size = 0
        elif dim_name == 'nv':
            # nv seems to be a special case of a dimension that is used in bounds (time_bounds, latitude_bounds, and
            # longitude_bounds). I'm unclear what it stands for, but it clearly represents the two edges of a grid cell
            # in time or space.
            dim_size = 2
        else:
            try:
                dim_size = bpch_obj.variables[dim_name].size
            except ValueError:
                print('Cannot find a variable corresponding to the dimension "{}". Will try skipping it.'.format(dim_name))

        nc_obj.createDimension(dim_name, size=dim_size)


def copy_nc_attrs(bpch_obj, nc_obj):
    """
    Copy attributes from the BPCH object to the netCDF object. This may be the root object, a variable, or any part of
    the hierarchy that implements the ncattrs() method on the BPCH side and the setncattr() method on the netCDF side.
    :param bpch_obj: the pseudo-netCDF object representing part of the BPCH file hierarchy
    :param nc_obj: the netCDF object representing part of the netCDF file hierarchy
    :return: none, modified nc_obj in place.
    """
    for attr in bpch_obj.ncattrs():
        nc_obj.setncattr(attr, getattr(bpch_obj, attr))


def copy_nc_vars(bpch_obj, nc_obj):
    for var_name, var in bpch_obj.variables.iteritems():
        nc_var = nc_obj.createVariable(var_name, var.dtype, dimensions=var.dimensions)
        copy_nc_attrs(var, nc_var)
        nc_var[:] = var[:]


def make_netcdf_file(bpch_file, save_dir, overwrite=True):
    if not os.path.isfile(bpch_file):
        raise IOError('BPCH file "{}" does not exist'.format(bpch_file))
    elif not os.path.isdir(save_dir):
        raise IOError('Save directory "{}" does not exist'.format(save_dir))

    ncdf_name = os.path.join(save_dir, os.path.basename(os.path.splitext(bpch_file)[0]) + '.nc')
    if not overwrite and os.path.isfile(ncdf_name):
        print('{} already exists, skipping'.format(ncdf_name))
        return


    # The BPCH objects are not compatible with the "with" keyword so revert to the try-finally method
    print('Transforming {} into {}'.format(bpch_file, ncdf_name))
    bpch_obj = bpch(bpch_file)
    nc_obj = ncdf.Dataset(ncdf_name, mode='w')
    try:
        copy_nc_dims(bpch_obj, nc_obj)
        copy_nc_attrs(bpch_obj, nc_obj)
        copy_nc_vars(bpch_obj, nc_obj)
        # Add an extra attribute indicating the original BPCH file
        nc_obj.setncattr('history', 'Translated from {} using {} in {} on {}'.format(
            bpch_file, os.path.basename(__file__), os.getcwd(), datetime.now()
        ))
    finally:
        bpch_obj.close()
        nc_obj.close()


def parse_args():
    parser = argparse.ArgumentParser(description='Convert GEOS-Chem ND51 binary punch files into netCDF files')
    parser.add_argument('-s', '--save-dir', default='.', help='Directory to save to (default is %(default)s)')
    parser.add_argument('-n', '--no-overwrite', action='store_true', help='Do not overwrite existing netCDF files')
    parser.add_argument('bpch_files', nargs='+', help='The BPCH files to turn into netCDF files')

    return parser.parse_args()


def main():
    check_version(numpy.__version__, '1.8.0rc1', 'The BPCH package may not be compatible with numpy versions > 1.8.0rc1')
    args = parse_args()
    for bpch_file in args.bpch_files:
        make_netcdf_file(bpch_file, args.save_dir, not args.no_overwrite)


if __name__ == '__main__':
    main()
