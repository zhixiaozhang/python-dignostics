# Python-AFWA offline diagnostic package for caculating enviromental conditions from WRFOUT
# Copy module_diag_functions.f90 and recalc_wrfout.py in the same working directory
# Run Linux command: "python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions" to add diagnostic functions to numpy.f2py
# Notes: afwa.diag_functions.diag_map and afwa.diag_functions.diag_row return the 2-D diagnostic map and row, respectively.
# When the model domain is large, diag_row is recommended to be used with proper parallelization.
# zhixiao.zhang@utah.edu, 2021/09/21

# from multiprocessing import Pool
import numpy as np
import xarray as xr
import glob, os, sys, time
import yaml
from netCDF4 import Dataset
from wrf import getvar
import numpy.f2py as f2py
import diag_functions as afwa
import dask
from dask.distributed import Client, LocalCluster

def add_vars(filename):
    """
    Calculates WRF diagnostics and adds to the WRF file.
    ----------
    filename: string
        Input WRF filename

    Returns
    ----------
    status: True
        Returns status = True if success.
    """

    ncfile = Dataset(filename, 'r+')
    p = getvar(ncfile, "p",units="Pa")
    tk = getvar(ncfile, "temp",units="K")
    rh = getvar(ncfile, "rh")
    hgt = getvar(ncfile, "z",units="m")
    ostat,cape,cin,lcl,lfc,el,lpl = afwa.diag_functions.diag_map(tk,rh,p,hgt,1,1)
    
    tmp=ncfile.createVariable('AFWA_CAPE_MU','float64',('south_north', 'west_east'))
    tmp[:]=cape
    tmp.long_name= 'Corrected Most Unstable CAPE'
    tmp.units='J kg^-1'
    tmp.description='Accumulated positive buoyant energy between LFC and EL'
    
    tmp=ncfile.createVariable('AFWA_CIN_MU','float64',('south_north', 'west_east'))
    tmp[:]=cin
    tmp.long_name= 'Corrected Most Unstable CIN'
    tmp.units='J kg^-1'
    tmp.description='Accumulated negative buoyant energy between LCL and LFC'
    
    tmp=ncfile.createVariable('AFWA_LCL','float64',('south_north', 'west_east'))
    tmp[:]=lcl
    tmp.long_name= 'Corrected LCL Height'
    tmp.units='m'
    tmp.description='Lifted condensation level (MSL) starting from the most unstable parcel'
    
    tmp=ncfile.createVariable('AFWA_ZLFC','float64',('south_north', 'west_east'))
    tmp[:]=lfc
    tmp.long_name= 'Corrected LFC Height'
    tmp.units='m'
    tmp.description='Level of Free Convection (MSL) starting from the most unstable parcel'
    
    tmp=ncfile.createVariable('AFWA_EL','float64',('south_north', 'west_east'))
    tmp[:]=el
    tmp.long_name= 'Corrected height of Equilibrium Level'
    tmp.units='m'
    tmp.description='Equilibrium level (MSL) lifted from the most unstable parcel'
    
    tmp=ncfile.createVariable('AFWA_LPL','float64',('south_north', 'west_east'))
    tmp[:]=lpl
    tmp.long_name= 'Corrected initial height of the most unstable parcel'
    tmp.units='m'
    tmp.description='Most unstable parcel initial level (MSL) featured by the max theta-e below 500 mb'
    
    ncfile.close()

    return True

def calc_diag(filename, outdir, outbasename):
    """
    Calculates WRF diagnostics and writes to a separate netCDF file.
    ----------
    filename: string
        Input WRF filename
    outdir: string
        Output file directory.
    outbasename: string
        Output file basename.

    Returns
    ----------
    status: True
        Returns status = True if success.
    """
    
    # Make output filename
    fname = os.path.basename(filename)
    fdate = fname.split('_')[2]
    ftime = fname.split('_')[3]
    outfilename = f'{outdir}{outbasename}{fdate}_{ftime}.nc'

    # Read data file
    ncfile = Dataset(filename, 'r')
    # Get coordinates
    times = getvar(ncfile, "times")
    XLAT = getvar(ncfile, "XLAT")
    XLONG = getvar(ncfile, "XLONG")
    # Get dimension names
    ydim, xdim = XLAT.dims

    # Get WRF variables
    p = getvar(ncfile, "p", units="Pa")
    tk = getvar(ncfile, "temp", units="K")
    rh = getvar(ncfile, "rh")
    hgt = getvar(ncfile, "z", units="m")
    # Call diagnostic function
    ostat,cape,cin,lcl,lfc,el,lpl = afwa.diag_functions.diag_map(tk, rh, p, hgt, 1, 1)

    # Replace default missing values with NaN
    mv = np.nanmin(cape)
    FillValue = np.nan
    cape[cape == mv] = FillValue
    cin[cin == mv] = FillValue
    lfc[lfc == mv] = FillValue
    el[el == mv] = FillValue

    # Remove attributes ('projection' in particular conflicts with Xarray)
    attrs_to_remove = ['FieldType', 'projection', 'MemoryOrder', 'stagger', 'coordinates', 'missing_value']
    for key in attrs_to_remove:
        XLONG.attrs.pop(key, None)
        XLAT.attrs.pop(key, None)

    # Define xarray dataset
    print(f'Writing output ...')
    varlist = {'XLAT': ([ydim, xdim], XLAT, XLAT.attrs), \
               'XLONG': ([ydim, xdim], XLONG, XLONG.attrs), \
               'AFWA_CAPE_MU': (['time', ydim, xdim], np.expand_dims(cape, axis=0)), \
               'AFWA_CIN_MU': (['time', ydim, xdim], np.expand_dims(cin, axis=0)), \
               'AFWA_LCL': (['time', ydim, xdim], np.expand_dims(lcl, axis=0)), \
               'AFWA_ZLFC': (['time', ydim, xdim], np.expand_dims(lfc, axis=0)), \
               'AFWA_EL': (['time', ydim, xdim], np.expand_dims(el, axis=0)), \
               'AFWA_LPL': (['time', ydim, xdim], np.expand_dims(lpl, axis=0)), \
              }
    coordlist = {'time': (['time'], times.expand_dims('time', axis=0)), \
                #  'latitude': (['latitude'], XLAT, XLAT.attrs), \
                #  'longitude': (['longitude'], XLONG, XLONG.attrs), \
                }
    attrlist = {'title': 'WRF AFWA offline diagnostics', \
                'created_on':time.ctime(time.time()), \
                }
    dsout = xr.Dataset(varlist, coords=coordlist, attrs=attrlist)

    # Assign variable attributes
    dsout['AFWA_CAPE_MU'].attrs['long_name'] = 'Corrected Most Unstable CAPE'
    dsout['AFWA_CAPE_MU'].attrs['units'] = 'J kg^-1'
    dsout['AFWA_CAPE_MU'].attrs['description'] = 'Accumulated positive buoyant energy between LFC and EL'

    dsout['AFWA_CIN_MU'].attrs['long_name'] = 'Corrected Most Unstable CIN'
    dsout['AFWA_CIN_MU'].attrs['units'] = 'J kg^-1'
    dsout['AFWA_CIN_MU'].attrs['description'] = 'Accumulated negative buoyant energy between LCL and LFC'

    dsout['AFWA_LCL'].attrs['long_name'] = 'Corrected LCL Height'
    dsout['AFWA_LCL'].attrs['units'] = 'm'
    dsout['AFWA_LCL'].attrs['description'] = 'Lifted condensation level (MSL) starting from the most unstable parcel'

    dsout['AFWA_ZLFC'].attrs['long_name'] = 'Corrected LFC Height'
    dsout['AFWA_ZLFC'].attrs['units'] = 'm'
    dsout['AFWA_ZLFC'].attrs['description'] = 'Level of Free Convection (MSL) starting from the most unstable parcel'

    dsout['AFWA_EL'].attrs['long_name'] = 'Corrected height of Equilibrium Level'
    dsout['AFWA_EL'].attrs['units'] = 'm'
    dsout['AFWA_EL'].attrs['description'] = 'Equilibrium level (MSL) lifted from the most unstable parcel'

    dsout['AFWA_LPL'].attrs['long_name'] = 'Corrected initial height of the most unstable parcel'
    dsout['AFWA_LPL'].attrs['units'] = 'm'
    dsout['AFWA_LPL'].attrs['description'] = 'Most unstable parcel initial level (MSL) featured by the max theta-e below 500 mb'

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    # Remove FillValue attribute on coordinates
    coord_dict = {'XLAT': {'zlib':True, 'dtype':'float32', '_FillValue':None},
                  'XLONG': {'zlib':True, 'dtype':'float32', '_FillValue':None}}
    encoding.update(coord_dict)
    # Write output netCDF file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', unlimited_dims='time', encoding=encoding)
    print('Output saved: ', outfilename)

    return True
        
def main():

    config_file = sys.argv[1]

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    indir = config['indir']
    outdir = config['outdir']
    file_basename = config['file_basename']
    outbasename = config['outbasename']
    run_parallel = config['run_parallel']
    n_workers = config['n_workers']
    threads_per_worker = config['threads_per_worker']

    # Find all WRF files
    folderlist = sorted(glob.glob(f'{indir}'))
    nfolders = len(folderlist)
    filelist=[]
    for ifolder in range(0, nfolders):
        tmplist = sorted(glob.glob(f'{folderlist[ifolder]}/{file_basename}'))
        filelist=np.append(filelist,tmplist)
    nfiles = len(filelist)
    
    results = []
    if run_parallel == 0:
        # Serial processing

        for ifile in range(0, nfiles):
            result = calc_diag(filelist[ifile], outdir, outbasename)
            results.append(result)

    elif run_parallel == 1:
        # Parallel processing

        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

        for ifile in range(0, nfiles):
            result = dask.delayed(calc_diag)(filelist[ifile], outdir, outbasename)
            results.append(result)

        # Trigger dask computation
        final_result = dask.compute(*results)

if __name__ == "__main__":
    main()

