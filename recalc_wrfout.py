# Python-AFWA offline diagnostic package for caculating enviromental conditions from WRFOUT
# Copy module_diag_functions.f90 and recalc_wrfout.py in the same working directory
# Run Linux command: python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions to add diagnostic functions to numpy.f2py
# Notes: afwa.diag_functions.diag_map and afwa.diag_functions.diag_row return the 2-D diagnostic map and row, respectively.
# When the model domain is large, diag_row is recommended to be used.
# zhixiao.zhang@utah.edu, 2021/09/21

from multiprocessing import Pool
import numpy as np
import glob
from netCDF4 import Dataset
from wrf import getvar
import numpy.f2py as f2py
import diag_functions as afwa

def add_vars(filepath):
    ncfile = Dataset(filepath,'r+')
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
        
def main():
    # set up number of workers for running in parallel
    n_workers = 32
    run_parallel = 1

    # paths
    indir = '/uufs/chpc.utah.edu/common/home/zipser-group2/zhixiao/wrf_output'
    folder_basename = '20*'
    file_basename = 'wrfout_d01_*00'

    # Find all WRF files
    folderlist = sorted(glob.glob(f'{indir}/{folder_basename}'))
    nfolders = len(folderlist)
    filelist=[]
    for ifolder in range(0, nfolders):
        tmplist = sorted(glob.glob(f'{folderlist[ifolder]}/{file_basename}'))
        filelist=np.append(filelist,tmplist)
    nfiles = len(filelist)
    
    if run_parallel == 0:
        for ifile in range(0, nfiles):
            add_vars(filelist[ifile])
    elif run_parallel == 1:
        pool = Pool(n_workers)
        pool.starmap(add_vars, zip(filelist))
        pool.close()

if __name__ == "__main__":
    main()

