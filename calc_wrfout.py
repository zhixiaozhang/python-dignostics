# Run on linux command first: python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions 


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
    
    ncfile['AFWA_CAPE_MU'][:] = cape
    ncfile['AFWA_CIN_MU'][:] = cin
    ncfile['AFWA_LCL'][:] = lcl
    ncfile['AFWA_ZLFC'][:]=lfc
    ncfile['AFWA_EL'][:]=el
    ncfile['AFWA_LPL'][:]=lpl
    
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

