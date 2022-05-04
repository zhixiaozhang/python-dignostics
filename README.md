# Python-AFWA Diagnostic Package

Calculate enviromental conditions offline using history output from Weather Research and Forecasting (WRF) model

Example1: recalc_wrfout.py Perform calculation and write them into a seperate file

Example2: calc_wrfout.py Perform calculation and attach them into the input file

# Compile Fortran Module: 

Copy module_diag_functions.f90 and recalc_wrfout.py calc_wrfout.py in the same working directory

Run Linux command: "python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions" to add diagnostic functions to numpy.f2py

afwa.diag_functions.diag_map computes 2-D maps of mucape, mucin, lcl, lfc, el, and lpl from input (z, y, x) temperature, RH, pressure, and ASL height.

afwa.diag_functions.diag_vector computes 1-D vector of mucape, mucin, lcl, lfc, el, and lpl from input (z, x) temperature, RH, pressure, and ASL height.

# Notes: 

afwa.diag_functions.diag_map and afwa.diag_functions.diag_row return the 2-D diagnostic map and row, respectively.

If the analysis domain is large, subrountine diag_row is recommended to be used with proper parallelization.

# Contact:

zhixiao.zhang@utah.edu, 2021/09/21
