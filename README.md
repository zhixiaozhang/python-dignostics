# Python-AFWA Diagnostic Package

Calculate enviromental conditions offline using the history output from the Weather Research and Forecasting (WRF) model

Example1: recalc_wrfout.py Calculate lifted parcel paramenters and save them in a seperate file

Example2: calc_wrfout.py Calculate lifted parcel paramenters and attach them to the input file

# Compile Fortran Module

Copy module_diag_functions.f90 and recalc_wrfout.py calc_wrfout.py in the same working directory

Run Linux command: "python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions" to add diagnostic functions to numpy.f2py

# Call Fortran Module in Python

afwa.diag_functions.diag_map computes 2-D maps of mucape, mucin, lcl, lfc, el, and lpl from input (z, y, x) temperature, RH, pressure, and ASL height

afwa.diag_functions.diag_row computes 1-D vector of mucape, mucin, lcl, lfc, el, and lpl from input (z, x) temperature, RH, pressure, and ASL height

afwa.diag_functions.diag_map or afwa.diag_functions.diag_row (tk,rh,p,hgt,1,1) calculates the enviroments with parcel lifted from the most unstable layer.

afwa.diag_functions.diag_map or afwa.diag_functions.diag_row (tk,rh,p,hgt,0,1) calculates the enviroments with parcel lifted from surface.

# Notes

afwa.diag_functions.diag_map and afwa.diag_functions.diag_row return the 2-D diagnostic map and row, respectively.

If the analysis domain is overly large, subrountine diag_row is recommended to be used with proper parallelization.

Code is modified from AFWA module in WRF model

# Contact

zhixiao.zhang@utah.edu, 2021/09/21

# Collaborators

Drs. Adam Varble and Zhe Feng
