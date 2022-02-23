# python-dignostics
Python-AFWA offline diagnostic package for caculating enviromental conditions using history output from Weather Research and Forecasting (WRF) model
  Example1: recalc_wrfout.py Perform calculation and write them into a seperate files
  Example2: calc_wrfout.py Perform calculation and attach them into input model output files

Compile fortran model and add it to numpy: 
  Copy module_diag_functions.f90 and recalc_wrfout.py calc_wrfout.py in the same working directory
  Run Linux command: "python -m numpy.f2py -c module_diag_functions.f90 -m diag_functions" to add diagnostic functions to numpy.f2py

Notes: 
  afwa.diag_functions.diag_map and afwa.diag_functions.diag_row return the 2-D diagnostic map and row, respectively.
  When the model domain is larger, subrountine diag_row is recommended to be used with proper parallelization.

Contact:
zhixiao.zhang@utah.edu, 2021/09/21
