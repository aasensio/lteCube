from netCDF4 import Dataset as nf

ff = nf('test.nc', 'r')
wave = ff.variables['lambda1']
res = ff.variables['stokes1'][:]
