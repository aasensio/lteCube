from netCDF4 import Dataset as nf

ff = nf('test.nc', 'r')
res=ff.variables['stokes1'][:]
