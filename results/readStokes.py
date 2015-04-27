from netCDF4 import Dataset as nf

ff = nf('test.out', 'r')
res=ff.variables['stokes1'][:]