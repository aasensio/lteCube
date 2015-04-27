#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
# Convert string to lower case up to the first occurrence of a separator
def lower_to_sep(string, separator='='):
	line=string.partition(separator)
	string=str(line[0]).lower()+str(line[1])+str(line[2])
	return string

from configobj import ConfigObj
import sys
import os
from subprocess import call
from ipdb import set_trace as stop

if (len(sys.argv) < 2):
	print "You need to give the configuration file"
	exit()

print "Using configuration file = "+sys.argv[1]

# Transform all keys to lowercase to avoid problems with
# upper/lower case
f = open(sys.argv[1],'r')
input_lines = f.readlines()
f.close()
input_lower = ['']
for l in input_lines:
	input_lower.append(lower_to_sep(l)) # Convert keys to lowercase

config = ConfigObj(input_lower)

file = open('conf.input','w')

# Write general information
if (config['general']['verbose']):
	file.write("1\n")
else:
	file.write("0\n")
file.write("'"+config['model']['logtau500']+"'\n")
file.write("'"+config['model']['electron density']+"'\n")
file.write("'"+config['model']['temperature']+"'\n")
file.write("'"+config['model']['field strength']+"'\n")
file.write("'"+config['model']['field inclination']+"'\n")
file.write("'"+config['model']['field azimuth']+"'\n")
file.write("'"+config['model']['los velocity']+"'\n")
file.write(config['general']['nx']+'  '+config['general']['ny']+'  '+config['general']['ndepths']+'\n')
file.write("'"+config['general']['file with output results']+"'\n")
file.write(config['general']['initial pixel']+"\n")
file.write(config['general']['number of pixels in chunk']+"\n")
file.write(config['general']['number of spectral regions']+"\n")
n_regions = int(config['general']['number of spectral regions'])

for i in range(n_regions):
	n_lines = int(config['region '+str(i+1)]['number of lines'])

	file.write(config['region '+str(i+1)]['number of lines']+"\n")
	file.write(config['region '+str(i+1)]['number of wavelengths']+"\n")
	file.write(config['region '+str(i+1)]['first wavelength']+"\n")
	file.write(config['region '+str(i+1)]['wavelength step']+"\n")

# Check whether we are defining the spectral lines individually
# or we are using a file with a linelist
	if (config['region '+str(i+1)].has_key('linelist filename')):
		print "Reading linelist file {0}".format(config['region '+str(i+1)]['linelist filename'])
		f = open(config['region '+str(i+1)]['linelist filename'],'r')
		spectralLines = f.readlines()
		f.close()
		for l in spectralLines:
			file.write(l)
	else:		
		for j in range(n_lines):
			file.write(config['region '+str(i+1)]['line '+str(j+1)]['central wavelength']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['energy lower [cm^-1]']+ ' '+
				config['region '+str(i+1)]['line '+str(j+1)]['loggf']+ ' '+
				config['region '+str(i+1)]['line '+str(j+1)]['atomic number']+ ' '+
				config['region '+str(i+1)]['line '+str(j+1)]['ionization stage']+ ' '+
				config['region '+str(i+1)]['line '+str(j+1)]['atomic weight']+ ' '+
				config['region '+str(i+1)]['line '+str(j+1)]['sigma0']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['alpha']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['lower level']['j']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['lower level']['g']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['upper level']['j']+' '+
				config['region '+str(i+1)]['line '+str(j+1)]['upper level']['g']+'\n')
	
file.close()

# Run the code
try:	
#	call(['/usr/lib64/mpich/bin/mpirun','-np',sys.argv[2],'./lteCube', 'conf.input'])
	call(['mpiexec','-n',sys.argv[2],'./lteCube', 'conf.input'])	
except:
	pass
#	os.remove('conf.input')