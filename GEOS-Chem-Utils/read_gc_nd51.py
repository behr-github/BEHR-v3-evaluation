# read_gc_nd51
# Python script to read a whole bunch of .bpch files output
# from GEOS-Chem into a single Matlab structure that will
# have a similar format as that produced by the Matlab routine
# "read_geos_output". These will still need a little cleanup
# using the Matlab function "cleanup_py_struct"
from __future__ import print_function
import os
import sys
import glob
from bpch import bpch
from datetime import date
import numpy as np
import scipy.io as sio
import math


##### INPUT #####

loaddir = raw_input('Enter the directory to load the files from: ')

if not os.path.isdir(loaddir):
    raise IOError('Given load directory does not exist')
    
file_prefix = raw_input('Enter the file prefix (part before the date). Uses ts_satellite by default: ')

if file_prefix is '':
    file_prefix = 'ts_satellite*.bpch'
elif file_prefix[-1] is '*':
    file_prefix += '.bpch'
elif file_prefix[-1] is not '*' and file_prefix[-5:] is not '.bpch':
    file_prefix += '*.bpch'

savedir = raw_input('Enter the directory to save the .mat file to. Same as load dir by default: ')

if savedir is '':
    savedir = loaddir
    
if not os.path.isdir(savedir):
    raise IOError('Given save directory does not exist')
    
savename = raw_input('Enter the name to give the save file, before the variable name: ')

if savename is '':
    savename = 'pyBPCHparse.mat'


#### MAIN LOOP ####
if os.getcwd() is not loaddir:
    os.chdir(loaddir)

F = glob.glob(file_prefix)
# Apparently globs do not necessarily return in alphabetical order
# which in this case is bad because alphabetical order translates to
# chronological order, which is rather important.
F.sort() 
nF = len(F)

# Number of days between Jan 1st, 1985 and Jan 1st, 0001
datedelta = date(1985,1,1) - date(1,1,1)

first_var = 1

tmpf = bpch(F[0])
gcvars = tmpf.variables.keys()
n = len(gcvars)
for i in xrange(n):
    if i % 3 is 0 or i==(n-1):
        print(gcvars[i])
    else:
        print('{0: <16}'.format(gcvars[i]), '\t',end='')

varname = raw_input('Enter the variables to read, separated by a space: ')
vars = varname.split()

for i in range(len(vars)):
    if vars[i] not in gcvars:
        print('{0} is not present in the file - enter the correct spelling or q to quit:'.format(vars[i]))
        while True:
            new_var = raw_input(">> ")
            if new_var is "q" or new_var is "Q":
                exit(0)
            elif new_var in gcvars:
                vars[i] = new_var
                break
            else:
                print('{0} is still not present in the file - enter the correct spelling or q to quit:'.format(new_var))



for varname in vars:
        i=0
        first_one = True
        matdict_array = []
        for fname in F:
            #Visual indicator of progress
            if i>0:
                # Delete the previous progress message to keep the console tidy.
                print('\b'*(len(prog_str)+1), end='')
            
            prog_str = 'Reading files: {0} of {1}'.format(i,len(F))
            print(prog_str, end='')
            sys.stdout.flush()

            f = bpch(fname)
            
            no2 = f.variables[varname]
            # This gets imported as time x level x lat x lon and we expect the reverse:
            no2 = no2.transpose([3,2,1,0])

            t = f.variables['time']
            dayssince1985 = math.floor(t[0]/24)
            datenum = datedelta.days + dayssince1985 + 367 # Matlab date numbers are from Jan 0000, which necessitates the + 367 correction

            if first_one:
                first_one = False
                sz_no2 = no2.shape
                dataBlock = np.empty((sz_no2[0],sz_no2[1],sz_no2[2],nF))
                tVec = np.empty((nF,1))
            
            dataBlock[:,:,:,i] = no2.squeeze()
            tVec[i] = datenum
            i+=1
            
        print('')

        dataUnit = no2.units
        fullName = no2.long_name
        fullCat = ''
        modelName = 'GEOS5_47L'
        lat = f.variables['latitude']
        latres = round((180.0 / lat.size)*2)/2
        lon = f.variables['longitude']
        lonres = round((360.0 / lon.size)*2)/2
        modelRes = [lonres,latres]
        dataScale = np.nan
        molMass = np.nan
        tEdge = np.empty((nF+1,1))
        tEdge[:] = np.nan

        matdict = {'dataBlock':dataBlock, 'dataUnit':dataUnit, 'fullName':fullName, 'fullCat':fullCat, 'tVec':tVec, 'modelName':modelName, 'modelRes':modelRes, 'dataScale':dataScale, 'molMass':molMass, 'tEdge':tEdge} 
        matdict_array.append(matdict)

        this_savename = "{0}-{1}_py.mat".format(savename, varname.replace("$_",""))
        print("Saving {0} in {1}".format(this_savename, savedir))
        sio.savemat(os.path.join(savedir,this_savename),{'SatData':matdict_array})

    
    
    
    
