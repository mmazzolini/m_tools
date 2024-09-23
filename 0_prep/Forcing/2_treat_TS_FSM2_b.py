# script for creating netcdf files ready for musa

import numpy as np
import pandas as pd
import xarray as xr
# system stuff
from glob import glob
import os
import sys
import socket
import datetime
#import matplotlib.pyplot as plt
#import pyproj
from multiprocessing import Process
from multiprocessing import Pool

# import custom modules
if socket.gethostname() == 'ra':
    basepath = '/home/marcomaz/snowdepth/'
else :    
    basepath = '/uio/hypatia/geofag-felles/projects/snowdepth/'
mytoolspath = basepath + 'code/tools/own_repos/snowdepth/tools'
sys.path.append(mytoolspath)

dempath = basepath+'data/Modeling/MuSA/DEM/Dem_full.nc'
topopath = basepath+'data/Modeling/TopoSCALE/'
forcingpath = basepath+'data/Modeling/MuSA/MuSA_input/Forcing_full/'

#define some functions

def ord_to_dt(ord):
    dt  = datetime.datetime.fromordinal(int(ord)) + datetime.timedelta(hours=int(ord%1*24)) - datetime.timedelta(days = 366)
    return dt


def get_indexes(xi,xf,yi,yf,ts):
    
    ei = np.where(ts.xg >= xi)[0][0]
    ef = np.where(ts.xg <= xf)[0][-1]
    
    ni = np.where(ts.yg >= yi)[0][-1]
    nf = np.where(ts.yg <= yf)[0][0]
    
    return ei,ef,ni,nf


#define resolution
res = 20

ts = xr.open_dataset(topopath+'TopoSCALE_Izas(2).nc')

dem = xr.open_dataset(dempath)

lookup_map = np.full((ts.xg.size,ts.yg.size),np.nan)

for i in np.arange(ts.xg.size):
    for j in np.arange(ts.yg.size):

        lookup_map[i,j] = ts.cn[(j)+(i)*(ts.yg.size)]

xi,xf,yi,yf = dem.easting.min().values,dem.easting.max().values+20,dem.northing.min().values-20,dem.northing.max().values

ei,ef,ni,nf = get_indexes(xi,xf,yi,yf,ts)

ts = ts.sel(easting=slice(ei,ef),northing=slice(nf,ni))

lookup_map = lookup_map[ei:ef,nf:ni]


datetimelist = []
for date in ts.t.data:
	datetimelist.append(ord_to_dt(date))

ts=ts.assign_coords({'t':datetimelist})

varlist = ['T','U','q','LW','SW','P','ps']

ts = ts.sel(t=slice("2019-09-01", "2020-08-31"))

times = ts.t

northings = ts.northing

eastings = ts.easting

x = ts.xg
y = ts.yg
dalist = []
for var in varlist:
    
    arr = np.full((ts.dims['easting'],ts.dims['northing'],ts.dims['t']),np.nan)
    
    for i in np.arange(x.size):
        for j in np.arange(y.size):
            
            c = lookup_map[i,j]-1
            arr[i,j,:] = ts[var][:,int(c)].astype(float)
    
    da = xr.DataArray(data = arr.transpose(), coords={'t':times,
    									  'northing':northings,
    									  'easting':eastings})
    da.name = ts[var].name
    da.attrs['long_name'] = ts[var].long_name
    dalist.append(da)

ds = xr.merge(dalist)

Tc = ds.T*0.01
q = ds.q*(10**-6)
ps= ds.ps*100

w=0.5*(1-np.sqrt(1-4*q))

e=0.5*ps*(-1+np.sqrt(1+4*w/0.622))
es = 610.94*np.exp(17.625*Tc/(243.04+Tc))
RH = e/es*100
RH.data[RH.data>100]=100


RH.name = 'RH'
RH = RH.assign_coords(coords={'t':times,
    						  'northing':northings,
    						  'easting':eastings})
RH.attrs['long_name'] = 'Relative Humidity [%]'
ds = xr.merge((ds,RH,ts.xg,ts.yg))
ds.attrs = {'name':'forcing for FSM musa'}

#I start saving
print('I start saving')

def save_day(label,group,forcingpath):
    group.to_netcdf(forcingpath+label.strftime("%Y-%m-%d")+'.nc')

pool = Pool(35)
pool.map(save_day,[(label,group,forcingpath) for label,group in ts.groupby('t.date')])

pool.close()
pool.join()

print('I saved')
    
"""
processes = [Process(target=save_day,args=(label,group,forcingpath,)) for label,group in ds.groupby('t.date')]

for process in processes:
    process.start()

for process in processes:
    process.join()

print('Done',flush=True)
"""

#how I did before
#for label,group in ds.groupby('t.date'):
#	group.to_netcdf(forcingpath+label.strftime("%Y-%m-%d")+'.nc')
