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

# import custom modules
if socket.gethostname() == 'ra':
    basepath = '/home/marcomaz/snowdepth/'
else :    
    basepath = '/uio/hypatia/geofag-felles/projects/snowdepth/'
mytoolspath = basepath + 'code/tools/own_repos/snowdepth/tools'
sys.path.append(mytoolspath)
wd = basepath+'marcomaz/Izas_op_2/'
#dempath = basepath+'data/MuSA_op/MuSA_input/DEM/DEM.nc'
topopath = basepath+'data/Modeling/TopoSCALE/'
dempath = basepath+'data/Modeling/MuSA/DEM/Dem_full.nc'


#define resolution
res = 20

ts = xr.open_dataset(topopath+'TopoSCALE_Izas(2).nc')

#code for making the DEM smaller

xi,xf,yi,yf = 705557.5,714700,4732077.5,4736300

#xi,xf,yi,yf = 709786,710000,4733668,4733902

def get_indexes(xi,xf,yi,yf,ts):
    ei = np.where(ts.xg >= xi)[0][0]
    ef = np.where(ts.xg <= xf)[0][-1]
    
    ni = np.where(ts.yg >= yi)[0][-1]
    nf = np.where(ts.yg <= yf)[0][0]
    
    return ei,ef,ni,nf

ei,ef,ni,nf = get_indexes(xi,xf,yi,yf,ts)

ts = ts.sel(easting=slice(ei,ef),northing=slice(nf,ni))


y = ts.northing

x = ts.easting

x = ts.xg
y = ts.yg

DEM = ts.Z.transpose()

DEM = DEM.assign_coords({'northing':y,'easting':x})

DEM=DEM.to_dataset()
DEM.to_netcdf(dempath)
