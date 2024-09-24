# PROCESS IZAS ACQUISITIONS FOR EVERY BEAM.
# 1) SELECT THE DATA
# 2) COREGISTER THE DEM
# 3) SAVE DATA WITH NEW DEM SAMPLED

# import packages

import pyarrow.feather as feather  # efficient data IO
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
#import pickle
import datetime
#import h5py
# import scipy.signal # ?
# work with geodata and rasters
import geopandas as gpd  # load shapefiles -- does not work on server
import pandas as pd
from shapely.geometry import MultiPoint, Polygon, LineString
from shapely.geometry import mapping  
from shapely import wkt
#import xarray as xr

#import pyproj
from pyproj import Transformer
#projstr= "+proj=utm +zone=33 +ellips=WGS84 +datum=WGS84 +units=m +no_defs"
#myProj = pyproj.Proj(projstr)
import h5py
import icepyx as ipx
import xdem
from xdem import coreg
import geoutils as gu

# system stuff
from glob import glob
import os
import sys
import socket
# ?? check whether that's needed / in the rigth place
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

# import custom modules
if socket.gethostname() == 'ra':
    basepath = '/home/marcomaz/snowdepth/'
else :    
    basepath = '/uio/hypatia/geofag-felles/projects/snowdepth/'
mytoolspath = basepath + 'code/tools/own_repos/snowdepth/tools'


sys.path.append(mytoolspath)
#import PointExtract
import QLtools
# replaces gda_lib_DT for all gdf functions - except for those requiring rasterio
import IC2tools#import PointExtract
# replaces gda_lib_DT for all gdf functions - except for those requiring rasterio

import wquantiles as wq

def wmedianf(group):
    return wq.quantile(group['dhcoreg'],group['yapc_score'],0.5)
    
# These are only if available with icesat2 / IC2 /xdem-dev environment
try:
    import contextily as ctx
    import getIC2data
    import rasterio
    # from rasterio.merge import merge # merge three dems
    import rasterio.plot
    from rasterio.mask import mask
    from rasterio import Affine  # or from affine import Affine
    import gda_lib_DT  # only used for the rasterio DEM extraction
except ImportError as e:
    print('Not all modules can be loaded in this environment:')
    print(os.environ['CONDA_DEFAULT_ENV'])
    print(e)
    pass


# enable kml support  
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'
supported_drivers['KML'] = 'rw'  # native gpd has no functionality to read kml
supported_drivers['KMZ'] = 'rw'  # native gpd has no functionality to read kml


wd = basepath+'marcomaz/Izas_op_2/'
forcing_path =  basepath+'data/Modeling/MuSA_dist/MuSA_input/Forcing/2020-01-01.nc'
datapath = wd+'processed_data/'
dempath = basepath+'data/DEM/'
figurefolder = wd+'/data_analysis/figures/'
fn3 = datapath+'ATL03_'

#y0=xr.open_dataset(forcing_path)['yg'].min().values
y0 = 4733697.5

wgs84 = 'EPSG:4326'
utm30n = 'EPSG:25830'

AOIname = 'Izas_Spain_sliderule_2m'
geoidfile = basepath+'/data/DEM/Datum_data/es_ign_egm08-rednap_30N.tif'
DEMpath = basepath+'data/DEM/SPAIN/IZAS/MDT02-ETRS89-HU30-0145-3-COB2_basin.tif'

slidepath = basepath+'data/ICESat-2/sliderule/ATL03/ATL03_yapc.feather'

coords_S = [-0.53, 42.659, -0.17, 42.84]

#plot settings
# 1ms corresponds to 7m, i.e. 10 bursts (of 0-8 photons each). That's usually OK for mountains (no veg)

windowsize = '2.910056ms'
windowlength = '21m'
ywind = 20
snowlims = [-2, 5]  # ylim for the snow depth transects
icptc = [[0.8, 0.5, 0.2]]  # icesat points
icmedc = [[0.6, 0.2, 0.1]]  # icesat median
icstdc = [[0.8, 0.8, 0.6]]  # icesat std fill




# open DEM & geoid

geoid = xdem.DEM(geoidfile)
DEM = xdem.DEM(DEMpath)
#slope = xdem.terrain.slope(DEM)
#aspect = xdem.terrain.aspect(DEM)

with open(wd + 'DATA/sl3.feather', 'rb') as f:
    sl = pd.DataFrame(feather.read_feather(f))
    
dates = sl.date.unique()
#dates = [20190204]

strong_beams = [10,20,30]

print(f'Iteration through the following dates: {dates}')
#iterate trough dates and strong beams
for date in dates:
    for beam in strong_beams:
        S03df = sl[((sl.date == date) & (sl.pb == beam))].copy()
        S_goodforcoreg= ((S03df.slope < 50) & (abs(S03df.dhnocoreg)<15))
        
        yapc_thresh = S03df[S_goodforcoreg].yapc_score.quantile(0.4)
        
        S_goodforcoreg = (S_goodforcoreg) & (S03df.yapc_score > yapc_thresh)
        
        if S_goodforcoreg.sum() > 100:
        
            print(date)
            
            #create folder:            
            wd_date = wd+str(date)+'_'+str(beam)+'/'
            if not os.path.exists(wd_date):
                os.makedirs(wd_date)
            fig,ax=plt.subplots()
            plt.hist(S03df[S_goodforcoreg].dhnocoreg, 100)
            ax.set_title('dh for coreg basis - '+AOIname)
            plt.savefig(wd_date+'dh_coregbasis_hist.png')
            
            
            C = S_goodforcoreg.sum()# check the histogram - done with abs(dh)<150
            print(f'for date {str(date)}, coregistering with {C} points')
            
            S03_coregpts = S03df[S_goodforcoreg.values].rename(columns={'x': 'lon', 'y': 'lat'})
            
            S03_coregpts['z'] = S03_coregpts['height']-S03_coregpts['geoid']
            
            try:
                nuth_kaab_S = coreg.NuthKaab(return_for_plot=True, weighting=True)
                nuth_kaab_S.fit_pts(S03_coregpts, DEM, verbose=True)
                nuth_kaab_S._plot_iterations(min_count=10, file_name=wd_date+'parameter_plot_')
            except:
                print("something went wrong in coregistration")
            
            else:
                IC2tools.savecoreginfo(wd_date, nuth_kaab_S)
                coregcoeff_S = nuth_kaab_S.to_matrix()
                
                DEM_c = DEM.copy()
                DEM_c.shift(-coregcoeff_S[0, 3], -coregcoeff_S[1, 3])
                DEM_c.save(wd_date+'coreg_2_dem.tif')        
    
                S03df['DEMcoreg'] = DEM_c.interp_points(np.array((S03df['x'].values, S03df['y'].values)).T, input_latlon=False, order=3)
                S03df['dhcoreg'] = S03df.height-S03df.geoid-S03df.DEMcoreg
                
                slope = xdem.terrain.slope(DEM_c)
                aspect= xdem.terrain.aspect(DEM_c)
    
                S03df['slope'] = slope.interp_points(np.vstack([np.array(S03df['x'].values), np.array(S03df['y'].values) ]).T, input_latlon=False, order=1)                                        
                S03df['aspect2'] = aspect.interp_points(np.vstack([np.array(S03df['x'].values), np.array(S03df['y'].values) ]).T, input_latlon=False, order=1)
    
                #SCATTER DH vs ASPECT
                fig, ax = plt.subplots()
                plt.scatter(S03df[S_goodforcoreg].aspect2,S03df[S_goodforcoreg].dhcoreg,
                            s=np.exp(S03df[S_goodforcoreg].yapc_score/100)/16,alpha=0.1)            
                ax.set_title('dh vs aspect coreg '+str(date))
                ax.set_xlim(0,360)
                ax.set_ylim(-4,8)
                plt.savefig(wd_date+'dh_vs_aspect_scatter.png')            
    
    
                #compute local stats
                
                S_goodforstats = (abs(S03df.dhcoreg) < 10)
                yapc_thresh = S03df[S_goodforstats].yapc_score.quantile(0.25)  
                S_goodforstats = (S_goodforstats) & (S03df.yapc_score > yapc_thresh)
                
                # resampled to window size.
                S03df['segment'] = ((S03df.y-y0-ywind/2)//ywind)
                y2= S03df[S_goodforstats].groupby(by='segment').segment.mean()*ywind+ywind/2+y0
                y = S03df[S_goodforstats].groupby(by='segment').y.mean()
                x = S03df[S_goodforstats].groupby(by='segment').x.mean()
                z = S03df[S_goodforstats].groupby(by='segment').DEMcoreg.median()-S03df[S_goodforstats].groupby(by='segment').geoid.median()
                
                S03df['time'] = S03df.index.values
                time = S03df[S_goodforstats].groupby(by='segment').time.mean()
                n = S03df[S_goodforstats].groupby(by='segment').dhcoreg.count()
    
                median = S03df[S_goodforstats].groupby(by='segment').apply(wmedianf)
                yapc_score = S03df[S_goodforstats].groupby(by='segment').yapc_score.mean()
                std = S03df[S_goodforstats].groupby(by='segment').dhcoreg.std()
                slope = S03df[S_goodforstats].groupby(by='segment').slope.mean()
                aspect = S03df[S_goodforstats].groupby(by='segment').aspect2.mean()
                
                
                #ICESAT-2 vs DEM and SNOWDEPTH TRANSECT
    
                miny = z.min()-1
                maxy = S03df[S_goodforstats].height.max()+2
                
                fig, ax = plt.subplots(2, 1, figsize=(15, 15), sharex=True)
                # vs DEM
                ax[0].scatter(S03df[S_goodforstats].y/1000, S03df[S_goodforstats].height-S03df[S_goodforstats]['geoid']
                              , s= np.exp(S03df[S_goodforstats].yapc_score/50)/16, c=icptc, label='photon height')
        
                ax[0].plot(S03df[S_goodforstats].y/1000, S03df[S_goodforstats]['DEMcoreg'], 'k.-',
                           markersize=0.1, lw=0.5, label='DEMcoreg')
                
                ax[0].set_ylim([miny, maxy])
                ax[0].set_title('ICESAT-2 strong beam vs DEM')
                ax[0].set_ylabel('Elevation masl [m]')
    
                #SNOWDEPTH
                
                ax[1].scatter(S03df[S_goodforstats].y/1000, S03df[S_goodforstats].dhcoreg,
                              s=0.15, c=icptc, label='ICESat-2 photons')
                ax[1].fill_between(y/1000, median-std,
                                   median+std, alpha=0.5, edgecolor=None,
                                   facecolor=icstdc, label='ICESat-2 stdev range', zorder=0)
                
                ax[1].plot(y/1000, median,'r.-',
                           markersize=1, lw=1, label=f"ICESat-2 median ({windowlength})")
     
                
                
                ax[1].axhline(0, lw=1, c='grey')  # the ground
                ax[1].set_xlabel('Northing utm30n [km]')
                ax[1].set_ylim(snowlims)
                # ax[1].set_xlim(lats)
                ax[1].legend()
                ax[1].set_ylabel('Snow depth [m]')
                ax[1].set_title('SNOW DEPTHs')
                fig.suptitle(f"Snow depths on {date}, {beam}",
                             fontsize='large', y=0.92)
                
                plt.savefig(wd_date+f'DEM&SD_{windowlength}2.png')
    
                #SCATTER DH vs SLOPE            
                fig, ax = plt.subplots()
                plt.scatter(slope,median,s=np.exp(yapc_score/50)/8,alpha=0.3)
                ax.set_title('dh vs slope - '+AOIname)
                ax.set_ylim(-1,5)
                plt.savefig(wd_date+'dh_vs_slope.png')
                
                #SCATTER DH vs ASPECT
                fig, ax = plt.subplots()
                plt.scatter(aspect,median,s=np.exp(yapc_score/25)/120,alpha=0.3)            
                ax.set_title('dh vs aspect coreg '+str(date))
                ax.set_xlim(0,360)
                ax.set_ylim(-1,5)
                plt.savefig(wd_date+'dh_vs_aspect_scatter_agg.png')
                
    
    
                
                plt.close('all')
                
                with open(wd_date + 'stats_ph.feather', 'wb') as f:
                        feather.write_feather(S03df[S_goodforstats], f)
                        
                #save local stats:
                stats = pd.DataFrame()
                
                stats['median_dh'] = median
                stats['y'] = y
                stats['y2'] = y2
                stats['x'] = x
                stats['z'] = z
                stats['yapc_score'] = yapc_score
                stats['std'] = std
                stats['date'] = date
                stats['beam'] = beam
                stats['n'] = n
                stats['slope'] = slope
                stats['aspect'] = aspect
                stats['time'] = time
                stats = stats.dropna()
                
                print(stats.shape)
                with open(wd_date + 'data.feather', 'wb') as file:
                        feather.write_feather(stats, file)
                
                stats.to_csv(wd_date+'data.csv')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
