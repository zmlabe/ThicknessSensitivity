"""
Plot test files of the forcings

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 1 February 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm
import datetime

### Define directories
directorydata = '/surtsey/ypeings/'
directoryfigure = '/home/zlabe/Desktop/testseaice/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculate forcing file SIT constant - %s----' % titletime)

### Read in data 
data = Dataset(directorydata + 'SST-SIC-SIT_lens_2051-2080_polar.nc')
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
sit = data.variables['ice_thick'][:]
data.close()

lons,lats = np.meshgrid(lon,lat)

#### Test Data
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

for i in range(sit.shape[0]):
    fig = plt.figure()
    ax = plt.subplot(111)
    
    m = Basemap(projection='ortho',lon_0=300,lat_0=90,resolution='l')
             
    var = sit[i]
              
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='darkgrey',linewidth=0.3)
    
    cs = m.contourf(lons,lats,var,np.arange(0,5.1,0.1),latlon=True,extend='max')
    cs1 = m.contour(lons,lats,lats,np.arange(66.6,67.6,1),linewidths=1,colors='r',
                    linestyles='--',latlon=True)
                    
    cs.set_cmap('cubehelix')
    m.fillcontinents(color='dimgrey')
    
    cbar = plt.colorbar(cs,extend='both')    
    
    cbar.set_label(r'\textbf{SIT (m)}')  
    ticks = np.arange(0,6,1)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(list(map(str,ticks)))     
    
    plt.savefig(directoryfigure + 'polar_testplot_%s.png' % i,dpi=300)