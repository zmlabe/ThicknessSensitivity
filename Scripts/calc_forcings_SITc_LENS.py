"""
Create forcing file of constant sea ice thickness at 2 m

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 16 August 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm
import datetime

### Define directories
directorydata = '/surtsey/zlabe/LENS/ForcingPerturb/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculate forcing file SIT constant - %s----' % titletime)

### Set all values to 2 m
### Used NCO by:
# ncap2 -s 'where(ice_thick>0) ice_thick=2;' SST-SIT_lens_CTL.nc test.nc
# ncap2 -s 'where(ice_thick>0) ice_thick=2;' SST-SIC-SIT_lens_2051-2080_FICT.nc  SST-SIC-SIT_lens_2051-2080_FIC.nc

### Read in data 
data = Dataset(directorydata + 'SST-SIC-SIT_lens_2051-2080_FIC.nc')
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
sit = data.variables['ice_thick'][:]
data.close()

lons,lats = np.meshgrid(lon,lat)

#### Test Data
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='ortho',lon_0=300,lat_0=90,resolution='l')
         
var = sit[2]
          
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=0.3)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)

cs = m.contourf(lons,lats,var,np.arange(0,3,0.1),latlon=True,extend='both')

m.fillcontinents(color='dimgrey')
                
cs.set_cmap('cubehelix')

cbar = plt.colorbar(cs,extend='both')    

cbar.set_label(r'\textbf{SIT (m)}')  
ticks = np.arange(0,5,1)
cbar.set_ticks(ticks)
cbar.set_ticklabels(list(map(str,ticks)))     

plt.savefig(directoryfigure + 'test_sitconstant.png',dpi=300)

print('Completed: Script done!')