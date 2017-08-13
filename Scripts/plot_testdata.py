"""
plot test data
 
Notes
-----
    Source : WACCM4
    Author : Zachary Labe
    Date   : 28 July 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/surt/simu/test/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot test- %s----' % titletime 
          
### Read in data
data = Dataset(directorydata + 'FIT1_11.nc')
tasf1 = data.variables['TS'][0,:,:]
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
data.close()

data = Dataset(directorydata + 'HIT1_11.nc')
tash1 = data.variables['TS'][0,:,:]
data.close()

data = Dataset(directorydata + 'FIT2_11.nc')
tasf2 = data.variables['TS'][0,:,:]
data.close()

data = Dataset(directorydata + 'HIT2_11.nc')
tash2 = data.variables['TS'][0,:,:]
data.close()

tasff1 = tasf1 - 273.15
tashh1 = tash1 - 273.15

tasff2 = tasf2 - 273.15
tashh2 = tash2 - 273.15

diff = ((tasff1 + tasff2)/2.) - ((tashh1 + tashh2)/2.)

var = diff

### Plot
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-20,20.1,0.1)
barlim = np.arange(-20,21,10)

fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m.drawcoastlines(color='dimgray',linewidth=1)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var,limit,extend='both')

cmap = ncm.cmap('hotcold_18lev')            
cs.set_cmap(cmap)     
            
cbar = plt.colorbar(cs,orientation='vertical',
                    extend='both',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='y', size=.01)

plt.savefig(directoryfigure + 'test_diff_waccm_11ens.png',dpi=300)