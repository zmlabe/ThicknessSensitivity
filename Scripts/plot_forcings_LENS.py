"""
Forcing files for SITperturb from LENS (1979-2005; SST, SIC, SIT)
(2060-2080; SIT)

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 17 July 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as c
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_var_LENS as LV
import read_SeaIceThick_LENS as lens

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
print '\n' '----Plotting forcing files - %s----' % titletime 

ensembles = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 2006
year2 = 2080
years = np.arange(year1,year2+1,1)

data = Dataset(directorydata + 'SST-SIT_lens_1976-2005.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sith = data.variables['ice_thick'][:,:,:]
data.close()

data = Dataset(directorydata + 'SST-SIT_lens_2051-2080.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sitf = data.variables['ice_thick'][:,:,:]
data.close()

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

sith[np.where(sith == 0)] = np.nan            
varh = (sith[0] + sith[1] + sith[-1])/3.

sitf[np.where(sitf == 0)] = np.nan            
varf = (sitf[0] + sitf[1] + sitf[-1])/3.

varh[np.where(varh == 2.)] = np.nan
varhtemp = varh.copy()
varhtemp[np.where(varhtemp > 0)] = 1.

varf = varf*varhtemp

varf[np.where(varf == 2.)] = np.nan

varc = varh.copy()
varc[np.where(varc >= 0)] = 2.

### Set limits for contours and colorbars
limsit = np.arange(0,7.1,0.1)
barlimsit = np.arange(0,8,1)
limdiff = np.arange(-3,3.1,0.5)
barlimdiff = np.arange(-3,4,1)

fig = plt.figure()
ax = plt.subplot(121)

m1 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')

var, lons_cyclic = addcyclic(varh, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m1(lon2d, lat2d)
   
m1.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)       
m1.drawcoastlines(color='darkgrey',linewidth=0.1)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
m1.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m1.contourf(x,y,var,limsit,extend='max')
cs1 = m1.contour(x,y,var,barlimsit,linewidths=0.1,colors='darkgrey',
                linestyles='-')

m1.fillcontinents(color='dimgray')

def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(30,cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N-15)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
        
cmap = colormapSIT()      
cs.set_cmap('cubehelix')

###########################################################################

ax = plt.subplot(122)

m2 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')

var, lons_cyclic = addcyclic(varc, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m2(lon2d, lat2d)
  
m2.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)        
m2.drawcoastlines(color='darkgrey',linewidth=0.1)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
m2.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m2.contourf(x,y,var,limsit,extend='max')
cs1 = m2.contour(x,y,var,barlimsit,linewidths=0.1,colors='darkgrey',
                linestyles='-')

cmap = colormapSIT()      
cs.set_cmap('cubehelix') 

cbar_ax = fig.add_axes([0.312,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

#m.fillcontinents(color='dimgray')

cbar.set_label(r'\textbf{Thickness [m]}',fontsize=11,color='dimgray')
cbar.set_ticks(barlimsit)
cbar.set_ticklabels(map(str,barlimsit)) 
cbar.ax.tick_params(axis='x', size=.01)

plt.text(-0.15,23,r'\textbf{HIT}',color='dimgray',fontsize=30)
plt.text(0.93,23,r'\textbf{CIT}',color='dimgray',fontsize=30)
#plt.text(1.05,21.8,r'\textbf{CONSTANT}',color='dimgray',fontsize=10)
plt.text(-0.75,13.5,r'\textbf{DJF}',color='dimgray',fontsize=30,
         rotation=90)

plt.savefig(directoryfigure + 'sit_comp.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
fig = plt.figure()
ax = plt.subplot(121)

m3 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')

var, lons_cyclic = addcyclic(varf, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m3(lon2d, lat2d)

m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)          
m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m3.drawcoastlines(color='darkgrey',linewidth=0.1)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
m3.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m3.contourf(x,y,var,limsit,extend='max')
cs1 = m3.contour(x,y,var,barlimsit,linewidths=0.1,colors='darkgrey',
                linestyles='-')

def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(30,cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N-15)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
           
cmap = colormapSIT()      
cs.set_cmap('cubehelix')

m3.fillcontinents(color='dimgray')

cbar = m3.colorbar(cs,location='bottom',pad = 0.2,extend='max',
                  drawedges=False)
ticks = barlimsit
labels = map(str,barlimsit)
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.set_label(r'\textbf{Thickness [m]}',fontsize=11,color='dimgray')
cbar.ax.tick_params(axis='x', size=.01)

###########################################################################
diffsit = varf - varh

ax = plt.subplot(122)

m4 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')

var, lons_cyclic = addcyclic(diffsit, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m4(lon2d, lat2d)
     
m4.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)     
m4.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m4.drawcoastlines(color='darkgrey',linewidth=0.1)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.1,color='k',fontsize=6)
m4.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m4.contourf(x,y,var,limdiff,extend='both')
cs1 = m4.contour(x,y,var,barlimdiff,linewidths=0.1,colors='darkgrey',
                linestyles='-')

cmap = ncm.cmap('amwg_blueyellowred')      
cs.set_cmap(cmap) 
cbar = m4.colorbar(cs,location='bottom',pad = 0.2,extend='max',
                  drawedges=False)
ticks = barlimdiff
labels = map(str,barlimdiff)
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.set_label(r'\textbf{Difference [m]}',fontsize=11,color='dimgray')
cbar.ax.tick_params(axis='x', size=.01)

m4.fillcontinents(color='dimgray')

plt.annotate(r'\textbf{FIT}',xy=(-0.82,1.1),
             xycoords='axes fraction',color='dimgray',fontsize=30,alpha=1) 
plt.annotate(r'\textbf{FIT--HIT}',xy=(0.2,1.1),
             xycoords='axes fraction',color='dimgray',fontsize=30,alpha=1) 
plt.annotate(r'\textbf{DJF}',xy=(-1.45,0.55),
             xycoords='axes fraction',color='dimgray',fontsize=30,alpha=1,
             rotation=90) 

plt.savefig(directoryfigure + 'sit_diff.png',dpi=300)

print 'Completed: Script done!'