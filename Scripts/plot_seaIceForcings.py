"""
Plot sea ice forcing files for WACCM4 experiments for concentration
and thickness. Data derives from the LENS ensemble mean.

Notes
-----
    Reference : Kay et al. [2014] - CESM Large Ensemble Project (LENS)
    Author : Zachary Labe
    Date   : 25 October 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cmocean
import datetime

### Define directories
directorydata = '/home/zlabe/surt/LENS/ForcingPerturb/' 
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting forcing files - %s----' % titletime)

### Read in HIT forcing file
data = Dataset(directorydata + 'SST-SIT_lens_1976-2005.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sith = data.variables['ice_thick'][:,:,:]
sich = data.variables['ice_cov'][:,:,:]
data.close()

### Read in FIT forcing file
data = Dataset(directorydata + 'SST-SIT_lens_2051-2080.nc')
sitf = data.variables['ice_thick'][:,:,:]
data.close()

### Read in FIC forcing file
data = Dataset(directorydata + 'SST-SIC-SIT_lens_2051-2080_FIC.nc')
sicf = data.variables['ice_cov'][:,:,:]
data.close()

print('Completed: Data read!')

### Create 2d array of latitude and longitude
lons2,lats2 = np.meshgrid(lons,lats)

### Average over DJF 
sith[np.where(sith == 0)] = np.nan            
varh = (sith[0] + sith[1] + sith[-1])/3.
sitf[np.where(sitf == 0)] = np.nan            
varf = (sitf[0] + sitf[1] + sitf[-1])/3.
sich[np.where(sich == 0)] = np.nan            
varch = ((sich[0] + sich[1] + sich[-1])/3.)*100 # convert SIC to 1-100%
sicf[np.where(sicf == 0)] = np.nan            
varcf = ((sicf[0] + sicf[1] + sicf[-1])/3.)*100 # convert SIC to 1-100%

### Use land/Arctic mask
varh[np.where(varh == 2.)] = np.nan
varhtemp = varh.copy()
varhtemp[np.where(varhtemp > 0)] = 1.
varf = varf*varhtemp
varf[np.where(varf == 2.)] = np.nan

### Set limits for contours and colorbars
limsit = np.arange(0,5.1,0.1)
barlimsit = np.arange(0,6,1)
limdiff = np.arange(-100,0.1,2)
barlimdiff = np.arange(-100,1,50)

limsic = np.arange(0,101,2)
barlimsic = np.arange(0,101,50)
limsicdiff = np.arange(-100,0.1,2)
barlimdiffsic = np.arange(-100,1,50)

print('Completed: Data processed!')
###############################################################################
###############################################################################
###############################################################################
### Create subplot of sea ice thickness and concentration fields
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

print('Completed: Begin plotting!')
###############################################################################
ax = plt.subplot(231)

m1 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
m1.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)       
m1.drawcoastlines(color='darkgrey',linewidth=0.1)

cst = m1.contourf(lons2,lats2,varh,limsit,extend='max',latlon=True)

#ax.annotate(r'\textbf{FIT}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center') 
m1.fillcontinents(color='dimgray')  
cst.set_cmap('cubehelix')

###############################################################################
ax = plt.subplot(232)

m2 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
m2.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)        
m2.drawcoastlines(color='darkgrey',linewidth=0.1)

cstt = m2.contourf(lons2,lats2,varf,limsit,extend='max',latlon=True)

m2.fillcontinents(color='dimgray')
#ax.annotate(r'\textbf{HIT}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center')    
cstt.set_cmap('cubehelix') 

###############################################################################
diffsit = ((varf - varh)/varh)*100. # calculate % change

ax = plt.subplot(233)

m3 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')    
m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)     
m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m3.drawcoastlines(color='darkgrey',linewidth=0.1)

cstd = m3.contourf(lons2,lats2,diffsit,limdiff,latlon=True)
cs1 = m2.contour(lons2,lats2,lats2,np.arange(66.6,67.6,1),linewidths=1,colors='r',
                linestyles='--',latlon=True)

cmap = cmocean.cm.ice   
cstd.set_cmap(cmap)
#ax.annotate(r'\textbf{FIT-HIT}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center') 
m3.fillcontinents(color='dimgray')

###############################################################################
ax = plt.subplot(234)

m1 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')  
m1.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)       
m1.drawcoastlines(color='darkgrey',linewidth=0.1)

csc = m1.contourf(lons2,lats2,varch,limsic,latlon=True)

m1.fillcontinents(color='dimgray')
#ax.annotate(r'\textbf{FIC}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center')         
cmap = cmocean.cm.dense    
csc.set_cmap(cmap)

###########################################################################
ax = plt.subplot(235)

m2 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')

m2.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)        
m2.drawcoastlines(color='darkgrey',linewidth=0.1)

cscc = m2.contourf(lons2,lats2,varcf,limsic,latlon=True)

#ax.annotate(r'\textbf{HIT}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center') 
m2.fillcontinents(color='dimgray')
cmap = cmocean.cm.dense    
cscc.set_cmap(cmap)

###############################################################################
diffsic = ((varcf - varch)/varch)*100. # calculate % change

ax = plt.subplot(236)

m3 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)     
m3.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m3.drawcoastlines(color='darkgrey',linewidth=0.1)

cscd = m3.contourf(lons2,lats2,diffsic,limsicdiff,latlon=True)
cs1 = m2.contour(lons2,lats2,lats2,np.arange(66.6,67.6,1),linewidths=1,colors='r',
                linestyles='--',latlon=True)

cmap = cmocean.cm.ice   
cscd.set_cmap(cmap)
#ax.annotate(r'\textbf{FIC-HIT}',xy=(0.5,0.82),
#             xycoords='axes fraction',color='k',fontsize=14,alpha=1,
#             ha='center') 
m3.fillcontinents(color='dimgray')

### Add text boxes
perc = r'$\bf{\%}$' 
plt.annotate(r'\textbf{THICKNESS}',xy=(-2.5,1.92),
             xycoords='axes fraction',color='k',fontsize=13,alpha=1,
             rotation=90,ha='right') 
plt.annotate(r'\textbf{CONCENTRATION}',xy=(-2.5,0.93),
             xycoords='axes fraction',color='k',fontsize=13,alpha=1,
             rotation=90,ha='right') 
plt.annotate(r'\textbf{H}',xy=(-1.885,2.3),
             xycoords='axes fraction',color='k',fontsize=20,alpha=1,
             ha='center') 
plt.annotate(r'\textbf{F}',xy=(-0.7,2.3),
             xycoords='axes fraction',color='k',fontsize=20,alpha=1,
             ha='center') 
plt.annotate(r'\textbf{F--H}',xy=(0.5,2.3),
             xycoords='axes fraction',color='k',fontsize=20,alpha=1,
             ha='center') 

### Add thickness colorbar
cbar_ax = fig.add_axes([0.28,0.505,0.2,0.02])                
cbar = fig.colorbar(cst,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{SIT [m]}',fontsize=10,color='dimgray',
                         labelpad=2)
cbar.set_ticks(barlimsit)
cbar.set_ticklabels(list(map(str,barlimsit)))
cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar.outline.set_edgecolor('dimgray')

### Add concentration colorbar
cbar_ax = fig.add_axes([0.28,0.09,0.2,0.02])                
cbar = fig.colorbar(csc,cax=cbar_ax,orientation='horizontal',
                    drawedges=False)
cbar.set_label(r'\textbf{SIC [%s]}' % perc,fontsize=10,
                          color='dimgray',labelpad=2)
cbar.set_ticks(barlimsic)
cbar.set_ticklabels(list(map(str,barlimsic)))
cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar.outline.set_edgecolor('dimgray')

### Add percent change in SIT/SIC colrobar
cbar_ax = fig.add_axes([0.685,0.09,0.2,0.02])                
cbar = fig.colorbar(cscd,cax=cbar_ax,orientation='horizontal',
                    drawedges=False)
cbar.set_ticks(barlimdiffsic)
cbar.set_ticklabels(list(map(str,barlimdiffsic)))
cbar.set_label(r'\textbf{[$\Delta$%s]}' % perc,fontsize=10,
                          color='dimgray',labelpad=2)
cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar.outline.set_edgecolor('dimgray')

print('Completed: Subplot done!')

### Save figure
plt.savefig(directoryfigure + 'seaice_forcingfiles.png',dpi=900)

print('Completed: Script done!')