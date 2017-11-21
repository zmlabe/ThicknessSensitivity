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
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
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

months = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',
           r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']

### Read in HIT forcing file
#data = Dataset(directorydata + 'SST-SIT_lens_1976-2005.nc')
#lons = data.variables['lon'][:]
#lats = data.variables['lat'][:]
#sith = data.variables['ice_thick'][:,:,:]
#sich = data.variables['ice_cov'][:,:,:]
#data.close()
#
#### Read in FIT forcing file
#data = Dataset(directorydata + 'SST-SIT_lens_2051-2080.nc')
#sitf = data.variables['ice_thick'][:,:,:]
#data.close()
#
#### Read in FIC forcing file
#data = Dataset(directorydata + 'SST-SIC-SIT_lens_2051-2080_FIC.nc')
#sicf = data.variables['ice_cov'][:,:,:]
#data.close()
#
#lons2,lats2 = np.meshgrid(lons,lats)
#
#print('Completed: Data read!')
#
#### Create 2d array of latitude and longitude
#lons2,lats2 = np.meshgrid(lons,lats)
#
#### Average over DJF 
#sith[np.where(sith == 0)] = np.nan            
#varh = sith
#sitf[np.where(sitf == 0)] = np.nan            
#varf = sitf
#
#sich[np.where(sich == 0)] = np.nan            
#varch = sich*100 # convert SIC to 1-100%
#sicf[np.where(sicf == 0)] = np.nan            
#varcf = sicf*100 # convert SIC to 1-100%
#
#### Use land/Arctic mask
#varh[np.where(varh == 2.)] = np.nan
#varhtemp = varh.copy()
#varhtemp[np.where(varhtemp > 0)] = 1.
#varf = varf*varhtemp
#varf[np.where(varf == 2.)] = np.nan
#
#### Calculate differences 
#diffsitq = varf - varh
#diffsicq = varcf - varch
#
#diffsit = np.append(diffsitq[9:],diffsitq[:3],axis=0)
#diffsic = np.append(diffsicq[9:],diffsicq[:3],axis=0)
#
#diffs = np.append(diffsit,diffsic,axis=0)
#
#print('Completed: Data processed!')
###############################################################################
###############################################################################
###############################################################################
### Create subplots of sea ice anomalies 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limsit = np.arange(-3,0.1,0.25)
barlimsit = np.arange(-3,1,3)

limsic = np.arange(-100,1,5)
barlimsic = np.arange(-100,1,100)

perc = r'$\bf{\%}$' 
    
fig = plt.figure()
for i in range(len(months)):
        
    var = diffs[i]

    ax1 = plt.subplot(2,6,i+1)
    m = Basemap(projection='npstere',boundinglat=51,lon_0=270,resolution='l',
                round =True,area_thresh=10000)
    
#    var, lons_cyclic = addcyclic(var, lons)
#    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
#    x, y = m(lon2d, lat2d)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    m.drawcoastlines(color='darkgrey',linewidth=0.1)
    m.fillcontinents(color='dimgray')
    
    if i < 6:
        limit = limsit
        barlim = barlimsit
        extendq = 'min'
    else:
        limit = limsic
        barlim = barlimsic
        extendq = 'max'
        
    cs = m.contourf(lons2,lats2,var,limit,extend=extendq,latlon=True)
    
    cmap = cmocean.cm.dense_r 
    cs.set_cmap(cmap)
    
    ### Add experiment text to subplot
    ax1.annotate(r'\textbf{%s}' % months[i],xy=(0,0),xytext=(0.865,0.90),
                 textcoords='axes fraction',color='k',fontsize=11,
                 rotation=320,ha='center',va='center')
    
    if i == 5:
        ### Add thickness colorbar
        cbar_ax = fig.add_axes([0.412,0.529,0.2,0.02])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                            extend='min',extendfrac=0.07,drawedges=False)
        cbar.set_label(r'\textbf{$\Delta$SIT [m]}',fontsize=9,color='dimgray',
                                 labelpad=-7)
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
        cbar.outline.set_edgecolor('dimgray')

    elif i == 11:
        ### Add concentration colorbar
        cbar_ax = fig.add_axes([0.412,0.2315,0.2,0.02])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                            drawedges=False)
        cbar.set_label(r'\textbf{$\Delta$SIC [%s]}' % perc,fontsize=9,
                                  color='dimgray',labelpad=-7)
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
        cbar.outline.set_edgecolor('dimgray')
    
plt.subplots_adjust(hspace=-0.35)
plt.subplots_adjust(wspace=0)
        
plt.savefig(directoryfigure + 'monthly_seaiceanoms.png',dpi=300)