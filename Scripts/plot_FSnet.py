"""
Plot net surface energy flux difference between HIT and FIT experiments. 
These are sea ice thickness perturbation experiments using WACCM4.

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting temperature - %s----' % titletime)

### Alott time series
year1 = 1960
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for energy flux terms
lat,lon,time,lev,fsh = MO.readExperi(directorydata,'RNET','HIT','surface')
lat,lon,time,lev,fsf = MO.readExperi(directorydata,'RNET','FIT','surface')

### Separate per periods (ON,DJ,FM)
fsh_on = np.nanmean(fsh[:,9:10,:,:],axis=1)
fsf_on = np.nanmean(fsf[:,9:10,:,:],axis=1)

fsh_dj,fsf_dj = UT.calcDecJan(fsh,fsf,lat,lon,'surface',1)

fsh_fm = np.nanmean(fsh[:,1:2,:,:],axis=1)
fsf_fm = np.nanmean(fsf[:,1:2,:,:],axis=1)

### Calculate period differenceds
diff_on = np.nanmean((fsf_on-fsh_on),axis=0)
diff_dj = np.nanmean((fsf_dj-fsh_dj),axis=0)
diff_fm = np.nanmean((fsf_fm-fsh_fm),axis=0)

### Calculate significance
stat_on,pvalue_on = UT.calc_indttest(fsh_on,fsf_on)
stat_dj,pvalue_dj = UT.calc_indttest(fsh_dj,fsf_dj)
stat_fm,pvalue_fm = UT.calc_indttest(fsh_fm,fsf_fm)

###########################################################################
###########################################################################
###########################################################################
### Plot net surface energy flux
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-50,50.1,1)
barlim = np.arange(-50,51,25)

fig = plt.figure()
ax1 = plt.subplot(131)

m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(diff_on, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvalue_onq,lons_cyclic = addcyclic(pvalue_on, lon)
pvalue_onq,lons_cyclic = shiftgrid(180.,pvalue_onq,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m.drawcoastlines(color='darkgray',linewidth=0.2)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var*-1,limit,extend='both')
cs1 = ax1.scatter(x,y,pvalue_onq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax1.annotate(r'\textbf{ON}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('nrl_sirkes')            
cs.set_cmap(cmap)   
m.fillcontinents(color='dimgray')

###########################################################################

ax2 = plt.subplot(132)

m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(diff_dj, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvalue_djq,lons_cyclic = addcyclic(pvalue_dj, lon)
pvalue_djq,lons_cyclic = shiftgrid(180.,pvalue_djq,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m.drawcoastlines(color='darkgray',linewidth=0.2)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var*-1,limit,extend='both')
cs1 = ax2.scatter(x,y,pvalue_djq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax2.annotate(r'\textbf{DJ}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('nrl_sirkes')            
cs.set_cmap(cmap)  
m.fillcontinents(color='dimgray')

###########################################################################

ax3 = plt.subplot(133)

m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(diff_fm, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvalue_fmq,lons_cyclic = addcyclic(pvalue_fm, lon)
pvalue_fmq,lons_cyclic = shiftgrid(180.,pvalue_fmq,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
m.drawcoastlines(color='darkgray',linewidth=0.2)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var*-1,limit,extend='both')
cs1 = ax3.scatter(x,y,pvalue_fmq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax3.annotate(r'\textbf{FM}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('nrl_sirkes')            
cs.set_cmap(cmap)  
m.fillcontinents(color='dimgray')

cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.01)

plt.savefig(directoryfigure + 'FSnet_diff.png',dpi=300)
print('Completed: Script done!')

