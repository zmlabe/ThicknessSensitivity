"""
Plot temperature comparisons between HIT and FIT experiments. These are 
sea ice thickness perturbation experiments using WACCM4.

Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
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
directoryfigure = '/home/zlabe/Desktop/TestPerturb/'
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
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for surface temperature data
lat,lon,time,lev,tash = MO.readExperi(directorydata,'T2M','HIT','surface')
lat,lon,time,lev,tasf = MO.readExperi(directorydata,'T2M','FIT','surface')

### Separate per periods (ON,DJ,FM)
tash_on = np.nanmean(tash[:,9:10,:,:],axis=1)
tasf_on = np.nanmean(tasf[:,9:10,:,:],axis=1)

tash_dj,tasf_dj = UT.calcDecJan(tash,tasf,lat,lon,'surface',1)

tash_fm = np.nanmean(tash[:,1:2,:,:],axis=1)
tasf_fm = np.nanmean(tasf[:,1:2,:,:],axis=1)

### Calculate period differenceds
diff_on = np.nanmean((tasf_on-tash_on),axis=0)
diff_dj = np.nanmean((tasf_dj-tash_dj),axis=0)
diff_fm = np.nanmean((tasf_fm-tash_fm),axis=0)
diff_onq = tasf_on-np.nanmean(tash_on,axis=0)
diff_djq = tasf_dj-np.nanmean(tash_dj,axis=0)
diff_fmq = tasf_fm-np.nanmean(tash_fm,axis=0)

### Calculate significance
stat_on,pvalue_on = UT.calc_indttest(tash_on,tasf_on)
stat_dj,pvalue_dj = UT.calc_indttest(tash_dj,tasf_dj)
stat_fm,pvalue_fm = UT.calc_indttest(tash_fm,tasf_fm)

###########################################################################
###########################################################################
###########################################################################
### Plot surface temperature
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-10,10.1,0.5)
barlim = np.arange(-10,11,5)

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
m.drawcoastlines(color='dimgray',linewidth=0.8)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var,limit,extend='both')
cs1 = ax1.scatter(x,y,pvalue_onq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax1.annotate(r'\textbf{ON}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)   

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
m.drawcoastlines(color='dimgray',linewidth=0.8)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var,limit,extend='both')
cs1 = ax2.scatter(x,y,pvalue_djq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax2.annotate(r'\textbf{DJ}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)  

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
m.drawcoastlines(color='dimgray',linewidth=0.8)
parallels = np.arange(-90,90,45)
meridians = np.arange(-180,180,60)
#m.drawparallels(parallels,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,True,True],
#                linewidth=0.6,color='dimgray',fontsize=6)
#m.drawlsmask(land_color='dimgray',ocean_color='mintcream')

cs = m.contourf(x,y,var,limit,extend='both')
cs1 = ax3.scatter(x,y,pvalue_fmq,color='k',marker='.',alpha=0.5,
                edgecolor='k',linewidth=0.2)

ax3.annotate(r'\textbf{FM}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)  

cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(wspace=0.01)

plt.savefig(directoryfigure + 'T2M_diff_FIT-HIT.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
#for i in xrange(diff_onq.shape[0]):
#    ax3 = plt.subplot(7,6,i+1)
#    
#    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
#                area_thresh=10000.)
#    
#    var, lons_cyclic = addcyclic(diff_onq[i], lon)
#    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
#    x, y = m(lon2d, lat2d)
#    
##    pvalue_fmq,lons_cyclic = addcyclic(pvalue_fm, lon)
##    pvalue_fmq,lons_cyclic = shiftgrid(180.,pvalue_fmq,lons_cyclic,start=False)
#              
#    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
#    m.drawcoastlines(color='dimgray',linewidth=0.2)
#    parallels = np.arange(-90,90,45)
#    meridians = np.arange(-180,180,60)
#    #m.drawparallels(parallels,labels=[True,True,True,True],
#    #                linewidth=0.6,color='dimgray',fontsize=6)
#    #m.drawmeridians(meridians,labels=[True,True,True,True],
#    #                linewidth=0.6,color='dimgray',fontsize=6)
#    #m.drawlsmask(land_color='dimgray',ocean_color='mintcream')
#    
#    cs = m.contourf(x,y,var,limit,extend='both')
##    cs1 = ax3.scatter(x,y,pvalue_fmq,color='k',marker='.',alpha=0.5,
##                    edgecolor='k',linewidth=0.2)
#    
#    cmap = ncm.cmap('NCV_blu_red')            
#    cs.set_cmap(cmap)  
#
#cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='max',extendfrac=0.07,drawedges=False)
#cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(map(str,barlim)) 
#cbar.ax.tick_params(axis='x', size=.01)
#
#plt.subplots_adjust(wspace=0.00)
#plt.subplots_adjust(hspace=0)
#
#plt.savefig(directoryfigure + 't2m_FIT-HIT.png',dpi=300)
print('Completed: Script done!')

