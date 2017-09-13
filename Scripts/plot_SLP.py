"""
Plot SLP comparisons between HIT and FIT experiments. These are 
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
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting sea level pressure - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for SLP data
lat,lon,time,lev,slph = MO.readExperi(directorydata,'SLP','HIT','surface')
lat,lon,time,lev,slpf = MO.readExperi(directorydata,'SLP','FIT','surface')

### Separate per periods (ON,DJ,FM)
slph_on = np.nanmean(slph[:,9:11,:,:],axis=1)
slpf_on = np.nanmean(slpf[:,9:11,:,:],axis=1)

slph_dj,slpf_dj = UT.calcDecJan(slph,slpf,lat,lon,'surface',1)

slph_fm = np.nanmean(slph[:,1:3,:,:],axis=1)
slpf_fm = np.nanmean(slpf[:,1:3,:,:],axis=1)

### Calculate period differenceds
diff_on = np.nanmean((slpf_on-slph_on),axis=0)
diff_dj = np.nanmean((slpf_dj-slph_dj),axis=0)
diff_fm = np.nanmean((slpf_fm-slph_fm),axis=0)
diff_onq = slpf_on-np.nanmean(slph_on,axis=0)
diff_djq = slpf_dj-np.nanmean(slph_dj,axis=0)
diff_fmq = slpf_fm-np.nanmean(slph_fm,axis=0)
    
stat_on,pvalue_on = UT.calc_indttest(slph_on,slpf_on)
stat_dj,pvalue_dj = UT.calc_indttest(slph_dj,slpf_dj)
stat_fm,pvalue_fm = UT.calc_indttest(slph_fm,slpf_fm)

###########################################################################
###########################################################################
###########################################################################
### Plot sea level pressure data
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-6,6.1,0.5)
barlim = np.arange(-6,7,3)

### Begin figure
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

cmap = ncm.cmap('nrl_sirkes')            
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

cmap = ncm.cmap('nrl_sirkes')            
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

cmap = ncm.cmap('nrl_sirkes')            
cs.set_cmap(cmap)  

cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.01)

plt.savefig(directoryfigure + 'SLP_diff_FIT-HIT.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Set limits for contours and colorbars
#limit = np.arange(-15,16.1,1)
#barlim = np.arange(-15,16,5)
#
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
#    cmap = ncm.cmap('nrl_sirkes')            
#    cs.set_cmap(cmap)  
#
#cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='max',extendfrac=0.07,drawedges=False)
#cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(map(str,barlim)) 
#cbar.ax.tick_params(axis='hPa', size=.01)
#
#plt.subplots_adjust(wspace=0.00)
#plt.subplots_adjust(hspace=0)
#
#plt.savefig(directoryfigure + 'slp_diffens_on.png',dpi=300)
print('Completed: Script done!')

