"""
Plot precipitation comparisons between HIT and FIT experiments. These are 
sea ice thickness perturbation experiments using WACCM4.

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
print('\n' '----Plotting Z50 - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for precipitation data
lat,lon,time,lev,prh = MO.readExperi(directorydata,'P','HIT','surface')
lat,lon,time,lev,prf = MO.readExperi(directorydata,'P','FIT','surface')

### Separate per periods (ON,DJ,FM)
prh_on = np.nanmean(prh[:,9:11,:,:],axis=1)
prf_on = np.nanmean(prf[:,9:11,:,:],axis=1)

prh_dj,prf_dj = UT.calcDecJan(prh,prf,lat,lon,'surface',1)

prh_fm = np.nanmean(prh[:,1:3,:,:],axis=1)
prf_fm = np.nanmean(prf[:,1:3,:,:],axis=1)

### Calculate period differenceds
diff_on = np.nanmean((prf_on-prh_on),axis=0)
diff_dj = np.nanmean((prf_dj-prh_dj),axis=0)
diff_fm = np.nanmean((prf_fm-prh_fm),axis=0)
diff_onq = prf_on-np.nanmean(prh_on,axis=0)
diff_djq = prf_dj-np.nanmean(prh_dj,axis=0)
diff_fmq = prf_fm-np.nanmean(prh_fm,axis=0)

### Calculate significance    
stat_on,pvalue_on = UT.calc_indttest(prh_on,prf_on)
stat_dj,pvalue_dj = UT.calc_indttest(prh_dj,prf_dj)
stat_fm,pvalue_fm = UT.calc_indttest(prh_fm,prf_fm)

###########################################################################
###########################################################################
###########################################################################
### Plot precipitation
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-2,2.1,0.05)
barlim = np.arange(-2,3,1) 

### Begin plot
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
cs1 = m.contourf(x,y,pvalue_onq,colors='None',hatches=['....'],
             linewidths=0.4)

ax1.annotate(r'\textbf{ON}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('precip4_diff_19lev')            
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
cs1 = m.contourf(x,y,pvalue_djq,colors='None',hatches=['....'],
             linewidths=0.4)

ax2.annotate(r'\textbf{DJ}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('precip4_diff_19lev')            
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
cs1 = m.contourf(x,y,pvalue_fmq,colors='None',hatches=['....'],
             linewidths=0.4)

ax3.annotate(r'\textbf{FM}',
            xy=(0, 0),xytext=(0.35,1.05),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cmap = ncm.cmap('precip4_diff_19lev')            
cs.set_cmap(cmap)  

cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.01)

plt.savefig(directoryfigure + 'PR_diff_FIT-HIT.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
#limit = np.arange(-5,5.1,0.1)
#barlim = np.arange(-5,6,5) 
#
#
#for i in xrange(diff_fmq.shape[0]):
#    ax3 = plt.subplot(7,6,i+1)
#    
#    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
#                area_thresh=10000.)
#    
#    var, lons_cyclic = addcyclic(diff_fmq[i], lon)
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
#    cmap = ncm.cmap('precip4_diff_19lev')            
#    cs.set_cmap(cmap)  
#
#cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='max',extendfrac=0.07,drawedges=False)
#cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray')
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(map(str,barlim)) 
#cbar.ax.tick_params(axis='x', size=.01)
#
#plt.subplots_adjust(wspace=0.00)
#plt.subplots_adjust(hspace=0)
#
#plt.savefig(directoryfigure + 'pr_diffens_fm.png',dpi=300)
print('Completed: Script done!')