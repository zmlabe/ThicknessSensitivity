"""
Plot Z30 comparisons between HIT and FIT experiments. These are 
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
year1 = 1960
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for 30 mb height data
lat,lon,time,lev,Z30h = MO.readExperi(directorydata,'Z30','HIT','surface')
lat,lon,time,lev,Z30f = MO.readExperi(directorydata,'Z30','FIT','surface')

### Separate per periods (ON,DJ,FM)
Z30h_on = np.nanmean(Z30h[:,9:10,:,:],axis=1)
Z30f_on = np.nanmean(Z30f[:,9:10,:,:],axis=1)

Z30h_dj,Z30f_dj = UT.calcDecJan(Z30h,Z30f,lat,lon,'surface',1)

Z30h_fm = np.nanmean(Z30h[:,1:2,:,:],axis=1)
Z30f_fm = np.nanmean(Z30f[:,1:2,:,:],axis=1)

### Calculate period differenceds
diff_on = np.nanmean((Z30f_on-Z30h_on),axis=0)
diff_dj = np.nanmean((Z30f_dj-Z30h_dj),axis=0)
diff_fm = np.nanmean((Z30f_fm-Z30h_fm),axis=0)
diff_onq = Z30f_on-np.nanmean(Z30h_on,axis=0)
diff_djq = Z30f_dj-np.nanmean(Z30h_dj,axis=0)
diff_fmq = Z30f_fm-np.nanmean(Z30h_fm,axis=0)

### Calculate significance    
stat_on,pvalue_on = UT.calc_indttest(Z30h_on,Z30f_on)
stat_dj,pvalue_dj = UT.calc_indttest(Z30h_dj,Z30f_dj)
stat_fm,pvalue_fm = UT.calc_indttest(Z30h_fm,Z30f_fm)

###########################################################################
###########################################################################
###########################################################################
### Plot 30 mb heights
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-100,100.1,1)
barlim = np.arange(-100,101,50) 

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
cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.01)

plt.savefig(directoryfigure + 'Z30_diff.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
for i in xrange(diff_fmq.shape[0]):
    ax3 = plt.subplot(7,6,i+1)
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(diff_fmq[i], lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
#    pvalue_fmq,lons_cyclic = addcyclic(pvalue_fm, lon)
#    pvalue_fmq,lons_cyclic = shiftgrid(180.,pvalue_fmq,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='dimgray',linewidth=0.2)
    parallels = np.arange(-90,90,45)
    meridians = np.arange(-180,180,60)
    #m.drawparallels(parallels,labels=[True,True,True,True],
    #                linewidth=0.6,color='dimgray',fontsize=6)
    #m.drawmeridians(meridians,labels=[True,True,True,True],
    #                linewidth=0.6,color='dimgray',fontsize=6)
    #m.drawlsmask(land_color='dimgray',ocean_color='mintcream')
    
    cs = m.contourf(x,y,var,limit,extend='both')
#    cs1 = ax3.scatter(x,y,pvalue_fmq,color='k',marker='.',alpha=0.5,
#                    edgecolor='k',linewidth=0.2)
    
    cmap = ncm.cmap('nrl_sirkes')            
    cs.set_cmap(cmap)  

cbar_ax = fig.add_axes([0.312,0.03,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0)

plt.savefig(directoryfigure + 'Z30_diffens_fm.png',dpi=300)
print 'Completed: Script done!'

