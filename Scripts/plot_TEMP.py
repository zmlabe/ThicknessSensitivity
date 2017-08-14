"""
Plot zonal wind comparisons between HIT and FIT experiments. These are 
sea ice thickness perturbation experiments using WACCM4.

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
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

### Call function for surface temperature data
#lat,lon,time,lev,th = MO.readExperi(directorydata,'TEMP','HIT','profile')
#lat,lon,time,lev,tf = MO.readExperi(directorydata,'TEMP','FIT','profile')
#
#### Separate per periods (ON,DJ,FM)
#th_on = np.nanmean(th[:,9:10,:,:,:],axis=1)
#tf_on = np.nanmean(tf[:,9:10,:,:,:],axis=1)
#
#th_dj,tf_dj = UT.calcDecJan(th,tf,lat,lon,'profile',lev.shape[0])
#
#th_fm = np.nanmean(th[:,1:2,:,:,:],axis=1)
#tf_fm = np.nanmean(tf[:,1:2,:,:,:],axis=1)
#
##### Calculate period differenceds
#diff_on = np.nanmean((tf_on-th_on),axis=0)
#diff_dj = np.nanmean((tf_dj-th_dj),axis=0)
#diff_fm = np.nanmean((tf_fm-th_fm),axis=0)
#
##### Calculate significance
#stat_on,pvalue_on = UT.calc_indttest(np.nanmean(th_on,axis=3),
#                                     np.nanmean(tf_on,axis=3))
#stat_dj,pvalue_dj = UT.calc_indttest(np.nanmean(th_dj,axis=3),
#                                     np.nanmean(tf_dj,axis=3))
#stat_fm,pvalue_fm = UT.calc_indttest(np.nanmean(th_fm,axis=3),
#                                     np.nanmean(tf_fm,axis=3))
#
#### Calculate zonal mean
#zdiff_on = np.nanmean((diff_on),axis=2)
#zdiff_dj = np.nanmean((diff_dj),axis=2)
#zdiff_fm = np.nanmean((diff_fm),axis=2)

############################################################################
############################################################################
############################################################################
##### Plot surface temperature
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-4,4.1,0.1)
barlim = np.arange(-4,5,1)
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
ax1 = plt.subplot(131)

ax1.spines['top'].set_color('darkgrey')
ax1.spines['right'].set_color('darkgrey')
ax1.spines['bottom'].set_color('darkgrey')
ax1.spines['left'].set_color('darkgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='darkgrey')
ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='darkgrey')    
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')


cs = plt.contourf(lat,lev,zdiff_on,limit,extend='both')
cs2 = plt.contour(lat,lev,zdiff_on,np.arange(0,1,1),
                  linewidths=0.6,colors='dimgrey')
cs1 = plt.scatter(latq,levq,pvalue_on,color='k',marker='.',alpha=0.7,
                edgecolor='k',linewidth=0.5)

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.xlim([0,90])
plt.ylim([1000,10])
plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap) 

ax1.annotate(r'\textbf{ON}',
            xy=(0, 0),xytext=(0.34,1.02),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

###########################################################################
ax2 = plt.subplot(132)

ax2.spines['top'].set_color('darkgrey')
ax2.spines['right'].set_color('darkgrey')
ax2.spines['bottom'].set_color('darkgrey')
ax2.spines['left'].set_color('darkgrey')
ax2.spines['left'].set_linewidth(2)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['right'].set_linewidth(2)
ax2.spines['top'].set_linewidth(2)
ax2.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='darkgrey')
ax2.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='darkgrey')    
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')

cs = plt.contourf(lat,lev,zdiff_dj,limit,extend='both')
cs2 = plt.contour(lat,lev,zdiff_dj,np.arange(0,1,1),
                  linewidths=0.6,colors='dimgrey')
cs1 = plt.scatter(latq,levq,pvalue_dj,color='k',marker='.',alpha=0.7,
                edgecolor='k',linewidth=0.5)

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.xlim([0,90])
plt.ylim([1000,10])
plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap) 

ax2.annotate(r'\textbf{DJ}',
            xy=(0, 0),xytext=(0.35,1.02),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

###########################################################################
ax3 = plt.subplot(133)

ax3.spines['top'].set_color('darkgrey')
ax3.spines['right'].set_color('darkgrey')
ax3.spines['bottom'].set_color('darkgrey')
ax3.spines['left'].set_color('darkgrey')
ax3.spines['left'].set_linewidth(2)
ax3.spines['bottom'].set_linewidth(2)
ax3.spines['right'].set_linewidth(2)
ax3.spines['top'].set_linewidth(2)
ax3.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='darkgrey')
ax3.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='darkgrey')    
ax3.xaxis.set_ticks_position('bottom')
ax3.yaxis.set_ticks_position('left')

cs = plt.contourf(lat,lev,zdiff_fm,limit,extend='both')
cs2 = plt.contour(lat,lev,zdiff_fm,np.arange(0,1,1),
                  linewidths=0.6,colors='dimgrey')
cs1 = plt.scatter(latq,levq,pvalue_fm,color='k',marker='.',alpha=0.7,
                edgecolor='k',linewidth=0.5)

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.xlim([0,90])
plt.ylim([1000,10])
plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap) 

ax3.annotate(r'\textbf{FM}',
            xy=(0, 0),xytext=(0.35,1.02),xycoords='axes fraction',
            fontsize=25,color='dimgrey',rotation=0)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(bottom=0.21)

plt.savefig(directoryfigure + 'T_diff.png',dpi=300)

