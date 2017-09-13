"""
Plot geopotential comparisons between HIT and FIT experiments. These are 
sea ice thickness perturbation experiments using WACCM4. This script is
for DAILY data.

Notes
-----
    Author : Zachary Labe
    Date   : 6 September 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput as DO
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
print('\n' '----Plotting Daily Geopotential - %s----' % titletime)

#### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for surface temperature data
lat,lon,time,lev,z_h = DO.readMeanExperi(directorydata,'GEOP',
                                        'HIT','profile')
lat,lon,time,lev,z_f = DO.readMeanExperi(directorydata,'GEOP',
                                        'FIT','profile')
                                        
#### Calculate significance
stat,pvalue = UT.calc_indttest(z_h,z_f)
                                                                    
### Calculate ensemble mean
z_diff= np.nanmean(z_f-z_h,axis=0)       

############################################################################
############################################################################
############################################################################
##### Plot temperature profile
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-150,150.1,15)
barlim = np.arange(-150,151,75)
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
timeq = np.arange(0,212,1)
timeqq,levq = np.meshgrid(timeq,lev)

fig = plt.figure()
ax1 = plt.subplot(111)

ax1.spines['top'].set_color('dimgrey')
ax1.spines['right'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='dimgrey')
ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='dimgrey')    
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')


cs = plt.contourf(timeq,lev,z_diff.transpose(),limit,extend='both')
                  
plt.contourf(timeqq,levq,pvalue.transpose(),colors='None',hatches=['////'],
             linewidth=5)                  
                  
#pvalue[np.where(np.isnan(pvalue))] = 0.0                  
#                  
#cs1 = plt.contour(timeqq,levq,pvalue.transpose(),np.arange(0,2,1),linestyles='-')

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')
plt.ylabel(r'\textbf{Pressure (hPa)}',color='dimgrey',fontsize=15,
                     labelpad=1)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()
plt.xlim([30,210])
plt.ylim([1000,10])

cmap = ncm.cmap('temp_diff_18lev')            
cs.set_cmap(cmap) 

#ax1.annotate(r'\textbf{ON}',
#            xy=(0, 0),xytext=(0.34,1.02),xycoords='axes fraction',
#            fontsize=25,color='dimgrey',rotation=0)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(bottom=0.21)

plt.savefig(directoryfigure + 'Z_daily_diff_FIT-HIT.png',dpi=300)
print('Completed: Script done!')                                