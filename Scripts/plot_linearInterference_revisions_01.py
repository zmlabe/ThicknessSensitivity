"""
Calculate linear interference. Plot created during revisions #1.

Notes
-----
    Author : Zachary Labe
    Date   : 1 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import cmocean
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
print('\n' '----Plotting linear interference - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

historicalforced = []
futureforced = []
lonss = []
varnames = ['GEOPxwave_all','GEOPxwave1','GEOPxwave2']
for v in range(len(varnames)):
    ### Call function for geopotential height data from reach run
    lat,lon1,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','profile')
    lat,lon1,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','profile')

    ### Modify lons 
    lon1[-1] = 360

    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT']
    experiments = [r'\textbf{FIT--HIT}']
    runs = [varhit,varfit]
    
    ### Separate per periods (m)
    var_m = [varhit[:,2,:,:],varfit[:,2,:,:]]
    
    ### Take at 60N 
    latq = np.where((lat>59) & (lat<61))[0]
    
    ### Compute comparisons for M - taken ensemble average at 60N (index 79)
    diff_FITHIT = np.nanmean(var_m[1] - var_m[0],axis=0)
    diffruns_m = diff_FITHIT[:,latq,:].squeeze()
    future = np.nanmean(var_m[1],axis=0)
    futureforcedq = np.nanmean(var_m[1][:,:,latq,:].squeeze(),axis=0)
    
    historicalforcedq = np.nanmean(var_m[0][:,:,latq,:].squeeze(),axis=0)
    
    historicalforced.append(historicalforcedq)
    futureforced.append(diffruns_m)
    lonss.append(lon1)
    
###########################################################################
###########################################################################
###########################################################################
#### Plot climatological waves
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])

fig = plt.figure()
for i in range(3):
    ax1 = plt.subplot(3,1,i+1)
    
    ### Calculate correlations
    corr = np.corrcoef(historicalforced[i].ravel(),futureforced[i].ravel())[0][1]
    
    lonq,levq = np.meshgrid(lonss[i],lev)
    
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
            
    cs = plt.contourf(lonq,levq,historicalforced[i],100,extend='both') 
    cs1 = plt.contour(lonq,levq,futureforced[i],np.arange(-100,101,5),
                      colors='k',linewidths=1.5) 
    
    waveq = ['Total Wave','Wave 1','Wave 2']
    plt.text(384,80,r'\textbf{%s}' % waveq[i],color='k',
             fontsize=8,rotation=0,ha='center',va='center')
    plt.text(384,140,r'\textbf{R=%s}' % str(corr)[:4],color='k',
         fontsize=8,rotation=0,ha='center',va='center')
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    xxlabels = ['0','60E','120E','180','120W','60W','0']
    
    if i==2:
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,361,60),xxlabels,fontsize=8)
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
    else:
        plt.ylim([1000,10])
        plt.xticks([])
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
    
    cmap = cmocean.cm.balance            
    cs.set_cmap(cmap) 

plt.savefig(directoryfigure + 'FITHIT_linearInterference_March.png',dpi=300)
print('Completed: Script done!')


print('Completed: Script done!')

