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
import calc_Utilities as UT
import read_MonthlyOutput as MO

### Define directories
directorydata = '/home/zlabe/surt/LENS/ForcingPerturb/' 
#directorydata = '/surtsey/zlabe/simu/'
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

lons2,lats2 = np.meshgrid(lons,lats)

print('Completed: Data read!')

### Create 2d array of latitude and longitude
lons2,lats2 = np.meshgrid(lons,lats)

### Average over DJF 
sith[np.where(sith <= 0.15)] = np.nan            
varh = sith
sitf[np.where(sitf <= 0.15)] = np.nan            
varf = sitf

sich[np.where(sich <= 0.15)] = np.nan            
varch = sich*100 # convert SIC to 1-100%
sicf[np.where(sicf <= 0.15)] = np.nan            
varcf = sicf*100 # convert SIC to 1-100%

### Use land/Arctic mask
varh[np.where(varh == 2.)] = np.nan
varhtemp = varh.copy()
varhtemp[np.where(varhtemp > 0)] = 1.
varf = varf*varhtemp
varf[np.where(varf == 2.)] = np.nan

#### Create 2d array of latitude and longitude
latq = np.where(lats>=65)[0]

lons2,lats2 = np.meshgrid(lons,lats[latq])

varf = varf[:,latq,:]
varh = varh[:,latq,:]
varcf = varcf[:,latq,:]
varch = varch[:,latq,:]

### Calculate months
varnf = np.append(varf[9:],varf[:3],axis=0)
varnh = np.append(varh[9:],varh[:3],axis=0)

varncf = np.append(varcf[9:],varcf[:3],axis=0)
varnch = np.append(varch[9:],varch[:3],axis=0)

### Calculate monthly averages
sitfmean = UT.calc_weightedAve(varnf,lats2)
sithmean = UT.calc_weightedAve(varnh,lats2)

sicfmean = UT.calc_weightedAve(varncf,lats2)
sichmean = UT.calc_weightedAve(varnch,lats2)

print('Completed: Data processed!')
###############################################################################
###############################################################################
###############################################################################
### Create subplots of sea ice anomalies 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
fig = plt.figure()
ax = plt.subplot(211)

plt.plot(sithmean,linewidth=3,marker='o',markersize=6,
         color=cmocean.cm.balance(0.2),label='Historical')
plt.plot(sitfmean,linewidth=3,marker='o',markersize=6,
         color=cmocean.cm.balance(0.8),label='Future')

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.yticks(np.arange(0,7,1),list(map(str,np.arange(0,7,1))))
plt.ylim([0,3])
xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
plt.xticks(np.arange(0,7,1),xlabels)
plt.xlim([0,5])

ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.45)

plt.ylabel(r'\textbf{Sea Ice Thickness (m)}',color='k',fontsize=8,labelpad=9)

#######################
#######################
#######################

ax = plt.subplot(212)

plt.plot(sichmean,linewidth=3,marker='o',markersize=6,
         color=cmocean.cm.balance(0.2),label='Historical')
plt.plot(sicfmean,linewidth=3,marker='o',markersize=6,
         color=cmocean.cm.balance(0.8),label='Future')

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.yticks(np.arange(0,101,10),list(map(str,np.arange(0,101,10))))
plt.ylim([30,90])
xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
plt.xticks(np.arange(0,7,1),xlabels)
plt.xlim([0,5])

plt.ylabel(r'\textbf{Sea Ice Concentration (\%)}',color='k',fontsize=8)

plt.legend(shadow=False,fontsize=6,loc='lower right',
           fancybox=True,frameon=False,ncol=1)
ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.45)

fig.subplots_adjust(hspace=0.4)
        
plt.savefig(directoryfigure + 'monthly_seaiceanoms.png_areamean.png',dpi=300)