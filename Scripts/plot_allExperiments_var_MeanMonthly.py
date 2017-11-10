"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for MEAN MONTHLY data for all variables.

Notes
-----
    Author : Zachary Labe
    Date   : 6 November 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MeanMonthlyOutput as DM
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/MeanMonthly/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Mean Monthly Data - %s----' % titletime)

#### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = ['LHFLX','SHFLX','FLNS']
runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--CIT}',
               r'\textbf{HIT--CIT}',r'\textbf{FIC--CIT}',
               r'\textbf{FICT--FIT}',r'\textbf{FICT--HIT}']

### Call functions for mean monthly data for polar cap
def readData(varnames):
    """
    Read in data for selected variables and calculate differences
    between experiments
    """

    lat,lon,time,lev,varhit = DM.readMeanExperi(directorydata,
                                                '%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,varfit = DM.readMeanExperi(directorydata,
                                                '%s' % varnames,
                                                'FIT','surface')
    lat,lon,time,lev,varcit = DM.readMeanExperi(directorydata,
                                                '%s' % varnames,
                                                'CIT','surface')
    lat,lon,time,lev,varfic = DM.readMeanExperi(directorydata,
                                                '%s' % varnames,
                                                'FIC','surface')
    lat,lon,time,lev,varfict = DM.readMeanExperi(directorydata,
                                                '%s' % varnames,
                                                'FICT','surface')
    
    ### Compare experiments
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Compute comparisons for experiments - take ensemble average
    diff_FITHIT = np.nanmean(varfit - varhit,axis=0)
    diff_FITCIT = np.nanmean(varfit - varcit,axis=0)
    diff_HITCIT = np.nanmean(varhit - varcit,axis=0)
    diff_FICCIT = np.nanmean(varfic - varcit,axis=0)
    diff_FICTFIT = np.nanmean(varfict - varfit,axis=0)
    diff_FICTHIT = np.nanmean(varfict - varhit,axis=0)
    diffruns = [diff_FITHIT,diff_FITCIT,diff_HITCIT,diff_FICCIT,
                              diff_FICTFIT,diff_FICTHIT]
    
    return diffruns,runs,lat,lon

### Call function to read data for selected variable
diffruns_lh,runs_lh,lat,lon = readData('LHFLX')
diffruns_sh,runs_sh,lat,lon = readData('SHFLX')
diffruns_long,runs_long,lat,lon = readData('FLNS')

### Create 2d array of latitude and longitude
lon2,lat2 = np.meshgrid(lon,lat)

###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
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
        
total_hitq = np.nanmean(runs_lh[0],axis=0) + np.nanmean(runs_sh[0],axis=0)
total_fitq = np.nanmean(runs_lh[1],axis=0) + np.nanmean(runs_sh[1],axis=0) 
total_citq = np.nanmean(runs_lh[2],axis=0) + np.nanmean(runs_sh[2],axis=0) 
total_ficq = np.nanmean(runs_lh[3],axis=0) + np.nanmean(runs_sh[3],axis=0) 
total_fictq = np.nanmean(runs_lh[4],axis=0) + np.nanmean(runs_sh[4],axis=0) 

total_hit = np.append(total_hitq[8:],total_hitq[:3])
total_fit = np.append(total_fitq[8:],total_fitq[:3])
total_cit = np.append(total_citq[8:],total_citq[:3])
total_fic = np.append(total_ficq[8:],total_ficq[:3])
total_fict = np.append(total_fictq[8:],total_fictq[:3])

long_hitq = np.nanmean(runs_long[0],axis=0)
long_fitq = np.nanmean(runs_long[1],axis=0)
long_citq = np.nanmean(runs_long[2],axis=0)
long_ficq = np.nanmean(runs_long[3],axis=0)
long_fictq = np.nanmean(runs_long[4],axis=0)

long_hit = np.append(long_hitq[8:],long_hitq[:3])
long_fit = np.append(long_fitq[8:],long_fitq[:3])
long_cit = np.append(long_citq[8:],long_citq[:3])
long_fic = np.append(long_ficq[8:],long_ficq[:3])
long_fict = np.append(long_fictq[8:],long_fictq[:3])

totallhsh = [total_hit,total_fit,total_cit,total_fic,total_fict]
longg = [long_hit,long_fit,long_cit,long_fic,long_fict]
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')

color=iter(cmocean.cm.thermal(np.linspace(0.1,0.9,len(totallhsh))))
for i in range(len(runnames)):
    c=next(color)
    plt.plot(totallhsh[i],linewidth=2.5,color=c,alpha=1,
             label = r'\textbf{%s}' % runnames[i],linestyle='-')
    
color=iter(cmocean.cm.thermal(np.linspace(0.1,0.9,len(totallhsh))))
for i in range(len(runnames)):
    c=next(color)
    plt.plot(longg[i],linewidth=1.5,color=c,alpha=1,linestyle='--')

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{Fluxes [W/m${^{2}}$]}',color='k',fontsize=13)

plt.yticks(np.arange(25,56,5),list(map(str,np.arange(25,56,5))))
plt.ylim([25,55])

xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
plt.xticks(np.arange(0,7,1),xlabels)
plt.xlim([0,6])

plt.savefig(directoryfigure + 'test.png',dpi=300)
