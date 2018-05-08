"""
Plots OCT-MAR spatial correlations -- test script so far!

Notes
-----
    Author : Zachary Labe
    Date   : 16 November 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput as MO
import cmocean
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting spatial correlations - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

months = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
varnames = ['U10','Z30','U300','Z500','SLP','T2M']

corrvar = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'CIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'FIC',r'CIT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIC--CIT}']
    runs = [varhit,varfit,varfict,varfic]
    
    ### Separate per months
    varmo_fit = np.append(varfit[:,9:,:,:],varfit[:,0:3,:,:],
             axis=1)
    varmo_hit = np.append(varhit[:,9:,:,:],varhit[:,0:3,:,:],
             axis=1)
    varmo_fict = np.append(varfict[:,9:,:,:],varfict[:,0:3,:,:],
              axis=1)
    varmo_fic = np.append(varfic[:,9:,:,:],varfic[:,0:3,:,:],
          axis=1)
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    diff_fithit = np.nanmean(varmo_fit - varmo_hit,axis=0)
    diff_fictfit = np.nanmean(varmo_fict - varmo_fic,axis=0) 
    
    ### Calculate significance
    pvalue_FITHITdjfq = np.empty((12,lat.shape[0],lon.shape[0]))
    pvalue_FICCITdjfq = np.empty((12,lat.shape[0],lon.shape[0]))
    for i in range(6):
        stat_FITHITdjf,pvalue_FITHITdjfq[i] = UT.calc_indttest(varmo_fit[:,i,:,:],
                                                           varmo_hit[:,i,:,:])
        stat_FICCITdjf,pvalue_FICCITdjfq[i] = UT.calc_indttest(varmo_fict[:,i,:,:],
                                                           varmo_fic[:,i,:,:])
    
    ### Create mask of significant values
    pvalue_FITHITdjf[np.where(np.isnan(pvalue_FITHITdjf))] = 0.0
    pvalue_FICCITdjf[np.where(np.isnan(pvalue_FICCITdjf))] = 0.0
        
    pvalue_FITHIT = pvalue_FITHITdjf
    pvalue_FICCIT = pvalue_FICCITdjf
    
    ### Keep only values significant in both SIT and SIC responses    
    diff_fithitq = diff_fithit * pvalue_FITHITdjf   
    diff_ficcitq = diff_fictfit * pvalue_FICCITdjf
    
    ### Calculate the correlations
    corrs = []
    for i in range(len(diff_fithit)):
        corrsq = UT.calc_spatialCorr(diff_fithitq[i],diff_ficcitq[i],
                                     lat,lon,'yes')
        corrs.append(corrsq)
    corrvar.append(corrs)
corrvar = np.asarray(corrvar)

##### Save file
np.savetxt(directorydata2 + 'patterncorr_sig.txt',corrvar.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FIC-CIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR, 95% significance]',
           newline='\n\n')

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


#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')
#
#plt.plot([0]*len(corrvar),linewidth=3,color='dimgrey',linestyle='--')
#
#color=iter(ncm.cmap('MPL_gnuplot2')(np.linspace(0,0.85,len(corrvar))))
#for i in range(len(corrvar)):
#    c=next(color)
#    corrvar[np.where(corrvar == 0.)] = np.nan
#    plt.plot(corrvar[i],linewidth=2.5,color=c,alpha=1,
#             label = r'\textbf{%s}' % varnames[i],linestyle='-',
#             marker='o',markersize=9)
#
#plt.legend(shadow=False,fontsize=9,loc='lower center',
#           fancybox=True,frameon=False,ncol=4)
#plt.ylabel(r'\textbf{Pattern Correlation [r]}',color='dimgrey',fontsize=13)
#
#plt.yticks(np.arange(-1,1.1,0.5),list(map(str,np.arange(-1,1.1,0.5))))
#plt.ylim([-1,1])
#
#xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
#plt.xticks(np.arange(0,6,1),xlabels)
#plt.xlim([0,5])
#
#ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.7)
#
#plt.savefig(directoryfigure + 'patterncorrs_monthly.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(111)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on')

cs = plt.pcolormesh(corrvar,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-1,vmax=1)

for i in range(corrvar.shape[0]):
    for j in range(corrvar.shape[1]):
        plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % corrvar[i,j],fontsize=7,
                 color='k',va='center',ha='center')

cs.set_cmap(cmocean.cm.curl)

ylabels = [r'\textbf{U10}',r'\textbf{Z30}',r'\textbf{U300}',r'\textbf{Z500}',
           r'\textbf{SLP}',r'\textbf{T2M}']
plt.yticks(np.arange(0.5,7.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=-2)

xlabels = [r'\textbf{OCT}',r'\textbf{NOV}',r'\textbf{DEC}',
           r'\textbf{JAN}',r'\textbf{FEB}',r'\textbf{MAR}']
plt.xticks(np.arange(0.5,6.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,6])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(-1,2,1)
labels = list(map(str,np.arange(-1,2,1)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Pattern Correlation [r]}',
               color='dimgrey',labelpad=3,fontsize=12)

plt.savefig(directoryfigure + 'patterncorrs_monthly_mesh_sig.png',dpi=900)