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
varnames = ['U10','Z30','U300','Z500','SLP','T2M','THICK','RNET']

corrvar = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FIC','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'FICT',r'FIC']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FICT--FIC}']
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
    diff_fictfit = np.nanmean(varmo_fict - varmo_fit,axis=0) 
    
    corrs = []
    for i in range(len(diff_fithit)):
        corrsq = UT.calc_spatialCorr(diff_fithit[i],diff_fictfit[i],
                                     lat,lon,'yes')
        corrs.append(corrsq)
    corrvar.append(corrs)
corrvar = np.asarray(corrvar)

#### Save file
np.savetxt(directorydata2 + 'patterncorr.txt',corrvar.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FICT-FIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]',newline='\n\n')

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

plt.plot([0]*len(corrvar),linewidth=3,color='dimgrey',linestyle='--')

color=iter(ncm.cmap('MPL_gnuplot2')(np.linspace(0,0.85,len(corrvar))))
for i in range(len(corrvar)):
    c=next(color)
    corrvar[np.where(corrvar == 0.)] = np.nan
    plt.plot(corrvar[i],linewidth=2.5,color=c,alpha=1,
             label = r'\textbf{%s}' % varnames[i],linestyle='-',
             marker='o',markersize=9)

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=4)
plt.ylabel(r'\textbf{Pattern Correlation [r]}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-1,1.1,0.5),list(map(str,np.arange(-1,1.1,0.5))))
plt.ylim([-1,1])

xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
plt.xticks(np.arange(0,6,1),xlabels)
plt.xlim([0,5])

ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.7)

plt.savefig(directoryfigure + 'patterncorrs_monthly.png',dpi=300)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#    
##    ### Calculate pearsonr correlation point by point
##    patternco = np.empty((len(months), diff_fithit.shape[2],
##                          diff_fithit.shape[3]))
##    patternpval = np.empty((len(months), diff_fithit.shape[2],
##                      diff_fithit.shape[3]))
##    for mo in range(len(months)):
##        for x in range(diff_fithit.shape[2]):
##            for y in range(diff_fithit.shape[3]):
##               patternco[mo,x,y],patternpval[mo,x,y] = sts.pearsonr(
##                                                       diff_fithit[:,mo,x,y],
##                                                       diff_fictfit[:,mo,x,y])
##               
##        print('Completed: pattern correlation for -- %s!' % months[mo])
##    
#################################################################################
#################################################################################
#################################################################################
##### Plot figure
##    plt.rc('text',usetex=True)
##    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
##    
##    ### Set limits for contours and colorbars
##    limit = np.arange(-1,1.05,0.05)
##    barlim = np.arange(-1,2,1)
##    
##    fig = plt.figure()
##    for i in range(len(months)):
##        
##        var = patternco[i]
###        pvar = pruns_djf[i]
##            
##        ax1 = plt.subplot(2,3,i+1)
##        m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
##                    area_thresh=10000.)
##        
##        var, lons_cyclic = addcyclic(var, lon)
##        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
##        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
##        x, y = m(lon2d, lat2d)
##        
###        pvar,lons_cyclic = addcyclic(pvar, lon)
###        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
##                  
##        m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
##        m.drawcoastlines(color='dimgray',linewidth=0.65)
##        
##        cs = m.contourf(x,y,var,limit,extend='both')
###        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
##        
##        cmap = cmocean.cm.balance  
##        cs.set_cmap(cmap)
##            
##        m.drawcoastlines(color='dimgray',linewidth=0.8)
##        
##        ### Add experiment text to subplot
##        ax1.annotate(r'\textbf{%s}' % months[i],xy=(0,0),xytext=(0.865,0.90),
##                     textcoords='axes fraction',color='dimgrey',fontsize=11,
##                     rotation=320,ha='center',va='center')
##    
##    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
##    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
##                        extend='max',extendfrac=0.07,drawedges=False)
##    cbar.set_label(r'\textbf{Pattern Correlation (r)}',fontsize=11,
##                             color='dimgray')
##    cbar.set_ticks(barlim)
##    cbar.set_ticklabels(list(map(str,barlim))) 
##    cbar.ax.tick_params(axis='x', size=.01)
##    cbar.outline.set_edgecolor('dimgrey')
##    
##    plt.subplots_adjust(wspace=0.01)
##    plt.subplots_adjust(hspace=0.01)
##    plt.subplots_adjust(bottom=0.15)
##    
##    plt.savefig(directoryfigure + '%s_patterncorr.png' % varnames[v],
##                dpi=300)