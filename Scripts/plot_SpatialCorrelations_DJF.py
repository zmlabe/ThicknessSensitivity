"""
Compute spatial correlation (r) between SIT and SIC responses for DJF

Notes
-----
    Author : Zachary Labe
    Date   : 21 February 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput as MO
import cmocean
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)

currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting DJF SIT-SIC spatial correlations - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

months = [r'DJF']
varnames = ['U10','Z30','U300','Z500','SLP','T2M','RNET']

corrvar = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'CIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'FIC',r'CIT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIC--CIT}']
    runs = [varhit,varfit,varfic,varcit]
    
    ### Separate per 2 month periods
    varmo_djf = np.empty((4,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):   
        varmo_djf[i],varmo_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',1)    
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    diff_fithit_djf = np.nanmean(varmo_djf[1] - varmo_djf[0],axis=0)
    diff_ficcit_djf = np.nanmean(varmo_djf[2] - varmo_djf[3],axis=0)
    
    ### Calculate spatial correlation
    corrs = UT.calc_spatialCorr(diff_fithit_djf,diff_ficcit_djf,lat,lon,'yes')
    corrvar.append(corrs)
corrvar = np.asarray(corrvar)
    
#### Save file
np.savetxt(directorydata2 + 'patterncorr_DJF.txt',corrvar.transpose(),
           delimiter=',',fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FIC-CIT to get the relative \n' \
           ' contributions of SIT and SIC [DJF]',newline='\n\n')

###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

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

vals,cor = np.meshgrid(np.arange(2),corrvar)

cs = plt.pcolormesh(cor,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=-1,vmax=1)

for i in range(cor.shape[0]):
    for j in range(cor.shape[1]):
        plt.text(j+0.5,i+0.5,r'\textbf{%+1.2f}' % cor[i,j],fontsize=6,
                 color='k',va='center',ha='center')

cs.set_cmap(cmocean.cm.curl)

ylabels = [r'\textbf{U10}',r'\textbf{Z30}',r'\textbf{U300}',r'\textbf{Z500}',
           r'\textbf{SLP}',r'\textbf{T2M}',r'\textbf{RNET}']
plt.yticks(np.arange(0.5,8.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=0.7)

xlabels = [r'\textbf{DJF}']
plt.xticks(np.arange(0.5,6.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,1])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(-1,2,1)
labels = list(map(str,np.arange(-1,2,1)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Pattern Correlation [R]}',
               color='dimgrey',labelpad=3,fontsize=12)

plt.savefig(directoryfigure + 'patterncorrs_DJF_mesh.png',dpi=300)