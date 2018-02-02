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
import calc_Utilities as UT
import cmocean

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

rootvar = []
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
    runnames = [r'HIT',r'FIT',r'FICT',r'FIC']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIC--CIT}']
    runs = [varhit,varfit,varfic,varcit]
    
    ### Separate per months
    varmo_fit = np.append(varfit[:,9:,:,:],varfit[:,0:3,:,:],
             axis=1)
    varmo_hit = np.append(varhit[:,9:,:,:],varhit[:,0:3,:,:],
             axis=1)
    varmo_fic = np.append(varfic[:,9:,:,:],varfic[:,0:3,:,:],
              axis=1)
    varmo_cit = np.append(varcit[:,9:,:,:],varcit[:,0:3,:,:],
          axis=1)
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    diff_fithit = np.nanmean(varmo_fit - varmo_hit,axis=0)
    diff_ficcit= np.nanmean(varmo_fic - varmo_cit,axis=0) 
    
    root = []
    for i in range(len(diff_fithit)):
        rootq = UT.calc_RMSE(diff_fithit[i],diff_ficcit[i],lat,lon,'yes')
        root.append(rootq)
    rootvar.append(root)
rootvar = np.asarray(rootvar)

##### Save file
np.savetxt(directorydata2 + 'rmse.txt',rootvar.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains root mean square error ' \
           '\n between FIT-HIT and FIC-CIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]',newline='\n\n')

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

cs = plt.pcolormesh(rootvar,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=0,vmax=60)

for i in range(rootvar.shape[0]):
    for j in range(rootvar.shape[1]):
        plt.text(j+0.5,i+0.5,r'\textbf{%1.2f}' % rootvar[i,j],fontsize=6,
                 color='tomato',va='center',ha='center')

cs.set_cmap(cmocean.cm.ice_r)

ylabels = [r'\textbf{U10}',r'\textbf{Z30}',r'\textbf{U300}',r'\textbf{Z500}',
           r'\textbf{SLP}',r'\textbf{T2M}',r'\textbf{THICK}',
           r'\textbf{RNET}']
plt.yticks(np.arange(0.5,8.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=0.7)

xlabels = [r'\textbf{OCT}',r'\textbf{NOV}',r'\textbf{DEC}',
           r'\textbf{JAN}',r'\textbf{FEB}',r'\textbf{MAR}']
plt.xticks(np.arange(0.5,6.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,6])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(0,61,30)
labels = list(map(str,np.arange(0,61,30)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Root Mean Square Error [RMSE]}',
               color='dimgrey',labelpad=3,fontsize=12)

plt.subplots_adjust(top=0.8)

plt.savefig(directoryfigure + 'rmse_monthly_mesh.png',dpi=300)

