"""
Plot figure 3 in manuscript for dynamical responses to sea ice loss in WACCM4
experiments [FIT-HIT, FIC-CIT, FICT-HIT]. Current variables include SLP,
Z500, and Z30. Time period includes December through February [DJF].

Notes
-----
    Author : Zachary Labe
    Date   : 4 April 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import cmocean
import itertools

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
print('\n' '----Plotting Poster Figure 1 - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnamesn = ['SLP','SLP']
experimentsn = [r'\textbf{$\Delta$SIT}',r'\textbf{$\Delta$SIC}',
               r'\textbf{$\Delta$NET}']
runnamesn = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
letters = ["a","b","c","d","e","f","g","h","i"]
def calcVarResp(varnames):
    ### Call function for variable data from each run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames,'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames,'FIT','surface')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                            '%s' % varnames,'CIT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                            '%s' % varnames,'FIC','surface')
    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames,'FICT','surface')
    
    ### Concatonate runs
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Separate per periods (ON,DJ,FM)
    var_djf = np.empty((5,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):
        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',1)   
    
    ### Compute climatology    
    climohit = np.nanmean(var_djf[1],axis=0)
    climocit = np.nanmean(var_djf[2],axis=0)
    climo_djf = [climohit,climocit,climohit]
    
    ### Compute comparisons for FM - taken ensemble average
    diff_FITHIT = np.nanmean(var_djf[1] - var_djf[0],axis=0)
    diff_FICCIT = np.nanmean(var_djf[3] - var_djf[2],axis=0)
    diff_FICTHIT = np.nanmean(var_djf[4] - var_djf[0],axis=0)
    diffruns_djf = [diff_FITHIT,diff_FICCIT]
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(var_djf[1],var_djf[0])
    stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(var_djf[3],var_djf[2])
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(var_djf[4],var_djf[0])
    pruns_djf = [pvalue_FITHIT,pvalue_FICCIT]
    
    return diffruns_djf,pruns_djf,climo_djf,lat,lon

### Call variables
diffslp,pslp,climoslp,lat,lon = calcVarResp('SLP')

## Create lists for plotting
variables = list(itertools.chain(*[diffslp]))
pvalues = list(itertools.chain(*[pslp]))
climos = list(itertools.chain(*[climoslp]))
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJ
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for v in range(len(variables)):
    ax = plt.subplot(1,2,v+1)
    
    ### Retrieve variables and pvalues
    var = variables[v]
    pvar = pvalues[v]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'SLP':
        limit = np.arange(-4,4.1,0.05)
        barlim = np.arange(-4,5,4)
    elif varnamesn[v] == 'Z500':
        limit = np.arange(-60,60.1,5)
        barlim = np.arange(-60,61,30) 
    elif varnamesn[v] == 'Z30':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,50) 
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    climoq,lons_cyclic = addcyclic(climos[v], lon)
    climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    
    cs = m.contourf(x,y,var,limit,extend='both')
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
    if varnamesn[v] == 'Z30': # the interval is 250 m 
        cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                        colors='k',linewidths=1.1,zorder=10)
    if varnamesn[v] == 'RNET':
        m.drawcoastlines(color='darkgray',linewidth=0.3)
        m.fillcontinents(color='dimgrey')
    else:
        m.drawcoastlines(color='dimgrey',linewidth=0.8)
    
    if varnamesn[v] == 'SLP':
        cmap = cmocean.cm.balance          
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'Z500':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'Z30':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    
    if v == 0:        
        ax.annotate(r'\textbf{SEA ICE}',xy=(0,0),xytext=(0.5,1.16),
                     textcoords='axes fraction',color='dimgrey',
                     fontsize=22,rotation=0,ha='center',va='center')
        ax.annotate(r'\textbf{\underline{THICKNESS}}',xy=(0,0),xytext=(0.5,1.045),
             textcoords='axes fraction',color='dimgrey',
             fontsize=22,rotation=0,ha='center',va='center')
    elif v == 1:        
        ax.annotate(r'\textbf{SEA ICE}',xy=(0,0),xytext=(0.5,1.16),
                     textcoords='axes fraction',color='dimgrey',
                     fontsize=22,rotation=0,ha='center',va='center')
        ax.annotate(r'\textbf{\underline{CONCENTRATION}}',xy=(0,0),xytext=(0.5,1.045),
                 textcoords='axes fraction',color='dimgrey',
                 fontsize=22,rotation=0,ha='center',va='center')
        
    ax.set_aspect('equal')
            
    ###########################################################################
    cbar_ax = fig.add_axes([0.304,0.11,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    if varnamesn[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'Z30':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'SLP':
        cbar.set_label(r'\textbf{Sea Level Pressure [hPa]}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'U10' or varnamesn[v] == 'U300':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnamesn[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnamesn[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnamesn[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
    elif varnamesn[v] == 'RNET':
        cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=11,color='dimgray') 
    
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))    
    cbar.ax.tick_params(axis='x', size=.001,labelcolor='dimgrey')
    cbar.outline.set_edgecolor('dimgrey')
    plt.tight_layout()
    
    plt.annotate(r'\textbf{December -- February}',xy=(0,0),xytext=(0.5,0.165),
         textcoords='figure fraction',color='dimgrey',
         fontsize=8,rotation=0,ha='center',va='center')
    
    
plt.savefig(directoryfigure + 'PosterFig1.png',dpi=1000)

print('Completed: Script done!')

