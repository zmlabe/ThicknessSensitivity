"""
Plot figure 3 in manuscript for dynamical responses to sea ice loss in WACCM4
experiments [FIT-HIT, FIC-CIT, FICT-HIT]. Current variables include SLP,
Z500, and Z30. Time period includes December through February [DJF].

Notes
-----
    Author : Zachary Labe
    Date   : 3 February 2018
"""

### Import modules
import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
from matplotlib import cbook
from matplotlib.colors import Normalize
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
print('\n' '----Plotting Figure 3 - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnamesn = np.repeat(['T2M','T925','T850'],2)
experimentsn = [r'\textbf{$\Delta$SIT}',r'\textbf{$\Delta$SIC}']
runnamesn = [r'HIT',r'FIT',r'CIT',r'FIC']
letters = ["a","b","c","d","e","f"]
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
    
    ### Concatonate runs
    runs = [varhit,varfit,varcit,varfic]
    
    ### Separate per periods (ON,DJ,FM)
    var_djf = np.empty((4,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):
        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',1)   
    
    ### Compute comparisons for FM - taken ensemble average
    diff_FITHIT = np.nanmean(var_djf[1] - var_djf[0],axis=0)
    diff_FICCIT = np.nanmean(var_djf[3] - var_djf[2],axis=0)
    diffruns_djf = [diff_FITHIT,diff_FICCIT]
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(var_djf[1],var_djf[0])
    stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(var_djf[3],var_djf[2])
    pruns_djf = [pvalue_FITHIT,pvalue_FICCIT]
    
    return diffruns_djf,pruns_djf,lat,lon

### Call variables
difft2m,pt2m,lat,lon = calcVarResp('T2M')
difft925,pt925,lat,lon = calcVarResp('T925')
difft850,pt850,lat,lon = calcVarResp('T850')

## Create lists for plotting
variables = list(itertools.chain(*[difft2m,difft925,difft850]))
pvalues = list(itertools.chain(*[pt2m,pt925,pt850]))
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJF

### Functions to center colormap on 0 - white    
class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint
norm = MidPointNorm(midpoint=0)

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for v in range(len(variables)):
    ax = plt.subplot(3,2,v+1)
    
    ### Retrieve variables and pvalues
    var = variables[v]
    pvar = pvalues[v]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'T2M':
        limit = np.arange(-5,15.1,0.25)
        barlim = np.arange(-5,16,5)
    elif varnamesn[v] == 'T925':
        limit = np.arange(-5,15.1,0.25)
        barlim = np.arange(-5,16,5)
    elif varnamesn[v] == 'T850':
        limit = np.arange(-5,15.1,0.25)
        barlim = np.arange(-5,16,5)
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    
    cs = m.contourf(x,y,var,limit,extend='max',norm=norm)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])

    m.drawcoastlines(color='dimgrey',linewidth=0.6)
    
    if varnamesn[v] == 'T2M':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap)    
    elif varnamesn[v] == 'T925':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap)    
    elif varnamesn[v] == 'T850':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap)    
            
    ### Add experiment text to subplot
    if any([v == 0,v == 2,v == 4]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.14,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=16,rotation=90,ha='center',va='center')
    if any([v == 0,v == 1]):
        ax.annotate(r'\textbf{%s}' % experimentsn[v],xy=(0,0),xytext=(0.5,1.09),
                     textcoords='axes fraction',color='k',
                     fontsize=13,rotation=0,ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=7)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 55:
        cbar_ax = fig.add_axes([0.74,0.65,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'T2M':
            cbar.set_label(r'\textbf{[T2M]$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)
        elif varnamesn[v] == 'T925':
            cbar.set_label(r'\textbf{[T2M]$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)
        elif varnamesn[v] == 'T850':
            cbar.set_label(r'\textbf{[T2M]$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 3:
        cbar_ax = fig.add_axes([0.365,0.09,0.3,0.023])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'T2M':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=13,color='dimgrey',labelpad=2.)
        elif varnamesn[v] == 'T925':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=13,color='dimgrey',labelpad=2.)
        elif varnamesn[v] == 'T850':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=13,color='dimgrey',labelpad=2.)     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=8,pad=2.1) 
        ticklabs = cbar.ax.get_xticklabels()
        cbar.ax.set_xticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='x', size=.0001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 55:
        cbar_ax = fig.add_axes([0.74,0.15,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'T2M':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)
        elif varnamesn[v] == 'T925':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)
        elif varnamesn[v] == 'T850':
            cbar.set_label(r'\textbf{$^\circ$C}',
                           fontsize=9,color='dimgrey',labelpad=1.2)        
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
fig.subplots_adjust(wspace=-0.6,hspace=0)
    
plt.savefig(directoryfigure + 'T_siberiaCooling_DJF_FITHIT.png',dpi=900)

print('Completed: Script done!')

