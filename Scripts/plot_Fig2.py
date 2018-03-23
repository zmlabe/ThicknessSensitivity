"""
Plot figure 2 in manuscript for dynamical responses to sea ice loss in WACCM4
experiments [FIT-HIT, FIC-CIT, FICT-HIT]. Current variables include T2M and
RNET. Time period includes December through February [DJF].

Notes
-----
    Author : Zachary Labe
    Date   : 4 February 2018
"""

### Import modules
import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from matplotlib import cbook
from matplotlib.colors import Normalize
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import read_MeanMonthlyOutput as DM
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
print('\n' '----Plotting Fig 2 - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Define constants
varnames = ['T2M']
experiments = [r'\textbf{$\Delta$SIT}',r'\textbf{$\Delta$SIC}',
               r'\textbf{$\Delta$NET}']
runnames = [r'HIT',r'FIT',r'HIT2',r'FICT2',r'FICT']

### Functions to read temperature
def readTemp(varnames):
    """
    Read in temperature data for selected variables and calculate differences
    between experiments
    """
    for v in range(len(varnames)):
        ### Call function for T2M data from reach run
        lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'HIT',
                                                'surface')
        lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'FIT',
                                                'surface')
        lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'CIT',
                                                'surface')
        lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'FIC',
                                                'surface')
        lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                                 '%s' % varnames[v],'FICT',
                                                 'surface')
        
        ### Create 2d array of latitude and longitude
        lon2,lat2 = np.meshgrid(lon,lat)
        
        ### Concatonate runs
        runs = [varhit,varfit,varcit,varfic,varfict]
        
        ### Separate per periods (DJF)
        var_djf = np.empty((5,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
        for i in range(len(runs)):
            var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,
                                                   'surface',1)    
        
        ### Compute comparisons for FM - taken ensemble average
        diff_FITHIT = np.nanmean(var_djf[1] - var_djf[0],axis=0)
        diff_FICCIT = np.nanmean(var_djf[3] - var_djf[2],axis=0)
        diff_FICTHIT = np.nanmean(var_djf[4] - var_djf[0],axis=0)
        diffruns_djf = [diff_FITHIT,diff_FICCIT,diff_FICTHIT]
        
        ### Calculate significance for FM
        stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(var_djf[1],var_djf[0])
        stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(var_djf[3],var_djf[2])
        stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(var_djf[4],var_djf[0])
        pruns_djf = [pvalue_FITHIT,pvalue_FICCIT,pvalue_FICTHIT]
        
    return diffruns_djf,pruns_djf,lat,lon

###############################################################################
###############################################################################
###############################################################################
# Function to read surface heat flux data
def readFlux(varnames):
    """
    Read in heat flux data for selected variables and calculate differences
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
    diff_FICCIT = np.nanmean(varfic - varcit,axis=0)
    diff_FICTHIT = np.nanmean(varfict - varhit,axis=0)
    diffruns = [diff_FITHIT,diff_FICCIT,diff_FICTHIT]
    
    return diffruns,runs,lat,lon

###########################################################################
###########################################################################
###########################################################################
# Read data for net surface energy budget
difftotallhshq = np.genfromtxt(directorydata2+'weightedsic_SHLH.txt',
                              skip_header=2,delimiter=',')
difftotallhsh = difftotallhshq.transpose()
temps,ptemps,lat,lon = readTemp(varnames)

### Create 2d array of latitude and longitude
lon2,lat2 = np.meshgrid(lon,lat)
    
###########################################################################
###########################################################################
###########################################################################
#### Plot surface temperature and rnet

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

### Begin plot
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
fig = plt.figure()

###############################################################################
### Set limits
var = temps[0]
pvar = ptemps[0]
       
limit = np.arange(-5,15.1,0.25)
barlim = np.arange(-5,16,5)
    
ax1 = plt.subplot(2,2,1)
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

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)    

m.drawcoastlines(color='dimgrey',linewidth=0.7)

ax1.annotate(r'%s' % experiments[0],xy=(0,0),xytext=(0.1,0.90),
             textcoords='axes fraction',color='k',fontsize=16,
             rotation=45,ha='center',va='center')
ax1.annotate(r'\textbf{[%s]}' % 'a',xy=(0,0),
        xytext=(0.89,0.9),xycoords='axes fraction',
        color='dimgrey',fontsize=7)

###############################################################################
### Set limits
var = temps[1]
pvar = ptemps[1]
       
limit = np.arange(-5,15.1,0.25)
barlim = np.arange(-5,16,5)
    
ax1 = plt.subplot(2,2,2)
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

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)    

m.drawcoastlines(color='dimgrey',linewidth=0.7)

ax1.annotate(r'%s' % experiments[1],xy=(0,0),xytext=(0.1,0.90),
             textcoords='axes fraction',color='k',fontsize=16,
             rotation=45,ha='center',va='center')
ax1.annotate(r'\textbf{[%s]}' % 'b',xy=(0,0),
        xytext=(0.89,0.9),xycoords='axes fraction',
        color='dimgrey',fontsize=7)

###############################################################################
### Set limits
var = temps[2]
pvar = ptemps[2]
       
limit = np.arange(-5,15.1,0.25)
barlim = np.arange(-5,16,5)
    
ax1 = plt.subplot(2,2,3)
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

cmap = ncm.cmap('NCV_blu_red')            
cs.set_cmap(cmap)    

m.drawcoastlines(color='dimgrey',linewidth=0.7)

ax1.annotate(r'%s' % experiments[2],xy=(0,0),xytext=(0.08,0.88),
             textcoords='axes fraction',color='k',fontsize=16,
             rotation=45,ha='center',va='center')
ax1.annotate(r'\textbf{[W/m${^{2}}$]}',xy=(0,0),xytext=(1.18,1),
     textcoords='axes fraction',color='dimgrey',fontsize=7,
     rotation=0,ha='center',va='center')
ax1.annotate(r'\textbf{[%s]}' % 'c',xy=(0,0),
        xytext=(0.89,0.9),xycoords='axes fraction',
        color='dimgrey',fontsize=7)

###############################################################################
ax = plt.axes([.543, .183, .24, .31]) 

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

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

color=iter(cmocean.cm.matter(np.linspace(0.3,1,len(difftotallhsh))))
for i in range(len(difftotallhsh)):
    c=next(color)
    plt.plot(difftotallhsh[i],linewidth=2,color=c,alpha=1,
             label = r'\textbf{%s}' % experiments[i],linestyle='-',
             marker='o',markersize=4)

plt.legend(shadow=False,fontsize=5,loc='lower left',
           fancybox=True,frameon=True,ncol=3,bbox_to_anchor=(0.05, 0.13),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4,
           edgecolor='dimgrey')

plt.yticks(np.arange(0,126,25),list(map(str,np.arange(0,126,25))),fontsize=6)
plt.ylim([0,100])

xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
plt.xticks(np.arange(0,7,1),xlabels,fontsize=6)
plt.xlim([0,6])

ax.annotate(r'\textbf{[%s]}' % 'd',xy=(0,0),
        xytext=(0.89,0.9),xycoords='axes fraction',
        color='dimgrey',fontsize=7)
    
cbar_ax = fig.add_axes([0.31,0.09,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{[T2M]$^\circ$C}',
               fontsize=13,color='dimgrey',labelpad=2)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelsize=8,pad=2.1)
cbar.outline.set_edgecolor('dimgrey')
    
fig.subplots_adjust(wspace=-0.4,hspace=0,bottom=0.15)
    
plt.savefig(directoryfigure + 'Fig2.png',dpi=600)
    
print('Completed: Script done!')

