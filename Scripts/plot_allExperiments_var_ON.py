"""
Plot temperature comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, CIT, FIC, FICT

Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import cmocean

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
print('\n' '----Plotting ON variable data - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnames = ['Z500','Z50','Z30','SLP','T2M','U10','RNET','P','THICK','U300',
            'SWE']
varnames = ['SWE']
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,tasfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,tascit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'CIT','surface')
    lat,lon,time,lev,tasfic = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,tasfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--CIT}',
               r'\textbf{HIT--CIT}',r'\textbf{FIC--CIT}',
               r'\textbf{FICT--FIT}',r'\textbf{FICT--HIT}']
    runs = [tashit,tasfit,tascit,tasfic,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    tas_on = np.empty((5,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
    tas_dj = np.empty((5,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
    tas_fm = np.empty((5,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
    for i in range(len(runs)):
        tas_on[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1)    
        tas_dj[i],tas_dj[i] = UT.calcDecJan(runs[i],runs[i],lat,lon,'surface',1)    
        tas_fm[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
    
    ### Compute comparisons for FM - taken ensemble average
    diff_FITHIT = np.nanmean(tas_on[1] - tas_on[0],axis=0)
    diff_FITCIT = np.nanmean(tas_on[1] - tas_on[2],axis=0)
    diff_HITCIT = np.nanmean(tas_on[0] - tas_on[2],axis=0)
    diff_FICCIT = np.nanmean(tas_on[3] - tas_on[2],axis=0)
    diff_FICTFIT = np.nanmean(tas_on[4] - tas_on[1],axis=0)
    diff_FICTHIT = np.nanmean(tas_on[4] - tas_on[0],axis=0)
    diffruns_on = np.asarray([diff_FITHIT,diff_FITCIT,diff_HITCIT,diff_FICCIT,
                              diff_FICTFIT,diff_FICTHIT])
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(tas_on[1],tas_on[0])
    stat_FITCIT,pvalue_FITCIT = UT.calc_indttest(tas_on[1],tas_on[2])
    stat_HITCIT,pvalue_HITCIT = UT.calc_indttest(tas_on[0],tas_on[2])
    stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(tas_on[3],tas_on[2])
    stat_FICTFIT,pvalue_FICTFIT = UT.calc_indttest(tas_on[4],tas_on[1])
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(tas_on[4],tas_on[0])
    pruns_on = np.asarray([pvalue_FITHIT,pvalue_FITCIT,pvalue_HITCIT,pvalue_FICCIT,
                           pvalue_FICTFIT,pvalue_FICTHIT])
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot various variables for ON
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'Z500':
        limit = np.arange(-60,60.1,1)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z50':
        limit = np.arange(-60,60.1,1)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z30':
        limit = np.arange(-60,60.1,1)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'SLP':
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,3)
    elif varnames[v] == 'U10' or varnames[v] == 'U300':
        limit = np.arange(-10,10.1,1)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'RNET':    
        limit = np.arange(-50,50.1,1)
        barlim = np.arange(-50,51,25)
    elif varnames[v] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,1) 
    elif varnames[v] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    
    fig = plt.figure()
    for i in range(len(diffruns_on)):
        var = diffruns_on[i]
        pvar = pruns_on[i]
        
        if varnames[v] == 'RNET':
            var = var*-1.
        
        ax1 = plt.subplot(2,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                     linewidths=0.4)
        
        if varnames[v] == 'T2M':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
            cmap = ncm.cmap('nrl_sirkes')           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z50':
            cmap = ncm.cmap('nrl_sirkes')     
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z30':
            cmap = ncm.cmap('nrl_sirkes')   
            cs.set_cmap(cmap)  
        elif varnames[v] == 'SLP':
            cmap = ncm.cmap('nrl_sirkes')           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U300':
            cmap = ncm.cmap('temp_diff_18lev')           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'RNET':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'P':
            cmap = ncm.cmap('precip4_diff_19lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'THICK':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)
        elif varnames[v] == 'SWE':
            cmap = cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
            
        if varnames[v] == 'RNET':
            m.drawcoastlines(color='darkgray',linewidth=0.3)
            m.fillcontinents(color='dimgrey')
        else:
            m.drawcoastlines(color='dimgray',linewidth=0.8)
            
        ### Add experiment text to subplot
        if i >= 4:
            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(0.88,0.885),
                         textcoords='axes fraction',color='k',
                         fontsize=11,rotation=319,ha='center',va='center')
        else:
            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(0.865,0.90),
                         textcoords='axes fraction',color='k',fontsize=11,
                         rotation=320,ha='center',va='center')
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    if varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z50':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z30':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'RNET':
        cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + 'allExperiments_ON_%s.png' % varnames[v],
                dpi=300)

print('Completed: Script done!')

