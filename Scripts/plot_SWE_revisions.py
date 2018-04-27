"""
Plot comparisons to understand contributions of SWE to cooling over Eurasia.
This was requested in paper revisions for GRL.

Notes
-----
    Author : Zachary Labe
    Date   : 27 April 2018
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
print('\n' '----Plotting SWE analysis - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['SWE']
runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
experiments = [r'\textbf{$\Delta$SIT}',r'\textbf{$\Delta$SIC}',
               r'\textbf{$\Delta$NET}']
period = 'DJF'
for v in range(len(varnames)):
    ### Call function for SWE data from each run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'CIT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    
    ### Concatonate runs
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Separate per periods (ON,DJ,FM)
    var_djf = np.empty((5,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):
        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',1)   
    
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
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot T2M
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    fig = plt.figure()
    for i in range(len(experiments)):
        var = diffruns_djf[i]
        pvar = pruns_djf[i]
        
        ### Set limits for contours and colorbars
        if varnames[v] == 'T2M':
            limit = np.arange(-10,10.1,0.5)
            barlim = np.arange(-10,11,5)
        elif varnames[v] == 'Z500':
            limit = np.arange(-60,60.1,1)
            barlim = np.arange(-60,61,30) 
        elif varnames[v] == 'Z30':
            limit = np.arange(-100,100.1,5)
            barlim = np.arange(-100,101,50) 
        elif varnames[v] == 'SLP':
            limit = np.arange(-6,6.1,0.5)
            barlim = np.arange(-6,7,3)
        elif varnames[v] == 'U10' or varnames[v] == 'U300':
            limit = np.arange(-10,10.1,1)
            barlim = np.arange(-10,11,5)
        elif varnames[v] == 'SWE':
            limit = np.arange(-25,25.1,1)
            barlim = np.arange(-25,26,25)
        elif varnames[v] == 'P':
            limit = np.arange(-2,2.1,0.05)
            barlim = np.arange(-2,3,1) 
        elif varnames[v] == 'THICK':
            limit = np.arange(-60,60.1,3)
            barlim = np.arange(-60,61,30)
        elif varnames[v] == 'EGR':
            limit = np.arange(-0.2,0.21,0.02)
            barlim = np.arange(-0.2,0.3,0.2)
        elif varnames[v] == 'RNET':    
            limit = np.arange(-50,50.1,1)
            barlim = np.arange(-50,51,25)
            
        ax1 = plt.subplot(1,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                    area_thresh=10000.)
        
        if varnames[v] == 'RNET':
            var = var*-1.
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
        if varnames[v] == 'Z30': # the interval is 250 m 
            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                            colors='k',linewidths=1.5,zorder=10)
        if varnames[v] == 'RNET':
            m.drawcoastlines(color='darkgray',linewidth=0.3)
            m.fillcontinents(color='dimgrey')
        else:
            m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'T2M':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
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
            cs.set_cmap(cmap) 
        elif varnames[v] == 'SWE':
            cmap = cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
        elif varnames[v] == 'P':
            cmap = ncm.cmap('precip4_diff_19lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'THICK':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)  
        elif varnames[v] == 'RNET':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap) 
        
        ### Add experiment text to subplot
        ax1.annotate(r'\textbf{%s}' % experiments[i],xy=(0,0),xytext=(0.5,1.1),
             textcoords='axes fraction',color='dimgrey',fontsize=18,
             rotation=0,ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    if varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z30':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
    elif varnames[v] == 'RNET':
        cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=11,color='dimgray') 
    
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))    
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    
    plt.savefig(directoryfigure + 'SWEexperiments_%s_%s.png' % (varnames[v],
                                                                     period),
                                                                     dpi=300)
    
print('Completed: Script done!')