"""
Plots DJF for climatological wave number X for WACCM4 experiments

Notes
-----
    Author : Zachary Labe
    Date   : 12 November 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT

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
print('\n' '----Plotting Climo Wave - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnames = ['Z300']
for v in range(len(varnames)):
    ### Call function for geopotential height data from reach run
    lat,lon1,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon1,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon1,time,lev,varcit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'CIT','surface')
    lat,lon1,time,lev,varfic = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIC','surface')
    lat,lon1,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIC--CIT}',
                   r'\textbf{FICT--HIT}']
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Separate per periods (DJF)
    var_djf = np.empty((5,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):
        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon1,'surface',1)    
    
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
    
    ### Read in wave number 
    lat,lon,time,lev,wavef = MO.readExperi(directorydata,
                                        '%sxwave1' % varnames[v],'FIT',
                                        'surface')
    lat,lon,time,lev,wavefc = MO.readExperi(directorydata,
                                        '%sxwave1' % varnames[v],'FIC',
                                        'surface')
    lat,lon,time,lev,wavefict = MO.readExperi(directorydata,
                                    '%sxwave1' % varnames[v],'FICT',
                                    'surface')
    
    wavef_djf,wavef_djf = UT.calcDecJanFeb(wavef,wavef,lat,lon,'surface',1) 
    wavefc_djf,wavefc_djf = UT.calcDecJanFeb(wavefc,wavefc,lat,lon,'surface',1) 
    wavefict_djf,wavefict_djf = UT.calcDecJanFeb(wavefict,wavefict,lat,lon,'surface',1) 
    
    climowavef = np.nanmean(wavef_djf,axis=0)
    climowavefc = np.nanmean(wavefc_djf,axis=0)
    climowavefict = np.nanmean(wavefict_djf,axis=0)
    wavelist = [climowavef,climowavefc,climowavefict]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot U
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    fig = plt.figure()
    for i in range(len(experiments)):
        var = diffruns_djf[i]
        pvar = pruns_djf[i]
        climo = wavelist[i]
        
        limit = np.arange(-60,61,5)
        barlim = np.arange(-60,61,30)
            
        ax1 = plt.subplot(1,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon1)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon1)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
#        climo,lons_cyclic = addcyclic(climo,lon)
#        climo,lons_cyclic = shiftgrid(180.,climo,lons_cyclic,start=False)
#        lon2dc, lat2dc = np.meshgrid(lons_cyclic, lat)
        lon2c,lat2c = np.meshgrid(lon, lat)
                  
        m.drawmapboundary(fill_color='white',color='w',linewidth=0.7)
        m.drawcoastlines(color='dimgray',linewidth=0.65)
        
        cs = m.contourf(x,y,var,limit,extend='both',alpha=0.7,antiliased=True)
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
        cs2 = m.contour(lon2c,lat2c,climo,np.arange(-200,201,60),
                        colors='k',linewidths=1.5,latlon=True,zorder=10)
        
        cmap = ncm.cmap('nrl_sirkes')   
        cs.set_cmap(cmap)  
            
        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        ### Add experiment text to subplot
        ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(0.865,0.90),
                     textcoords='axes fraction',color='k',fontsize=11,
                     rotation=320,ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.23,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=True)
    cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    
    plt.subplots_adjust(wspace=0.01)
    
    plt.savefig(directoryfigure + '%s_climowave1.png' % varnames[v],
                dpi=300)
    
print('Completed: Script done!')

