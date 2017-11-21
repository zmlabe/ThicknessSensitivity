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
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT']
    experiments = [r'\textbf{FIT--HIT}']
    runs = [varhit,varfit]
    
    ### Separate per periods (Feb,Mar)
    varh_f = runs[0][:,1,:,:]
    varh_m = runs[0][:,2,:,:]
    
    varf_f = runs[1][:,1,:,:]
    varf_m = runs[1][:,2,:,:]
    
    ### Compute comparisons for FM - taken ensemble average
    diff_feb = np.nanmean(varf_f - varh_f,axis=0)
    diff_mar = np.nanmean(varf_m - varh_m,axis=0)
    
    ### Calculate significance for FM
    stat_feb,pvalue_feb = UT.calc_indttest(varfit[:,1,:,:],varhit[:,1,:,:])
    stat_mar,pvalue_mar = UT.calc_indttest(varfit[:,1,:,:],varhit[:,1,:,:])
    
    ### Read in wave number 
    lat,lon,time,lev,waveh= MO.readExperi(directorydata,
                                        '%sxwave1' % varnames[v],'HIT',
                                        'surface')
#    lat,lon,time,lev,wavef = MO.readExperi(directorydata,
#                                        '%sxwave1' % varnames[v],'FIT',
#                                        'surface')
    
    climowaveh_feb = np.nanmean(waveh[:,1,:,:],axis=0)
    climowaveh_mar = np.nanmean(waveh[:,2,:,:],axis=0)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot climatological wave
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    fig = plt.figure()
    ax = plt.subplot(121)
    
    var = diff_feb
    pvar = pvalue_feb
    climo = climowaveh_feb
    
    limit = np.arange(-60,61,5)
    barlim = np.arange(-60,61,30)
        
    m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon1)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    lon2c,lat2c = np.meshgrid(lon, lat)
              
    m.drawmapboundary(fill_color='white',color='w',linewidth=0.7)
    m.drawcoastlines(color='dimgray',linewidth=0.65)
    
    cs = m.contourf(x,y,var,limit,extend='both',alpha=0.7,antiliased=True)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
    cs2 = m.contour(lon2c,lat2c,climo,np.arange(-200,201,50),
                    colors='k',linewidths=1.5,latlon=True,zorder=10)
    
    cmap = ncm.cmap('nrl_sirkes')   
    cs.set_cmap(cmap)  
        
    m.drawcoastlines(color='dimgray',linewidth=0.8)
    
    ### Add experiment text to subplot
    ax.annotate(r'\textbf{FEB}',xy=(0,0),xytext=(0.5,1.1),
                 textcoords='axes fraction',color='dimgray',fontsize=23,
                 rotation=0,ha='center',va='center')
    ax.annotate(r'\textbf{FIT--HIT}',xy=(0,0),xytext=(-0.1,0.5),
             textcoords='axes fraction',color='dimgray',fontsize=23,
             rotation=90,ha='center',va='center')
    
    ###########################################################################
    ax = plt.subplot(122)
    
    var = diff_mar
    pvar = pvalue_mar
    climo = climowaveh_mar
    
    limit = np.arange(-60,61,5)
    barlim = np.arange(-60,61,30)
        
    m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon1)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    lon2c,lat2c = np.meshgrid(lon, lat)
              
    m.drawmapboundary(fill_color='white',color='w',linewidth=0.7)
    m.drawcoastlines(color='dimgray',linewidth=0.65)
    
    cs = m.contourf(x,y,var,limit,extend='both',alpha=0.7,antiliased=True)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
    cs2 = m.contour(lon2c,lat2c,climo,np.arange(-200,201,50),
                    colors='k',linewidths=1.5,latlon=True,zorder=10)
    
    cmap = ncm.cmap('nrl_sirkes')   
    cs.set_cmap(cmap)  
        
    m.drawcoastlines(color='dimgray',linewidth=0.8)
    
    ### Add text
    ax.annotate(r'\textbf{MAR}',xy=(0,0),xytext=(0.5,1.1),
             textcoords='axes fraction',color='dimgray',fontsize=23,
             rotation=0,ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.13,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=True)
    cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgray')
    
    plt.subplots_adjust(wspace=0.01)
    
    plt.savefig(directoryfigure + '%s_climowave1.png' % varnames[v],
                dpi=300)
    
###############################################################################
###############################################################################
###############################################################################
varnames = ['Z30']
for v in range(len(varnames)):
    ### Call function for geopotential height data from reach run
    lat,lon1,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon1,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT']
    experiments = [r'\textbf{FIT--HIT}']
    runs = [varhit,varfit]
    
    ### Separate per periods (Feb,Mar)
    varh_f = runs[0][:,1,:,:]
    varh_m = runs[0][:,2,:,:]
    
    varf_f = runs[1][:,1,:,:]
    varf_m = runs[1][:,2,:,:]
    
    ### Compute comparisons for FM - taken ensemble average
    diff_feb = np.nanmean(varf_f - varh_f,axis=0)
    diff_mar = np.nanmean(varf_m - varh_m,axis=0)
    
    ### Calculate significance for FM
    stat_feb,pvalue_feb = UT.calc_indttest(varfit[:,1,:,:],varhit[:,1,:,:])
    stat_mar,pvalue_mar = UT.calc_indttest(varfit[:,2,:,:],varhit[:,2,:,:])
    
    ### Calculate climatology    
    climowaveh_feb = np.nanmean(varhit[:,1,:,:],axis=0)
    climowaveh_mar = np.nanmean(varhit[:,2,:,:],axis=0)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot Z30
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    fig = plt.figure()
    ax = plt.subplot(121)
    
    var = diff_feb
    pvar = pvalue_feb
    climo = climowaveh_feb
    
    limit = np.arange(-100,101,5)
    barlim = np.arange(-100,101,50)
        
    m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon1)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    climo,lons_cyclic = addcyclic(climowaveh_feb, lon1)
    climo,lons_cyclic = shiftgrid(180.,climo,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='dimgray',linewidth=0.65)
    
    cs = m.contourf(x,y,var,limit,extend='both',alpha=1)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
    cs2 = m.contour(x,y,climo,np.arange(21900,23500,250),
                    colors='k',linewidths=1.5,zorder=10)
    
    cmap = ncm.cmap('nrl_sirkes')   
    cs.set_cmap(cmap)  
        
    m.drawcoastlines(color='dimgray',linewidth=0.8)
    
    ### Add experiment text to subplot
    ax.annotate(r'\textbf{FEB}',xy=(0,0),xytext=(0.5,1.1),
                 textcoords='axes fraction',color='dimgray',fontsize=23,
                 rotation=0,ha='center',va='center')
    ax.annotate(r'\textbf{FIT--HIT}',xy=(0,0),xytext=(-0.1,0.5),
             textcoords='axes fraction',color='dimgray',fontsize=23,
             rotation=90,ha='center',va='center')
    
    ###########################################################################
    ax = plt.subplot(122)
    
    var = diff_mar
    pvar = pvalue_mar
    climo = climowaveh_mar
    
    limit = np.arange(-125,126,5)
    barlim = np.arange(-125,126,125)
        
    m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon1)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon1)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    climo,lons_cyclic = addcyclic(climowaveh_mar, lon1)
    climo,lons_cyclic = shiftgrid(180.,climo,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='dimgray',linewidth=0.65)
    
    cs = m.contourf(x,y,var,limit,extend='both',alpha=1)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
    cs2 = m.contour(x,y,climo,np.arange(21900,23500,250),
                    colors='k',linewidths=1.5,zorder=10)
    
    cmap = ncm.cmap('nrl_sirkes')   
    cs.set_cmap(cmap)  
        
    m.drawcoastlines(color='dimgray',linewidth=0.8)
    
    ### Add text
    ax.annotate(r'\textbf{MAR}',xy=(0,0),xytext=(0.5,1.1),
             textcoords='axes fraction',color='dimgray',fontsize=23,
             rotation=0,ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.13,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=True)
    cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgray')
    
    plt.subplots_adjust(wspace=0.01)
    
    plt.savefig(directoryfigure + '%s_FIT-HIT.png' % varnames[v],
                dpi=300)    
    
print('Completed: Script done!')

