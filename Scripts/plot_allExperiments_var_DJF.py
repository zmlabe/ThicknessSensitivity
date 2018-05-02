"""
Plots DJF for various parameters for WACCM4 experiments

Notes
-----
    Author : Zachary Labe
    Date   : 7 November 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/vertical/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting temperature - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnames = ['U','TEMP','GEOP']
varnames = ['TEMP']
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','profile')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','profile')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'CIT','profile')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIC','profile')
    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--CIT}',
                   r'\textbf{HIT--CIT}',r'\textbf{FIC--CIT}',
                   r'\textbf{FICT--FIT}',r'\textbf{FICT--HIT}']
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Separate per periods (DJF)
    var_djf = np.empty((5,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3],
                        varhit.shape[4]))
    for i in range(len(runs)):
        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'profile',17)    
    
    ### Compute comparisons for FM - taken ensemble average
    diff_FITHIT = np.nanmean(var_djf[1] - var_djf[0],axis=0)
    diff_FITCIT = np.nanmean(var_djf[1] - var_djf[2],axis=0)
    diff_HITCIT = np.nanmean(var_djf[0] - var_djf[2],axis=0)
    diff_FICCIT = np.nanmean(var_djf[3] - var_djf[2],axis=0)
    diff_FICTFIT = np.nanmean(var_djf[4] - var_djf[1],axis=0)
    diff_FICTHIT = np.nanmean(var_djf[4] - var_djf[0],axis=0)
    diffruns_djf = np.asarray([diff_FITHIT,diff_FITCIT,diff_HITCIT,diff_FICCIT,
                              diff_FICTFIT,diff_FICTHIT])
        
    ### Calculate zonal mean
    zdiff_FITHIT = np.nanmean(diff_FITHIT,axis=2)
    zdiff_FITCIT = np.nanmean(diff_FITCIT,axis=2)
    zdiff_HITCIT = np.nanmean(diff_HITCIT,axis=2)
    zdiff_FICCIT = np.nanmean(diff_FICCIT,axis=2)
    zdiff_FICTFIT = np.nanmean(diff_FICTFIT,axis=2)
    zdiff_FICTHIT = np.nanmean(diff_FICTHIT,axis=2)
    zdiffruns_djf = np.asarray([zdiff_FITHIT,zdiff_FITCIT,zdiff_HITCIT,
                                zdiff_FICCIT,zdiff_FICTFIT,zdiff_FICTHIT])
    
    ## Calculate climo
    zclimo_hit = np.apply_over_axes(np.nanmean,var_djf[0],(0,3)).squeeze()
    zclimo_fit = np.apply_over_axes(np.nanmean,var_djf[1],(0,3)).squeeze()
    zclimo_cit = np.apply_over_axes(np.nanmean,var_djf[2],(0,3)).squeeze()
    zclimo_fic = np.apply_over_axes(np.nanmean,var_djf[3],(0,3)).squeeze()
    zclimo_fict = np.apply_over_axes(np.nanmean,var_djf[4],(0,3)).squeeze()
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(np.nanmean(var_djf[1],axis=3),
                                                 np.nanmean(var_djf[0],axis=3))
    stat_FITCIT,pvalue_FITCIT = UT.calc_indttest(np.nanmean(var_djf[1],axis=3),
                                                 np.nanmean(var_djf[2],axis=3))
    stat_HITCIT,pvalue_HITCIT = UT.calc_indttest(np.nanmean(var_djf[0],axis=3),
                                                 np.nanmean(var_djf[2],axis=3))
    stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(np.nanmean(var_djf[3],axis=3),
                                                 np.nanmean(var_djf[2],axis=3))
    stat_FICTFIT,pvalue_FICTFIT = UT.calc_indttest(np.nanmean(var_djf[4],axis=3),
                                                 np.nanmean(var_djf[1],axis=3))
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(np.nanmean(var_djf[4],axis=3),
                                                 np.nanmean(var_djf[0],axis=3))
    pruns_djf = np.asarray([pvalue_FITHIT,pvalue_FITCIT,pvalue_HITCIT,
                           pvalue_FICCIT,pvalue_FICTFIT,pvalue_FICTHIT])
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot U
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) z
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'U':
        limit = np.arange(-3,3.1,0.1)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'TEMP':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,1)
    elif varnames[v] == 'GEOP':
        limit = np.arange(-60,61,1)
        barlim = np.arange(-60,61,30)
        
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    latq,levq = np.meshgrid(lat,lev)
    
    fig = plt.figure()
    ax1 = plt.subplot(131)
    
    ax1.spines['top'].set_color('dimgrey')
    ax1.spines['right'].set_color('dimgrey')
    ax1.spines['bottom'].set_color('dimgrey')
    ax1.spines['left'].set_color('dimgrey')
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
    ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')    
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    
    cs = plt.contourf(lat,lev,zdiff_FITHIT,limit,extend='both')
    
    if varnames[v] == 'U': 
        cs2 = plt.contour(lat,lev,zclimo_hit,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')
    plt.contourf(latq,levq,pvalue_FITHIT,colors='None',hatches=['////'],
                 linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    
    plt.ylabel(r'\textbf{Pressure (hPa)',fontsize=13)
    
    if varnames[v] == 'U':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'TEMP':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'GEOP':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax1.annotate(r'\textbf{$\Delta$SIT}',
                xy=(0, 0),xytext=(0.5,1.06),xycoords='axes fraction',
                fontsize=19,color='dimgrey',rotation=0,ha='center',
                va='center')
    
    #############################################################################
    ax2 = plt.subplot(132)
    
    ax2.spines['top'].set_color('dimgrey')
    ax2.spines['right'].set_color('dimgrey')
    ax2.spines['bottom'].set_color('dimgrey')
    ax2.spines['left'].set_color('dimgrey')
    ax2.spines['left'].set_linewidth(2)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['right'].set_linewidth(2)
    ax2.spines['top'].set_linewidth(2)
    ax2.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
    ax2.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')    
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    
    cs = plt.contourf(lat,lev,zdiff_FICCIT,limit,extend='both')
    if varnames[v] == 'U':
        cs2 = plt.contour(lat,lev,zclimo_cit,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')
    plt.contourf(latq,levq,pvalue_FICCIT,colors='None',hatches=['////'],
                 linewidth=5) 
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax2.annotate(r'\textbf{$\Delta$SIC}',
                xy=(0, 0),xytext=(0.5,1.06),xycoords='axes fraction',
                fontsize=19,color='dimgrey',rotation=0,ha='center',
                va='center')
    
    if varnames[v] == 'U':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'TEMP':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'GEOP':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    
    ############################################################################
    ax3 = plt.subplot(133)
    
    ax3.spines['top'].set_color('dimgrey')
    ax3.spines['right'].set_color('dimgrey')
    ax3.spines['bottom'].set_color('dimgrey')
    ax3.spines['left'].set_color('dimgrey')
    ax3.spines['left'].set_linewidth(2)
    ax3.spines['bottom'].set_linewidth(2)
    ax3.spines['right'].set_linewidth(2)
    ax3.spines['top'].set_linewidth(2)
    ax3.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
    ax3.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')    
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    
    cs = plt.contourf(lat,lev,zdiff_FICTHIT,limit,extend='both')
    if varnames[v] == 'U':        
        cs2 = plt.contour(lat,lev,zclimo_hit,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')
    plt.contourf(latq,levq,pvalue_FICTHIT,colors='None',hatches=['////'],
                 linewidth=5) 
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    
    if varnames[v] == 'U':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'TEMP':
        cmap = ncm.cmap('NCV_blu_red')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'GEOP':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax3.annotate(r'\textbf{$\Delta$NET}',
                xy=(0, 0),xytext=(0.5,1.06),xycoords='axes fraction',
                fontsize=19,color='dimgrey',rotation=0,ha='center',
                va='center')
    
    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    elif varnames[v] == 'TEMP':
        cbar.set_label(r'\textbf{[T]$^\circ$C}',fontsize=11,color='dimgray',
                                 labelpad=0)
    elif varnames[v] == 'GEOP':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    
    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(bottom=0.21)
    
    plt.savefig(directoryfigure + 'DJF_vertical_%s.png' % varnames[v],dpi=300)
print('Completed: Script done!')

