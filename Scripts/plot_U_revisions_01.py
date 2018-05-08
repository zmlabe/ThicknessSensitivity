"""
Plots tests of zonal wind for regions compared with net sea ice effect. Plot
was created during revisions #1.

Notes
-----
    Author : Zachary Labe
    Date   : 1 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
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
print('\n' '----Plotting U tests- %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnames = ['U']
for v in range(len(varnames)):
#    ### Call function for surface temperature data from reach run
#    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
#                                            '%s' % varnames[v],'CIT','profile')
#    lat,lon,time,lev,varfpol = MO.readExperi(directorydata,
#                                            '%s' % varnames[v],'FPOL','profile')
#    lat,lon,time,lev,varfsub = MO.readExperi(directorydata,
#                                             '%s' % varnames[v],'FSUB','profile')
#    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
#                                            '%s' % varnames[v],'FICT','profile')
#    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
#                                             '%s' % varnames[v],'HIT','profile')
#    
#    ### Create 2d array of latitude and longitude
#    lon2,lat2 = np.meshgrid(lon,lat)
#    
#    ### Concatonate runs
#    runnames = [r'CIT',r'FPOL',r'FSUB',r'FICT',r'HIT']
#    experiments = [r'\textbf{FPOL--CIT}',r'\textbf{FSUB--CIT}',
#                   r'\textbf{FICT--HIT']
#    runs = [varcit,varfpol,varfsub,varfict,varhit]
#    
#    ### Separate per periods (DJF)
#    var_djf = np.empty((5,varcit.shape[0]-1,varcit.shape[2],varcit.shape[3],
#                        varcit.shape[4]))
#    for i in range(len(runs)):
#        var_djf[i],var_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
#                                              lon,'profile',17)   
#       
#    ### Compute comparisons for FM - taken ensemble average
#    diff_FSUBCIT = np.nanmean(var_djf[1] - var_djf[0],axis=0)
#    diff_FPOLCIT = np.nanmean(var_djf[2] - var_djf[0],axis=0)
#    diff_FICTHIT = np.nanmean(var_djf[3] - var_djf[4],axis=0)
#    diffruns_djf = [diff_FSUBCIT,diff_FPOLCIT,diff_FICTHIT]
#        
#    ### Calculate zonal mean
#    zdiff_FSUBCIT = np.nanmean(diff_FSUBCIT,axis=2)
#    zdiff_FPOLCIT = np.nanmean(diff_FPOLCIT,axis=2)
#    zdiff_FICTHIT = np.nanmean(diff_FICTHIT,axis=2)
#    zdiffruns_djf = np.asarray([zdiff_FSUBCIT,zdiff_FPOLCIT,zdiff_FICTHIT])
#    
#    ### Calculate climo
#    zclimo_cit = np.apply_over_axes(np.nanmean,var_djf[0],(0,3)).squeeze()
#    zclimo_FSUB = np.apply_over_axes(np.nanmean,var_djf[1],(0,3)).squeeze()
#    zclimo_FPOL = np.apply_over_axes(np.nanmean,var_djf[2],(0,3)).squeeze()
#    zclimo_hit = np.apply_over_axes(np.nanmean,var_djf[4],(0,3)).squeeze()
#    
#    ### Calculate signicance for DJF
#    stat_FSUBCIT,pvalue_FSUBCIT = UT.calc_indttest(np.nanmean(var_djf[1],axis=3),
#                                                 np.nanmean(var_djf[0],axis=3))
#    stat_FPOLCIT,pvalue_FPOLCIT = UT.calc_indttest(np.nanmean(var_djf[2],axis=3),
#                                                 np.nanmean(var_djf[0],axis=3))
#    stat_FICTHIT,pvalue_FICTHIT= UT.calc_indttest(np.nanmean(var_djf[3],axis=3),
#                                                 np.nanmean(var_djf[4],axis=3))
#    pruns_djf = np.asarray([pvalue_FSUBCIT,pvalue_FPOLCIT,pvalue_FICTHIT])
    
    ###########################################################################
    ###########################6################################################
    ###########################################################################
    #### Plot U
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'U':
        limit = np.arange(-3,3.1,0.1)
        barlim = np.arange(-3,4,1)
        
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    latq,levq = np.meshgrid(lat,lev)
    
    fig = plt.figure()
    ax1 = plt.subplot(132)
    
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
    
    
    cs = plt.contourf(lat,lev,zdiff_FSUBCIT+zdiff_FPOLCIT,limit,extend='both')
    
    if varnames[v] == 'U': 
        cs2 = plt.contour(lat,lev,zclimo_cit,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')
#    plt.contourf(latq,levq,pvalue_FICTHIT,colors='None',hatches=['////'],
#                 linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    
    if varnames[v] == 'U' or varnames[v] == 'V':
        cmap = ncm.cmap('NCV_blu_red')          
        cs.set_cmap(cmap) 
    elif varnames[v] == 'TEMP':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'GEOP':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'EGR':
        cmap = cmocean.cm.curl           
        cs.set_cmap(cmap) 
    
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax1.annotate(r'\textbf{$\Delta$$\bf{\sum}$Regions',
                xy=(0, 0),xytext=(0.2,1.02),xycoords='axes fraction',
                fontsize=10,color='dimgrey',rotation=0)

    ###########################################################################
    ax3 = plt.subplot(131)
    
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
    
    if varnames[v] == 'U' or varnames[v] == 'V':
        cmap = ncm.cmap('NCV_blu_red')              
        cs.set_cmap(cmap) 
        
    plt.ylabel(r'\textbf{Pressure (hPa)',fontsize=12)
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax3.annotate(r'\textbf{$\Delta$NET}',
                xy=(0, 0),xytext=(0.3,1.02),xycoords='axes fraction',
                fontsize=18,color='dimgrey',rotation=0)
    
    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{[U] m/s}',fontsize=11,color='dimgray',
                                 labelpad=0)
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.ax.tick_params(labelsize=8)
    
    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(bottom=0.21)
    
    ax1 = plt.subplot(133)
    
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
    
    
    cs = plt.contourf(lat,lev,(zdiff_FICTHIT - (zdiff_FSUBCIT+zdiff_FPOLCIT)),limit,extend='both')
    
    if varnames[v] == 'U': 
        cs2 = plt.contour(lat,lev,zclimo_cit,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')
#    plt.contourf(latq,levq,pvalue_FSUBCIT,colors='None',hatches=['////'],
#                 linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    
    if varnames[v] == 'U' or varnames[v] == 'V':
        cmap = ncm.cmap('NCV_blu_red')          
        cs.set_cmap(cmap) 
    elif varnames[v] == 'TEMP':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'GEOP':
        cmap = ncm.cmap('temp_diff_18lev')            
        cs.set_cmap(cmap) 
    elif varnames[v] == 'EGR':
        cmap = cmocean.cm.curl           
        cs.set_cmap(cmap) 
    
    plt.xlabel(r'\textbf{Latitude ($^{\circ}$N)',fontsize=8,labelpad=0)
    ax1.annotate(r'\textbf{$\Delta$NET -- $\Delta$$\sum$Regions}',
                xy=(0, 0),xytext=(0.0,1.02),xycoords='axes fraction',
                fontsize=10,color='dimgrey',rotation=0)
    
    plt.savefig(directoryfigure + 'U_revisionsTest_%s.png' % varnames[v],dpi=900)
print('Completed: Script done!')

