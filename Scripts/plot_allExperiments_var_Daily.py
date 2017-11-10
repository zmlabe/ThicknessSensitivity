"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data for all variables.

Notes
-----
    Author : Zachary Labe
    Date   : 6 September 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput as DO
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/Daily/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Daily Geopotential - %s----' % titletime)

#### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = ['GEOP','TEMP','U','V']
runnames = [r'HIT',r'FIT',r'CIT',r'FIC',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--CIT}',
               r'\textbf{HIT--CIT}',r'\textbf{FIC--CIT}',
               r'\textbf{FICT--FIT}',r'\textbf{FICT--HIT}']

### Call functions for variable profile data for polar cap
for v in range(len(varnames)):
    lat,lon,time,lev,varhit = DO.readMeanExperi(directorydata,
                                                '%s' % varnames[v],
                                                'HIT','profile')
    lat,lon,time,lev,varfit = DO.readMeanExperi(directorydata,
                                                '%s' % varnames[v],
                                                'FIT','profile')
    lat,lon,time,lev,varcit = DO.readMeanExperi(directorydata,
                                                '%s' % varnames[v],
                                                'CIT','profile')
    lat,lon,time,lev,varfic = DO.readMeanExperi(directorydata,
                                                '%s' % varnames[v],
                                                'FIC','profile')
    lat,lon,time,lev,varfict = DO.readMeanExperi(directorydata,
                                                '%s' % varnames[v],
                                                'FICT','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Compare experiments
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--CIT}',
                   r'\textbf{HIT--CIT}',r'\textbf{FIC--CIT}',
                   r'\textbf{FICT--FIT}',r'\textbf{FICT--HIT}']
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Compute comparisons for experiments - take ensemble average
    diff_FITHIT = np.nanmean(varfit - varhit,axis=0)
    diff_FITCIT = np.nanmean(varfit - varcit,axis=0)
    diff_HITCIT = np.nanmean(varhit - varcit,axis=0)
    diff_FICCIT = np.nanmean(varfic - varcit,axis=0)
    diff_FICTFIT = np.nanmean(varfict - varfit,axis=0)
    diff_FICTHIT = np.nanmean(varfict - varhit,axis=0)
    diffruns = np.asarray([diff_FITHIT,diff_FITCIT,diff_HITCIT,diff_FICCIT,
                              diff_FICTFIT,diff_FICTHIT])
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(varfit,varhit)
    stat_FITCIT,pvalue_FITCIT = UT.calc_indttest(varfit,varcit)
    stat_HITCIT,pvalue_HITCIT = UT.calc_indttest(varhit,varcit)
    stat_FICCIT,pvalue_FICCIT = UT.calc_indttest(varfic,varcit)
    stat_FICTFIT,pvalue_FICTFIT = UT.calc_indttest(varfict,varfit)
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(varfict,varhit)
    pruns = np.asarray([pvalue_FITHIT,pvalue_FITCIT,pvalue_HITCIT,
                           pvalue_FICCIT,pvalue_FICTFIT,pvalue_FICTHIT])
                                                 
    ############################################################################
    ############################################################################
    ############################################################################
    ##### Plot daily profile with height for selected variable
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'GEOP':
        limit = np.arange(-150,150.1,15)
        barlim = np.arange(-150,151,75)
    elif varnames[v] == 'TEMP':
        limit = np.arange(-3,3.1,0.2)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'U':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'V':
        limit = np.arange(-0.3,0.305,0.01)
        barlim = np.arange(-0.3,0.31,0.15)
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    timeq = np.arange(0,212,1)
    timeqq,levq = np.meshgrid(timeq,lev)
    
    fig = plt.figure()
    for i in range(len(experiments)):
        ax1 = plt.subplot(2,3,i+1)
        
        var = diffruns[i]
        pvar = pruns[i]
        
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
        
        cs = plt.contourf(timeq,lev,var.transpose(),limit,extend='both')
                          
        plt.contourf(timeqq,levq,pvar.transpose(),colors='None',
                     hatches=['////'])                  
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
        plt.xlim([30,210])
        plt.ylim([1000,10])
        
        if varnames[v] == 'GEOP':
            cmap = ncm.cmap('temp_diff_18lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'TEMP':
            cmap = ncm.cmap('NCV_blu_red')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'U':
            cmap = ncm.cmap('temp_diff_18lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'V':
            cmap = ncm.cmap('temp_diff_18lev')            
            cs.set_cmap(cmap) 

    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray',labelpad=1)
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    
    ### Add experiments
    plt.annotate(r'%s' % experiments[0],xy=(0,0),xytext=(0.233,0.903),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    plt.annotate(r'%s' % experiments[1],xy=(0,0),xytext=(0.515,0.903),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    plt.annotate(r'%s' % experiments[2],xy=(0,0),xytext=(0.795,0.903),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    plt.annotate(r'%s' % experiments[3],xy=(0,0),xytext=(0.233,0.505),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    plt.annotate(r'%s' % experiments[4],xy=(0,0),xytext=(0.515,0.505),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    plt.annotate(r'%s' % experiments[5],xy=(0,0),xytext=(0.795,0.505),
                 textcoords='figure fraction',color='dimgrey',
                 fontsize=16,ha='center',va='center')
    
    ### Add y-label
    plt.annotate(r'\textbf{Pressure (hPa)}',xy=(0,0),xytext=(0.05,0.54),
             textcoords='figure fraction',color='k',
             fontsize=16,ha='center',va='center',rotation=90)
    
    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(hspace=0.35)
    plt.subplots_adjust(bottom=0.19)
    
    plt.savefig(directoryfigure + 'allExperiments_%s_daily.png' % varnames[v],
                dpi=300)
    print('Completed: Script done!')                                