"""
Plot figure 4 in manuscript for zonal wind responses to sea ice loss in WACCM4
experiments [FIT-HIT, FICT-HIT].Time period includes all months.

Notes
-----
    Author : Zachary Labe
    Date   : 3 February 2018
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
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Fig. 4 - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

varnames = ['U']
for v in range(len(varnames)):
#    ### Call function for surface temperature data from reach run
#    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
#                                            '%s' % varnames[v],'HIT','profile')
#    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
#                                            '%s' % varnames[v],'FIT','profile')
#    lat,lon,time,lev,varfict = MO.readExperi(directorydata,
#                                             '%s' % varnames[v],'FICT','profile')
#    
#    ### Create 2d array of latitude and longitude
#    lon2,lat2 = np.meshgrid(lon,lat)
#    
#    ### Concatonate runs
#    runnames = [r'HIT',r'FIT',r'FICT']
#    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FICT--HIT}']
#    runs = [varhit,varfit,varfict]
#    
#    ### Separate per months
#    varmo_fit = np.append(varfit[:,9:,:,:,:],varfit[:,0:3,:,:,:],
#             axis=1)
#    varmo_hit = np.append(varhit[:,9:,:,:,:],varhit[:,0:3,:,:,:],
#             axis=1)
#    varmo_fict = np.append(varfict[:,9:,:,:,:],varfict[:,0:3,:,:,:],
#              axis=1)
#    
#    ### Compute comparisons for FM - taken ensemble average
#    diff_FITHIT = np.nanmean(varmo_fit - varmo_hit,axis=0)
#    diff_FICTHIT = np.nanmean(varmo_fict - varmo_hit,axis=0)
#    diffruns_djf = [diff_FITHIT,diff_FICTHIT]
#        
#    ### Calculate zonal mean
#    zdiff_FITHIT = np.nanmean(diff_FITHIT,axis=3)
#    zdiff_FICTHIT = np.nanmean(diff_FICTHIT,axis=3)
#    zdiffruns = np.append(zdiff_FITHIT,zdiff_FICTHIT,axis=0)
#    
#    ## Calculate climo
#    zclimo_hit1 = np.nanmean(varhit,axis=0)
#    zclimo_hit2 = np.nanmean(zclimo_hit1,axis=3)
#    zclimoq = np.append(zclimo_hit2[9:,:,:],zclimo_hit2[0:3,:,:],axis=0)
#    zclimo = np.append(zclimoq,zclimoq,axis=0)
#    
#    ### Calculate significance for each month
#    stat_FITHIT = np.empty((12,len(lev),len(lat)))
#    stat_FICTHIT = np.empty((12,len(lev),len(lat)))
#    pvalue_FITHIT = np.empty((12,len(lev),len(lat)))
#    pvalue_FICTHIT = np.empty((12,len(lev),len(lat)))
#    for i in range(12):
#        stat_FITHIT[i],pvalue_FITHIT[i] = UT.calc_indttest(np.nanmean(varfit[:,i,:,:,:],axis=3),
#                                                     np.nanmean(varhit[:,i,:,:,:],axis=3))
#        stat_FICTHIT[i],pvalue_FICTHIT[i] = UT.calc_indttest(np.nanmean(varfict[:,i,:,:,:],axis=3),
#                                                     np.nanmean(varhit[:,i,:,:,:],axis=3))
#    pruns_FITHIT = np.append(pvalue_FITHIT[9:],pvalue_FITHIT[0:3],axis=0)
#    pruns_FICTHIT = np.append(pvalue_FICTHIT[9:],pvalue_FICTHIT[0:3],axis=0)
#    
#    pruns = np.append(pruns_FITHIT,pruns_FICTHIT,axis=0)
    
    ###########################################################################
    ###########################################################################
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
    for i in range(12):
        ax1 = plt.subplot(2,6,i+1)
    
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
                
        cs = plt.contourf(lat,lev,zdiffruns[i],limit,extend='both')
        
        if varnames[v] == 'U': 
            cs2 = plt.contour(lat,lev,zclimo[i],np.arange(-20,101,5),
                              linewidths=0.5,colors='dimgrey')
            
        plt.contourf(latq,levq,pruns[i],colors='None',hatches=['////'],
                     linewidth=5)   
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=7)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=7)
        plt.minorticks_off()
        
        plt.xlim([0,90])
        plt.ylim([1000,10])
        
        if i==1 or i==2 or i==3 or i==4 or i==5 or i==7 or i==8 or i==9 or i==10 or i==11:
            ax1.tick_params(labelleft='off') 
        if i < 6:
            ax1.tick_params(labelbottom='off') 
            
        if varnames[v] == 'U':
            cmap = ncm.cmap('NCV_blu_red')            
            cs.set_cmap(cmap) 
    
        labelmonths = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.08),xycoords='axes fraction',
                        fontsize=17,color='dimgrey',rotation=0,
                        ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{[U] m/s}',fontsize=11,color='dimgrey',
                                 labelpad=0)
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
          
    plt.annotate(r'\textbf{$\Delta$SIT}',
            xy=(0, 0),xytext=(0.055,0.73),xycoords='figure fraction',
            fontsize=21,color='k',rotation=90,
            ha='center',va='center')        
    plt.annotate(r'\textbf{$\Delta$NET}',
        xy=(0, 0),xytext=(0.055,0.36),xycoords='figure fraction',
        fontsize=21,color='k',rotation=90,
        ha='center',va='center')  
    plt.annotate(r'\textbf{Latitude ($^{\circ}$N)',
        xy=(0, 0),xytext=(0.515,0.15),xycoords='figure fraction',
        fontsize=8,color='k',rotation=0,
        ha='center',va='center')  

    plt.subplots_adjust(hspace=0.2)
    plt.subplots_adjust(bottom=0.21)
    
plt.savefig(directoryfigure + 'Fig4.png',dpi=600)
print('Completed: Script done!')

