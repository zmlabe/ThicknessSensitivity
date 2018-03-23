"""
Plots SEP-MAR for fpol-cit, fsub-cit for various parameters in the WACCM4
experiments.

Notes
-----
    Author : Zachary Labe
    Date   : 20 February 2018
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
directoryfigure = '/home/zlabe/Desktop/verticalMonthly/'
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

varnames = ['U','TEMP','GEOP','EGR']
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'CIT','profile')
    lat,lon,time,lev,varfpol = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FPOL','profile')
    lat,lon,time,lev,varfsub = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FSUB','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'CIT',r'FPOL',r'FSUB']
    experiments = [r'\textbf{FPOL--CIT}',r'\textbf{FSUB--CIT}']
    runs = [varcit,varfpol,varfsub]
    
    ### Separate per months
    varmo_fpol = np.append(varfpol[:,9:,:,:,:],varfpol[:,0:3,:,:,:],
             axis=1)
    varmo_cit = np.append(varcit[:,9:,:,:,:],varcit[:,0:3,:,:,:],
             axis=1)
    varmo_fsub = np.append(varfsub[:,9:,:,:,:],varfsub[:,0:3,:,:,:],
              axis=1)
    
    ### Compute comparisons for FM - taken ensemble average
    diff_fpolcit = np.nanmean(varmo_fpol - varmo_cit,axis=0)
    diff_fsubcit = np.nanmean(varmo_fsub - varmo_cit,axis=0)
    diffruns_djf = [diff_fpolcit,diff_fsubcit]
        
    ### Calculate zonal mean
    zdiff_fpolcit = np.nanmean(diff_fpolcit,axis=3)
    zdiff_fsubcit = np.nanmean(diff_fsubcit,axis=3)
    zdiffruns = np.append(zdiff_fpolcit,zdiff_fsubcit,axis=0)
    
    ## Calculate climo
    zclimo_cit1 = np.nanmean(varcit,axis=0)
    zclimo_cit2 = np.nanmean(zclimo_cit1,axis=3)
    zclimoq = np.append(zclimo_cit2[9:,:,:],zclimo_cit2[0:3,:,:],axis=0)
    zclimo = np.append(zclimoq,zclimoq,axis=0)
    
    ### Calculate significance for each month
    stat_fpolcit = np.empty((12,len(lev),len(lat)))
    stat_fsubcit = np.empty((12,len(lev),len(lat)))
    pvalue_fpolcit = np.empty((12,len(lev),len(lat)))
    pvalue_fsubcit = np.empty((12,len(lev),len(lat)))
    for i in range(12):
        stat_fpolcit[i],pvalue_fpolcit[i] = UT.calc_indttest(np.nanmean(varfpol[:,i,:,:,:],axis=3),
                                                     np.nanmean(varcit[:,i,:,:,:],axis=3))
        stat_fsubcit[i],pvalue_fsubcit[i] = UT.calc_indttest(np.nanmean(varfsub[:,i,:,:,:],axis=3),
                                                     np.nanmean(varcit[:,i,:,:,:],axis=3))
    pruns_fpolcit = np.append(pvalue_fpolcit[9:],pvalue_fpolcit[0:3],axis=0)
    pruns_fsubcit = np.append(pvalue_fsubcit[9:],pvalue_fsubcit[0:3],axis=0)
    
    pruns = np.append(pruns_fpolcit,pruns_fsubcit,axis=0)
    
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
    elif varnames[v] == 'TEMP':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,1)
    elif varnames[v] == 'GEOP':
        limit = np.arange(-60,61,1)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.08,0.081,0.005)
        barlim = np.arange(-0.08,0.09,0.04)
        
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    latq,levq = np.meshgrid(lat,lev)
    
    fig = plt.figure()
    for i in range(12):
        ax1 = plt.subplot(2,6,i+1)
        
        clmq = i
    
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
                              linewidths=0.3,colors='dimgrey')
        if varnames[v] == 'EGR': 
            cs2 = plt.contour(lat,lev,zclimo[i],np.arange(-0.8,0.9,0.2),
                              linewidths=0.3,colors='dimgrey')
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
        elif varnames[v] == 'TEMP':
            cmap = ncm.cmap('NCV_blu_red')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'GEOP':
            cmap = ncm.cmap('temp_diff_18lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.curl           
            cs.set_cmap(cmap) 
    
        labelmonths = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.08),xycoords='axes fraction',
                        fontsize=17,color='dimgrey',rotation=0,
                        ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    elif varnames[v] == 'TEMP':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
    elif varnames[v] == 'GEOP':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    
            
    plt.annotate(r'\textbf{FPOL--HIT2}',
            xy=(0, 0),xytext=(0.045,0.73),xycoords='figure fraction',
            fontsize=17,color='k',rotation=90,
            ha='center',va='center')        
    plt.annotate(r'\textbf{FSUB--HIT2}',
        xy=(0, 0),xytext=(0.045,0.36),xycoords='figure fraction',
        fontsize=17,color='k',rotation=90,
        ha='center',va='center')  
#    plt.annotate(r'\textbf{Latitude [$^\circ$N]}',
#        xy=(0, 0),xytext=(0.5,0.1),xycoords='figure fraction',
#        fontsize=17,color='k',rotation=0,
#        ha='center',va='center') 
    plt.subplots_adjust(hspace=0.2)
    plt.subplots_adjust(bottom=0.21)
    
    plt.savefig(directoryfigure + 'mo_vertical_region_%s.png' % varnames[v],dpi=300)
print('Completed: Script done!')

