"""
Compute ratio (%) between SIT and NET responses for DJF

Notes
-----
    Author : Zachary Labe
    Date   : 24 April 2018 [Revisions #1]
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput as MO
import cmocean
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)

currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting DJF SIT-SIC ratio - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

months = [r'DJF']
varnames = ['U10','Z30','U300','Z500','SLP','T2M','RNET']

ratiovar = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'HIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'FICT',r'HIT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FICT--HIT}']
    runs = [varhit,varfit,varfic,varcit]
    
    ### Separate per 2 month periods
    varmo_djf = np.empty((4,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):   
        varmo_djf[i],varmo_djf[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',1)    
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    diff_fithit_djf = np.nanmean(varmo_djf[1] - varmo_djf[0],axis=0)
    diff_ficcit_djf = np.nanmean(varmo_djf[2] - varmo_djf[3],axis=0)
    
    ### Calculate significance 
    stat_FITHITdjf,pvalue_FITHITdjf = UT.calc_indttest(varmo_djf[1],varmo_djf[0])
    stat_FICCITdjf,pvalue_FICCITdjf = UT.calc_indttest(varmo_djf[2],varmo_djf[3])
    
    ### Create mask of significant values
    pvalue_FITHITdjf[np.where(np.isnan(pvalue_FITHITdjf))] = 0.0
    pvalue_FICCITdjf[np.where(np.isnan(pvalue_FICCITdjf))] = 0.0
        
    pvalue_FITHIT = [pvalue_FITHITdjf]
    pvalue_FICCIT = [pvalue_FICCITdjf]
    
    ### Create mask of shared significant values
    mask = np.asarray(pvalue_FITHIT) * np.asarray(pvalue_FICCIT)
    
    ### Slice out lats below different regions
    slicearea = False
    if slicearea == True:
        if varnames[v] == 'U300':
            latq = np.where((lat>=50) & (lat<=70))[0]
            latqq = lat[latq]
        elif varnames[v] == 'Z500':
            latq = np.where((lat>=65) & (lat<=90))[0]
            latqq = lat[latq]
        else:
            latq = np.where(lat>40)[0]
            latqq = lat[latq]
        
    elif slicearea == False:
        ### Slice out lats below 40
        latq = np.where(lat>40)[0]
        latqq = lat[latq]
    else:
        print(ValueError('Wrong region sliced!'))
    
    ### Create 2nd meshgrid with lats > 40N
    lonnew,latnew=np.meshgrid(lon,latqq)
    
    ### Create mask for DJF
    mask = mask[:,latq,:]
    
    ### Keep only values significant in both SIT and SIC responses    
    diff_fithit_djfq = diff_fithit_djf[latq,:] * pvalue_FITHITdjf[latq,:]
    
    diff_ficcit_djfq = diff_ficcit_djf[latq,:] * pvalue_FICCITdjf[latq,:]
    
    fithit = [diff_fithit_djfq]
    ficcit = [diff_ficcit_djfq]
    
    def calc_iceRatio(varx,vary,maske,up,down):
        """
        Compute relative % difference
        """
        print('\n>>> Using calc_iceRatio function!')
        
        ### Mask extremes
        if maske == True:
            print('MASKING EXTREMES!')
            
            varxup = np.nanpercentile(varx,up)
            varxdo = np.nanpercentile(varx,down)
            
            varyup = np.nanpercentile(vary,up)
            varydo = np.nanpercentile(vary,down)
            
            print(varxup,varxdo)
            print(varyup,varydo)
            
            varx[np.where((varx >= varxup) | (varx <= varxdo))] = np.nan
            vary[np.where((vary >= varyup) | (vary <= varydo))] = np.nan
        
        percchange = (abs(varx)/abs(vary)) * 100.
        
        ### Test if real values
        if np.isnan(percchange).all() == True:
            percchange[np.where(np.isnan(percchange))] = 0.0
        if percchange > 500:
            percchange = 0.0
                
        print('*Completed: Finished calc_iceRatio function!')
        return percchange,varx,vary
    
    fithitave = np.empty((3))
    ficcitave = np.empty((3))
    for i in range(len(fithit)):
        fithitave[i] = UT.calc_weightedAve(abs(fithit[i]),latnew)
        ficcitave[i] = UT.calc_weightedAve(abs(ficcit[i]),latnew)

    ratio = []
    for i in range(len(fithit)):
        percchangeq,varx,vary = calc_iceRatio(fithitave[i],ficcitave[i],False,95,5)
        
        ratio.append(percchangeq)
    ratiovar.append(ratio)
meanratiovar = np.asarray(ratiovar).squeeze()
    
#### Save file
np.savetxt(directorydata2 + 'sitNETratio_DJF_R1area.txt',np.round(meanratiovar.transpose(),1),delimiter=',',
           fmt='%3.1f',header='  '.join(varnames)+'\n',
           footer='\n File contains ratio values of relative contributions' \
           '\n between FIT-HIT and FICT-HIT to get the relative \n' \
           ' contributions of SIT and SIC [bimonth, DJF]',newline='\n\n')

###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on')

vals,rat = np.meshgrid(np.arange(2),meanratiovar)

cs = plt.pcolormesh(rat,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=0,vmax=50)

for i in range(rat.shape[0]):
    for j in range(rat.shape[1]):
        plt.text(j+0.5,i+0.5,r'\textbf{%3.1f}' % rat[i,j],fontsize=6,
                 color='k',va='center',ha='center')

cs.set_cmap(cmocean.cm.thermal_r)

ylabels = [r'\textbf{U10}',r'\textbf{Z30}',r'\textbf{U300}',r'\textbf{Z500}',
           r'\textbf{SLP}',r'\textbf{T2M}',r'\textbf{RNET}']
plt.yticks(np.arange(0.5,7.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=0.7)

xlabels = [r'\textbf{DJF}']
plt.xticks(np.arange(0.5,4.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,1])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50)
ticks = np.arange(0,51,50)
labels = list(map(str,np.arange(0,51,50)))
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')
cbar.set_label(r'\textbf{Ratio [\%]}',
               color='dimgrey',labelpad=3,fontsize=12)

plt.subplots_adjust(top=0.8)

plt.savefig(directoryfigure + 'SITNET_ratio_mesh_DJFarea.png',dpi=300)