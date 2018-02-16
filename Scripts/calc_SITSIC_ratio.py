"""
Compute ratio (%) between SIT and SIC responses

Notes
-----
    Author : Zachary Labe
    Date   : 15 February 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput as MO
import cmocean
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
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
print('\n' '----Plotting SIT-SIC ratio - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

months = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
varnames = ['U10','Z30','U300','Z500','SLP','T2M','RNET']
#varnames = ['SLP']

ratiovar = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'CIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FIT',r'FIC',r'CIT']
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIC--CIT}']
    runs = [varhit,varfit,varfic,varcit]
    
    ### Separate per 2 month periods
    varmo_on = np.empty((4,varhit.shape[0],varhit.shape[2],varhit.shape[3]))
    varmo_dj = np.empty((4,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
    varmo_fm = np.empty((4,varhit.shape[0],varhit.shape[2],varhit.shape[3]))
    for i in range(len(runs)):
        varmo_on[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1)    
        varmo_dj[i],varmo_dj[i] = UT.calcDecJan(runs[i],runs[i],lat,lon,'surface',1)    
        varmo_fm[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    diff_fithit_on = np.nanmean(varmo_on[1] - varmo_on[0],axis=0)
    diff_ficcit_on = np.nanmean(varmo_on[2] - varmo_on[3],axis=0)
    
    diff_fithit_dj = np.nanmean(varmo_dj[1] - varmo_dj[0],axis=0)
    diff_ficcit_dj = np.nanmean(varmo_dj[2] - varmo_dj[3],axis=0)
    
    diff_fithit_fm = np.nanmean(varmo_fm[1] - varmo_fm[0],axis=0)
    diff_ficcit_fm = np.nanmean(varmo_fm[2] - varmo_fm[3],axis=0)
    
    ### Calculate significance 
    stat_FITHITon,pvalue_FITHITon = UT.calc_indttest(varmo_on[1],varmo_on[0])
    stat_FICCITon,pvalue_FICCITon = UT.calc_indttest(varmo_on[2],varmo_on[3])

    stat_FITHITdj,pvalue_FITHITdj = UT.calc_indttest(varmo_dj[1],varmo_dj[0])
    stat_FICCITdj,pvalue_FICCITdj = UT.calc_indttest(varmo_dj[2],varmo_dj[3])

    stat_FITHITfm,pvalue_FITHITfm = UT.calc_indttest(varmo_fm[1],varmo_fm[0])
    stat_FICCITfm,pvalue_FICCITfm = UT.calc_indttest(varmo_fm[2],varmo_fm[3])
    
    ### Create mask of significant values
    pvalue_FITHITon[np.where(np.isnan(pvalue_FITHITon))] = 0.0
    pvalue_FICCITon[np.where(np.isnan(pvalue_FICCITon))] = 0.0

    pvalue_FITHITdj[np.where(np.isnan(pvalue_FITHITdj))] = 0.0
    pvalue_FICCITdj[np.where(np.isnan(pvalue_FICCITdj))] = 0.0

    pvalue_FITHITfm[np.where(np.isnan(pvalue_FITHITfm))] = 0.0
    pvalue_FICCITfm[np.where(np.isnan(pvalue_FICCITfm))] = 0.0
        
    pvalue_FITHIT = [pvalue_FITHITon,pvalue_FITHITdj,pvalue_FITHITfm]
    pvalue_FICCIT = [pvalue_FICCITon,pvalue_FICCITdj,pvalue_FICCITfm]
    
    ### Create mask of shared significant values
    mask = np.asarray(pvalue_FITHIT) * np.asarray(pvalue_FICCIT)
    
    ### Slice out lats below 40
    latq = np.where(lat>40)[0]
    latqq = lat[latq]
    
    ### Create 2nd meshgrid with lats > 40N
    lonnew,latnew=np.meshgrid(lon,latqq)
    
    ### Create mask for ON, DJ, FM
    mask = mask[:,latq,:]
    
    ### Keep only values significant in both SIT and SIC responses
#    diff_fithit_onq = diff_fithit_on[latq,:] * mask[0,:,:]
#    diff_fithit_djq = diff_fithit_dj[latq,:] * mask[1,:,:]
#    diff_fithit_fmq = diff_fithit_fm[latq,:] * mask[2,:,:]
#    
#    diff_ficcit_onq = diff_ficcit_on[latq,:] * mask[0,:,:]
#    diff_ficcit_djq = diff_ficcit_dj[latq,:] * mask[1,:,:]
#    diff_ficcit_fmq = diff_ficcit_fm[latq,:] * mask[2,:,:]
    
    diff_fithit_onq = diff_fithit_on[latq,:] * pvalue_FITHITon[latq,:]
    diff_fithit_djq = diff_fithit_dj[latq,:] * pvalue_FITHITdj[latq,:]
    diff_fithit_fmq = diff_fithit_fm[latq,:] * pvalue_FITHITfm[latq,:]
    
    diff_ficcit_onq = diff_ficcit_on[latq,:] * pvalue_FICCITon[latq,:]
    diff_ficcit_djq = diff_ficcit_dj[latq,:] * pvalue_FICCITdj[latq,:]
    diff_ficcit_fmq = diff_ficcit_fm[latq,:] * pvalue_FICCITfm[latq,:]
    
    ### Change 0 to nan as to no affect the averaging
#    diff_fithit_onq[np.where(diff_fithit_onq == 0.0)] = np.nan
#    diff_fithit_djq[np.where(diff_fithit_djq == 0.0)] = np.nan
#    diff_fithit_fmq[np.where(diff_fithit_fmq == 0.0)] = np.nan
#    
#    diff_ficcit_onq[np.where(diff_ficcit_onq == 0.0)] = np.nan
#    diff_ficcit_djq[np.where(diff_ficcit_djq == 0.0)] = np.nan
#    diff_ficcit_fmq[np.where(diff_ficcit_fmq == 0.0)] = np.nan
    
    fithit = [diff_fithit_onq,diff_fithit_djq,diff_fithit_fmq]
    ficcit = [diff_ficcit_onq,diff_ficcit_djq,diff_ficcit_fmq]
    
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
#        fithit[i][np.where(fithit[0] == 0.0)] = np.nan
#        ficcit[i][np.where(ficcit[0] == 0.0)] = np.nan
        fithitave[i] = UT.calc_weightedAve(abs(fithit[i]),latnew)
        ficcitave[i] = UT.calc_weightedAve(abs(ficcit[i]),latnew)

    ratio = []
    for i in range(len(fithit)):
        percchangeq,varx,vary = calc_iceRatio(fithitave[i],ficcitave[i],False,95,5)
        
        ratio.append(percchangeq)
    ratiovar.append(ratio)
meanratiovar = np.asarray(ratiovar).squeeze()
#ratiovar[np.where(np.isnan(ratiovar))] = 0.0
#meanratiovar = UT.calc_weightedAve(ratiovar[:,:,:,:],latnew)

varyy = abs(fithit[0])
fig = plt.figure()
ax = plt.subplot(111)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=0.3)
parallels = np.arange(-90,90,30)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
cs = m.contourf(lonnew,latnew,varyy[:,:],55,latlon=True,extend='both')             
cs.set_cmap(cmocean.cm.thermal)
cbar = plt.colorbar(cs,extend='both')    
plt.savefig(directoryfigure + 'test_ratio.png',dpi=300)
    
#### Save file
np.savetxt(directorydata2 + 'sicsitratio.txt',meanratiovar.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains ratio values of relative contributions' \
           '\n between FIT-HIT and FIC-CIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]',newline='\n\n')

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

cs = plt.pcolormesh(meanratiovar,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=0,vmax=50)

for i in range(meanratiovar.shape[0]):
    for j in range(meanratiovar.shape[1]):
        plt.text(j+0.5,i+0.5,r'\textbf{%3.1f}' % meanratiovar[i,j],fontsize=6,
                 color='r',va='center',ha='center')

cs.set_cmap(cmocean.cm.tempo)

ylabels = [r'\textbf{U10}',r'\textbf{Z30}',r'\textbf{U300}',r'\textbf{Z500}',
           r'\textbf{SLP}',r'\textbf{T2M}',r'\textbf{RNET}']
plt.yticks(np.arange(0.5,7.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=0.7)

xlabels = [r'\textbf{ON}',r'\textbf{DJ}',r'\textbf{FM}']
plt.xticks(np.arange(0.5,4.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,3])

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

plt.savefig(directoryfigure + 'SITSIC_ratio_mesh.png',dpi=300)