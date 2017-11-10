"""
Script calculates NAO index - currently a TEST script

Notes
-----
    Author : Zachary Labe
    Date   : 6 September 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput as DO
import calc_Utilities as UT
from eofs.standard import Eof
import scipy.stats as sts
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/nao/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating NAO Index - %s----' % titletime)

#### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call function for Z500 data (daily)
lat,lon,time,lev,z500_h = DO.readMeanExperi(directorydata,'Z500',
                                        'CIT','surface')
lat,lon,time,lev,z500_f = DO.readMeanExperi(directorydata,'Z500',
                                        'FICT','surface')
                                        
##### Calculate significance
#stat,pvalue = UT.calc_indttest(TEMP_h,TEMP_f)
                                                                    
### Calculate ensemble mean
z5_diffq = z500_f-z500_h

### Calculate for DJFM
z5_diff = z5_diffq[:,91:,:,:]     

#### Slice over (20-90N) and (90W-40E)
latq = np.where((lat>=20) & (lat<=90))[0]
latnao = lat[latq]

lonnew = np.mod(lon, 360.0) - 180.0
lonq = np.where((lonnew>=-90) & (lonnew<=40))[0]
lonnao = lonnew[lonq]

z5_diffn = z5_diff[:,:,latq,:]
z5_diffnao = z5_diffn[:,:,:,lonq]

z5n_h = np.nanmean(z500_h[:,91:,latq,:],axis=0)
z5nao_h = z5n_h[:,:,lonq]

### Calculate NAO
# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(latnao)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(z5nao_h, weights=wgts)

# Retrieve the leading EOF, expressed as the covariance between the leading PC
# time series and the input SLP anomalies at each grid point.
eof1 = solver.eofsAsCovariance(neofs=1).squeeze()
pc1 = solver.pcs(npcs=1, pcscaling=1).squeeze()

### Calculate NAO index
def NAOIndex(anomz5,eofpattern):
    """
    Calculate NAO index by regressing Z500 onto the EOF1 pattern
    """
    print('\n>>> Using NAO Index function!')       
    
    nao = np.empty((anomz5.shape[0],anomz5.shape[1]))
    for i in range(anomz5.shape[0]):
        print('Regressing ensemble ---> %s!' % (i+1))
        for j in range(anomz5.shape[1]):
            varx = np.ravel(anomz5[i,j,:,:])
            vary = np.ravel(eofpattern[:,:])
            mask = np.isfinite(varx) & np.isfinite(vary)     
            
            nao[i,j],intercept,r,p_value,std_err = sts.stats.linregress(
                                                                  varx[mask],
                                                                  vary[mask]) 
    print('*Completed: finished with NAO function!')
    return nao
            
### Calculate NAO index
naoindex = NAOIndex(z5_diffnao,eof1)
pc1 = (naoindex+np.mean(naoindex))/np.std(naoindex)
pc1 = np.nanmean(pc1,axis=0)
    
#### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)
    
varf = eof1[:,:]    

m = Basemap(projection='ortho',lon_0=-20,lat_0=60,resolution='l',
            area_thresh=10000.)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgray',linewidth=0.8)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0,color='k',fontsize=6)
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

# Make the plot continuous
barlim = np.arange(-8,9,4)
values = np.arange(-8,8.1,0.25)

lon2,lat2 = np.meshgrid(lonnao,latnao)

cs = m.contourf(lon2,lat2,varf,45,
                extend='both',latlon=True)
cs1 = m.contour(lon2,lat2,varf,
                linewidths=0.1,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap) 

cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)                    
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(labelsize=8)   

plt.savefig(directoryfigure + 'testeof1_FICT.png',dpi=300)

############################################################################
############################################################################
############################################################################
#### Plot NAO index
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')

ax.yaxis.grid(zorder=1,color='darkgrey',alpha=0.35)

zeroline = [0]*122

pc1_masked = np.ma.masked_less_equal(pc1, 0)

plt.bar(np.arange(len(pc1)),pc1,color='tomato',edgecolor='tomato',zorder=9) 
plt.bar(np.arange(len(pc1)),pc1_masked.filled(np.nan),color='dodgerblue',
        edgecolor='dodgerblue',zorder=10)

#plt.plot(pc1,linewidth=2.5,color='dodgerblue',alpha=1,
#        linestyle='-')

#plt.legend(shadow=False,fontsize=9,loc='lower center',
#           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{NAO Index (Z500)}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,1),list(map(str,np.arange(-5,6,1))),fontsize=9)
plt.ylim([-3,3])

xlabels = [r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,121,30),xlabels,fontsize=9)
plt.xlim([0,121])

plt.savefig(directoryfigure + 'NAO_Index_FIT-HIT.png',dpi=300)