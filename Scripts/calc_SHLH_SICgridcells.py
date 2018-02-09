"""
Calculate weighted turbulent heat flux (sensible + latent) for areas
over the control with a minimum of 10% sea ice concentration.

Notes
-----
    Author : Zachary Labe
    Date   : 9 February 2018
"""

### Import modules
import numpy as np
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculate weighted turbulent fluxes - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Define constants
runnames = [r'HIT',r'FIT',r'HIT2',r'FICT2',r'FICT']
experiments = ['FIT--HIT','FICT2--HIT2','FICT--HIT']

### Read in SIC data
lat,lon,time,lev,sic = MO.readExperi(directorydata,'SIC','HIT','surface')

### Find where ice is < 15% (values are 0 to 100 in sic array)
sicq = sic[5,:,:,:].copy()
sicq[np.where(sicq < 10)] = 0.0
sicq[np.where((sicq >= 10) & (sicq <= 100))] = 1.
sicn = np.append(sicq[8:],sicq[:3],axis=0)

###############################################################################
###############################################################################
###############################################################################
# Function to read surface heat flux data
def readFlux():
    """
    Read in heat flux data for selected variables and calculate differences
    between experiments
    """
    ### Call function for latent heat flux
    lat,lon,time,lev,lhhit = MO.readExperi(directorydata,
                                            'LHFLX','HIT','surface')
    lat,lon,time,lev,lhfit = MO.readExperi(directorydata,
                                            'LHFLX','FIT','surface')
    lat,lon,time,lev,lhcit = MO.readExperi(directorydata,
                                            'LHFLX','CIT','surface')
    lat,lon,time,lev,lhfic = MO.readExperi(directorydata,
                                            'LHFLX','FIC','surface')
    lat,lon,time,lev,lhfict = MO.readExperi(directorydata,
                                             'LHFLX','FICT','surface')
    ### Call function for sensible heat flux
    lat,lon,time,lev,shhit = MO.readExperi(directorydata,
                                            'SHFLX','HIT','surface')
    lat,lon,time,lev,shfit = MO.readExperi(directorydata,
                                            'SHFLX','FIT','surface')
    lat,lon,time,lev,shcit = MO.readExperi(directorydata,
                                            'SHFLX','CIT','surface')
    lat,lon,time,lev,shfic = MO.readExperi(directorydata,
                                            'SHFLX','FIC','surface')
    lat,lon,time,lev,shfict = MO.readExperi(directorydata,
                                            'SHFLX','FICT','surface')
    ### Calculate turbulent heat fluxes
    varhit = lhhit + shhit
    varfit = lhfit + shfit
    varcit = lhcit + shcit
    varfic = lhfic + shfic
    varfict = lhfict + shfict
        
    ### Compare experiments
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Compute comparisons for experiments - take ensemble average
    diff_FITHIT = np.nanmean(varfit - varhit,axis=0)
    diff_FICCIT = np.nanmean(varfic - varcit,axis=0)
    diff_FICTHIT = np.nanmean(varfict - varhit,axis=0)
    diffruns = [diff_FITHIT,diff_FICCIT,diff_FICTHIT]
    
    return diffruns,runs,lat,lon

### Call function to read data for selected variable
diffruns_rnet,runs_rnet,lat,lon = readFlux()

difftotal_FITHITq = diffruns_rnet[0] + diffruns_rnet[0]
difftotal_FICCITq = diffruns_rnet[1] + diffruns_rnet[1]
difftotal_FICTHITq = diffruns_rnet[2] + diffruns_rnet[2]

difftotal_FITHIT = np.append(difftotal_FITHITq[8:],difftotal_FITHITq[:3],axis=0)
difftotal_FICCIT = np.append(difftotal_FICCITq[8:],difftotal_FICCITq[:3],axis=0)
difftotal_FICTHIT = np.append(difftotal_FICTHITq[8:],difftotal_FICTHITq[:3],axis=0)
difftotallhsh = [difftotal_FITHIT,difftotal_FICCIT,difftotal_FICTHIT]

#### Take average above 40N
latq = np.where(lat > 40)[0]
latslice = lat[latq]
lon2,lat2 = np.meshgrid(lon,latslice)

### Mask out values not over SIC grid cells
rnetvals = []
for i in range(len(difftotallhsh)):
    rnetvalsq = difftotallhsh[i] * sicn
    rnetvalsq[np.where(rnetvalsq == 0.0)] = np.nan
    rnetvalsq = rnetvalsq[:,latq,:]
    
    rnetvals.append(rnetvalsq)
    
### Calculated weighted average 
weightedrnet = np.empty((len(rnetvals),sicn.shape[0]))
for i in range(len(rnetvals)):
    weightedrnet[i,:] = UT.calc_weightedAve(rnetvals[i],lat2)
    
#### Create files for rnet
np.savetxt(directorydata2 + 'weightedsic_SHLH.txt',weightedrnet.transpose(),
           delimiter=',',header='  '.join(experiments)+'\n',
       footer='\n File contains net turbulet energy flux response' \
       '\n which are weighted above 40N for SIC cells >10% \n' \
       ' in all months of the year',newline='\n\n')

print('Completed: Script done!')

fig = plt.figure()
ax = plt.subplot(111)
m = Basemap(projection='ortho',lon_0=300,lat_0=90,resolution='l')        
#var = sicn[4,latq,:]
a=rnetvals[2]
var = a[-1]
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=0.3)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
cs = m.contourf(lon2,lat2,var,55,latlon=True,extend='both')             
cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,extend='both')    
plt.savefig(directoryfigure + 'test_sic.png',dpi=300)

