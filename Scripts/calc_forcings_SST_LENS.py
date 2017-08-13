"""
Forcing files for SITperturb from LENS (1979-2005; SST)

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 19 July 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as c
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm
import datetime

### Define directories
directorydata1 = '/surtsey/zlabe/LENS/SST/'
directorydata3 = '/surtsey/zlabe/LENS/ForcingPerturb/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate forcing files SST - %s----' % titletime 

ensembles = map(str,np.arange(2,10,1)) + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 2006
year2 = 2080
years = np.arange(year1,year2+1,1)

### Pick years
yearmin = 2051
yearmax = 2080
yearq = (np.where((years >= yearmin) & (years <= yearmax))[0])*12
          
### Read in functions
ssth = []
for i in xrange(len(ensembles)):
    filename = directorydata1 + 'SST_1920_2005_%03d.nc' % int(ensembles[i])
    data = Dataset(filename)
    lats = data.variables['ULAT'][:]
    lons = data.variables['ULONG'][:]
    ssthq = data.variables['SST'][yearq[0]:,0,:,:]
    ssth.append(ssthq)
    print('Completed: Read ensemble #%s!') % ensembles[i]
    data.close()
        
    del ssthq
    
sstf = []
for i in xrange(len(ensembles)):    
    if i > 31:
        filename = directorydata1 + 'SST_2006_2100_%03d.nc' % int(ensembles[i])
        data = Dataset(filename)
        lats = data.variables['ULAT'][:]
        lons = data.variables['ULONG'][:]
        sstfq = data.variables['SST'][yearq[0]:-20*12,0,:,:]
        sstf.append(sstfq)
        print('Completed: Read ensemble #%s!') % ensembles[i]
        data.close()
    else:
        filename = directorydata1 + 'SST_2006_2080_%03d.nc' % int(ensembles[i])
        data = Dataset(filename)
        lats = data.variables['ULAT'][:]
        lons = data.variables['ULONG'][:]
        sstfq = data.variables['SST'][yearq[0]:,0,:,:]
        sstf.append(sstfq)
        print('Completed: Read ensemble #%s!') % ensembles[i]
        data.close()
    
sstf = np.asarray(sstf)

### Average composite for years
sstave = np.mean(sstf,axis=0)
del sstf
sstfn = np.reshape(sstave,(30,12,lats.shape[0],lats.shape[1]))
sstfnn = np.mean(sstfn,axis=0)

print '\n Completed: Average over years %s - %s!' % (yearmin,yearmax)


def netcdfLENS(lats,lons,var,varqq,directory):
    print '\n>>> Using netcdf4LENS function!'
    
    name = 'lens_comps_%s_20062080.nc' % varqq
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'LENS %s interpolated on 1x1 grid' % varqq
    
    ### Dimensions
    ncfile.createDimension('months',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable(varqq,'f4',('months','lat','lon'))
    
    ### Units
    if varqq == 'sst':
        varns.units = 'K'
    elif varqq == 'sic':
        varns.units = 'fraction'
    elif varqq == 'sit':
        varns.units = 'm'
    ncfile.title = 'LENS %s' % varqq
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NCAR LENS'
    ncfile.references = 'Kay et al. [2013]'
    
    ### Data
    months[:] = np.arange(1,12+1,1)
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print '*Completed: Created netCDF4 File!'

#netcdfLENS(lats,lons,sstfnn,'sst',directorydata3)

data = Dataset(directorydata3 + 'lens_comps_sst_20062080.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sit = data.variables['sst'][:]
data.close()
data = Dataset(directorydata3 + 'lens_comps_sst_19762005.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sit2 = data.variables['sst'][:]
data.close()

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

#m = Basemap(projection='robin',lon_0=0,resolution='l')
m = Basemap(projection='ortho',lon_0=300,lat_0=-25,resolution='l')
         
var = sit[9] - sit2[9]
          
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=0.3)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,var,80,latlon=True,extend='both')
cs1 = m.contour(lons,lats,var,50,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('GMT_ocean')        
cs.set_cmap(cmap)

cbar = plt.colorbar(cs,extend='both')    
#cbar.set_label(r'\textbf{SIT}')  
#ticks = np.arange(0,8,1)
#cbar.set_ticks(ticks)
#cbar.set_ticklabels(map(str,ticks))     

plt.savefig(directoryfigure + 'testss_sstf.png',dpi=300)

print 'Completed: Script done!'