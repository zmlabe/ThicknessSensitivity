"""
Forcing files for SITperturb from LENS (1979-2005; SST, SIC, SIT)
(2060-2080; SIT)

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 11 July 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as c
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm
import datetime
import read_var_LENS as LV
import read_SeaIceThick_LENS as lens

### Define directories
directorydata1 = '/home/zlabe/Surtsey3/CESM_large_ensemble/' 
directorydata2 = '/home/zlabe/Surtsey3/'
directorydata3 = '/surtsey/zlabe/LENS/ForcingPerturb/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate forcing files - %s----' % titletime 

ensembles = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 2006
year2 = 2080
years = np.arange(year1,year2+1,1)
          
### Read in functions
#sst,lats,lons = LV.readLENSEnsemble(directorydata1,'SST') # until 2080
#sic,lats,lons = LV.readLENSEnsemble(directorydata1,'SIC') # until 2080
#sit,lats,lons = lens.readLENSEnsemble(directorydata2,'None','historical')
#sit,lats,lons = lens.readLENSEnsemble(directorydata2,'None','rcp85')

### Pick years
yearmin = 2051
yearmax = 2080
yearq = np.where((years >= yearmin) & (years <= yearmax))[0]

### Average composite for years
#sstn = np.nanmean(sst[:,yearq,:,:,:],axis=1)
#sicn = np.nanmean(sic[:,yearq,:,:,:],axis=1)
#sitn = np.mean(sit[:,yearq,:,:,:],axis=1)

#sst = None
#sic = None
#sit = None

print '\n Completed: Average over years %s - %s!' % (yearmin,yearmax)

### Average over ensembles
#sst_ens = np.nanmean(sstn,axis=0)
#sic_ens = np.nanmean(sicn,axis=0)
#sit_ens = np.mean(sitn,axis=0)

#del sstn
#del sicn
#del sitn

#sit_ens[np.where(sit_ens > 12)] = np.nan 

print 'Completed: Average over all ensembles!'

def netcdfLENS(lats,lons,var,varqq,directory):
    print '\n>>> Using netcdf4LENS function!'
    
    name = 'lens_comp_%s_20512080.nc' % varqq
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'LENS %s interpolated on 1x1 grid' % varqq
    
    ### Dimensions
    ncfile.createDimension('months',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat'))
    longitude = ncfile.createVariable('lon','f4',('lon'))
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

#netcdfLENS(lats,lons,sst_ens,'sst',directorydata3)
#netcdfLENS(lats,lons,sic_ens,'sic',directorydata3)
#netcdfLENS(lats,lons,sit_ens,'sit',directorydata3)

data = Dataset(directorydata3 + 'lens_comp_sit_20512080.nc')
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
sit = data.variables['sit'][:]
data.close()

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='robin',lon_0=0,resolution='l')
m = Basemap(projection='npstere',boundinglat=67,lon_0=270,resolution='l',round =True)

sit[np.where(sit == 0)] = np.nan            
var = sit[9]

lons2,lats2 = np.meshgrid(lons,lats)
          
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=0.3)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,True,True],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons2,lats2,var,80,latlon=True,extend='both')
cs1 = m.contour(lons2,lats2,var,50,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)

def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(30,cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N-15)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
        
cmap = ncm.cmap('GMT_ocean')   
#cmap = colormapSIT()      
cs.set_cmap(cmap)

cbar = plt.colorbar(cs,extend='both')    
cbar.set_label(r'\textbf{SIT}')  
#ticks = np.arange(0,8,1)
#cbar.set_ticks(ticks)
#cbar.set_ticklabels(map(str,ticks))     

plt.savefig(directoryfigure + 'test_sit.png',dpi=300)

print 'Completed: Script done!'