"""
Plot figure 2 in manuscript for dynamical responses to sea ice loss in WACCM4
experiments [FIT-HIT, FIC-CIT, FICT-HIT]. Current variables include T2M and
RNET. Time period includes December through February [DJF].

Notes
-----
    Author : Zachary Labe
    Date   : 4 February 2018
"""

### Import modules
import numpy as np
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import matplotlib.pyplot as plt

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
print('\n' '----Plotting Fig 2 - %s----' % titletime)

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
def readFlux(varnames):
    """
    Read in heat flux data for selected variables and calculate differences
    between experiments
    """
    for v in range(len(varnames)):
        ### Call function for surface temperature data from reach run
        lat,lon,time,lev,varhit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'HIT','surface')
        lat,lon,time,lev,varfit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'FIT','surface')
        lat,lon,time,lev,varcit = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'CIT','surface')
        lat,lon,time,lev,varfic = MO.readExperi(directorydata,
                                                '%s' % varnames[v],'FIC','surface')
        lat,lon,time,lev,varfict = MO.readExperi(directorydata,
                                                 '%s' % varnames[v],'FICT','surface')
        
    ### Compare experiments
    runs = [varhit,varfit,varcit,varfic,varfict]
    
    ### Compute comparisons for experiments - take ensemble average
    diff_FITHIT = np.nanmean(varfit - varhit,axis=0)*-1
    diff_FICCIT = np.nanmean(varfic - varcit,axis=0)*-1
    diff_FICTHIT = np.nanmean(varfict - varhit,axis=0)*-1
    diffruns = [diff_FITHIT,diff_FICCIT,diff_FICTHIT]
    
    return diffruns,runs,lat,lon

### Call function to read data for selected variable
diffruns_rnet,runs_rnet,lat,lon = readFlux(['RNET'])

difftotal_FITHITq = diffruns_rnet[0] + diffruns_rnet[0]
difftotal_FICCITq = diffruns_rnet[1] + diffruns_rnet[1]
difftotal_FICTHITq = diffruns_rnet[2] + diffruns_rnet[2]

difftotal_FITHIT = np.append(difftotal_FITHITq[8:],difftotal_FITHITq[:3],axis=0)
difftotal_FICCIT = np.append(difftotal_FICCITq[8:],difftotal_FICCITq[:3],axis=0)
difftotal_FICTHIT = np.append(difftotal_FICTHITq[8:],difftotal_FICTHITq[:3],axis=0)
difftotallhsh = [difftotal_FITHIT,difftotal_FICCIT,difftotal_FICTHIT]

### Take average above 40N
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
    
### Create files for rnet
np.savetxt(directorydata2 + 'weightedsic_rnets.txt',weightedrnet.transpose(),
           delimiter=',',header='  '.join(experiments)+'\n',
       footer='\n File contains net surface energy flux response' \
       '\n which are weighted above 40N for SIC cells >10% \n' \
       ' in all months of the year',newline='\n\n')

print('Completed: Script done!')

