"""
Script reads in monthly latitude data from WACCM4 experiments 
(CIT,HIT,FIT,FICT,FIT)
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    readExperi(directory,varid,experi,level)
"""

def readExperi(directory,varid,experi,level):
    """
    Function reads monthly latitude data from WACCM4 simulations

    Parameters
    ----------
    directory : string
        working directory for stored WACCM4 experiments (remote server)
    varid : string
        variable name to read
    experi : string
        experiment name (CIT or HIT or FIT or FIC or FICT)
    level : string
        Height of variable (surface or profile)
        

    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    time : 1d numpy array
        standard time (days since 1870-1-1, 00:00:00)
    var : 3d numpy array or 4d numpy array 
        [year,month,lat] or [year,month,level,lat]

    Usage
    -----
    lat,lon,time,lev,var = readExperi(directory,varid,experi,level)
    """
    print('\n>>> Using readExperi function! \n')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Call files
    totaldirectory = directory + experi + '/monthly/'
    filename = totaldirectory + varid + '_1900-2000.nc'
    
    if any([experi == 'FPOL',experi == 'FSUB']):
        directory = '/home/zlabe/green/simu/'
        totaldirectory = directory + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
    
    ### Read in Data
    if level == 'surface': # 3d variables
        data = Dataset(filename,'r')
        varq = data.variables['%s' % varid][:,:,:,0]
        data.close()
        
        dataq = Dataset(totaldirectory + 'T2M_1900-2000.nc')
        time = dataq.variables['time'][:]
        lev = 'surface'
        lat = dataq.variables['latitude'][:]
        lon = dataq.variables['longitude'][:]
        dataq.close()
    elif level == 'profile': # 4d variables
        data = Dataset(filename,'r')
        varq = data.variables['%s' % varid][:,:,:,0]
        data.close()
        
        dataq = Dataset(totaldirectory + 'TEMP_1900-2000.nc')
        time = dataq.variables['time'][:]
        lev = dataq.variables['level'][:]
        lat = dataq.variables['latitude'][:]
        lon = dataq.variables['longitude'][:]
        dataq.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:4],varid))
    
    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(int(varq.shape[0]/12),months,
                              int(lat.shape[0])))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(int(varq.shape[0]/12),months,int(lev.shape[0]),
                      int(lat.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')

    print('\n*Completed: Finished readExperi function!')
    return lat,lon,time,lev,var

### Test function -- no need to use    
#directory = '/surtsey/zlabe/simu/'
#varid = 'EPZ'
#experi = 'HIT'
#level = 'profile'
#    
#lat,lon,time,lev,var = readExperi(directory,varid,experi,level)
