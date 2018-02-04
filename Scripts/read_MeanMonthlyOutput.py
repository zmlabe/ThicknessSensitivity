"""
Script reads in MEAN monthly data from WACCM4 experiments (CIT,HIT,FIT)
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    readMeanExperi(directory,varid,experi,level)
"""

def readMeanExperi(directory,varid,experi,level):
    """
    Function reads monthly data from WACCM4 simulations

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
    var : 2d numpy array or 3d numpy array 
        [year,month] or [year,month,level]

    Usage
    -----
    lat,lon,time,lev,var = readMeanExperi(directory,varid,experi,level)
    """
    print('\n>>> Using readMeanExperi function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Call files
    totaldirectory = directory + experi + '/monthly/'
    filename = totaldirectory + varid + '_mean.nc'
    
    ### Read in Data
    if level == 'surface': # 1d variables
        data = Dataset(filename,'r')
        time = data.variables['time'][:]
        lev = 'surface'
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 2d variables
        data = Dataset(filename,'r')
        time = data.variables['time'][:]
        lev = data.variables['level'][:]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:4],varid))
    
    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 2d variables
        var = np.reshape(varq,((varq.shape[0]//12),months))
    elif level == 'profile': # 3d variables
        var = np.reshape(varq,((varq.shape[0]//12),months,int(lev.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')

    print('*Completed: Finished readExperi function!')
    return lat,lon,time,lev,var

#### Test function -- no need to use    
#directory = '/surtsey/zlabe/simu/'
#varid = 'LHFLX'
#experi = 'FIT'
#level = 'surface'
#lat,lon,time,lev,var = readMeanExperi(directory,varid,experi,level)
