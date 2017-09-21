"""
Script reads LENS data for selected variables
 
Notes
-----
    Author : Zachary Labe
    Date   : 28 November 2016
    
Usage
-----
    lats,lons,var = readLENS(directory,varq)
"""
    
def readLENSEnsemble(directory,varq):
    """
    Function reads LENS ensembles netCDF4 data array

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    varq : string
        variable from LENS

    Returns
    -------
    lats : 1d array
        latitudes
    lons : 1d array
        longitudes
    varq : 5d array [ens,year,month,lat,lon]
        selected variable

    Usage
    -----
    lats,lons,var = readLENS(directory,varq)
    """
    
    print('\n>>> Using readLENS function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ens = ['02','03','04','05','06','07','08','09'] + \
        list(map(str,np.arange(10,36,1))) + list(map(str,np.arange(101,106,1)))
    
    ### Modify directory
    directory = directory + '%s/' % (varq)
    
    if varq == 'SST':
        varn = np.empty((len(ens),75*12,384,320)) # 96 for all
        for i in range(len(ens)):
            if int(ens[i]) > 33:
                filename = '%s_2006_2100_0%s.nc' % (varq,ens[i])
                
                data = Dataset(directory + filename)
                lats = data.variables['ULAT'][:]
                lons = data.variables['ULONG'][:]
                varn[i,:,:,:] = np.squeeze(data.variables['%s' % varq][:-240,:,:]) # -2080
                data.close()
            else:
                filename = '%s_2006_2080_0%s.nc' % (varq,ens[i])
                
                data = Dataset(directory + filename)
                lats = data.variables['ULAT'][:]
                lons = data.variables['ULONG'][:]
                varn[i,:,:,:] = np.squeeze(data.variables['%s' % varq][:,:,:]) # -2080
                data.close()
            
            if int(ens[i]) > 100:
                filename = '%s_2006_2100_%s.nc' % (varq,ens[i])
                
                data = Dataset(directory + filename)
                lats = data.variables['ULAT'][:]
                lons = data.variables['ULONG'][:]
                varn[i,:,:,:] = np.squeeze(data.variables['%s' % varq][:-240,:,:]) # -2080
                data.close()
            
            print('Completed: Read LENS Ensemble #%s - %s!' % (ens[i],varq))
        
    else:
        varn = np.empty((len(ens),161*12,96,144)) # 96 for all
        for i in range(len(ens)):
            filename = '%s_0%s_1920_2080.nc' % (varq,ens[i])
            
            if int(ens[i]) > 100:
                filename = '%s_%s_1920_2100.nc' % (varq,ens[i])
                
            print(directory +filename)
                
            data = Dataset(directory + filename)
            lats = data.variables['latitude'][:]
            lons = data.variables['longitude'][:]
            varn[i,:,:,:] = data.variables['%s' % varq][:-240,:,:] # -2080
            data.close()
            
            print('Completed: Read LENS Ensemble #%s - %s!' % (ens[i],varq))
                
    var = np.reshape(varn,(len(ens),int(varn.shape[1]/12),12,
                           int(lats.shape[0]),int(lons.shape[0])))
    var = np.squeeze(np.asarray(var))
    
    ### Modify Units
    if varq == 'SLP':
        var = var/100. #Pa to hPa
    elif varq == 'T2M' or varq == 'T':
        var = var - 273.15 #K to C
    elif varq == 'SIT':
        var[np.where(var < 0)] = np.nan
        var[np.where(var > 12)] = np.nan
        
    ### Missing values    
    var[np.where(var <= -9999)] = np.nan
        
    print('*Completed: Read %s data!' % varq)
    
    return var,lats,lons
    
#var,lats,lons = readLENSEnsemble('/home/zlabe/Surtsey3/CESM_large_ensemble/','SST')