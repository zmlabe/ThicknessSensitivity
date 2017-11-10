"""
Functions are useful untilities for SITperturb experiments
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    [1] calcDecJan(varx,vary,lat,lon,level,levsq)
    [2] calcDecJanFeb(varx,vary,lat,lon,level,levsq)
    [3] calc_indttest(varx,vary)
    [4] calc_weightedAve(var,lats)
"""

def calcDecJan(varx,vary,lat,lon,level,levsq):
    """
    Function calculates average for December-January

    Parameters
    ----------
    varx : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    vary : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    level : string
        Height of variable (surface or profile)
    levsq : integer
        number of levels
        
    Returns
    -------
    varx_dj : 3d array or 4d array
        [year,lat,lon] or [year,lev,lat,lon]
    vary_dj : 3d array
        [year,lat,lon] or [year,lev,lat,lon]

    Usage
    -----
    varx_dj,vary_dj = calcDecJan(varx,vary,lat,lon,level,levsq)
    """
    print('\n>>> Using calcDecJan function!')
    
    ### Import modules
    import numpy as np
    
    ### Reshape for 3d variables
    if level == 'surface':    
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12),
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12),
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_dj = np.empty((varx.shape[0]-1,lat.shape[0],lon.shape[0]))
        vary_dj = np.empty((vary.shape[0]-1,lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djappendh = np.append(varxravel[11+i,:,:],varxravel[12+i,:,:])
            djappendf = np.append(varyravel[11+i,:,:],varyravel[12+i,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)
    ### Reshape for 4d variables
    elif level == 'profile':
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12.),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12.),levsq,
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_dj = np.empty((int(varx.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        vary_dj = np.empty((int(vary.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djappendh = np.append(varxravel[11+i,:,:,:],
                                  varxravel[12+i,:,:,:])
            djappendf = np.append(varyravel[11+i,:,:,:],
                                  varyravel[12+i,:,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                               
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
                                
    print('Completed: Organized data by months (ON,DJ,FM)!')

    print('*Completed: Finished calcDecJan function!')
    return varx_dj,vary_dj

def calcDecJanFeb(varx,vary,lat,lon,level,levsq):
    """
    Function calculates average for December-January-February

    Parameters
    ----------
    varx : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    vary : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    level : string
        Height of variable (surface or profile)
    levsq : integer
        number of levels
        
    Returns
    -------
    varx_djf : 3d array or 4d array
        [year,lat,lon] or [year,lev,lat,lon]
    vary_djf : 3d array
        [year,lat,lon] or [year,lev,lat,lon]

    Usage
    -----
    varx_djf,vary_djf = calcDecJanFeb(varx,vary,lat,lon,level,levsq)
    """
    print('\n>>> Using calcDecJan function!')
    
    ### Import modules
    import numpy as np
    
    ### Reshape for 3d variables
    if level == 'surface':    
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12),
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12),
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_djf = np.empty((varx.shape[0]-1,lat.shape[0],lon.shape[0]))
        vary_djf = np.empty((vary.shape[0]-1,lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djfappendh1 = np.append(varxravel[11+i,:,:],varxravel[12+i,:,:])
            djfappendf1 = np.append(varyravel[11+i,:,:],varyravel[12+i,:,:])  
            djfappendh = np.append(djfappendh1,varxravel[13+i,:,:])
            djfappendf = np.append(djfappendf1,varyravel[13+i,:,:]) 
            varx_djf[counter,:,:] = np.nanmean(np.reshape(djfappendh,
                                    (3,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)                   
            vary_djf[counter,:,:] = np.nanmean(np.reshape(djfappendf,
                                    (3,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)
    ### Reshape for 4d variables
    elif level == 'profile':
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12.),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12.),levsq,
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_djf = np.empty((int(varx.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        vary_djf = np.empty((int(vary.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djfappendh1 = np.append(varxravel[11+i,:,:,:],
                                  varxravel[12+i,:,:,:])
            djfappendf1 = np.append(varyravel[11+i,:,:,:],
                                  varyravel[12+i,:,:,:]) 
            djfappendh = np.append(djfappendh1,
                                  varxravel[13+i,:,:,:])
            djfappendf = np.append(djfappendf1,
                                  varyravel[13+i,:,:,:])  
            varx_djf[counter,:,:] = np.nanmean(np.reshape(djfappendh,
                                    (3,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                   
            vary_djf[counter,:,:] = np.nanmean(np.reshape(djfappendf,
                                    (3,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                               
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
                                
    print('Completed: Organized data by months (DJF)!')

    print('*Completed: Finished calcDecJanFeb function!')
    return varx_djf,vary_djf
    
def calc_indttest(varx,vary):
    """
    Function calculates statistical difference for 2 independent
    sample t-test

    Parameters
    ----------
    varx : 3d array
    vary : 3d array
    
    Returns
    -------
    stat = calculated t-statistic
    pvalue = two-tailed p-value

    Usage
    -----
    stat,pvalue = calc_ttest(varx,vary)
    """
    print('\n>>> Using calc_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_ind(varx,vary,nan_policy='omit')
    
    ### Significant at 95% confidence level
    pvalue[np.where(pvalue >= 0.05)] = np.nan
    pvalue[np.where(pvalue < 0.05)] = 1.
    
    print('*Completed: Finished calc_ttest function!')
    return stat,pvalue

def calc_weightedAve(var,lats):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    
    Parameters
    ----------
    var : 5d,4d,3d array of a gridded variable
    lats : 2d array of latitudes
    
    Returns
    -------
    meanvar : weighted average for 3d,2d,1d array

    Usage
    -----
    meanvar = calc_weightedAve(var,lats)
    """
    print('\n>>> Using calc_weightedAve function!')
    
    ### Import modules
    import numpy as np
    
    ### Calculate weighted average for various dimensional arrays
    if var.ndim == 5:
        meanvar = np.empty((var.shape[0],var.shape[1],var.shape[2]))
        for ens in range(var.shape[0]):
            for i in range(var.shape[1]):
                for j in range(var.shape[2]):
                    varq = var[ens,i,j,:,:]
                    mask = np.isfinite(varq) & np.isfinite(lats)
                    varmask = varq[mask]
                    areamask = np.cos(np.deg2rad(lats[mask]))
                    meanvar[ens,i,j] = np.nansum(varmask*areamask) \
                                        /np.sum(areamask)  
    elif var.ndim == 4:
        meanvar = np.empty((var.shape[0],var.shape[1]))
        for i in range(var.shape[0]):
            for j in range(var.shape[1]):
                varq = var[i,j,:,:]
                mask = np.isfinite(varq) & np.isfinite(lats)
                varmask = varq[mask]
                areamask = np.cos(np.deg2rad(lats[mask]))
                meanvar[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
    elif var.ndim == 3:
        meanvar = np.empty((var.shape[0],var.shape[1]))
        for i in range(var.shape[0]):
            varq = var[i,:,:]
            mask = np.isfinite(varq) & np.isfinite(lats)
            varmask = varq[mask]
            areamask = np.cos(np.deg2rad(lats[mask]))
            meanvar[i] = np.nansum(varmask*areamask)/np.sum(areamask)
    else:
        ValueError('Variable has the wrong dimensions!')
     
    print('Completed: Weighted variable average!')
    
    print('*Completed: Finished calc_weightedAve function!')
    return meanvar
