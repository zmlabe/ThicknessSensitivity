"""
Functions are useful untilities for SITperturb experiments
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    [1] calcDecJan(varx,vary,lat,lon,level,levsq)
    [2] calc_indttest(varx,vary)
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
                           (varx.shape[0]*12.,
                            lat.shape[0],lon.shape[0]))
        varyravel = np.reshape(vary.copy(),
                               (vary.shape[0]*12.,
                                lat.shape[0],lon.shape[0])) 
                               
        varx_dj = np.empty((varx.shape[0]-1,lat.shape[0],lon.shape[0]))
        vary_dj = np.empty((vary.shape[0]-1,lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i/12
            djappendh = np.append(varxravel[11+i,:,:],varxravel[12+i,:,:])
            djappendf = np.append(varyravel[11+i,:,:],varyravel[12+i,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,lat.shape[0],lon.shape[0])),axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,lat.shape[0],lon.shape[0])),axis=0)
    ### Reshape for 4d variables
    elif level == 'profile':
        varxravel = np.reshape(varx.copy(),
                           (varx.shape[0]*12.,levsq,
                            lat.shape[0],lon.shape[0]))
        varyravel = np.reshape(vary.copy(),
                               (vary.shape[0]*12.,levsq,
                                lat.shape[0],lon.shape[0])) 
                               
        varx_dj = np.empty((varx.shape[0]-1,levsq,
                            lat.shape[0],lon.shape[0]))
        vary_dj = np.empty((vary.shape[0]-1,levsq,
                            lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i/12
            djappendh = np.append(varxravel[11+i,:,:,:],
                                  varxravel[12+i,:,:,:])
            djappendf = np.append(varyravel[11+i,:,:,:],
                                  varyravel[12+i,:,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,levsq,lat.shape[0],
                                     lon.shape[0])),axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,levsq,lat.shape[0],
                                     lon.shape[0])),axis=0)                               
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
                                
    print('Completed: Organized data by months (ON,DJ,FM)!')

    print('*Completed: Finished calcDecJan function!')
    return varx_dj,vary_dj
    
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
