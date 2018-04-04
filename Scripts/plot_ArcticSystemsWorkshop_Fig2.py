"""
Plot figure 3 in manuscript for dynamical responses to sea ice loss in WACCM4
experiments [FIT-HIT, FIC-CIT, FICT-HIT]. Current variables include SLP,
Z500, and Z30. Time period includes December through February [DJF].

Notes
-----
    Author : Zachary Labe
    Date   : 4 April 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import cmocean
import itertools

### Directory and time
directorydata = '/home/zlabe/Documents/Projects/Tests/SIV_animate/Data/'                   
directoryfigure = '/home/zlabe/Desktop/'              

now = datetime.datetime.now()
currentmn = str(now.month-1)
currentdy = str(now.day)
currentyr = str(now.year)
years = np.arange(1979,2018,1)

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Poster Figure 2 - %s----' % titletime)

### Read data
years,j,f,d = np.genfromtxt(directorydata + 'monthly_piomas.txt',
                           unpack=True,delimiter='',usecols=[0,1,2,12])

siv = (j[1:] + f[1:] + d[:-1])/3

### Plot Figure
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
ax = plt.subplot()

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

ax.yaxis.grid(zorder=1,color='darkgrey',alpha=1,linewidth=0.4)

plt.plot(years[1:],siv,color=cmocean.cm.balance(0.78),linewidth=3.5,marker='o',markersize=7,
         label=r'\textbf{PIOMAS v2.1 [Zhang and Rothrock, 2003]}')

plt.xticks(np.arange(1980,2021,10),list(map(str,np.arange(1980,2021,10))),
           fontsize=13,color='dimgrey')
plt.yticks(np.arange(14,29,2),list(map(str,np.arange(14,29,2))),fontsize=13,
           color='dimgrey')

plt.ylabel(r'\textbf{VOLUME [$\times$1000 km$^{3}$]}',
                     color='k',fontsize=16)
plt.title(r'\textbf{DEC-FEB : ARCTIC SEA ICE}',color='K',fontsize=27)

le = plt.legend(shadow=False,fontsize=8,loc='upper center',
           bbox_to_anchor=(0.27, 0.07),fancybox=True,frameon=False,ncol=1)
for text in le.get_texts():
    text.set_color('dimgrey') 

plt.xlim([1980,2020])
plt.ylim([14,28])

plt.savefig(directoryfigure + 'PosterFig2.png',dpi=1000)