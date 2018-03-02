import sys
sys.path.append("../utils")

import numpy as np
from netCDF4 import Dataset, num2date
from climatology_utils import *
import matplotlib

# ETOPO file
tfile = '../data/ETOPO2v2c_f4.nc'

# Path to binned swh data
path = '../data/case_all_sat_binned_swh_1992-08-23_2016-08-23.nc'
nc = Dataset(path, 'r')

time = num2date(nc.variables['time'][:], nc.variables['time'].units)

data = nc.variables['swhcor'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
lon_grid, lat_grid = np.meshgrid(lon,lat)

lon0 = -125.4
lat0 = 35.4

hs0 = closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, data)

monthly_hs = group_by_month(time, hs0)

thrs = [2,2.5, 3,4,5]
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
percent_hs = {}
for m in months:
    percent_hs[m] = []

percent_hs = {}
for th in thrs:
    percent_hs[str(th)] = []

for th in thrs:
    for m in range(12):
        N = monthly_hs['data'][m].compressed().size
        percent_hs[str(th)].append(round(100.*sum(monthly_hs['data'][m].compressed() > th)/N,1))

month_number = np.arange(1,13)
mk = ['-o', '-^', '-d', '-s', '-*'] 
msizes = 7*np.ones(5,np.int)
msizes[-1] = 12

plt.figure()
for i, th in enumerate(thrs):
     perc = np.array(percent_hs[str(th)])
     plt.plot(month_number, perc, mk[i],markersize=msizes[i], color='k',linewidth=.5, label = str(th)+" m")
labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']
plt.xticks(np.arange(1,13), labels, fontsize=20)
plt.yticks(np.arange(0,100,20),fontsize=20)
plt.ylabel('Percent of Days [$\%$]', fontsize=20, labelpad = 14)
matplotlib.rcParams['legend.handlelength'] = 0
plt.grid(linestyle='dotted')
legend = plt.legend(loc = 'best', fancybox=True, frameon =True,edgecolor='k', fontsize=14, numpoints=1)
legend.get_frame().set_linewidth(1)
plt.show()
