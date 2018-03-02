
import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec

from climatology_utils import *

ncw = Dataset('../data/wnd_19940823_20120823.nc', 'r')
time_wind = num2date(ncw.variables['time'][:], ncw.variables['time'].units)
lat = ncw.variables['latitude'][:]
lon = ncw.variables['longitude'][:]+360
lon_grid, lat_grid = np.meshgrid(lon, lat)
stats = 'mean'

# ETOPO file
tfile = '/Users/bia/work/data/etopo/ETOPO2v2c_f4.nc'
# Creating Mask for topography
etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z', sub = 1, lon_grid=lon, lat_grid=lat)
depth = -100

uwnd = mask_etopo(bath=etopo_interp, var=ncw.variables['uwnd'][:], depth=depth)
vwnd = mask_etopo(bath=etopo_interp, var=ncw.variables['vwnd'][:], depth=depth)

# Calval
lon0, lat0 = 360-125.4, 35.4
u0 = closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, uwnd)
v0 = closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, vwnd)
w0 = np.ma.sqrt( u0**2 + v0**2)

monthly_calval_series = monthly_average(time_wind,w0)
monthly_calval = group_by_month(monthly_calval_series['time'],np.ma.array(monthly_calval_series['mean']))
std_calval = np.array(monthly_calval['std'])
N_calval = np.array(monthly_calval['N'])
stderr_calval =std_calval/np.sqrt(N_calval-1)

um = np.array(group_by_month(time_wind, u0)[stats])
vm = np.array(group_by_month(time_wind, v0)[stats])
mod = (um**2 +vm**2)**.5
theta0 = np.arctan2(vm,um)*180/np.pi 
az = trig2az(theta0)

###############################################################

#CDIP
cdip_path = '../data/CDIP/*.nc'
cdip_files = sorted(glob.glob(cdip_path))
buoys={'latitude':[], 'longitude':[], 'id':[]}
for f in cdip_files:
    stn = os.path.basename(f)[:3]
    nc = Dataset(f,'r')
    buoys['latitude'].append(nc.variables['metaStationLatitude'][0])
    buoys['longitude'].append(nc.variables['metaStationLongitude'][0]+360)
    buoys['id'].append(stn)

ind = np.argsort(buoys['latitude'])[::-1]

markers = ['-o','-s','-^','-*','-p', '-d','-v']
msizes = 7*np.ones(7,np.int)
msizes[3] = 12
colors = ['w','w','w','w','w','w', 'w', 'w', 'w', 'w']
mews = 2*np.ones(7, np.int)
leter = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
props = dict(boxstyle='round', facecolor='white')
labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']

#############################################################

fig = plt.figure(figsize=(25,12))
gs0 = gridspec.GridSpec(2,1)
gs1 = gridspec.GridSpecFromSubplotSpec(2, 4, height_ratios=[4, 1], subplot_spec=gs0[0], hspace=0)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 4, height_ratios=[4, 1], subplot_spec=gs0[1],hspace=0)

ax = plt.subplot(gs1[0])
ax.plot(monthly_calval['month'],monthly_calval['mean'], '-D', color='r', markersize=8, linewidth =.5, label = 'calval')
ax.fill_between(monthly_calval['month'],monthly_calval['mean'] - stderr_calval,monthly_calval['mean'] + stderr_calval, facecolor='r', alpha=0.4, zorder=2)
ax.axhline(y=np.mean(monthly_calval['mean']), ls='--', color='k', linewidth=1)
plt.xlim([1,15])
l = plt.legend(loc= 'upper right',fontsize =18, numpoints=1, frameon=True)
l.legendPatch.set_alpha(1)
matplotlib.rcParams['legend.handlelength'] = 0
plt.xticks(np.arange(1,13))
plt.grid(linestyle='dotted')
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.linspace(np.min(monthly_calval['mean']), np.max(monthly_calval['mean']),5), fontsize=20)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.text(0.05, 0.14, leter[0], transform=ax.transAxes, fontsize=22, verticalalignment='top', bbox=props)
ax1 = plt.subplot(gs1[4], sharex = ax)
ax1.quiver(np.array(monthly_calval['month']), 0, um/mod, vm/mod,scale = 10)
plt.xlim([0,13])
plt.xticks(np.arange(1,13), labels, fontsize=20)
labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']
plt.grid(linestyle='dotted')
plt.setp(ax1.get_yticklabels(), visible=False)
numbers = [1,2,3,0,1,2,3]

for j, i in enumerate(ind):
    n = numbers[j]
    stn = buoys['id'][i]
    lon0 = buoys['longitude'][i]
    lat0 = buoys['latitude'][i]
    u = closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, uwnd)
    v = closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, vwnd)
    w = np.ma.sqrt( u**2 + v**2)
    um = np.array(group_by_month(time_wind, u)[stats])
    vm = np.array(group_by_month(time_wind, v)[stats])
    mod = (um**2 +vm**2)**.5
    monthly_wind_series = monthly_average(date_time=time_wind, data=w)
    monthly_wind =  group_by_month(date_time=monthly_wind_series['time'], data= np.ma.array(monthly_wind_series['mean']))
    N = np.array(monthly_wind['N'])
    std = monthly_wind['std']
    stderr =std/np.sqrt(N-1)
    if j<=2:
        ax = plt.subplot(gs1[n])
        ax1 = plt.subplot(gs1[n+4], sharex = ax)
    else:
        ax = plt.subplot(gs2[n])
        ax1 = plt.subplot(gs2[n+4], sharex = ax)
    ax.plot(monthly_wind['month'], monthly_wind['mean'], '%s' %markers[i], color='k', markerfacecolor= colors[i], markersize=msizes[i], markeredgecolor='black',mew=mews[i], linewidth =.5, label = '%s' %stn)
    ax.fill_between(monthly_wind['month'], monthly_wind['mean'] - stderr,monthly_wind['mean'] + stderr, facecolor='k', alpha=0.4, zorder=2)
    ax.axhline(y=np.mean(monthly_wind['mean']), ls='--', color='k', linewidth=1)
    plt.xlim([1,15])
    l = ax.legend(loc= 'upper right',fontsize =18, numpoints=1, frameon=True)
    l.legendPatch.set_alpha(1)
    matplotlib.rcParams['legend.handlelength'] = 0
    plt.yticks(np.linspace(np.min(monthly_wind['mean']), np.max(monthly_wind['mean']),5), fontsize=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xticks(np.arange(1,13), labels, fontsize=20)
    ax.grid(linestyle='dotted')
    plt.setp(ax.get_xticklabels(), visible=False)
    #ytick = np.linspace(np.min(monthly_wind['mean']), np.max(monthly_wind['mean']),5)
    #plt.yticks(ytick, fontsize=20)
    ax.text(0.05, 0.14, leter[j+1], transform=ax.transAxes, fontsize=22, verticalalignment='top', bbox=props)
    ax1.quiver(np.array(monthly_calval['month']), 0, um/mod, vm/mod,scale = 10)
    plt.xlim([0,13])
    plt.xticks(np.arange(1,13), labels, fontsize=20)
    ax1.set_yticks([-0.05,0, 0.05])
    ax1.grid(linestyle='dotted')
    plt.setp(ax1.get_yticklabels(), visible=False)
fig.text(0.07, 0.5, 'U$_{10}$ [m/s]', va='center', rotation='vertical', fontsize=22)
matplotlib.rcParams['ytick.labelsize'] = 20
plt.show()
