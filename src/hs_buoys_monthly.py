
import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec

from climatology_utils import *

# ETOPO file
tfile = '../data/ETOPO2v2c_f4.nc'

# SWOT calval coordinates
lon_calval, lat_calval = 360-125.4, 35.4
##########################################################################
#  Load  WW3 significant wave height and direction at the calval site 
##########################################################################

# Loads Hs 
path_hs_ww3 = '../data/hs_19940823_20120823.nc'
nchs = Dataset(path_hs_ww3, 'r')
ww3_time = num2date(nchs.variables['time'][:], nchs.variables['time'].units)
xw, yw = np.meshgrid(nchs.variables['longitude'][:]+360, nchs.variables['latitude'][:])

lon0, lat0 = closest_grid_point(lon_calval, lat_calval, xw, yw)
hs_ww3 = nchs.variables['hs'][:]
ww3_hs0 = closest_grid_to_buoy(lon0, lat0, xw, yw, hs_ww3)

# Loads peak direction
path_dp_ww3 = '../data/dp_19940823_20120823.nc'
ncdp = Dataset(path_dp_ww3, 'r')
dp_ww3 = ncdp.variables['dp'][:]
dp0 = closest_grid_to_buoy(lon0, lat0, xw, yw, dp_ww3)

ww3_hs_series = monthly_average(ww3_time,ww3_hs0)
ww3_hs_calval = group_by_month(ww3_hs_series['time'],np.ma.array(ww3_hs_series['mean']))
ww3_std_calval = np.array(ww3_hs_calval['std'])
ww3_N_calval = np.array(ww3_hs_calval['N'])
ww3_stderr_calval = ww3_std_calval/np.sqrt(ww3_N_calval-1)

dp_series = monthly_average(ww3_time,dp0)
dp_calval = group_by_month(dp_series['time'],np.ma.array(dp_series['mean']))

dir_calval = np.pi*direction_from_to(az2trig(np.array(dp_calval['mean'])))/180

xc = np.cos(dir_calval)
yc = np.sin(dir_calval)

##############################################################################
# Load altimetry data at the Calval
##############################################################################
path_sat = '../data/eastern_pacific_all_sat_binned_swh_1994-08-23_2012-08-23.nc'
time, x, y, swh = read_sat(path_sat, tfile)
hs0 = closest_grid_to_buoy(lon0, lat0, x, y, swh)

monthly_calval_series = monthly_average(time,hs0)
monthly_calval = group_by_month(monthly_calval_series['time'],np.ma.array(monthly_calval_series['mean']))
std_calval = np.array(monthly_calval['std'])
N_calval = np.array(monthly_calval['N'])
stderr_calval =std_calval/np.sqrt(N_calval-1)

##########################################################################
# Load CDIP buoys lat and lon
#########################################################################
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

###########################################################################
# Plotting
##########################################################################
markers = ['-o','-s','-^','-*','-p', '-d','-v']
msizes = 7*np.ones(7,np.int)
msizes[3] = 12
colors = ['w','w','w','w','w','w', 'w', 'w', 'w', 'w']
mews = 2*np.ones(7, np.int)
leter = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
props = dict(boxstyle='round', facecolor='white')
#--------------------------------------------------------------------------
fig = plt.figure(figsize=(25,12))
gs0 = gridspec.GridSpec(2,1)
gs1 = gridspec.GridSpecFromSubplotSpec(2, 4, height_ratios=[4, 1], subplot_spec=gs0[0], hspace=0)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 4, height_ratios=[4, 1], subplot_spec=gs0[1],hspace=0)
ax = plt.subplot(gs1[0])
ax.plot(monthly_calval['month'],monthly_calval['mean'], '-D', color='r',\
        markersize=8, linewidth =.5, label = 'calval')
ax.fill_between(monthly_calval['month'],monthly_calval['mean'] - stderr_calval,\
        monthly_calval['mean'] + stderr_calval, facecolor='r', alpha=0.4, zorder=2)
ax.axhline(y=np.mean(monthly_calval['mean']), ls='--', color='red', linewidth=1)
ax.plot(ww3_hs_calval['month'],ww3_hs_calval['mean'], '-D', color='darkred',\
        markersize=4, linewidth =.5, label = 'WW3 calval')
ax.fill_between(ww3_hs_calval['month'],ww3_hs_calval['mean'] - \
        ww3_stderr_calval, ww3_hs_calval['mean'] + ww3_stderr_calval,\
        facecolor='darkred', alpha=0.4, zorder=2)
ax.axhline(y=np.mean(ww3_hs_calval['mean']), ls='--', color='darkred',\
        linewidth=1)
plt.xlim([1,15])
l = plt.legend(loc= 'upper right',fontsize =18, numpoints=1, frameon=True)
l.legendPatch.set_alpha(1)
matplotlib.rcParams['legend.handlelength'] = 0
plt.xticks(np.arange(1,13))
plt.grid(linestyle='dotted')
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.linspace(np.min(monthly_calval['mean']), np.max(\
    monthly_calval['mean']),5), fontsize=20)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.text(0.05, 0.14, leter[0], transform=ax.transAxes, fontsize=22,\
        verticalalignment='top', bbox=props)
ax1 = plt.subplot(gs1[4], sharex = ax)
ax1.quiver(np.array(monthly_calval['month']), 0, xc, yc,scale = 10)
plt.xlim([0,14])
labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']
plt.xticks(np.arange(1,13), labels, fontsize=20)
ax1.set_yticks([-0.05,0, 0.05])
ax1.grid(linestyle='dotted')
plt.setp(ax1.get_yticklabels(), visible=False)
labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']
plt.xticks(np.arange(1,13), labels, fontsize=20)
plt.grid(linestyle='dotted')
plt.setp(ax1.get_yticklabels(), visible=False)
numbers = [1,2,3,0,1,2,3]

for j, i in enumerate(ind):
    n = numbers[j]
    stn = buoys['id'][i]
    lon0 = buoys['longitude'][i]
    lat0 = buoys['latitude'][i]
    time, lat, lon, depth, buoyname, Hs, Dp, Tp = readCDIP(stn, path = 'url')
    monthly_date_series = monthly_average(date_time=time, data=Hs)
    monthly_data = group_by_month(monthly_date_series['time'],\
            np.ma.array(monthly_date_series['mean']))
    monthly_dir_series = monthly_average(date_time=time, data=Dp)
    monthly_dir = group_by_month(monthly_dir_series['time'],\
            np.ma.array(monthly_dir_series['mean']))

    dp= np.pi*direction_from_to(az2trig(np.array(monthly_dir['mean'])))/180.
    dx = np.cos(dp)
    dy = np.sin(dp)

    N = np.array(monthly_data['N'])
    std = monthly_data['std']
    stderr =std/np.sqrt(N-1)

    if j<=2:
        ax = plt.subplot(gs1[n])
        ax1 = plt.subplot(gs1[n+4], sharex = ax)
    else:
        ax = plt.subplot(gs2[n])
        ax1 = plt.subplot(gs2[n+4], sharex = ax)
    ax.plot(monthly_data['month'], monthly_data['mean'], '%s' %markers[i],\
            color='k', markerfacecolor= colors[i], markersize=msizes[i],\
            markeredgecolor='black',mew=mews[i], linewidth =.5, label='%s' %stn)
    ax.fill_between(monthly_data['month'], monthly_data['mean'] - stderr,\
            monthly_data['mean'] + stderr, facecolor='k', alpha=0.4, zorder=2)
    ax.axhline(y=np.mean(monthly_data['mean']), ls='--', color='k', linewidth=1)
    plt.xlim([1,15])
    l = ax.legend(loc= 'upper right',fontsize =18, numpoints=1, frameon=True)
    l.legendPatch.set_alpha(1)
    matplotlib.rcParams['legend.handlelength'] = 0
    plt.yticks(np.linspace(np.min(monthly_data['mean']), np.max(\
            monthly_data['mean']),5), fontsize=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xticks(np.arange(1,13), labels, fontsize=20)
    ax.grid(linestyle='dotted')
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.text(0.05, 0.14, leter[j+1], transform=ax.transAxes, fontsize=22,\
            verticalalignment='top', bbox=props)
    ax1.quiver(np.array(monthly_calval['month']), 0, dx, dy,scale = 10)
    plt.xlim([0,13])
    plt.xticks(np.arange(1,13), labels, fontsize=20)
    labels = ['','Feb','','Apr','','Jun','','Aug','','Oct','','Dec']
    ax1.set_yticks([-0.05,0, 0.05])
    ax1.grid(linestyle='dotted')
    plt.setp(ax1.get_yticklabels(), visible=False)
fig.text(0.07, 0.5, 'Hs [m]', va='center', rotation='vertical', fontsize=22)
matplotlib.rcParams['ytick.labelsize'] = 20
plt.show()
