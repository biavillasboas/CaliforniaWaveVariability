import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
from climatology_utils import *
from scipy.ndimage.filters import gaussian_filter


# ETOPO file
tfile = '../data/etopo/ETOPO2v2c_f4.nc'

# Path to binned altimetry data
path = '../data/eastern_pacific_all_sat_binned_swh_1994-08-23_2012-08-23.nc' 

nc = Dataset(path, 'r')
time = num2date(nc.variables['time'][:], nc.variables['time'].units)
initial_time = time[0].date()
final_time = time[-1].date()
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

# Creating Mask for topography
etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z',\
                        sub = 1, lon_grid=lon, lat_grid=lat)

# mask shalower than 'depth'
depth = -100
swhcor = mask_etopo(bath=etopo_interp, var=nc.variables['swhcor'][:], depth=depth)

# computes monthly averages
monthly_swhcor = group_by_month(date_time=time, data=swhcor)

# Reads SWOT Calval ground-track
data = np.loadtxt('../data/ephem_calval_june2015_sph.txt')
lon_track = data[:,1]
lat_track = data[:,2]

# -------------- Define CASE domain limits ----------------
lat_case = [22,47]
lat0 = lat_case[0]+(lat_case[-1]-lat_case[0])/2
lon_case = [220, 255]
lon0 = lon_case[0]+(lon_case[-1]-lon_case[0])/2
# -------------------------------------------------

# Define map and plot variables
mon = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
xc, yc = np.meshgrid(lon, lat)
lin = np.linspace(1.1,4.6,60 )
cmap = 'YlGnBu_r'
map = Basemap(llcrnrlon=lon_case[0],llcrnrlat=lat_case[0],
        urcrnrlon=lon_case[-1],urcrnrlat=lat_case[-1],resolution='i',
        projection='lcc',lat_1=lat_case[0],lat_2=lat_case[-1],
        lat_0=lat0, lon_0=lon0)
xm, ym = map(xc,yc)
xs, ys = map(360-125.4, 35.4) # SWOT CalVal site
xt, yt = map(lon_track,lat_track) # SWOT ground-track
parallels = np.round(np.arange(lat[0],lat[-1]-10,10))
meridians = np.round(np.arange(lon[0]-10,lon[-1],15))
cs = map.contourf(xm,ym, monthly_swhcor['mean'][-1], lin, cmap=cmap)
levels = [2,2.5,3,3.5, 4]
#------------------------------------------------
# Plot 

fig=plt.figure(figsize=(14,11))
for i in np.arange(-1,11):
    ax = fig.add_subplot(3,4,i+2)
    data = gaussian_filter(monthly_swhcor['mean'][i],.8)
    c = map.contour(xm,ym, data, levels, colors='.1')
    plt.clabel(c, c.levels, inline=True, manual=True, fmt = '%.1f', fontsize=18)
    cf = map.contourf(xm,ym, monthly_swhcor['mean'][i], lin, cmap=cmap)
    for cont in cf.collections:
            cont.set_edgecolor("face")
    map.plot(xt , yt, '--', color='.8', linewidth = 1.5)
    map.plot(xs,ys, 'Dr', markersize=12, mew=2, markeredgecolor='black')
    map.drawparallels(parallels,linewidth=.5,labels=[0,0,0,0],fontsize=22, color = '0.5')
    map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,0],fontsize=22, color = '0.6')
    map.fillcontinents(color='0.5', lake_color='0.5')
    if i in [-1,3,7]:
        map.drawparallels(parallels,linewidth=.5,labels=[1,0,0,0],fontsize=22, color = '0.5')
    if i in [7,8,9,10]:
        map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,1],fontsize=22, color = '0.6')
    plt.title(mon[i], fontsize=22)
plt.tight_layout()
cax = plt.axes([0.3, 0.0001, 0.5, 0.01])
cbar = plt.colorbar(cs, cax=cax, orientation = 'horizontal', format = '%.1f')
cbar.set_label('Hs [m]',fontsize=22)
cbar.ax.tick_params(labelsize=22)
plt.subplots_adjust(hspace= 0.05, bottom=0.05)
plt.show()
