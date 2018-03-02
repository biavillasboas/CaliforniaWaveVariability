import numpy as np

from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt

from climatology_utils import *
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter

ncw = Dataset('../data/wnd_19940823_20120823.nc', 'r')
ncf = Dataset('../data/fp_19940823_20120823.nc', 'r')

lat = ncw.variables['latitude'][:]
lat0 = lat[0]+(lat[-1]-lat[0])/2 
lon = ncw.variables['longitude'][:]+360
lon0 = lon[0]+(lon[-1]-lon[0])/2
time = num2date(ncw.variables['time'][:], ncw.variables['time'].units)

# ETOPO file
tfile = '/Users/bia/work/data/etopo/ETOPO2v2c_f4.nc'
# Creating Mask for topography
etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z',\
                                sub = 1, lon_grid=lon, lat_grid=lat)

fp = ncf.variables['fp'][:]
fp = np.ma.masked_where(fp<0, fp)

u = ncw.variables['uwnd'][:]
v = ncw.variables['vwnd'][:]
w = np.ma.sqrt( u**2 + v**2)

w_angle = np.arctan2(v,u)*180./np.pi
wdir = trig2az(w_angle)

g = 9.8
cp = g/(2*np.pi*fp)

A = cp/w

thrs = 1.2
N = np.zeros(A.shape)
N[A<=thrs] = 1

depth = -100
age = mask_etopo(bath=etopo_interp, var=A, depth=depth)
number = mask_etopo(bath=etopo_interp, var=N, depth=depth)

monthly_age =  group_by_month(date_time=time, data=age)
monthly_number =  group_by_month(date_time=time, data=number)

#################################################################
lat_case = [22,47]
lat0 = lat_case[0]+(lat_case[-1]-lat_case[0])/2
lon_case = [220, 255]
lon0 = lon_case[0]+(lon_case[-1]-lon_case[0])/2

mon = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
xc, yc = np.meshgrid(lon, lat)
map = Basemap(llcrnrlon=lon_case[0],llcrnrlat=lat_case[0],
        urcrnrlon=lon_case[-1],urcrnrlat=lat_case[-1],resolution='i',
        projection='lcc',lat_1=lat_case[0],lat_2=lat_case[-1],
        lat_0=lat0, lon_0=lon0)
xm, ym = map(xc,yc)
parallels = np.round(np.arange(lat[0],lat[-1]-10,10))
meridians = np.round(np.linspace(lon[0], lon[-1],3))

# SWOT Calval tracks
data = np.loadtxt('../data/ephem_calval_june2015_sph.txt')
lon_track = data[:,1]
lat_track = data[:,2]
xs, ys = map(360-125.4, 35.4)
xt, yt = map(lon_track,lat_track)
########################################################
Fraction of days of wind seas
########################################################
cmap = 'inferno'
lin = np.linspace(0,.5,10)
cs = map.contourf(xm,ym, monthly_number['mean'][-1], lin, cmap=cmap)
fig=plt.figure(figsize=(14,11), dpi=100)
for i in np.arange(-1,11):
    ax = fig.add_subplot(3,4,i+2)
    cnt = map.contourf(xm,ym, monthly_number['mean'][i], lin, cmap=cmap)
    map.plot(xt , yt, '--', color='.8', linewidth = 1.5)
    map.plot(xs,ys, 'Dr', markersize=12, mew=2, markeredgecolor='black')
    map.fillcontinents(color='0.5', lake_color='0.5')
    map.drawparallels(parallels,linewidth=.5,labels=[0,0,0,0],color = '0.5')
    map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,0],color = '0.6')
    if i in [-1,3,7]:
        map.drawparallels(parallels,linewidth=.5,labels=[1,0,0,0],fontsize=22, color = '0.5')
    if i in [7,8,9,10]:
        map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,1],fontsize=22, color = '0.6')
    plt.title(mon[i], fontsize=22)
plt.tight_layout()
plt.subplots_adjust(hspace= 0.05, bottom=0.05)
# rect = [left, bottom, width, height]
cax = plt.axes([0.25, 0.0001, 0.5, 0.01])
cbar = plt.colorbar(cs, cax=cax, orientation = 'horizontal', format = '%.1f')
cbar.set_label('fraction of wind-waves',fontsize=22)
cbar.set_ticks(np.linspace(0,.5,6), update_ticks=True)
cbar.set_ticklabels(np.linspace(0,.5,6), update_ticks=True)
cbar.ax.tick_params(labelsize=22)
plt.show()
