
import numpy as np

from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from climatology_utils import *
import cmocean

###########################################################

ncw = Dataset('../data/wnd_19940823_20120823.nc', 'r')

lat = ncw.variables['latitude'][:]
lat0 = lat[0]+(lat[-1]-lat[0])/2
lon = ncw.variables['longitude'][:]+360
lon0 = lon[0]+(lon[-1]-lon[0])/2
time = num2date(ncw.variables['time'][:], ncw.variables['time'].units)

u = ncw.variables['uwnd'][:]
v = ncw.variables['vwnd'][:]
w = np.ma.sqrt( u**2 + v**2)

monthly_wind =  group_by_month(date_time=time, data=w)
monthly_u = group_by_month(date_time=time, data=u)
monthly_v = group_by_month(date_time=time, data=v)

lat_case = [22,47]
lat0 = lat_case[0]+(lat_case[-1]-lat_case[0])/2
lon_case = [220, 255]
lon0 = lon_case[0]+(lon_case[-1]-lon_case[0])/2

# SWOT Calval tracks
data = np.loadtxt('../data/ephem_calval_june2015_sph.txt')
lon_track = data[:,1]
lat_track = data[:,2]

########################################################
sub=5
mon = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
xc, yc = np.meshgrid(lon, lat)
cmap = 'viridis'
map = Basemap(llcrnrlon=lon_case[0],llcrnrlat=lat_case[0],
        urcrnrlon=lon_case[-1],urcrnrlat=lat_case[-1],resolution='i',
        projection='lcc',lat_1=lat_case[0],lat_2=lat_case[-1],
        lat_0=lat0, lon_0=lon0)
xm, ym = map(xc,yc)
xs, ys = map(360-125.4, 35.4)
xt, yt = map(lon_track,lat_track)
parallels = np.round(np.arange(lat[0],lat[-1]-10,10))
meridians = np.round(np.linspace(220, 240,3))
levels = [7,9]
lin = np.linspace(4,10.1,60)
cs = map.contourf(xm,ym, monthly_wind['mean'][-1], lin, cmap=cmap)
fig=plt.figure(figsize=(14,11))
for i in np.arange(-1,11):
    ax = fig.add_subplot(3,4,i+2)
    if i in [5]:
        c = map.contour(xm,ym, monthly_wind['mean'][i], levels, colors='1')
        plt.clabel(c, c.levels, inline=True, manual=False, fmt = '%.1f', fontsize=20)
    pl = map.contourf(xm,ym, monthly_wind['mean'][i], lin, cmap=cmap)
    for cont in pl.collections:
        cont.set_edgecolor("face")
    
    mod = (monthly_u['mean'][i][::sub,::sub]**2 + monthly_v['mean'][i][::sub,::sub]**2)**0.5 
    map.quiver(xm[::sub,::sub],ym[::sub,::sub], monthly_u['mean'][i][::sub,::sub]/mod, monthly_v['mean'][i][::sub,::sub]/mod, scale=15, width=0.005, headwidth=4)
    map.plot(xt , yt, '--', color='.8', linewidth = 1.5)
    map.plot(xs,ys, 'Dr', markersize=10, mew=2, markeredgecolor='black')
    
    map.fillcontinents(color='0.5', lake_color='0.5')
    map.drawparallels(parallels,linewidth=.5,labels=[0,0,0,0], fontsize=22, color = '0.5')
    map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,0],fontsize=22, color = '0.6')
    if i in [-1,3,7]:
        map.drawparallels(parallels,linewidth=.5,labels=[1,0,0,0],fontsize=22, color = '0.5')
    if i in [7,8,9,10]:
        map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,1],fontsize=22, color = '0.6')
    plt.title(mon[i], fontsize=22)
plt.tight_layout()
# rect = [left, bottom, width, height]
cax = plt.axes([0.25, 0.0001, 0.5, 0.01])
cbar = plt.colorbar(cs, cax=cax, orientation = 'horizontal', format = '%.1f')
cbar.set_label('Wind Speed [m/s]',fontsize=22)
cbar.ax.tick_params(labelsize=22)
plt.subplots_adjust(hspace= 0.05, bottom=0.05)
plt.show()
