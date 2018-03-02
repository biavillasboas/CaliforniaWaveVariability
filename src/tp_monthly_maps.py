import numpy as np

from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from climatology_utils import *
import cmocean


nc = Dataset('../data/fp_19940823_20120823.nc', 'r')

lat = nc.variables['latitude'][:]
lat0 = lat[0]+(lat[-1]-lat[0])/2
lon = nc.variables['longitude'][:]+360
lon0 = lon[0]+(lon[-1]-lon[0])/2
time = num2date(nc.variables['time'][:], nc.variables['time'].units)

# ETOPO file
tfile = '../data/ETOPO2v2c_f4.nc'
# Creating Mask for topography
etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z',\
                                        sub = 1, lon_grid=lon, lat_grid=lat)

depth = -100
tp = 1./mask_etopo(bath=etopo_interp, var=nc.variables['fp'][:], depth=depth)
tp.mask[tp<0] = True

monthly_tp =  group_by_month(date_time=time, data=tp)

# --------- CASE -----------------------
lat_case = [22,47]
lat0 = lat_case[0]+(lat_case[-1]-lat_case[0])/2
lon_case = [220, 255]
lon0 = lon_case[0]+(lon_case[-1]-lon_case[0])/2


# SWOT Calval tracks
data = np.loadtxt('../data/ephem_calval_june2015_sph.txt')
lon_track = data[:,1]
lat_track = data[:,2]

mon = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
xc, yc = np.meshgrid(lon, lat)
map = Basemap(llcrnrlon=lon_case[0],llcrnrlat=lat_case[0],
        urcrnrlon=lon_case[-1],urcrnrlat=lat_case[-1],resolution='i',
        projection='lcc',lat_1=lat_case[0],lat_2=lat_case[-1],
        lat_0=lat0, lon_0=lon0)
xm, ym = map(xc,yc)
xs, ys = map(360-125.4, 35.4)
xt, yt = map(lon_track,lat_track)
parallels = np.round(np.arange(lat[0],lat[-1]-10,10))
meridians = np.round(np.linspace(220, 240,3))

#------------ contours----------------------------
stat = 'mean'
lin = np.linspace(8,15,12)
cmap = cmocean.cm.tempo

levels = [10,11,12]
cs = map.contourf(xm,ym, monthly_tp[stat][-1], lin, cmap=cmap)
fig=plt.figure(figsize=(14,11), dpi=100)
for i in np.arange(-1,11):
    ax = fig.add_subplot(3,4,i+2)
    pl = map.contourf(xm,ym, monthly_tp[stat][i], lin, cmap=cmap)
    #data = gaussian_filter(monthly_tp[stat][i],.8)
    #c = map.contour(xm,ym, data, levels, colors='.1')
    #plt.clabel(c, c.levels, inline=True, manual=False, fmt = '%.1f', fontsize=18)
    map.fillcontinents(color='0.5', lake_color='0.5')
    map.drawparallels(parallels,linewidth=.5,labels=[0,0,0,0],color = '0.5')
    map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,0],color = '0.6')
    map.plot(xt , yt, '--', color='.8', linewidth = 1.5)
    map.plot(xs,ys, 'Dr', markersize=10, mew=2, markeredgecolor='black')
    if i in [-1,3,7]:
        map.drawparallels(parallels,linewidth=.5,labels=[1,0,0,0],fontsize=22, color = '0.5')
    if i in [7,8,9,10]:
        map.drawmeridians(meridians,linewidth=.5,labels=[0,0,0,1],fontsize=22, color = '0.6')
    plt.title(mon[i], fontsize=22)
#plt.tight_layout(pad=6,h_pad=2,w_pad=6)
plt.tight_layout()
# rect = [left, bottom, width, height]
cax = plt.axes([0.25, 0.001, 0.5, 0.01])
cbar = plt.colorbar(cs, cax=cax, orientation = 'horizontal', format = '%i')
cbar.set_label('Tp [s]',fontsize=22)
cbar.ax.tick_params(labelsize=22)
cbar.set_ticks(range(8,16), update_ticks=True)
plt.subplots_adjust(hspace= 0.05, bottom=0.05)
plt.show()
