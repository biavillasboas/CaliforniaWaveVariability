import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
from climatology_utils import *
import colorcet as cc


nc = Dataset('../data/dp_19940823_20120823.nc', 'r')

lat = nc.variables['latitude'][:]
lat0 = lat[0]+(lat[-1]-lat[0])/2
lon = nc.variables['longitude'][:]+360
lon0 = lon[0]+(lon[-1]-lon[0])/2
time = num2date(nc.variables['time'][:], nc.variables['time'].units)


# ETOPO file
tfile = '/Users/bia/work/data/etopo/ETOPO2v2c_f4.nc'
# Creating Mask for topography
etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z',\
                                        sub = 1, lon_grid=lon, lat_grid=lat)

depth = -100
dp = mask_etopo(bath=etopo_interp, var=nc.variables['dp'][:], depth=depth)

monthly_dp =  group_by_month(date_time=time, data=dp)

#------------ contours----------------------------
stats = 'mean'
directions = np.pi*direction_from_to(az2trig(dp))/180
xds = np.cos(directions)
yds = np.sin(directions)

xd = np.array(group_by_month(date_time=time, data=xds)['mean'])
yd = np.array(group_by_month(date_time=time, data=yds)['mean'])

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
#cmap = 'viridis'
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
stats = 'mean'

cmap = 'cet_cyclic_mrybm_35_75_c68'
lin = np.linspace(0,360,24)
sub = 3

cs = map.contourf(xm,ym, monthly_dp[stats][-1], lin, cmap=cmap)
fig=plt.figure(figsize=(14,11), dpi=100)
for i in np.arange(-1,11):
    ax = fig.add_subplot(3,4,i+2)
    pl = map.contourf(xm,ym, monthly_dp[stats][i], lin, cmap=cmap)
    mod = (xd[i][::sub,::sub]**2 + yd[i][::sub,::sub]**2)**0.5
    map.quiver(xm[::sub,::sub],ym[::sub,::sub], xd[i][::sub,::sub]/mod, yd[i][::sub,::sub]/mod, scale=15, width=0.005, headwidth=3)
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
plt.tight_layout()
cax = plt.axes([0.25, 0.001, 0.5, 0.01])
cbar = plt.colorbar(cs, cax=cax, orientation = 'horizontal', format = '%i')
cbar.set_label('Dp [degrees]',fontsize=22)
cbar.ax.tick_params(labelsize=22)
cbar.set_ticks([0,45,90,135,180,225,270,315,360], update_ticks=True)
plt.show()

