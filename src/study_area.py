###########################################################
# Plots Figure 1 of Villas Boas et al. 2017
##########################################################

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import cm

import os, glob

# load etopo
tfile = '/Users/bia/work/data/etopo/ETOPO1_Ice_g_gmt4.grd'
nc_topo = Dataset(tfile, 'r')

lon = nc_topo.variables['x'][:] +180
lat = nc_topo.variables['y'][:]

ind_lat = (lat>=10) & (lat<=70)
ind_lon = (lon>=190) & (lon<=290)

lon = lon[ind_lon]
lat = lat[ind_lat]
xx, yy = np.meshgrid(lon, lat)

ind = np.where(nc_topo.variables['x'][:]<0)[0][-1]

# unwraps batymehtry for lon from 0 to 360
z_lat = nc_topo.variables['z'][ind_lat]
topo = np.ma.masked_all(z_lat.shape)
topo[:,:ind+1] = z_lat[:,ind:]
topo[:,ind+1:] = z_lat[:,:ind]

z = topo[:,ind_lon]
z.mask[z>0] = True

# CASE domain boundaries
lat_case = [22,47]
lat0 = lat_case[0]+(lat_case[-1]-lat_case[0])/2
lon_case = [220, 255]
lon0 = lon_case[0]+(lon_case[-1]-lon_case[0])/2

# Creates map
map = Basemap(llcrnrlon=lon_case[0],llcrnrlat=lat_case[0],
        urcrnrlon=lon_case[-1],urcrnrlat=lat_case[-1],resolution='i',
        projection='lcc',lat_1=lat_case[0],lat_2=lat_case[-1],
        lat_0=lat0, lon_0=lon0)

xm, ym = map(xx,yy)

# SWOT calval site coordinates
xs, ys = map(360-125.4, 35.4)

# SWOT Calval tracks
data = np.loadtxt('../data/ephem_calval_june2015_sph.txt')
lon_track = data[:,1]
lat_track = data[:,2]
xt, yt = map(lon_track,lat_track)

# Pt Conception and Pt Arena coordinates for text label
x_pa, y_pa = map(237, 38.9)
x_pc, y_pc= map(240.2, 34.8)

############################################################################3
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
markers = ['o','s','^','*', 'p', 'd', 'v']
msizes = 8*np.ones(7,np.int)
msizes[3] = 12
bbox_props = dict(boxstyle="round", facecolor='white')

lin = np.round(np.arange(-6300,100,100))
ticks = [-6000, -4000,-2000,-500,0]
parallels = np.round(np.arange(25,50,5))
meridians = np.arange(200, 255,5)

fig = plt.figure(figsize=(12, 10))
cs = map.contourf(xm,ym, z,lin, cmap=cm.Blues_r)
map.shadedrelief()
map.drawparallels(parallels,labels=[1,0,0,0],fontsize=18, color = '0.3')
map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=18, color = '0.3')
map.plot(xt , yt, '--', color='.8', linewidth = 1.5)
map.plot(xs,ys, 'Dr', markersize=12, mew=2, markeredgecolor='black', label = 'calval')
for i in ind:
    lat0 = buoys['latitude'][i]
    lon0 = buoys['longitude'][i]
    stn = buoys['id'][i]
    xb, yb = map(lon0,lat0)
    mk = markers[i]
    map.plot(xb, yb, mk ,color = 'w', markersize=msizes[i],
            markeredgecolor='black',mew=2, linestyle="None",label=stn)
plt.legend(loc= 'upper left', numpoints=1, fancybox=True, frameon =True,
        shadow=True, framealpha=.7, fontsize=18)
plt.text(x_pa, y_pa, 'Pt. Arena',fontsize=14, bbox=bbox_props)
plt.text(x_pc, y_pc, 'Pt. Conception',fontsize=14, bbox=bbox_props)
# rect = [left, bottom, width, height]
cax = plt.axes([0.25, 0.02, 0.5, 0.01])
cbar=plt.colorbar(cs, cax=cax,orientation='horizontal')
cbar.set_ticks(ticks)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Depth [m]',fontsize=16)
plt.show()
