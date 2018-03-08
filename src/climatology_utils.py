"""
"""
import numpy as np
import glob   
import os
from numpy import pi

from netCDF4 import Dataset, num2date
import datetime
from datetime import date, timedelta
from scipy import interpolate
from scipy.linalg import inv
from scipy.io import loadmat
import matplotlib.pyplot as plt


def closest_grid_to_buoy(lon0, lat0, lon_grid, lat_grid, data):

    ind_lat, ind_lon = \
        np.unravel_index(((lat_grid-lat0)**2 + (lon_grid-lon0)**2).argmin(), lat_grid.shape)
    if len(data.shape)==3:
        data_point = np.squeeze(data[:,ind_lat,ind_lon])
    if len(data.shape)==2:
        data_point = np.squeeze(data[ind_lat,ind_lon])
    return data_point

def closest_grid_point(lon0, lat0, lon_grid, lat_grid):

    ind_lat, ind_lon = \
        np.unravel_index(((lat_grid-lat0)**2 + (lon_grid-lon0)**2).argmin(), lat_grid.shape)

    return lon_grid[ind_lat, ind_lon], lat_grid[ind_lat,ind_lon]

def trig2az(theta):

    """
    Trigonometric angle to azimuth. Meteorological convention.

    """
    az0 = np.ma.masked_all(theta.shape)
    idx1 = (90>=theta) & (theta>=-90)
    idx2 = (90<theta) & (theta<=180)
    idx3 = (theta<-90) & (theta>=-180)
    az0[idx1] = abs(theta[idx1] - 90)
    az0[idx2] = (90 - theta[idx2]) + 360
    az0[idx3] = abs(theta[idx3]) + 90
    az = az0.copy()
    az[az0<=180] = az0[az0<=180] + 180
    az[az0>180] = az0[az0>180] - 180
    return az


def az2trig(az):
    assert ((az<=360) & (az>=0)).all(), "Azimuth out of range"
    theta = np.ma.masked_all(az.shape)
    ind1 = (az>=0.) & (az<=90.)
    ind2 = (az>90.) & (az<=360.)
    theta[ind1] = 90. - az[ind1]
    theta[ind2] =  90. - az[ind2] + 360.
    return theta


def direction_from_to(theta):
    direction = np.ma.masked_all(theta.shape)
    ind1 = (theta>=180.)
    ind2 = (theta<180.)
    direction[ind1] = theta[ind1] - 180
    direction[ind2] = theta[ind2] + 180
    return direction

def save_etopo(tfile, lon_topo, lat_topo, z_topo, sub, output):

    etopo = Dataset(tfile, 'r')
    lat = etopo.variables[lat_topo][::sub]
    lon = etopo.variables[lon_topo][::sub]
    z = etopo.variables[z_topo][::sub,::sub]
    if ((lon >= -180.).all() and (lon <= 180.).all()):
        lon0 = np.where(lon>=0)[0][0]
        lon = lon + 180
        topo = z.copy()
        topo[:,:lon0] = z[:,lon0:]
        topo[:,lon0:] = z[:,:lon0]
    else: 
        topo = z.copy()

    nc = Dataset(output, 'w', format='NETCDF4')
    lat_dim = nc.createDimension('latitude', lat.size)
    lon_dim = nc.createDimension('longitude', lon.size)

    vars={}
    vars['latitude'] = nc.createVariable('latitude',etopo.variables[lat_topo].dtype.str,('latitude',))
    vars['longitude'] = nc.createVariable('longitude',etopo.variables[lat_topo].dtype.str,('longitude',))
    vars['z'] = nc.createVariable('z', etopo.variables[z_topo].dtype.str, ('latitude','longitude',))

    attributes = etopo.variables['x'].ncattrs()
    try:
        attributes.remove('_FillValue')
    except:
        pass
    for a in attributes:
        setattr(vars['longitude'], a, getattr(etopo.variables['x'], a))

    attributes = etopo.variables['y'].ncattrs()
    try:
        attributes.remove('_FillValue')
    except:
        pass
    for a in attributes:
        setattr(vars['latitude'], a, getattr(etopo.variables['y'], a))

    attributes = etopo.variables['z'].ncattrs()
    try:
        attributes.remove('_FillValue')
    except:
        pass
    for a in attributes:
        setattr(vars['z'], a, getattr(etopo.variables['z'], a))


    vars['latitude'][:] = lat
    vars['longitude'][:] = lon
    vars['z'][:] = topo

    nc.close()

    return lat, lon, topo


def interp_etopo(tfile, lon_topo, lat_topo, z_topo, sub, lon_grid, lat_grid):

    etopo = Dataset(tfile, 'r')
    lats = etopo.variables[lat_topo][::sub]
    lons = etopo.variables[lon_topo][::sub]
    z = etopo.variables[z_topo][::sub,::sub]
    if ((lons >= -180.).all() and (lons <= 180.).all()):
        lon0 = np.where(lons>=0)[0][0]
        lons = lons + 180
        topo = z.copy()
        tmp = z[0,lon0:].shape[0]
        topo[:,:tmp] = z[:,lon0:]
        topo[:,tmp:] = z[:,:lon0]
    else: 
        topo = z.copy()

    grid = interpolate.interp2d(lons, lats, topo, kind='linear')
    lon, lat = lon_grid, lat_grid
    etopo_interp = grid(lon, lat)
    return etopo_interp

def mask_etopo(bath, var, depth):
    rs = bath.reshape(1, var.shape[1], var.shape[2])
    topo_3d = np.repeat(rs, var.shape[0],axis=0)
    mask_var = np.ma.masked_where(topo_3d>depth, var)
    return mask_var
    
def group_by_month(date_time, data):

    years = np.array([y.year for y in date_time])
    months = np.array([m.month for m in date_time])

    monthly_data = {}
    monthly_data = {'month': [], 'data':[],'mean': [], 'median': [],'std': [],'N': [] }
    for m in range(1,13):
        ind = months==m
        monthly_data['month'].append(m)
        tmp = np.ma.array(data[ind])
        monthly_data['data'].append(tmp)
        monthly_data['mean'].append(np.ma.mean(tmp, axis=0))
        monthly_data['median'].append(np.ma.median(tmp, axis=0))
        monthly_data['std'].append(np.ma.std(tmp, axis=0))
        monthly_data['N'].append(np.size(tmp.compressed()))

    return monthly_data

def monthly_average(date_time, data):

    assert isinstance(date_time, np.ndarray), 'date_time should be an array'

    years = np.array([y.year for y in date_time])
    months = np.array([m.month for m in date_time])

    monthly_data = {}
    monthly_data = {'time': [], 'mean': [], 'median': [],'std': [],'N': [] }
    for year in np.unique(years):
        for month in np.unique(months):
            ind_year = years == year 
            ind_month = months == month
            ind = ind_year*ind_month
            try:
                tmp = data[ind]
                time = date_time[ind]
                delta_t = datetime.timedelta(seconds=np.mean([(t-time[0]).total_seconds() for t in time]))
                monthly_data['time'].append(time[0] + delta_t)
                monthly_data['mean'].append(np.ma.mean(tmp, axis=0))
                monthly_data['median'].append(np.ma.median(tmp, axis=0))
                monthly_data['std'].append(np.ma.std(tmp, axis=0))
                monthly_data['N'].append(np.size(tmp.compressed()))
            except: pass

    return monthly_data

def time_average(date_time, data, time_step):

    t_size = np.float(data.shape[0])/time_step +1
    new_data = []
    new_time = []
    
    for n,i in enumerate(np.arange(0,data.shape[0]-time_step,time_step)):
        new_data.append(np.ma.mean(data[i:i+time_step], axis=0))
        time = date_time[i:i+time_step]
        delta_t = datetime.timedelta(seconds=np.mean([(t-time[0]).total_seconds() for t in time]))
        new_time.append(time[0] + delta_t)

    new_time = np.ma.array(new_time)
    new_data = np.ma.array(new_data)

    return new_data, new_time



def plot_annual_cycle(monthly_data, stderr, figtitle, figname):

    plt.figure(figsize=(12, 6))
    plt.plot(monthly_data['month'], monthly_data['mean'], 'b', label = 'Monthly mean')
    plt.fill_between(monthly_data['month'], monthly_data['mean'] - stderr,monthly_data['mean'] + stderr, facecolor='blue', alpha=0.3)
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    plt.xticks(np.arange(1,13), labels)
    plt.ylabel('Hs_corr [m]', fontsize =14, labelpad=14)
    plt.xlabel('Time [month]', fontsize =14, labelpad=14)
    plt.legend(fontsize =12)
    plt.title(figtitle)
    plt.savefig(figname, bbox_inches = 'tight', dpi=300)
    plt.show()

def plot_2cycles(monthly_data, stderr, figtitle, figname):
    
    t = np.arange(1,25)
    plt.figure(figsize=(12, 6))
    plt.plot(t, 2*monthly_data['mean'], 'b', label = 'Monthly mean')
    plt.fill_between(t, 2*monthly_data['mean'] - stderr,2*monthly_data['mean'] + stderr, facecolor='blue', alpha=0.3)
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    plt.xticks(np.arange(1,25), 2*labels)
    plt.ylabel('Hs_corr [m]', fontsize =14, labelpad=14)
    plt.xlabel('Time [month]', fontsize =14, labelpad=14)
    plt.legend(fontsize =12)
    plt.title(figtitle)
    plt.savefig(figname, bbox_inches = 'tight', dpi=300)
    plt.show()


def OverlapingDates(dates1, dates2):

    dates1 = np.asanyarray(dates1)
    dates2 = np.asanyarray(dates2)
    ovlp_d1 = []
    ind1 = []
    for dn, d in enumerate(dates1):
        if d in dates2:
            ovlp_d1.append(d)
            ind1.append(dn)

    ovlp_d2 = []
    ind2 = []
    for dn, d in enumerate(dates2):
        if d in dates1:
            ovlp_d2.append(d)
            ind2.append(dn)

    if (not ind1) & (not ind2):
        print "no overlaping dates"
    else:
        assert (dates1[ind1]==dates2[ind2]).all()

    return ind1, ind2


# CDIP Archived Dataset URL

def readCDIP(stn, path = 'url'):

    if path == 'url':
        data_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/' + stn + 'p1/' + stn + 'p1_historic.nc'
    else:
        data_url = path

    nc = Dataset(data_url)
    t = nc.variables['waveTime'][:]
    Hs = np.ma.array(nc.variables['waveHs'][:])
    Tp = np.ma.array(nc.variables['waveTp'][:])
    Dp = np.ma.array(nc.variables['waveDp'][:])

    time = num2date(t[:], nc.variables['waveTime'].units)
    lat = nc.variables['metaStationLatitude'][:]
    lon = nc.variables['metaStationLongitude'][:]
    depth = nc.variables['metaWaterDepth'][:]
    buoyname = " ".join(nc.variables['metaStationName'][:-45])

    return time, lat, lon, depth, buoyname, Hs, Dp, Tp

def regularly_spaced_daily(time, data):
    t0 = time[0].date()
    tf = time[-1].date()
    delta_days = (tf-t0).days + 1
    time_days = np.array( [ t0 + datetime.timedelta(days=d) for d in range(delta_days)] )
    Ndays = time_days.size
    date = np.array([t.date() for t in time])
    
    if len(data.shape)==3:
        Ny, Nx = data.shape[1:]
        hs_series = np.ma.masked_all([Ndays, Ny, Nx])
    
    elif len(data.shape)==2:
        Nx = data.shape[1]
        hs_series = np.ma.masked_all([Ndays,Nx])
    
    else:
        hs_series = np.ma.masked_all(Ndays)

    for i, d in enumerate(time_days):
        ind = date == d
        if ind.any():
            hs_series[i] = np.ma.mean(data[ind], axis=0)

    return time_days, hs_series

def interp_map(lon_in, lat_in, data, sub, lon_grid, lat_grid):

    etopo = Dataset(tfile, 'r')
    lats = etopo.variables[lat_in][::sub]
    lons = etopo.variables[lon_t][::sub]
    z = etopo.variables[data][::sub,::sub]
    if ((lons >= -180.).all() and (lons <= 180.).all()):
        lon0 = np.where(lons>=0)[0][0]
        lons = lons + 180
        new = data.copy()
        tmp = data[0,lon0:].shape[0]
        new[:,:tmp] = data[:,lon0:]
        new[:,tmp:] = data[:,:lon0]
    else:
        new = data.copy()

    grid = interpolate.interp2d(lons, lats, new, kind='linear')
    lon, lat = lon_grid, lat_grid
    data_interp = grid(lon, lat)
    return data_interp

def read_sat(path_sat, tfile):

    nc_sat = Dataset(path_sat, 'r')
    time_sat = num2date(nc_sat.variables['time'][:], nc_sat.variables['time'].units)
    time_sat = np.array([t.date() for t in time_sat])

    lat_sat = nc_sat.variables['lat'][:]
    lon_sat = nc_sat.variables['lon'][:]
    swh_sat = nc_sat.variables['swhcor'][:]

    etopo_interp = interp_etopo(tfile=tfile, lon_topo='longitude', lat_topo = 'latitude', z_topo = 'z',\
sub = 1, lon_grid=lon_sat, lat_grid=lat_sat)

    depth = -100
    xx_sat, yy_sat = np.meshgrid(lon_sat, lat_sat)
    ind = etopo_interp>=depth
    x_sat = np.ma.masked_array(xx_sat,mask=ind)
    y_sat = np.ma.array(yy_sat,mask=ind)

    return time_sat, x_sat, y_sat, swh_sat

def readCDIP(stn, path = 'url'):

    if path == 'url':
        data_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/' + stn + 'p1/' + stn + 'p1_historic.nc'
    else:
        data_url = path

    nc = Dataset(data_url)
    t = nc.variables['waveTime'][:]
    Hs = np.ma.array(nc.variables['waveHs'][:])
    Tp = np.ma.array(nc.variables['waveTp'][:])
    Dp = np.ma.array(nc.variables['waveDp'][:])

    time = num2date(t[:], nc.variables['waveTime'].units)
    lat = nc.variables['metaStationLatitude'][:]
    lon = nc.variables['metaStationLongitude'][:]
    depth = nc.variables['metaWaterDepth'][:]
    buoyname = " ".join(nc.variables['metaStationName'][:-45])

    return time, lat, lon, depth, buoyname, Hs, Dp, Tp
