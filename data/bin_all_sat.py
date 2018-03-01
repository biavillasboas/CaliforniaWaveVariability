""" 
"""

import numpy as np
from datetime import datetime, timedelta
import os
import glob
import copy

from netCDF4 import Dataset, num2date, date2num

#           chose satellite from list: 
# sat = ['al', 'e1', ' e2',  'en',  'g2',  'h2',  'j1',  'j2',  'tp']

sat = 'all_sat'
datapath = '../data/ifremer_wave_merge_11.1/*'
years = sorted(glob.glob(datapath))

# Size of the regular grid in degrees
regular_grid = 1.
lon_size = int(360./regular_grid)
lat_size = int(180./regular_grid)

grid_lat = np.arange(-90,90,regular_grid)
grid_lon = np.arange(0,360,regular_grid)

# initializing the outputs
data_base = {}
# Here *_tmp stands for the variable and *_tmp2, to the variable squared
data_base['data_tmp'] = np.ma.zeros([lat_size, lon_size])
data_base['data_cor'] = np.ma.zeros([lat_size, lon_size])
data_base['N'] = np.ma.zeros([lat_size, lon_size])

for year in years:
    files = sorted(glob.glob(year + '/*.nc'))
    for f in files:
        nc0 = Dataset(f, 'r')
        swh = nc0.variables['swh'][:]
        swhcor = nc0.variables['swhcor'][:]
        lat = nc0.variables['lat'][:]
        lon = nc0.variables['lon'][:]
        lon[lon<0] = lon[lon<0]+360
        date_times =  num2date(nc0.variables['time'][:], nc0.variables['time'].units)
        days = np.array([t.date() for t in date_times])
        # One item at dict data for each day
        # Now process each measurement.
        day0 = np.unique(days)[0]
        dname = day0.isoformat()
        for i in range(len(days)):
            indlat = np.abs(grid_lat -lat[i]).argmin()
            indlon = np.abs(grid_lon -lon[i]).argmin()
            data_base['data_tmp'][indlat, indlon] += swh[i]
            data_base['data_cor'][indlat, indlon] += swhcor[i]
            data_base['N'][indlat, indlon] += 1

        print "Saving file for %s" %dname
        output = '../data/products/binned_swh_ifremer_1x1/'+sat+'/daily_files/binned_swh_ifremer_%s_%s.nc' % \
            (sat, np.unique(days)[0].strftime('%Y%m%d'))

        assert not os.path.isfile(output), \
            "Already exists: %s" % output

        nc = Dataset(output, 'w', format='NETCDF4')
        time_dim = nc.createDimension('time', 1)
        lat_dim = nc.createDimension('latitude', lat_size)
        lon_dim = nc.createDimension('longitude', lon_size)

        vars={}
        vars['time'] = nc.createVariable('time',
                nc0.variables['time'].dtype.str, ('time',))
        vars['lat'] = nc.createVariable('lat',
                nc0.variables['lat'].dtype.str, ('latitude',))
        vars['lon'] = nc.createVariable('lon',
                '<i4', ('longitude',))
        vars['swh'] = nc.createVariable('swh',
                nc0.variables['swh'].dtype.str, ('latitude','longitude',))
        vars['swhcor'] = nc.createVariable('swhcor',
                nc0.variables['swhcor'].dtype.str, ('latitude','longitude',))
        
        for my_variable in vars.keys():
            # Coping the attributes
            try:
                attributes = nc0.variables[my_variable].ncattrs()
                attributes.remove('_FillValue')
                for a in attributes:
                    setattr(vars[my_variable], a,
                            getattr(nc0.variables[my_variable], a))
            
            except: 
                attributes = nc0.variables[my_variable].ncattrs()
                for a in attributes:
                    setattr(vars[my_variable], a,
                            getattr(nc0.variables[my_variable], a))

        setattr(vars['lon'], 'valid_range', np.array([0,36000.]))
        
        vars['N'] = nc.createVariable('N',
                nc0.variables['swh'].dtype.str, ('latitude','longitude',))

        setattr(vars['N'], 'units', 'dimensionless' )
        setattr(vars['N'], 'long_name', 'Number of data points in each bin' )

        N = np.ma.masked_where(data_base['N']==0, data_base['N'])
        vars['lat'][:] = grid_lat
        vars['lon'][:] = grid_lon
        vars['N'][:] = N
        t0 = datetime(day0.year,day0.month,day0.day,0,0)
        vars['time'][:] = date2num(t0, vars['time'].units)
        vars['swh'][:] = np.ma.masked_where(data_base['N']==0, data_base['data_tmp'])/N
        vars['swhcor'][:] = np.ma.masked_where(data_base['N']==0, data_base['data_cor'])/N

        nc.close()

        # Remove what was already saved to keep memory under control
        print("Finished with %s. I'll remove it from the buffer" % dname)
        data_base = {}
        data_base['data_tmp'] = np.ma.zeros([lat_size, lon_size])
        data_base['data_cor'] = np.ma.zeros([lat_size, lon_size])
        data_base['N'] = np.ma.zeros([lat_size, lon_size])
