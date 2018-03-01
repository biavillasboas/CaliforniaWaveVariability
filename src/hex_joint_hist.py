import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset, num2date
from climatology_utils import *
from matplotlib.ticker import FormatStrFormatter

import seaborn as sns

stn = '029'
path_cdip = '../data/CDIP/%sp1_historic.nc'%stn
time_buoy, lat0, lon0, local_depth, buoyname, Hs, Dp, Tp = readCDIP(stn, path_cdip)

winter_months = [12,1,2]
summer_months = [5,6,7]

ind_summer = np.array([d.month in summer_months for d in time_buoy])
ind_winter = np.array([d.month in winter_months for d in time_buoy])

hs_summer = Hs[ind_summer]
tp_summer = Tp[ind_summer]
dp_summer = Dp[ind_summer]

hs_winter = Hs[ind_winter]
tp_winter = Tp[ind_winter]
dp_winter = Dp[ind_winter]

def joint_plots(data1, data2, name1, name2, figpath, color, letter):
    with sns.axes_style("ticks"):
        g = sns.jointplot(data1, data2, kind="hex", color = color, joint_kws={'gridsize':20}, marginal_kws={'bins':20,'color':'white'}, stat_func=None)
        g.ax_marg_x.hist(data1, 20, histtype='bar', color=color, edgecolor='k')
        g.ax_marg_y.hist(data2, 20, histtype='bar', orientation='horizontal',color=color, edgecolor='k')
        props = dict(boxstyle='round', facecolor='white', edgecolor='black', linewidth=1)
        g.fig.text(0.14, 0.17, letter, fontsize=26, verticalalignment='top', bbox=props)
        g.ax_joint.tick_params(labelsize=20)
        g.ax_joint.yaxis.set_major_formatter(FormatStrFormatter("%i"))
        g.ax_joint.xaxis.set_major_formatter(FormatStrFormatter("%i"))
        cax = g.fig.add_axes([.65, .2, .01, .3])
        cbar = plt.colorbar(cax=cax, format='%.1f')
        ma = cbar.vmax
        ticks = np.linspace(0,ma,7)
        labels = np.round(ticks/1e3,1)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=16)
        x0, x1 = g.ax_joint.get_xlim()
        g.set_axis_labels('%s' %name1, ' %s' %name2, fontsize=20,labelpad=14)
        g.savefig(figpath, dpi=300, bbox_inches='tight')
        plt.show()

# Summer plots
joint_plots(hs_summer.compressed(), tp_summer.compressed(), 'Hs [m]', 'Tp [s]', 'summer_hstp_hex.jpg', color='salmon', letter='e')
joint_plots(hs_summer.compressed(), dp_summer.compressed(), 'Hs [m]', 'Dp [degrees]', 'summer_hsdp_hex.jpg', color='steelblue', letter='d')
joint_plots(tp_summer.compressed(), dp_summer.compressed(), 'Tp [s]', 'Dp [degrees]', 'summer_tpdp_hex.jpg', color='grey', letter='f')


# Winter plots
joint_plots(hs_winter.compressed(), tp_winter.compressed(), 'Hs [m]', 'Tp [s]', '/Users/bia/Desktop/winter_hstp_hex.jpg', color='salmon', letter='b')
joint_plots(hs_winter.compressed(), dp_winter.compressed(), 'Hs [m]', 'Dp [degrees]', '/Users/bia/Desktop/winter_hsdp_hex.jpg', color='steelblue', letter='a')
joint_plots(tp_winter.compressed(), dp_winter.compressed(), 'Tp [s]', 'Dp [degrees]', '/Users/bia/Desktop/winter_tpdp_hex.jpg', color='grey', letter='c')
