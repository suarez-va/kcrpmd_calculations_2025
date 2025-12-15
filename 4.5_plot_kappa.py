"""
Plotting script for kappa transmission coefficient from Part 4

"""

import os
import matplotlib.pyplot as plt   # plots
import numpy as np

time = np.loadtxt("kappa_data/time.txt")
kappa_avg = np.loadtxt("kappa_data/kappa_avg.txt")
kappa_se = np.loadtxt("kappa_data/kappa_se.txt")
indices = (np.arange(1,6) * np.floor(len(time) / 6)).astype(int)

gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
fig_kw={'figsize':(5.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1)) = plt.subplots(1,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.set_xlabel("time (a.u.)", fontsize = 15)
ax1.set_ylabel(r"$\kappa$(t)", fontsize = 15)
ax1.plot(time, kappa_avg, color='k', linewidth=2)
#ax1.errorbar(time[indices], kappa_avg[indices], kappa_se[indices], fmt='none', ecolor='k', elinewidth=1, zorder=5)
ax1.errorbar(time[indices], kappa_avg[indices], np.absolute(kappa_se[indices]), fmt='none', ecolor='k', elinewidth=1, zorder=5)

plt.show()

plt.savefig('transmission_coefficient', bbox_inches='tight')

