import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
mpl.rcParams['mathtext.fontset'] = 'cm'
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
cm = plt.get_cmap('Set1')
cm2 = plt.get_cmap('Set2')

mpl.rcParams.update({'errorbar.capsize': 2})

n = np.array([0.5, 1.5, 3.5, 7.])
t_min = np.array([0, 1, 2, 3])
t_max = np.array([1, 2, 3, 3])

fig,ax = plt.subplots(figsize=(4,2.8))

#plt.step(t, n, where='post')
plt.plot(t_min, n, marker = 'o', linestyle = '', color='k', markersize=4)
plt.hlines(n, t_min, t_max, linestyle = ':', linewidth=0.5)

plt.ylim(-1, 7.5)
plt.xlim(-0.1, 4)

ax.set_yticks(np.arange(-1, 8, 1))


ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_ylabel('Number of observed deficits', fontsize = 8)
ax.set_xlabel('Age', fontsize = 8)

ax.text(0.03, 0.13, r'$n_{t-1}$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, color='k',fontsize = 12, zorder=1000000)

ax.text(0.279, 0.25, r'$n_{t}$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, color='k',fontsize = 12, zorder=1000000)

ax.text(0.53, 0.485, r'$n_{t+1}$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, color='k',fontsize = 12, zorder=1000000)

ax.text(0.775, 0.89, r'$n_{t+2}$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, color='k',fontsize = 12, zorder=1000000)

plt.tight_layout()
plt.savefig('deficit_counts.pdf')



#n = np.array([0.5, 2, 4., 6., 8.])
n = np.array([0.5, 1.5, 3.5, 7.])
t_min = np.array([0, 1, 2, 3, 4])
t_max = np.array([1, 2, 3, 4, 5])

n_d = np.array([2., 2.5, 3.75])
n_r = np.array([1, 0.5, 0.25])

fig,ax = plt.subplots(figsize=(4,3.2))

#plt.step(t, n, where='post')
#plt.plot(t_min, n, marker = 'o', linestyle = '')
#plt.plot(t_max[:3], n_r, marker = 'o', linestyle = '', color = cm(2))
plt.plot(t_max[:3], n_r, marker = 'o', linestyle = '', color = '#009E73')
#plt.hlines(n, t_min, t_max, linestyle = ':')

plt.ylim(-1, 7)
plt.xlim(-0.1, 4)



ax.set_xticklabels([])
ax.set_yticklabels([])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_ylabel('Number repaired', fontsize = 18)
ax.set_xlabel('Age', fontsize = 18)

plt.tight_layout()
#plt.savefig('repair_counts.pdf')

fig,ax = plt.subplots(figsize=(4,3.2))

#plt.step(t, n, where='post')
#plt.plot(t_min, n, marker = 'o', linestyle = '')
#plt.plot(t_max[:3], n_d, marker = 'o', linestyle = '', color = cm(0))
plt.plot(t_max[:3], n_d, marker = 'o', linestyle = '', color = '#CC79A7')
#plt.plot(t_max[:2], n_r, marker = 's', linestyle = '')
#plt.plot(t_max[:2], n_d-n_r+np.array([0,2]), marker = '^', linestyle = '')
#plt.hlines(n, t_min, t_max, linestyle = ':')

plt.ylim(-1, 7)
plt.xlim(-0.1, 4)

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_ylabel('Number damaged', fontsize = 18)
ax.set_xlabel('Age', fontsize = 18)

plt.tight_layout()
#plt.savefig('n_d.pdf')



fig,ax = plt.subplots(2,1, figsize=(2,2.8))
ax=ax.flatten()

#plt.step(t, n, where='post')
#plt.plot(t_min, n, marker = 'o', linestyle = '')
ax[0].plot(t_max[:3], n_r, marker = 'o', linestyle = '', color = '#009E73', markersize=4)
#plt.hlines(n, t_min, t_max, linestyle = ':')

#plt.ylim(-1, 7)
#plt.xlim(-0.1, 4)

ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))

ax[0].set_xticklabels([])
ax[0].set_yticklabels([])

ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)

ax[0].set_ylabel('Number repaired', fontsize = 8)
#ax[0].set_xlabel('Age', fontsize = 18)

#plt.step(t, n, where='post')
#plt.plot(t_min, n, marker = 'o', linestyle = '')
ax[1].plot(t_max[:3], n_d, marker = 'o', linestyle = '', color = '#CC79A7', markersize=4)
#plt.plot(t_max[:2], n_r, marker = 's', linestyle = '')
#plt.plot(t_max[:2], n_d-n_r+np.array([0,2]), marker = '^', linestyle = '')
#plt.hlines(n, t_min, t_max, linestyle = ':')

#ax[1].set_ylim(-1, 7)
#ax[1].set_xlim(-0.1, 4)

ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))

ax[1].set_xticklabels([])
ax[1].set_yticklabels([])

ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)

ax[1].set_ylabel('Number damaged', fontsize = 8)
ax[1].set_xlabel('Age', fontsize = 8)

plt.tight_layout()
plt.savefig('repairdamage_counts.pdf')
