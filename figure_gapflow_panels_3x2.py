'''
    Raul Valenzuela
    raul.valenzuela@colorado.edu

'''


import seaborn as sns
import gapflow
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import timedelta
from rv_utilities import discrete_cmap

from matplotlib import rcParams

rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 15
rcParams['mathtext.default'] = 'sf'

cmap = discrete_cmap(7, base_cmap='Set1')
color1 = cmap(0)
color2 = cmap(1)

# import os
# homed = os.path.expanduser('~')
homed = '/Users/raulvalenzuela/Data'


def get_subdates(gf, tta):

    bool_list = list()

    "mesowest obs have 55 minutes offset "
    offset = timedelta(minutes=55)

    for date in gf.index:
        if date-offset in tta:
            bool_list.append(True)
        else:
            bool_list.append(False)

    return bool_list


''' creates plot with seaborn style '''
with sns.axes_style("white"):
    sns.set_style('ticks',
                  {'xtick.direction': u'in',
                   'ytick.direction': u'in'}
                  )
    scale = 1.2
    fig = plt.figure(figsize=(6*scale, 8*scale))

    gs0 = gridspec.GridSpec(3, 2)

    axes = range(5)
    axes[0] = plt.subplot(gs0[0], gid='(a) 12-14Jan03')
    axes[1] = plt.subplot(gs0[1], gid='(b) 21-23Jan03')
    axes[2] = plt.subplot(gs0[2], gid='(c) 09Jan04')
    axes[3] = plt.subplot(gs0[3], gid='(d) 02Feb04')
    axes[4] = plt.subplot(gs0[4], gid='(e) 16-18Feb04')

# params = dict(wdir_surf=130, wdir_wprof=170,
#               rain_czd=0.25, nhours=2)

params = dict(wdir_thres=150,
              wdir_layer=[0, 500],
              rain_czd=0.25,
              nhours=2
              )
try:
    tta
except NameError:
    tta = gapflow.get_tta_dates([2003, 2004], params)

cases = [8, 9, 11, 12, 13]
wprof_hgt = 0.5  # [km] first gate to top layer
for n, c in enumerate(cases):
    gf = gapflow.run(c,
                     plot_theory=True,
                     grid=False,
                     wprof_hgt=wprof_hgt,
                     ax=axes[n],
                     homedir=homed,
                     color_surf=(0, 0, 0, 0),
                     color_wp=(0, 0, 0, 0)
                     )
    subdates = get_subdates(gf, tta.index)
    sub = gf.loc[subdates]
    x1, y1 = sub.pdiff.values, sub.ucomp.values
    x2, y2 = sub.pdiff.values, sub.wp_ucomp.values
    cl_surf = (0, 0, 0.5)
    cl_wprf = (0, 0.5, 0)
    h1 = axes[n].scatter(x1, y1, s=90, marker='x',
                         color=cl_surf, lw=2)
    h2 = axes[n].scatter(x2, y2, s=100, marker='+',
                         color=cl_wprf, lw=2)

    axes[n].text(0.05, 0.95, axes[n].get_gid(), size=14, va='top',
                 weight='bold', transform=axes[n].transAxes,
                 backgroundcolor='w', clip_on=True)

h3 = axes[0].lines[0]
h4 = axes[0].lines[1]
axes[2].legend([h1, h2, h3, h4],
               ['Surface',
                'Wprof (mean 0.1-{:3.1f} km MSL)'.format(wprof_hgt),
                'M95',
                r'$\pm$50% M95 $C_D$'],
               scatterpoints=1,
               loc=(0.01, 0.4))
axes[2].set_ylabel('BBY zonal wind $[m s^{-1}]$', ha='center')
axes[4].set_xlabel('Pressure difference, BBY-SCK [hPa]', ha='left')
axes[0].set_xticklabels([])
axes[1].set_xticklabels([])
axes[2].set_xticklabels([])
axes[1].set_yticklabels([])
axes[3].set_yticklabels([])
plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.95)

# plt.show()

place = '/Users/raulvalenzuela/Documents/'
fname = 'fig_gap_flow.png'
plt.savefig(place+fname,
            dpi=300, format='png',papertype='letter',
            bbox_inches='tight')
