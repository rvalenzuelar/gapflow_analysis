
import seaborn as sns
import gapflow
import matplotlib.pyplot as plt
#import Windprof2 as wp
import xpol_tta_analysis as xta
import numpy as np
import matplotlib.gridspec as gridspec
from datetime import timedelta

from matplotlib import rcParams
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 15

#import os
#homed = os.path.expanduser('~')
homed = '/localdata'


def get_subdates(gf,tta):
    
    subdates=[]
    for date in gf.index:
        if date in tta:
            subdates.append(True)
        else:
            subdates.append(False)

    return subdates

''' creates plot with seaborn style '''
with sns.axes_style("white"):
    plt.figure(figsize=(8, 8))
    
    gs0 = gridspec.GridSpec(2, 2)

    axes = range(4)
    axes[0] = plt.subplot(gs0[0],gid='(a) 12-14Jan03')
    axes[1] = plt.subplot(gs0[1],gid='(b) 21-23Jan03')
    axes[2] = plt.subplot(gs0[2],gid='(c) 02Feb03')
    axes[3] = plt.subplot(gs0[3],gid='(d) 16-18Feb04')

params = dict(wdir_surf=130,wdir_wprof=170,
              rain_czd=0.25,nhours=2)

try:
    tta
except NameError:
    tta = np.array(xta.get_tta_dates([2003,2004],params))
    tta = tta + timedelta(minutes=55)
    

cases=[8,9,12,11]
for n,c in enumerate(cases):
    gf = gapflow.run(c, plot=True, grid=False,
                     ax=axes[n], homedir=homed,
                     color_surf=(0, 0, 0, 0),
                     color_wp=(0, 0, 0, 0))
    
#    tta = wp.get_tta_times(case=str(n+8), homedir=homed)
#    sub = gf.loc[tta + timedelta(minutes=55)]
    
    
    subdates = get_subdates(gf,tta)
    sub = gf.loc[subdates]
    x1, y1 = sub.pdiff.values, sub.ucomp.values
    x2, y2 = sub.pdiff.values, sub.wp_ucomp.values
    h1=axes[n].scatter(x1, y1, s=90, marker='x', color=(0, 0, 0.5),lw=2)
    h2=axes[n].scatter(x2, y2, s=100, marker='+', color=(0, 0.5, 0),lw=2)
    axes[n].text(0.05,0.95,axes[n].get_gid(),size=14,va='top',
            weight='bold',transform=axes[n].transAxes,
            backgroundcolor='w',clip_on=True)
    
h3=axes[0].lines[0]
h4=axes[0].lines[1]
axes[0].legend([h1,h2,h3,h4],
                ['Surface','Wprof','M95',r'$\pm$50% M95 $C_D$'],
                scatterpoints=1)
axes[0].set_ylabel('BBY zonal wind [m s-1]',ha='right')
axes[2].set_xlabel('Pressure difference, BBY-SCK [hPa]',ha='left')
axes[0].set_xticks([])
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[3].set_yticks([])
plt.subplots_adjust(wspace=0.1, hspace=0.1)

plt.show()

#fname='/home/raul/Desktop/gap_flow.png'
#plt.savefig(fname, dpi=300, format='png',papertype='letter',
#            bbox_inches='tight')