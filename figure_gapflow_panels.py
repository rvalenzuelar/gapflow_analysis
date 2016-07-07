
import seaborn as sns
import gapflow
import matplotlib.pyplot as plt
import Windprof2 as wp
import xpol_tta_analysis as xta
import numpy as np
import os
from datetime import timedelta

homed = os.path.expanduser('~')
# homed = '/localdata'


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
    fig, axes = plt.subplots(2, 3, figsize=(8.5, 6),
                             sharex=True, sharey=True)
axes = axes.flatten()

try:
    tta
except NameError:
    tta = np.array(xta.get_tta_dates([2003,2004]))
    tta = tta + timedelta(minutes=55)
    
n = 0
while n < 6:
    gf = gapflow.run(n+8, plot=True, grid=False, ax=axes[n], homedir=homed,
                     color_surf=(0, 0, 0, 0),color_wp=(0, 0, 0, 0))
    
#    tta = wp.get_tta_times(case=str(n+8), homedir=homed)
#    sub = gf.loc[tta + timedelta(minutes=55)]
    
    
    subdates=get_subdates(gf,tta)
    sub = gf.loc[subdates]
    x1, y1 = sub.pdiff.values, sub.ucomp.values
    x2, y2 = sub.pdiff.values, sub.wp_ucomp.values
    h1=axes[n].scatter(x1, y1, s=90, marker='x', color=(0, 0, 0.5),lw=2)
    h2=axes[n].scatter(x2, y2, s=100, marker='+', color=(0, 0.5, 0),lw=2)
    n += 1
    
axes[0].legend([h1,h2],['Surface','Wprof'],scatterpoints=1)
axes[0].set_ylabel('BBY zonal wind [m s-1]')
axes[4].set_xlabel('Pressure difference, BBY-SCK [hPa]')

plt.subplots_adjust(wspace=0, hspace=0,
                    left=0.1, right=0.98,
                    bottom=0.1, top=0.95)

plt.show()
