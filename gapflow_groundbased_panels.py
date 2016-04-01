import gapflow
import matplotlib.pyplot as plt
import Windprof2 as wp
import os
from datetime import timedelta

homed = os.path.expanduser('~')

fig, axes = plt.subplots(2, 3, figsize=(8.5, 6), sharex=True, sharey=True)
axes = axes.flatten()

n = 8
while n < 14:
    gf = gapflow.run(n, plot=True, grid=False, ax=axes[n-8])
    tta = wp.get_tta_times(case=str(n), homedir=homed)
    sub = gf.loc[tta+timedelta(minutes=55)]
    x1, y1 = sub.pdiff.values, sub.ucomp.values
    x2, y2 = sub.pdiff.values, sub.wp_ucomp.values
    axes[n-8].scatter(x1, y1, s=100, marker='x', color=(0, 0, 0.5))
    axes[n-8].scatter(x2, y2, s=100, marker='+', color=(0, 0.5, 0))
    n += 1
axes[0].legend(scatterpoints=1)
axes[0].set_ylabel('BBY zonal wind [m s-1]')
axes[4].set_xlabel('Pressure difference, BBY-SCK [hPa]')
# axes[7].set_visible(False)
# axes[8].set_visible(False)
plt.subplots_adjust(wspace=0, hspace=0,
                    left=0.1, right=0.98,
                    bottom=0.1, top=0.95)
