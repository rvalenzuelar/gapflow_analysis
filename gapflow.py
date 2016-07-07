'''
 Analysis of gap flows as appear in
 Valenzuela & Kingmisll (2014)

 Raul Valenzuela
 raul.valenzuela@colorado.edu
 January, 2016
'''

import Meteoframes as mf
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.patches as patches
import pandas as pd
#import os
import Thermodyn as tm

from datetime import datetime
from matplotlib.path import Path


def run(case, plot=False, grid=True, ax=None, homedir=None,
        color_surf=(0, 0, 0.5, 0.5),color_wp=(0, 0.5, 0, 0.5)):

    f = get_filenames(case,homedir)
    meso = mf.parse_mesowest_excel(f[0])

    if len(f[1]) > 1:
        '  more than one day of obs '
        surf = mf.parse_surface(f[1][0])
        for ff in f[1][1:]:
            surf = surf.append(mf.parse_surface(ff))
    else:
        ' only one day '
        surf = mf.parse_surface(f[1][0])

    ''' resample to 1min so we can find
    mesowest index '''
    surf = surf.resample('1T').interpolate()

    t = get_times(case)
    mpress = meso.loc[t[0]: t[1]]['PMSL'].values
    mesoidx = meso.loc[t[0]: t[1]].index

    if case in [1, 2]:
        bias = 9
    else:
        bias = 0

    spress = surf.loc[mesoidx]['press'].values - bias
    swspd = surf.loc[mesoidx]['wspd'].values
    swdir = surf.loc[mesoidx]['wdir'].values
    ucomp = -swspd*np.sin(np.radians(swdir))

    # blh = get_BLH(case)
    blh = 500
    massPa, massU = mass_eq(air_density=1.24, BLH=blh)
    pressDiff = spress - mpress
    d = {'ucomp': ucomp, 'wspd': swspd, 'wdir': swdir,
         'pdiff': pressDiff, 'Bpress': spress, 'Kpress': mpress}
    gapflow = pd.DataFrame(data=d, index=mesoidx)

    ' removes rows with NaN'
    gapflow = gapflow[np.isfinite(gapflow['Kpress'])]

    ' add wind profiler data at target altitude'
    out = get_windprof(case, gapflow_time=gapflow.index,
                       target_hgt_km=0.5,homedir=homedir)
    wp_wspd, wp_wdir = out
    wp_ucomp = -wp_wspd*np.sin(np.radians(wp_wdir))
    gapflow['wp_ws'] = wp_wspd
    gapflow['wp_wd'] = wp_wdir
    gapflow['wp_ucomp'] = wp_ucomp

    path = make_polygon(massPa, massU)

    gapflow = check_polygon(gapflow, path)
    gapflow['gapflow'] = ((gapflow.poly is True) & (gapflow.wdir <= 120))
#    sub = gapflow[(gapflow.poly is True) & (gapflow.wdir <= 120)]

    if plot:
#        timetxt = 'Case {} {}\nBeg: {} UTC\nEnd: {} UTC'
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 7))
            ax.set_xlabel('Pressure difference, BBY-SCK [hPa]')
            ax.set_ylabel('BBY zonal wind [m s-1]')
        ax.scatter(gapflow['pdiff'], gapflow['ucomp'],
                   color=color_surf, label='surf')
        ax.scatter(gapflow['pdiff'], gapflow['wp_ucomp'],
                   color=color_wp, label='wp')
        # ax.scatter(sub['pdiff'], sub['ucomp'], color='r')
        ax.plot(massPa/100, massU[0], marker=None)
        ax.plot(massPa/100, massU[1], linestyle='--', color='r')
        ax.plot(massPa/100, massU[2], linestyle='--', color='r')
        if grid:
            ax.grid(True)
        ax.set_xlim([-12, 1])
        ax.set_ylim([-20, 15])
        ini = mesoidx[0]
        end = mesoidx[-1]
        date = ini.strftime('%b-%Y ')
        beg = ini.strftime('%d')
        end = end.strftime('%d')
#        ax.text(0.03, 0.76, timetxt.format(str(case).zfill(2),
#                                           date, beg, end),
        if beg == end:
            ax.text(0.03, 0.85,'{} {}'.format(beg,date),       
                    fontsize=12,
                    transform=ax.transAxes)
        else:
            ax.text(0.03, 0.85,'{}-{} {}'.format(beg,end,date),       
                    fontsize=12,
                    transform=ax.transAxes)

    return gapflow


def get_windprof(case, gapflow_time=None, target_hgt_km=None,homedir=None):
    import Windprof2 as wp
#    import os
    from scipy.interpolate import interp1d
    from datetime import timedelta

#    homedir = os.path.expanduser('~')
    out = wp.make_arrays(resolution='coarse',
                         surface=False, case=str(case), period=False,
                         homedir=homedir)
    wspd, wdir, time, hgt = out
    time = np.array(time)

    ' match surface obs timestamp'
    time2 = time - timedelta(minutes=5)

    ' get correspoing target time index '
    index_dict = dict((value, idx) for idx, value in enumerate(time2))
    target_idx = [index_dict[x] for x in gapflow_time]

    ' wp time equivalent to gapflow time '
    target_time = time2[target_idx] + timedelta(minutes=5)

    wspd_target = []
    wdir_target = []
    for tt in target_time:
        idx = np.where(time == tt)[0]
        ws = np.squeeze(wspd[:, idx])
        wd = np.squeeze(wdir[:, idx])
        ' interpolate at target altitude '
        fws = interp1d(hgt, ws)
        fwd = interp1d(hgt, wd)
        new_ws = fws(target_hgt_km)
        new_wd = fwd(target_hgt_km)
        wspd_target.append(new_ws)
        wdir_target.append(new_wd)

    return np.array(wspd_target), np.array(wdir_target)


def check_polygon(df, path):

    x = df['pdiff'].values
    y = df['ucomp'].values

    coords = zip(x, y)
    poly = []
    for c in coords:
        poly.append(path.contains_point(c))

    df['poly'] = pd.Series(poly, index=df.index)

    return df


def make_polygon(X, Y):

    x = np.concatenate((X, X[::-1]))/100.
    y = np.concatenate((Y[2], Y[1][::-1]))

    vertices = zip(x, y)
    npoints = len(vertices)

    codes = [Path.MOVETO]+[Path.LINETO]*(npoints-2)+[Path.CLOSEPOLY]

    path = Path(vertices, codes)
    # patch=patches.PathPatch(path,facecolor='orange')

    return path


def get_BLH(case):
    '''
    retrieve a boundary layer height based on
    a subjective analysis of vertical directional
    shear in wind profiler
    '''
    case = str(case)
    blh = {'1': 100,
           '2': 100,
           '3': 500,
           '4': 500,
           '5': 500,
           '6': 500,
           '7': 500,
           '8': 100,
           '9': 200,
           '10': 150,
           '11': 500,
           '12': 200,
           '13': 450,
           '14': 500}
    return blh[case]


def get_filenames(casenum,basedir):

#    basedir = os.path.expanduser('~')
    b = basedir+'/SURFACE'
    surf_files = {1:   [b+'/case01/KSCK.xls', [b+'/case01/bby98018.met', b+'/case01/bby98019.met']],
                  2:   [b+'/case02/KSCK.xls', [b+'/case02/bby98026.met', b+'/case02/bby98027.met']],
                  3:   [b+'/case03/KSCK.xls', [b+'/case03/bby01023.met', b+'/case03/bby01024.met']],
                  4:   [b+'/case04/KSCK.xls', [b+'/case04/bby01025.met', b+'/case04/bby01026.met']],
                  5:   [b+'/case05/KSCK.xls', [b+'/case05/bby01040.met', b+'/case05/bby01041.met']],
                  6:   [b+'/case06/KSCK.xls', [b+'/case06/bby01042.met']],
                  7:   [b+'/case07/KSCK.xls', [b+'/case07/bby01048.met']],
                  8:   [b+'/case08/KSCK.xls', [b+'/case08/bby03012.met', b+'/case08/bby03013.met', b+'/case08/bby03014.met']],
                  9:   [b+'/case09/KSCK.xls', [b+'/case09/bby03021.met', b+'/case09/bby03022.met', b+'/case09/bby03023.met']],
                  10: [b+'/case10/KSCK.xls', [b+'/case10/bby03046.met', b+'/case10/bby03047.met']],
                  11: [b+'/case11/KSCK.xls', [b+'/case11/bby04009.met']],
                  12: [b+'/case12/KSCK.xls', [b+'/case12/bby04033.met']],
                  13: [b+'/case13/KSCK.xls', [b+'/case13/bby04047.met', b+'/case13/bby04048.met', b+'/case13/bby04049.met']],
                  14: [b+'/case14/KSCK.xls', [b+'/case14/bby04056.met']]
                  }

    return surf_files[casenum]


def get_times(casenum):

    slice_times = {1: [datetime(1998, 1, 18, 0, 56), datetime(1998, 1, 18, 23, 56)],
                   2: [datetime(1998, 1, 26, 0, 56), datetime(1998, 1, 27, 3, 56)],
                   3: [datetime(2001, 1, 23, 0, 0), datetime(2001, 1, 25, 0, 0)],
                   4: [datetime(2001, 1, 25, 0, 0), datetime(2001, 1, 27, 0, 0)],
                   5: [datetime(2001, 2, 9, 0, 0), datetime(2001, 2, 11, 0, 0)],
                   6: [datetime(2001, 2, 11, 0, 0), datetime(2001, 2, 12, 0, 0)],
                   7: [datetime(2001, 2, 17, 0, 0), datetime(2001, 2, 18, 0, 0)],
                   8: [datetime(2003, 1, 12, 0, 0), datetime(2003, 1, 15, 0, 0)],
                   9: [datetime(2003, 1, 21, 0, 0), datetime(2003, 1, 24, 0, 0)],
                   10: [datetime(2003, 2, 15, 0, 0), datetime(2003, 2, 17, 0, 0)],
                   11: [datetime(2004, 1, 9, 0, 0), datetime(2004, 1, 10, 0, 0)],
                   12: [datetime(2004, 2, 2, 0, 0), datetime(2004, 2, 3, 0, 0)],
                   13: [datetime(2004, 2, 16, 0, 0), datetime(2004, 2, 19, 0, 0)],
                   # 13: [datetime(2004, 2, 16, 0, 0), datetime(2004, 2, 17, 6, 0)],
                   14: [datetime(2004, 2, 25, 0, 0), datetime(2004, 2, 26, 0, 0)]}

    return slice_times[casenum]


def mass_eq(air_density=1.24, BLH=500):
    ' Gap Flow based on Mass et al (1995) MWR '

    H = BLH  # [m] height of well-mixed boundary layer
    BLcoeff = 2.8  # boundary layer coef (Deardorff 1972)
    npoints = 100
    delPa = np.linspace(0, -1200, npoints)  # [Pa]
    delX = 100000  # [m] distance btwn gap entrance and BBY
    rho = air_density  # [kg m-3] average air density
    PGF = -(1/rho)*(delPa/delX)

    Cd = np.array([7.5e-3, 0.5*7.5e-3, 1.5*7.5e-3])  # drag coefficient
    umass = []
    for c in Cd:
        K = BLcoeff*c/H
        U2 = (PGF/K) * (1-np.exp(-2*K*delX))
        umass.append(-np.real(np.sqrt(U2)))

    return [delPa, umass]


def air_density():

    Rd = 287.  # [J K-1 kg-1]

    # density values do not vary significantly
    Tv = tm.virtual_temperature(C=stemp, mixing_ratio=smixr/1000.)+273.15
    air_density1 = (spress*100.)/(Rd*Tv)

    Tv = tm.virtual_temperature(C=mtemp, mixing_ratio=smixr/1000.)+273.15
    air_density2 = (mpress*100.)/(Rd*Tv)
