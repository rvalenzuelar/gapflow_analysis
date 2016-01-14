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
import pandas as pd
import sys

from datetime import datetime

basedir = '/home/rvalenzuela/SURFACE'

def main(argv):

	case=int(argv)
	f = get_filenames(case)
	meso = mf.parse_mesowest_excel(f[0])

	if len(f[1])>1:
		'  more than one day of obs '
		surf1 = mf.parse_surface(f[1][0])
		surf2 = mf.parse_surface(f[1][1])
		surf=surf1.append(surf2)
	else:
		surf = mf.parse_surface(f[1][0])
	' upsample to 1min'
	surf = surf.resample('1T').interpolate()


	t = get_times(case)
	mpress = meso.loc[ t[0]: t[1] ]['PMSL'].values
	mesoidx = meso.loc[ t[0]: t[1] ].index

	spress = surf.loc[mesoidx]['press'].values
	swspd = surf.loc[mesoidx]['wspd'].values
	swdir = surf.loc[mesoidx]['wdir'].values
	ucomp=-swspd*np.sin(np.radians(swdir))

	# print surf.loc[mesoidx]
 	# print meso.loc[ t[0]: t[1] ]

	pressDiff = spress - mpress
	[massPa, massU]  = mass_eq(100)

	fig,ax=plt.subplots(figsize=(8,7))
	ax.scatter(pressDiff, ucomp)
	ax.plot(massPa/100,massU[0],marker=None)
	ax.plot(massPa/100,massU[1],linestyle='--',color='r')
	ax.plot(massPa/100,massU[2],linestyle='--',color='r')
	ax.grid(True)
	ax.set_xlabel('Pressure difference, BBY-SCK [hPa]')
	ax.set_ylabel('BBY zonal wind [m s-1]')
	ax.set_xlim([-10, 10])
	ax.set_ylim([-20, 15])
	ini=mesoidx[0]
	end=mesoidx[-1]
	titletimes=ini.strftime('%Y: ')+ini.strftime('%b-%d %H to ')+end.strftime('%b-%d %H UTC ')
	plt.suptitle('Gap Flow Analysis Case'+str(case).zfill(2)+'\n'+titletimes)
	plt.show(block=False)


def get_filenames(casenum):

	b=basedir
	surf_files={1:   [b+'/case01/KSCK.xls', [b+'/case01/bby98018.met', b+'/case01/bby98019.met']],
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

	slice_times = {  1: [datetime(1998,1,18,0,56), datetime(1998,1,18,23,56)], 
					  2: [datetime(1998,1,26,0,56), datetime(1998,1,27,3,56)], 
					  3: [datetime(2001,1,23,0,0), datetime(2001,1,25,0,0)], 
					  4: [datetime(2001,1,25,0,0), datetime(2001,1,27,0,0)],
					  5: [datetime(2001,2,9,0,0), datetime(2001,2,11,0,0)],
					  6: [datetime(2001,2,11,0,0), datetime(2001,2,12,0,0)],
					  7: [datetime(2001,2,17,0,0), datetime(2001,2,18,0,0)],
					  8: [datetime(2003,1,12,0,0), datetime(2003,1,15,0,0)],
					  9: [datetime(2003,1,21,0,0), datetime(2003,1,24,0,0)],
					10: [datetime(2003,2,15,0,0), datetime(2003,2,17,0,0)],
					11: [datetime(2004,1,9,0,0), datetime(2004,1,10,0,0)],
					12: [datetime(2004,2,2,0,0), datetime(2004,2,3,0,0)],
					13: [datetime(2004,2,16,0,0), datetime(2004,2,19,0,0)],
					14: [datetime(2004,2,25,0,0), datetime(2004,2,26,0,0)]}

	return slice_times[casenum]

def mass_eq(npoints):

	H = 500 # [m]
	BLcoeff = 2.8 # boundary layer coef (Deardorff 1972)
	delPa = np.linspace(450,-900,npoints) # [Pa]
	delX = 100000 # [m] distance btwn gap entrance and BBY
	rho = 1.24 # [kg m-3] average air density
	PGF = -(1/rho)*(delPa/delX)

	Cd = [7.5e-3, 0.5*7.5e-3, 1.5*7.5e-3] # drag coefficient
	umass=[]
	for c in Cd:
		K = BLcoeff*c/H 
		U2 = (PGF/K)* (1-np.exp(-2*K*delX))
		umass.append(-np.real(np.sqrt(U2)))

	return [delPa, umass]

main(sys.argv[1])