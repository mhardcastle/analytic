# Add flux densities

from astropy.table import Table
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from constants import *

t=Table.read('source-table.fits')
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

t['flux']=(1e26*t['l150_rest']*(1+t['z'])/(4*np.pi*cosmo.luminosity_distance(t['z']).to(u.m)**2.0)).value
# factor 2 because D is the half-size
t['LAS']=2*np.sin(t['theta'])*(180*60*60/np.pi)*(t['D']/cosmo.angular_diameter_distance(t['z']).to(u.m)).value
# consistency with LOFAR catalogue
t['Length']=2*np.sin(t['theta'])*t['D']/kpc
t.write('source-table.fits',overwrite=True)
