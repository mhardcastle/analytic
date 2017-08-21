# Add flux densities

from astropy.table import Table
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

t=Table.read('source-table.fits')
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

t['flux']=(1e26*t['l150']*(1+t['z'])/(4*np.pi*cosmo.luminosity_distance(t['z']).to(u.m)**2.0)).value
t['LAS']=(180*60*60)*(t['D']/cosmo.angular_diameter_distance(t['z']).to(u.m)).value/(2*np.pi)
t.write('source-table.fits',overwrite=True)
