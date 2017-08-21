from astropy.table import Table

rms=1e-4
beam=6
fluxcut=1e-2
t = Table.read('source-table.fits')

filter = t['live']
filter &= t['flux']>fluxcut
filter &= t['flux']*5*beam**2.0/(t['LAS']**2.0) > 3*rms

t2 = t[filter]

t2.write('obs-table.fits',overwrite=True)
