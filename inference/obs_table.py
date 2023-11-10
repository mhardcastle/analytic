from astropy.table import Table

rms=70e-6
beam=6
fluxcut=0.5e-3
t = Table.read('source-table.fits')

filter = t['live']
filter &= t['flux']>fluxcut
# try to define the actual observational SB limit
fluxlim=(0.3*rms/40.79)*t['LAS']**2.0

filter &=t['flux']>fluxlim

#filter &= t['flux']*beam**2.0/(t['LAS']**2.0) > 3*rms

t2 = t[filter]

t2.write('obs-table.fits',overwrite=True)
