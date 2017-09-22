# analytic

This code runs the analytic models described by Hardcastle (2017 -- to
be submitted). Code used to generate the plots from the paper is
provided and can be used for some additional examples.

## Setup instructions

The code relies on a synchrotron library, pysynch. Clone pysynch from
https://github.com/mhardcastle/pysynch.git
and do `python setup.py install`.

Once pysynch is installed the `analytic` code can be used simply by
making sure that the user's `PYTHONPATH` points to this directory. It
does not have its own installation script at the moment.

## Basic use

All code needs to contain the line

```from solver import Evolve_RG```

`Evolve_RG` is the object that evolves a radio galaxy in a particular
environment. Methods of `Evolve_RG` are used to carry out the
operations the user wants.

You may find it useful to do

```
from constants import *
```

which defines some constants and units.

Initialize the environment:

```
env = Evolve_RG('universal',M500=2e13)
```

Define a set of times over which you wish to solve for the dynamics,
then solve for them:

```
tmax = 200
tv = np.logspace(-5,np.log10(tmax),100)*Myr
Q = 1e38
env.solve(Q,tv)
```

The more time intervals you have, the more accurate (but also slower)
the solution.

After `env.solve`, results are stored as attributes of `env` in numpy
arrays matching the time array `tv`:

* `env.tv` are the time values supplied (s)
* `env.R` is lobe length (m)
* `env.Rp` is lobe width (m)
* `env.m1` is forward Mach number
* `env.mp1` is transverse Mach number 
* `env.vl` is lobe volume (m^3)
* `env.vt` is total volume (m^3)

You may want to calculate synchrotron emission:

```
env.findsynch(150e6)
env.findcorrection((150e6,1.4e9))
```

After this some more attributes are set:

* `env.nu_ref` is the reference frequency supplied (Hz)
* `env.B` is the magnetic field strength (T)
* `env.synch` is the uncorrected synchrotron luminosity at the
  reference frequency (W/Hz)
* `env.freqs` are the frequencies at which corrections were requested (Hz)
* `env.corrs` is an array of corrections with size (number of times, number of freqs).
* `env.corr_synch` is an array of synchrotron luminosities after
  correction with the same size as `env.corrs` (W/Hz)

At any time you may save the output

```
env.save('output.pickle')
```

Later on you can load it again:

```
new_env = Evolve_RG.load('output.pickle')
```

## Tests

The code in `test/test_plots.py` tests most of the key routines. It
takes a few minutes to run.
