To redo the plotting for the paper:

setenv AP /home/mjh/git/analytic/paper
setenv PYTHONPATH /home/mjh/git/analytic

# Fig 1 Pressure profiles

python $AP/prprof.py
pdfcrop prprof.pdf

# Fig 2 Example environments

rm example*.pickle
python $AP/example_environments.py
python $AP/example_env_plots.py
pdfcrop example_env.pdf

# Fig 3 HK13 comparison

python $AP/hk13beta.py
python $AP/plot_pluto_hk13.py
pdfcrop hk13-compare.pdf

# Fig 4 Synchrotron vs time
python $AP/example_envs_highz.py
python $AP/example_envs_no_adiabatic.py
python $AP/example-spectra-plots.py
python $AP/example-spectra-vary-plots.py
pdfcrop example-spectra.pdf
pdfcrop example-spectra-vary.pdf

# Fig 5 remnant

python $AP/example_remnant.py
python $AP/example-tstop-spectra-plots.py
pdfcrop example-tstop-spectra.pdf

# Fig 6 Inverse-Compton LC

python $AP/plot_remnant_ic.py
pdfcrop remnant_ic.pdf

# Fig 7 PDD

on cluster
qsub -t 0-129 $AP/pdd.qsub
then copy over and
cd pdd
python $AP/pdd_plots.py
python $AP/pdd_plots.py highz
pdfcrop lllr.pdf
pdfcrop spix.pdf
pdfcrop highz-lllr.pdf
pdfcrop highz-spix.pdf
cd ..

# Fig 8 histograms
on cluster:
python $AP/make_distribution.py
qsub -t 0-10000 $AP/distribution.qsub
rinse and repeat if necessary
python $AP/update_distribution.py # get key quantities from runs
python $AP/flux_distribution.py # add fluxes
python $AP/obs_table.py

then

python $AP/plot_hists.py
pdfcrop sample-hists.pdf

# Fig 9 remnant fraction
python $AP/remnant_fraction.py
pdfcrop remnant-fraction.pdf

# Fig 10 radio power/jet power

python $AP/plot_ql150.py
python $AP/plot_ql150_lowz.py
pdfcrop plot_ql150.pdf
pdfcrop plot_ql150_lowz.pdf

# Fig 11 radio power/environment

python $AP/plot_l150m500_lowz.py
pdfcrop plot_150m500.pdf

# Fig age-size

python $AP/plot_agesize.py
pdfcrop agesize.pdf
