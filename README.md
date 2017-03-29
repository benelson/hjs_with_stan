# hjs_with_stan


Hot Jupiter population inference using a statistical modeling package [Stan](http://mc-stan.org/).


Dependencies
=======

Besides a standard Python installation with numpy, matplotlib, etc., you will also need...

-[daft](http://daft-pgm.org/), for probabilistic graphical modeling

-[pystan](https://pystan.readthedocs.io/en/latest/), for doing MCMC with Stan (in addition to any of Stan's dependencies)


Things of interest
=======

nb_pgm.ipynb: notebook to generate the probabilistic graphical model figures.

nb_nea_figures.ipynb: notebook to generate the paper figures based on planets queried from the NASA Exoplanet Archive.

nb_eod_figures.ipynb: notebook to generate the paper figures based on planets queried from exoplanets.org.

nb_model_comparison.ipynb: notebook to generate model comparison results.


Data and MCMC runs can be found in the *test* directories. *test_nea* uses NASA Exoplanet Archive data and *test_eod* uses exoplanets.org data. Runs based on a single or combination of datasets are listed under the *test_\<nea or eod\>/sample_\<dataset name\>* directories.

Stan models are listed under the *.stan* file extension.


Citations
=======
Coming soon...
