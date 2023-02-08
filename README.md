# Plots for "Tidally induced spiral arm wraps encoded in phase space" and others

This repository contains the code for some of the plots in the paper "[Tidally induced spiral arm wraps encoded in phase space](https://www.aanda.org/articles/aa/full_html/2022/12/aa44064-22/aa44064-22.html)" (see file `fourier_plots.ipynb`) and other plots related to galactic dynamics.

## Notebooks

* `compute_actions.ipynb`: Computes actions from positions and velocities, and studies the correlation between action and age.
* `cumuls.ipynb`: Compares Lz-VR curves for clusters of different ages, using a statistical test.
* `decompositions_actions.ipynb`: Functional PCA applied to the curves in relation to the actions.
* `decomposition_tests.ipynb`: Tests of functional PCA with `skfda`.
* `emd.ipynb`: Similar plots to the ones in `fourier_plots.ipynb` but using Empirical Mode Decomposition and the Hilbert-Huang transform instead of the Fourier transform.
* `fourier_discarded_plots.ipynb`: Discarded plots for "[Tidally induced spiral arms]('https://www.aanda.org/articles/aa/full_html/2022/12/aa44064-22/aa44064-22.html')".
* `fourier_2d.ipynb`: Tests of 2-D Fourier decomposition
* `fourier_plots.ipynb`: Plots for "[Tidally induced spiral arms]('https://www.aanda.org/articles/aa/full_html/2022/12/aa44064-22/aa44064-22.html')".
* `rotation_curves.ipynb`: Fitting parameters of rotation curves.
* `z_data_exploration.ipynb`: Some plots with models and Gaia data involving the three dimensions.
