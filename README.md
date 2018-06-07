# stellar_spectra
Python tools for dealing with spectra in near-infrared exoplanet imaging data.

This repository contains various functions for loading and fitting spectra. It was developed for use with exoplanet imaging data. For example:

1. Flux-scaling a model stellar spectrum to various photometric datapoints.

2. Using the model stellar spectrum to convert a set of exoplanet contrast measurements to fluxes.

3. Fitting a set of exoplanet fluxes to a set of models.

The code is split into three main parts: modules for loading spectral models (bt_nextgen.py, tremblin.py and morley_2012.py), a script to calculate the zeropoint flux of a given filter (calc_zerpoint.py) and a module that contains more general tools (spectral_tools.py).

## Loading Spectral Models

* bt_nextgen.py : For loading various spectra from the Lyon and Exeter groups (BT-Nextgen, BT-Settl, BT-Dusty and BT-COND have been tested to work). Should be able to load ascii files and fits files with the read_ascii_file and read_fits_file functions.

* tremblin.py : For loading spectra from Tremblin+2015 (http://cdsads.u-strasbg.fr/abs/2015ApJ...804L..17T).

* morley_2012.py : For loading spectra from Morley+2012 (http://cdsads.u-strasbg.fr/abs/2012ApJ...756..172M), and extended to work on the Morley+2014 spectra (http://cdsads.u-strasbg.fr/abs/2014ApJ...787...78M).

## spectral_tools.py

* Can bin spectra to different resolutions, or calculate the predicted flux for a given filter using bin_spectrum, bin_spectrum_filter_curve, bin_spectrum_gaussian, bin_spectrum_all.

* Scale a model spectrum in flux to match some photometric points using scale_spectral_model_to_flux.

* Calculate the goodness-of-fit for a two spectra being the same using spectra spectrum_chi2 (e.g. a model and a dataset), and calculate the best-fit flux ratio between two spectra using fit_to_flux.

* Flux calibrate spectra from the splat / SpeX Prism library (https://github.com/aburgasser/splat) using calibrate_spex_spectrum. This makes sure it is in the correct units and scaled to the right distance for useful comparison with data.


Copyright (C) 2018 Anthony Cheetham

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.