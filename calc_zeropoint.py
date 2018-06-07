#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 09:07:49 2017

@author: cheetham
"""

# Calculate the zero-point of a certain filter

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.table import Table
import spectral_tools
from astropy import units

#######
# Load the filter
#filter_file = './filters/cit_j.txt'
#filter_file = './filters/johnson_m.txt'
#filter_file = './filters/2mass_j.txt'
filter_file = './filters/NACO_Mp.dat'
#filter_file = './filters/NIRC2_Lp.dat'
#filter_file = './filters/SPHERE_IRDIS_D_H23.dat'
#filter_file = './filters/SPHERE_IRDIS_D_K12.dat'
#filter_file = './filters/SPHERE_IRDIS_B_Y.dat'
filter_data = np.loadtxt(filter_file,unpack=True)

if filter_data.shape[0] == 3:
#    filter_wavs,filter_curve,temp = filter_data
    filter_wavs,temp,filter_curve = filter_data
    filter_wavs *= 1e-3
else:
    filter_wavs,filter_curve = filter_data

#filter_wavs*=1e-3
#######


# Load the Vega spectrum
t = Table.read('alpha_lyr_stis_008.fits')
vega_wavs_angs = t['WAVELENGTH'] # in Angstroms
vega_flux = t['FLUX'] # Flambda, units ?

flux_conversion = (units.erg/units.second/(units.cm**2)/units.angstrom).to (units.watt/(units.m**2)/units.micron)
vega_flux *= flux_conversion

# Convert to microns
vega_wavs = vega_wavs_angs*1e-4

# Integrate the spectrum
vega_filter_flux = spectral_tools.bin_spectrum_filter_curve(vega_wavs,vega_flux,filter_wavs,
                     filter_curve)
print('Using filter '+filter_file)
print('Flux of vega is '+str(vega_filter_flux))

#######

# Plot Vega and the filter
plt.figure(0)
#plt.clf()
plt.plot(vega_wavs,vega_flux)
plt.ylim(0.,0.1e-7)
yl = plt.ylim()
plt.plot(filter_wavs,0.95*filter_curve*yl[1]/max(filter_curve))
plt.xlim(0.9,np.max([3.0,np.max(filter_wavs)]))
plt.xlim(0.9,5.5)