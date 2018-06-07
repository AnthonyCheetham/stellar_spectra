#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:16:30 2017

@author: cheetham
"""

import numpy as np
import astropy.io.fits as pf
from astropy.table import Table
import stellar_spectra
from stellar_spectra import spectral_tools
import matplotlib.pyplot as plt

n_samples=10000 # for using Monte Carlo method to calculate the uncertainties

# Model of stellar spectrum
stellar_spec_dir='/Volumes/Cheetham2/spectra/BT-NextGen_M-0.5_a+0.0_hot/'
stellar_spec_file='lte088-4.5-0.5a+0.0.BT-NextGen.7.bz2' # Need Teff=8800K, log(g)=? and [M/H]=0

# SED files
sed_file = '/Users/cheetham/code/hip_65426/stellar_sed/'+stellar_spec_file.replace('.bz2','')+'_test.fits'
vosa_sed_file = '/Users/cheetham/code/hip_65426/stellar_sed/edited_sed.fits'

# Spectrum in contrast, that will be converted to flux
ifs_contrast_file = '/Users/cheetham/code/hip_65426/ifs_contrast_file.dat'

# Locations of the sphere filters
filter_dir = './filters/'
sphere_h_filter_file=filter_dir+'SPHERE_IRDIS_D_H23.dat'
sphere_k_filter_file=filter_dir+'SPHERE_IRDIS_D_K12.dat'

# File to save the output
save_file = './test.dat'

refit_stellar_model = False # For scaling a model of the stellar spectrum
convert_photometry = True # For converting the planet contrast into flux
save_photometry = False

##########
# Get the stellar SED
##########
if refit_stellar_model:

    # Photometry of the star from a few sources
    phot_2mass_j = 6.826
    phot_2mass_j_uncert = 0.018
    phot_2mass_h = 6.853
    phot_2mass_h_uncert = 0.048
    phot_2mass_k = 6.770
    phot_2mass_k_uncert = 0.017
    
    hip_mag=7.0427
    hip_mag_uncert= 0.0011
    
    phot_tycho_b = 7.111
    phot_tycho_b_uncert = 0.015
    phot_tycho_v = 6.996
    phot_tycho_v_uncert = 0.010
    
    # WISE
    phot_wise_w1=6.774
    phot_wise_w1_uncert=0.047
    phot_wise_w2=6.794
    phot_wise_w2_uncert=0.020
    phot_wise_w3=6.783
    phot_wise_w3_uncert=0.016
    phot_wise_w4=6.60
    phot_wise_w4_uncert=0.056
    
    # Load the stellar spectral template
    if stellar_spec_file.endswith('.fits'):
        output = stellar_spectra.bt_nextgen.read_fits_file(stellar_spec_dir+stellar_spec_file,
                                                  start_wav = 0.1,end_wav = 100,ext_no=1)
    else:
        start_row=105000
        output=stellar_spectra.bt_nextgen.read_ascii_file(stellar_spec_dir+stellar_spec_file,
                    start_row=start_row,end_row=start_row+94000,include_dilution_factor=False)
    
    spectral_wavs,spectral_flux=output
    
    
    photometric_mags = [phot_2mass_j,phot_2mass_h,phot_2mass_k,phot_tycho_b,
                        phot_tycho_v,phot_wise_w1,phot_wise_w2,phot_wise_w3,
                        phot_wise_w4]
    photometric_mags_uncerts = [phot_2mass_j_uncert,phot_2mass_h_uncert,
                                phot_2mass_k_uncert,phot_tycho_b_uncert,
                                phot_tycho_v_uncert,phot_wise_w1_uncert,
                                phot_wise_w2_uncert,phot_wise_w3_uncert,
                                phot_wise_w4_uncert]
    photometry_sources = ['2mass_j','2mass_h','2mass_k','tycho_b','tycho_v',
                          'wise_w1','wise_w2','wise_w3','wise_w4']
    
    # Now flux calibrate the spectrum
    calibrated_sed = spectral_tools.scale_spectral_model_to_flux(spectral_wavs,
            spectral_flux,photometric_mags,photometric_mags_uncerts,photometry_sources,
            n_montecarlo=1000,plot=True)
    
    # Fix the plot a bit
    plt.figure(0)
    plt.ylim(1e-17,1e-9)
    plt.xlim(1e-1,1e2)
    
    # Save it out:
    t = Table(data=[spectral_wavs,calibrated_sed])
    t.write(sed_file,format='fits',overwrite=True)
    


if convert_photometry:

    ##########
    # Prepare the companion spectrum in contrast
    ##########
    # We need to fill these arrays with the data we have
    data_wavs = [] # Central wavelengths of each filter. e.g. data_wavs[ix] = 3.8
    methods = [] # The type of the filter be 'filter' or 'Gaussian' for each filter. e.g. methods[ix] = 'filter'
    spectral_data = [] # the data for each filter (FWHM or filter curve). e.g. spectral_data[ix]=[filter_wavs,filter_curve] or spectral_data=0.8

    comp_flux_ratio = [] # The flux ratio for each filter
    comp_flux_ratio_uncerts = [] # the uncertainty in the flux ratio

    # some SPHERE contrasts in different bands
    h2_crat = 3.501e-5
    h3_crat = 4.894e-5
    k1_crat = 1.022e-4
    k2_crat = 1.406e-4
    h2_crat_err = 1.66e-6
    h3_crat_err = 2.653e-6
    k1_crat_err = 2.843e-6
    k2_crat_err = 4.217e-6

    # Make a list with all of the information
    data_wavs.extend([1.5888,1.6671,2.110,2.251]) # microns
    methods.extend(['filter','filter','filter','filter'])
    
    sphere_h_filter_wavs,filter_h2,filter_h3 = np.loadtxt(sphere_h_filter_file,unpack=True)
    sphere_h_filter_wavs *= 1e-3 # convert to microns
    
    sphere_k_filter_wavs,filter_k1,filter_k2 = np.loadtxt(sphere_k_filter_file,unpack=True)
    sphere_k_filter_wavs *= 1e-3 # convert to microns
    
    spectral_data.extend([[sphere_h_filter_wavs,filter_h2],
                          [sphere_h_filter_wavs,filter_h3],
                          [sphere_k_filter_wavs,filter_k1],
                          [sphere_k_filter_wavs,filter_k2]])

    comp_flux_ratio.extend([h2_crat,h3_crat,k1_crat,k2_crat])
    comp_flux_ratio_uncerts.extend([h2_crat_err,h3_crat_err,k1_crat_err,k2_crat_err])

    # Add some IFS data
    ifs_wavs,ifs_crats,ifs_crat_uncerts = np.loadtxt(ifs_contrast_file)
    n_ifs = len(ifs_wavs)
    data_wavs.extend(ifs_wavs)
    methods.extend(np.repeat('gaussian',n_ifs))
    ifs_fwhm = 2*(ifs_wavs[1]-ifs_wavs[0])
    spectral_data.extend(np.repeat(ifs_fwhm,n_ifs))

    comp_flux_ratio.extend(ifs_crats)
    comp_flux_ratio_uncerts.extend(ifs_crat_uncerts)

    #######
    # Load the star SED from the file made above
    #######
    t = Table.read(sed_file)
    sed_wavs = []
    sed_flux = []
    for row in t:
        sed_wavs.append(row[0])
        sed_flux.append(row[1])
    sed_wavs = np.array(sed_wavs)
    sed_flux = np.array(sed_flux)
    
    # Remove the non-interesting part
    good_ix = (sed_wavs> 0.7) & (sed_wavs < 7.0)
    sed_flux = sed_flux[good_ix]
    sed_wavs = sed_wavs[good_ix] #microns
    
    ##########
    # Turn the stellar SED into a measurement in each filter
    ##########
    stellar_flux = spectral_tools.bin_spectrum_all(sed_wavs,sed_flux,
                                 data_wavs,spectral_data,methods)    
    
    ##########
    # Multiply by the contrast to get the companion spectrum
    ##########        
    comp_fluxes = np.array(stellar_flux) * comp_flux_ratio
    comp_flux_uncerts = np.array(stellar_flux) * comp_flux_ratio_uncerts
    

    # Plot the result
    plt.figure(2)
    plt.clf()
    plt.errorbar(data_wavs,comp_fluxes,yerr=comp_flux_uncerts,fmt=',')
    
    # Write it to the output file
    if save_photometry:
            
        header='Wavelength(um)  Flux(Wm-2um-1)  Flux_uncertainty'
        save_data = np.array([data_wavs,comp_fluxes,comp_flux_uncerts]).transpose()
        np.savetxt(save_file,save_data,header=header,fmt='%16.10e')