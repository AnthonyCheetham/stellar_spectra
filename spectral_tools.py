# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:27:26 2016

Module containing an assortment of tools needed to use the spectra
(e.g. for fitting, )

@author: cheetham
"""

import sys,os,gzip,bz2
from astropy import units
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate,optimize
from astropy import units


##################

def bin_spectrum(wavs_in,flux_in,wav_start,n_wav,delta_wav):
    ''' Bins a high res spectrum to a low res version by simply taking the mean
    of the fluxes contained in each wavelength bin
    wav_start is the middle of the first bin
    n_wav is the number of bins
    delta_wav is the width of each bin
    '''
    
    # Initialise some arrays
    wavs_out = [] # the output wavelengths
    flux_out = [] # the output flux

    # Loop through the wavelength channels
    for wav_ix in range(n_wav):
        
        wav = wav_start + wav_ix*delta_wav
        minwav = wav - delta_wav/2.
        maxwav = wav + delta_wav/2.
        wavs_out.append(wav)
        
        # Find the relevant wavelengths of the input array
        relevant_ix = (wavs_in > minwav) & (wavs_in < maxwav)
        
        # Take the mean of the fluxes over this range
        flux = np.mean(flux_in[relevant_ix])
        flux_out.append(flux)
    
    # Turn into numpy arrays and return
    wavs_out = np.array(wavs_out)
    flux_out = np.array(flux_out)
    
    return [wavs_out,flux_out]
    
##################

def bin_spectrum_gaussian(wavs_in,flux_in,wavs_out,spectral_fwhms):
    ''' Bins a high res spectrum to a low res version by convolving with a Gaussian
    wavs_in, flux_in : the wavelengths and fluxes of the input spectrum
    wavs_out : the output wavelengths
    spectral_fwhms : the FWHMs of the spectral PSF
    '''
    
    flux_out = []
    
    if np.size(wavs_out) ==1:
        wavs_out =[wavs_out]
    
    # Loop through the wavelength channels
    for ix,wav in enumerate(wavs_out):
                
        # Convert the FWHM into the parameter c
        if np.size(spectral_fwhms) >1:
            c= spectral_fwhms[ix]/2.35482
        else:
            c= spectral_fwhms/2.35482
        gauss = np.exp(- (wavs_in - wav)**2 / (2*c**2)) / (c*np.sqrt(2*np.pi))
        gauss /= np.sum(gauss)
        
        # Take the mean of the fluxes over this range
        flux = np.nansum(flux_in*gauss)
        flux_out.append(flux)
    
    # Turn it into a numpy array and return
    flux_out = np.array(flux_out)
    
    return flux_out

##################

def bin_spectrum_filter_curve(wavs_in,flux_in,filter_wavs,filter_curve):
    ''' Bins a high res spectrum to a low res version using a known filter
        transmission function
    wavs_in, flux_in : the wavelengths and fluxes of the input spectrum
    wavs_out : the output wavelengths
    spectral_fwhms : the FWHMs of the spectral PSF
    '''
    
    # Make an interpolation function for the filter curve
    filter_interp = interpolate.interp1d(filter_wavs,filter_curve,fill_value=0,
                                         bounds_error=False)
    
    # Put the filter curve on the same grid as the input array
    filter_func = filter_interp(wavs_in)
    
    # Make sure the filter curve has unit area
    filter_func = filter_func/np.sum(filter_func)
        
    # Sum them together
    flux = np.nansum(flux_in*filter_func)
    
    return flux

##################

def bin_spectrum_all(wavs_in,flux_in,wavs_out,spectral_fwhms,methods):
    ''' This is a wrapper for the gaussian and filter spectral binning
    functions, that allows for a combination of each function for different
    wavelength ranges.
        wavs_in, flux_in : the wavelengths and fluxes of the input spectrum
        wavs_out : the output wavelengths
        spectral_fwhms : the FWHMs of the spectral PSF or a list = [filter_wavs,filter_transmission]
        method: an array of length = len(spectral_fwhms) where each entry
            is either 'Gaussian' or 'filter', so we know which function to use
    '''
    
    # output array for the fluxes
    flux_out = []
    
    # Loop through wavelenths
    for ix,wav in enumerate(wavs_out):
        
        this_method = methods[ix]
        
        if this_method.lower() == 'gaussian':
            flux = bin_spectrum_gaussian(wavs_in,flux_in,wav,spectral_fwhms[ix])
            
        elif this_method.lower() == 'filter':
            filter_wavs,filter_curve = spectral_fwhms[ix]
            flux = bin_spectrum_filter_curve(wavs_in,flux_in,filter_wavs,filter_curve)
        else:
            print('Unknown spectral binning method: '+str(this_method))
            raise Exception
        
        flux_out.append(flux)
    flux_out = np.array(flux_out)
    flux_out = flux_out.ravel()
    
    return flux_out
        

##################

def scale_spectral_flux(wavs_in,flux_in,filter_file,mag,mag_zeroflux, interp_order=3):
    '''
    Calculates the scaling factor for a synthetic spectrum using the flux of the star.
    wavs_in, flux_in : the wavelengths and fluxes of the spectrum
    filter_file : the file containing the filter curve for the photometric standard.
        It is assumed that the first column/row is the wavelengths and the second is the throughput
        It is also assumed that the units of wavs_in and the wavelength column are the same
    mag : the magnitude of the star
    mag_zeroflux : the zero-point flux of the photometric standard
    '''
    
    # Load the 2MASS filter curve
    # From Cohen et al 2003:
    # " these curves are designed to be integrated directly over stellar 
    #   spectra in Fλ form, in order to calculate synthetic photometric 
    #   magnitudes. "
    filter_curve=np.loadtxt(filter_file)
    # In case the data are stored in columns in the file, we need to transpose it
    if filter_curve.shape[1] < filter_curve.shape[0]:
        filter_curve = filter_curve.transpose()
        
    # Does it need to be normalised? I guess not...
    # Normalise it so that the total throughput is 1. Use the trapezoidal rule
#    filter_curve[1] /= np.trapz(filter_curve[1],filter_curve[0])
    # And multiply it by the central wavelength to get the right units back...
#    filter_curve[1] /= np.sum(filter_curve[0]*filter_curve[1])
#    filter_curve[1] /= np.sum(filter_curve[1])
    
    # Make an interpolation function
    filter_curve_interp=interpolate.interp1d(filter_curve[0],filter_curve[1],
                     kind=interp_order,fill_value=0.,bounds_error=False)
                     
    # Calculate the expected J band star flux (in units of energy / area / wavelength)
    star_flux = 10**(-mag/2.5)*mag_zeroflux
    
    # Calculate the J band flux that the input model has
    # Then divide the expected by the model to get the scaling factor
    model_flux = np.trapz(flux_in*filter_curve_interp(wavs_in),wavs_in)
#    model_flux = np.sum(flux_in*filter_curve_interp(wavs_in))
    flux_scale_factor = star_flux / model_flux
    
    # To check, calculate the J band flux of the scaled spectrum  in the
    #  same way and convert it back to magnitudes
    mag_out = -2.5*np.log10(np.trapz(flux_in*filter_curve_interp(wavs_in)*flux_scale_factor,wavs_in)/mag_zeroflux)
#    mag_out = -2.5*np.log10(np.sum(flux_in*filter_curve_interp(wavs_in)*flux_scale_factor)/mag_zeroflux)
    print 'Output and expected mag:',mag_out,mag
    
    return flux_scale_factor
    

#def spectrum_chi2(flux_factor,flux = 0., flux_uncerts=1.,
#                      template_spectrum=0.):
#    ''' Function for scipy.optimize.minimize to minimize when
#    fitting two spectra to each other with a single unknown parameter:
#    the flux scaling between them.'''
#    resids = (flux) - (template_spectrum*flux_factor)
#    chi2 = np.sum((resids /flux_uncerts)**2)
#    return chi2

##################

def spectrum_chi2(flux_factor,flux = 0., flux_uncerts=1.,
                      template_spectrum=0.,spectral_fwhms=1.):
    ''' Function for scipy.optimize.minimize to minimize when
    fitting two spectra to each other with a single unknown parameter:
    the flux scaling between them.
    This version weights each measurement by the spectral FWHM,
    so that broadband measurements contribute more to the fit
    '''
    resids = (flux) - (template_spectrum*flux_factor)
    weights = spectral_fwhms / np.nansum(spectral_fwhms)
    chi2 = np.nansum(weights*(resids /flux_uncerts)**2)
    return chi2

##################

def fit_to_flux(template_spectrum,flux,flux_uncerts,spectral_fwhms,return_chi2=False):
    ''' 
    Scale a template spectrum to match a set of flux datapoints with uncertainties.
    
    template_spectrum: a set of fluxes to be scaled
    flux: a set of measured fluxes
    flux_uncerts: the uncertainties on the measured fluxes
    spectram_fwhms: the FWHM of the filter used to weight each measurement
    
    Each element of template_spectrum must be equivalent to the same element in 
    flux and flux_uncerts (i.e. same filter and units)
    
    '''
    x0 = [np.nanmax(flux)/np.nanmax(template_spectrum)]
    args = (flux,flux_uncerts,template_spectrum,spectral_fwhms)
    fitting_result = optimize.minimize(spectrum_chi2,x0,args,method='Nelder-Mead')
    flux_factor = fitting_result['x']
    
    if return_chi2:        
        return flux_factor,fitting_result['fun']
    else:
        return flux_factor



##################

def scale_spectral_model_to_flux(spectral_wavs,spectral_flux,photometric_mags,
                                 photometric_mags_uncerts,photometry_sources,
                                 n_montecarlo=1000,plot=True):
    ''' Fits a spectral model to a set of observed photometric points.
    For each photometric point you have to define an uncertainty and a source.
    For each source, a filter curve and zeropoint must be hard-coded here.
        Filter curves must be 2-column files with Wavelength (um), Transmission (fraction) as columns.
        Filter curve files must be named the same as their source name +'.txt'
        Zeropoints must be in F_lambda (W m**-2 um**-1)
        
    n_montecarlo = number of realisations used to convert the uncertainties in magnitudes
            into uncertainties in flux
    
    '''
    module_dir = os.path.dirname(os.path.abspath(__file__))+os.sep
    
    # Convert everything into fluxes and get the information needed for bin_spectrum_all
    wavs_out = []
    spectral_fwhms = []
    spectral_data = []
    methods = []
    star_fluxes = []
    star_fluxes_uncerts = []
    
    for ix,source in enumerate(photometry_sources):
        
        source = source.lower()
        
        #### 2MASS #####
        # Taken from Cohen+2013 
        # i.e. http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
        if source == '2mass_j':
            zeropoint = 3.129e-13 #± 5.464E-15 (W cm**-2 um**-1)
            zeropoint *= units.watt / units.cm**2 / units.micron
            filter_fwhm = 0.162 # um
        elif source == '2mass_h':
            zeropoint = 1.133e-13 #± 2.212E-15 (W cm**-2 um**-1)
            zeropoint *= units.watt / units.cm**2 / units.micron
            filter_fwhm = 0.251
        elif source == '2mass_k':
            zeropoint = 4.283e-14 #± 8.053E-16 (W cm**-2 um**-1)
            zeropoint *= units.watt / (units.cm**2) / units.micron
            filter_fwhm = 0.262
            
        #### Hipparcos #####
        
        #### CIT #####
        # These are actually fake filter curves that I made up. 
        # They are top-hat functions with the right bandwidth
        elif source == 'cit_j':
            zeropoint = 3.13429253067e-09
            zeropoint *= units.watt / (units.m**2) / units.micron
            filter_fwhm = 0.240
            


        #### TYCHO #####
        # Taken from the SVO website
        # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=TYCHO/TYCHO.V
        elif source == 'tycho_v':
            zeropoint = 3.984e-9 
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 0.1665
        elif source == 'tycho_b':
            zeropoint = 6.589e-9 
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 0.1455
        
        #### WISE #####
        #  Also taken from the SVO website
        elif source == 'wise_w1':
            zeropoint = 8.238e-12
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 0.6626 # um
        elif source == 'wise_w2':
            zeropoint = 2.431e-12
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 1.1073 # um
        elif source == 'wise_w3':
            zeropoint = 6.570e-14
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 5.505571 # um
        elif source == 'wise_w4':
            zeropoint = 4.995e-15
            zeropoint *= units.erg/(units.cm**2)/units.s/units.angstrom
            filter_fwhm = 4.1016 # um
        
        #### Unknown #####
        else:
            raise IOError('Unknown photometric source: '+str(source))
            
        ##################
        
        # Convert the flux
        zeropoint = zeropoint.to( units.watt / (units.m**2) / units.micron )
        mags_mc = np.random.normal(loc=photometric_mags[ix],
                      scale=photometric_mags_uncerts[ix],size=n_montecarlo)
        star_flux_mc = 10**(-mags_mc/2.5)*zeropoint.value
        
        star_fluxes.append(np.mean(star_flux_mc))
        star_fluxes_uncerts.append(np.std(star_flux_mc))
        
        # Load the filter curve
        filter_file = module_dir+'filters'+os.sep+source+'.txt'
        filter_wavs,filter_curve = np.loadtxt(filter_file,unpack=True)
        # Add the info needed for the fitting        
        wavs_out.append(np.sum(filter_wavs*filter_curve/np.sum(filter_curve)))
        spectral_fwhms.append(filter_fwhm)
        spectral_data.append([filter_wavs,filter_curve])
        methods.append('filter')
    
    
    # Now bin the model spectrum to give predictions for each filter
    model_flux = bin_spectrum_all(spectral_wavs,spectral_flux,wavs_out,spectral_data,methods)
    
    # Now calculate the flux factor needed to scale the model spectrum to the datapoints
    flux_factor = fit_to_flux(model_flux,star_fluxes,star_fluxes_uncerts,spectral_fwhms)
    
    calibrated_spectrum = spectral_flux*flux_factor
    
    if plot:
        plt.figure(0)
        plt.clf()
        plt.plot(spectral_wavs,calibrated_spectrum,'b')
        plt.errorbar(wavs_out,star_fluxes,yerr=star_fluxes_uncerts,xerr=spectral_fwhms,fmt='ro')
        
        #and the predictions
        plt.plot(wavs_out,model_flux*flux_factor,'^')
        
        plt.xlabel('Wavelength (um)')
        plt.ylabel('Flux (W m-2 um-1)')
        
        plt.xscale('log')
        plt.yscale('log')
        
        plt.figure(3)
        plt.clf()
        plt.plot(wavs_out,(star_fluxes-(model_flux*flux_factor))/star_fluxes_uncerts,'x')
        plt.xlabel('Wavelength (um)')
        plt.ylabel('Residuals (sigma)')
        chi2 = np.sum(((np.array(star_fluxes) - model_flux*flux_factor)/star_fluxes_uncerts)**2)
        print('Chi2: '+str(chi2))
        print('RedChi2: '+str(chi2/(len(star_fluxes)-1)))
    
    return calibrated_spectrum
        
##################
    
def calibrate_spex_spectrum(sp,targ_dist = 41.7,silent=False,
                output_units = (units.W / (units.micron *  units.m**2))):
    ''' Input a spectrum from splat, and this outputs the spectrum
    properly scaled to the distance of the target and in the right units
    
    output_units must be an astropy.units type
    
    The output has the format [wavelengths,flux,is_ok]
    where is_ok is True/False depending on whether the calibration
    was successful.
    '''
    
    is_ok = True
    
    # First, flux calibrate the spectrum to its correct H band mag
    bd_hmag = sp.h_2mass
    if (bd_hmag == '') or (bd_hmag == '--'):
        bd_hmag = 0
        is_ok=False
    else:
        bd_hmag = np.float(bd_hmag)
    
    # Check that it actually has a H band spectrum
    if (sp.wave.value.min() < 1.6) & (sp.wave.value.max() > 1.6):
        # We only need to calibrate if we didn't already do it
        if 'Flux calibrated with 2MASS H filter to an apparent magnitude of '+str(bd_hmag) not in sp.history:
            sp.fluxCalibrate('2MASS H',bd_hmag)
    else:
        is_ok = False
    
    # Now scale the flux by the distance
    bd_dist =sp.distance
    if (bd_dist == '') or (bd_dist == '--'):
        bd_dist = 10.
        is_ok = False
    else:
        bd_dist = np.float(bd_dist)
    distance_factor = (bd_dist/targ_dist)**2
    
    # And convert the units
    unit_conv_factor = sp.flux.unit.to(output_units)
    
    flux_out = sp.flux.value * distance_factor * unit_conv_factor
    
    if (is_ok == False) and (silent == False):
        print('Warning! calibrate_spex_spectrum couldnt calibrate the flux.')
        print('  This is likely due to a missing distance or H magnitude')
    
    return [sp.wave.value,flux_out,is_ok]

