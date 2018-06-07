# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:26:50 2016

Module for reading and manipulating the Morley et al. (2012) T/Y dwarf
spectral models and Morley et al. (2014) Y dwarf models

@author: cheetham
"""

import glob
import matplotlib.pyplot as plt

import gzip,bz2,os,sys
from astropy import units,constants
import numpy as np

def read_file(spectrum_file,start_row=2,morley_2014 = False,
              include_dilution_factor = False,radius=1.,distance=10.):
    """
    Reads in a Morley et al. spectrum and outputs the wavelengths in microns 
    flux in W / m**2 / um
    start_row allows for reading only a subset of the file to save time
    (and to skip headers).

    Note that the spectra are NOT corrected for the dilution factor by default
     i.e. (radius/distance)**2 
     
     radius in Jupiter radii
     distance in pc (default 10pc)
     
     if we want to read in the morley_2014 models, the first column is actually
     the frequency in Hz, so we need to convert that
    """
    
    # Make sure start_row is an integer!
    start_row = np.int(start_row)
    
    # Load the spectrum using numpy
    spec = np.loadtxt(spectrum_file,skiprows=start_row)
    wavs = spec[:,0]
    flux = spec[:,1]
    
    if morley_2014:
        wavs  = constants.c.value / wavs # lambda = c/f
        # Convert to microns
        wavs *= 1e6

        
    # Reverse the direction so the wavelengths increase
    wavs = wavs[::-1]
    flux = flux[::-1]
    
    # Convert the spectra from erg/s/cm**2/Hz into W/m**2/um
    # Do this step by step to make sure it is correct
    # First, the erg/s/cm**2 to W/m**2
    conv_factor1 = (units.erg / units.s / (units.cm**2)).to(units.W / (units.m**2))
    # And the factor for 1./Hz into 1./m (F_lambda = F_nu * c /lambda**2)
    conv_factor2 = constants.c.value / (wavs*1e-6)**2
    # And now convert 1/m into 1/um
    conv_factor3 = ((units.m**(-1)).to(units.micron**(-1)))
    
#    print type(conv_factor1), type(conv_factor2),type(conv_factor3)
    flux_si = flux * conv_factor1 * conv_factor2 * conv_factor3
    
    if include_dilution_factor:
        dilution_factor = (radius * units.jupiterRad.to(units.pc) / distance)**2
        flux_si*=dilution_factor

    
    return [wavs,flux_si]

def build_grid(spectra_directory, min_wav_um = 0.9, max_wav_um = 2.4,morley_2014 = False):
    ''' Build a grid from the Morley models in a given directory
    Will load the files, combine them and then make an interpolation 
    function in Teff, g, 
    
    set morley_2014 = True to read the morley et al 2014 files
    
    '''
    # Make sure it ends in / 
    if spectra_directory[-1] != os.sep:
        spectra_directory+=os.sep
    
    # Find all of the files
    files = glob.glob(spectra_directory+'sp*')
    
    # Arrays to store the parameters
    all_teffs = []
    all_gs = []
    all_fseds = []
    
    # Arrays to store the data
    all_wavs = []
    all_fluxes = []
    
    # Loop through them and load the data
    for ix in range(len(files)):
        
        name = files[ix]
        # Get the params from the file name (a bit hacked but it works)
        cleanname = name.replace(spectra_directory,'')
        cleanname = cleanname.strip('sp_t')
        cleanname = cleanname.split('_')[0]
        teff = np.int(cleanname.split('g')[0])
        g = np.int(cleanname.split('g')[1].split('f')[0])
        fsed = np.int(cleanname.split('f')[1].split('h')[0]) # this h allows it to handle the Morley 2014 models as well
        
        g *=100 # the filename has g in SI, but logg is usually in cgs
        
        all_teffs.append(teff)
        all_gs.append(g)
        all_fseds.append(fsed)
        
        # Download the actual data and put it in an array
        wavs,flux_si = read_file(name,start_row = 30e3,morley_2014=morley_2014)
        
        relevant_ix = (wavs > min_wav_um) & (wavs < max_wav_um)
        wavs = wavs[relevant_ix]
        flux_si = flux_si[relevant_ix]
        
        all_wavs.append(wavs)
        all_fluxes.append(flux_si)
        
        if (ix % 10) ==9:
            print('Done '+str(ix+1)+' of '+str(len(files)))
            
    # Turn them into numpy arrays
    all_teffs = np.array(all_teffs)
    all_gs = np.array(all_gs)
    all_fseds = np.array(all_fseds)
    all_wavs = np.array(all_wavs)
    all_fluxes = np.array(all_fluxes)
    
    # return everything
    return [all_teffs,all_gs,all_fseds,all_wavs,all_fluxes]

#def rearrange_grid():
#    # Find the unique parameters used for the grid
#    teffs = np.unique(all_teffs)
#    gs = np.unique(all_gs)
#    fseds = np.unique(all_fseds)
#    
#    n_teff = len(teffs)
#    n_gs = len(gs)
#    n_fseds = len(fseds)
#    n_wavs = all_wavs.shape[1]
#    
#    # Rearrange it into a grid of models
#    fluxes = np.zeros((n_teff,n_gs,n_fseds,n_wavs))


# temp code
#spectra_directory ='/Volumes/Cheetham2/spectra/sulfideclouds/'
#temp = glob.glob(spectra_directory+'sp*')
#wavs,flux_si = read_file(temp[0],start_row=20e3)

#plt.clf()
#plt.plot(wavs,flux_si)