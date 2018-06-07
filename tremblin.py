# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:29:34 2016

Module to read in the Tremblin 2015 Brown Dwarf model spectra
http://adsabs.harvard.edu/abs/2015ApJ...804L..17T

@author: cheetham
"""

import gzip,bz2,glob,os
from astropy import units
from astropy.io import fits
import numpy as np


def read_file(spectrum_file,min_wav=0.5,max_wav = 2.5,distance=False,
              radius=False):
    ''' Read in a single file with the format of the Tremblin 2015
    brown dwarf models.
    The first column is wavenumber in angstroms
    The second column is F(nubar) in erg/s/cm
    radius is the radius of the companion in Jupiter radii
    
    These are accurate for R=0.1R_sun at 10pc
    You can provide a distance (in pc) or a radius (in R_sun) to
    scale these correctly
    '''
    wavenum,flux_nubar_cgs = np.loadtxt(spectrum_file,unpack=True)
    
    # Now convert everything to SI and microns,
#    F (lambda) = F (nubar) * nubar ** 2
    flux_cgs = flux_nubar_cgs * wavenum**2
    
    # Convert the units
    input_units = units.erg / (units.s * units.cm**3)
    output_units = (units.W / (units.micron *  units.m**2))
    unit_conv = input_units.to(output_units)
    
    flux_si = flux_cgs * unit_conv

    # Convert wavelength to microns    
    wavelength_cm = 1./wavenum
    wavelength_micron = wavelength_cm * (units.cm.to(units.micron))
    
    # Just return the relevant part
    good_ix = (wavelength_micron < max_wav) & (wavelength_micron > min_wav)
    wavelength_micron = wavelength_micron[good_ix]
    flux_si = flux_si[good_ix]
    
    if distance:
        flux_si *= (10./distance)**2
    if radius:
        radius_sol = radius * units.jupiterRad.to(units.solRad)
        flux_si *= (radius_sol/0.1)**2
    
    return wavelength_micron,flux_si

def build_grid(spectra_directory,min_wav=0.5,max_wav=2.5):
    '''
    Iterates through a directory of Tremblin models and puts them all
    together into a single array, with the relevant parameters
    '''

    # Find all of the files    
    files = glob.glob(spectra_directory+'spec*')
    
    # Arrays to store the info
    all_wavs = []
    all_fluxes = []
    all_teffs = []
    all_loggs = []
    all_kzzs = []
    all_metallicities = []
    all_gammas = []
    
    # Loop through them
    for ix,spectrum_file in enumerate(files):
        
        # Load it
        wavs,flux = read_file(spectrum_file,min_wav=min_wav,max_wav=max_wav)
        
        # Get the important info from the header
        cleanname = spectrum_file.replace(spectra_directory,'')
        # Remove prefix and file extension
        cleanname = cleanname.replace('spec_lr','')
        cleanname = cleanname.replace('.dat','')
        # Remove the symbols for each quantity
        cleanname = cleanname.replace('t','').replace('g','').replace('k','')
        # Get the sign of the metallicity
        cleanname = cleanname.replace('p','+').replace('m','-')
        # Split it
        cleanname = cleanname.split('_')
        
        teff,logg,kzz,metallicity,gamma = cleanname[1:]
        
        # Store the info
        all_teffs.append(teff)
        all_loggs.append(logg)
        all_kzzs.append(kzz)
        all_metallicities.append(metallicity)
        all_gammas.append(gamma)
        all_wavs.append(wavs)
        all_fluxes.append(flux)
    
    # Make them arrays
    all_teffs = np.array(all_teffs,dtype=np.float)
    all_loggs = np.array(all_loggs,dtype=np.float)
    all_kzzs = np.array(all_kzzs,dtype=np.float)
    all_metallicities = np.array(all_metallicities,dtype=np.float)
    all_gammas = np.array(all_gammas,dtype=np.float)
    all_wavs = np.array(all_wavs)
    all_fluxes = np.array(all_fluxes)
        
    # Return everything except the kzzs and gammas, since they are constant...
    return all_wavs,all_fluxes,all_teffs,all_loggs,all_metallicities
    

# Test it
#spectra_directory = '/Volumes/Cheetham2/spectra/Tremblin/k0_p0_g0/'
#files = glob.glob(spectra_directory + 'spec*.dat')
#spectrum_file = files[0]
#
#wavelength_micron,flux_si = read_file(files[0])

#plt.figure(3)
#plt.clf()  
#plt.plot(wavelength_micron,flux_si)

#grid = build_grid(spectra_directory)
#all_wavs,all_fluxes,all_teffs,all_loggs,all_metallicities = grid