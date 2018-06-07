# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:14:37 2016

Module for loading and manipulating the data from the BT-NextGen stellar spectra

@author: cheetham
"""


import gzip,bz2
import glob,os,time
from astropy import units
from astropy.io import fits
import numpy as np



def read_ascii_file(spectrum_file,start_row=0,end_row=-1,blackbody=False,stellar_radius=1.0,
              distance=1.,include_dilution_factor=True):
    """
    Reads in a BT-NextGen spectrum and outputs the wavelengths in microns and 
    flux in W / m**2 / um
    start_row and end_row allow for reading only a subset of the file to save time
    (and to skip headers). But note that this way you cant read in the last row.
    blackbody : Set to true to return the blackbody spectrum instead
    stellar_radius : Stellar radius in solar radii
    distance : distance in pc
    """
    
    if spectrum_file.endswith('.bz2'):
        with bz2.BZ2File(spectrum_file,'r') as myf:
            x=myf.read()
    elif spectrum_file.endswith('.gz'):
        with gzip.open(spectrum_file,'r') as myf:
            x=myf.read()
    elif spectrum_file.endswith('.xz'):
        
        from backports import lzma # This can be tricky to install...
        x = lzma.open(spectrum_file).read()
        
    else:
        raise IOError("Unrecognized file type: "+str(spectrum_file))
        
    start_row = np.int(start_row)
    end_row = np.int(end_row)
    
    data=x.split('\n')
    output=[]
    all_output=[]
        
    for ix,row in enumerate(data[start_row:end_row]):
        # remove duplicate white space, not sure how to do this easily
        temp=row.strip()
        for rep in range(50):
            temp=temp.replace('  ',' ')
        
        # And change Ds to Es for the exponential
        temp=temp.replace('D','e')
    
        # Split it based on white space and append the results to the output array    
        if temp != '':
            split=temp.split(' ')
            all_output.append(split)
            output.append(np.float64(split[0:3]))
    output=np.array(output) # make it a numpy array so it is 2D
    if output.ndim != 2:
        raise Exception("Couldnt read the spectra file!"+spectrum_file)
    
    wavs=output[:,0]*1e-4 # turn to microns
    flux=output[:,1]
    
    if blackbody:
        flux = output[:,2]
    
    DF=-8.  # For "all most recent models"
#    DF=-28.9007901434 # For NextGen T > 5000K
#    DF=-26.9007901434 # For NextGen T < 5000K
    
    # convert flux to ergs/sec/cm**2/A
    flux=10**(flux+DF) # this is the conversion eqn from the website (https://phoenix.ens-lyon.fr/Grids/FORMAT)
    
    # convert to W/m**2/um
    # 1 erg/s = 1e-7 W
    # 1 cm**-2 = 1e4 m**-2
    # 1 A**-1 = 1e4 um**-1
    flux_si = flux * 1e-7 * 1e4 * 1e4 # 
#    flux2=flux*units.erg/units.s/(units.cm**2)/units.angstrom
#    flux_si=flux2.to(units.watt/(units.m**2)/units.micron)
#    print 'median log flux:',np.median(np.log10(flux_si))
    
    # Scale by the distance and stellar radii
    # * (radius/distance)^2 in same units.
    rsol_pc = 2.25461e-8
    if include_dilution_factor:
        dilution_factor = (stellar_radius * rsol_pc / distance)**2
        flux_si*=dilution_factor

    
    return [wavs,flux_si]

def read_fits_file(spectrum_file,start_wav = 0.9, end_wav = 2.3,ext_no =1,
                      return_header = False):
    ''' Read in a BT-SETTL spectrum file in fits format
    ext_no : the extension number for the HDUlist
    start_wav and end_wav : the wavelength range in microns to return
    
    Returns:
        wavs,flux : the wavelengths and fluxes in microns and W / (m2 um)
    '''
    
    # Load the file
    hdulist = fits.open(spectrum_file)
    x=hdulist[ext_no].data
    wavs = x['Wavelength']
    flux = x['Flux']
    
    header = hdulist[ext_no].header
    
    # Get rid of the unneeded rows
    good_rows = (wavs > start_wav) & (wavs < end_wav)
    wavs = wavs[good_rows]
    flux = flux[good_rows]
    
    # Check that the units are already ok
    if x.columns[1].unit != 'W / (m2 um)':
        print(' Warning! Units of the BT-SETTL file are not being corrected!')

    if return_header:
        return wavs,flux,header
    else:
        return wavs,flux

def build_btsettl_grid(spectra_directory,start_wav = 0.9, end_wav = 2.3,
               ext_no =1,start_row = 179799, end_row = 334799,filestr='lte*'):
                   
    ''' Loop through a directory of BT-SETTL spectra files in fits format and
    combine them all into a grid of models
    if nextgen_format = True, it will use read_file instead of read_btsettl_file
        Note that the ASS2009 grids go down to 400K and are read by read_file,
        while the other grids stop at ~1200K and are read by read_btsettl_file
    '''
    
    # Make sure it ends in / 
    if spectra_directory[-1] != os.sep:
        spectra_directory+=os.sep
    
    # Find all of the files
    files = glob.glob(spectra_directory+filestr)
    
    # Get the arrays ready 
    all_wavs = []
    all_fluxes = []
    all_teffs = []
    all_loggs = []
    all_metallicities = []
    
    # For the nextgen format files, we will get all wavelengths in the 
    #  first loop and then look only at the relevant ones after that
#    start_row=0
#    end_row=-1
    
    # Loop through them
    for ix,spectrum_file in enumerate(files):
        
        # Is it FITS or ASCII?
        is_fits = '.fits' in spectrum_file
        
        # Load everything
        if is_fits:
            wavs,flux,header = read_fits_file(spectrum_file,start_wav=start_wav,
                     end_wav=end_wav,ext_no=ext_no, return_header=True)
                         
            # Get the parameters for the model
            teff = header['PHXTEFF']
            logg = header['PHXLOGG']
            if 'PHXFE_H' in header.keys():
                metallicity = header['PHXFE_H']
            else:
                metallicity = 0
            
        else:
            wavs, flux = read_ascii_file(spectrum_file,start_row=start_row,
              end_row=end_row,include_dilution_factor=False)
            
            if ix ==1:
                print 'Loaded the section of the grid from:',np.min(wavs),'-',np.max(wavs),'um'
              
              
            # If it is the first loop, find which rows we need to use next time
#            if ix ==-1:
#                start_row = np.argmin(np.abs(wavs-start_wav))
#                end_row = np.argmin(np.abs(wavs-end_wav))
#                
#                wavs = wavs[start_row:end_row]
#                flux = flux[start_row:end_row]
            
            # Or get rid of the bad wavelengths manually
            good_ix = (wavs > start_wav) & (wavs < end_wav)
            wavs = wavs[good_ix]            
            flux = flux[good_ix]
            # And sort it just in case
            sort_ix = np.argsort(wavs)
            wavs = wavs[sort_ix]
            flux = flux[sort_ix]
            
            # Now get the parameters from the filename
            cleanname = spectrum_file.replace(spectra_directory,'')
            cleanname = cleanname.replace('lte','').split('.BT')[0].replace('a','')
            # This is a bit hacky but it should work. Blame the stupid filenames
            cleanname=cleanname.replace('-',' -')
            cleanname=cleanname.replace('+',' +')
            arr=np.array(cleanname.split(' '),dtype=np.float)
            teff, logg, metallicity = arr[0:3]
            teff*=100 # the filename has T/100
        
        # Save everything
        all_wavs.append(wavs)
        all_fluxes.append(flux)
        all_teffs.append(teff)
        all_loggs.append(logg)
        all_metallicities.append(metallicity)
        
        # Print the progress
        if (ix % 10) == 9:
            print('Done '+str(ix+1)+' of '+str(len(files)))
    
    # Turn them into arrays
    all_wavs = np.array(all_wavs)
    all_fluxes = np.array(all_fluxes)
    all_teffs = np.array(all_teffs)
    all_loggs = np.array(all_loggs)
    all_metallicities = np.array(all_metallicities)
    
    return all_wavs,all_fluxes,all_teffs,all_loggs,all_metallicities
    