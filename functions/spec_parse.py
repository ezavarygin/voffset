#!/usr/bin/env python

import numpy as np
from barak.convolve import convolve_constant_dv
from astropy.io import fits


class Up_parse:
    """
    The class is to parse a UVES_popler output file.
    
    It provides wavelength, flux, sigma_error, mask (valid pixels) arrays.
    """
    def __init__(self, path_to_fits,off_set=0.0):
        shift = 1. + off_set/299792.458 # (1+v/c)
        self.fits_opened = fits.open(path_to_fits)
        self.file_name = path_to_fits.split('/')[-1]
        self.wave   = np.array([self.pix2wave(pix)*shift for pix in range(self.fits_opened[0].data[0].size)])
        self.flux   = self.fits_opened[0].data[0]
        self.error  = self.fits_opened[0].data[1]
        self.disp   = self.fits_opened[0].header['UP_DISP']
        self.mask   = np.array([True if status == 1.0 else False for status in self.fits_opened[0].data[4]])
#        self.status = self.fits_opened[0].data[4] # clipping status in Up
        self.fits_opened.close()
#        self.weight = np.array([1.0 if mask else 1E-10 for mask in self.mask]) # no working in spline?
        self.length = len(self.flux)
        

    def pix2wave(self,pixel):
        """
        Convert a pixel number into a wavelength.

        It is based on the info from the header.
        """
        log_lin = self.fits_opened[0].header['DC-FLAG']
        if log_lin == 1:
            wl_0     = 10**self.fits_opened[0].header['CRVAL1']
            CRPIX1   = self.fits_opened[0].header['CRPIX1']
            log_disp = self.fits_opened[0].header['CD1_1']
            w_i=wl_0*10**(((pixel+1)-CRPIX1)*log_disp)
        else:
            sys.exit('wavelength scale should be log-linear!')
        return w_i

    def fix_clip(self):
        """
        Makes flux values in clipped pixels to be np.nan.
        """
        for i in range(self.length):
            if not self.mask[i]:
                self.flux[i] = np.nan

    def fix_clip_err_nan(self):
        """
        Makes error values in clipped pixels to be np.nan.
        """
        for i in range(self.length):
            if not self.mask[i]:
                self.error[i] = np.nan

    def fix_clip_err_med(self):
        """
        Makes error values in clipped pixels to be mean.
        """
        for i in range(self.length):
            if not self.mask[i]:
                self.error[i] = 1.0
        median = np.median(self.error)
        for i in range(self.length):
            if not self.mask[i]:
                self.error[i] = median

    def convolve(self, res):
        """
        Convolves the spectrum using the Barak package.
        
        res is a path to the file with Gaussian kernel specified with '--fwhm'.
        """
        if res == None:
            self.flux = convolve_constant_dv(self.wave,self.flux,wa_dv=self.wave,npix=4)
        else:
            self.flux = convolve_constant_dv(self.wave,self.flux,vfwhm=self.spec_fwhm_file(res)[self.file_name])

    def spec_fwhm_file(self,res):
        """
        Creates a dictionary with "spectrum name : fwhm".
        
        It takes fwhms from the res file specified with '--fwhm'.
        """
        with open(res,'r') as f:
            res_split = f.read().split('\n')
        fwhm = dict()
        for line in res_split:
            fwhm[line.split()[0]] = float(line.split()[1])
        return fwhm

    def get_mask(self,an_spec):
        """
        Defines an array of valid pixels.
        
        Make a projection of bad pixels in the anchor spectrum onto the one being under analysis.
        The pixel numbers are in frame of the analysed specrum.
        The output array therefore takes into account bad pixels of both anchor and analysed spectra.
        """
        spec_mask = self.mask
        pix2pix = self.pix2an_pix(an_spec)
        an_length = an_spec.length
        if an_spec.file_name == self.file_name:
            pass
        else:
            for i in range(self.length):
                if spec_mask[i]:
                    if pix2pix[i] < 0 or pix2pix[i] >= an_length:
                        spec_mask[i] = False
                    else:
                        if an_spec.mask[pix2pix[i]]:
                            pass
                        else:
                            spec_mask[i] = False
                else:
                    pass
        self.mask = spec_mask

    def pix2an_pix(self,an_spec):
        """
        Defines a relation between pixels in an_spec and the analysed spectrum.
        
        Defines which pixel numbers in an_spec the pixels in the analysed spectrumn corresponds to.
        Returns array of len(self.flux) with corresponding pixel indices of the anchor spectrum.
        """
        an_wl_0     = 10**an_spec.fits_opened[0].header['CRVAL1']
        an_CRPIX1   = an_spec.fits_opened[0].header['CRPIX1']
        an_log_disp = an_spec.fits_opened[0].header['CD1_1']
        an_index = []
        for i,wave in enumerate(self.wave):
            an_index.append(int(round(np.log10(wave/an_wl_0)/an_log_disp + an_CRPIX1 - 1.)))
        return np.array(an_index)
