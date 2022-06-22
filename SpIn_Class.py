from matplotlib import lines
import numpy as np
from astropy import units as u
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from pyparsing import line
from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import find_lines_derivative
from specutils.manipulation import extract_region
from astropy.nddata import StdDevUncertainty
import warnings
from specutils.fitting import fit_continuum


class Spectra():

    def __init__(self, wavelength, flux, error = False):
        
        '''
        Initilizes the spectra class
        INPUT: 
            wavelength has to be dimensionless and have values in angstroms
            flux also has to be dimensionless and normalized around 0
            error assumes standard deviation uncertainty, defaults to No error
        
        Creates a Spectra object
        '''
        #self.wavelength = wavelength
        #self.flux = flux
        self.voigt_params = {}
        self.sub_spectrum = {}
        #In order to find lines, NaN values can't exist
        nan_vals=np.isnan(flux)
        flux=flux[~nan_vals]
        wavelength=wavelength[~nan_vals]

        if error:
            spec_error = StdDevUncertainty(error)
            self.spectrum = Spectrum1D(spectral_axis= wavelength * u.AA, 
                                       flux = flux * u.dimensionless_unscaled, 
                                       uncertainty = spec_error) 
        else:
            self.spectrum = Spectrum1D(spectral_axis= wavelength * u.AA, 
                                       flux = flux * u.dimensionless_unscaled)
            
    def plotting_spec(self):

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (7, 5), constrained_layout = True)
        ax.plot(self.sub_spectrum.wavelength, self.sub_spectrum.flux, color = 'black', alpha = 0.5)
        ax.plot(self.sub_spectrum.wavelength, self.voigt_params(self.sub_spectrum.wavelength), color = 'blueviolet', lw = 4)

        return ax

    def Fitting(self, line_wavelength, window):
        '''
        Function to fit a Voigt profile to line center provided and window region

        INPUT:
            line_wavelength has to have units of Angstroms
            window also has units of Angstroms
        '''
        #check if wavelength is within spectral region 
        if line_wavelength.value > self.spectrum.wavelength.value.min() and line_wavelength.value < self.spectrum.wavelength.value.max():
            #Setting the sub region
            sub_region = SpectralRegion(line_wavelength - window, line_wavelength + window)
        
            #extracting the sub_region above
            self.sub_spectrum[f'{line_wavelength.value}'] = extract_region(self.spectrum, sub_region)

            voigt_model = models.Voigt1D(x_0 = line_wavelength, bounds = {'x_0': (line_wavelength.value - 0.5, line_wavelength.value + 0.5)})
            voigt_fitter = fitting.LevMarLSQFitter()

            self.voigt_params[f'{line_wavelength.value}'] = voigt_fitter(model = voigt_model, 
                                         x = self.sub_spectrum.wavelength, 
                                         y = self.sub_spectrum.flux)
        else:
            print("Given wavelength is not within the spectral range.")

    def line_finding(self, threshold=1, sub = False):
        '''
        Function to find extreme absorption or emission lines
        INPUT:
            threshold is used for setting significance value to find emission/absorption
        '''
        with warnings.catch_warnings():  # Ignore warnings
            warnings.simplefilter('ignore')
            if ~sub:
                lines = find_lines_derivative(self.spectrum, flux_threshold=threshold)
            else: 
                lines = find_lines_derivative(self.sub_spectrum, flux_threshold=threshold)

        return lines

    def line_loop(self, lines):
        '''
        '''
        for line in lines['line_center']:
            sub_region = SpectralRegion(line - window * 3, line + window * 3)
            self.sub_spectrum[f'{line.value}'] = extract_region(self.spectrum, sub_region)
            self.continuum_fitting()
            self.Fitting(line, sub_region, line)

    def continuum_fitting(self):
        '''
        Function to fit the local continuum around the located absorption or emission line
        '''
        #add sigma clipping in the future for threshold
        sig_lines = self.line_finding(threshold = 0.5, sub = True)





