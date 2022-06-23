from matplotlib import lines
import numpy as np
from astropy import units as u
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import find_lines_derivative
from specutils.manipulation import extract_region
from astropy.nddata import StdDevUncertainty
import warnings
from scipy.signal import savgol_filter
from scipy.signal import medfilt
from astropy.stats import SigmaClip
from scipy import interpolate

class Spectra():
    """Initializes the Spectra Class

        Attributes:
            voight_params (function): Function of the best fit profile of the sub_spectrum.
            spectrum (specutils.Spectrum1D): specutils inherited class. Build a Spectrum1D object based 
            off of wavelength (array) and flux (array). Error (optional) defaults to False.
            sub_spectrum (specutils.Spectrum1D): A zoomed-in region cutout of spectrum based on window size.
            continuum_divide (specutils.Spectrum1D): A continuum-divided spectrum.
            f: the interpolated continuum function
            mask: This is the boolean mask that masks the values that got sigma clipped
    """    

    def __init__(self, wavelength, flux, error = False):        

        #self.wavelength = wavelength
        #self.flux = flux
        self.voigt_params = None
        self.sub_spectrum = None
        self.continuum_divide = None
        self.f = None 
        self.mask = None 
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
        """Plotting function

        Returns:
            axes.Axes object: Plots the sub_spectrum and the best fit Voight profile. 
        """        

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (7, 5), constrained_layout = True, facecolor = 'white')
        ax.plot(self.sub_spectrum.wavelength, self.sub_spectrum.flux, color = 'black')
        ax.plot(self.sub_spectrum.wavelength, self.voigt_params(self.sub_spectrum.wavelength), color = 'blueviolet', lw = 3)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')

        return ax

    def Fitting(self, line_wavelength, window, limit_bound = 10 * u.AA):
        """Function to fit a Voigt profile to line center provided and window region

        Args:
            line_wavelength (float): Has to have units of Angstroms.
            window (float): Has to have units of Angstroms.
            limit_bound (float, optional): Changes line center bounds to fit within a narrower region. Defaults to 10*u.AA.
        """                

        #check if wavelength is within spectral region 
        if line_wavelength.value > self.spectrum.wavelength.value.min() and line_wavelength.value < self.spectrum.wavelength.value.max():
            #Setting the sub region
            sub_region = SpectralRegion(line_wavelength - window, line_wavelength + window)
        
            #extracting the sub_region above
            if self.continuum_divide:
                self.sub_spectrum = extract_region(self.continuum_divide, sub_region)
            else:
                self.sub_spectrum = extract_region(self.spectrum, sub_region)

            voigt_model = models.Voigt1D(x_0 = line_wavelength, bounds = {'x_0': (line_wavelength.value - limit_bound.value, line_wavelength.value + limit_bound.value)})
            voigt_fitter = fitting.LevMarLSQFitter()

            self.voigt_params = voigt_fitter(model = voigt_model, 
                                         x = self.sub_spectrum.wavelength, 
                                         y = self.sub_spectrum.flux)
        else:
            print("Given wavelength is not within the spectral range.")

    def line_finding(self, threshold, continuum = True):
        """Function to find extreme absorption or emission lines.

        Args:
            threshold (float): The treshold a flux value must be above or below to be considered an emission or absorption feature.
            continuum (bool, optional): Set to False if you didn't use the continuum_fit function on the spectra. Defaults to True.

        Returns:
            lines (QTable): Table of emission and absorption lines. Column headers: Line center (line_center), line type (line_type) and 
            index of line center (line_center_index) are stored for each line.
        """                        

        with warnings.catch_warnings():  # Ignore warnings
            warnings.simplefilter('ignore')
            if continuum:
                lines = find_lines_derivative(self.continuum_divide, threshold)
            else:
                lines = find_lines_derivative(self.spectrum, threshold)

        return lines

    def continuum_fit(self, window = 11):
        """Fits the continuum

        Args:
            window (int): A user defined odd integer that defines the window length to perform the smoothing of the spectrum
        
        """        

        sigclip = SigmaClip(sigma = 1.5)

        self.mask = sigclip(data = self.spectrum.flux, masked = True)

        smoothed_spectrum = savgol_filter(x = self.spectrum.flux[~self.mask.mask], window_length = window, polyorder = 3)
        
        self.f = interpolate.interp1d(self.spectrum.wavelength[~self.mask.mask], smoothed_spectrum, fill_value='extrapolate')
        
        flux_divide = (self.spectrum.flux / self.f(self.spectrum.wavelength)) - 1
        
        self.continuum_divide = Spectrum1D(spectral_axis= self.spectrum.wavelength, flux = flux_divide)

    def plotting_continuum(self):    
        """Plots a side by side comparison of the original spectrum with the continuum
           
           Returns:
               fig (matplotlib fig object): a matplotlib  fig object that is defaulted to a figsize of (14, 5)
               axes: A list of axes (axes[0], axes[1]) where axes[0] holds the spectrum and the continuum fit plot and axes[1] hold the continuum divided spectrum plot
        """ 
        #makes a 1 row and 2 column fig, ax variables 
        fig, axes = plt.subplots(1, 2, figsize = (14,5), dpi = 100, facecolor = 'white')
        
        #plotting the spectrum, masked_spectrum and continuum
        axes[0].plot(self.spectrum.wavelength, self.spectrum.flux, color = 'black')
        axes[0].plot(self.spectrum.wavelength[~self.mask.mask], self.spectrum.flux[~self.mask.mask], color = 'red', alpha = .7, label = 'sigma clip')
        axes[0].plot(self.spectrum.wavelength, self.f(self.spectrum.wavelength), color = 'blue', label = 'continuum guess')
        axes[0].legend()
        axes[0].set_title('Original Spec')
        axes[0].set_ylabel('Flux')
        axes[0].set_xlabel('Wavelength')
    
        #plotting the continuum divided spectrum
        axes[1].plot(self.continuum_divide.wavelength, self.continuum_divide.flux, color = 'k')
        axes[1].set_title('Continuum Divided')
        axes[1].set_ylabel('Flux')
        axes[1].set_xlabel('Wavelength')

        return fig, axes
