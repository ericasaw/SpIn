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

class Spectra():
    """Initializes the Spectra Class
    """    

    def __init__(self, wavelength, flux, error = False):        
        """Creates a Spectra object

        Args:
            wavelength (array):  Has to be dimensionless and have values in angstroms
            flux (array): Has to be dimensionless
            error (bool, optional): Assumes standard deviation uncertainty. Defaults to False.
        """

        #self.wavelength = wavelength
        #self.flux = flux
        self.voigt_params = None
        self.sub_spectrum = None
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
            _type_: _description_
        """        

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (7, 5), constrained_layout = True)
        ax.plot(self.sub_spectrum.wavelength, self.sub_spectrum.flux, color = 'black', alpha = 0.5)
        ax.plot(self.sub_spectrum.wavelength, self.voigt_params(self.sub_spectrum.wavelength), color = 'blueviolet', lw = 4)

        return ax

    def Fitting(self, line_wavelength, window):        
        """Function to fit a Voigt profile to line center provided and window region

        Args:
            line_wavelength (float): Has to have units of Angstroms.
            window (float): Has to have units of Angstroms.
        """        

        #check if wavelength is within spectral region 
        if line_wavelength.value > self.spectrum.wavelength.value.min() and line_wavelength.value < self.spectrum.wavelength.value.max():
            #Setting the sub region
            sub_region = SpectralRegion(line_wavelength - window, line_wavelength + window)
        
            #extracting the sub_region above
            self.sub_spectrum = extract_region(self.spectrum, sub_region)

            voigt_model = models.Voigt1D(x_0 = line_wavelength, bounds = {'x_0': (line_wavelength.value - 0.5, line_wavelength.value + 0.5)})
            voigt_fitter = fitting.LevMarLSQFitter()

            self.voigt_params = voigt_fitter(model = voigt_model, 
                                         x = self.sub_spectrum.wavelength, 
                                         y = self.sub_spectrum.flux)
        else:
            print("Given wavelength is not within the spectral range.")

    def line_finding(self):        
        """Function to find extreme absorption or emission lines

        Returns:
            _type_: _description_
        """        

        with warnings.catch_warnings():  # Ignore warnings
            warnings.simplefilter('ignore')
            lines = find_lines_derivative(self.spectrum, flux_threshold=0.75)

        return lines

    def continuum_fit(self):

        """Fits the continuum
        """       

        #smoothed_spectrum = savgol_filter(x = self.spectrum.flux.value, window_length=1301, polyorder = 3)
        #medfilt_smoothed = medfilt(volume = self.spectrum.flux.value, kernel_size=151)
        
        sigclip = SigmaClip(sigma = 1.5)

        mask = sigclip(data = self.spectrum.flux, masked = True)

        smoothed_continuum = savgol_filter(x = self.spectrum.flux[~mask.mask], window_length=401, polyorder = 3)
        
        plt.plot(self.spectrum.wavelength, self.spectrum.flux, color = 'black', alpha = .7)
        plt.plot(self.spectrum.wavelength[~mask.mask], self.spectrum.flux[~mask.mask], color = 'red', alpha = .7)
        plt.plot(self.spectrum.wavelength[~mask.mask], smoothed_continuum, color = 'blue', alpha = .5)
        #plt.plot(self.spectrum.wavelength, medfilt_smoothed, color = 'blue', alpha = .5)
        plt.show()
