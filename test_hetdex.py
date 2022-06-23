from cgi import test
from SpIn_Class import Spectra
from astropy import units as u
from matplotlib.pylab import plt
import numpy as np

test_wavelength, test_spec = np.loadtxt('test_data/HDR_LAE.txt', unpack = True)

our_spec = Spectra(test_wavelength, test_spec)
#our_spec.continuum_fit(window = 11)
#lines = our_spec.line_finding(1, continuum=False)
#print(lines[lines['line_type'] == 'emission'])
#our_spec.Fitting(4334.98  * u.AA, 10 * u.AA)#, 5 * u.AA)
#fig, ax = our_spec.plotting_continuum()
#plt.show()
#our_spec.plotting_spec()
#plt.show()
fig = plt.figure(figsize = (7, 5), facecolor='white')
plt.plot(test_wavelength, test_spec, color = 'black')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.show()