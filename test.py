from SpIn_Class import Spectra
from muler.igrins import IGRINSSpectrum
from muler.igrins import IGRINSSpectrumList
from astropy import units as u
from matplotlib.pylab import plt


test_spec = IGRINSSpectrumList.read("SDCH_20201202_0059.spec_a0v.fits").trim_edges()[14]

our_spec = Spectra(test_spec.wavelength.value, test_spec.flux.value)
our_spec.continuum_fit(window = 401)
lines = our_spec.line_finding(0.1)
print(lines)
our_spec.Fitting(16100.62780311558  * u.AA, 20 * u.AA, 5 * u.AA)
fig, ax = our_spec.plotting_continuum()
plt.show()