from SpIn_Class import Spectra
from muler.igrins import IGRINSSpectrum
from muler.igrins import IGRINSSpectrumList
from astropy import units as u
from matplotlib.pylab import plt


test_spec = IGRINSSpectrumList.read("SDCH_20201202_0059.spec_a0v.fits").trim_edges()[11]

our_spec = Spectra(test_spec.wavelength.value, test_spec.flux.value)
#our_spec.Fitting(16493.3 * u.AA, 10 * u.AA)
#ax = our_spec.plotting_spec()

#plt.show()

#lines= our_spec.line_fitting()
#print(lines)
our_spec.continuum_fit()
