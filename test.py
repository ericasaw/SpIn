from SpIn_Class import Spectra
from muler.igrins import IGRINSSpectrum
from astropy import units as u
from matplotlib.pylab import plt


test_spec = IGRINSSpectrum(file="SDCH_20201202_0059.spec_a0v.fits")

our_spec = Spectra(test_spec.wavelength.value, test_spec.normalize().flux.value - 1)
#our_spec.Fitting(16493.3 * u.AA, 10 * u.AA)
#ax = our_spec.plotting_spec()

#plt.show()

lines= our_spec.line_finding()
print(lines)

