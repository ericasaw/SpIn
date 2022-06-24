from SpIn_Class import Spectra
import numpy as np

def test_hetdexSpec():

    #reads in the hetdex Spectra
    test_wavelength, test_spec = np.loadtxt('test_data/HDR_LAE.txt', unpack = True)
    our_spec = Spectra(test_wavelength, test_spec)
    
    lines = our_spec.line_finding(1, continuum=False)
    
    emission_lines = lines[lines['line_type'] == 'emission']['line_center'].value
    #print(emission_lines)
    #expected Line Test
    expected_line = 4334.98 #angstroms

    #expected_line in emission_lines
    assert (abs(expected_line - emission_lines) < 1.5).any()

    print('No Errors :D!!')


    #our_spec.Fitting(4334.98  * u.AA, 10 * u.AA)#, 5 * u.AA)
    #fig, ax = our_spec.plotting_continuum()
    #plt.show()
    #our_spec.plotting_spec()
    #plt.show()

test_hetdexSpec()
    