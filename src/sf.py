# sf.py
from scipy.interpolate import interp1d
import numpy as np

class sf:
    '''
    Atomic scattering factors.
    '''

    def __init__(atom, pressure, temperature):
        self.atom = atom
        self.pressure = pressure 
    def xuv_index(self, eV):
        '''
        Based on atomic scattering factors from LBNL.

        returns a function for the index of refraction for a given photon energy.

        '''

        nrg, f1, f2 = np.genfromtxt('sf/' + self.atom + '.nff', dtype = float, skip_header = 1, 
                                    usecols = (0,1,2), delimiter = '\t', unpack = True)

        re = 2.8179 * 10 ** -15 #classical electron radius
        kb = 1.3806488 * 10 ** -23 #Boltzmann constant
        f1_interp = interp1d(nrg, f1)
        f2_interp = interp1d(nrg, f2)
        wl = 1240 / eV * 10 ** -9
        dens = self.pressure/kb/self.temperature #density
        return 1 - re * wl ** 2 / 2 / np.pi * dens * (f1_interp(eV) + 1j * f2_interp(eV))

    def xuv_absorption(self, eV):
        '''
        Based on atomic scattering factors from LBNL

        returns the absorption crossection for a given photon energy.
        '''
        nrg, f2 = np.genfromtxt('sf/' + self.atom + '.nff', dtype = float, skip_header = 1, 
                                    usecols = (0,2), delimiter = '\t', unpack = True)
        re = 2.8179 * 10 ** -15
        f2_interp = interp1d(nrg, f2)
        wl = 1240 / eV * 10 **-9
        return 2 * wl * re * f2_interp(eV)