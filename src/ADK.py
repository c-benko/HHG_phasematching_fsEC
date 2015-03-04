import numpy as np
# My imports
from atom import *
from laser import *
from gas import *


class ADK:

    '''
    Containts two adk calculations, the standard one and the modified one by Tong and Lin.

    In general, they will return the ionization fraction using the cycle averaged rate 
    integrated over the pulse. 

    It also contains the steady-state ionization fraction inside a fsEC. See comments.
    '''

    def __init__(self, Atom='Xe', Intensity = .5, Pulse_FWHM = 120e-15):
        self.Atom = Atom
        self.Intensity = Intensity
        self.Pulse_FWHM = Pulse_FWHM

        # loaded quantities
        self.at = atom(self.Atom, 0 , 0, 0)

    def wbar_adk(self, inten):
        '''
        Standard ADK ionization rate calculation. See "Fundamentals of Attosecond Optics".

        returns in units of [1/s]
        '''

        

        p = self.at.adk_params()
        F = np.sqrt(inten / 355.0) + 1e-24

        w_adk = p['abs_Cnl_sq'] * p['G_lm'] * p['Ip'] / 27.2 * (2 * p['F0'] / F) ** (
            2 * p['n_star'] - p['ang_m'] - 1) * np.exp(- 2 * p['F0'] / 3 / F)

        w_bar_adk = np.sqrt(3 * F / np.pi / p['F0']) * w_adk

        return w_bar_adk * 41.341 * 10 ** 15

    def ion_frac_adk(self):
        '''
        integrate over cycle averaged rate, return ionized fraction. 
        '''

        las = laser(self.Pulse_FWHM, self.Intensity)
        tspan = np.linspace(0, 2 * las.pulse_FWHM, 1000)
        dt = (tspan[2] - tspan[1])
        w = []
        for i in range(len(tspan)):
            w.append(self.wbar_adk(las.pulse(tspan[i])))
        return 1 - np.exp(-1 * sum(w) * dt)

    def wbar_adk_TL(self, inten):
        '''
        Implements the over the barrier ionization correction by 
        Tong and Lin J. Phys. B: At. Mol. Opt. Phys. 38 (2005) 2593–2600
        '''
                
        p = self.at.adk_params()
        F = np.sqrt(inten / 355.0) + 1e-24
        kappa = np.sqrt(2 * p['Ip'] / 27.2)

        w_adk = p['abs_Cnl_sq'] * p['G_lm'] * p['Ip'] / 27.2 * (2 * p['F0'] / F) ** (
            2 * p['n_star'] - p['ang_m'] - 1) * np.exp(- 2 * p['F0'] / 3 / F)

        w_bar_adk = np.sqrt(
            3 * F / np.pi / p['F0']) * np.exp(-1 * p['alpha'] / (p['Ip'] / 27.2) * F / kappa ** 3) * w_adk

        return w_bar_adk * 41.341 * 10 ** 15

    def ion_frac_adk_TL(self):
        '''
        Integrate over cycle averaged rate, return ionized fraction. 
        implements the over the barrier ionization correction by 
        Tong and Lin J. Phys. B: At. Mol. Opt. Phys. 38 (2005) 2593–2600
        '''
        las = laser(self.Pulse_FWHM, self.Intensity)
        tspan = np.linspace(0, 2 * las.pulse_FWHM, 1000)
        dt = (tspan[2] - tspan[1])
        w = []
        for i in range(len(tspan)):
            w.append(self.wbar_adk_TL(las.pulse(tspan[i])))
        return 1 - np.exp(-1 * sum(w) * dt)

    def steady_state(self, frac, kp):
        '''
        returns the ionization fraction inside the fsEC. 
        See Allison et al. PRL 2011.

        For our experiment, kp = 0.2 - 0.3

        '''
        return frac / (frac + kp)


