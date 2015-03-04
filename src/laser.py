import numpy as np

class laser:
    '''
    The laser class, very basic. Takes intensity and pulse width, contains function form of pulse.

    Takes pulse_FWHM in [s], Int in [10**14 W cm**-2], lam is wavelength in [nm], spot is w0 in [um].
    pulse(self, t, strength, sigma) outputs a sin**2 pulse for the normalized strength and sigma.

    It also contains a few functions for general gaussian beam properties. 

    '''
    def __init__(self, pulse_FWHM = 100e-15 , Int = 1, lam = 1.07e-6, spot = 17e-6):
        self.pulse_FWHM = pulse_FWHM
        self.Int = Int
        self.lam = lam
        self.spot = spot
        self.zr = np.pi*self.spot ** 2 / self.lam
        self.p0 = self.Int * np.pi / 2 * self.spot ** 2

    def pulse(self, t):
        '''
        Gaussian laser pulse.
        '''
        p = self.Int * np.exp(-2.77 * (t / self.pulse_FWHM ) ** 2 )
        
        return p

    def pulse_s2(self, t):
        '''
        IMPORTANT:
        only valid for t >= - self.pulse_FWHM  and  t <= self.pulse_FWHM
        '''
        p = self.Int * np.sin( np.pi * (t - self.pulse_FWHM) / 2 / self.pulse_FWHM ) ** 2
        
        return p

    def int_z(self, z):
        '''
        returns intensity as function of distance from focus. 
        '''
        return 2 * self.p0 / np.pi / self.spot ** 2 / (1 + (z / self.zr) ** 2 )

    def d_int_z(self, z):
        '''
        returns dI / dz
        '''
        return -2 * self.p0 / np.pi / self.spot ** 2 * z / (self.zr * (1 + (z/self.zr) ** 2 ) ** 2 )

