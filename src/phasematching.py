# phasematching.py
# Absorption limited case inside and fsEC
# at each value of the intensity, the effect of the cavity is implemented by
# calculating the effect of ionization on "squashing" the pulse.
# The reduced intensity is used to calculation the HHG yield.
# The phase shifts and blue shifts are calculated as a function of single pass
# and multiplied by the finesse to account for the cavity.

# my imports
from atom import *
from gas import *
from laser import *
from cavity import *
from ADK import *

class phase_matching:
    '''
    Contains a calculation of phase matching conditions in the 
    absorption-limited case.

    References:
    E. Constant et al. Phys. Rev. Lett.  82, 1668–1671 (1999).
    T. K. Allison. UC Berkeley Thesis (2010).
    T. K. Allison et al. Phys. Rev. Lett. 107, 183903 (2011).
    S. Hadrich et al. Nat. Photon. 8, 779–783 (2014).

    '''
    def __init__(self, Atom = 'Xe', q = 17, Power = 20 , Pulse_FWHM = 120e-15, 
        Spot = 17e-6, Lam = 1070e-9, Pressure = .1, Length = 0.1e-3, Temperature = 100, 
        Z = 0, IC = .008, Loss = 0.0, SS = 'on', BS = 'on'):

        ## input quantities
        self.Atom = Atom
        self.q = q
        self.Power = Power      
        self.Pulse_FWHM = Pulse_FWHM
        self.Spot = Spot 
        self.Lam = Lam
        self.Pressure = Pressure
        self.Length = Length 
        self.Temperature = Temperature
        self.Z = Z
        self.IC = IC
        self.Loss = Loss
        self.SS = SS
        self.BS = BS

        ## constants
        self.re = 2.8179 * 10 ** -15 #classical electron radius
        self.c = 3 * 10 ** 8 
        self.frep = 154 * 10 ** 6
        self.omega0 = self.c * 2 * np.pi / self.Lam

        ## generated quantities used in functions below
        # atom
        self.at = atom(self.Atom, self.Lam , self.Pressure , self.Temperature)
        
        # abs cross section
        self.cross = self.at.xuv_absorption(self.q * 1240 / (self.Lam * 10 ** 9))
        
        # Ip
        self.Ip = self.at.adk_params()['Ip']
        
        # cavity
        self.cav = cavity(self.IC, self.Loss)
        self.Finesse = self.cav.finesse(self.IC + self.Loss)
        self.BuildUp = self.cav.build_up(self.IC, self.Loss)
        self.Intensity = 2 * self.BuildUp * self.Power / (np.pi * self.Pulse_FWHM * self.frep * self.Spot ** 2) * 10 ** -18

        # gas jet
        self.targ = gas(self.Pressure, self.Temperature, self.Length)
        
        # adk class
        self.adk = ADK(self.Atom, self.Intensity, self.Pulse_FWHM)
        
        # gas
        self.g = gas(self.Pressure, self.Temperature, self.Length)

        # power reduction effects
        if self.BS == 'on':
            self.pred = (self.power_reduction()**2 * self.Finesse / np.pi)
            self.red_Intensity = (1 - self.pred) * self.Intensity
            self.red_adk = ADK(self.Atom, self.red_Intensity, self.Pulse_FWHM)
            self.red_las = laser(self.Pulse_FWHM, self.red_Intensity, self.Lam, self.Spot)
        elif self.BS == 'off':
            self.pred = 0
            self.red_Intensity = self.Intensity
            self.red_adk = ADK(self.Atom, self.red_Intensity, self.Pulse_FWHM)
            self.red_las = laser(self.Pulse_FWHM, self.red_Intensity, self.Lam, self.Spot)
## HHG Parameter Calculations
    def dipole(self, inten):
        '''
        Returns the dipole. 
        '''
         
        Up = 9.33 * inten * (self.Lam * 10 ** 6)**2 
        harm = self.q * 1240 / (self.Lam * 10 ** 9) 
        Icut =  (harm - self.Ip) / (3.14 * 9.33 * (self.Lam * 10 ** 6) ** 2)
        if inten < Icut:
            Aq = ( inten / Icut ) ** 10.6
        else:
            Aq = ( inten / Icut ) ** 4.6
        return Aq

    def abs_length(self, eta):
        '''
        Returns the absorption length in [m].
        '''

        return 1 / ( self.cross *  self.g.density * (1-eta))

    def abs_cav(self):
        '''
        Returns the absorption length in [m] for an XUV beam propagating
        in a residual background pressure in the vacuum chamber.

        The propagation distance is assumed to be 1 M. 

        1 atm of backing pressure is approximately 1.3e-5 to 1.3e-6 atm
        of residual gas.
        '''

        g = gas(self.Pressure * 1.3e-6, 300, 1)
        return np.exp(-1 * self.cross *  g.density * 1 ) 


    def coh_length(self, eta):
        '''
        Returns the coherence length in [m].
        '''
        

        # index difference
        dn, eta_crit = self.at.eta_crit(self.q * 1240 / (self.Lam * 10 ** 9))

        # rayleigh range
        zr = self.red_las.zr

        # wavevector mismatch
        dk_atomic = self.q * 2 * np.pi / self.Lam * dn
        dk_plasma = -1*self.q * self.targ.density * self.re * self.Lam * eta
        dk_gouy = -1*self.q / zr / (1 + self.Z **2 / zr ** 2 )
        dk_dip = 2 * self.at.alpha1 * self.red_Intensity / zr * (self.Z / zr) / (1 + self.Z **2 / zr ** 2 )
    
        return np.pi / abs(dk_atomic + dk_plasma + dk_gouy + dk_dip)

## Cavity effects calculations
    def nonlinear_phaseshift(self, eta):
        '''
        Calculates the phase shift from ionization. Steady state portion is included and not subtracted off
        as in the Allison et al. result.

        See T.K.Allison et al. PRL 2011
        '''

        return -1 * self.re * self.Lam * self.g.density * self.Length * eta

    def power_reduction(self):
        '''
        Calculates the power reduction from ionization. dE/E where E is the field envelope.
        The effect should be small. alpha2 is a factor to crudely account for spatial averaging
        over the focus. 

        This factor can be thought of as additional loss to the cavity as a function of intensity.
        It can be inserted as loss to calculate the new finesse. Doing so reproduces our experimental
        results in Allison et al.

        See T.K.Allison et al. PRL 2011 for more information.
        '''

        alpha2 = .2
        if self.SS == 'on':
            eta = self.adk.steady_state(self.adk.ion_frac_adk(), self.at.kp())
        else:
            eta = self.adk.steady_state(self.adk.ion_frac_adk(), 1e6)

        return (alpha2 * self.Ip * 1.6 * 10 ** -19 * 2 * np.pi * 
                self.g.density * self.Length * (1 - eta) * 
                self.adk.wbar_adk(self.Intensity) / (self.Intensity * 10 **18))

    def blue_shift(self, phi, dt):
        '''
        Calculate the central frequency shift from ionization. The multi-pass effect is incorporated
        by the F/pi factor in front of the length. 

        See Dutin et al. Opt. Lett. 2010
        '''

        return (-1 *  np.diff(phi) / dt ) / self.omega0

## Harmonic Yield Calculation        
    def harmonic_yield(self):
        '''
        Returns the harmonic yield as a function of time. 
        This is required to perform the integrated signal. 
        
        Needs to output index as a function of time as well for the pump.
        Add this. 

        '''

        # time span
        num = 100
        tspan = np.linspace(-2 * self.Pulse_FWHM, 2 * self.Pulse_FWHM, num)
        dt = tspan[2] - tspan[1]

        # parameters
        Lmed = self.targ.length
        Dens = self.targ.density

        # bins for outputs
        wbar = np.array(np.zeros(num))
        harm_yield = np.array(np.zeros(num))
        eta = np.array(np.zeros(num))
        Lcoh = np.array(np.zeros(num))
        phi = np.array(np.zeros(num))
        freq = np.array(np.zeros(num-1))

        # Steady-state ionization fraction
        if self.SS == 'on':
            eta_ss = self.red_adk.steady_state(self.red_adk.ion_frac_adk(), self.at.kp())
        else:
            eta_ss = self.red_adk.steady_state(self.red_adk.ion_frac_adk(), 1e6)

        for i in range(num):
            wbar[i] = self.red_adk.wbar_adk_TL(self.red_las.pulse(tspan[i]))
            eta[i] =  1 - (1 - eta_ss) * np.exp(-1 * sum(wbar) * dt)
            Labs = self.abs_length(eta[i]) 
            Lcoh[i] = self.coh_length(eta[i])
            phi[i] = self.nonlinear_phaseshift(eta[i]) * self.Finesse / np.pi
             

            fac1 = 4 * Labs ** 2 * Dens ** 2 * (1 - eta[i]) ** 2  * self.dipole( self.red_las.pulse(tspan[i]) )
            fac2 = 1 / (1 + 4 * np.pi ** 2 * Labs ** 2 / Lcoh[i] ** 2)
            fac3 = 1 + np.exp(-Lmed / Labs) - 2 * np.cos(np.pi * Lmed / Lcoh[i]) * np.exp(-Lmed / 2 / Labs)
            harm_yield[i] =  fac1 * fac2 * fac3

        freq = self.blue_shift(phi, dt)

        return tspan, self.red_las.pulse(tspan), Lcoh, eta, harm_yield, num, phi, freq, Labs
    
    def int_harmonic_yield(self):
        '''
        Returns the integrated harmonic yield. 
        This is required to perform the integrated signal. 
        '''

        tspan, pulse, Lcoh, eta, harm_yield, num, phi, freq, Labs = self.harmonic_yield()
        dt = tspan[2] - tspan[1]
        return sum(harm_yield) * dt * self.abs_cav(), Lcoh[int(num/2)], eta[-1], Labs, self.targ.length, phi[0], (phi[-1]-phi[0]), max(freq), self.pred, self.BuildUp

    




