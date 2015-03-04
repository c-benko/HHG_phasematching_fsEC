# atom class

from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import sys, os
here = os.path.dirname(os.path.abspath(__file__))

class atom:
    '''
    The atom class contains relevent parameters for determining the ionization
    rate in a strong field. Parameters come from 'Fundamentals of Attosecond Optics' and 
    use the ADK formalism. It also contains information for calculation XUV dispersion and XUV 
    photoabsorption. It contains the index of refraction for the driving laser wavelenth.

    Some parameters are stored in other classes but passed in during init. See comments of __init__.
    '''

    def __init__(self, Atom , Lam , Pressure , Temperature ):
        self.Atom = Atom

        #get from laser class, defaults given. 
        self.Lam = Lam
        
        #get from gas class, defaults given.
        self.Pressure = Pressure
        self.Temperature = Temperature

        # intensity dependent phase
        self.alpha1 = 2
        self.alpha2 = 22

        # loaded using numpy 
        nrg, f1, f2 = np.genfromtxt(here + '/sf/' + self.Atom + '.txt', dtype = float, skip_header = 1, 
                                    usecols = (0,1,2), delimiter = '\t', unpack = True)

        # load using pandas
        # d = pd.read_csv('sf/' + self.Atom + '.txt', skiprows = 1, delimiter = '\t')
        # nrg = d.values[:,0]
        # f1 = d.values[:,1]
        # f2 = d.values[:,2]

        self.nrg = nrg
        self.f1 = f1
        self.f2 = f2

    def adk_params(self):
        '''
        ADK parameters. See 'Fundamentals of Attosecond Optics'

        return [F0, n_star, l_star, ang_l, ang_m, abs_Cnl_sq, G_lm, Ip]
        '''
        F0 = {'Xe': 0.84187, 'Kr': 1.04375, 'Ar': 1.24665, 'Ne': 1.99547, 'He': 2.42946}
        n_star = {'Xe': 1.05906, 'Kr': 0.98583, 'Ar': 0.92915, 'Ne': 0.7943, 'He': 0.74387}
        l_star = {'Xe': 0.05906, 'Kr': 0.98583, 'Ar': 0.92915, 'Ne': 0.7943, 'He': 0.74387}
        ang_l = {'Xe': 1, 'Kr': 1, 'Ar': 1, 'Ne': 1, 'He': 0}
        ang_m = {'Xe': 0, 'Kr': 0, 'Ar': 0, 'Ne': 0, 'He': 0}
        abs_Cnl_sq = {'Xe': 3.88241, 'Kr': 4.02548,'Ar': 4.11564, 'Ne': 4.24355, 'He': 4.25575}
        G_lm = {'Xe': 3, 'Kr': 3,'Ar': 3, 'Ne': 3, 'He': 3}
        Ip =  {'Xe': 12.129, 'Kr': 13.99,'Ar': 15.759, 'Ne': 21.564, 'He': 24.587}
        alpha = {'Xe': 9, 'Kr': 9,'Ar': 9, 'Ne': 9, 'He': 7}

        return {'F0':F0[self.Atom], 'n_star': n_star[self.Atom], 'l_star': l_star[self.Atom],
                'ang_l': ang_l[self.Atom], 'ang_m': ang_m[self.Atom], 'abs_Cnl_sq':abs_Cnl_sq[self.Atom],
                'G_lm': G_lm[self.Atom], 'Ip':Ip[self.Atom], 'alpha': alpha[self.Atom]}

    def xuv_index(self, eV):
        '''
        Based on atomic scattering factors from LBNL.

        returns a function for the index of refraction for a given photon energy.

        '''
        re = 2.8179 * 10 ** -15 #classical electron radius
        kb = 1.3806488 * 10 ** -23 #Boltzmann constant
        f1_interp = interp1d(self.nrg, self.f1)
        f2_interp = interp1d(self.nrg, self.f2)
        wl = 1240 / eV * 10 ** -9
        dens = self.Pressure/kb/self.Temperature #density
        return 1 - re * wl ** 2 / 2 / np.pi * dens * (f1_interp(eV) + 1j * f2_interp(eV))

    def xuv_absorption(self, eV):
        '''
        Based on atomic scattering factors from LBNL

        returns the absorption crossection for a given photon energy.
        '''
        re = 2.8179 * 10 ** -15
        f2_interp = interp1d(self.nrg, self.f2)
        wl = 1240 / eV * 10 **-9
        return 2 * wl * re * f2_interp(eV)



    def drive_index(self):
        '''
        Based on Börzsönyi APPLIED OPTICS / Vol. 47, No. 27 / 20 September 2008

        returns the index of refraction of the driving laser for a given wavelength,
        pressure and temperature.
        '''
        B1 = {'Xe': 103701.61 * 10 ** -8, 'Kr': 26102.88 * 10 ** -8, 'Ar': 20332.29 * 10 ** -8, 'Ne': 9154.48 * 10 ** -8, 'He': 4977.77 * 10 ** -8}
        C1 = {'Xe': 12.75 * 10 ** -6, 'Kr': 2.01 * 10 ** -6, 'Ar': 206.12  * 10 ** -6, 'Ne': 656.97 * 10 ** -6, 'He': 28.54 * 10 ** -6}
        B2 = {'Xe': 31228.61 * 10 ** -8, 'Kr': 56946.82 * 10 ** -8, 'Ar': 8.066 * 10 ** -8, 'Ne': 4018.63 * 10 ** -8, 'He': 1856.94 * 10 ** -8}
        C2 = {'Xe': 0.561 * 10 ** -3, 'Kr': 10.043 * 10 ** -3, 'Ar': 1.24665 * 10 ** -3, 'Ne': 5.728 * 10 ** -3, 'He': 7.760 * 10 ** -3}
        wl = self.Lam * 10 ** 6
        return np.sqrt( 1 + ( self.Pressure * 273 / self.Temperature ) * (
                        B1[self.Atom] * wl ** 2 / (wl ** 2 - C1[self.Atom] )  +
                        B2[self.Atom] * wl ** 2 / (wl ** 2 - C2[self.Atom] ) ) )        

    def drive_index2(self):
        '''
        base on The Refractive Indices and Verdet Constants of the Inert Gases, 
        Proc. R. Soc. Lond. A 1960 259, doi: 10.1098/rspa.1960.0237

        '''
        A = {'Xe': 1.366e-3, 'Kr': 8.377e-4 , 'Ar': 5.547e-4 , 'Ne': 1.335e-4 , 'He': 6.927e-5 }
        B1 = {'Xe': 9.02e5, 'Kr': 6.7e5, 'Ar': 5.15e5, 'Ne': 2.24e5, 'He': 2.24e5}
        B2 = {'Xe': 1.81e12, 'Kr': 8.84e11, 'Ar': 4.19e11, 'Ne': 8.09e10, 'He': 5.94e10}
        B3 = {'Xe': 4.89e18, 'Kr': 1.49e18, 'Ar': 4.09e17, 'Ne': 3.56e16, 'He': 1.72e16}
        B4 = {'Xe': 1.45e25, 'Kr': 2.74e24, 'Ar': 4.32e23, 'Ne': 0, 'He': 0}
        B5 = {'Xe': 4.34e31, 'Kr': 5.10e30, 'Ar': 0, 'Ne': 0, 'He': 0}
        wl = self.Lam * 10 ** 10
        return np.sqrt( 1 + A[self.Atom] * (1 + B1[self.Atom] / wl ** 2 + B2[self.Atom] / wl ** 4 + B3[self.Atom] / wl ** 6 + B4[self.Atom] / wl ** 8 + B5[self.Atom] / wl ** 10))

    def eta_crit(self, eV):
        '''
        Critical ionization fraction.
        '''
        re = 2.8179 * 10 ** -15 #classical electron radius
        kb = 1.3806488 * 10 ** -23 #Boltzmann constant
        Natm = 1.013 * 10 ** 5 / kb / self.Temperature

        dn = np.real(self.drive_index() - self.xuv_index(eV) )
        eta_crit = 1 / (1 + Natm * re * self.Lam ** 2 / 2 / np.pi / dn)
        return dn, eta_crit

    def kp(self):
        '''
        Decay rate, see Allison et al. PRL (2011)
        '''
        kp = {'Xe': .08, 'Kr': .2 , 'Ar': .3 , 'Ne': .4 , 'He': .5 }
        return kp[self.Atom]





