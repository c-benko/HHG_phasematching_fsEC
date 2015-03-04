# Heyl's equation for phase matching pressure versus spot size.
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

# my imports
from atom import *

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator

plt.close('all')

lam = 1070.0e-9 #wavelength
Length = 150e-6 # target length
P0 = 1.0 # 1 bar
temp = 200.0 # in K
q = 17.0 # harmonic order
eV = 17 * 1240 / (lam * 10**9) # harmonic photon energy
w0 = np.linspace(9e-6, 100e-6, 100) #spot size
finesse = 400

at = atom('Xe', lam, P0, temp)
dn, etc = at.eta_crit(eV)

def pm_pressure(eta, lam, P0, temp, eV, w0):
    dn, etc = at.eta_crit(eV)
    pressure = lam ** 2 / (2 * np.pi ** 2 * w0 ** 2 * dn * (1 - eta / etc))
    re = 2.8179 * 10 ** -15 #classical electron radius
    kb = 1.3806488 * 10 ** -23 #Boltzmann constant
    phase = -1 * re * lam * pressure*10**5 / kb / temp * Length * eta
    return pressure, phase

fig, ax = plt.subplots(2, sharex = True)
pressure1, phase1 = pm_pressure(0.1*etc, lam, P0, temp, eV, w0)
pressure2, phase2 = pm_pressure(0.5*etc, lam, P0, temp, eV, w0)
pressure3, phase3 = pm_pressure(0.9*etc, lam, P0, temp, eV, w0)
ax[0].loglog(w0*10**6, pressure1, label = '0.1 x $\eta_{crit}$')
ax[0].loglog(w0*10**6, pressure2, label = '0.5 x $\eta_{crit}$')
ax[0].loglog(w0*10**6, pressure3, label = '0.9 x $\eta_{crit}$')
ax[0].set_ylabel('Pressure [bar]')
ax[0].legend(loc = 'lower left')
ax[0].grid(which = 'minor')
ax[0].set_title('Phase matching pressure')
ax[1].set_title('$\eta_{crit}=$' + "{0:.2f}".format(round(etc,2)) + ', Finesse = ' + str(finesse) + ', 150 $\mu$m Target.')

ax[1].semilogx(w0*10**6, phase1*finesse/np.pi, label = '0.1 x $\eta_{crit}$')
ax[1].semilogx(w0*10**6, phase2*finesse/np.pi, label = '0.5 x $\eta_{crit}$')
ax[1].semilogx(w0*10**6, phase3*finesse/np.pi, label = '0.9 x $\eta_{crit}$')
ax[1].set_xlabel('Beam Radius [$\mu$m]')
ax[1].set_ylabel('Phase [rad]')
ax[1].legend(loc = 'lower left')
ax[1].set_xlim(min(w0*10**6), max(w0*10**6))
ax[1].set_ylim(-40,0)
ax[1].yaxis.set_major_locator(LinearLocator(5))
ax[1].grid(which = 'minor')
ax[0].text(-0.1, 1.15, 'a)', transform=ax[0].transAxes,
      fontsize=24, fontweight='bold', va='top', ha='right')
ax[1].text(-0.1, 1.15, 'b)', transform=ax[1].transAxes,
      fontsize=24, fontweight='bold', va='top', ha='right')
# plt.show() 

w0 = np.linspace(10e-6, 40e-6, 50)
finesse = np.linspace(100,500,50)
phase = np.zeros((len(finesse), len(w0)))
phase2 = np.zeros((len(finesse), len(w0)))


for i in range(len(finesse)):
    for j in range(len(w0)):
        pressure3, phase3 = pm_pressure(0.1*etc, lam, P0, temp, eV, w0[j])
        phase[i,j] = -1*phase3*finesse[i]/np.pi

fig, ax = plt.subplots(2, sharex = True)
ext = [min(w0)*10**6, max(w0)*10**6,min(finesse),max(finesse) ]
ax[0].yaxis.set_major_locator(LinearLocator(5))
image = ax[0].imshow(phase, extent = ext ,
                    aspect='auto', origin='lower',interpolation = 'bicubic')
im = plt.colorbar(image, ax=ax[0])
levels = [np.pi]
CS = ax[0].contour(phase, levels, hold='on', colors = 'grey',
        origin='lower', extent=ext)
CS.collections[0].set_label('$\pi$ rad')
im.set_label('Phase [rad]')
# ax[0].set_xlabel('Beam Radius [$\mu$m]')
ax[0].set_ylabel('Finesse')
ax[0].set_title('$\eta$ = 0.1 x $\eta_{crit}$')
ax[0].legend(loc = 'lower right')


for i in range(len(finesse)):
    for j in range(len(w0)):
        pressure3, phase3 = pm_pressure(0.5*etc, lam, P0, temp, eV, w0[j])
        phase2[i,j] = -1*phase3*finesse[i]/np.pi

ext = [min(w0)*10**6, max(w0)*10**6,min(finesse),max(finesse) ]
image = ax[1].imshow(phase2, extent = ext ,
                    aspect='auto', origin='lower',interpolation = 'bicubic')
im = plt.colorbar(image, ax=ax[1])
levels = [np.pi]
CS = ax[1].contour(phase2, levels, hold='on', colors = 'grey',
        origin='lower', extent=ext, linewidth = 3)
CS.collections[0].set_label('$\pi$ rad')
im.set_label('Phase [rad]')
ax[1].yaxis.set_major_locator(LinearLocator(5))
ax[1].set_xlabel('Beam Radius [$\mu$m]')
ax[1].set_ylabel('Finesse')
ax[1].set_title('$\eta$ = 0.5 x $\eta_{crit}$')
ax[1].legend(loc = 'lower right')

ax[0].text(-0.1, 1.15, 'a)', transform=ax[0].transAxes,
      fontsize=24, fontweight='bold', va='top', ha='right')
ax[1].text(-0.1, 1.15, 'b)', transform=ax[1].transAxes,
      fontsize=24, fontweight='bold', va='top', ha='right')
plt.show()

