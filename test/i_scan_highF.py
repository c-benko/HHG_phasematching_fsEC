# i_scan.py
import sys, os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

from phasematching import *
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator



num = 100
iscan = np.linspace(3, 30, num)*.8
pd = 120e-15
spot = 17e-6

harm = np.array(np.zeros(num))
Lcoh = np.array(np.zeros(num))
ion = np.array(np.zeros(num))
Labs = np.array(np.zeros(num))
Lmed = np.array(np.zeros(num))
phiSS = np.array(np.zeros(num))
dphi = np.array(np.zeros(num))
freq = np.array(np.zeros(num))
pred = np.array(np.zeros(num))
buildup = np.array(np.zeros(num))

for i in range(num):
    sim = phase_matching('Xe', 17,  iscan[i] ,  pd, spot, 1070e-9,  .005, 0.10e-3,  300, 0, .008,0, 'on')#.01
    harm[i], Lcoh[i], ion[i], Labs[i], Lmed[i], phiSS[i], dphi[i], freq[i], pred[i], buildup[i] = sim.int_harmonic_yield()


fig, ax = plt.subplots(2, 3, figsize = (15,8))
matplotlib.rcParams.update({'font.size': 16})

ax[0,0].plot(iscan  , harm, 'k-', linewidth = 2)
ax[0,0].set_xlabel('Injected power [W]')
ax[0,0].set_ylabel('Harmonic Yield [arb.]')

ax[0,1].plot(iscan  , Lcoh * 10 ** 3, 'k-', linewidth = 2, label = 'Coherence Length')
ax[0,1].plot(iscan  , Lmed * 10 ** 3, 'r-', linewidth = 2, label = 'Medium Length')
ax[0,1].plot(iscan  , Labs * 10 ** 3, 'b-', linewidth = 2, label = 'Absorption Length')
ax[0,1].set_xlabel('Injected power [W]')
ax[0,1].set_ylabel('Length [mm]')
ax[0,1].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[0,1].legend(loc = 'center right', prop={'size':12})

ax[0,2].plot(iscan  , buildup*iscan*(1-pred), 'k-', linewidth = 2, label = 'Reduced')
ax[0,2].plot(iscan  , buildup*iscan, 'r-', linewidth = 2, label = 'No Gas')
ax[0,2].yaxis.set_major_locator(LinearLocator(5))
ax2 = ax[0,2].twinx()
ax[0,2].set_ylim(min(buildup*iscan*(1-pred)),max(buildup*iscan*(1-pred)))
ax2.plot(iscan  , (2*buildup*iscan*(1-pred)/(np.pi*154e6 * pd * (spot)**2))*10**-18, 'k-')
ax2.set_ylabel('Peak Intensity [10$^{14}$ W cm$^{-2}$]')
ax2.set_ylim(min((2*buildup*iscan*(1-pred)/(np.pi*154e6 * pd * (spot)**2))*10**-18),max((2*buildup*iscan*(1-pred)/(np.pi*154e6 * pd * (spot)**2))*10**-18))
ax2.yaxis.set_major_locator(LinearLocator(5))
ax2.grid()

ax[0,2].legend(loc = 2)
ax[0,2].set_xlabel('Injected power [W]')
ax[0,2].set_ylabel('IC Power [W]')

ax[1,0].plot(iscan  , ion, 'm-', linewidth = 2, label = 'Ionization Fraction')
# ylim(.01,1)
ax[1,0].set_xlabel('Injected power [W]')
ax[1,0].set_ylabel('Ionization Fraction')

ax[1,1].plot(iscan  , Lcoh / Labs, 'k-', linewidth = 2, label = 'L_coh / L_abs ( > 3)')
ax[1,1].plot(iscan  , Lmed / Labs, 'r-', linewidth = 2, label = 'L_med / L_abs ( > 5)')
ax[1,1].set_xlabel('Injected power [W]')
ax[1,1].set_ylabel('Ratio')
ax[1,1].legend(prop={'size':12})

ax[1,2].plot(iscan, phiSS, 'k-', linewidth = 2, label = 'Steady State Phase')
ax[1,2].yaxis.set_major_locator(LinearLocator(5))

ax2 = ax[1,2].twinx()
ax2.plot(iscan, dphi, 'r-', linewidth = 2, label = 'Phase Change')
ax2.set_ylabel('Phase Shift over pulse [rad]', color = 'r')

ax[1,2].set_xlabel('Injected power [W]')
ax[1,2].set_ylabel('Steady State Phase Shift [rad]')
ax2.yaxis.set_major_locator(LinearLocator(5))
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.grid()

fig.set_tight_layout(True)
plt.show()
