# testrun.py
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))
from phasematching import *
import matplotlib.pyplot as plt
from pylab import *

sim = phase_matching('Ar', 100,  88 ,  45e-15, 40e-6, 1070e-9,  1, 20e-3,  300, 0, .015, 0, 'on')
# sim = phase_matching('Xe', 13,  35 ,  300e-15, 22e-6, 1070e-9,  .1, .3e-3,  200, 0, .015, 0, 'on')
close('all')
t, y1 , y2, y3, y4, y5, y6, y7,y8 = sim.harmonic_yield()

fig, ax1 = plt.subplots()
ax1.plot(t, y1/max(y1),'b-', label = 'Pulse: ' + str(max(y1)))
ax1.plot(t, y3,'b--', label = 'Ionization Fraction')
ax1.plot(t, y4/max(y4),'b.-', label = 'Harmonic Yield')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Ionization fraction/normalized intensity', color='b')
plt.legend( loc = 1)
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
plot(t, y2 * 10 ** 3,'r-', label = 'L coh [mm]' )
ax2.set_ylabel('Coherence Length [mm]', color='r')
ax2.grid(False)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.legend(loc = 3)


fig, ax = subplots()
ax.plot(t, y6)

fig, ax = subplots()
ax.plot(t[:-1], y7)
plt.show()
