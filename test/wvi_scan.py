import sys, os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

from phasematching import *
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator



num = 20
wscan = np.linspace(13e-6, 40e-6, num)
iscan = np.linspace(6, 80, num)

harm = np.array((np.zeros((num, num))))
Lcoh = np.array((np.zeros((num, num))))
ion = np.array((np.zeros((num, num))))
Labs = np.array((np.zeros((num, num))))
Lmed = np.array((np.zeros((num, num))))
phiSS = np.array((np.zeros((num, num))))
dphi = np.array((np.zeros((num, num))))
freq = np.array((np.zeros((num, num))))
pred = np.array((np.zeros((num, num))))
buildup = np.array((np.zeros((num, num))))

for i in range(num):
    for j in range(num):
        sim = phase_matching('Xe', 17,  iscan[i]*.8 ,  120e-15, wscan[j], 1070e-9,  .5, 0.7e-3,  200, 0, .015,0, 'on')
        harm[i,j], Lcoh[i,j], ion[i,j], Labs[i,j], Lmed[i,j], phiSS[i,j], dphi[i,j], freq[i,j], pred[i,j], buildup[i,j] = sim.int_harmonic_yield()

f, ax = plt.subplots()
image = ax.imshow(harm, extent = [min(iscan),max(iscan),min(wscan)*10**6, max(wscan)*10**6 ],
                    aspect='auto', origin='lower',interpolation = 'bicubic')#,vmin=0, vmax=9e28)
plt.colorbar(image, ax=ax)
# ax.set_ylim(min(pscan),max(pscan))
# ax.set_xlim(min(iscan), max(iscan))


plt.show()