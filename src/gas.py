# gas.py

import numpy as np

class gas:
    '''
    Basic gas class that contains parameters for the gas jet.

    Can be extended later to include jet profiles.

    pressure is in [atm]. temperature is in [K], length in [m].
    '''

    def __init__(self, pressure, temperature, length):
        self.pressure = pressure
        self.temperature = temperature
        self.length = length
        self.density = self.pressure * 1.013 * 10 ** 5 / (1.3806488 * 10 ** -23) / self.temperature
        self.pl = self.pressure
     