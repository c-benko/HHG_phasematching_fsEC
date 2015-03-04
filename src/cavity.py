import numpy as np

class cavity:
    def __init__(self, IC, Loss):
        self.IC = IC
        self.Loss = Loss

    def finesse(self, loss):
        '''
        returns the finesse as a function of loss. loss includes input coupler.
        '''
        return 2 * np.pi / loss

    def build_up(self, Ti , loss):
        '''
        returns to the power enhancement

        Ti is transmission of input coupler.
        loss includes input coupler.
        '''
        return Ti * self.finesse(Ti + loss) ** 2 / np.pi ** 2