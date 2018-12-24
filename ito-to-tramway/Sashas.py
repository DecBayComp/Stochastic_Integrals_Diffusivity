"""
A class that treats dx and dy instead of x and y in translocations
"""


import numpy as np

from tramway.helper import Translocations


class Sashas(Translocations):
    def _extract_space(self):
        return np.asarray(self.origins[['dx', 'dy']])

    def _extract_time(self):
        return np.asarray(self.origins['dt'])

    @property
    def space_cols(self):
        return ['x', 'y']

    @space_cols.setter
    def space_cols(self, cs):
        pass

    @property
    def time_col(self):
        return 't'

    @time_col.setter
    def time_col(self, cs):
        pass
