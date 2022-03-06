##Wallis, 2022
#Please visit openamoc.info or investigate the README file to learn more
from typing import AnyStr, Callable
from collections import defaultdict

import numpy as np

#functions for calling variables
class State(np.ndarray):
    _variables_ = []

    @classmethod
    def from_dict(cls, dictlike):
        return np.asarray([dictlike[varbl]
            for varbl in cls._variables_]).view(cls)

    def __getitem__(self, key):
        if isinstance(key, str):
            key = self._variables_.index(key)
        return np.ndarray.__getitem__(self, key)

    def __str__(self):
        return ''.join((type(self).__name__, '(',
            ', '.join('{0:s}:{1:f}'.format(varbl, self[varbl])
                for varbl in self._variables_),
            ')'))

#functions for shaping variable change
class BalanceBase:
    def __init__(self):
        self.fluxes = list()
        self.inertias = defaultdict(lambda : 1.0)

    def set_inertia(self, **kwargs):
        for varbl, inertia in kwargs.items():
            self.inertias[varbl] = inertia

    def add_flux(self, sink:AnyStr=None, source:AnyStr=None):
        def wrapper(func):
            self.fluxes.append((func, source, sink))
            return func
        return wrapper

    def tend(self, state:State):
        convs = defaultdict(float)
        for flux, source, sink in self.fluxes:
            conv = flux(state)
            if source is not None:
                convs[source] -= conv
            if sink is not None:
                convs[sink] += conv
        result = defaultdict(float)
        for varbl, conv in convs.items():
            result[varbl] = convs[varbl]/self.inertias[varbl]
        return type(state).from_dict(result)

#functions for integration and iteration
class OdeSolver:
    def __init__(self, tend:Callable, dt:float, scheme:AnyStr):
        self.tend = tend
        self.dt = dt
        self.scheme = getattr(self, scheme)

    @staticmethod
    def rk4(state, dt, tend):
        while True:
            k1 = dt*tend(state)
            k2 = dt*tend(state+k1/2)
            k3 = dt*tend(state+k2/2)
            k4 = dt*tend(state+k3)
            state = state + (k1+2*(k2+k3)+k4)/6
            yield state

    def iter_states(self, state:State):
        return self.scheme(state, self.dt, self.tend)

class BalanceModel(BalanceBase, OdeSolver):
    def __init__(self, dt:float, scheme:AnyStr='rk4'):
        BalanceBase.__init__(self)
        OdeSolver.__init__(self, tend=self.tend, dt=dt, scheme=scheme)
