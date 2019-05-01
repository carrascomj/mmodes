#!/usr/bin/python3

# Fixed-step ODE solver classes. Although scipy.ode provide variable-step more
# reliable algorithms, sometimes it's not feashible to use them.

# Author: Jorge Carrasco Muriel
# e-mail: jorge.cmuriel@alumnos.upm.es
# Date of first version: 15/03/2019

import numpy as np

class ODESolver():
    ''' Major parent ODE solver class '''
    def __init__(self, f, mod, dt = 0.01):
        if not callable(f):
            raise NotIntegratorError('f is %s, not a function' % type(f))
        self.f = lambda u, t: np.asarray(f(u, t, mod))
        self.dt = dt

    def set_initial_value(self, y, t):
        ''' Sets initial condition of the system, as in scipy.ode '''
        self.y = np.asarray(y)
        self.t = t

class FEA(ODESolver):
    ''' Forward Euler Approach '''
    def integrate(self, maxT, step):
        dt, y, t, f, = self.dt, self.y, self.t, self.f
        self.y = y + dt*f(y, t)
        self.t = t + dt

class RungeKutta4(ODESolver):
    ''' 4th-order Runge Kutta method '''
    def integrate(self, maxT, step):
        dt, y, t, f, = self.dt, self.y, self.t, self.f
        dt2 = dt/2.0
        K1 = dt*f(y, t)
        K2 = dt*f(y + 0.5*K1, t + dt2)
        K3 = dt*f(y + 0.5*K2, t + dt2)
        K4 = dt*f(y + K3, t + dt)
        self.y = y + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        self.t = t + dt
