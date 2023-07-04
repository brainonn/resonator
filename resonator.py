# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from matplotlib import pyplot as plt
import numpy as np

class Resonator:
    def __init__(self, fs, fp, Q, C0, m=None):
        self.fs = fs
        self.fp = fp
        self.m = m if m is not None else (fp ** 2 - fs ** 2) / fs ** 2
        self.Q = Q
        self.C0 = C0
    
    @classmethod
    def fromrlc(cls, R, L, C, C0):
        fs = 1 / (2 * np.pi * np.sqrt(L * C))
        fp = 1 / (2 * np.pi * np.sqrt(L * C * C0 / (C + C0)))
        Q = 1 / R * np.sqrt(L / C)
        return cls(fs, fp, Q, C0)
    
    def phase(self, f):
        return np.arctan(self.Q * self.fs / (self.m * f) 
                         * ((f / (self.Q * self.fs)) ** 2 
                            + (f ** 2 / self.fs ** 2 - 1) 
                            * (f ** 2 / self.fs ** 2 
                               - self.fp ** 2 / self.fs ** 2)))
    
    def r(self, f):
        return self.m / (self.fs * self.Q * self.C0) \
            * (1 / ((f / self.Q / self.fs) ** 2 
               + (self.fp ** 2 / self.fs ** 2 - f ** 2 / self.fs ** 2) ** 2))
    
    def x(self, f):
        return - 1 / (f * self.C0) * (((f / self.Q / self.fs) ** 2 
                                       + (f ** 2 / self.fs ** 2 - 1) 
                                       * (f ** 2 / self.fs ** 2 
                                          - self.fp ** 2 / self.fs ** 2)) 
                                      / ((f / self.Q / self.fs) ** 2 
                                      + (self.fp ** 2 / self.fs ** 2 
                                         - f ** 2 / self.fs ** 2) ** 2))
    def impedance(self, f):
        return self.r(f) + 1j * self.x(f)
    
    def g(self, f):
        return self.r(f) / np.abs(self.impedance(f)) ** 2
    def b(self, f):
        return - self.x(f) / np.abs(self.impedance(f)) ** 2
    def admittance(self, f):
        return 1 / self.impedance(f)

def x(C, f):
    return -1 / (2 * np.pi * f * C)

def fbd(f, C1, C2):
    return 1j * x(C2, f) / (r.r(f) + 1j * (r.x(f) + x(C2, f)))

def A(f, C1, C2, g=1):
    return -g * 1j * x(C1, f) * (r.r(f) + 1j * (r.x(f) + x(C2, f))) \
        / (r.r(f) + 1j * (r.x(f) + x(C1, f) + x(C2, f)))

if __name__ == '__main__':
    fs = 1055561
    fp = 1055721
    Q = 16251.907728114395
    C0 = 7.60393485503531e-12
    C1 = 36e-12
    C2 = 36e-12
    r = Resonator(fs, fp, Q, C0)
    freqs = np.arange(0.999 * fs, 1.001 * fp) 
    plt.plot(freqs, r.x(freqs))
    plt.axvline(fs, color='red')
    plt.axvline(fp, color='black')
