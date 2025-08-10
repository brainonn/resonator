# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

class Resonator:
    def __init__(self, fs, Q, C0, fp=None, m=None):
        """
        Parameters
        ----------
        fs : float, Hz
            Series resonance frequency.
        Q : float
            Quality factor.
        C0 : float, F
            Static capacitance.
        fp : float, Hz
            Parallel resonance frequency. If None will be calcultaed from
            capacitance ratio `m`. The default is None.
        m : float, optional
            Capacitance ratio. If None will be calculated from resonance
            frequencies. The default is None.

        Returns
        -------
        None.

        """
        if fp is None and m is None:
            raise ValueError('You should specify either paraller resonance'
                             'frequency or capacitance ratio')
        self.fs = fs
        self.fp = fp if fp is not None else fs * np.sqrt(1 + m)
        self.m = m if m is not None else (fp ** 2 - fs ** 2) / fs ** 2
        self.Q = Q
        self.C0 = C0
    
    @classmethod
    def fromrlc(cls, R, L, C, C0):
        """
        Initiation form RLC equivalent parameters

        Parameters
        ----------
        R : float, Ohm
            Equivalent resistance.
        L : float, H
            Equivalent inductance.
        C : float, F
            Dynamic capacitance.
        C0 : float, F
            Static capacitance.

        Returns
        -------
        Resonator
            Class member with parameters derived from RLC.

        """
        fs = 1 / (2 * np.pi * np.sqrt(L * C))
        fp = 1 / (2 * np.pi * np.sqrt(L * C * C0 / (C + C0)))
        Q = 1 / R * np.sqrt(L / C)
        return cls(fs, Q, C0, fp=fp)
    
    def phase(self, f):
        """
        

        Parameters
        ----------
        f : float, Hz
            Frequency to get phase at.

        Returns
        -------
        float
            phase at frequency.

        """
        return np.arctan(self.Q * self.fs / (self.m * f) 
                         * ((f / (self.Q * self.fs)) ** 2 
                            + (f ** 2 / self.fs ** 2 - 1) 
                            * (f ** 2 / self.fs ** 2 
                               - self.fp ** 2 / self.fs ** 2)))
    
    def r(self, f):
        """

        Parameters
        ----------
        f : float, Hz
            Frequency to get real part of impedance at.

        Returns
        -------
        float, Ohm
            Equivivalent resistance.

        """
        return self.m / (self.fs * self.Q * self.C0) \
            * (1 / ((f / self.Q / self.fs) ** 2 
               + (self.fp ** 2 / self.fs ** 2 - f ** 2 / self.fs ** 2) ** 2))
    
    def x(self, f):
        """
        

        Parameters
        ----------
        f : float, Hz
            Frequency to get imaginary part of impedance at.

        Returns
        -------
        float, Ohm
            Reactance.

        """
        return - 1 / (f * self.C0) * (((f / self.Q / self.fs) ** 2 
                                       + (f ** 2 / self.fs ** 2 - 1) 
                                       * (f ** 2 / self.fs ** 2 
                                          - self.fp ** 2 / self.fs ** 2)) 
                                      / ((f / self.Q / self.fs) ** 2 
                                      + (self.fp ** 2 / self.fs ** 2 
                                         - f ** 2 / self.fs ** 2) ** 2))
    def impedance(self, f):
        """
        

        Parameters
        ----------
        f : float, Hz
            Frequency to get impedance at.

        Returns
        -------
        float, Ohm
            Impedance.

        """
        return self.r(f) + 1j * self.x(f)
    
    def g(self, f):

        """
        

        Parameters
        ----------
        f : float, Hz
            Frequency to get impedance at.

        Returns
        -------
        float, siemens
            Conductance.

        """
        return self.r(f) / np.abs(self.impedance(f))**2
    
    def b(self, f):
        """

        Parameters
        ----------
        f : float, Hz
            Frequency to get imaginary part of impedance at.

        Returns
        -------
        float, Ohm
            Susceptance.


        """
        return - self.x(f) / np.abs(self.impedance(f))**2
    

    def admittance(self, f):
        """
        

        Parameters
        ----------
        f : float, Hz
            Frequency to get impedance at.

        Returns
        -------
        float, Ohm
            Admittance.

        """
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
    C1 = 160e-12 #at amplifier output
    C2 = 36e-12 #at amplifier input
    Ro = 1e3
    g = 25e-3
    r = Resonator(fs, Q, C0, fp=fp)
    freqs = np.arange(0.999 * fs, 1.001 * fp)
    
    x1 = lambda f: x(C1, f)
    x2 = lambda f: x(C2, f)

    
    G = lambda f: -g * Ro * A(f, C1, C2, g=-1) / (Ro + A(f, C1, C2, g=-1))
    beta = lambda f: fbd(f, C1, C2)
    d = lambda f: r.x(f) + (1 + r.r(f) / Ro) * x1(f) + x2(f)
    
    
    fr2 = root_scalar(d, x0=fs, x1=fp)
    
    if fr2.converged:
        print(fr2.root)
    
    Gb = lambda f: g * x1(f) * x2(f)\
        / (r.r(f) - x1(f) / Ro * (r.x(f) + x2(f)) + 
           1j * (r.x(f) + x1(f) + x2(f) + x1(f) * r.r(f) / Ro))
    
    plt.plot(freqs, np.angle(Gb(freqs)))

    plt.axvline(fs, color='green')
    plt.axvline(fr2.root, color='red')
    plt.axvline(fp, color='black')
