
from math import *
import scipy.integrate as integrate
import numpy as np
import scipy.optimize as op
import mpmath as mp
import scipy.special as special
import csv

class CosmoBasis(object):
    # Initializer fixes constant parameters
    def __init__(self, H_0, O_M, O_V, T_CMB, besseltable):
        # Hubble constant H_0 [km s^-1 Mpc^-1]
        self.H_0 = H_0
        # Hubble parameter h [ dimensionless ]
        self.h = self.H_0 / 100.0
        # Relative Matter density
        self.O_M = O_M
        # Relative Vacuum density
        self.O_V = O_V
        # Fixing the relative radiation density in all cases
        self.O_R = 4.165E-1/self.H_0**2
        # Fixing c [m/s]
        self.c = 299792458.0
        # Total density
        self.O_tot = self.O_M + self.O_V + self.O_R
        # Hubble distance [Mpc]
        self.D_H = self.c / (1000.0 * H_0)
        # Hubble time
        self.t_H = 1.0 / H_0
        # Boltzmann Constant [m^2 kg s^-2 K^-1]
        self.k_b = 1.3806488 * 10**(-23)
        # Baryon mass (mass of neutron) [kg]
        self.m_b = 1.674927351 * 10**(-27)
        # Electron mass [kg]
        self.m_e = 9.10938291 * 10**(-31)
        # Charge of electron [Coulomb]
        self.e = 1.60217657 * 10**(-19)

        # Planck Constant h [m^2 kg s^-1]
        self.h_planck = 6.62606957 * 10**(-34)
        # Gravitational constant [m^3 kg^-1 s^2]
        self.G = 6.67384 * 10**(-11)
        # T_* = hc/k_b Lambda_21cm [K]
        self.T_star = self.h_planck * self.c / (self.k_b * 0.21)
        # CMB temperature [K]
        self.T_CMB = T_CMB
        # T_gamma is usually set to T_CMB [K]
        self.T_gamma = self.T_CMB 
        # Spontaneous decay-rate of the spin-flip transition [s^-1]
        self.A_10 = 2.85 * 10**(-15)
        # Baryon density
        self.O_b = 0.044
        # Logarithmic tilt of the the spectrum of fluctuations
        self.n_s = 0.95
        # Variance of matter fluctuations today smoothed on a scale of 8 h^-1 Mpc
        self.sigma_8 = 0.8
           
        # Reionization regime
        self.delta_z = 4.0
        self.z_rei = 10.0
        self.z_CMB = 1100

        #TODO: get a good bias factor
        # bias factor
        self.b_bias = 1
        # k_eq scalles at horizon crossing [Mpc^-1]
        self.k_eq = 0.073 * self.O_M * self.h**2
        
        #this includes a matrix of all the values in the bessel function table
        self.besseltable = self.read_bessel_table(besseltable)

#######################################################################
    # This function reads in a table of bessel function values and
    # retruns it as a matrix.
    def read_bessel_table(self, bt):

        file = open(bt, 'r')
        data = csv.reader(file, delimiter=' ')
        bessel_table = [row for row in data ]
        
        self.bessel_lmin = int(bessel_table[0][0])
        self.bessel_lmax = int(bessel_table[0][1])
        self.bessel_xmin = float(bessel_table[0][2])
        self.bessel_xmax = float(bessel_table[0][3]) 
        self.bessel_xstep = float(bessel_table[0][4]) 
           
        bessel_table.pop(0)
        return bessel_table

    # Interpolation function for Spherical bessel functions.
    # This works as intended and is quite fast. Requires a precalculated 
    # grid of bessel values. Use code in Bessel/ to generate this list.
    def sphbess_interp(self, l, x):
        n = (x-self.bessel_xmin)/self.bessel_xstep
        y0 = float(self.besseltable[l][int(n)])
        y1 = float(self.besseltable[l][int(n) + 1])

        res = (y1-y0)*(n-int(n)) + y0
        return res
    
    # Helper function to compute spherical bessel functions
    def sphbess(self, l, x):
        if x > 10*l*l:
            bess = 1/x * sin(x - l*pi/2)
        else:
            bess = special.sph_jn(l,x)[0][l]
            if isnan(bess):
                bess = 1/x * sin(x - l*pi/2)
        return bess
    def new_sphbess(self, l, x):
        bess = special.sph_jn(l,x)[0][l]
        if isnan(bess):
            bess = 1/x * mp.sin(x - l*pi/2)
        return bess

    # Helper function to compute spherical bessel functions fast
    def fsphbess(self, n, x):
        return sqrt(pi/(2*x)) * mp.fp.besselj(n+0.5, x)

    # Helper function to compute spherical bessel functions fast
    def mp_sphbess(self, n, x):
        return sqrt(pi/(2*x)) * mp.fp.besselj(n+0.5, x)

    # Helper function, denominator for various integrals
    # E(z) = H(z)/H_0 
    def E(self, z):
        return sqrt(self.O_V + self.O_R * (1+z)**4 +
                self.O_M * (1+z)**3 + (1-self.O_tot) * (1+z)**2)

    # Helper function, Z = int (dz/E, 0, z)
    def Z(self, z):
        integral = integrate.quad(lambda x: 1.0/self.E(x), 0, z)
        return integral[0]
    def Z_opt(self, z):
        integral = mp.fp.quad(lambda x: 1.0/self.E(x), [0, z])
        return integral

    # Curvature term in the metric:
    #   sinh x, x or sin x 
    def S_k(self, x):
        if (self.O_V + self.O_M < 1.0):
            return sinh(x) 
        elif (self.O_V + self.O_M == 1.0):
            return x
        else:
            return sin(x)

    # Simple Simpson integration method. Requires even number of steps.
    def integrate_simps(self, f, a, b, steps):
        if steps % 2:
            raise ValueError("There must be an even number of steps, received n = %d" % steps)
        stepsize = float(b-a)/float(steps)
        res = f(a) + f(b)

        for n in range(1, steps, 2):
            res += 4 * f(a + n * stepsize)
        for n in range(2, steps - 1, 2):
            res += 2 * f(a + n * stepsize)

        return res * stepsize / 3.0
