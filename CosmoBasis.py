from math import *
import scipy.integrate as integrate
import csv
import sys
sys.path.append('camb4py-master/build/lib.linux-x86_64-2.7/camb4py/')
import camb4py
from fortranfile import FortranFile
sys.path.append('Bessel/Bessels_from_CAMB/')
import bessels # This is the fortran version
import bessel # This is the python version

class CosmoBasis(object):
    # Initializer fixes constant parameters
    def __init__(self, params):
        
        # store the initial parameters as the fiducial parameters.
        self.fiducial_params = params
        # check whether enough parameters were given
        # adds fiducial values for those not given
        self.check_params(params)
        
        #### Physical Constants ####
        # Fixing c [m/s]
        self.c = 299792458.0
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
        
        # generates useful parameters from input
        self.generate_params(self.fiducial_params)
        
        # CMB temperature [K]
        self.T_CMB = self.fiducial_params["T_CMB"]
        # T_gamma is usually set to T_CMB [K]
        self.T_gamma = self.T_CMB 
        # Spontaneous decay-rate of the spin-flip transition [s^-1]
        self.A_10 = 2.85 * 10**(-15)
        # Logarithmic tilt of the the spectrum of fluctuations
        self.n_s = 0.95
        # Variance of matter fluctuations today smoothed on a scale of 8 h^-1 Mpc
        self.sigma_8 = 0.8
           
        # Reionization regime
        self.delta_z_rei = 4.0
        self.z_rei = 10.0
        self.z_CMB = 1100

        #TODO: get a good bias factor
        # beta factor
        self.beta = 0.7
        # bias factor
        self.b_bias = self.O_M**0.6 / self.beta
        # k_eq scales at horizon crossing [Mpc^-1]
        self.k_eq = 0.073 * self.O_M * self.h**2
        
        #create CAMB object and results.
        self.camb = camb4py.load('CAMB/camb', defaults='CAMB/params.ini')
        
        print "CosmoBasis initialized"
#######################################################################

    # This function takes the list of cosmological parameters and translates them to parameters
    # which will be used directly in Cosmological computations.
    def generate_params(self, params):
        # Hubble constant H_0 [km s^-1 Mpc^-1]
        self.H_0 = params["hubble"]
        # Hubble parameter h [ dimensionless ]
        self.h = self.H_0 / 100.0
        # Baryon density
        self.O_b = params["ombh2"] / self.h**2
        # CDM density
        self.O_cdm = params["omch2"] / self.h**2
        # Neutrino density
        self.O_nu = params["omnuh2"] / self.h**2
        # TODO: IMPORTANT, find out if we should ignore O_R
        #       or how to include it
        # Fixing the relative radiation density in all cases
        self.O_R = 4.165E-1/self.H_0**2
        # Curvature term
        self.O_k = params["omk"] + self.O_R
        # Relative Matter density
        self.O_M = self.O_b + self.O_cdm + self.O_nu
        # total Omega
        self.O_tot = 1.0 - self.O_k
        # Relative Vacuum density
        self.O_V = self.O_tot - self.O_M
        # Hubble distance [Mpc]
        self.D_H = self.c / (1000.0 * self.H_0)
        # Hubble time
        self.t_H = 1.0 / self.H_0
        
        return None
    
    # This function checks if all necessary parameters have been added to the dictionary.
    # If not, it adds values.
    def check_params(self, params):
        if "ombh2" not in params.keys():
            self.fiducial_params["ombh2"] = 0.0226
            print "ombh2 not defined! Fiducial value assumed"
    
        if "omch2" not in params.keys():
            self.fiducial_params["omch2"] = 0.112
            print "omch2 not defined! Fiducial value assumed"

        if "omnuh2" not in params.keys():
            self.fiducial_params["omnuh2"] = 0.00064
            print "omnuh2 not defined! Fiducial value assumed"

        if "omk" not in params.keys():
            self.fiducial_params["omk"] = 0.0
            print "omk not defined! Fiducial value assumed"

        if "hubble" not in params.keys():
            self.fiducial_params["hubble"] = 70.0
            print "hubble not defined! Fiducial value assumed"

        if "T_CMB" not in params.keys():
            self.fiducial_params["T_CMB"] = 2.7255
            print "T_CMB not defined! Fiducial value assumed"

        if "zmin" not in params.keys():
            self.fiducial_params["zmin"] = 7.0
            print "zmin not defined! Fiducial value assumed"

        if "zmax" not in params.keys():
            self.fiducial_params["zmax"] = 9.0
            print "zmax not defined! Fiducial value assumed"

        if "zsteps" not in params.keys():
            self.fiducial_params["zsteps"] = 100
            print "zsteps not defined! Fiducial value assumed"

        if "Pk_steps" not in params.keys():
            self.fiducial_params["Pk_steps"] = 3
            print "Pk_steps not defined! Fiducial value assumed"


        return None

    # This function uses the subroutine bjl from CAMB to efficiently compute spherical bessel functions. 
    # If for some reason f2py is not working, then comment first line here and uncomment second.
    # This will use a Python translation of the Fortran code, but this is slower!
    def sphbess_camb(self, l, x):
        #return bessels.bjl(l, x)
        return bessel.bjl(l, x)

    # Helper function, denominator for various integrals
    # E(z) = H(z)/H_0 
    def E(self, z):
        return sqrt(self.O_V + self.O_R * (1+z)**4 +
                self.O_M * (1+z)**3 + self.O_k * (1+z)**2)

    # Helper function, Z = int (dz/E, 0, z)
    def Z(self, z):
        integral = integrate.quad(lambda x: 1.0/self.E(x), 0, z)
        return integral[0]
    
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

    def mpc_to_m(self, x):
        conv_factor = 3.0856776 * 10**22
        return conv_factor * x
    def m_to_mpc(self, x):
        conv_factor = 3.0856776 * 10**22
        return x / conv_factor 


