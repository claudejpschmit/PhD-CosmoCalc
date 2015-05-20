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
import copy

class CosmoBasis(object):
    # Initializer fixes constant parameters
    def __init__(self, params):
        
        # store the initial parameters as the fiducial parameters.
        self.fiducial_params = copy.deepcopy(params)
        # check whether enough parameters were given
        # adds fiducial values for those not given
        self.check_params(params)
        self.current_params = copy.deepcopy(self.fiducial_params)

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

    # Note on neutrinos:
    # Here we add the neutrino density to the matter density as well as to the radiation density. 
    # For high redshift, neutrinos are relativistic and add to the radiation density, also at 
    # high redshifts, radiation is dominating so having the neutrinos in has a great effect,
    # since O_nu is of the same order as O_gamma.
    # At late times O_R becomes very unimportant so adding O_nu in has a very minor effect.
    # As for adding O_nu to O_M, O_b/O_nu = 100 and O_CDM/O_nu = 1000 so adding it in or not 
    # should not have a great effect at any epoch.

    def generate_params(self, params):
        # CMB temperature [K]
        p = copy.deepcopy(params)
        self.T_CMB = p["T_CMB"]
        # T_gamma is usually set to T_CMB [K]
        self.T_gamma = self.T_CMB
        # Hubble constant H_0 [km s^-1 Mpc^-1]
        self.H_0 = p["hubble"]
        # Hubble parameter h [ dimensionless ]
        self.h = self.H_0 / 100.0
        # Baryon density
        self.O_b = p["ombh2"] / self.h**2
        # CDM density
        self.O_cdm = p["omch2"] / self.h**2
        # Neutrino density - non-relativistic
        self.O_nu = p["omnuh2"] / self.h**2
        # Relatve Photon density
        self.O_gamma =  pi**2 * (self.T_CMB/11605.0)**4 / (15 * 8.098 * 10**(-11) * self.h**2)
        # Relativistic Neutrinos: rho_nu = 3* 7/8 * (4/11)^(4/3) * rho_gamma
        self.O_nu_rel = self.O_gamma *3.0*7.0/8.0*(4.0/11.0)**(4.0/3.0)
        # Fixing the relative radiation density in all cases
        self.O_R = self.O_gamma + self.O_nu_rel
        # Curvature term
        self.O_k = p["omk"]
        # Relative Matter density
        self.O_M = self.O_b + self.O_cdm + self.O_nu 
        # total Omega
        self.O_tot = 1 - self.O_k
        # Relative Vacuum density
        self.O_V = self.O_tot - self.O_M - self.O_R
        # Hubble distance [Mpc]
        self.D_H = self.c / (1000.0 * self.H_0)
        # Hubble time
        self.t_H = 1.0 / self.H_0

        self.current_params = p
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

        # TODO: Changed this from 100 to 1000
        if "zsteps" not in params.keys():
            self.fiducial_params["zsteps"] = 100
            print "zsteps not defined! Fiducial value assumed"

        if "Pk_steps" not in params.keys():
            self.fiducial_params["Pk_steps"] = 3
            print "Pk_steps not defined! Fiducial value assumed"

        if "k_steps" not in params.keys():
            self.fiducial_params["k_steps"] = 40000
            print "zsteps not defined! Fiducial value assumed"

        return None

    # This function uses the subroutine bjl from CAMB to efficiently compute spherical bessel functions. 
    # If for some reason f2py is not working, then comment first line here and uncomment second.
    # This will use a Python translation of the Fortran code, but this is slower!
    def sphbess_camb(self, l, x):
        # The fortran version
        #return bessels.bjl(l, x)
        # The python version
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
        if (self.O_tot < 1.0):
            return sinh(x) 
        elif (self.O_tot == 1.0):
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


    def params_to_planck15(self, params):
        params["T_CMB"] = 2.7255
        params["ombh2"] = 0.02237
        params["omch2"] = 0.1187
        params["hubble"] = 67.74
        params["tau"] = 0.068
        params["tilt"] = 0.9678
        params["omk"] = 0.0
        params["omnuh2"] = 0.0025
        params["sigma_8"] = 0.8159


