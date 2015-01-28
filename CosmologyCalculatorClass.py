######################################################################
# This file contains a calculator for basic cosmological calculations #
# similar to Ned Wrigths' calulator:                                  # 
# http://www.astro.ucla.edu/~wright/CosmoCalc.html                    #
# or http://arxiv.org/pdf/astro-ph/0609593v2.pdf                      # 
# aka A Cosmology Calculator for the World Wide Web                   #
# See also, http://arxiv.org/pdf/astro-ph/9905116v4.pdf               #
# for a review of the formulas used.                                  #
# Or Distance measures in Cosmology, David W. Hogg                    #
#######################################################################

from math import *
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.special as special
from sympy import Symbol
import numpy as np

class CosmoCalc(object):
    # Initializer fixes constant parameters
    def __init__(self, H_0, O_M, O_V, T_CMB):
        # Hubble constant H_0 [km s^-1 Mpc^-1]
        self.H_0 = H_0
        # Hubble parameter h [ dimensionless ]
        self.h = self.H_0 / 100
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
#######################################################################
    
    # Helper function, denominator for various integrals
    # E(z) = H(z)/H_0 
    def E(self, z):
        return sqrt(self.O_V + self.O_R * (1+z)**4 +
                self.O_M * (1+z)**3 + (1-self.O_tot) * (1+z)**2)

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
    
    # Hubble Time [s * Mpc/km]: 
    #   t_H = 1/H_0
    def hubble_time(self):
        return self.t_H

    # Hubble Distance [Mpc]: 
    #   D_H = c/H_0
    def hubble_dist(self):
        return self.D_H

    # Comoving distance (line of sight) [Mpc], 
    #   aka. Comoving radial distance (D_now):
    #   D_C = D_H int(dz / E(z),0,z) = D_H * Z = cZ/H_0
    def comoving_radial_dist(self, z):
        return self.hubble_dist() * self.Z(z)
    def D_C(self, z):
        return self.comoving_radial_dist(z)
    def D_now(self, z):
        return self.comoving_radial_dist(z)


    # Comoving distance (transverse) [Mpc],
    #   aka. Proper motion distance:
    #   D_M = D_H/sqrt(1-O_tot) S_k(sqrt(1-O_tot) D_C/D_H)
    def comoving_dist_transverse(self, z):
        O_k = 1-self.O_M-self.O_V
        return self.D_H / sqrt(abs(1-O_k)) * \
                self.S_k(sqrt(abs(1-O_k)) * \
                self.D_C(z)/self.D_H)
    def D_M(self, z):
        return self.comoving_dist_transverse(z)

    # Angular diameter distance [Mpc]:
    #   D_A = D_M / (1+z)
    # 
    # Second Form:
    # Angular diam dist between 2 objects at redshifts z & z2 [Mpc]
    #   formula only holds for O_tot <= 1:
    #   D_A12 = 1/(1+z_2) * (D_M2 sqrt(1+(1-O_tot) D_M1^2/D_H^2) -
    #                       D_M1 sqrt(1+(1-O_tot) D_M2^2/D_H^2) )
    def angular_diam_dist(self, z, z2 = None):
        O_k = 1-self.O_M-self.O_V
        root = sqrt(abs(1-self.O_tot))
        
        if z2 is None:
            return self.D_H * self.S_k((root) * self.Z(z))/ \
                    ((1+z) * root)
        elif (self.O_V + self.O_M <= 1.0):
            return 1.0 / (1.0 + z2) * \
                    ( self.D_M(z2) * sqrt(1 + (1 - O_k) *\
                    self.D_M(z)**2 / self.D_H**2) - \
                    self.D_M(z) * sqrt(1 + (1 - O_k) *\
                    self.D_M(z2)**2 / self.D_H**2) )
        else:
            print "Error: D_A12 formula invalid for O_tot > 1.0"
            return -1
    def D_A(self, z, z2 = None):
        return self.angular_diam_dist(z, z2)

    # Luminosity distance [Mpc]:
    #   D_L = (1 + z)^2 * D_A
    def luminosity_dist(self, z):
        return self.D_A(z) * (1+z)**2
    def D_L(self, z):
        return self.luminosity_dist(z)

    # Comoving Volume [Mpc^3]:
    #   V_C = int D_H (1+z)^2 * D_A^2 / E(z) dOmega dz
    def comoving_volume(self, z):
        vol = 4*pi*self.D_H
        integral = integrate.quad(lambda x: (1+x)**2 * \
                self.D_A(x)**2/self.E(x), 0, z)
        return vol*integral[0]
    def V_C(self, z):
        return self.comoving_volume(z)
    
    # Age of the Universe at redshift z [s * Mpc/km]:
    #   t(z) = t_H int(dz/((1+z) E(z)), z, inf)
    #def age_of_universe(self, z):
    #    age = integrate.quad(lambda a: 1.0/self.D(a), 0.0, 1.0/(1.0+z))
    #   return self.convert_to_gy(age[0]/H_0)
    def age_of_universe(self, z):
        age = integrate.quad(lambda x: 1.0/((1+x)*self.E(x)) \
                , z, np.inf)
        return age[0] * self.t_H

    # Light travel time [s * Mpc/km]:
    #   ltt = t(0) - t(z)
    def light_travel_time(self, z):
        return self.age_of_universe(0) - self.age_of_universe(z)

    # Distance based on light travel time [Mpc]
    # D_ltt = c * (t(0) - t(z))
    def D_ltt(self, z):
        return self.light_travel_time(z) * self.c / 1000.0
    
    # Hubble constant with redshift [km * s^-1 * Mpc^-1]
    def H(self, z):
        return self.H_0 * self.E(z)
    # Hubble constant in SI units [s^-1]
    def H_SI(self, z):
        return self.H(z) * 1000 / (3.08567758 * 10E16)
    # critical density [kg/m^3]
    def rho_crit(self, z):
        return 3.0 * self.H_SI(z)**2 / (8.0 * pi * self.G)
    # matter density [kg/m^3]
    def rho_M(self, z):
        rho_M_0 = self.O_M * self.rho_crit(0)
        return rho_M_0 * (1+z)**3
    # Radiation density [kg/m^3]
    def rho_R(self, z):
        rho_R_0 = self.O_R * self.rho_crit(0)
        return rho_R_0 * (1+z)**4
    # Vacuum density [kg/m^3] (constant in z)
    def rho_V(self, z):
        return self.O_V * self.rho_crit(0)
    # Relative Matter density with redshift
    def Omega_M(self, z):
        #return self.rho_M(z) / self.rho_crit(z)
        return self.O_M * (1+z)**3 / self.E(z)**2
    # Relative Radiation density with redshift
    def Omega_R(self, z):
        #return self.rho_R(z) / self.rho_crit(z)
        return self.O_R * (1+z)**4 / self.E(z)**2
    # Relative Vacuum density with redshift
    def Omega_V(self, z):
        #return self.rho_V(z) / self.rho_crit(z)
        return self.O_V / self.E(z)**2

    # Estimate for the number of Baryons in the Universe
    def num_baryons(self):
        num = 4.0/3.0 * pi * (self.c/self.H_SI(0))**3 * self.rho_crit(0)
        return num / self.m_b
    
    # Neutral fraction x_HI
    def x_HI(self, z, delta_z, z_reionization):
        return 0.5 * (tanh((z-z_reionization) / delta_z) + 1)

    ##################### Thermodynamics  ########################
    
    # Temperature of the universe in terms of redshift z [K]
    def T(self, z):
        return self.T_CMB * (1+z)
    # Total hydrogen density, n_H + n_p [m^-3]
    def n_H_tot (self, z):
        return 1.6 * (1+z)**3

    # Number densities [m^-3]
    # baryon density
    def n_b(self, z):
        return self.O_b * self.rho_crit(z) / self.m_b
    # hydrogen density
    # TODO: this has a coefficient that is only approximately 1.
    def n_H(self, z):
        coeff = 1.0 
        return coeff * self.n_b(z)
    # (free) proton density
    def n_p(self, z):
        return self.x_HI(z, self.delta_z, self.z_rei)
    # (free) electron density
    def n_e(self, z):
        return self.n_p(z)

    ####################################################################
    ####### Two point Correlation of brightness tem. fluctuations ######
    ####################################################################

    ############### first we ignore redshift space distortions #########
   
    # Calculates the Tb correlation without redshift space distortions 
    # between redshifts z1 & z2.
    # use eg. z1 = 6 & z2 = 10 
    def corr_Tb_no_distortions(self, l, k1, k2, z_low, z_high):
        #TODO: figure out integration limits kkk
        integral = integrate.quad(lambda k: k**2 * self.P_delta(k) *\
                    self.M(l, k1, k, z_low, z_high) *\
                    self.M(l, k2, k, z_low, z_high), 0.001, 10)
        return integral[0]
    
    # M_l(k,k') = 2b/pi * int dr r^2 deltaT_b(r) j_l(kr) j_l(k'r) []
    def M(self, l, k1, k2, z_low, z_high):
        r_low = self.comoving_radial_dist(z_low)
        r_high = self.comoving_radial_dist(z_high)

        integral = integrate.quad(lambda r: r**2 * self.delta_Tb_bar(r) *\
                    special.sph_jn(l, k1*r)[0][l-1] *\
                    special.sph_jn(l, k2*r)[0][l-1] *\
                    self.P_growth(r), r_low, r_high)
        return 2*self.b_bias/pi*integral[0]

    # Mean Brightness Temperature fluctuations at distance r (comoving) [K]
    def delta_Tb_bar(self, r):
        #TODO: find sensible constant
        constant = 1
        return constant

    # we seperate the power spectrum P(k,a) into a growing mode and P(k)
    # P(k,a) = P_delta(k) * P_growth(a)
    def P_growth(self, r):
        #convert r to a
        # here we use an approximation for the universe to be an Einstein de Sitter
        # universe. This makes inverting r(z) easier and analytically doable. 
        # In general, a numerical inversion is necessary!
        # The factor of 1000 is to balance units.
        z = (2*self.c / (2*self.c - 1000 * self.H_0 * r)) - 1

        res = self.D1(z) / self.D1(0)
        return res**2

    def D1(self, z):
        # need to relate a to z
        # have done this
        prefactor = 5 * self.O_M / 2 * self.E(z) 
        # x = a' 
        integral = integrate.quad(lambda x: (1+x) / self.E(x)**3, z, np.inf)
        
        return prefactor * integral[0]
        
    # This is the part of the power spectrum that only depends on the scale k
    # Units: []
    def P_delta(self, k):
        amplitude = 1
        x = k / self.k_eq
        transfer_function = self.transfer(x)**2

        P = amplitude * transfer_function 
        return P 

    # BBKS transfer function, Dodelson (7.70)
    def transfer(self, x):
        res = np.log(1 + 0.171 * x) / (0.171 * x)
        bracket = 1+ 0.284 * x + (1.18 * x)**2 + (0.399 * x)**3 + (0.490 * x)**4
        res = res * bracket**(-0.25)
        return res
