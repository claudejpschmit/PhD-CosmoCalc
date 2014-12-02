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
from sympy import Symbol
import numpy as np
import matplotlib.pyplot as plt

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

        # Ratio of spin degeneracy factors of the 
        # 1S triplet and singlet state
        self.ratio_g1g0 = 3
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
        # T_
        # Spontaneous decay-rate of the spin-flip transition [s^-1]
        self.A_10 = 2.85 * 10**(-15)
        # Baryon density
        self.O_b = 0.044
        # Logarithmic tilt of the the spectrum of fluctuations
        self.n_s = 0.95
        # Variance of matter fluctuations today smoothed on a scale of 8 h^-1 Mpc
        self.sigma_8 = 0.8
        
        # Collisional Coupling scattering rates [cm^3 s^-1]
        # between H and H 
        kappa_HH = {'1' : 1.38 * 1E-13,
                 '2' : 1.43 * 1E-13, 
                 '5' : 4.65 * 1E-13,
                 '10' : 2.88 * 1E-12, 
                 '20' : 1.78 * 1E-11,
                 '50' : 6.86 * 1E-11, 
                 '100' : 1.19 * 1E-10,
                 '200' : 1.75 * 1E-10, 
                 '500' : 2.66 * 1E-10,
                 '1000' : 3.33 * 1E-10, 
                 '2000' : 0,
                 '3000' : 0, 
                 '5000' : 0,
                 '7000' : 0, 
                 '10000' : 0,
                 '15000' : 0, 
                 '20000' : 0}  
        x = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        y = [kappa_HH[str(i)] for i in x]
        self.kappa_HH = interpolate.interp1d(x, y)

        # between H and p
        kappa_Hp = {'1' : 0.4028,
                         '2' : 0.4517, 
                         '5' : 0.4301,
                         '10' : 0.3699, 
                         '20' : 0.3172,
                         '50' : 0.3047, 
                         '100' : 0.3379,
                         '200' : 0.4043, 
                         '500' : 0.5471,
                         '1000' : 0.7051, 
                         '2000' : 0.9167,
                         '3000' : 1.070, 
                         '5000' : 1.301,
                         '7000' : 1.480, 
                         '10000' : 1.695,
                         '15000' : 1.975, 
                         '20000' : 2.201 }   
        x = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,
                3000, 5000, 7000, 10000, 15000, 20000]
        y = [kappa_Hp[str(i)] for i in x]
        self.kappa_Hp = interpolate.interp1d(x, y)


        # between H and e
        kappa_He = {'1' : 0.239,
                         '2' : 0.337, 
                         '5' : 0.503,
                         '10' : 0.746, 
                         '20' : 1.05,
                         '50' : 1.63, 
                         '100' : 2.26,
                         '200' : 3.11, 
                         '500' : 4.59,
                         '1000' : 5.92, 
                         '2000' : 7.15,
                         '3000' : 7.71, 
                         '5000' : 8.17,
                         '7000' : 8.32, 
                         '10000' : 8.37,
                         '15000' : 8.29, 
                         '20000' : 8.11 }  
        y = [kappa_He[str(i)] for i in x]
        self.kappa_He = interpolate.interp1d(x, y)

 
        # Einstein Coefficients
        self.B_10 = self.A_10 * 0.21**3 / (2 * self.h_planck * self.c)
        self.B_01 = self.ratio_g1g0 * self.B_10
               
        # Oscillator strength of the Lyalpha transition
        self.f_alpha = 0.4162
        # TODO: figure out S_alpha
        self.S_alpha = 1.0

        # Reionization regime
        self.delta_z = 4.0
        self.z_rei = 10.0
        self.z_CMB = 1100
#######################################################################
    
    # Helper function, denominator for various integrals
    # E(z) 
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

    ####################### 21cm Stuff  ##########################

    # Oth order brighness temperature [mK]
    def T_b_0th_Order(self, z):
        norm = 27 * self.O_b * self.h**2 / 0.023 * sqrt(0.15 / (10 * self.O_M * self.h**2))
        spinTemp = (self.T_S(z) - self.T(z))/self.T_S(z)
        return norm * self.x_HI(z, self.delta_z, self.z_rei) * (1+z)**0.5 * spinTemp 

    #### Spin temperature Calculation ####
    #### This is according to Chapter 9 of:
    #### ned.ipac.caltech.edu/level5/Sept06/Loeb/Loeb_contents.html

    def C_10(self, T_k):
        # TODO: this is a function of the gas temperature?!
        #       and what is n_H ? we have n_h prop (1+z)^3
        # ??? self.C_10 = 4.0 / 3.0 * self.kappa_H * n_H ? 
        pass
    
    def C_01(self, T_k):
        # TODO: ??? not defined
        pass
    # TODO: the color temperature of the surrounding bath of radio photons
    def T_alpha(self, z):
        return self.T_k(z)
    # TODO: the gas kinetic temperature
    def T_k(self, z):
        return self.T_CMB * (1 + self.z_CMB) * ((1+z)/(1+self.z_CMB))**2 + 0.5 * 10**2 * (tanh(z-self.z_CMB/self.z_CMB) + 1)


    def T_S (self, z):
        num = 1 + self.x_alpha(z) + self.x_c(z)
        den = 1.0/self.T_gamma + self.x_alpha(z)/self.T_alpha(z) + self.x_c(z)/self.T_k(z)
        return num / den

    ########
    # Total Collisional Coupling Coefficient
    def x_c(self, z):
        norm = self.T_star / (self.T_gamma * self.A_10)
        res = self.n_H(z) * self.kappa_HH(self.T_k(z))\
                + self.n_p(z) * self.kappa_Hp(self.T_k(z)) \
                + self.n_e(z) * self.kappa_He(self.T_k(z))
        return res * norm

    # Wouthuysen-Field effect coupling
    def x_alpha(self, z):
        num = 16 * pi**2 * self.T_star * self.e**2 * self.f_alpha * \
                self.S_alpha * self.J_alpha(z)
        den = 27 * self.A_10 * self.T_gamma * self.m_e * self.c
        return num / den

    def J_alpha(self, z):
        prefactor = 1.165 * 10**12 * (1+z)
        return  0.5 * prefactor * (tanh(z-self.z_CMB/self.z_CMB) + 1)

    ############################ Plotting ############################

    def plot_distances(self):
        normalization = self.H_0/self.c*1000
        x = [i/100.0 for i in range(0,1000)]
        y1 = [self.D_A(i) * normalization for i in x]
        y2 = [self.D_L(i) * normalization for i in x]
        y3 = [self.D_now(i) * normalization for i in x]
        y4 = [self.D_ltt(i) * normalization for i in x]

        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'$D_A$')
        plot2 = plt.loglog(x, y2, basex = 10, basey = 10, label = r'$D_L$')
        plot3 = plt.loglog(x, y3, basex = 10, basey = 10, label = r'$D_{now}$')
        plot4 = plt.loglog(x, y4, basex = 10, basey = 10, label = r'$D_{ltt}$')
        plt.legend(loc = 'upper left')
        plt.title(r'Cosmological Distances vs Redshift for $\Omega_{tot} =$ %s' % (self.O_M + self.O_V))
        plt.xlabel('Redshift z')
        plt.ylabel(r'$H_0 D / c$')
        plt.ylim(0.01)
        plt.grid(True)
        
        plt.savefig('distances.png')
        plt.show()

    def plot_densities_rho(self, zmax, step):
        maxrange = zmax * step
        x = [float(i)/step for i in range(0,maxrange)]
        y1 = [self.rho_M(z) for z in x]
        y2 = [self.rho_R(z) for z in x]
        y3 = [self.rho_V(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$\rho_M$')
        plot2 = plt.plot(x, y2, label = r'$\rho_R$')
        plot3 = plt.plot(x, y3, label = r'$\rho_V$')
        plt.legend(loc = 'upper left')
        plt.title(r'Densities vs Redshift' )
        plt.xlabel('Redshift z')
        plt.ylabel(r'$\rho$')
        plt.grid(True)

        plt.savefig('densities.png')
        plt.show()
    
    def plot_densities_Omega(self, zmax, step):
        maxrange = zmax * step
        x = [float(i)/step for i in range(0, maxrange)]
        y1 = [self.Omega_M(z) for z in x]
        y2 = [self.Omega_R(z) for z in x]
        y3 = [self.Omega_V(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$\Omega_M$')
        plot2 = plt.plot(x, y2, label = r'$\Omega_R$')
        plot3 = plt.plot(x, y3, label = r'$\Omega_V$')
        plt.legend(loc = 'upper left')
        plt.title(r'Relative Densities vs Redshift' )
        plt.xlabel('Redshift z')
        plt.ylabel(r'$\Omega$')
        plt.grid(True)
        plt.xscale('log')

        plt.savefig('relative_densities.png')
        plt.show()

    def plot_H(self, zmax, step):
        maxrange = zmax * step
        x = [float(i)/step for i in range(0, maxrange)]
        y1 = [self.H(z)/100.0 for z in x]
        plot1 = plt.plot(x, y1, label = r'$h$')
        plt.legend(loc = 'upper left')
        plt.title(r'Hubble parameter $h$ vs Redshift' )
        plt.xlabel('Redshift $z$')
        plt.ylabel(r'$h$')
        plt.grid(True)

        plt.savefig('hubble_parameter.png')
        plt.show()
        
    def plot_x_HI(self, zmax, step):
        maxrange = zmax * step
        x = [float(i)/step for i in range(0, maxrange)]
        y1 = [self.x_HI(z, 0.5, 10) for z in x]
        plot1 = plt.plot(x, y1, label = r'$x_{HI}$')
        plt.legend(loc = 'upper left')
        plt.title(r'Neutral fraction $x_{HI}$ vs Redshift' )
        plt.xlabel(r'Redshift $z$')
        plt.ylabel(r'$h$')
        plt.grid(True)

        plt.savefig('x_HI.png')
        plt.show()

    
    def plot_kappa_HH(self, T_k_min, T_k_max, step):
        points = {'1' : 1.38 * 10**(-13),
                 '2' : 1.43 * 10**(-13), 
                 '5' : 4.65 * 10**(-13),
                 '10' : 2.88 * 10**(-12), 
                 '20' : 1.78 * 10**(-11),
                 '50' : 6.86 * 10**(-11), 
                 '100' : 1.19 * 10**(-10),
                 '200' : 1.75 * 10**(-10), 
                 '500' : 2.66 * 10**(-10),
                 '1000' : 3.33 * 10**(-10), 
                 }  
        x2 = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        y2 = [points[str(i)] for i in x2]
        stepsize = (T_k_max - T_k_min) / float(step)
        x = [T_k_min + float(i) * stepsize for i in range(0, step)]
        y1 = [self.kappa_HH(T_k) for T_k in x]
        #for y in y1:
        #    print y
        plot1 = plt.plot(x, y1, label = r'$\kappa_{1-0}^{HH}$')
        plot2 = plt.plot(x2, y2, 'ro',label = r'data')
        plt.legend(loc = 'upper left')
        plt.title(r'$\kappa_{1-0}^{HH}$ vs $T_k$' )
        plt.xlabel(r'$T_k$')
        plt.ylabel(r'$\kappa_{1-0}^{HH}$')
        plt.grid(True)

        plt.savefig('kappa_HH.png')
        plt.show()

    def plot_kappa_Hp(self, T_k_min, T_k_max, step):
        points = {'1' : 0.4028,
                         '2' : 0.4517, 
                         '5' : 0.4301,
                         '10' : 0.3699, 
                         '20' : 0.3172,
                         '50' : 0.3047, 
                         '100' : 0.3379,
                         '200' : 0.4043, 
                         '500' : 0.5471,
                         '1000' : 0.7051, 
                         '2000' : 0.9167,
                         '3000' : 1.070, 
                         '5000' : 1.301,
                         '7000' : 1.480, 
                         '10000' : 1.695,
                         '15000' : 1.975, 
                         '20000' : 2.201 }   
        x2 = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,
                3000, 5000, 7000, 10000, 15000, 20000]
        y2 = [points[str(i)] for i in x2]
        stepsize = (T_k_max - T_k_min) / float(step)
        x = [T_k_min + float(i) * stepsize for i in range(0, step)]
        y1 = [self.kappa_Hp(T_k) for T_k in x]
        #for y in y1:
        #    print y
        plot1 = plt.plot(x, y1, label = r'$\kappa_{1-0}^{Hp}$')
        plot2 = plt.plot(x2, y2, 'ro',label = r'data')
        plt.legend(loc = 'upper left')
        plt.title(r'$\kappa_{1-0}^{Hp}$ vs $T_k$' )
        plt.xlabel(r'$T_k$')
        plt.ylabel(r'$\kappa_{1-0}^{Hp}$')
        plt.grid(True)
        
        plt.savefig('kappa_Hp.png')
        plt.show()

    def plot_kappa_He(self, T_k_min, T_k_max, step):
        points = {'1' : 0.239,
                         '2' : 0.337, 
                         '5' : 0.503,
                         '10' : 0.746, 
                         '20' : 1.05,
                         '50' : 1.63, 
                         '100' : 2.26,
                         '200' : 3.11, 
                         '500' : 4.59,
                         '1000' : 5.92, 
                         '2000' : 7.15,
                         '3000' : 7.71, 
                         '5000' : 8.17,
                         '7000' : 8.32, 
                         '10000' : 8.37,
                         '15000' : 8.29, 
                         '20000' : 8.11 }  
        x2 = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,
                3000, 5000, 7000, 10000, 15000, 20000]
        y2 = [points[str(i)] for i in x2]
        stepsize = (T_k_max - T_k_min) / float(step)
        x = [T_k_min + float(i) * stepsize for i in range(0, step)]
        y1 = [self.kappa_He(T_k) for T_k in x]
        #for y in y1:
        #    print y
        plot1 = plt.plot(x, y1, label = r'$\kappa_{1-0}^{He}$')
        plot2 = plt.plot(x2, y2, 'ro',label = r'data')
        plt.legend(loc = 'upper left')
        plt.title(r'$\kappa_{1-0}^{He}$ vs $T_k$' )
        plt.xlabel(r'$T_k$')
        plt.ylabel(r'$\kappa_{1-0}^{He}$')
        plt.grid(True)
        
        plt.savefig('kappa_He.png')
        plt.show()

    def plot_T_b(self, zmin, zmax, step):
        stepsize = (zmax - zmin) / float(step)
        x = [zmin + float(i) * stepsize for i in range(0, step)]
        y1 = [self.T_b_0th_Order(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$T_b$')
        plt.legend(loc = 'upper left')
        plt.title(r'$T_b$ vs $z$' )
        plt.xlabel(r'$z$')
        plt.ylabel(r'$T_b$')
        plt.grid(True)
        
        plt.savefig('T_b.png')
        plt.show()

    def plot_T_k(self, zmin, zmax, step):
        stepsize = (zmax - zmin) / float(step)
        x = [zmin + float(i) * stepsize for i in range(0, step)]
        y1 = [self.T_k(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$T_k$')
        plt.legend(loc = 'upper left')
        plt.title(r'$T_k$ vs $z$' )
        plt.xlabel(r'$z$')
        plt.ylabel(r'$T_k$')
        plt.grid(True)
        
        plt.savefig('T_k.png')
        plt.show()
