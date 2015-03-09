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

from CosmoBasis import CosmoBasis

from math import *
import scipy.integrate as integrate
import numpy as np
import mpmath as mp
#from skmonaco import mcquad

class CosmoCalc(CosmoBasis):
    
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
                    ( self.D_M(z2) * (1 + sqrt(1 - O_k) *\
                    self.D_M(z)**2 / self.D_H**2) - \
                    self.D_M(z) * (1 + sqrt(1 - O_k) *\
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
    # This has units [MPc^6 K^2]
    def corr_Tb(self, l, k1, k2, z_low, z_high, k_low, k_high):
        #TODO: figure out integration limits kkk
        integral = self.integrate_simps(lambda k: k**2 *\
                    self.P_delta(k, units_k = 'mpc-1', units_P = 'mpc3') *\
                    self.M(l, k1, k, z_low, z_high) *\
                    self.M(l, k2, k, z_low, z_high), k_low, k_high, 1000)
        return integral

#########################################################################

    def Tb_integrand(self, l, k1, k2, k3, z1, z2):
        res = k3**2*self.P_delta(k3, units_k = 'mpc-1', units_P = 'mpc3') *\
                    self.M_integrand(l, k1, k3, z1) *\
                    self.M_integrand(l, k2, k3, z2)
  
        return res


##########################################################################
########################################################################
# this function uses simple Simpson integration and interpolates spherical Bessels.

    def M(self, l, k1, k2, z_low, z_high):
        def integrand(z):
            r = self.D_C(z)
            integral = integrate.quad(lambda x: (1+x) / self.E(x)**3, z, np.inf)

            return r**2 * self.delta_Tb_bar(r) * self.sphbess_interp(l,k1*r) *\
                    self.sphbess_interp(l,k2*r)  * integral[0]
        
        integral = self.integrate_simps(lambda z: integrand(z), z_low, z_high, 100)

        B = integrate.quad(lambda x: (1+x) / self.E(x)**3, 0, np.inf) 
        prefactor = 2*self.b_bias*self.c/(B[0]*pi*self.H_0*1000)
        res = prefactor * integral  
        return res


#########################################################################


########################################################################
# this function uses scipy integration

    def M_scipy(self, l, k1, k2, z_low, z_high):
        def integrand(z):
            r = self.D_C(z)
            integral = integrate.quad(lambda x: (1+x) / self.E(x)**3, z, np.inf)

            return r**2 * self.delta_Tb_bar(r) * self.sphbess_interp(l,k1*r) * self.sphbess_interp(l,k2*r)* integral[0]
        
        integral = integrate.quad(lambda z: integrand(z), z_low, z_high)

        B = integrate.quad(lambda x: (1+x) / self.E(x)**3, 0, np.inf) 
        prefactor = 2*self.b_bias*self.c/(B[0] *pi*self.H_0*1000)
        res = prefactor * integral[0]  
        return res


#########################################################################
#########################################################################
# this function uses mpmath integration

    def M_mp(self, l, k1, k2, z_low, z_high):
        def integrand(z):
            r = self.D_C(z)
            integral = integrate.quad(lambda x: (1+x) / self.E(x)**3, z, np.inf)

            return r**2 * self.delta_Tb_bar(r) * self.sphbess_interp(l,k1*r) * self.sphbess_interp(l,k2*r)* integral[0] 

        integral = mp.quad(lambda z: integrand(z), [z_low, z_high])

        B = integrate.quad(lambda x: (1+x) / self.E(x)**3, 0, np.inf) 
        prefactor = 2*self.b_bias*self.c/(B[0] *pi*self.H_0*1000)
        res = prefactor * integral  
        return res


#########################################################################
########################################################################
# this function uses scipy integration and no growth function

    def M_scipy_ng(self, l, k1, k2, z_low, z_high):
        def integrand(z):
            r = self.D_C(z)
            
            return r**2 * self.delta_Tb_bar(r) * self.sphbess_interp(l,k1*r) * self.sphbess_interp(l,k2*r) / self.E(z)         
        
        integral = integrate.romberg(lambda z: integrand(z), z_low, z_high, divmax = 30)

        
        prefactor = 2*self.b_bias*self.c/(pi*self.H_0*1000)
        res = prefactor * integral  
        return res


#########################################################################
#########################################################################
# this function uses mpmath integration and no growth function

    def M_mp_ng(self, l, k1, k2, z_low, z_high):
        def integrand(z):
            r = self.D_C(z)
            return r**2 * self.delta_Tb_bar(r) * self.sphbess_interp(l,k1*r) * self.sphbess_interp(l,k2*r) / self.E(z)  

        integral = mp.quad(lambda z: integrand(z), [z_low, z_high])
        
        prefactor = 2*self.b_bias*self.c/(pi*self.H_0*1000)
        res = prefactor * integral  
        return res


#########################################################################

    # Mean Brightness Temperature fluctuations at distance r (comoving) [K]
    def delta_Tb_bar(self, r):
        #TODO: find sensible constant
        constant = 1
        return constant


    # we seperate the power spectrum P(k,a) into a growing mode and P(k)
    # P(k,a) = P_delta(k) * P_growth(a)
    def P_growth(self, z):
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
    # Units: 
    # output:units_P = default:        P in [h^-3 * MPc^3]
    #        units_P = mpc3 or Mpc3:   P in [MPc^3]
    # input: units_k = default:        k in [h Mpc^-1]
    #        units_k = mpc-1 or Mpc-1: k in [MPc^-1]
    
    def P_delta(self, k, units_k = 'default', units_P = 'default'):
        
        if units_P == 'default':

            if units_k == 'default':
                keq = 0.073 * self.O_M * self.h 
                k_factor = k
            elif units_k == 'Mpc-1' or units_k == 'mpc-1':
                keq = 0.073 * self.O_M * self.h**2
                k_factor = k / self.h
        elif units_P == 'Mpc3' or units_P == 'mpc3':
            if units_k == 'default':
                keq = 0.073 * self.O_M * self.h 
                k_factor = k / self.h**3
            elif units_k == 'Mpc-1' or units_k == 'mpc-1':
                keq = 0.073 * self.O_M * self.h**2
                k_factor = k / self.h**4
        
        if self.Omega_M == 1:
            delta_H = 1.9 * 10**(-5)
        else:
            delta_H = 4.5 * 10**(-5)
        A = 2*pi**2 * (self.c/1000)**4 * delta_H**2 / 100**4
        x = k / keq
        transfer_function_sq = self.transfer(x)**2
        P = A * k_factor * transfer_function_sq

        return P 

    # BBKS transfer function, Dodelson (7.70)
    def transfer(self, x):
        res = np.log(1 + 0.171 * x) / (0.171 * x)
        bracket = 1+ 0.284 * x + (1.18 * x)**2 + (0.399 * x)**3 + (0.490 * x)**4
        res = res * bracket**(-0.25)
        return res
    
    # Calculates the second  and third 
    # term in the full Tb correlation 
    # between redshifts z1 & z2.
    # use eg. z1 = 6 & z2 = 10 
    # This has units []
    def corr_Tb_term23(self, l, k1, k2, z_low, z_high, k_low, k_high):
        #TODO: figure out integration limits kkk
        integral = integrate.quad(lambda k: k**(-2) *\
                    self.P_delta(k, units_k = 'mpc-1', units_P = 'mpc3') *\
                    self.M(l, k1, k, z_low, z_high) *\
                    self.N(l, k2, k, z_low, z_high, k_low, k_high), k_low, k_high)
        beta = 1

        return beta*integral[0]
    
    # Calculates the fourth 
    # term in the full Tb correlation 
    # between redshifts z1 & z2.
    # use eg. z1 = 6 & z2 = 10 
    # This has units []
    def corr_Tb_term4(self, l, k1, k2, z_low, z_high, k_low, k_high):
        #TODO: figure out integration limits kkk
        integral = integrate.quad(lambda k: k**(-6) *\
                    self.P_delta(k, units_k = 'mpc-1', units_P = 'mpc3') *\
                    self.N(l, k1, k, z_low, z_high, k_low, k_high) *\
                    self.N(l, k2, k, z_low, z_high, k_low, k_high), k_low, k_high)
        beta = 1
        return beta**2 * integral[0]

    def N_bar(self, l, k1,k2, z_low, z_high, k_low, k_high):
        def integrand(z):
            r = self.D_C(z)
            pref = 1.0/(self.E(z)*1000.0*self.H_0)**2
            sums = r**2 * self.I(l-1, l-1, k1, k2, r) - r * (l+1)/ k1 * self.I(l-1, l, k1, k2, r) - r * (l+1) / k2 * self.I(l, l-1, k1, k2, r) + (l+1)**2 / (k1*k2) * self.I(l, l, k1, k2, r)
            return pref * sums

        integral = mp.quad(lambda z: integrand(z), [z_low, z_high])


    def I(self, l1, l2, k1, k2, z):
        
        r = self.D_C(z)
        prefactor = 2*self.b_bias*self.c/pi
        res = self.delta_Tb_bar(r) * self.sphbess_interp(l1,k1*r) * self.sphbess_interp(l2,k2*r)
        #integral = integrate.quad(lambda x: (1+x) / self.E(x)**3, z, np.inf)
        
        return prefactor * res
