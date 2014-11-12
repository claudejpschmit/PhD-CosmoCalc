# This file contains a calculator for basic cosmological calculations 
# similar to Ned Wrigths' calulator: 
# http://www.astro.ucla.edu/~wright/CosmoCalc.html 
# See also, http://arxiv.org/pdf/astro-ph/9905116v4.pdf 
# for a review of the formulas used. 
# Or Distance measures in Cosmology, David W. Hogg

from math import *
import sys
import scipy.integrate as integrate
import numpy

class CosmoCalc(object):
    # Initializer
    def __init__(self, H_0, O_M, O_V):
        self.H_0 = H_0
        self.O_M = O_M
        self.O_V = O_V
        # Fixing the radiation density in all cases
        self.O_R = 4.165E-1/self.H_0**2
        # Fixing c [m/s]
        self.c = 299792458.0
        self.O_tot = self.O_M + self.O_V + self.O_R
        # Hubble distance [Mpc]
        self.D_H = self.c / (1000.0 * H_0)
        # Hubble time
        self.t_H = 1.0 / H_0
    
 ############################################
    
    # Helper function, denominator for various integrals
    def E(self, z):
        return sqrt(self.O_V + self.O_R * (1+z)**4 +
                self.O_M * (1+z)**3 + (1-self.O_tot) * (1+z)**2)

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
    #   D_C = D_H int(dz / E(z),0,z)
    def comoving_radial_dist(self, z):
        integral = integrate.quad(lambda x: 1.0/self.E(x), 0, z)
        return self.hubble_dist() * integral[0]
    def D_C(self, z):
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
    def ang_diam2(self, z):
        pass    
    def angular_diam_dist(self, z, z2 = None):
        O_k = 1-self.O_M-self.O_V
        if z2 is None:
            return self.D_M(z) / (1.0 + z)
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
    def comoving_volume(self, z):
        O_k = 1-self.O_M-self.O_V
        if (self.O_V + self.O_M < 1.0):
            return 4.0 * pi * self.D_H**3 / (2 * (1-O_k)) * \
                   (self.D_M(z)/self.D_H * \
                   sqrt(1 + (1 - O_k) * self.D_M(z)**2 / \
                   self.D_H**2) - 1 / sqrt(abs(1-O_k)) * \
                   asinh(sqrt(abs(1-O_k)) * self.D_M(z) / \
                   self.D_H))
        elif (self.O_V + self.O_M == 1.0):
            return 4.0 * pi * self.D_M(z)**3 / 3.0 
        else:
            return 4.0 * pi * self.D_H**3 / (2 * abs(1-O_k)) * \
                   (self.D_M(z)/self.D_H * \
                   sqrt(1 + (1 - O_k) * self.D_M(z)**2 / \
                   self.D_H**2) - 1 / sqrt(abs(1-O_k)) * \
                   asin(sqrt(abs(1-O_k)) * self.D_M(z) / \
                   self.D_H))
    def V_C(self, z):
        return self.comoving_volume(z)
    def com_vol2(self,z):
        vol = 4*pi*self.D_H
        integral = integrate.quad(lambda x: (1+x)**2 * \
                self.D_A(x)**2/self.E(x), 0, z)
        return vol*integral[0]

    # Comoving volume element:
    # dV_C = D_H (1+z)^2 * D_A^2 / E(z) dOmega dz
    ##########

    # Age of the Universe at redshift z [s * Mpc/km]:
    #   t(z) = t_H int(dz/((1+z) E(z)), z, inf)
    #def age_of_universe(self, z):
    #    age = integrate.quad(lambda a: 1.0/self.D(a), 0.0, 1.0/(1.0+z))
     #   return self.convert_to_gy(age[0]/H_0)
    def age_of_universe(self, z):
        age = integrate.quad(lambda x: 1.0/((1+x)*self.E(x)), z, numpy.inf)
        return age[0] * self.t_H

    # Light travel time [s * Mpc/km]:
    #   ltt = t(0) - t(z)
    def light_travel_time(self, z):
        return self.age_of_universe(0) - self.age_of_universe(z)

############### Input/Output ###################

def convert_to_gy(age):
    return age * 10**10 * 3.08568 / (365.25 * 24 * 3600) 
# Setting default parameters
H_0 = 70
O_M = 0.3
z = 3
O_V = 0.7
default = raw_input("Do you wish to use default parameters? (y/n): ")

# User defined parameters
if default.lower() == "n" or default.lower() == "no":
    H_0 = float(raw_input("Insert Hubble Constant H_0: "))
    O_M = float(raw_input("Insert Matter density Omega_M: "))
    z = float(raw_input("Insert redshift z: "))
    O_V = float(raw_input("Insert Vacuum density Omega_V: "))

# Output
calc = CosmoCalc(H_0, O_M, O_V)
print "For a Universe with H0 = %s, Omega_M = %s, Omega_V = %s and z = %s:" % (calc.H_0, calc.O_M, calc.O_V, z)
print "It is now %s Gyr since the Big Bang." % \
        convert_to_gy(calc.age_of_universe(0))
print "The age at redshift z was %s Gyr." % \
        convert_to_gy(calc.age_of_universe(z))
print "The light travel time was %s Gyr." % \
        convert_to_gy(calc.light_travel_time(z))
print "The comoving radial distance is %s MPc." % \
        calc.comoving_radial_dist(z)
vol = calc.comoving_volume(z) / 10**9
vol2 = calc.com_vol2(z) / 10**9
print "The comoving volume within redshift z is %s Gpc^3." % \
        vol
print "The comoving volume using method 2 is %s Gpc^3." % \
        vol2

print "The angular size distance D_A is %s MPc." % \
        calc.angular_diam_dist(z)
print "The luminosity distance D_L is %s MPc." % \
        calc.luminosity_dist(z)
