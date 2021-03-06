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
import scipy.interpolate as interpolate
import numpy as np
from sys import getsizeof
import copy

class CosmoCalc(CosmoBasis):
   
    def __init__(self, params):
        CosmoBasis.__init__(self, params)

        # Creating a list of Comoving distances between z_low & z_high to speed up M_l integration
        # Also creating list of growth functions
        self.zmin_Ml = self.fiducial_params["zmin"]
        self.zmax_Ml = self.fiducial_params["zmax"]
        self.zsteps_Ml = self.fiducial_params["zsteps"] 
        self.stepsize_Ml = (self.zmax_Ml - self.zmin_Ml)/float(self.zsteps_Ml)
        self.Pk_steps = self.fiducial_params["Pk_steps"]
        self.k_steps = self.fiducial_params["k_steps"]

        # define the fiducial model r values.
        # We define the q_Ml to start out equal to the fiducial values and then 
        # set r_Ml to be the fiducial comoving distances as well.
        self.q_Ml = []
        self.update_q()
        self.r_Ml = copy.deepcopy(self.q_Ml)

        self.H_f = []
        for n in range(0,self.zsteps_Ml + 1):
            z = self.zmin_Ml + n * self.stepsize_Ml
            self.H_f.append(self.H(z))
         

        # factor of 1000 in order to get the units right later on
        # ie Ml wil have units [Mpc^3 * Mpc^(3/2) mK]
        self.prefactor_Ml = 2*self.b_bias*self.c/pi

        # Initialize a list of parameters that have already been used. 
        # Also, initialise a list of PK interpolator functions.
        # The index should relate both lists.
        self.Pk_params_used = []
        self.Pk_interpolators_used = []
        
        #defining default interpolator
        # This calculates the matter power spectrum in units [h^(-3)Mpc^3]
        #self.Pk_interp = self.Pk_update_interpolator(self.params, self.zmin_Ml, self.zmax_Ml, 3)
        self.Pk_update_interpolator(self.fiducial_params)

        # Let's user know that initialisation is done
        print "CosmoCalc initialized"

#############################################################################

    # Method to update q(z), using the current parameters in the Calculator.
    def update_q(self):
        # This function has to be called after updating the parameters.
               
        self.q_Ml = []
        for n in range(0,self.zsteps_Ml + 1):
            z = self.zmin_Ml + n * self.stepsize_Ml
            self.q_Ml.append(self.D_C(z))
        return None 
    
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
        #O_k = 1-self.O_M-self.O_V
        if (self.O_k == 0.0):
            return self.D_C(z)
        else:
            return self.D_H / sqrt(abs(self.O_k)) * \
                self.S_k(sqrt(abs(self.O_k)) * \
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
        #O_k = 1-self.O_M-self.O_V
        root = sqrt(abs(self.O_k))
        print root 
        if z2 is None:
            if (self.O_tot == 1.0):
                return self.D_H  * self.Z(z)/ (1+z)
            else:
                return self.D_H * self.S_k((root) * self.Z(z))/ \
                    ((1+z) * root)
        elif (self.O_tot <= 1.0):
            return 1.0 / (1.0 + z2) * \
                    ( self.D_M(z2) * sqrt(1 + self.O_k *\
                    self.D_M(z)**2 / self.D_H**2) - \
                    self.D_M(z) * sqrt(1 + self.O_k *\
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
    
    # Basic Cosmological parameter functions.
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
        return self.x_HI(z)
    # (free) electron density
    def n_e(self, z):
        return self.n_p(z)

    ################## Power Spectrum Calculations #####################
    #TODO: Also, check for units, which ones will we need? May have to introduce norm. factors.
    # This function takes a set of Cosmological parameters and uses CAMB 
    # to compute the Power spectrum for those. The results are then used to update the interpolation
    # function for the power spectrum P(k,z).
    # Important! The resulting function takes input k values in units of [h Mpc^-1] and gives results
    #   in units of [h^-3 Mpc^3]
    def Pk_update_interpolator(self, params):

        #loop through a list of all sets of parameters used thus far
        #if the set of parameters has already been used, set f to the 
        #associated interpolation function. If not, calculate new
        #function and add it to the list of interpolated functions and return it.
        #also add the set of parameters to the list of used parameters.
        
        # This function now includes optimization so to avoid unecessary
        # CAMB calls for parameter values that have already been calculated in the past.
        
        index = 0
        while (index <= len(self.Pk_params_used) - 1):
            if (self.Pk_params_used[index] != params):
                index += 1
            else:
                break

        if index == len(self.Pk_params_used):
            #add new Pk interpolator
            table = []
            table_z = []
            zi = self.zmin_Ml
            zf = self.zmax_Ml
            nsteps = self.Pk_steps
            z_stepsize = float(zf-zi)/float(nsteps)

            for n in range(0, nsteps + 1):
                # create array with all z values
                z = zi + n * z_stepsize
                table_z.append(z)

                # generate p(k,z) from camb
                camb_result_dict = self.camb(**{'ombh2':params["ombh2"], 'omch2':params["omch2"],\
                        'omnuh2':params["omnuh2"], 'omk':params["omk"], 'hubble':params["hubble"],\
                        'transfer_redshift(1)':z})
                res = camb_result_dict["transfer_matterpower"]
           
                table.append(res)
        
            #number of kvalues
            nkvals = len(table[0])
            table_k = []
            for n in range(0,nkvals):
                #store all kvalues into table_k
                table_k.append(table[0][n][0])

            table_Pkz = []
            for m in range(0,nsteps+1):
                row = []
                for n in range(0,nkvals):
                    row.append(table[m][n][1])
                table_Pkz.append(row)

            f = interpolate.interp2d(table_k, table_z, table_Pkz, kind = 'cubic')
           
            # Now the parameters and associated interpolation function are stored 
            # for later use. This should not be a memory problem, one such function
            # has a size of 64bytes.

            # Note on copy.deepcopy: This is necessary as Python only ever passes 
            # references to the memory address, so if we only put params, and then later
            # call this function with a different set of parameters, Pk_params_used points
            # to the parameters passed to the function and has essentially forgotten about
            # the parameters that came before. ie. without the copy statement we append
            # the memory address rather than the value found at that address.
            self.Pk_params_used.append(copy.deepcopy(params))
            self.Pk_interpolators_used.append(copy.deepcopy(f))
            #print "new interpolator created"
        else:
            f = self.Pk_interpolators_used[index]
            #print "old Pk used"



        # finally set the interpolator to what we need.

        self.Pk_interp = f
        #self.Pk_interp = interpolate.interp2d(table_k, table_z, table_Pkz, kind = 'cubic')

        return None

    # ## The following defines methods to calculate the Power Spectrum from scratch ##
    def Pkz_calc(self, k, z):
        return self.P_growth(z) * self.P_delta(k)
    
    # we seperate the power spectrum P(k,a) into a growing mode and P(k)
    # P(k,z) = P_delta(k) * P_growth(z)
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
        
        if self.O_M == 1:
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
        bracket = 1 + 0.284 * x + (1.18 * x)**2 + (0.399 * x)**3 + (0.490 * x)**4
        res = res * bracket**(-0.25)
        return res 

    ####################################################################
    ####### Two point Correlation of brightness tem. fluctuations ######
    ####################################################################

    ############### first we ignore redshift space distortions #########
   
    # Calculates the Tb correlation without redshift space distortions 
    # between redshifts z1 & z2.
    # use eg. z1 = 6 & z2 = 10 
    # This has units [MPc^6 mK^2]
    def corr_Tb(self, l, k1, k2, k_low, k_high):
        integral = self.integrate_simps(lambda k: k**2 *\
                    self.M(l, k1, k) *\
                    self.M(l, k2, k), k_low, k_high, self.k_steps)
        return integral
    
    ############### Then, we include redshift space distortions #########
   
    # Calculates the Tb correlation with redshift space distortions 
    # between redshifts z1 & z2.
    # use eg. z1 = 6 & z2 = 10 
    # This has units [MPc^6 mK^2]
    def corr_Tb_rsd(self, l, k1, k2, k_low, k_high):
        def integrand(k):
            m1 = self.M(l, k1, k)
            n1 = self.N_bar(l, k1, k)
            if k1 == k2:
                m2 = m1
                n2 = n1
            else:
                m2 = self.M(l, k2, k)
                n2 = self.N_bar(l, k2, k)

            bb = self.b_bias * self.beta
            bb2 = bb**2
            res = k**2 * m1 * m2 + bb * k * (m1 * n2 + n1 * m2) + bb2 * n1 * n2
            return res

        integral = self.integrate_simps(lambda k: integrand(k), k_low, k_high, self.k_steps)
        return integral

    def Cl(self, l, k1, k2, k_low, k_high):
        return self.corr_Tb_rsd(l, k1, k2, k_low, k_high)

    ########################################################################
    
    ############### Then we define the M_l and N_l matrices. ###############

    # This computes Ml in units [Mpc^3 * Mpc^{3/2} * mK]
    def M(self, l, k1, k2):

        def integrand(z):
            n_old = (z - self.zmin_Ml)/self.stepsize_Ml
            if abs(n_old-int(n_old)) > 0.5:
                n = int(n_old)+1
            else:
                n = int(n_old)
            r = self.r_Ml[n]
            q = self.q_Ml[n]
          
            # Note that we are dividing by the fiducial H(z).
            # The factor of 1000 is there to change from km/s/Mpc to m/s/Mpc
            return r**2 * self.delta_Tb_bar(z) * self.sphbess_camb(l,k1*r) *\
                    self.sphbess_camb(l,k2*q)  * sqrt(self.Pk_interp(k2*self.h,z)/self.h**3) /\
                    (self.H_f[n]*1000.0)
        
        integral = self.integrate_simps(lambda z: integrand(z), self.zmin_Ml, self.zmax_Ml, self.zsteps_Ml)

        res = self.prefactor_Ml * integral  
        return res

    # This computes Nl_bar in units [Mpc^2 * Mpc^{3/2} * mK]
    def N_bar(self, l, k1,k2):
        
        def integrand(z):
            n_old = (z - self.zmin_Ml)/self.stepsize_Ml
            if abs(n_old-int(n_old)) > 0.5:
                n = int(n_old)+1
            else:
                n = int(n_old)
            r = self.r_Ml[n]
            q = self.q_Ml[n]     
           
            # Note that we are dividing by the fiducial H(z).
            # The factor of 1000 is there to change from km/s/Mpc to m/s/Mpc
            pref = 1.0/(self.H_f[n]*1000.0*(1.0 + z)) * self.prefactor_Ml

            pk = sqrt(self.Pk_interp(k2*self.h, z)/self.h**3)
            dtb = self.delta_Tb_bar(z)
            pkdtb = pk * dtb
            jl1r = self.sphbess_camb(l - 1, k1 * r)
            jl2r = self.sphbess_camb(l, k1 * r)
            jl1q = self.sphbess_camb(l - 1, k2 * q)
            jl2q = self.sphbess_camb(l, k2 * q)

            sums = k1 * r * jl1r * jl1q -\
                   k1 * r * (l+1) / (k2 * q) * jl1r * jl2q -\
                   (l+1) * jl2r * jl1q +\
                   (l+1)**2 / (k2 * q) * jl2r * jl2q 
            
            return pref * r * pkdtb * sums
        
        integral = self.integrate_simps(lambda z: integrand(z), self.zmin_Ml, self.zmax_Ml, self.zsteps_Ml)

        return integral
    
       
    #########################################################################
    # Mean Brightness Temperature fluctuations at distance r (comoving) [mK]
    # This does not include the peculiar velocity in form of the velocity gradient, as that 
    # is separately introduced into the formalism
    def delta_Tb_bar(self, z):
        constant_A = 27*self.O_b * self.h**2 / 0.023 * sqrt(0.015/(self.O_M * self.h**2))
        x_HI = self.x_HI(z)
        T_S = self.T_S(z)
        T_g = self.T(z)

        return constant_A * x_HI * (T_S - T_g)/T_S * sqrt(1+z)
        #return constant_A * x_HI * sqrt(1+z)

    def T_S(self, z):
        Ti = 10.0 #something
        Tf = 500.0 #something
        rate = 2.0 #something
        zi = self.z_rei
        zf = zi - self.delta_z_rei #z at end of reionization
        deltaT = abs(Ti-Tf)

        return deltaT * (1/pi * atan(rate * (-(z - (zi + zf)/2.0)))+0.5) + Ti

    def x_HI(self, z):
        rate = 2.0 #something
        zi = self.z_rei
        zf = zi - self.delta_z_rei #z at end of reionization

        return 1/pi * atan(rate * (z - (zi + zf)/2.0))+0.5
    
    def T_K(self, z):
        zd = 200.0 #something
        if z > zd:
            res = self.T(z)
        else:
            Td = self.T(zd) 
            Tf = 1000.0 #something
            rate = 1.0 #something
            z_on = self.z_rei
            z0 = (2*z_on - 4)/2.0
            tanh_term = (0.5*(tanh(rate * (-(z - z0)))+1)) * Tf 
            res = Td * (1+z)**2/(1+zd)**2 + tanh_term
        
        return res
