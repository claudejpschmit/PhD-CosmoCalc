#######################################################################
# This file contains a plotting class for the Cosmology Calcuator     #
#######################################################################

from CosmologyCalculatorClass import CosmoCalc
import matplotlib.pyplot as plt
import numpy as np
from math import *

class CosmoPlot(CosmoCalc):
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
        y4 = [self.Omega_V(z)+self.Omega_R(z)+self.Omega_M(z) for z in x] 
        plot1 = plt.plot(x, y1, label = r'$\Omega_M$')
        plot2 = plt.plot(x, y2, label = r'$\Omega_R$')
        plot3 = plt.plot(x, y3, label = r'$\Omega_V$')
        plot4 = plt.plot(x, y4, label = r'Sum')
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
        
    def plot_P(self, kmin, kmax, step, units_k = 'default', units_P = 'default'):
        stepsize = (kmax - kmin) / float(step)
        
        x = [kmin + float(i) * stepsize for i in range(0, step)]
        y1 = [self.P_delta(k, units_k, units_P) for k in x]
        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'P(k)')
        
        if units_k == 'default':
            plt.xlabel(r'$k$ $(h \, Mpc^{-1})$')
        elif units_k == 'mpc-1' or units_k == 'Mpc-1':
            plt.xlabel(r'$k$ $(Mpc^{-1})$')
        
        if units_P == 'default':
            plt.ylabel(r'$P(k)$ $(h^{-3} Mpc^3)$')
        elif units_P == 'mpc3' or units_P == 'Mpc3':
            plt.ylabel(r'$P(k)$ $(Mpc^3)$')


        if units_k == 'default' and units_P == 'default':
            xcomp, ycomp = np.loadtxt('compare.dat', unpack = True)
            plot2 = plt.loglog(xcomp, ycomp, basex = 10, basey = 10, label = r'P(k) from iCosmo')
            xcomp2, ycomp2 = np.loadtxt('test_matterpower.dat', unpack = True)
            plot3 = plt.loglog(xcomp2, ycomp2, basex = 10, basey = 10, label = r'P(k) from CAMB')

        plt.legend(loc = 'upper right')
        plt.title(r'Matter Power Spectrum' )
        plt.grid(True)
        plt.xlim([10**(-3), 10])
        plt.ylim([1, 100000])

        
        plt.savefig('powerspectrum.png')
        plt.show()

    def plot_P_growth(self, z_low, z_high, step):
        stepsize = (z_high - z_low) / float(step)
        x = [z_low + float(i) * stepsize for i in range(0, step)]
        y1 = [self.P_growth(z) for z in x]
        y2 = [e**(-z) for z in x]
        
        plot1 = plt.plot(x, y1, label = r'$P_{growth}$')
        plot2 = plt.plot(x, y2, label = r'exponential')
        
        plt.legend(loc = 'upper left')
        plt.title(r'$P_{growth}$ vs $z$' )
        plt.xlabel(r'$z$')
        plt.ylabel(r'$P_{growth}$')
        plt.grid(True)
        
        plt.savefig('P_growth.png')
        plt.show()

    def plot_bessel_camb(self, l, xmin, xmax, stepsize):
        
        nsteps = int((xmax-xmin)/float(stepsize))
        x = [xmin + n*stepsize for n in range(0, nsteps)]
        y = [self.sphbess_camb(l,x1) for x1 in x]

        plot1 = plt.plot(x, y, label = r'$j_l(x)$')
        
        plt.legend(loc = 'upper left')
        plt.title(r'$j_l(x)$ vs $x$' )
        plt.xlabel(r'$x$')
        plt.ylabel(r'$j_l(x)$')
        plt.grid(True)
        plt.show()
    
    def plot_P_camb(self, kmin, kmax, nsteps):
        stepsize = (kmax - kmin) / float(nsteps)
        x = [kmin + float(i) * stepsize for i in range(0, nsteps)]
        y1 = [self.Pk_interp(k,0) for k in x]
        plot1 = plt.scatter(x, y1, label = r'$P(k,z=z_0)$',color='red')
        y2 = [self.Pk_interp(k,0+0.5) for k in x]
        plot2 = plt.scatter(x, y2, label = r'$P(k,z=z_0+0.5)$',color='blue')
        plt.xlabel(r'$k$ $(h \, Mpc^{-1})$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylabel(r'$P(k)$ $(h^{-3} Mpc^3)$')

        plt.legend(loc = 'upper right')
        plt.title(r'Matter Power Spectrum' )
        plt.grid(True)
        plt.xlim([10**(-3), 10])
        plt.ylim([1, 100000])
        
        plt.savefig('camb_powerspectrum.png')
        plt.show()
    
    def plot_P_camb_5points(self, kmin, kmax, nsteps, param_key, variation_size):
        stepsize = (kmax - kmin) / float(nsteps)
        x = [kmin + float(i) * stepsize for i in range(0, nsteps)]
        
        params = self.fiducial_params
        p0 = params[param_key]

        y1 = [self.Pk_interp(k,0) for k in x]
        plot1 = plt.plot(x, y1, label = r'$P(k,z=z_0) at x$',color='red')

        params[param_key] = p0 - variation_size
        self.Pk_update_interpolator(params)       
        y2 = [self.Pk_interp(k,0) for k in x]
        plot2 = plt.plot(x, y2, label = r'$P(k,z=z_0) at x - h$',color='blue')
        
        params[param_key] = p0 - 2*variation_size
        self.Pk_update_interpolator(params)       
        y3 = [self.Pk_interp(k,0) for k in x]
        plot3 = plt.plot(x, y3, label = r'$P(k,z=z_0) at x - 2h$',color='black')
       
        params[param_key] = p0 +2* variation_size
        self.Pk_update_interpolator(params)       
        y4 = [self.Pk_interp(k,0) for k in x]
        plot4 = plt.plot(x, y4, label = r'$P(k,z=z_0) at x + 2h$',color='green')

        params[param_key] = p0 + variation_size
        self.Pk_update_interpolator(params)       
        y5 = [self.Pk_interp(k,0) for k in x]
        plot5 = plt.plot(x, y5, label = r'$P(k,z=z_0) at x + h$',color='yellow')


        plt.xlabel(r'$k$ $(h \, Mpc^{-1})$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylabel(r'$P(k)$ $(h^{-3} Mpc^3)$')

        plt.legend(loc = 'upper right')
        plt.title(r'Matter Power Spectrum' )
        plt.grid(True)
        plt.xlim([10**(-3), 10])
        plt.ylim([1, 100000])
        
        #plt.savefig('camb_powerspectrum.png')
        plt.show()
   
    def plot_dTb(self, zmin, zmax, nsteps):
        stepsize = (zmax - zmin) / float(nsteps)
        x = [zmin + float(i) * stepsize for i in range(0, nsteps)]
        y1 = [self.delta_Tb_bar(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$\delta\bar{T}_b(z)$')
        plt.xlabel(r'$z$ ')
        plt.ylabel(r'$\delta\bar{T}_b(z)$')

        plt.legend(loc = 'upper right')
        plt.title(r'$\delta\bar{T}_b(z)$ vs. $z$' )
        plt.grid(True)
        plt.show()
        
    def plot_xHI(self, zmin, zmax, nsteps):
        stepsize = (zmax - zmin) / float(nsteps)
        x = [zmin + float(i) * stepsize for i in range(0, nsteps)]
        y1 = [self.x_HI(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$x_{HI}(z)$')
        plt.xlabel(r'$z$ ')
        plt.ylabel(r'$x_{HI}(z)$')

        plt.legend(loc = 'upper right')
        plt.title(r'$x_{HI}(z)$ vs. $z$' )
        plt.grid(True)
        plt.show()
    def plot_Ts(self, zmin, zmax, nsteps):
        stepsize = (zmax - zmin) / float(nsteps)
        x = [zmin + float(i) * stepsize for i in range(0, nsteps)]
        y1 = [self.T_S(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$T_{S}(z)$')
        plt.xlabel(r'$z$ ')
        plt.ylabel(r'$T_{S}(z)$')

        plt.legend(loc = 'upper right')
        plt.title(r'$T_{S}(z)$ vs. $z$' )
        plt.grid(True)
        plt.show()
        
    def plot_Tk(self, zmin, zmax, nsteps):
        stepsize = (zmax - zmin) / float(nsteps)
        x = [zmin + float(i) * stepsize for i in range(0, nsteps)]
        y1 = [self.T_K(z) for z in x]
        y2 = [self.T(z) for z in x]
        plot1 = plt.plot(x, y1, label = r'$T_{K}(z)$')
        plot2 = plt.plot(x, y2, label = r'$T_{CMB}(z)$')
        plt.xlabel(r'$z$ ')
        plt.ylabel(r'$T_{K}(z)$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

        plt.legend(loc = 'upper right')
        plt.title(r'$T_{K}(z)$ vs. $z$' )
        plt.grid(True)
        plt.show()
