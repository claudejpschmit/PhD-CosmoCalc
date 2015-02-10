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

    def plot_corr_Tb_no_distortions(self, l_min, l_max, k1, k2, z_low, z_high, step):
        stepsize = (l_max - l_min) / float(step)
        x = [l_min + float(i) * stepsize for i in range(0, step)]
        y1 = [self.corr_Tb_no_distortions(l, k1, k2, z_low, z_high) for l in x]
        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'$\langle \delta T_{b,l}(k)\delta T_{b,l}(k\')\rangle$ vs $l$ for $k = 0.01$ and $k\' = 0.1$')
        plt.legend(loc = 'upper left')
        plt.title(r'$\langle \delta T_{b,l}(k)\delta T_{b,l}(k\')\rangle$ vs $l$' )
        plt.xlabel(r'$l$')
        plt.ylabel(r'$\langle \delta T_{b,l}(k)\delta T_{b,l}(k\')\rangle$')
        plt.grid(True)
        plt.ylim([10, 1000])
        
        plt.savefig('correlation.png')
        plt.show()
    
    
    def plot_Ml(self, l, k_fixed, k2_low, k2_high, z_low, z_high, step):
        stepsize = (k2_high - k2_low) / float(step)
        x = [k2_low + float(i) * stepsize for i in range(0, step)]
        y1 = [self.M(l, k_fixed, k2, z_low, z_high) for k2 in x]
        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'$M_l(k, k\')$')
        plt.legend(loc = 'upper left')
        plt.title(r'$M_l(k, k\')$ vs $k\'$ for $l = 5$' )
        plt.xlabel(r'$k\'$')
        plt.ylabel(r'$M_l(k, k\')$')
        plt.grid(True)
        
        plt.savefig('M_l.png')
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

