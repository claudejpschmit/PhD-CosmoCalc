#######################################################################
# This file contains a plotting class for the Cosmology Calcuator     #
#######################################################################

from CosmologyCalculatorClass import CosmoCalc
import matplotlib.pyplot as plt

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

    def plot_P(self, kmin, kmax, step):
        stepsize = (kmax - kmin) / float(step)
        x = [kmin + float(i) * stepsize for i in range(0, step)]
        y1 = [self.P_delta(k) for k in x]
        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'P(k) vs k')
        plt.legend(loc = 'upper left')
        plt.title(r'test' )
        plt.xlabel(r'$k$')
        plt.ylabel(r'$P(k)$')
        plt.grid(True)
        
        plt.savefig('powerspectrum.png')
        plt.show()
