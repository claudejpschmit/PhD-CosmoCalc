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
        plot1 = plt.loglog(x, y1, basex = 10, basey = 10, label = r'$T_k$')
        plt.legend(loc = 'upper left')
        plt.title(r'$T_k$ vs $z$' )
        plt.xlabel(r'$z$')
        plt.ylabel(r'$T_k$')
        plt.grid(True)
        plt.ylim([10, 1000])
        
        plt.savefig('T_k.png')
        plt.show()
