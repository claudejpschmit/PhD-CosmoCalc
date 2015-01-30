##################################################################
# This file handles plotting for the Cosmology Calculator class. #
##################################################################

from CosmologyPlotterClass import CosmoPlot
import argparse

############## Parsing input ##############

descr = 'This program uses the Cosmology Calculator Class \
         to plot various cosmological results including ...'

parser = argparse.ArgumentParser(description=descr)
parser.add_argument('--version', action='version', 
        version='%(prog)s v0.1')
parser.add_argument('--H_0', metavar = 'H_0', 
        type = float, default = 70.0,
        help = 'Hubble Constant [km/s/Mpc]')
parser.add_argument('--O_M', metavar = 'O_M', 
        type = float, default = 0.3, help = 'Matter density')
parser.add_argument('--O_V', metavar = 'O_V', 
        type = float, default = 0.7, help = 'Vacuum density')
parser.add_argument('--z', metavar = 'z', 
        type = float, default = 3, help = 'redshift')
parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature')


args = parser.parse_args()
z = args.z

################# Output ################## 

plotter = CosmoPlot(args.H_0, args.O_M, args.O_V, args.T_CMB)
#plotter.plot_distances()
#plotter.plot_densities_rho(5000, 100)
#plotter.plot_densities_Omega(10000, 100)
#plotter.plot_H(1000, 100)
#plotter.plot_x_HI(50, 100)
#plotter.plot_kappa_HH(1, 1000, 1000)
#plotter.plot_kappa_Hp(1, 20000, 1000)
#plotter.plot_kappa_He(1, 20000, 1000)
#plotter.plot_T_b(5, 200, 1000)
#plotter.plot_T_k(5, 1000, 1000)
#plotter.plot_T_CMB(0,1000, 1000)
#plotter.plot_age(200,1100, 2000)

#### this has errors
#plotter.plot_corr_Tb_no_distortions(1, 1000, 0.01, 0.1, 6, 10, 1000)

#### wrong because amplitude isn't right yet.
plotter.plot_P(0.001,10, 10000, units_k = 'default', units_P = 'default')
