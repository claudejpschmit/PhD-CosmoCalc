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
        help = 'Hubble Constant [km/s/Mpc], default is H_0 = 70.0')
parser.add_argument('--O_M', metavar = 'O_M', 
        type = float, default = 0.3, help = 'Matter density, default is O_M = 0.3')
parser.add_argument('--O_V', metavar = 'O_V', 
        type = float, default = 0.7, help = 'Vacuum density, default is O_V = 0.7')
parser.add_argument('--O_b', metavar = 'O_b',
        type = float, default = 0.046, help = 'Baryon density, default is O_b = 0.046')
parser.add_argument('--z', metavar = 'z', 
        type = float, default = 3, help = 'redshift, default is z = 3')
parser.add_argument('--z_low', metavar = 'z_low', 
        type = float, default = 7, help = 'lower bound for z integration, default is z_low = 7')
parser.add_argument('--z_high', metavar = 'z_high', 
        type = float, default = 9, help = 'top bound for z integration, default is z_high = 9')
parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature, default is T_CMB = 2.75')

args = parser.parse_args()
z = args.z

################# Output ################## 
# generate parameters
h = args.H_0 / 100.0
ombh2 = args.O_b * h**2
omch2 = (args.O_M - args.O_b) * h**2
omk = 1.0 - args.O_M - args.O_V
omnuh2 = 0.00064
params = {"ombh2":ombh2, "omch2":omch2, "omnuh2":omnuh2, "omk":omk, "hubble":args.H_0, "zmin":args.z_low, "zmax":args.z_high, "T_CMB":args.T_CMB}
# initialize plotter
plotter = CosmoPlot(params)
plotter.plot_P_camb_5points(0.001, 10, 10000, "ombh2", 0.005)


#plotter.plot_dTb(5, 20,100)
#plotter.plot_xHI(0, 20, 100)
#plotter.plot_Ts(0,20,100)
#plotter.plot_Tk(0,1000,1000)
#plotter.plot_P_camb(0.001, 10, 10000)
#plotter.plot_distances()
#plotter.plot_bessel_camb(1, 0, 100, 0.1)
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
#plotter.plot_P_growth(0,10,100)
#plotter.plot_P(0.001,10, 10000, units_k = 'default', units_P = 'default')
