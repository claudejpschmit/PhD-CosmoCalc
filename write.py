##################################################################
# This file handles plotting for the Cosmology Calculator class. #
##################################################################

from CosmologyWriterClass import CosmoWrite
import argparse

############## Parsing input ##############

descr = 'This program uses the Cosmology Calculator Class \
        to write the function M_l to a file. This program \
        gives the values of M_l(k1,k2) for a specified l & k1 \
        over a range of k2. The following parameters should \
        be defined: l, k_fixed, k_low, k_high, steps'

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
parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature, default is T_CMB = 2.75')
parser.add_argument('--l', metavar = 'l', 
        type = int, default = 5, help = 'Spherical Bessel index, default is l = 5')
parser.add_argument('--k_fixed', metavar = 'k_fixed', 
        type = float, default = 0.1, help = 'k1, default is k_fixed = 0.1')
parser.add_argument('--k2_low', metavar = 'k2_low', 
        type = float, default = 0.01, help = 'lower bound for k2, default is k2_low = 0.01')
parser.add_argument('--k2_high', metavar = 'k2_high', 
        type = float, default = 0.5, help = 'top bound for k2, default is k2_high = 0.5')
parser.add_argument('--steps', metavar = 'steps', 
        type = int, default = 10000, help = 'Number of steps between k_low and k_high, default is steps = 10000')
parser.add_argument('--stepsize', metavar = 'stepsize', 
        type = float, default = 0, help = 'Stepsize. Non-zero value leads to overwriting steps variable, default is H_0 = 70.0')
parser.add_argument('--z_low', metavar = 'z_low', 
        type = float, default = 7, help = 'lower bound for z integration, default is z_low = 7')
parser.add_argument('--z_high', metavar = 'z_high', 
        type = float, default = 9, help = 'top bound for z integration, default is z_high = 9')
parser.add_argument('--method', metavar = 'method', 
        type = int, default = 2, help = 'Integration method, either 0 for scipy or 1 for mp or 2 for simpson, default is method = 2')

args = parser.parse_args()

################# Output ################## 
# generate parameters
h = args.H_0 / 100.0
ombh2 = args.O_b * h**2
omch2 = (args.O_M - args.O_b) * h**2
omk = 1.0 - args.O_M - args.O_V
omnuh2 = 0.00064
params = {"ombh2":ombh2, "omch2":omch2, "omnuh2":omnuh2, "omk":omk, "hubble":args.H_0}
# initialize writer
writer = CosmoWrite(params, args.z_low, args.z_high, args.T_CMB)

if args.method == 0:
    writer.calculate_Ml_scipy(args.l, args.k_fixed, args.k2_low, args.k2_high, args.z_low, args.z_high, args.steps, args.stepsize)
elif args.method == 1:
    writer.calculate_Ml_mp(args.l, args.k_fixed, args.k2_low, args.k2_high, args.z_low, args.z_high, args.steps, args.stepsize)
elif args.method == 2:
    writer.calculate_Ml(args.l, args.k_fixed, args.k2_low, args.k2_high, args.steps, args.stepsize)
