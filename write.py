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
        help = 'Hubble Constant [km/s/Mpc]')
parser.add_argument('--O_M', metavar = 'O_M', 
        type = float, default = 0.3, help = 'Matter density')
parser.add_argument('--O_V', metavar = 'O_V', 
        type = float, default = 0.7, help = 'Vacuum density')
parser.add_argument('--z', metavar = 'z', 
        type = float, default = 3, help = 'redshift')
parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature')
parser.add_argument('--l', metavar = 'l', 
        type = float, default = 5, help = 'Spherical Bessel index')
parser.add_argument('--k_fixed', metavar = 'k_fixed', 
        type = float, default = 0.1, help = 'k1')
parser.add_argument('--k2_low', metavar = 'k2_low', 
        type = float, default = 0.01, help = 'lower bound for k2')
parser.add_argument('--k2_high', metavar = 'k2_high', 
        type = float, default = 0.5, help = 'top bound for k2')
parser.add_argument('--steps', metavar = 'steps', 
        type = float, default = 10000, help = 'Number of steps between k_low and k_high')
parser.add_argument('--stepsize', metavar = 'stepsize', 
        type = float, default = 0, help = 'Stepsize. Non-zero value leads to averwriting steps variable')
parser.add_argument('--z_low', metavar = 'z_low', 
        type = float, default = 7, help = 'lower bound for z integration')
parser.add_argument('--z_high', metavar = 'z_high', 
        type = float, default = 9, help = 'top bound for z integration')


args = parser.parse_args()

################# Output ################## 

writer = CosmoWrite(args.H_0, args.O_M, args.O_V, args.T_CMB)

#writer.calculate_Ml(args.l, args.k_fixed, args.k2_low, args.k2_high, args.z_low, args.z_high, int(args.steps), args.stepsize)
writer.calculate_Ml_opt(args.l, args.k_fixed, args.k2_low, args.k2_high, args.z_low, args.z_high, int(args.steps), args.stepsize)
