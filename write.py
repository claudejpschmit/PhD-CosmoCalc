##################################################################
# This file handles plotting for the Cosmology Calculator class. #
##################################################################

from CosmologyWriterClass import CosmoWrite
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

writer = CosmoWrite(args.H_0, args.O_M, args.O_V, args.T_CMB)

writer.calculate_Ml(l=5, k_fixed=0.1, k2_low=0.2, k2_high=1, z_low=7, z_high=9, step=10)
