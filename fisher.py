###########################################
# This file handles user input and output #
# for the Fisher Class                    #
###########################################

from FisherClass import Fisher
import argparse
import math

############## Parsing input ##############

descr = 'The following program uses the FisherClass and thus the\
         CosmoCalculator Class to perform Fisher matrix calculations.\n \
         The following parameters can be entered manually, H_0, O_M, O_b,\
         z and T_CMB. Default values are assumed if nothing is entered.'

parser = argparse.ArgumentParser(description=descr)
parser.add_argument('--version', action='version', 
        version='%(prog)s v0.1')
parser.add_argument('--H_0', metavar = 'H_0', 
        type = float, default = 70.0,
        help = 'Hubble Constant [km/s/Mpc], default is H_0 = 70.0')
parser.add_argument('--O_M', metavar = 'O_M', 
        type = float, default = 0.3, help = 'Matter density, default is O_M = 0.3')
parser.add_argument('--O_b', metavar = 'O_b',
        type = float, default = 0.046, help = 'Baryon density, default is O_b = 0.046')
parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature, default is T_CMB = 2.75')
parser.add_argument('--z', metavar = 'z', 
        type = float, default = 3, help = 'redshift, default is z = 3')

args = parser.parse_args()
z = args.z

################# Output ################## 

# generate parameters
h = args.H_0 / 100.0
omnuh2 = 0.00064
ombh2 = args.O_b * h**2
omch2 = (args.O_M - args.O_b) * h**2 - omnuh2
omk = 0.0
omrh2 = 4.165E-1/10000.0
params = {"ombh2":ombh2, "omch2":omch2, "omnuh2":omnuh2, "omk":omk, "hubble":args.H_0}
# initializing the Fisher object
fisher_obj = Fisher(params)


fisher_obj.write_logder("ombh2", 0.0226, 0.00000001, 0.000001, 100, "_very_low_stepsize")
#fisher_obj.compute_Cl(10)
