###########################################
# This file handles user input and output #
# for the Cosmology Calculator Class      #
###########################################

from CosmologyCalculatorClass import CosmoCalc
import argparse
import cProfile as profile

############## Parsing input ##############

descr = 'This program uses the Cosmology Calculator Class \
         to calculate various cosmological results including \
         the age of the universe and various cosmological distances.\n \
         The following parameters can be entered manually, H_0, O_M, O_V,\
         z and T_CMB. Default values are assumed if nothing is entered.'

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
parser.add_argument('--z', metavar = 'z', 
        type = float, default = 3, help = 'redshift, default is z = 3')
parser.add_argument('--z_low', metavar = 'z_low', 
        type = float, default = 7, help = 'lower bound for z integration, default is z_low = 7')
parser.add_argument('--z_high', metavar = 'z_high', 
        type = float, default = 9, help = 'top bound for z integration, default is z_high = 9')

parser.add_argument('--T_CMB', metavar = 'T_CMB', 
        type = float, default = 2.75, help = 'CMB temperature, default is T_CMB = 2.75')
parser.add_argument('--p', metavar = 'profiling', 
        type = float, default = False, help = 'Turns on (1) or off (0) program profiling')
parser.add_argument('--bessel', metavar = 'bessel', 
        type = str, default = 'Bessel/CAMB/JL_unformatted.bin', help = 'Filename of spherical bessel table, default is path = Bessel/CAMB/JL_unformatted.bin')

args = parser.parse_args()
z = args.z

################# Output ################## 

# TODO: Maybe include this in the class
def convert_to_gy(age):
    return age * 10**10 * 3.08568 / (365.25 * 24 * 3600)

calc = CosmoCalc(args.H_0, args.O_M, args.O_V, args.z_low, args.z_high, args.T_CMB, args.bessel)


if not args.p:
    print "\nFor a Universe with H0 = %s, Omega_M = %s, Omega_V = %s, z = %s and T_CMB = %s K:\n" % (calc.H_0, calc.O_M, calc.O_V, z, calc.T_CMB)
    print "It is now %s Gyr since the Big Bang." % \
        convert_to_gy(calc.age_of_universe(0))
    print "The age at redshift z was %s Gyr." % \
        convert_to_gy(calc.age_of_universe(z))
    print "The light travel time was %s Gyr." % \
        convert_to_gy(calc.light_travel_time(z))
    print "The comoving radial distance is %s MPc." % \
        calc.comoving_radial_dist(z)
    vol = calc.comoving_volume(z) / 10**9
    print "The comoving volume within redshift z is %s Gpc^3." % \
        vol
    print "The angular size distance D_A is %s MPc." % \
        calc.angular_diam_dist(z)
    print "The luminosity distance D_L is %s MPc." % \
        calc.luminosity_dist(z)
    print "The number of baryons in the Universe is %s" % \
        calc.num_baryons()

    
    print "#################################"
    print calc.corr_Tb(10, 0.1, 0.2, 0.01, 1)
    #print calc.M_mp(3, 0.01, 1.0, 7, 9)
    print calc.D_C(args.z_low)
    print calc.D_C(args.z_high)
    #print calc.M_scipy(3, 0.01, 1.0, 7, 9)
    #print calc.sphbess_camb(3000,10000)

else:
    profile.run('calc.corr_Tb(10, 0.1, 0.2, 0.01, 1)', filename="Tb.profile")

    #profile.run('calc.M(3, 0.1, 1.0, 7, 9)', filename="interpol.profile")
    #profile.run('calc.M_mp(3, 0.1, 1.0, 7, 9)', filename="mpquad.profile")
    #profile.run('calc.M_scipy(3, 0.1, 1.0, 7, 9)', filename="scipyquad.profile")
