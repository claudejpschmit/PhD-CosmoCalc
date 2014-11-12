from CosmologyCalculatorClass import CosmoCalc

def convert_to_gy(age):
    return age * 10**10 * 3.08568 / (365.25 * 24 * 3600) 
# Setting default parameters
H_0 = 70
O_M = 0.3
z = 3
O_V = 0.7
default = raw_input("Do you wish to use default parameters? (y/n): ")

# User defined parameters
if default.lower() == "n" or default.lower() == "no":
    H_0 = float(raw_input("Insert Hubble Constant H_0: "))
    O_M = float(raw_input("Insert Matter density Omega_M: "))
    z = float(raw_input("Insert redshift z: "))
    O_V = float(raw_input("Insert Vacuum density Omega_V: "))

# Output
calc = CosmoCalc(H_0, O_M, O_V)
print "For a Universe with H0 = %s, Omega_M = %s, Omega_V = %s and z = %s:\n" % (calc.H_0, calc.O_M, calc.O_V, z)
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
