#######################################################################
# This file contains a writing class for the Cosmology Calcuator,     #
# which writes data calculated by the calculator to a file.           #
#######################################################################

from CosmologyCalculatorClass import CosmoCalc
import numpy as np

class CosmoWrite(CosmoCalc):
    ############################ Writing ############################
 
    def calculate_Ml(self, l, k_fixed, k2_low, k2_high, z_low, z_high, step):
        filename = 'M_'+str(l)+'('+str(k_fixed)+',['+str(k2_low)+','+str(k2_high)+']).'+str(step)+'steps.dat'
        stepsize = (k2_high - k2_low) / float(step)
        x = [k2_low + float(i) * stepsize for i in range(0, step)]
        g = open(filename, 'w')
        g.close()
        
        for k2 in x:
            res = str(self.M(l, k_fixed, k2, z_low, z_high))
            f = open(filename, 'a')
            f.write(str(k2)+" "+res+"\n")
            f.close()   
