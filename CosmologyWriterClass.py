#######################################################################
# This file contains a writing class for the Cosmology Calcuator,     #
# which writes data calculated by the calculator to a file.           #
#######################################################################

from CosmologyCalculatorClass import CosmoCalc
import numpy as np
import os

class CosmoWrite(CosmoCalc):
    ############################ Writing ############################
    def calculate_Ml(self, l, k_fixed, k2_low, k2_high, step, stepsize = 0):

        if stepsize == 0:
            stepsize = (k2_high - k2_low) / float(step)
        else:
            step = int((k2_high - k2_low) / stepsize)

        
        dir_path = 'output/'
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        
        filename = 'M_'+str(int(l))+'('+str(k_fixed)+',['+str(k2_low)+','+str(k2_high)+']).'+str(step)+'steps.dat'
        x = [k2_low + float(i) * stepsize for i in range(0, step)]
        g = open(os.path.join(dir_path, filename), 'w')
        g.close()
        
        for k2 in x:
            res = str(self.M(l, k_fixed, k2))
            f = open(os.path.join(dir_path, filename), 'a')
            f.write(str(k2)+" "+res+"\n")
            f.close()

    def calculate_Nl(self, l, k_fixed, k2_low, k2_high, step, stepsize = 0):

        if stepsize == 0:
            stepsize = (k2_high - k2_low) / float(step)
        else:
            step = int((k2_high - k2_low) / stepsize)

        
        dir_path = 'output/'
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        
        filename = 'N_'+str(int(l))+'('+str(k_fixed)+',['+str(k2_low)+','+str(k2_high)+']).'+str(step)+'steps.dat'
        x = [k2_low + float(i) * stepsize for i in range(0, step)]
        g = open(os.path.join(dir_path, filename), 'w')
        g.close()
        
        for k2 in x:
            res = str(self.N_bar(l, k_fixed, k2))
            f = open(os.path.join(dir_path, filename), 'a')
            f.write(str(k2)+" "+res+"\n")
            f.close()
