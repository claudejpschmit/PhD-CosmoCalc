from CosmoCalculatorClass import CosmoCalc

class Fisher(object):
    def __init__(self, H_0, O_M, O_V, z_low_integration, z_high_integration, T_CMB):
        
        # setup the working calculator object.
        self.calc = CosmoCalc(H_0, O_M, O_V, z_low_integration, z_high_integration, T_CMB)
        
        #defining default parameters:
        # params[0] = ombh2
        # params[1] = omch2
        # params[2] = omnuh2
        # params[3] = omk
        # params[4] = hubble
        self.params = {"ombh2":0.0226, "omch2":0.112, "omnuh2":0.00064, "omk":0.0, "hubble":H_0}

 ###################################################################################################

    def update_Model(params):
        pass
    
    def compute_Cl():
        pass 

    def derivative():
        pass

