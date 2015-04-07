from CosmoCalculatorClass import CosmoCalc

class Fisher(object):
    def __init__(self, params, z_low_integration, z_high_integration, T_CMB):
        
        # setup the working calculator object.
        self.calc = CosmoCalc(params, z_low_integration, z_high_integration, T_CMB)
        
        #defining initial parameters.
        self.params = params
        self.fiducial_params = params

    def __call__():
        pass
 ###################################################################################################

    def update_Model(params):
        # update q's,
        # update update P(k,z)
        
        # I need to have all my functions taking parameters into account.
        # so translate the cosmo parameters into what is used in those equations.
        # then it should be sufficient to initialize cosmocalc only with a list of parameters.
        # or possibly a dictionary.
        # Maybe have a function generate_params(params) which does this.
        pass
    
    def compute_Cl(self, l, k1, k2):
        pass 

    # computes matrix element for Fl.
    def compute_Fl(self, l, param1, param2):
        pass
    
    # this function returns a Fisher matrix element. F(param1, param2)
    def F(self, param1, param2):
        lmax = 10
        sum = 0
        for l in range(0,lmax):
            sum += (2 * l + 1) * self.compute_Fl(l, param1, param2)
        return sum

    # Five-point stensil method
    # TODO: change what the derivative is with respect to.
    def derivative(func, parameter, stepsize):
        x = parameter
        h = stepsize
        num = -func(x + 2 * h) + 8 * func(x + h) - 8 * func(x - h) + func(x - 2 * h)
        return num / float(12 * h)
