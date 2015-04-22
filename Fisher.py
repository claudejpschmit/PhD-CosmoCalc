
from math import *
from CosmologyCalculatorClass import CosmoCalc

class Fisher(object):
    def __init__(self, params):
        
        # setup the working calculator object.
        self.calc = CosmoCalc(params)
        
        # defining initial parameters.
        self.params = self.calc.fiducial_params
        self.fiducial_params = self.calc.fiducial_params
        # defining parameter variation matrix
        # This determines by how much each parameter will be varied during each step.
        self.var_params = {}
        self.model_params_keys = ["ombh2", "omch2", "omnuh2", "omk", "hubble"]
        for key in self.model_params_keys:
            # varying parameters by 1% of the fiducial value
            if self.params[key] == 0.0:
                self.var_params[key] = 0.0001
            else:
                self.var_params[key] = self.params[key]/100.0

        self.Cl = []
        self.Cl_inv = []

        self.kmin = 0.001
        self.kmax = 5
        ksteps = 100
        kstepsize = (self.kmax - self.kmin)/float(ksteps)
        self.krange = []
        for n in range(0,ksteps+1):
            self.krange.append(self.kmin + n * kstepsize)
        print "Fisher initialized"
    
    def __call__(self, param_key1, param_key2):
        '''
        while (index = something):
                # generate new parameters for the run
            new_params = update_parameters()
                # update the model using the new parameters
            update_model(new_params)
                # use these parameters to find the fisher matrix
            result(index) = compute_Fl_matrix()
        '''
        pass
 ###################################################################################################

    def update_Model(self, new_params):
        # update q's,
        # update update P(k,z)
        
        # I need to have all my functions taking parameters into account.
        # so translate the cosmo parameters into what is used in those equations.
        # then it should be sufficient to initialize cosmocalc only with a list of parameters.
        # or possibly a dictionary.
        # Maybe have a function generate_params(params) which does this.
        self.calc.generate_params(new_params)
        self.calc.update_q()
        self.calc.Pk_update_interpolator(new_params)
        # At this point the calc object should be updated with the new parameters and can be used 
        # to calculate anything in terms of these new parameters.

        return None
    
    def compute_Cl(self, l):
        # update self.Cl matrix
        # for a certain range of k values compute cl's and append to matrix
        self.Cl = []
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            self.Cl.append(row) 
        return None
    
    def compute_Cl_inv(self):
        self.Cl_inv = np.linalg.inv(self.Cl)
        return None
    
    # This calculates a single matrix element of the derivative matrix.
    def Cl_derivative(self, l, param_key, k1, k2):
        
        h = self.var_params[param_key]
        print ("calculation with h = ",h)
        # storing the original parameter
        x = self.params[param_key]

        self.params[param_key] = x + 2*h
        self.update_Model(self.params)
        f1 = self.calc.Cl(l, k1, k2, self.kmin, self.kmax)           
        
        self.params[param_key] = x + h
        self.update_Model(self.params)
        f2 = self.calc.Cl(l, k1, k2, self.kmin, self.kmax)

        self.params[param_key] = x - h
        self.update_Model(self.params)
        f3 = self.calc.Cl(l, k1, k2, self.kmin, self.kmax)

        self.params[param_key] = x - 2*h
        self.update_Model(self.params)
        f4 = self.calc.Cl(l, k1, k2, self.kmin, self.kmax)

        # reset the model to what it was before. This may not be necessary depending
        # on how the program will later operate.
        self.params[param_key] = x
        #self.update_Model(self.params)

        # actually calculating the derivative matrix
        num = -f1 + 8 * f2 - 8 * f3 + f4
        result =  num / (12.0 * h)

        return result

    # This calculates a single matrix element of the derivative matrix.
    def Cl_loglog_derivative(self, l, param_key, k1, k2):
        
        h = self.var_params[param_key]
        print ("calculation with h = ",h)
        # storing the original parameter
        x = self.params[param_key]

        self.params[param_key] = x + 2*h
        self.update_Model(self.params)
        f1 = log(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))           
        
        self.params[param_key] = x + h
        self.update_Model(self.params)
        f2 = log(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))

        self.params[param_key] = x - h
        self.update_Model(self.params)
        f3 = log(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))

        self.params[param_key] = x - 2*h
        self.update_Model(self.params)
        f4 = log(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))

        # reset the model to what it was before. This may not be necessary depending
        # on how the program will later operate.
        self.params[param_key] = x
        #self.update_Model(self.params)

        # actually calculating the derivative matrix
        num = -f1 + 8 * f2 - 8 * f3 + f4
        result =  num / (12.0 * h)

        return x*result



    def Cl_derivative_matrix(self, l, param_key):
        
        h = self.var_params[param_key]
        # storing the original parameter
        x = self.params[param_key]

        # here we're generating f(x+2h) etc matrices for the various kvalues.
        f1matrix = []
        f2matrix = []
        f3matrix = []
        f4matrix = []

        self.params[param_key] = x + 2*h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f1matrix.append(row)
        
        self.params[param_key] = x + h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f2matrix.append(row)
        
        self.params[param_key] = x - h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f3matrix.append(row)
        
        self.params[param_key] = x - 2*h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f4matrix.append(row)
        
        # reset the parameter to what it was before. 
        # Updating the model may not be necessary as this is done at the beginning
        # of this routine anyways.
        # Edit: Updating the model to a past model should be very cheap now, so
        #       one might as well just do it.
        self.params[param_key] = x
        self.update_Model(self.params)

        # actually calculating the derivative matrix
        res = []
        for k1 in self.krange:
            row = []
            for k2 in self.krange:        
                # Five-point stensil method      
                num = -f1matrix[k1][k2] + 8 * f2matrix[k1][k2] - 8 * f3matrix[k1][k2] + f4matrix[k1][k2]
                result =  num / (12.0 * h)

                row.append(result)
                
            res.append(row)

        return res

    # computes matrix element for Fl wrt 2 parameter names.
    def compute_Fl(self, l, param_key1, param_key2):
        Cl_alpha = self.Cl_derivative_matrix(l,param_key1)
        if param_key1 == param_key2:
            Cl_beta = Cl_alpha
        else:
            Cl_beta = self.Cl_derivative_matrix(l,param_key2)
        
        self.compute_Cl(l)
        self.compute_Cl_inv()
        
        product = np.dot(Cl_alpha, self.Cl_inv)
        product = np.dot(product, Cl_beta)
        product = np.dot(product, self.Cl_inv)

        return 0.5 * np.trace(product)
        
    # this function returns a Fisher matrix element. F(param1, param2)
    def F(self, param_key1, param_key2):
        lmax = 1
        sum = 0
        for l in range(0,lmax + 1):
            sum += (2 * l + 1) * self.compute_Fl(l, param_key1, param_key2)
        return sum

