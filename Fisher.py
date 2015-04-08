from CosmologyCalculatorClass import CosmoCalc

class Fisher(object):
    def __init__(self, params, z_low_integration, z_high_integration, T_CMB):
        
        # setup the working calculator object.
        self.calc = CosmoCalc(params, z_low_integration, z_high_integration, T_CMB)
        
        # defining initial parameters.
        self.params = params
        self.fiducial_params = params
        # defining parameter variation matrix
        # This determines by how much each parameter will be varied during each step.
        self.var_params = {}
        for key in params:
            # varying parameters by 1% of the fiducial value
            self.var_params[key] = params[key]/100.0

        self.Cl = []
        self.Cl_inv = []

        self.kmin = 1
        self.kmax = 10
        ksteps = 1
        kstepsize = (self.kmax - self.kmin)/ksteps
        self.krange = []
        for n in range(0,ksteps+1):
            self.krange.append(self.kmin + n * kstepsize)
        print("K range", self.krange)
        print "Fisher initialized"
    def __call__(self, param1, param2):
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

    def Cl_derivative_matrix(self, l, param_index):
        
        h = self.var_params[param_index]
        # storing the original parameter
        x = self.params[param_index]

        # here we're generating f(x+2h) etc matrices for the various kvalues.
        f1matrix = []
        f2matrix = []
        f3matrix = []
        f4matrix = []

        self.params[param_index] = x + 2*h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f1matrix.append(row)
        
        self.params[param_index] = x + h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f2matrix.append(row)
        
        self.params[param_index] = x - h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f3matrix.append(row)
        
        self.params[param_index] = x - 2*h
        self.update_Model(self.params)
        for k1 in self.krange:
            row = []
            for k2 in self.krange:
                row.append(self.calc.Cl(l, k1, k2, self.kmin, self.kmax))
            f4matrix.append(row)
        
        # reset the model to what it was before. This may not be necessary depending
        # on how the program will later operate.
        self.params[param_index] = x
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
    def compute_Fl(self, l, param_index1, param_index2):
        Cl_alpha = self.Cl_derivative_matrix(l,param_index1)
        Cl_beta = self.Cl_derivative_matrix(l,param_index2)

        # The inverse has to be calculated only once!! Put it somewhere else.
        #Cl_inverse = np.inv(self.Cl_matrix)
        
        product = np.dot(Cl_alpha, self.Cl_inverse)
        product = np.dot(product, Cl_beta)
        product = np.dot(product, self.Cl_inverse)

        return 0.5 * np.trace(product)
        
    # this function returns a Fisher matrix element. F(param1, param2)
    def F(self, param1, param2):
        lmax = 1
        sum = 0
        for l in range(0,lmax + 1):
            self.compute_Cl(l)
            self.compute_Cl_inv(l)
            sum += (2 * l + 1) * self.compute_Fl(l, param1, param2)
        return sum

