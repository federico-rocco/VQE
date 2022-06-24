# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:39:07 2022

@author: cosmo
"""
import numpy as np

class Minimizer:
    def __init__(self, 
                method,
                max_iter=200, # Minimizer iterations.
                max_eval=200,  # Funtion evaluations.
                tol=1e-08,
                disp=True,
                adapt=False):
                            
        #super().__init__()
        self.max_iter = max_iter
        self.max_eval = max_eval
        self.tol = tol
        self.method = method
        self.disp = disp
        self.gradient = None

        # Set scipy.minimize options
        self.bounds = None
        max_eval_str = 'maxfev' # Different string notation for different methods
        #if self.method == 'L-BFGS-B':
            #self.bounds = [(0,2*np.pi) for i in theta]
            #max_eval_str = 'maxfun'
        self.options = {'disp':self.disp,
                        'maxiter': self.max_iter,
                        max_eval_str: self.max_eval}
        if method.lower() == 'cobyla':
            del self.options[max_eval_str]
        if method.lower() == 'nelder-mead' and adapt == True:
            self.options['adaptive'] = True



    def __call__(self,theta):
        self.bounds = [(0,2*np.pi) for i in theta]
        from scipy.optimize import minimize, differential_evolution
        result = minimize(self.loss_function,
                          theta,
                          #bounds=self.bounds,
                          method=self.method,
                          options=self.options,
                          jac=self.gradient,
                          tol=self.tol)
        params = result.x
        return params


    def set_loss_function(self,loss_function):
        self.loss_function = loss_function

