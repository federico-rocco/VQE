# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:39:07 2022

@author: cosmo
"""
import numpy as np
from qiskit.algorithms.optimizers import SPSA, COBYLA, CG, SLSQP 

class Minimizer:
    def __init__(self, 
                method,
                max_iter=200, # Minimizer iterations.
                tol=1e-08,
                disp=True):

        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.disp = disp
        self.gradient = None
        self.scipy = False

        if self.method == 'spsa':
            self.opt = SPSA(maxiter=self.max_iter)
        elif self.method == 'cobyla':
            self.opt = COBYLA(maxiter=self.max_iter, tol=self.tol, disp=self.disp)
        elif self.method == 'cg':
            self.opt = CG(maxiter=self.max_iter, tol=self.tol, disp=self.disp)
        elif self.method == 'slsqp':
            self.opt = SLSQP(maxiter=self.max_iter, tol=self.tol, disp=self.disp)
        else:
            self.scipy = True
        



    def __call__(self, loss_function, theta):               
            
        if self.scipy == True:
            
            from scipy.optimize import differential_evolution
            self.bounds = [(0,2*np.pi) for i in theta]
            result = differential_evolution(loss_function,
                              bounds=self.bounds,
                              maxiter=self.max_iter,
                              tol=self.tol,
                              disp=self.disp)
                
        else:
            result = self.opt.minimize(loss_function, theta)
            
        params = result.x
        return params

