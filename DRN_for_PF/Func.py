#!/usr/bin/env python
# coding: utf-8

import math
import numpy as np
import pandas as pd
from sklearn import preprocessing, neighbors
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import pickle
import matplotlib.pyplot as plt
import os
import scipy.optimize as opt
import abc

class Func(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def pdf(self):
        pass

    @abc.abstractmethod
    def mean(self):
        pass

    @abc.abstractmethod
    def resolution(self):
        pass

    @abc.abstractmethod
    def popt_guesses(self):
        pass

    @classmethod
    def fit(self, bin_centers, bin_heights, prev_guess=None):
        try:
            p0 = self.popt_guesses(bin_centers, bin_heights, prev_guess)
            popt, pcov = opt.curve_fit(self.pdf, bin_centers, bin_heights, p0=p0, bounds = self.bounds, sigma = np.sqrt(bin_heights))
            perr = np.sqrt(np.diag(pcov))
            
            mean, meanerr = self.mean(popt, perr)
            res, reserr = self.resolution(popt, perr)
            
            plotx = np.linspace(min(bin_centers), max(bin_centers), 100)
            plt.plot(plotx, self.pdf(plotx, *popt))
            
            chisq = self.chisquare(bin_centers, bin_heights, popt, self.nparam)
            
            return mean, meanerr, res, reserr, chisq
        except:
            return -1, -1, -1, -1, -1
    @classmethod
    def chisquare(self, bin_centers, bin_heights, popt, nparam):
         #compute chisq
        chisq = 0
        for i in range(len(bin_centers)):
            x = bin_centers[i]
            y = bin_heights[i]
            fx = self.pdf(x,*popt)
            r = fx-y
            chisq+=(r*r)/fx
        return chisq[0]/( len(bin_centers) -nparam - 1)
