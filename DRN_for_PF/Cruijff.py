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

from Func import Func

class Cruijff(Func):
    nparam = 6
    bounds = [(0.0, -np.inf,0.0,0.0,0.0,0.0), (np.inf,np.inf,np.inf,np.inf,np.inf,np.inf)]
    
    @classmethod
    def pdf(self, x, A, mean, sigmaL, alphaL, sigmaR, alphaR):
        x = np.asarray(x)
        if x.ndim==0:
            x  = np.array([x])

        answer = np.zeros(x.shape)

        left = x < mean
        right = np.logical_not(left)

        top = np.square(x-mean)

        bottomL = 2* (sigmaL*sigmaL + alphaL*top[left])
        bottomR = 2* (sigmaR*sigmaR + alphaR*top[right])

        answer[left] = A * np.exp(-np.divide(top[left],bottomL))
        answer[right] = A * np.exp(-np.divide(top[right],bottomR))
        return answer

    @classmethod
    def mean(self, popt,perr):
        return popt[1], perr[1]

    @classmethod
    def resolution(self, popt,perr):
        res = (popt[2] + popt[4])/(2 * popt[1])

        #error propegation
        dsl = perr[2]/popt[2]
        dsr = perr[4]/popt[4]
        dm = perr[1]/popt[1]
        dres = res*np.sqrt(dsl*dsl + dsr*dsr + dm*dm)

        return res, dres

    @classmethod
    def popt_guesses(self, bin_centers, bin_heights, prev_guess):
        mean_guess = bin_centers[np.argmax(bin_heights)]
        sigma_guess = np.sqrt(np.average((bin_centers-mean_guess)**2, weights=bin_heights))
        sigmaL_guess = sigma_guess
        sigmaR_guess = sigma_guess
        if prev_guess is None:
            alphaL_guess = 7e-2
            alphaR_guess = 7e-2
        else:
            alphaL_guess = prev_guess[3]
            alphaR_guess = prev_guess[5]

        A_guess = np.max(bin_heights) / self.pdf(0,1,0,sigmaL_guess, alphaL_guess, sigmaR_guess, alphaR_guess)
        A_guess = A_guess[0]
        p0 = [A_guess, mean_guess, sigmaL_guess, alphaL_guess, sigmaR_guess, alphaR_guess]
        return p0   
