#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from Cruijff import Cruijff

matplotlib.use("AGG")

transforms = {
    'trueE' : lambda pred, rawE, trueE: pred/trueE,
    'ratio' : lambda pred, rawE, trueE : pred*rawE/trueE,
    'logratio' : lambda pred, rawE, trueE : np.exp(pred)*rawE/trueE,
    'ratiolog' : lambda pred, rawE, trueE : np.exp(pred*np.log(rawE))/trueE,
    'ratioflip' : lambda pred, rawE, trueE : rawE/pred/trueE,
    'logratioflip' : lambda pred, rawE, trueE : rawE/np.exp(pred)/trueE,
    'ratiologflip' : lambda pred, rawE, trueE : np.exp(np.log(rawE)/pred)/trueE,
    'logratiolog' : lambda pred, rawE, trueE : transforms['ratiolog'](np.exp(pred), rawE, trueE),
    'logratiologflip' : lambda pred, rawE, trueE : transforms['ratiologflip'](np.exp(pred), rawE, trueE),
    'log' : lambda pred, rawE, trueE : np.exp(pred)/trueE,
    'BDT' : lambda pred, rawE, trueE: pred
}

class Hist:
    @classmethod
    def find_core(self, histdata, p):
        '''
        Find the central p*100% of the data
        Defined to be the minimum interval that contains p% of the data.
        '''
        bottom = 0
        top = bottom+p
        i = 0
        lens = []
        step = 0.005
        while(top<1):
            rangebounds = np.quantile(histdata,[bottom,top])
            lens += [rangebounds[1]-rangebounds[0]]
            i+=1
            bottom = step*i
            top = bottom+p
        min_i = np.argmin(lens)
        bottom = step*min_i
        top = bottom+p
        rangebounds = np.quantile(histdata,[bottom,top])

        #return 0.8, 1.1
        return rangebounds[0], rangebounds[1]

    @classmethod
    def make_histdata(self, y_pred, rawE, trueE, binvar, binmin, binmax, target):
        '''
        Turn y_pred into correctly filters E_pred/E_true values

        @param y_pred: predicted y values
        @param rawE: raw energy
        @param trueE: true energy

        @param binvar: the variable that we are binning w.r.t
        @param binmin: the minimum value of binvar
        @param binmax: the maximum value of binvar

        @param target: the name of the target (eg 'logratioflip')
        '''
        transform = transforms[target]

        histdata = transform(y_pred, rawE, trueE)

        if binvar is not None and binmin is not None and binmax is not None: 
            filtmask = np.logical_and(binvar<binmax, binvar>=binmin)
            print("\tSelected %d events"%np.sum(filtmask))
            center = np.mean(binvar[filtmask])
            histdata = histdata[filtmask]
        else:
            center = None

        return histdata, center

    @classmethod
    def do_histogram(self, y_pred, rawE, trueE, 
            target, 
            prefix,
            binvar = None, binmin = None, binmax = None, binname=None, 
            func=Cruijff, coresize = 0.95, watermark=True):        
        '''
        Do a histogram

        @param y_pred: prediction y values
        @param rawE: raw energy
        @param trueE: true energy

        @param target: the name of the target (eg 'logratioflip')

        @param prefix: the prefix the name all the output files with

        @param binvar: the variable that we are binning w.r.t
        @param binmin: the minimum value of binvar
        @param binmax: the maximum value of binvar
        @param binnname: the name of the binned variable, formatted nicely to put on the plot

        @param func: the function to fit with
        @param coresize: the coresize (0 < coresize <= 1)
        '''
        #clear the plot
        plt.cla()
        plt.clf()
        
        if binvar is not None:
            print("Binning %s from %f to %f"%(binname,binmin,binmax))

        histdata, center = self.make_histdata(y_pred, rawE, trueE, 
                binvar, binmin, binmax, 
                target)

        statmean = np.mean(histdata)
        statvar = np.var(histdata)
        statcount = len(histdata)
        print("\tstats:")
        print('\t\tmean: %.3f'%statmean)
        print('\t\tvariance: %.3f'%statvar)
        print('\t\tcount: %d'%statcount)

        if statcount<8:
            return -1, -1, -1, -1, -1, -1

        rangemin, rangemax = self.find_core(histdata, coresize)
        rangemask = np.logical_and(histdata>rangemin, histdata<rangemax)
        histdata = histdata[rangemask]
        print("\tCore region of %d points"%len(histdata))

        nbins = int(np.cbrt(len(histdata)))
        
        title = ''
        outname = ""
        if binvar is not None and binmin is not None and binmax is not None and binname is not None:
            title+='\n'
            title+=r'$%.3f < $'%binmin
            title+=binname
            title+=r' $<%.3f$'%binmax

            outname += '_%s%.3f-%.3f'%(binname, binmin, binmax)

        bin_heights, bin_borders, _ = plt.hist(histdata, bins=nbins, label='histogram')
        bin_centers = (bin_borders[:-1] + bin_borders[1:]) / 2
        
        mean, meanerr, res, reserr, chisq = func.fit(bin_centers, bin_heights)

        plt.xlabel(r'$E_{pred}/E_{true}$')
        plt.ylabel('counts')

        if watermark:
            plt.text(0.68, 0.90, 'CMS', transform=plt.gca().transAxes, fontsize=16, color='black', fontweight='bold')
            plt.text(0.68, 0.85, 'Work in progress', transform=plt.gca().transAxes, fontsize=12, color='black', style='italic')

        title = r'$\mu=%.3f,\ \sigma/\mu=%.3f\%%,\ \chi^2 = %.3f$' %(mean, res*100,chisq) + title
        plt.title(title)
        plt.grid(True)

        textstr='Stat Mean: %.3f\nStat Var: %.5f\nN: %d'%(statmean, statvar, statcount)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=props)

        histoutname = prefix + '_HIST' + outname + '.png'
        plt.savefig(histoutname, format='png',bbox_inches='tight')

        plt.clf()
        plt.cla()

        data = {'Left_Edge' : bin_borders[:-1],
                'Right_Edge' : bin_borders[1:],
                'Height' : bin_heights}
        df = pd.DataFrame.from_dict(data)

        dataoutname = prefix + '_HISTDATA' + outname + '.csv'
        df.to_csv(dataoutname)
        
        print('\tMean =',mean)
        print('\tRes =',res)

        return mean, meanerr, res, reserr, chisq, center

