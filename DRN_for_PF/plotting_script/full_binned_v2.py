#import pandas as pd #dataframes etc                                                                                                       
#import matplotlib.pyplot as plt #plotting                                                                                                 
import pickle
import awkward as ak
import numpy as np
import os, sys
#import seaborn as sns
#import pickle
#import numpy as np

folder = sys.argv[1]
print("input predictions: ",folder)
folder2 = sys.argv[2]
print("input predictions: ",folder2)
outfileName= sys.argv[3]
out_fname = '%s/%s'%(folder,outfileName)
print("output file: ",out_fname)
inpickle_folder =sys.argv[4]
print("input pickle files are picked from ",inpickle_folder)
inpickle2_folder =sys.argv[5]
print("input pickle files are picked from ",inpickle2_folder)

pred_v2 ="%s/pred_tb.pickle"%folder
predPickle = open(pred_v2, "rb")
preds_ratio1= np.asarray(pickle.load(predPickle))
########## other pickle files ############### 
pred_v2 ="%s/pred_tb.pickle"%folder2
predPickle = open(pred_v2, "rb")
preds_ratio2 = np.asarray(pickle.load(predPickle))
preds_ratio=ak.concatenate((preds_ratio1,preds_ratio2), axis=0)
############################################

#print(preds_ratio[preds_ratio>3])
#preds_ratio[preds_ratio>3] = 3
#print(preds_ratio[preds_ratio>3])

#trueEn ="%s/ratioflip_target.pickle"%inpickle_folder
trueEn ="%s/trueE_target.pickle"%inpickle_folder
trueEnPickle = open(trueEn,"rb")
trueEn_ratio1 = np.asarray(pickle.load(trueEnPickle))
RechitEn ="%s/recHitEn.pickle"%inpickle_folder
RechitEnPickle = open(RechitEn,"rb")
RechitEn_pkl1 =pickle.load(RechitEnPickle)

########## other pickle files ###############
trueEn ="%s/trueE_target.pickle"%inpickle2_folder
trueEnPickle = open(trueEn,"rb")
trueEn_ratio2 = np.asarray(pickle.load(trueEnPickle))
RechitEn ="%s/recHitEn.pickle"%inpickle2_folder
RechitEnPickle = open(RechitEn,"rb")
RechitEn_pkl2 =pickle.load(RechitEnPickle)
trueEn_ratio=ak.concatenate((trueEn_ratio1,trueEn_ratio2), axis=0)
#RechitEn_pkl=ak.concatenate((RechitEn_pkl1,RechitEn_pkl2), axis=0)
############################################


# hit_z ="/home/rusack/shared/pickles/HGCAL_TestBeam/Test_alps/2To3M/Hit_Z.pickle"
# hit_zPickle = open(hit_z,"rb")
# z =pickle.load(hit_zPickle)
# frac =  ((z<54)*0.0105) + (np.logical_and(z>54, z<154)*0.0789) + ((z>154)*0.0316)

rawE1 = ak.sum(RechitEn_pkl1, axis=1)
rawE2 = ak.sum(RechitEn_pkl2, axis=1)
rawE=ak.concatenate((rawE1,rawE2), axis=0)

print("==================================")
print("Prediction ratio = ",preds_ratio[0])
print("rawE/trueE = ",trueEn_ratio[0])
print("rawE = ",rawE[0])
print("==================================")
print(preds_ratio[0],trueEn_ratio[0],rawE[0])
#preds_trueEn =rawE*(1/preds_ratio) 
#trueEn_pkl = rawE*(1/trueEn_ratio)
preds_trueEn = preds_ratio#rawE*(1/preds_ratio)
trueEn_pkl = trueEn_ratio#rawE*(1/trueEn_ratio)

valid_idx_file ="%s/all_valididx.pickle"%inpickle_folder
train_idx_file ="%s/all_trainidx.pickle"%inpickle_folder
valid_idx_f = open(valid_idx_file,"rb")
valid_idx1 = np.asarray(pickle.load(valid_idx_f))
print("Length of validation array : ",len(valid_idx1))
train_idx_f = open(train_idx_file,"rb")
train_idx1 = np.asarray(pickle.load(train_idx_f))
print("Length of trained array : ",len(train_idx1))

print("################### idx for other file #########")
valid_idx_file ="%s/all_valididx.pickle"%inpickle2_folder
train_idx_file ="%s/all_trainidx.pickle"%inpickle2_folder
valid_idx_f = open(valid_idx_file,"rb")
valid_idx2 = np.asarray(pickle.load(valid_idx_f))
print("Length of validation array : ",len(valid_idx2))
train_idx_f = open(train_idx_file,"rb")
train_idx2 = np.asarray(pickle.load(train_idx_f))
print("Length of trained array : ",len(train_idx2))
train_idx=ak.concatenate((train_idx1,train_idx2), axis=0)
valid_idx=ak.concatenate((valid_idx1,valid_idx2), axis=0)
print("Added Length of validation array : ",len(valid_idx))
print("Added Length of trained array : ",len(train_idx))
print("################################################")

print("Predicted energy (trained) = ",preds_trueEn[train_idx[1]])
print("True energy (trained) = ",trueEn_pkl[train_idx[1]])
print("Predicted energy (validated) = ",preds_trueEn[valid_idx[1]])
print("True energy (validated) = ",trueEn_pkl[valid_idx[1]])
bin_range =[2, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 114, 124, 134, 144, 154, 164, 174, 184, 194]
#bin_range = np.arange(10,200,4)
print(bin_range)
print("Total bins : ", len(bin_range))

import ROOT, array

fout= ROOT.TFile(out_fname, 'RECREATE')# ROOT.TFile("./aggre2_maxlr6e04_75ep/constlr_4feat/ep100/hist_85bins_updateRelwt_ratio_Constlr_maxlr1e04_2aggr_ep100_4feat_NSM.root", 'RECREATE')
#import ROOT
hist_pred_Valid=[]
hist_true_Valid=[]
hist_pred_Train=[]
hist_true_Train=[]
hist_predTrue_Valid=[]
hist_norm_predTrue_Valid=[]
hist_predTrue_Train=[]
hist_norm_predTrue_Train=[]
hist_pred_Tbdata=[]
hist_true_Tbdata=[]
hist_predTrue_Tbdata=[]
hist_norm_predTrue_Tbdata=[]
Energy=[20,50,80,100,120,200,250,300]
M=35 # number of histograms

#Train_norm_pred_trueEn=ROOT.TGraph('Train_norm_pred_trueEn', '(pred vs true) in GeV', 200, 0, 200)
#Valid_norm_pred_trueEn=ROOT.TGraph('Valid_norm_pred_trueEn', '(pred vs true) in GeV', 200, 0, 200)
#Train_norm_pred_trueEn_reso=ROOT.TGraph('Train_norm_pred_trueEn_reso', '(pred vs true) in GeV', 200, 0, 200)
#Valid_norm_pred_trueEn_reso=ROOT.TGraph('Valid_norm_pred_trueEn_reso', '(pred vs true) in GeV', 200, 0, 200)

for i_hist in range(len(bin_range)-1):
    print(i_hist," : ",bin_range[i_hist])
    if(bin_range[i_hist]<10):
        xhigh_pred = 2.0*bin_range[i_hist]
    elif(bin_range[i_hist]<100):
        xhigh_pred = 4.0*bin_range[i_hist]
    else: #(bin_range[i_hist]>=100 and bin_range[i_hist]<200):):
        xhigh_pred= 10.0*bin_range[i_hist]

    xhigh_true= 2.0*bin_range[i_hist]

    xhigh_diff= 20
    xlow_diff= -20
    xhigh_norm= 5
    xlow_norm= -5
    
    name1='TrueEn_%i_to_%i' %(bin_range[i_hist],bin_range[i_hist+1])#,u[i_hist],v[i_hist],typee[i_hist])
    hist_pred_Valid.append(ROOT.TH1F('Valid_Predi_%s' % name1, """:"Predicted energy in GeV":""", 500, 0, xhigh_pred))
    hist_true_Valid.append(ROOT.TH1F('Valid_trueEn_%s' % name1, """:"true Beam energy in GeV":""", 500, 0,xhigh_true ))
    hist_pred_Train.append(ROOT.TH1F('Train_Predi_%s' % name1, """:"Predicted energy in GeV":""", 500, 0, xhigh_pred))
    hist_true_Train.append(ROOT.TH1F('Train_trueEn_%s' % name1, """:"true Beam energy in GeV":""", 500, 0, xhigh_true))
    hist_predTrue_Valid.append(ROOT.TH1F('Valid_Diff_Predi_%s' % name1, """:"Predicted -true in GeV":""", 500, xlow_diff, xhigh_diff))
    hist_norm_predTrue_Valid.append(ROOT.TH1F('Valid_norm_pred_trueEn_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))
    hist_predTrue_Train.append(ROOT.TH1F('Train_Diff_Predi_%s' % name1, """:"Predicted -true in GeV":""", 500, xlow_diff, xhigh_diff))
    hist_norm_predTrue_Train.append(ROOT.TH1F('Train_norm_pred_trueEn_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))

valid_predEn_all=[]
valid_trueEn_all=[]
# for ibin in range(len(bin_range)-1):
#     #print(ibin)
#     if(ibin==0):
#         bin_range[0]=9.0
    #print(ibin, bin_range[ibin])
for i in range(len(valid_idx)):
    valid_trueEn=(trueEn_pkl[valid_idx[i]])
    valid_predEn=(preds_trueEn[valid_idx[i]])
    trueEn=np.empty(len(bin_range)-1,dtype='float')
    predEn=np.empty(len(bin_range)-1,dtype='float')
    for ibin in range(len(bin_range)-1):
        if(ibin==0):
            inext=5
        else:
            inext=4
        #f(ibin<len(bin_range)-1):
        if(valid_trueEn>=bin_range[ibin] and valid_trueEn <=bin_range[ibin+1]):
            #print(bin_range[ibin],bin_range[ibin+1])
            diff= valid_trueEn - valid_predEn
            norm = diff/valid_trueEn
            hist_pred_Valid[ibin].Fill(valid_predEn)
            hist_true_Valid[ibin].Fill(valid_trueEn)
            hist_predTrue_Valid[ibin].Fill(diff)
            hist_norm_predTrue_Valid[ibin].Fill(norm)
train_predEn_all=[]
train_trueEn_all=[]
# for ibin in range(len(bin_range)-1):
#     #print(ibin)
#     if(ibin==0):
#         bin_range[0]=9.0
    #print(ibin, bin_range[ibin])
for i in range(len(train_idx)):
    train_trueEn=(trueEn_pkl[train_idx[i]])
    train_predEn=(preds_trueEn[train_idx[i]])
    trueEn=np.empty(len(bin_range)-1,dtype='float')
    predEn=np.empty(len(bin_range)-1,dtype='float')
#    diff= train_trueEn - train_predEn
#    norm = diff/train_trueEn
#    Train_norm_pred_trueEn.Fill(train_trueEn,norm)
    for ibin in range(len(bin_range)-1):
        if(ibin==0):
            inext=5
        else:
            inext=4
        #if(ibin<len(bin_range)-1):
        if(train_trueEn>=bin_range[ibin] and train_trueEn <=bin_range[ibin+1]):
            diff= train_trueEn - train_predEn
            norm = diff/train_trueEn
            hist_pred_Train[ibin].Fill(train_predEn)
            hist_true_Train[ibin].Fill(train_trueEn)
            hist_predTrue_Train[ibin].Fill(diff)
            hist_norm_predTrue_Train[ibin].Fill(norm)

fun1 = ROOT.TF1('fun1', 'gaus', -0.5, 0.8 )
fun2 = ROOT.TF1('fun2', 'gaus', -0.5, 0.8 )

n=len(bin_range)-1
x=array.array("d", [0]*n)
rms_tr= array.array("d", [0]*n)
mean_va= array.array("d", [0]*n)
rms_va= array.array("d", [0]*n)
mean_tr= array.array("d", [0]*n)

for ibin in range(len(bin_range)-1):
    x[ibin]=bin_range[ibin]
    hist_norm_predTrue_Train[ibin].Fit(fun1,"R")
    mean_tr[ibin]=fun1.GetParameter(1)
    rms_tr[ibin]=fun1.GetParameter(2)/(1.0 + min(0.0, fun1.GetParameter(1)))
    hist_norm_predTrue_Valid[ibin].Fit(fun2,"R")
    mean_va[ibin]=fun2.GetParameter(1)
    rms_va[ibin]=fun2.GetParameter(2)/(1.0 + min(0.0, fun2.GetParameter(1)))



Train_norm_pred_trueEn=ROOT.TGraph(len(x),x,mean_tr)
Train_norm_pred_trueEn_reso=ROOT.TGraph(len(x),x,rms_tr)
Valid_norm_pred_trueEn=ROOT.TGraph(len(x),x,mean_va)
Valid_norm_pred_trueEn_reso=ROOT.TGraph(len(x),x,rms_va)
fout.cd()
for i in range(len(bin_range)-1):
    hist_pred_Valid[i].Write()
    hist_true_Valid[i].Write()
    hist_pred_Train[i].Write()
    hist_true_Train[i].Write()
    hist_predTrue_Valid[i].Write()
    hist_norm_predTrue_Valid[i].Write()
    hist_predTrue_Train[i].Write()
    hist_norm_predTrue_Train[i].Write()


Train_norm_pred_trueEn.Write("Train_norm_pred_trueEn")
Train_norm_pred_trueEn_reso.Write("Train_norm_pred_trueEn_reso")
Valid_norm_pred_trueEn.Write("Valid_norm_pred_trueEn")
Valid_norm_pred_trueEn_reso.Write("Valid_norm_pred_trueEn_reso")

fout.Close()
