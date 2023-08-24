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
outfileName= sys.argv[2]
out_fname = '%s/%s'%(folder,outfileName)
print("output file: ",out_fname)
inpickle_folder =sys.argv[3]
print("input pickle files are picked from ",inpickle_folder)

pred_v2 ="%s/pred_tb.pickle"%folder
predPickle = open(pred_v2, "rb")
preds_ratio = np.asarray(pickle.load(predPickle))
#print(preds_ratio[preds_ratio>3])
#preds_ratio[preds_ratio>3] = 3
#print(preds_ratio[preds_ratio>3])

#trueEn ="%s/ratioflip_target.pickle"%inpickle_folder
trueEn ="%s/trueE_target.pickle"%inpickle_folder
trueEnPickle = open(trueEn,"rb")
trueEn_ratio = np.asarray(pickle.load(trueEnPickle))
RechitEn ="%s/recHitEn.pickle"%inpickle_folder
RechitEnPickle = open(RechitEn,"rb")
RechitEn_pkl =pickle.load(RechitEnPickle)
Eta="%s/eta.pickle"%inpickle_folder
EtaPickle = open(Eta,"rb")
eta_pkl=np.asarray(pickle.load(EtaPickle))

Ecal="%s/ecal.pickle"%inpickle_folder
ecalPickle = open(Ecal,"rb")
ecal_pkl=np.asarray(pickle.load(ecalPickle))
# hit_z ="/home/rusack/shared/pickles/HGCAL_TestBeam/Test_alps/2To3M/Hit_Z.pickle"
# hit_zPickle = open(hit_z,"rb")
# z =pickle.load(hit_zPickle)
# frac =  ((z<54)*0.0105) + (np.logical_and(z>54, z<154)*0.0789) + ((z>154)*0.0316)

rawE = ak.sum(RechitEn_pkl, axis=1)
print("==================================")
print("Prediction ratio = ",preds_ratio[0])
print("rawE/trueE = ",trueEn_ratio[0])
print("rawE = ",rawE[0])
print("eta = ",eta_pkl[0])
print("==================================")
print(preds_ratio[0],trueEn_ratio[0],rawE[0])

########## target: true E / raw E                                                                                                                                          
preds_trueEn =rawE*(preds_ratio)
trueEn_pkl = rawE*(trueEn_ratio)
########## target: raw E / true E                                                                                                                                          
#preds_trueEn = rawE*(1/preds_ratio)                                                                                                                                       
#trueEn_pkl = rawE*(1/trueEn_ratio)                                                                                                                                        
########## target: true E                                                                                                                                                  
#preds_trueEn = preds_ratio                                                                                                                                                
#trueEn_pkl=trueEn_ratio            

print("Length of prediction array : ",len(preds_trueEn))
print("Length of true Energy array : ",len(trueEn_pkl))
print("Length of eta array : ",len(eta_pkl))

valid_idx_file ="%s/all_valididx.pickle"%inpickle_folder
train_idx_file ="%s/all_trainidx.pickle"%inpickle_folder
valid_idx_f = open(valid_idx_file,"rb")
valid_idx = np.asarray(pickle.load(valid_idx_f))
print("Length of validation array : ",len(valid_idx))
print("Length of EH validation array : ",len(valid_idx[ecal_pkl[valid_idx]>0]))

train_idx_f = open(train_idx_file,"rb")
train_idx = np.asarray(pickle.load(train_idx_f))
print("Length of trained array : ",len(train_idx))
print("Length of EH trained array : ",len(train_idx[ecal_pkl[train_idx]>0]))

print("Predicted energy (trained) = ",preds_trueEn[train_idx[1]])
print("True energy (trained) = ",trueEn_pkl[train_idx[1]])
print("Predicted energy (validated) = ",preds_trueEn[valid_idx[1]])
print("True energy (validated) = ",trueEn_pkl[valid_idx[1]])

#bin_range=[0,0.05,0.10.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.75]
bin_range = np.arange(0.00,3.0,0.04)
print(bin_range)
print("Total bins : ", len(bin_range))
Ebin_range=[0, 2, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 114, 124, 134, 144, 154, 164, 174, 184, 194, 200]
print(Ebin_range)
print("Total bins : ", len(Ebin_range))

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
#Train_trueE_eta_response=ROOT.TH2F(
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
    xhigh_norm= 3
    xlow_norm= -3
    
    name1='eta_%f_to_%f' %(bin_range[i_hist],bin_range[i_hist+1])#,u[i_hist],v[i_hist],typee[i_hist])
    hist_pred_Valid.append(ROOT.TH1F('Valid_Predi_%s' % name1, """:"Predicted energy in GeV":""", 500, 0, xhigh_pred))
    hist_true_Valid.append(ROOT.TH1F('Valid_trueEn_%s' % name1, """:"true Beam energy in GeV":""", 500, 0,xhigh_true ))
    hist_pred_Train.append(ROOT.TH1F('Train_Predi_%s' % name1, """:"Predicted energy in GeV":""", 500, 0, xhigh_pred))
    hist_true_Train.append(ROOT.TH1F('Train_trueEn_%s' % name1, """:"true Beam energy in GeV":""", 500, 0, xhigh_true))
    hist_predTrue_Valid.append(ROOT.TH1F('Valid_Diff_Predi_%s' % name1, """:"Predicted -true in GeV":""", 500, xlow_diff, xhigh_diff))
    hist_norm_predTrue_Valid.append(ROOT.TH1F('Valid_norm_pred_trueEn_%s' % name1, """:"(pred-true)/true in GeV":""", 450, xlow_norm, xhigh_norm))
    hist_predTrue_Train.append(ROOT.TH1F('Train_Diff_Predi_%s' % name1, """:"Predicted -true in GeV":""", 500, xlow_diff, xhigh_diff))
    hist_norm_predTrue_Train.append(ROOT.TH1F('Train_norm_pred_trueEn_%s' % name1, """:"(pred-true)/true in GeV":""", 450, xlow_norm, xhigh_norm))
#    for e_hist in range(len(Ebin_range)-1):
        
valid_predEn_all=[]
valid_trueEn_all=[]
# for ibin in range(len(bin_range)-1):
#     #print(ibin)
#     if(ibin==0):
#         bin_range[0]=9.0
    #print(ibin, bin_range[ibin])

count= array.array("d", [0]*(len(bin_range)-1))
i=0
ibin=0
for i in range(len(valid_idx)):
    if(ecal_pkl[valid_idx[i]]==0):
        valid_trueEn=(trueEn_pkl[valid_idx[i]])
        valid_predEn=(preds_trueEn[valid_idx[i]])
        trueEn=np.empty(len(bin_range)-1,dtype='float')
        predEn=np.empty(len(bin_range)-1,dtype='float')
        
        for ibin in range(len(bin_range)-1):
            if(abs(eta_pkl[valid_idx[i]])>=bin_range[ibin] and abs(eta_pkl[valid_idx[i]]) <=bin_range[ibin+1]):
#            if(bin_range[ibin]==2):
                count[ibin]=count[ibin]+1
                diff= valid_trueEn - valid_predEn
                norm = diff/valid_trueEn
                hist_pred_Valid[ibin].Fill(valid_predEn)
                hist_true_Valid[ibin].Fill(valid_trueEn)
                hist_predTrue_Valid[ibin].Fill(diff)
                hist_norm_predTrue_Valid[ibin].Fill(norm)

                    
summ=0
for ibin in range(len(bin_range)-1):
    print("Total events in ",bin_range[ibin],"-",bin_range[ibin+1]," GeV is ",count[ibin])
    summ=summ+count[ibin]
print("Total events in Validation = ",summ)

train_predEn_all=[]
train_trueEn_all=[]
count2= array.array("d", [0]*(len(bin_range)-1))
print("length of train_trueEn = ", len(trueEn_pkl), " and length of predE = ",len(preds_trueEn))
i=0
ibin=0
for i in range(len(train_idx)):
    if(ecal_pkl[train_idx[i]]==0):
        train_trueEn=(trueEn_pkl[train_idx[i]])
        train_predEn=(preds_trueEn[train_idx[i]])
        trueEn=np.empty(len(bin_range)-1,dtype='float')
        predEn=np.empty(len(bin_range)-1,dtype='float')
        for ibin in range(len(bin_range)-1):
            if(abs(eta_pkl[train_idx[i]])>=bin_range[ibin] and abs(eta_pkl[train_idx[i]]) <=bin_range[ibin+1]):
                diff= train_trueEn - train_predEn
                norm = diff/train_trueEn
                count2[ibin]=count2[ibin]+1
                hist_pred_Train[ibin].Fill(train_predEn)
                hist_true_Train[ibin].Fill(train_trueEn)
                hist_predTrue_Train[ibin].Fill(diff)
                hist_norm_predTrue_Train[ibin].Fill(norm)
            

summ=0
ibin=0
for ibin in range(len(bin_range)-1):
    print("Total events in eta : ",bin_range[ibin],"-",bin_range[ibin+1]," is ",count2[ibin])
    summ=summ+count2[ibin]
print("Total events in Training = ",summ)

fun1 = ROOT.TF1('fun1', 'gaus', -0.35, 0.4 )
fun2 = ROOT.TF1('fun2', 'gaus', -0.35, 0.4 )

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
