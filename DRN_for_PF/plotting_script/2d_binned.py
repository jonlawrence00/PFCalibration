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


rawE = ak.sum(RechitEn_pkl, axis=1)
print("==================================")
print("Prediction ratio = ",preds_ratio[0])
print("rawE/trueE = ",trueEn_ratio[0])
print("rawE = ",rawE[0])
print("eta = ",eta_pkl[0])
print("==================================")
print(preds_ratio[0],trueEn_ratio[0],rawE[0])
#preds_trueEn =rawE*(1/preds_ratio) 
#trueEn_pkl = rawE*(1/trueEn_ratio)
preds_trueEn = preds_ratio#rawE*(1/preds_ratio)
trueEn_pkl = trueEn_ratio#rawE*(1/trueEn_ratio)
print("Length of prediction array : ",len(preds_trueEn))
print("Length of true Energy array : ",len(trueEn_pkl))
print("Length of eta array : ",len(eta_pkl))


valid_idx_file ="%s/all_valididx.pickle"%inpickle_folder
train_idx_file ="%s/all_trainidx.pickle"%inpickle_folder
valid_idx_f = open(valid_idx_file,"rb")
valid_idx = np.asarray(pickle.load(valid_idx_f))
print("Length of validation array : ",len(valid_idx))

train_idx_f = open(train_idx_file,"rb")
train_idx = np.asarray(pickle.load(train_idx_f))
print("Length of trained array : ",len(train_idx))

print("Predicted energy (trained) = ",preds_trueEn[train_idx[1]])
print("True energy (trained) = ",trueEn_pkl[train_idx[1]])
print("Predicted energy (validated) = ",preds_trueEn[valid_idx[1]])
print("True energy (validated) = ",trueEn_pkl[valid_idx[1]])

#bin_range=[0,0.05,0.10.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.75]
#bin_range = np.arange(0,3,0.5)
#bin_range =[0,1.55,2.5,2.75]
bin_range =[0,1.55]
print(bin_range)
print("Total bins : ", len(bin_range))
#Ebin_range = np.arange(0,204,4)
#Ebin_range =[0, 2, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 114, 124, 134, 144, 154, 164, 174, 184, 194, 200]
#Ebin_range2 =[2, 8, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104,108, 124, 134, 144, 154, 164, 174, 184, 194, 204, 214, 224]
Ebin_range =[100,104]
Ebin_range2 =[108,124]
print(Ebin_range)
print("Total bins : ", len(Ebin_range))

import ROOT, array

fout= ROOT.TFile(out_fname, 'RECREATE')# ROOT.TFile("./aggre2_maxlr6e04_75ep/constlr_4feat/ep100/hist_85bins_updateRelwt_ratio_Constlr_maxlr1e04_2aggr_ep100_4feat_NSM.root", 'RECREATE')
#import ROOT
hist_norm_predTrue_Valid=[]
hist_norm_predTrue_Train=[]
Energy=[20,50,80,100,120,200,250,300]
M=35 # number of histograms
#Train_trueE_eta_response=ROOT.TH2F(
hist2D_norm_Eta_predTrue_Train=ROOT.TProfile2D('hist2D_norm_Eta_predTrue_Train', """:"(pred-true)/true in GeV":""", 5, 0,2.5,20,0,200)
hist2D_norm_Eta_predTrue_Valid=ROOT.TProfile2D('hist2D_norm_Eta_predTrue_Valid', """:"(pred-true)/true in GeV":""", 5, 0,2.5,20,0,200)
for i_hist in range(len(bin_range)-1):
#    print(i_hist," : ",bin_range[i_hist])
    xhigh_norm= 5
    xlow_norm= -5
    tmp=[]
    tmp2=[]
    for e_hist in range(len(Ebin_range)-1):
        #print(i_hist," : ",bin_range[i_hist]," :::: ",e_hist," : ",Ebin_range[e_hist])
        name1='eta_%f_%f_E_%i_%i' %(bin_range[i_hist],bin_range[i_hist+1],Ebin_range[e_hist],Ebin_range[e_hist+1])#,u[i_hist],v[i_hist],typee[i_hist])
        tmp.append(ROOT.TH1F('Valid_norm_pred_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))       
        tmp2.append(ROOT.TH1F('Train_norm_pred_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))  
    hist_norm_predTrue_Valid.append(tmp)
    hist_norm_predTrue_Train.append(tmp2)
    
    
#hist_norm_predTrue_Valid[i_hist].append(ROOT.TH1F('Valid_norm_pred_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))
        #hist_norm_predTrue_Train[i_hist].append(ROOT.TH1F('Train_norm_pred_%s' % name1, """:"(pred-true)/true in GeV":""", 500, xlow_norm, xhigh_norm))

        
valid_predEn_all=[]
valid_trueEn_all=[]

count=0 #array.array("d", [0]*(len(bin_range)-1))
i=0
ibin=0
for i in range(len(valid_idx)):
    if(ecal_pkl[valid_idx[i]]>0):
        valid_trueEn=(trueEn_pkl[valid_idx[i]])
        valid_predEn=(preds_trueEn[valid_idx[i]])
        trueEn=np.empty(len(bin_range)-1,dtype='float')
        predEn=np.empty(len(bin_range)-1,dtype='float')
        for ibin in range(len(bin_range)-1):
            if(abs(eta_pkl[valid_idx[i]])>=bin_range[ibin] and abs(eta_pkl[valid_idx[i]]) <=bin_range[ibin+1]):
                for Ebin in range(len(Ebin_range)-1):
                    if(abs(trueEn_pkl[valid_idx[i]])>=Ebin_range[Ebin] and abs(trueEn_pkl[valid_idx[i]]) <=Ebin_range2[Ebin]):
                        diff= valid_trueEn - valid_predEn
                        norm = diff/valid_trueEn
                        hist_norm_predTrue_Valid[ibin][Ebin].Fill(norm)
                        count=count+1
            

print("Total Validate entries in "+str(bin_range[ibin])+" - "+str(bin_range[ibin+1])+" : "+str(count))

train_predEn_all=[]
train_trueEn_all=[]

#print("length of train_trueEn = ", len(trueEn_pkl), " and length of predE = ",len(preds_trueEn))
count=0
for i in range(len(train_idx)):
    if(ecal_pkl[train_idx[i]]>0):
        train_trueEn=(trueEn_pkl[train_idx[i]])
        train_predEn=(preds_trueEn[train_idx[i]])
        trueEn=np.empty(len(bin_range)-1,dtype='float')
        predEn=np.empty(len(bin_range)-1,dtype='float')
        for ibin in range(len(bin_range)-1):
            if(abs(eta_pkl[train_idx[i]])>=bin_range[ibin] and abs(eta_pkl[train_idx[i]]) <=bin_range[ibin+1]):
                for Ebin in range(len(Ebin_range)-1):
                    if(abs(trueEn_pkl[train_idx[i]])>=Ebin_range[Ebin] and abs(trueEn_pkl[train_idx[i]]) <=Ebin_range2[Ebin]):
                        diff= train_trueEn - train_predEn
                        norm = diff/train_trueEn
                        hist_norm_predTrue_Train[ibin][Ebin].Fill(norm)
                        count=count+1
                
print("Total Trained entries in "+str(bin_range[ibin])+" - "+str(bin_range[ibin+1])+" : "+str(count))


fun1 = ROOT.TF1('fun1', 'gaus', -0.35, 0.4 )
fun2 = ROOT.TF1('fun2', 'gaus', -0.35, 0.4 )

for ibin in range(len(bin_range)-1):
    for Ebin in range(len(Ebin_range)-1):
        hist_norm_predTrue_Train[ibin][Ebin].Fit(fun1,"R")
        mean=fun1.GetParameter(1)
        rms=fun1.GetParameter(2)/(1.0 + min(0.0, fun1.GetParameter(1)))
        hist2D_norm_Eta_predTrue_Valid.Fill(bin_range[ibin],Ebin_range[Ebin],mean)
        hist_norm_predTrue_Valid[ibin][Ebin].Fit(fun2,"R")
        mean=fun2.GetParameter(1)
        rms=fun2.GetParameter(2)/(1.0 + min(0.0, fun2.GetParameter(1)))
        hist2D_norm_Eta_predTrue_Train.Fill(bin_range[ibin],Ebin_range[Ebin],mean)

n=len(bin_range)-1
x=array.array("d", [0]*n)
rms_tr= array.array("d", [0]*n)
mean_va= array.array("d", [0]*n)
rms_va= array.array("d", [0]*n)
mean_tr= array.array("d", [0]*n)
Train_norm_pred_trueEn=[]
Train_norm_pred_trueEn_reso=[]
Valid_norm_pred_trueEn=[]
Valid_norm_pred_trueEn_reso=[]

for ibin in range(len(bin_range)-1):
    n=len(Ebin_range)-1
    x=array.array("d", [0]*n)
    rms_tr= array.array("d", [0]*n)
    mean_va= array.array("d", [0]*n)
    rms_va= array.array("d", [0]*n)
    mean_tr= array.array("d", [0]*n)
    for Ebin in range(len(Ebin_range)-1):
        x[Ebin]=Ebin_range[Ebin]
        hist_norm_predTrue_Train[ibin][Ebin].Fit(fun1,"R")
        mean_tr[Ebin]=fun1.GetParameter(1)
        rms_tr[Ebin]=fun1.GetParameter(2)/(1.0 + min(0.0, fun1.GetParameter(1)))
        hist_norm_predTrue_Valid[ibin][Ebin].Fit(fun2,"R")
        mean_va[Ebin]=fun2.GetParameter(1)
        rms_va[Ebin]=fun2.GetParameter(2)/(1.0 + min(0.0, fun2.GetParameter(1)))

    Train_norm_pred_trueEn.append(ROOT.TGraph(len(x),x,mean_tr))
    Train_norm_pred_trueEn_reso.append(ROOT.TGraph(len(x),x,rms_tr))
    Valid_norm_pred_trueEn.append(ROOT.TGraph(len(x),x,mean_va))
    Valid_norm_pred_trueEn_reso.append(ROOT.TGraph(len(x),x,rms_va))

    #name2='response_%f_%f' %(bin_range[ibin],bin_range[ibin+1])
    #name3='reso_%f_%f' %(bin_range[ibin],bin_range[ibin+1])
    #Train_norm_pred_trueEn.SetName('Train_%s'%name2)
    #Train_norm_pred_trueEn_reso.SetName('Train_%s'%name3)
    #Valid_norm_pred_trueEn.SetName('Valid_%s'%name2)
    #Valid_norm_pred_trueEn_reso.SetName('Valid_'%name3)

hist2D_norm_Eta_predTrue_Train.GetZaxis().SetRangeUser(-0.5,0.5)
hist2D_norm_Eta_predTrue_Valid.GetZaxis().SetRangeUser(-0.5,0.5)

fout.cd()
for i in range(len(bin_range)-1):
    name2='response_%f_%f' %(bin_range[i],bin_range[i+1])
    name3='reso_%f_%f' %(bin_range[i],bin_range[i+1])
    Train_norm_pred_trueEn[i].Write('Train_%s'%name2)
    Train_norm_pred_trueEn_reso[i].Write('Train_%s'%name3)
    Valid_norm_pred_trueEn[i].Write('Valid_%s'%name2)
    Valid_norm_pred_trueEn_reso[i].Write('Valid_%s'%name3)
    for j in range(len(Ebin_range)-1):
        hist_norm_predTrue_Valid[i][j].Write()
        hist_norm_predTrue_Train[i][j].Write()

hist2D_norm_Eta_predTrue_Valid.Write();
hist2D_norm_Eta_predTrue_Train.Write();

fout.Close()
