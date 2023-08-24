#define skim_cxx
#include "skim.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>

void skim::Loop()
{
  TString s_data="full";
  //  TString s_data="trans";
  //  TString s_data="barrel";
  //    TString s_data="ec_in";
  //   TString s_data="ec_out";
  TString ext="neg";

  //   In a ROOT session, you can do:
//      root> .L skim.C
//      root> skim t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
 
  cout << "nentries " << nentries << endl;
  TFile *f;
  
  //   TFile *f = new TFile("step3_eta2pt5to2pt75_pos.root","RECREATE");
  if(s_data.Contains("barrel"))f = new TFile("step3_eta1pt55.root","RECREATE");
  else if(s_data.Contains("ec_in")) f = new TFile("step3_eta1pt55to2pt5_"+ext+".root","RECREATE");
  else if(s_data.Contains("ec_out")) f = new TFile("step3_eta2pt5to2pt75_"+ext+".root","RECREATE");
  else if(s_data.Contains("trans")) f = new TFile("step3_eta1to1pt5_"+ext+".root","RECREATE");
  //  else f = new TFile("tmp.root","RECREATE");
  else f = new TFile("step3_2_350_eta3.root","RECREATE");
  //  else  f = new TFile("step3_eta2pt75_"+ext+".root","RECREATE");

  TH2D* h_recrawvstrue    = new TH2D("h_recrawvstrue","Rechit raw energy (x axis) wrt True energy (y axis)",800,0,400,800,0,400);
  TH2D* h_recrawvsraw    = new TH2D("h_recrawvsraw","Rechit raw energy (x axis) wrt Raw cluster energy (y axis)",800,0,400,800,0,400);
  TH2D* h_recXvsrecY    = new TH2D("h_recXvsrecY","Rechit (x axis) wrt Rechit (y axis)",10000,-300,300,10000,-300,300);
  TH2D* h_recXvsrecY_oneevt    = new TH2D("h_recXvsrecY_oneevt","Rechit (x axis) wrt Rechit (y axis) ",1200,-300,300,1200,-300,300);
  TH2D* h_recZvsrecY    = new TH2D("h_recZvsrecY","Rechit (z axis) wrt Rechit (y axis)",20000,-600,600,10000,-300,300);
  TH2D* h_recEtavsrecY    = new TH2D("h_recEtavsrecY","Rechit eta wrt Rechit (y axis)",100,-3,3,10000,-300,300);
  TH2D* h_recXvsrecY_ecal    = new TH2D("h_recXvsrecY_ecal","Rechit (x axis) wrt Rechit (y axis)",10000,-300,300,10000,-300,300);
  TH2D* h_recXvsrecY_hcal    = new TH2D("h_recXvsrecY_hcal","Rechit (x axis) wrt Rechit (y axis)",10000,-300,300,10000,-300,300);

  //  TFile *f = new TFile("step3_eta1to1pt5.root","RECREATE");  
  f->cd();
  TTree *newtree = fChain->CloneTree(0);
  float ecalEn=0.0;
  float	hcalEn=0.0;
  int EH=0;
  int H=0;

  fChain->SetBranchStatus("*",0); 

  fChain->SetBranchStatus("HcalRechit_posx",1);
  fChain->SetBranchStatus("HcalRechit_posy",1);
  fChain->SetBranchStatus("HcalRechit_posz",1);
  fChain->SetBranchStatus("HcalRechit_E",1);
  fChain->SetBranchStatus("HcalPFclustereta",1);
  fChain->SetBranchStatus("HcalRechit_depth",1);

  fChain->SetBranchStatus("EcalRechit_posx",1);
  fChain->SetBranchStatus("EcalRechit_posy",1);
  fChain->SetBranchStatus("EcalRechit_posz",1);
  fChain->SetBranchStatus("EcalRechit_E",1);
  fChain->SetBranchStatus("EcalPFclustereta",1);
  fChain->SetBranchStatus("EcalRechit_depth",1);
  fChain->SetBranchStatus("true",1);
  fChain->SetBranchStatus("eta",1);
  fChain->SetBranchStatus("phi",1);
  fChain->SetBranchStatus("ecal",1);
  fChain->SetBranchStatus("hcal",1);
  fChain->SetBranchStatus("Eecal",1);
  fChain->SetBranchStatus("Ehcal",1);
  fChain->SetBranchStatus("dr",1);
  fChain->SetBranchStatus("pfcID",1);
  fChain->SetBranchStatus("p", 1);
  newtree->Branch("ecalEn",&ecalEn,"ecalEn/F");
  newtree->Branch("hcalEn",&hcalEn,"hcalEn/F");

  //s->Branch("ecal",&ecal_,"ecal/F");

  int var=0;
  fChain->SetBranchStatus("HcalRechit_totE",1);
  vector<float> ecalEnergies,  hcalEnergies;
  int oneevt=0;

  Long64_t nbytes = 0, nb = 0, nSurvived=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //   for (Long64_t jentry=0; jentry<10000000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    //      if (abs(eta) >= 1.5 || abs(eta) < 1.0) continue;
    //cout<<jentry<<" : abs(eta) = "<<abs(eta)<<" , X hit = "<<EcalRechit_posx->size()+HcalRechit_posx->size()<<endl;


    float e = 0.0 , h = 0.0;
    if(pfcID->size()> 0) {
      for(unsigned ii = 0; ii < pfcID->size(); ii++) {
	if (pfcID->at(ii) == 5 && dr->at(ii) < 0.4) {
	  e += Eecal->at(ii);
	  h += Ehcal->at(ii);
	}
      }
      ecalEnergies.push_back(e);
      hcalEnergies.push_back(h);
      ecalEn=e;
      hcalEn=h;
    }
    else {
      ecalEnergies.push_back(ecal);
      hcalEnergies.push_back(hcal);
      ecalEn=ecal;
      hcalEn=hcal;
    }

      
    if(s_data.Contains("ec_out")){
      if(ext=="pos"){
	if (!(eta>2.5 && eta<=2.75)) continue;
	if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)>2.5 && HcalPFclustereta->at(0)<= 2.75)) continue;
	if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)>2.5 && EcalPFclustereta->at(0)<= 2.75)) continue;
      }
      else if(ext=="neg"){
	if (!(eta<-2.5 && eta>=-2.75)) continue;
	if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)<-2.5 && HcalPFclustereta->at(0)>=-2.75)) continue;
	if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)<-2.5 && EcalPFclustereta->at(0)>=-2.75)) continue;
      }
    }
    else if(s_data.Contains("ec_in")){
      if(ext=="pos"){
	if (!(eta>1.55 && eta<=2.5)) continue;
	if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)>1.55 && HcalPFclustereta->at(0)<=2.5)) continue;
	if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)>1.55 && EcalPFclustereta->at(0)<=2.5)) continue;
      }
      else if(ext=="neg"){
	if (!(eta<-1.55 && eta>=-2.5)) continue;
	if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)<-1.55 && HcalPFclustereta->at(0)>=-2.5)) continue;
	if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)<-1.55 && EcalPFclustereta->at(0)>=-2.5)) continue;	
      }
    }
    else if(s_data.Contains("trans")){
      if(ext=="pos"){
	if (!(eta>1 && eta<=1.5)) continue;
        if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)>1 && HcalPFclustereta->at(0)<=1.5)) continue;
        if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)>1 && EcalPFclustereta->at(0)<=1.5)) continue;
      }
      else if(ext=="neg"){
        if (!(eta<-1 && eta>=-1.5)) continue;
        if(HcalRechit_posx->size()>0 && !(HcalPFclustereta->at(0)<-1 && HcalPFclustereta->at(0)>=-1.5)) continue;
        if(EcalRechit_posx->size()>0 && !(EcalPFclustereta->at(0)<-1 && EcalPFclustereta->at(0)>=-1.5)) continue;
      }
    }

    else if(s_data.Contains("barrel")){
      if (abs(eta)>1.55) continue;
      if(HcalRechit_posx->size()>0 && abs(HcalPFclustereta->at(0))>1.55) continue;
      if(EcalRechit_posx->size()>0 && abs(EcalPFclustereta->at(0))>1.55) continue;
    }    
    else{
      if (abs(eta)>2.75) continue;
      if(HcalRechit_posx->size()>0 && abs(HcalPFclustereta->at(0))>2.75) continue;
      if(EcalRechit_posx->size()>0 && abs(EcalPFclustereta->at(0))>2.75) continue;
    }

      double TotalHcalenergy=0, TotalEcalenergy=0;
      for(int i=0; i<HcalRechit_posx->size();i++)
	{
	  TotalHcalenergy+=HcalRechit_E->at(i);
	}
      for(int i=0; i<EcalRechit_posx->size();i++)
	{
	  //	  cout<<i<<" : ECAL rechit x axis position : "<<EcalRechit_posx->at(i)<<" : ECAL rechit y axis position : "<<EcalRechit_posy->at(i)<<endl;
	  TotalEcalenergy+=EcalRechit_E->at(i);
	}
      float ecal_=ecalEnergies[jentry];
      float hcal_=hcalEnergies[jentry];

      if((TotalHcalenergy+TotalEcalenergy)==0) {
	cout<<jentry<<" : HCAL rechit x axis position : "<<HcalRechit_posx->size()<<" : HCAL rechit y axis position : "<<HcalRechit_posy->size()<<endl;
	cout<<jentry<<" : ECAL rechit x axis position : "<<EcalRechit_posx->size()<<" : ECAL rechit y axis position : "<<EcalRechit_posy->size()<<endl;       
	continue;
      }


      if(!(true_>4 && true_<350)) continue;
      if(!((ecal_+hcal_)>4 && (ecal_+hcal_)<350)) continue;
      for(int i=0; i<HcalRechit_posx->size();i++)
        {
	  
	  //	  if(EcalRechit_posx->size()==0 && oneevt==0 && eta==0) {
	  //  cout<<jentry<<": (x,y,z,E) == ("<<HcalRechit_posx->at(i)<<", "<<HcalRechit_posy->at(i)<<", "<<HcalRechit_posz->at(i)<<", "<<TotalHcalenergy<<")"<<endl;
	  //	  }
	  
	  if(jentry==1109) h_recXvsrecY_oneevt->Fill(HcalRechit_posx->at(i),HcalRechit_posy->at(i));
          h_recXvsrecY->Fill(HcalRechit_posx->at(i),HcalRechit_posy->at(i));
          h_recXvsrecY_hcal->Fill(HcalRechit_posx->at(i),HcalRechit_posy->at(i));
          h_recZvsrecY->Fill(HcalRechit_posz->at(i),HcalRechit_posy->at(i));
	  h_recEtavsrecY->Fill(HcalPFclustereta->at(i),HcalRechit_posy->at(i));
        }
      if(EcalRechit_posx->size()==0 && HcalRechit_posx->size()>0 && eta==0) oneevt++;
      for(int i=0; i<EcalRechit_posx->size();i++)
        {
          h_recXvsrecY->Fill(EcalRechit_posx->at(i),EcalRechit_posy->at(i));
          h_recXvsrecY_ecal->Fill(EcalRechit_posx->at(i),EcalRechit_posy->at(i));
          h_recZvsrecY->Fill(EcalRechit_posz->at(i),EcalRechit_posy->at(i));
          h_recEtavsrecY->Fill(EcalPFclustereta->at(i),EcalRechit_posy->at(i));
        }
      
      /*      
      if(pfcID->size()>0){
	cout<<jentry<<": Ecal energy : "<<ecal_<<", ECAL rechit energy : "<<TotalEcalenergy<<" , Hcal energy : "<<hcal_<<", HCAL rechit energy : "<<TotalHcalenergy<<endl;
      }
      */
      h_recrawvstrue->Fill(TotalEcalenergy+TotalHcalenergy,true_);      
      h_recrawvsraw->Fill(TotalEcalenergy+TotalHcalenergy,ecal_+hcal_);
      
      if(ecal_>0)
	EH++;
      else
	H++;
      //h_recrawvstrue->Draw("colz");
      nSurvived++;

      newtree->Fill();
  }
  gStyle->SetOptStat(0);
  h_recrawvstrue->GetXaxis()->SetRangeUser(4,350);
  h_recrawvstrue->GetYaxis()->SetRangeUser(4,350);
  h_recrawvsraw->GetXaxis()->SetRangeUser(4,350);
  h_recrawvsraw->GetYaxis()->SetRangeUser(4,350);
  h_recrawvstrue->GetXaxis()->SetTitle("E^{raw}_{rechit} (GeV)");
  h_recrawvstrue->GetYaxis()->SetTitle("E^{true} (GeV)");
  h_recrawvsraw->GetXaxis()->SetTitle("E^{raw}_{rechit} (GeV)");
  h_recrawvsraw->GetYaxis()->SetTitle("E^{raw}_{cluster} (GeV)");

    
  h_recZvsrecY->GetXaxis()->SetTitle("Rechit Z position (cm)");
  h_recXvsrecY->GetXaxis()->SetTitle("Rechit X position (cm)");
  h_recZvsrecY->GetYaxis()->SetTitle("Rechit Y position (cm)");
  h_recXvsrecY->GetYaxis()->SetTitle("Rechit Y position (cm)");
  h_recXvsrecY_ecal->GetYaxis()->SetTitle("ECAL Rechit Y position (cm)");
  h_recXvsrecY_hcal->GetYaxis()->SetTitle("HCAL Rechit Y position (cm)");
  h_recXvsrecY_ecal->GetXaxis()->SetTitle("ECAL Rechit X position (cm)");
  h_recXvsrecY_hcal->GetXaxis()->SetTitle("HCAL Rechit X position (cm)");
  h_recEtavsrecY->GetXaxis()->SetTitle("Rechit eta");
  h_recEtavsrecY->GetYaxis()->SetTitle("Rechit Y position (cm)");
  h_recXvsrecY_oneevt->GetXaxis()->SetTitle("Rechit X position (cm)");
  h_recXvsrecY_oneevt->GetYaxis()->SetTitle("Rechit Y position (cm)");
  h_recXvsrecY_oneevt->GetZaxis()->SetTitle("Rechit Z position (cm)");
    
  newtree->Write("s");
  h_recXvsrecY->Write();
  h_recZvsrecY->Write();
  h_recrawvstrue->Write();
  h_recrawvsraw->Write();
  h_recEtavsrecY->Write();
  h_recXvsrecY_ecal->Write();
  h_recXvsrecY_hcal->Write();
  h_recXvsrecY_oneevt->Write();

  cout<<"No. of entries survived: "<<nSurvived<<endl;
  cout<<"No. of EH hadrons : "<<EH<<endl;
  cout<<"No. of H hadrons : "<<H<<endl;
  
}


