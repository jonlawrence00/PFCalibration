//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 16 11:59:35 2023 by ROOT version 6.16/00
// from TTree s/ PFCalibration
// found on file: step3.root
//////////////////////////////////////////////////////////

#ifndef skim_h
#define skim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class skim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types

  Float_t  true_;
  Float_t         p;
  Float_t         ecal;
   Float_t         hcal;
   Float_t         ho;
   Float_t         eta;
   Float_t         phi;
   Int_t           charge;
   vector<float>   *dr;
   vector<float>   *Eecal;
   vector<float>   *Ehcal;
   vector<float>   *pfcID;
   vector<int>     *pfcs;
   vector<float>   *correcal;
   vector<float>   *corrhcal;
   Float_t         Ccorrecal;
   Float_t         Ccorrhcal;
   vector<float>   *EcalRechit_posx;
   vector<float>   *EcalRechit_posy;
   vector<float>   *EcalRechit_posz;
   vector<float>   *EcalRechit_E;
   vector<float>   *EcalRechit_dr;
   vector<float>   *EcalRechit_eta;
   vector<float>   *EcalRechit_phi;
   Float_t         EcalRechit_totE;
   vector<float>   *EcalRechit_depth;
   vector<float>   *EcalPFclustereta;
   vector<float>   *HcalRechit_posx;
   vector<float>   *HcalRechit_posy;
   vector<float>   *HcalRechit_posz;
   vector<float>   *HcalRechit_E;
   vector<float>   *HcalRechit_dr;
   vector<float>   *HcalRechit_eta;
   vector<float>   *HcalRechit_phi;
   Float_t         HcalRechit_totE;
   vector<float>   *HcalRechit_depth;
   vector<float>   *HcalPFclustereta;
   ULong64_t       run;
   ULong64_t       evt;
   ULong64_t       lumiBlock;
   ULong64_t       time;

   // List of branches
   TBranch        *b_true_;   //!
   TBranch        *b_p;   //!
   TBranch        *b_ecal;   //!
   TBranch        *b_hcal;   //!
   TBranch        *b_ho;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_dr;   //!
   TBranch        *b_Eecal;   //!
   TBranch        *b_Ehcal;   //!
   TBranch        *b_pfcID;   //!
   TBranch        *b_pfcs;   //!
   TBranch        *b_correcal;   //!
   TBranch        *b_corrhcal;   //!
   TBranch        *b_Ccorrecal;   //!
   TBranch        *b_Ccorrhcal;   //!
   TBranch        *b_EcalRechit_posx;   //!
   TBranch        *b_EcalRechit_posy;   //!
   TBranch        *b_EcalRechit_posz;   //!
   TBranch        *b_EcalRechit_E;   //!
   TBranch        *b_EcalRechit_dr;   //!
   TBranch        *b_EcalRechit_eta;   //!
   TBranch        *b_EcalRechit_phi;   //!
   TBranch        *b_EcalRechit_totE;   //!
   TBranch        *b_EcalRechit_depth;   //!
   TBranch        *b_EcalPFclustereta;   //!
   TBranch        *b_HcalRechit_posx;   //!
   TBranch        *b_HcalRechit_posy;   //!
   TBranch        *b_HcalRechit_posz;   //!
   TBranch        *b_HcalRechit_E;   //!
   TBranch        *b_HcalRechit_dr;   //!
   TBranch        *b_HcalRechit_eta;   //!
   TBranch        *b_HcalRechit_phi;   //!
   TBranch        *b_HcalRechit_totE;   //!
   TBranch        *b_HcalRechit_depth;   //!
   TBranch        *b_HcalPFclustereta;   //!
   TBranch        *b_orun;   //!

   skim(TTree *tree=0);
   virtual ~skim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
  virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef skim_cxx
skim::skim(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../rootfiles/step3_rechit_1264.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("../rootfiles/step3_rechit_1264.root");
	//	 f = new TFile("../../Run3/rootfile/PGun_step3_RECO_1100_2021_Run3.root");
      }
      f->GetObject("s",tree);

   }
   Init(tree);
}

skim::~skim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skim::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dr = 0;
   Eecal = 0;
   Ehcal = 0;
   pfcID = 0;
   pfcs = 0;
   correcal = 0;
   corrhcal = 0;
   EcalRechit_posx = 0;
   EcalRechit_posy = 0;
   EcalRechit_posz = 0;
   EcalRechit_E = 0;
   EcalRechit_dr = 0;
   EcalRechit_eta = 0;
   EcalRechit_phi = 0;
   EcalRechit_depth = 0;
   EcalPFclustereta = 0;
   HcalRechit_posx = 0;
   HcalRechit_posy = 0;
   HcalRechit_posz = 0;
   HcalRechit_E = 0;
   HcalRechit_dr = 0;
   HcalRechit_eta = 0;
   HcalRechit_phi = 0;
   HcalRechit_depth = 0;
   HcalPFclustereta = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("true", &true_, &b_true_);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("ecal", &ecal, &b_ecal);
   fChain->SetBranchAddress("hcal", &hcal, &b_hcal);
   fChain->SetBranchAddress("ho", &ho, &b_ho);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("dr", &dr, &b_dr);
   fChain->SetBranchAddress("Eecal", &Eecal, &b_Eecal);
   fChain->SetBranchAddress("Ehcal", &Ehcal, &b_Ehcal);
   fChain->SetBranchAddress("pfcID", &pfcID, &b_pfcID);
   fChain->SetBranchAddress("pfcs", &pfcs, &b_pfcs);
   fChain->SetBranchAddress("correcal", &correcal, &b_correcal);
   fChain->SetBranchAddress("corrhcal", &corrhcal, &b_corrhcal);
   fChain->SetBranchAddress("Ccorrecal", &Ccorrecal, &b_Ccorrecal);
   fChain->SetBranchAddress("Ccorrhcal", &Ccorrhcal, &b_Ccorrhcal);
   fChain->SetBranchAddress("EcalRechit_posx", &EcalRechit_posx, &b_EcalRechit_posx);
   fChain->SetBranchAddress("EcalRechit_posy", &EcalRechit_posy, &b_EcalRechit_posy);
   fChain->SetBranchAddress("EcalRechit_posz", &EcalRechit_posz, &b_EcalRechit_posz);
   fChain->SetBranchAddress("EcalRechit_E", &EcalRechit_E, &b_EcalRechit_E);
   fChain->SetBranchAddress("EcalRechit_dr", &EcalRechit_dr, &b_EcalRechit_dr);
   fChain->SetBranchAddress("EcalRechit_eta", &EcalRechit_eta, &b_EcalRechit_eta);
   fChain->SetBranchAddress("EcalRechit_phi", &EcalRechit_phi, &b_EcalRechit_phi);
   fChain->SetBranchAddress("EcalRechit_totE", &EcalRechit_totE, &b_EcalRechit_totE);
   fChain->SetBranchAddress("EcalRechit_depth", &EcalRechit_depth, &b_EcalRechit_depth);
   fChain->SetBranchAddress("EcalPFclustereta", &EcalPFclustereta, &b_EcalPFclustereta);
   fChain->SetBranchAddress("HcalRechit_posx", &HcalRechit_posx, &b_HcalRechit_posx);
   fChain->SetBranchAddress("HcalRechit_posy", &HcalRechit_posy, &b_HcalRechit_posy);
   fChain->SetBranchAddress("HcalRechit_posz", &HcalRechit_posz, &b_HcalRechit_posz);
   fChain->SetBranchAddress("HcalRechit_E", &HcalRechit_E, &b_HcalRechit_E);
   fChain->SetBranchAddress("HcalRechit_dr", &HcalRechit_dr, &b_HcalRechit_dr);
   fChain->SetBranchAddress("HcalRechit_eta", &HcalRechit_eta, &b_HcalRechit_eta);
   fChain->SetBranchAddress("HcalRechit_phi", &HcalRechit_phi, &b_HcalRechit_phi);
   fChain->SetBranchAddress("HcalRechit_totE", &HcalRechit_totE, &b_HcalRechit_totE);
   fChain->SetBranchAddress("HcalRechit_depth", &HcalRechit_depth, &b_HcalRechit_depth);
   fChain->SetBranchAddress("HcalPFclustereta", &HcalPFclustereta, &b_HcalPFclustereta);
   fChain->SetBranchAddress("run", &run, &b_orun);
   fChain->SetBranchAddress("evt", &evt, &b_orun);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_orun);
   fChain->SetBranchAddress("time", &time, &b_orun);
   Notify();
}

Bool_t skim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skim_cxx
