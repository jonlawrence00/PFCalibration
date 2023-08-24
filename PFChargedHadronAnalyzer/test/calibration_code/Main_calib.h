#ifndef MAIN_UL_H
#define MAIN_UL_H

#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include "TH2F.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <string>
#include <iostream>
#include <math.h>

using namespace std;

double sigC_ = 5.;
// unsigned sampleRangeHigh = 200;
unsigned sampleRangeHigh = 500;


//threshold
/////////////MM
double aEH = 3.5;
double aE = 3.5;

// double aEH = 4.5;
// double aE = 4.5;

double aH = 2.5;//3.0;

double aEHe = 3.5;
double aEe = 3.5;

//spandey Feb 23 2018
// double aEHe = 4.5;
// double aEe = 4.5;

double aHe = 2.5;
/////////////MM


/*
//threshold
/////////////SP
double aEH = 3.675;
double aE = 3.675;
double aH = 2.625;//3.0;

double aEHe = 3.675;
double aEe = 3.675;
double aHe = 2.625;
/////////////SP
*/


double lBs = 0.25; //0.25 B //0.5 E
double mBs=1.; //1 B
double hBs=5.0; //10 40 B 50 E

double RBE = 4; //rebinning for alpha-beta curves

//Eta factor : ecal*factor + hcal
double factorB = 1.3;//Optimal //A factor put in by hand to make the eta dependence agree
double factorE = 1.3; //better. Should fit for the value later. 
//double factorE = 1.7; //better. Should fit for the value later.  //shubham

vector<double> sigmas;
TF1* faBarrel;
TF1* fbBarrel;
TF1* fcBarrel;
TF1* faEtaBarrel;
TF1* fbEtaBarrel;


TF1* faBarrel52x;
TF1* fbBarrel52x;
TF1* fcBarrel52x;
TF1* faEtaBarrel52x;
TF1* fbEtaBarrel52x;

void LoadOldThresholds() {
  aEH = 3.5; aE = 3.5; aH = 2.5; //3.5 2.5
  aEHe = 3.5; aEe = 3.5; aHe = 2.5;
}

//spandey
void LoadNewThresholds() {
  aEH = 3.8; aE = 3.8; aH = 3.0; //3.5 2.5
  aEHe = 3.8; aEe = 3.8; aHe = 3.0;
}



// Double_t median1(TH1 *h1) {
//    //compute the median for 1-d histogram h1
//    Int_t nbins = h1->GetXaxis()->GetNbins();
//    Double_t *x = new Double_t[nbins];
//    Double_t *y = new Double_t[nbins];
//    for (Int_t i=0;i<nbins;i++) {
//       x[i] = h1->GetXaxis()->GetBinCenter(i+1);
//       y[i] = h1->GetBinContent(i+1);
//    }
//    Double_t median = TMath::Median(nbins,x,y);
//    delete [] x;
//    delete [] y;
//    return median;
// } 

void InitBarrelAlpha() {

  // faBarrel = new TF1("faBarrel","[0]+((([1]+([2]/(x^[8])))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]*(x/[5])))))",1.,1000.);
  // fbBarrel = new TF1("fbBarrel","[0]+([1]+[2]/x^1.1)*exp(-x/[3])-[4]*exp(-x^[5]*x/2000.)",1.,1000.);
  // fcBarrel = new TF1("fcBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  // faEtaBarrel = new TF1("faEtaBarrel","[0]+[1]*exp(-x^[3]/[2])+([4]/x^[5])*exp([6]*x)",1.,1000.);
  // fbEtaBarrel = new TF1("fbEtaBarrel","[0]+([1]*x+[2])*exp(-(x+[4])/[3])",1.,1000.);
  // faBarrel->SetParameter(0,1.04193);
  // fbBarrel->SetParameter(0,0.992287);
  // fcBarrel->SetParameter(0,0.978144);
  // faEtaBarrel->SetParameter(0,-0.00327152);
  // fbEtaBarrel->SetParameter(0,-6.58476);
  // faBarrel->SetParameter(1,16.2881);
  // fbBarrel->SetParameter(1,0.216829);
  // fcBarrel->SetParameter(1,0.0451738);
  // faEtaBarrel->SetParameter(1,-0.00327152);
  // fbEtaBarrel->SetParameter(1,0.0273361);
  // faBarrel->SetParameter(2,-14.9307);
  // fbBarrel->SetParameter(2,-5.0668);
  // fcBarrel->SetParameter(2,-0.49014);
  // fbEtaBarrel->SetParameter(2,12.4812);
  // faBarrel->SetParameter(3,0.80573);
  // fbBarrel->SetParameter(3,42.8035);
  // fcBarrel->SetParameter(3,8.84333);
  // fbEtaBarrel->SetParameter(3,505.069);
  // faBarrel->SetParameter(4,3.06364);
  // fbBarrel->SetParameter(4,0.0832237);
  // fcBarrel->SetParameter(4,1.18357);
  // faEtaBarrel->SetParameter(4,-0.00327152);
  // fbEtaBarrel->SetParameter(4,318.064);
  // faBarrel->SetParameter(5,9.23088);
  // fbBarrel->SetParameter(5,0.921352);
  // fcBarrel->SetParameter(5,12.6626);
  // faBarrel->SetParameter(6,0.174875);
  // faBarrel->SetParameter(7,0.522519);
  // faBarrel->SetParameter(8,0.0568171);

  // Barrel (fit made with |eta| < 1.2)
    faBarrel = new TF1("faBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
    fbBarrel = new TF1("fbBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
    fcBarrel = new TF1("fcBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
    faEtaBarrel = new TF1("faEtaBarrel","[0]+[1]*exp(-x/[2])",1.,1000.);
    fbEtaBarrel = new TF1("fbEtaBarrel","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
    faBarrel->SetParameter(0,1.15665);
    fbBarrel->SetParameter(0,0.994603);
    fcBarrel->SetParameter(0,0.956544);
    faEtaBarrel->SetParameter(0,0.014664);
    fbEtaBarrel->SetParameter(0,0.00975451);
    faBarrel->SetParameter(1,0.165627);
    fbBarrel->SetParameter(1,0.13632);
    fcBarrel->SetParameter(1,0.0857207);
    faEtaBarrel->SetParameter(1,-0.0426776);
    fbEtaBarrel->SetParameter(1,0.102247);
    faBarrel->SetParameter(2,0.827718);
    fbBarrel->SetParameter(2,-0.758013);
    fcBarrel->SetParameter(2,-0.44347);
    faEtaBarrel->SetParameter(2,431.054);
    fbEtaBarrel->SetParameter(2,436.21);
    faBarrel->SetParameter(3,231.339);
    fbBarrel->SetParameter(3,183.627);
    fcBarrel->SetParameter(3,63.3479);
    faBarrel->SetParameter(4,2.45332);
    fbBarrel->SetParameter(4,1);
    fcBarrel->SetParameter(4,1.24174);
    faBarrel->SetParameter(5,29.6603);
    fbBarrel->SetParameter(5,39.6784);
    fcBarrel->SetParameter(5,12.322);





  faBarrel->SetLineColor(kBlue+1);
  faBarrel->SetLineStyle(2);
  faBarrel->SetLineWidth(2);
  fbBarrel->SetLineColor(kBlue+1);
  fbBarrel->SetLineStyle(2);
  fbBarrel->SetLineWidth(2);
  fcBarrel->SetLineColor(kBlue+1);
  fcBarrel->SetLineStyle(2);
  fcBarrel->SetLineWidth(2);
  faEtaBarrel->SetLineColor(kBlue+1);
  faEtaBarrel->SetLineStyle(2);
  faEtaBarrel->SetLineWidth(2);
  fbEtaBarrel->SetLineColor(kBlue+1);
  fbEtaBarrel->SetLineStyle(2);
  fbEtaBarrel->SetLineWidth(2);

  faBarrel->SetTitle("Ecal coef 4XY");
  fbBarrel->SetTitle("Hcal coef 4XY");
  faEtaBarrel->SetTitle("#alpha coef 4XY");
  fbEtaBarrel->SetTitle("#beta coef 4XY");



 // faBarrel52x = new TF1("faBarrel52X","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
 //  fbBarrel52x = new TF1("fbBarrel52X","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
 //  fcBarrel52x = new TF1("fcBarrel52X","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
 //  faEtaBarrel52x = new TF1("faEtaBarrel52X","[0]+[1]*exp(-x/[2])",1.,1000.);
 //  fbEtaBarrel52x = new TF1("fbEtaBarrel52X","[0]+[1]*exp(-x/[2])",1.,1000.);
 //  faBarrel52x->SetParameter(0,1.30727);
 //  fbBarrel52x->SetParameter(0,0.998609);
 //  fcBarrel52x->SetParameter(0,0.961909);
 //  faEtaBarrel52x->SetParameter(0,256.375);
 //  fbEtaBarrel52x->SetParameter(0,0.0236908);
 //  faBarrel52x->SetParameter(1,0.297017);
 //  fbBarrel52x->SetParameter(1,0.185606);
 //  fcBarrel52x->SetParameter(1,2.48019);
 //  faEtaBarrel52x->SetParameter(1,-256.412);
 //  fbEtaBarrel52x->SetParameter(1,0.144532);
 //  faBarrel52x->SetParameter(2,0.176997);
 //  fbBarrel52x->SetParameter(2,-2.4317);
 //  fcBarrel52x->SetParameter(2,-6.18656);
 //  faEtaBarrel52x->SetParameter(2,1.02648e+06);
 //  fbEtaBarrel52x->SetParameter(2,122.074);
 //  faBarrel52x->SetParameter(3,71.526);
 //  fbBarrel52x->SetParameter(3,11.5238);
 //  fcBarrel52x->SetParameter(3,7.32584);
 //  faBarrel52x->SetParameter(4,2.16376);
 //  fbBarrel52x->SetParameter(4,1);
 //  fcBarrel52x->SetParameter(4,0.263953);
 //  faBarrel52x->SetParameter(5,31.9415);
 //  fbBarrel52x->SetParameter(5,9.98657);
 //  fcBarrel52x->SetParameter(5,387.311);
 
 faBarrel52x = new TF1("faBarrel52x","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fbBarrel52x = new TF1("fbBarrel52x","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fcBarrel52x = new TF1("fcBarrel52x","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  faEtaBarrel52x = new TF1("faEtaBarrel52x","[0]+[1]*exp(-x/[2])",1.,1000.);
  fbEtaBarrel52x = new TF1("fbEtaBarrel52x","[0]+[1]*exp(-x/[2])",1.,1000.);
  faBarrel52x->SetParameter(0,1.36776);
  fbBarrel52x->SetParameter(0,1.00174);
  fcBarrel52x->SetParameter(0,0.958822);
  faEtaBarrel52x->SetParameter(0,-0.0235949);
  fbEtaBarrel52x->SetParameter(0,-20.0611);
  faBarrel52x->SetParameter(1,0.401817);
  fbBarrel52x->SetParameter(1,-0.575566);
  fcBarrel52x->SetParameter(1,-20.7743);
  faEtaBarrel52x->SetParameter(1,0.145271);
  fbEtaBarrel52x->SetParameter(1,20.1811);
  faBarrel52x->SetParameter(2,-0.272423);
  fbBarrel52x->SetParameter(2,0.511645);
  fcBarrel52x->SetParameter(2,22.8622);
  faEtaBarrel52x->SetParameter(2,3.76578);
  fbEtaBarrel52x->SetParameter(2,733806);
  faBarrel52x->SetParameter(3,39.0306);
  fbBarrel52x->SetParameter(3,11.098);
  fcBarrel52x->SetParameter(3,1.29814);
  faBarrel52x->SetParameter(4,2.30666);
  fbBarrel52x->SetParameter(4,1);
  fcBarrel52x->SetParameter(4,0.223834);
  faBarrel52x->SetParameter(5,27.6117);
  fbBarrel52x->SetParameter(5,34.2481);
  fcBarrel52x->SetParameter(5,84.6003);


  faBarrel52x->SetLineColor(kOrange-3);
  faBarrel52x->SetLineStyle(7);
  faBarrel52x->SetLineWidth(2);
  fbBarrel52x->SetLineColor(kOrange-3);
  fbBarrel52x->SetLineStyle(7);
  fbBarrel52x->SetLineWidth(2);
  fcBarrel52x->SetLineColor(kOrange-3);
  fcBarrel52x->SetLineStyle(7);
  fcBarrel52x->SetLineWidth(2);
  faEtaBarrel52x->SetLineColor(kOrange-3);
  faEtaBarrel52x->SetLineStyle(7);
  faEtaBarrel52x->SetLineWidth(2);
  fbEtaBarrel52x->SetLineColor(kOrange-3);
  fbEtaBarrel52x->SetLineStyle(7);
  fbEtaBarrel52x->SetLineWidth(2);

  faBarrel52x->SetTitle("Ecal coef 52X");
  fbBarrel52x->SetTitle("Hcal coef 52X");
  faEtaBarrel52x->SetTitle("#alpha coef 52X");
  fbEtaBarrel52x->SetTitle("#beta coef 52X");




}



///////////////////////////////////////////////////////////////////////////////
//Container class that holds all the coefficients for a particular Etrue bin
///////////////////////////////////////////////////////////////////////////////

class ABC
{
   private:

      vector<double> ETrueEnergies_; //Variables that fall in this ETrue bin, 
      vector<double> ecalEnergies_;  //which is defined by binLowEdge and 
      vector<double> hcalEnergies_;  //binHighEdge. sigmaEcalHcal is simply the
      vector<double> etas_;          //uncertainty from the ecal and hcal(not
      vector<double> sigmaEcalHcal_; //uncertainty in the calibration consts)
      double binLowEdge_;
      double binHighEdge_;
      double etaMinFit_;     //From what etas to fit the B and C constants.
      double etaMaxFit_;     
      double etaMinEtaFit_;  //From what etas to fit the Alpha and Beta 
      double etaMaxEtaFit_;  //constants.
      bool isBarrel_;
      

      double ETrueAverage_; //Average and RMS of the ETrue values in the bin.
      double ETrueRMS_;
      double a_;  //Constant values.
      double b_;
      double c_;

      double sigmaB_;  //Uncertainty in the constants
      double sigmaC_;

   public:
      ABC(double binLowEdge, double binHighEdge, bool isBarrel); 

      bool addEntry(double ETrueEnergy, double ecalEnergy, double hcalEnergy,
                    double eta);  //Adds an event to the ETrue bin
      double getBinLowEdge();
      double getBinHighEdge();
      bool isBarrel();  //Checks if it is a barrel-type constant storage
      bool isEmpty();   //Checks if its empty
      bool isEmptyInFitRange();   //Checks if its empty in eta fit range
      unsigned getSize();        //Returns the various stored variables in the
      double getETrueAverage();  //ABC object
      double getETrueRMS();
      double getA();
      double getB();
      double getC();
      double getSigmaB();
      double getSigmaC();
       
      double getETrue(unsigned i);
      double getEcal(unsigned i); //Returns b*ecal for entry i
      double getHcal(unsigned i); //Returns c*hcal for entry i
      double getEta(unsigned i);
      double getNEntries();

      void computeETrueAverage();  //Computes the various calibration constants
      void computeETrueRMS();      //and other stored elements in the object.
      void computeA(double a);     //Right now the constant "a" is not computed
      void computeB();             //but rather just set.
      void computeC();
      bool computeBC();
      void clear();
};

ABC::ABC(double binLowEdge, double binHighEdge, 
	 bool _isBarrel) 

{
  binLowEdge_ = binLowEdge;
  binHighEdge_ = binHighEdge;
  isBarrel_ = _isBarrel;

   if(isBarrel_)
   {
      etaMinFit_ = 0.0;
      etaMaxFit_ = 1.6;//1.0
      etaMinEtaFit_ = 0.0;
      etaMaxEtaFit_ = 1.6;//1.3
   }
   else
   {
      etaMinFit_ = 1.6;
      etaMaxFit_ = 2.2; //FIXME 2.2
      etaMinEtaFit_ = 1.6;
      etaMaxEtaFit_ = 2.8; //FIXME 2.8
   }
   
   a_ = 0;
   b_ = 0;
   c_ = 0;
   ETrueAverage_ = 0;
   ETrueRMS_ = 0;
   sigmaB_ = 0;
   sigmaC_ = 0;
}
bool ABC::addEntry(double ETrue, double ecalEnergy, double hcalEnergy, double eta)
{

   
   double sigmaEcalHcal;

   if(isBarrel_) 
      sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max(ecalEnergy + 
                                                           hcalEnergy, 1.0)));
   else
      sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max(ecalEnergy + 
                                                           hcalEnergy, 1.0)));

   if( fabs(ecalEnergy + hcalEnergy + a_ - ETrue) > sigC_*sigmaEcalHcal || 
      (ecalEnergy + hcalEnergy) < 0.5 || ETrue < 1.0 || ETrue < binLowEdge_ || 
      ETrue> binHighEdge_ || fabs(eta) > etaMaxFit_ || fabs(eta) < etaMinFit_ )
      return false;


   ETrueEnergies_.push_back(ETrue);
   ecalEnergies_.push_back(ecalEnergy);
   hcalEnergies_.push_back(hcalEnergy);
   etas_.push_back(eta);
   sigmaEcalHcal_.push_back(sigmaEcalHcal);

   // if(ecalEnergy==0 && isBarrel_)
   //   cout<<hcalEnergy<<"  "<<fabs(ecalEnergy + hcalEnergy + a_ - ETrue)<<"   "<<sigmaEcalHcal<<endl;

   return true;
}

double ABC::getBinLowEdge() {return binLowEdge_;}
double ABC::getBinHighEdge() {return binHighEdge_;}
double ABC::getETrueAverage() {return ETrueAverage_;}
double ABC::getETrueRMS() {return ETrueRMS_;}
double ABC::getA() {return a_;}
double ABC::getB() {return b_;}
double ABC::getC() {return c_;}
double ABC::getSigmaB() {return sigmaB_;}
double ABC::getSigmaC() {return sigmaC_;}

double ABC::getETrue(unsigned i) {return ETrueEnergies_[i];}
double ABC::getEcal(unsigned i) {return ecalEnergies_[i];}
double ABC::getHcal(unsigned i) {return hcalEnergies_[i];}
double ABC::getEta(unsigned i) {return etas_[i];}
double ABC::getNEntries(){return ETrueEnergies_.size();}
bool ABC::isBarrel() {return isBarrel_;}
bool ABC::isEmpty() 
{
   if(ETrueEnergies_.size() == 0) return true;
   else return false;
}
bool ABC::isEmptyInFitRange() 
{
  //  cout<<"etas_.size():"<<etas_.size()<<endl;
   for(unsigned i = 0; i < etas_.size(); i++)
   {
      if(fabs(etas_[i]) < etaMaxFit_ && fabs(etas_[i]) > etaMinFit_)
         return false;
   }
   return true;
}
unsigned ABC::getSize()
{
   return ETrueEnergies_.size();
}
void ABC::computeETrueAverage()
{
   double totalETrue = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
     if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }   

     totalETrue += ETrueEnergies_[i];
   }
   
   ETrueAverage_ = totalETrue/(ETrueEnergies_.size() - numberSkipped);
} 
void ABC::computeETrueRMS()
{
   double totalETrueSquared = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i<ETrueEnergies_.size(); i++)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }
       
      totalETrueSquared += ETrueEnergies_[i]*ETrueEnergies_[i];
   }
   
   ETrueRMS_ = sqrt(
      (totalETrueSquared/(ETrueEnergies_.size() - numberSkipped) -
       ETrueAverage_*ETrueAverage_)/(ETrueEnergies_.size() - numberSkipped)); 
} 
void ABC::computeA(double a) 
{
   a_ = a; 
}
void ABC::computeB() 
{
   double totalEcalSquared = 0;
   double totalEMinusATimesEcal = 0;
   
   for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;

      totalEcalSquared += 2* ecalEnergies_[i]*ecalEnergies_[i]/
      (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      totalEMinusATimesEcal += 2*(ETrueEnergies_[i] - a_)*ecalEnergies_[i]/
      (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
   }
   
   b_ = totalEMinusATimesEcal/totalEcalSquared;
   sigmaB_ = sqrt(1/totalEcalSquared);
}
void ABC::computeC() 
{
   double totalHcalSquared = 0;
   double totalEMinusATimesHcal = 0;
   
   //cout<<" size C : "<<ETrueEnergies_.size()<<endl;

   for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
     totalHcalSquared += 2*hcalEnergies_[i]*hcalEnergies_[i]/
       (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
     totalEMinusATimesHcal += 2*(ETrueEnergies_[i] - a_)*hcalEnergies_[i]/
       (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
   }

   c_ = totalEMinusATimesHcal/totalHcalSquared;
   //cout<<totalEMinusATimesHcal<<"   "<<totalHcalSquared<<"   "<<c_<<endl;
   //cout<<"totalEMinusATimesHcal: "<<totalEMinusATimesHcal<<",   totalHcalSquared: "<<totalHcalSquared<<",   c_: "<<c_<<endl;

   // //FIXME MM
   // TH1F* cplot=new TH1F("cplot","",3000,-1, 5); 
   // float median=c_;
   // for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   //   {
   //     cplot->Fill( (ETrueEnergies_[i] - a_)/hcalEnergies_[i] );
   //   }

   // //  cout<<ETrueEnergies_[0]<<"  "<<" c_ "<<c_<<"   "<<(c_ + median1(cplot))/2.<<"   "<<median1(cplot)<<endl;
   // // if(ETrueEnergies_[0]>20)
   // //   c_ = (c_ + median1(cplot))/2.;
   // if(ETrueEnergies_[0]<=20)
   //   c_ = median1(cplot);

   // delete cplot;
   
   //================================


   sigmaC_ = sqrt(1/totalHcalSquared);
}
bool ABC::computeBC()
{
   ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs;
   ROOT::Math::SVector<double, 2> consts;
   ROOT::Math::SVector<double, 2> values;
   bool isInverted;
 
   coeffs(0, 0) = 0;
   coeffs(0, 1) = 0;
   coeffs(1, 0) = 0;
   coeffs(1, 1) = 0;
   consts(0) = 0;
   consts(1) = 0;
   //Create the matrix that will be inverted and the vector which will multiply
   //that matrix to find the b and c calibration constants.

   //cout<<" size AB : "<<ETrueEnergies_.size()<<endl;

   for(unsigned i = 0; i < ETrueEnergies_.size(); ++i)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;
      
      // if(i<10)
      // 	cout<<ETrueEnergies_[i]<<"   "<<ecalEnergies_[i]<<"   "<<hcalEnergies_[i]<<"   "<<etas_[i]<<"  "<<sigmaEcalHcal_[i]<<"   "<<a_;


      coeffs(0, 0) += 2*ecalEnergies_[i]*ecalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(0, 1) += 2*ecalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(1, 0) += 2*ecalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(1, 1) += 2*hcalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      consts(0) += 2*(ETrueEnergies_[i] - a_)*ecalEnergies_[i] /
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      consts(1) += 2*(ETrueEnergies_[i] - a_)*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

      // if(i<10)
      // 	cout<<" ---> "<<coeffs(0, 0)<<"  "<<coeffs(0, 0)
      // 	    <<"  "<<coeffs(0, 0)<<"  "<<coeffs(0, 0)
      // 	    <<"  "<<consts(0)<<"  "<<consts(1)<<endl;

   }
   isInverted = coeffs.Invert();
   // coeffs.Print(cout); cout<<endl;
   // consts.Print(cout);


   if(isInverted && sqrt(coeffs(0,0)) <  100000) //Make sure it inverted successfully (i.e. det != 0)
   {

      values = coeffs*consts;
      
      b_ = values(0);
      c_ = values(1);
      sigmaB_ = sqrt(coeffs(0,0));
      sigmaC_ = sqrt(coeffs(1,1));

      //cout<<" b: "<<b_<<"  c: "<<c_<<endl;

      // //FIXME MM
      // TH1F* cplot=new TH1F("cplot","",3000,-1, 5); 
      //   for(unsigned i = 0; i < ETrueEnergies_.size(); ++i)
      // 	  {
      // 	    if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;
	    
      // 	    cplot->Fill( (ETrueEnergies_[i] - a_ - b_*ecalEnergies_[i])/hcalEnergies_[i] ); 
      // 	  }

      // 	//	cout<<ETrueEnergies_[0]<<"  "<<" c_ "<<c_<<"   "<<(c_ + median1(cplot))/2.<<"   "<<median1(cplot)<<endl;
	
      // 	c_ = (c_ + median1(cplot))/2.;
      // 	//	c_ = median1(cplot);

      // 	delete cplot;
	//================================
	
      return true;


   }
   else return false;
}

void ABC::clear()
{
   ETrueEnergies_.clear();
   ecalEnergies_.clear();
   hcalEnergies_.clear();
   etas_.clear();
   
   delete &binLowEdge_;
   delete &binHighEdge_;
   delete &isBarrel_;
   delete &ETrueAverage_;
   delete &ETrueRMS_;
   delete &a_;
   delete &b_;
   delete &c_;
   delete &sigmaB_;
   delete &sigmaC_;
}





class AlphaBeta
{
private:
  vector<double> ETrueEnergies_;
  vector<double> ecalEnergies_;
  vector<double> hcalEnergies_;
  vector<double> etas_;
  vector<double> sigmaEcalHcal_;
  vector<double> a_;
      

  double binLowEdge_;
  double binHighEdge_;
  double etaMinFit_;
  double etaMaxFit_;
  bool isBarrel_;

  double ETrueAverage_;
  double ETrueRMS_;
  double alpha_;
  double gamma_;
  double beta_;
  double sigmaAlpha_;
  double sigmaGamma_;
  double sigmaBeta_;

public:
  AlphaBeta(double binLowEdge, double binHighEdge, bool isBarrel);
  bool addEntry(double ETrueEnergy, double ecalEnergies, double correctedHcal, double eta);  //Adds an event to the ETrue bin
  double getBinLowEdge();
      double getBinHighEdge();
      bool isBarrel();  //Checks if it is a barrel-type constant storage
      bool isEmpty();   //Checks if its empty
      bool isEmptyInFitRange();   //Checks if its empty in eta fit range
      unsigned getSize();        //Returns the various stored variables in the
      double getETrue(unsigned i);
      double getEcal(unsigned i);
      double getHcal(unsigned i);

      double getETrueAverage();  //AlphaBeta object
      double getETrueRMS();
      double getAlpha();
      double getBeta();
  double getGamma();
      double getSigmaAlpha();
      double getSigmaBeta();
  double getSigmaGamma();

      void correctEcal(unsigned i, double b);
      void correctHcal(unsigned i, double c);
      void computeSigmaEcalHcal();
      void computeETrueAverage();  //Computes the various calibration constants
      void computeETrueRMS();      //and other stored elements in the object.
      bool computeAlphaBeta();
  bool computeAlphaBetaGamma();
      void clear();
      
};

AlphaBeta::AlphaBeta(double binLowEdge, double binHighEdge, bool _isBarrel) 

{
   binLowEdge_ = binLowEdge;
   binHighEdge_ = binHighEdge;
   isBarrel_ = _isBarrel;

   if(isBarrel_)
   {
      etaMinFit_ = 0.0;
      etaMaxFit_ =1.2;// 1.5;//1.2
   }
   else
   {
     etaMinFit_ = 1.6; //FIXME 1.6
     etaMaxFit_ = 2.8; //FIXME 2.8
   }
   
   alpha_ = 0;
   beta_ = 0;
   ETrueAverage_ = 0;
   ETrueRMS_ = 0;
   sigmaAlpha_ = 0;
   sigmaBeta_ = 0;
}
bool AlphaBeta::addEntry(double ETrue, double ecal, 
                         double hcal, double eta)
{
   double a = 0;

   if(isBarrel_) 
   {
    
      if(ecal > 0)
         a = aEH;
      else
         a = aH;
   }
   else
   {
      if(ecal > 0)
         a = aEHe;
      else
         a = aHe;
   }
   if((ecal + hcal) < 0.5 || ETrue < 1.0 || 
      ETrue < binLowEdge_ || ETrue > binHighEdge_ || fabs(eta) > etaMaxFit_ || 
      fabs(eta) < etaMinFit_) return false;

   a_.push_back(a);
   ETrueEnergies_.push_back(ETrue);
   ecalEnergies_.push_back(ecal);
   hcalEnergies_.push_back(hcal);
   etas_.push_back(eta);

   return true;
}

double AlphaBeta::getBinLowEdge() {return binLowEdge_;}
double AlphaBeta::getBinHighEdge() {return binHighEdge_;}
double AlphaBeta::getETrueAverage() {return ETrueAverage_;}
double AlphaBeta::getETrueRMS() {return ETrueRMS_;}
double AlphaBeta::getAlpha() {return alpha_;}
double AlphaBeta::getBeta() {return beta_;}
double AlphaBeta::getGamma() {return gamma_;}
double AlphaBeta::getSigmaAlpha() {return sigmaAlpha_;}
double AlphaBeta::getSigmaBeta() {return sigmaBeta_;}
double AlphaBeta::getSigmaGamma() {return sigmaGamma_;}
bool AlphaBeta::isBarrel() {return isBarrel_;}
bool AlphaBeta::isEmpty() 
{
   if(ETrueEnergies_.size() == 0) return true;
   else return false;
}
bool AlphaBeta::isEmptyInFitRange() 
{
   for(unsigned i = 0; i < etas_.size(); i++)
   {
      if(fabs(etas_[i]) < etaMaxFit_ && fabs(etas_[i]) > etaMinFit_)
      {
         return false;
      }
   }
   return true;
}
unsigned AlphaBeta::getSize()
{
   return ETrueEnergies_.size();
}
double AlphaBeta::getETrue(unsigned i){return ETrueEnergies_[i];}
double AlphaBeta::getEcal(unsigned i){return ecalEnergies_[i];}
double AlphaBeta::getHcal(unsigned i){return hcalEnergies_[i];}
void AlphaBeta::correctEcal(unsigned i, double b)
{
   ecalEnergies_[i] = b* ecalEnergies_[i];
}
void AlphaBeta::correctHcal(unsigned i, double c)
{
   hcalEnergies_[i] = c* hcalEnergies_[i];
}
void AlphaBeta::computeSigmaEcalHcal()
{
   double sigmaEcalHcal;
   double ecal;
   double hcal;
   //double a;
   vector<bool> erase;
   

   for(unsigned entry = 0; entry < ETrueEnergies_.size(); entry++)
   {
      ecal = ecalEnergies_[entry];
      hcal = hcalEnergies_[entry];
      //a = a_[entry];
      erase.push_back(false);
      
      if(isBarrel_)
         sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*
                              (std::max(ecal + hcal, 1.0)));
      else
         sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*
                              (std::max(ecal + hcal, 1.0)));
      
      sigmaEcalHcal_.push_back(sigmaEcalHcal);
         
   }


}
void AlphaBeta::computeETrueAverage()
{
   double totalETrue = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
     if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }   

     totalETrue += ETrueEnergies_[i];
   }
   
   ETrueAverage_ = totalETrue/(ETrueEnergies_.size() - numberSkipped);
} 
void AlphaBeta::computeETrueRMS()
{
   double totalETrueSquared = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i<ETrueEnergies_.size(); i++)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }
       
      totalETrueSquared += ETrueEnergies_[i]*ETrueEnergies_[i];
   }
   
   ETrueRMS_ = sqrt(
      (totalETrueSquared/(ETrueEnergies_.size() - numberSkipped) -
       ETrueAverage_*ETrueAverage_)/(ETrueEnergies_.size() - numberSkipped)); 
} 



bool AlphaBeta::computeAlphaBeta()
{
   ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs;
   ROOT::Math::SVector<double, 2> consts;
   ROOT::Math::SVector<double, 2> values;
   bool isInverted;
   vector<double> etaPow;
  

   coeffs(0, 0) = 0;
   coeffs(0, 1) = 0;
   coeffs(1, 0) = 0;
   coeffs(1, 1) = 0;
   consts(0) = 0;
   consts(1) = 0;
   int N=0;

   if(isBarrel_) 
   {
      for(unsigned i = 0; i < etas_.size(); i++) 
         etaPow.push_back(etas_[i]*etas_[i]);
      

      for(unsigned i = 0; i < etas_.size(); ++i)
      { 
	if( fabs(ecalEnergies_[i] + hcalEnergies_[i] - ETrueEnergies_[i] + a_[i]) 
             > sigC_*sigmaEcalHcal_[i]) continue;   
         if(hcalEnergies_[i] == 0)continue;
         
	 if(ecalEnergies_[i]>0 ) {
         coeffs(0, 0) += 2.0*
	   ecalEnergies_[i]*ecalEnergies_[i]/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(0, 1) += 2.0*etaPow[i]*
	   ecalEnergies_[i]*ecalEnergies_[i]/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 0) += 2.0*etaPow[i]*
	   ecalEnergies_[i]*ecalEnergies_[i]/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 1) += 2.0*etaPow[i]*etaPow[i]*
	   ecalEnergies_[i]*ecalEnergies_[i]/
	   (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         consts(0) += 2.0*(ETrueEnergies_[i] -
			   a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*ecalEnergies_[i]/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);  
         
         consts(1) += 2.0*etaPow[i]*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*ecalEnergies_[i]/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
	 }
	 else {
	   coeffs(0, 0) += 2.0*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(0, 1) += 2.0*etaPow[i]*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 0) += 2.0*etaPow[i]*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 1) += 2.0*etaPow[i]*etaPow[i]*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         consts(0) += 2.0*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);  
         
         consts(1) += 2.0*etaPow[i]*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorB*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
	 }

    }
   }
   else
   {

     for(unsigned i = 0; i < etas_.size(); i++) {
       //etaPow.push_back((fabs(etas_[i]))*(fabs(etas_[i]) ) *  (fabs(etas_[i]))*(fabs(etas_[i]) ) );
       if (fabs(etas_[i]) > 2.5) {
       etaPow.push_back((fabs(etas_[i]) - 1.5)*(fabs(etas_[i]) - 1.5)  
			*(fabs(etas_[i]) - 1.5)*(fabs(etas_[i]) - 1.5));
       }
       else if (fabs(etas_[i]) < 2.5){
	 //etaPow.push_back((fabs(etas_[i]) - 1.5));
	 etaPow.push_back(0.103-0.053*(fabs(etas_[i]) - 2.0));
       }
	 
     }
     
      for(unsigned i = 0; i < etas_.size(); ++i)
      {  
         if ( fabs(ecalEnergies_[i] + hcalEnergies_[i] - ETrueEnergies_[i] + a_[i]) > sigC_*sigmaEcalHcal_[i]) continue;             
         if(hcalEnergies_[i] == 0)continue;
         
	 N++;

         coeffs(0, 0) += 2.0*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(0, 1) += 2.0*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 0) += 2.0*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 1) += 2.0*etaPow[i]*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         consts(0) += 2.0*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);  
         
         consts(1) += 2.0*etaPow[i]*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      }
      // coeffs.Print(cout);
      // cout<<"N = "<<N<<endl;
   }

   //Create the matrix that will be inverted and the vector which will multiply
   //that matrix to find the alpha and beta calibration constants.
  

   isInverted = coeffs.Invert();
   
   if(isInverted && sqrt(coeffs(0,0)) < 100000)
   {
      values = coeffs*consts;
      
      alpha_ = values(0);
      beta_ = values(1);
      sigmaAlpha_ = sqrt(coeffs(0,0));
      sigmaBeta_ = sqrt(coeffs(1,1));

      return true;
   }
   else
   {
      return false;
   }
} 


bool AlphaBeta::computeAlphaBetaGamma()
{
   ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs;
   ROOT::Math::SVector<double, 3> consts;
   ROOT::Math::SVector<double, 3> values;
   bool isInverted;
   vector<double> etaPow2;
   vector<double> etaPow;

   coeffs(0, 0) = 0;
   coeffs(0, 1) = 0;
   coeffs(0, 2) = 0;
   coeffs(1, 0) = 0;
   coeffs(1, 1) = 0;
   coeffs(1, 2) = 0;
   coeffs(2, 0) = 0;
   coeffs(2, 1) = 0;
   coeffs(2, 2) = 0;
   consts(0) = 0;
   consts(1) = 0;
   consts(2) = 0;

   int N=0;

   for(unsigned i = 0; i < etas_.size(); i++)  {
       // etaPow.push_back(etas_[i]*etas_[i]);
     etaPow.push_back(fabs(etas_[i]));
     etaPow2.push_back(fabs(etas_[i])*fabs(etas_[i]) );
      
   }

      for(unsigned i = 0; i < etas_.size(); ++i)
      {  
         if ( fabs(ecalEnergies_[i] + hcalEnergies_[i] - ETrueEnergies_[i] + a_[i]) > sigC_*sigmaEcalHcal_[i]) continue;             
         if(hcalEnergies_[i] == 0)continue;
         
	 N++;
	 
         coeffs(0, 0) += 2.0*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
	 coeffs(0, 1) += 2.0*etaPow[i]*
	   (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
	   (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
	   (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

         coeffs(0, 2) += 2.0*etaPow2[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 0) += 2.0*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

	 coeffs(1, 1) += 2.0*etaPow2[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(1, 2) += 2.0*etaPow2[i]*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

	 coeffs(2, 0) += 2.0*etaPow2[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

	 coeffs(2, 1) += 2.0*etaPow2[i]*etaPow[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
         
         coeffs(2, 2) += 2.0*etaPow2[i]*etaPow2[i]*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
	 

         consts(0) += 2.0*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);  
         
	 consts(1) += 2.0*etaPow[i]*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

         consts(2) += 2.0*etaPow2[i]*
            (ETrueEnergies_[i] - a_[i] - ecalEnergies_[i] - 
             hcalEnergies_[i])*
            (factorE*ecalEnergies_[i] + hcalEnergies_[i])/
            (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      }
   

   //Create the matrix that will be inverted and the vector which will multiply
   //that matrix to find the alpha and beta calibration constants.
      //    cout<<" N= "<<N<<endl;
      // coeffs.Print(cout);
      if(N==0) return false;
      isInverted = coeffs.Invert();
      // coeffs.Print(cout); cout<<endl;
      // consts.Print(cout);
      
   if(isInverted && sqrt(coeffs(0,0)) < 100000)
   {
      values = coeffs*consts;
      
      alpha_ = values(0);
      gamma_ = values(1);
      beta_ = values(2);

      //      cout<<endl<<alpha_<<"   "<<gamma_<<"   "<<beta_<<endl;

      sigmaAlpha_ = sqrt(coeffs(0,0));
      sigmaGamma_ = sqrt(coeffs(1,1));
      sigmaBeta_ = sqrt(coeffs(2,2));
     

      return true;
   }
   else
   {
      return false;
   }
} 





void AlphaBeta::clear()
{
   ETrueEnergies_.clear();
   ecalEnergies_.clear();
   hcalEnergies_.clear();
   etas_.clear();
   
   delete &binLowEdge_;
   delete &binHighEdge_;
   delete &isBarrel_;
   delete &ETrueAverage_;
   delete &ETrueRMS_;
   delete &a_;
   delete &alpha_;
   delete &beta_;
   delete &sigmaAlpha_;
   delete &sigmaBeta_;
}

///////////////////////////////////////////////////////////////////////////////
//The class that holds the calibration information over all values of ETrue. 
//It takes in a collection of ABC objects and fits each calibration 
//constant to a function. These functions are then used to find the calibrated
//energy.
///////////////////////////////////////////////////////////////////////////////
class Calibration
{
   private:
      double ETrueMax_;
      bool isBarrel_;
      
      vector<double> ETrueMeansABC_;
      vector<double> ETrueRMSsABC_;
      vector<double> ETrueMeansAlphaBeta_;
      vector<double> ETrueRMSsAlphaBeta_;

      vector<double> as_;
      vector<double> bs_;
      vector<double> cs_;
      vector<double> alphas_;
      vector<double> betas_;
  vector<double> gammas_;
 
      vector<double> sigmaBs_;
      vector<double> sigmaCs_;
      vector<double> sigmaAlphas_;
      vector<double> sigmaBetas_;
  vector<double> sigmaGammas_;
      
      TGraph *graphA_ ;
      TGraphErrors *graphB_;
      TGraphErrors *graphC_;
      TGraphErrors *graphAlpha_;
      TGraphErrors *graphBeta_ ;
  TGraphErrors *graphGamma_ ;

      TGraph *graphBError_;
      TGraph *graphCError_;
      TGraph *graphAlphaError_;
      TGraph *graphBetaError_ ;
      TF1* functionA_;
      TF1* functionB_;
      TF1* functionC_;
      TF1* functionAlpha_;
      TF1* functionBeta_;
  TF1* functionGamma_;
      
   public:
      Calibration();
      Calibration(double ETrueMax, bool isBarrel);
      //Adds a graph point by hand, i.e. putting in the individual values.
      void addGraphPoints(double ETrueAverage, double ETrueRMS,
                          double a, double b, double sigmaB, double C,
                          double sigmaC, double alpha, double sigmaAlpha,
                          double beta, double sigmaBeta,
			  double gamma, double sigmaGamma ); 
      //Adds a graph point by taking apart an ABC object.
      void addGraphPoints(ABC* abc);
      void addGraphPoints( AlphaBeta* alphabeta);
  //Creates the graphs after the points have been added.
  void initializeGraphs(string option);  
  double getETrueMax();      //Returns the various objects  that the 
  TGraphErrors* getGraph();  //calibration class holds.
  TF1* getFunctionA();
  TF1* getFunctionB();
  TF1* getFunctionC();
  TF1* getFunctionAlpha();
  TF1* getFunctionBeta();   
  TF1* getFunctionGamma();   
  //Returns calibrated energy without any eta dependence.
  double getCalibratedEnergy(double ETrue, double ecalEnergy, 
			     double hcalEnergy); 
  //Returns calibrated energy with eta dependence.
  double getCalibratedEnergy(double ETrue, double ecalEnergy, 
			     double hcalEnergy, double eta);
  //and with advanced pol dependency
  double getCalibratedEnergyWithGamma(double ETrue, double ecalEnergy, 
				      double hcalEnergy, double eta);
    
    void setETrueMax(double ETrueMax);
  bool fitAsToFunction(TF1 *functionA);  //Fits the functions to their 
  bool fitAsToFunction();                //graph points. One that takes in
  bool fitBsToFunction(TF1 *functionB);  //a function which then fits it to
  bool fitBsToFunction();                //that function, and one that is
      bool fitCsToFunction(TF1 *functionC);  //used to simply improve the fit.
      bool fitCsToFunction();
      bool fitAlphasToFunction(TF1 *functionAlpha);
      bool fitAlphasToFunction();
      bool fitBetasToFunction(TF1 *functionBetas);
      bool fitBetasToFunction();
  bool fitGammasToFunction(TF1 *functionGammas);
  bool fitGammasToFunction();
  
  bool setAlphasToFunction(TF1 *functionAlpha);
  bool setBetasToFunction(TF1 *functionBeta);
  
      void drawCoeffGraph(string graph, string tag);  //Makes and draws a graph of a 
      void drawSigmaGraph(string graph);  //coefficient. Takes in as an 
                                          //argument a, b, c, alpha, or beta.
      void printBs();
      void printCs();
      void printBetas();
      void printAlphas();
  void printGammas();
};

Calibration::Calibration()
{
   ETrueMax_ = 1000;
   isBarrel_ = true;
}
Calibration::Calibration(double ETrueMax, bool isBarrel) 
          
{
   ETrueMax_ = ETrueMax;
   isBarrel_ = isBarrel;
}
void Calibration::addGraphPoints(double ETrueAverage, double ETrueRMS,
                                 double a, double b, double sigmaB, double c,
                                 double sigmaC, double alpha, 
                                 double sigmaAlpha, double beta,
				 double sigmaBeta, double gamma,
				 double sigmaGamma)
{
   ETrueMeansABC_.push_back(ETrueAverage);
   ETrueRMSsABC_.push_back(ETrueRMS);
   as_.push_back(a);
   bs_.push_back(b);
   cs_.push_back(c);
   alphas_.push_back(alpha);
   betas_.push_back(beta);
   gammas_.push_back(gamma);
   

   sigmaBs_.push_back(sigmaB);
   sigmaCs_.push_back(sigmaC);
   sigmaAlphas_.push_back(sigmaAlpha);
   sigmaBetas_.push_back(sigmaBeta);
   sigmaGammas_.push_back(sigmaGamma);

}
void Calibration::addGraphPoints(ABC* abc)
{
   if(abc->isEmpty() || (abc->getSigmaC() == 0 && abc->getC() == 0)) return;

   ETrueMeansABC_.push_back(abc->getETrueAverage());
   ETrueRMSsABC_.push_back(abc->getETrueRMS());
   as_.push_back(abc->getA());
   bs_.push_back(abc->getB());
   cs_.push_back(abc->getC());


   sigmaBs_.push_back(abc->getSigmaB());
   sigmaCs_.push_back(abc->getSigmaC());
}

void Calibration::addGraphPoints(AlphaBeta* alphabeta)
{
   if(alphabeta->isEmpty()|| 
      (alphabeta->getSigmaBeta() == 0 && alphabeta->getBeta() == 0)) return;


   ETrueMeansAlphaBeta_.push_back(alphabeta->getETrueAverage());
   ETrueRMSsAlphaBeta_.push_back(alphabeta->getETrueRMS());
   alphas_.push_back(alphabeta->getAlpha());
   betas_.push_back(alphabeta->getBeta());
   gammas_.push_back(alphabeta->getGamma());
   sigmaAlphas_.push_back(alphabeta->getSigmaAlpha());
   sigmaBetas_.push_back(alphabeta->getSigmaBeta());
   sigmaGammas_.push_back(alphabeta->getSigmaGamma());

}
void Calibration::initializeGraphs(string option)
{
   vector<double> x;
   vector<double> sigmaX;
   vector<double> y;
   vector<double> sigmaY;
   
   if(option == "abc" || option == "ABC" || option == "all")
   {
      for(unsigned i = 0; i < ETrueMeansABC_.size(); i++)
      {
         if( bs_[i] == 0 && sigmaBs_[i] == 0.0) continue;
         
         x.push_back(ETrueMeansABC_[i]);
         sigmaX.push_back(ETrueRMSsABC_[i]);
         y.push_back(bs_[i]);
         sigmaY.push_back(sigmaBs_[i]);
         
      }
      
      graphB_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                 &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();
      
      for(unsigned i = 0; i < ETrueMeansABC_.size(); i++)
      {
         if( cs_[i] == 0 && sigmaCs_[i] == 0.0) continue;
         
         x.push_back(ETrueMeansABC_[i]);
         sigmaX.push_back(ETrueRMSsABC_[i]);
         y.push_back(cs_[i]);
         sigmaY.push_back(sigmaCs_[i]);
         
      }
      
      graphC_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                 &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();
   }
   if(option == "AlphaBeta" || option == "alphabeta" || 
      option == "alphaBeta" || option == "ALPHABETA" || option == "all")
   {      
      for(unsigned i = 0; i < ETrueMeansAlphaBeta_.size(); i++)
      {
         if( alphas_[i] == 0 && sigmaAlphas_[i] == 0.0) continue;
         
         x.push_back(ETrueMeansAlphaBeta_[i]);
         sigmaX.push_back(ETrueRMSsAlphaBeta_[i]);
         y.push_back(alphas_[i]);
         sigmaY.push_back(sigmaAlphas_[i]);
         
      }

      graphAlpha_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                  &sigmaY[0]);
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();
      
      for(unsigned i = 0; i < ETrueMeansAlphaBeta_.size(); i++)
      {
         if( betas_[i] == 0 && sigmaBetas_[i] == 0.0) continue;

         x.push_back(ETrueMeansAlphaBeta_[i]);
         sigmaX.push_back(ETrueRMSsAlphaBeta_[i]);
         y.push_back(betas_[i]);
         sigmaY.push_back(sigmaBetas_[i]);
         
      }
   
      graphBeta_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                    &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();

  for(unsigned i = 0; i < ETrueMeansAlphaBeta_.size(); i++)
      {
         if( betas_[i] == 0 && sigmaGammas_[i] == 0.0) continue;

         x.push_back(ETrueMeansAlphaBeta_[i]);
         sigmaX.push_back(ETrueRMSsAlphaBeta_[i]);
         y.push_back(gammas_[i]);
         sigmaY.push_back(sigmaGammas_[i]);
         
      }
   
      graphGamma_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                    &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();



   }
   
}
double Calibration::getETrueMax() {return ETrueMax_;}
TGraphErrors* Calibration::getGraph() {return graphB_;}
TF1* Calibration::getFunctionA() {return functionA_;}
TF1* Calibration::getFunctionB() {return functionB_;}
TF1* Calibration::getFunctionC() {return functionC_;}
TF1* Calibration::getFunctionAlpha() {return functionAlpha_;}
TF1* Calibration::getFunctionBeta() {return functionBeta_;}
TF1* Calibration::getFunctionGamma() {return functionGamma_;}

double Calibration::getCalibratedEnergyWithGamma(double ETrue, double ecalEnergy, 
						 double hcalEnergy, double eta)
{
   // double etaPow;
   // double etaPow2;

  //spandey change 3 oct, 2017
   double etaPow = 0;
   double etaPow2 = 0;

   double factor_;
   double a = functionA_->Eval(ETrue);
   double b = functionB_->Eval(ETrue);
   double c = functionC_->Eval(ETrue);
   double alpha = functionAlpha_->Eval(ETrue);
   double beta = functionBeta_->Eval(ETrue);
   double gamma = functionGamma_->Eval(ETrue);
   double counterAlpha = 0;
   double counterBeta = 0;

   if(isBarrel_) 
   {
      etaPow = eta*eta;
      factor_ = factorB;
      if(ecalEnergy>0) {
	counterAlpha = alpha;
	counterBeta = beta;
      }
      gamma = 0;
   }
   else 
   {
     etaPow = fabs(eta);
     etaPow2 = fabs(eta)*fabs(eta) ;
     // etaPow = eta*eta;
     factor_ = factorE;
   }

   // if(fabs(eta)>1.5 && ecalEnergy==0)
   //   cout<<alpha<<"  "<<beta<<"   "<<ETrue<<endl;
   
   return a + (1.0 + alpha + gamma*etaPow + factor_*beta*etaPow2)*b*ecalEnergy + 
     (1.0 + alpha + gamma*etaPow + beta*etaPow2 - counterAlpha - counterBeta*etaPow)*
     c*hcalEnergy;
}



void Calibration::setETrueMax(double ETrueMax){ETrueMax_ = ETrueMax;}


bool Calibration::fitAsToFunction(TF1 *functionA)
{
   functionA_ = functionA;
//   graphA_->Fit(functionA_->GetName(), "Q", "", 0, ETrueMax_);
   return true;
}
bool Calibration::fitAsToFunction()
{
   graphA_->Fit(functionA_->GetName(), "Q", "", 2., ETrueMax_);
   return true;
}

bool Calibration::fitBsToFunction(TF1 *functionB)
{
   functionB_ = functionB;
   graphB_->Fit(functionB_->GetName(), "Q", "", 2., ETrueMax_);
   return true;
}
bool Calibration::fitBsToFunction()
{
   graphB_->Fit(functionB_->GetName(), "Q", "", 2., ETrueMax_);
   return true;
}
bool Calibration::fitCsToFunction(TF1 *functionC)
{
   functionC_ = functionC;
   graphC_->Fit(functionC_->GetName(), "Q", "", 2., ETrueMax_);

   return true;
}
bool Calibration::fitCsToFunction()
{
   graphC_->Fit(functionC_->GetName(), "Q", "", 2., ETrueMax_);
   return true;
}

bool Calibration::setAlphasToFunction(TF1 *functionAlpha)
{
   functionAlpha_ = functionAlpha;
   return true;
}

bool Calibration::fitAlphasToFunction(TF1 *functionAlpha)
{
   functionAlpha_ = functionAlpha;
   graphAlpha_->Fit(functionAlpha_->GetName(), "Q", "", 2., ETrueMax_); //2.
   return true;
}
bool Calibration::fitAlphasToFunction()
{
  graphAlpha_->Fit(functionAlpha_->GetName(), "Q", "", 2., ETrueMax_); //2.
   return true;
}

bool Calibration::setBetasToFunction(TF1 *functionBeta)
{
  functionBeta_ = functionBeta;
   return true;
}

bool Calibration::fitBetasToFunction(TF1 *functionBeta)
{
   functionBeta_ = functionBeta;
   graphBeta_->Fit(functionBeta_->GetName(), "Q", "", 2.0, ETrueMax_);
   return true;
}
bool Calibration::fitBetasToFunction()
{
   graphBeta_->Fit(functionBeta_->GetName(), "Q", "", 2.0, ETrueMax_);
   return true;
}

bool Calibration::fitGammasToFunction(TF1 *functionGamma)
{
   functionGamma_ = functionGamma;
   graphGamma_->Fit(functionGamma_->GetName(), "Q", "", 2.0, ETrueMax_ );
   return true;
}
bool Calibration::fitGammasToFunction()
{
   graphGamma_->Fit(functionGamma_->GetName(), "Q", "", 2.0, ETrueMax_);
   return true;
}



void Calibration::drawCoeffGraph(string graph, string tag)
{
  

  cout<<" Tag = "<<tag<<endl;
  cout<<" Graph = "<<graph<<endl;
  TString fileName="resp_reso_"+graph+"_"+tag+".root";
  //char* fileName = new char[1000];
   string saveString;
   //TCanvas* canvas = new TCanvas( (graph+tag).c_str() , graph.c_str(), 1200,600 );
   TCanvas* canvas = new TCanvas( (graph+tag).c_str() , graph.c_str(), 500,300 );
   //   sprintf(fileName,"resp_reso_%s_%s.root",graph,tag);
   TFile* file3=new TFile(fileName,"recreate");
   file3->cd();
   TH2F* histo = new TH2F("histoCG", "", sampleRangeHigh, 0, sampleRangeHigh, sampleRangeHigh,  -3.0, 3.0); 

   canvas->cd();
   canvas->SetLogx();
   histo->SetStats(0);
   histo->GetXaxis()->SetTitle("True Energy(GeV)");
   histo->GetXaxis()->SetTitle("True Energy(GeV)");
   histo->GetYaxis()->SetTitle("Value");
   histo->Draw();
   
   gPad->SetGridx();
   gPad->SetGridy();
   
   // TLegend* leg=new TLegend(0.53,0.23,0.91,0.44);
   // leg->SetFillColor(0); leg->SetShadowColor(0);
   // TLine* line=new TLine(0,0,0,0);
   if(graph == "a" || graph == "A") 
     {
       histo->SetTitle("A Parameter ");
       graphB_->SetTitle("A vs True Energy");
       histo->GetYaxis()->SetTitle( ("A coefficient ("+tag+")").c_str() );
       graphB_->SetMarkerStyle(22);
       graphB_->SetMarkerSize(1);
       graphB_->SetFillColor(0);

       graphB_->Draw("P");
       // if(tag=="EH") {
       // 	faBarrel->Draw("Lsame+");
       //   	faBarrel52x->Draw("Lsame+");
       // }

       //  leg->AddEntry(
     TLegend *leg=new TLegend(0.30,0.25,0.90,0.35);
     // leg->AddEntry(histo,"[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))","");
     leg->AddEntry(histo,"[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3]))))","");//for UL 2016 ec
     // leg->AddEntry(histo,"[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))","");//for UL 2016 barrel 
     leg->SetTextAlign(32);
     leg->SetTextSize(0.04);
     leg->Draw();
     saveString = "ACoefficient" + tag + ".gif";
     canvas->SaveAs(saveString.c_str());
     saveString = "ACoefficient" + tag + ".png";
     canvas->SaveAs(saveString.c_str());
     graphB_->Write();

     }
   if(graph == "b" || graph == "B") 
     {
       histo->SetTitle("B Parameter");
       graphC_->SetTitle("B vs True Energy");
       histo->GetYaxis()->SetTitle( ("B coefficient ("+tag+")").c_str() );
       graphC_->SetMarkerStyle(22);
       graphC_->SetMarkerSize(1);
       graphC_->SetFillColor(0);

      graphC_->Draw("P");
      // if(tag=="EH") {
      // 	faBarrel->Draw("Lsame+");
	//   	faBarrel52x->Draw("Lsame+");
      //}

      //  leg->AddEntry(
     TLegend *leg=new TLegend(0.30,0.25,0.90,0.35);
     leg->AddEntry(histo,"[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))","");//for UL2016 barrel(EH)
     //     leg->AddEntry(histo,"[0]+([4]*(x-[5])*exp(-(x*[7])))+(([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))",""); //for UL2016 endcap(EH)
     leg->SetTextAlign(32);
     leg->SetTextSize(0.04);
     leg->Draw();
     saveString = "BCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
     saveString = "BCoefficient" + tag + ".png";
      canvas->SaveAs(saveString.c_str());

      graphC_->Write();

     }
   else if(graph == "c" || graph == "C") 
   {
      histo->SetTitle("C parameter ");
      graphC_->SetTitle("C vs True Energy");
      histo->GetYaxis()->SetTitle( ("C coefficient ("+tag+")").c_str() );
      graphC_->SetMarkerStyle(22);
      graphC_->SetMarkerSize(1);
      //graphC_->SetMarkerColor(1);
      graphC_->SetFillColor(0);
      
      graphC_->Draw("P");

     TLegend *leg=new TLegend(0.30,0.25,0.90,0.35);
     //     leg->AddEntry(histo,"[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",""); //for UL 2017/2016 barrel
     leg->AddEntry(histo,"[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",""); //for UL 2017/2016 endcap
     leg->SetTextAlign(32);
     leg->SetTextSize(0.04);
     leg->Draw();
      // if(tag=="EH") {
      // 	fbBarrel->Draw("Lsame+");
      // 	//   	fbBarrel52x->Draw("Lsame+");
      // }
      // else {
      // 	fcBarrel->Draw("Lsame+");
      // 	//    	fcBarrel52x->Draw("Lsame+");
      // }

      saveString = "CCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
      saveString = "CCoefficient" + tag + ".png";
      canvas->SaveAs(saveString.c_str());

      graphC_->Write();

   }
   else if(graph == "alpha" || graph == "Alpha") 
   {
      histo->SetTitle("Alpha parameter ");
      graphAlpha_->SetTitle("Alpha vs True Energy");
      histo->GetYaxis()->SetTitle( ("#alpha coefficient ("+tag+")").c_str() );
      graphAlpha_->SetMarkerStyle(22);
      graphAlpha_->SetMarkerSize(1);
      graphAlpha_->SetMarkerColor(2);
      graphAlpha_->SetFillColor(0);

      graphAlpha_->Draw("P");
      faEtaBarrel->Draw("Lsame+");
      //   faEtaBarrel52x->Draw("Lsame+");
     TLegend *leg=new TLegend(0.30,0.75,0.85,0.85);
     //leg->AddEntry(graphAlpha_,"[0]+[1]*exp(-x/[2])",""); //for UL2017 endcap
     //leg->AddEntry(graphAlpha_,"[0]+[1]*x^[3]*exp(-x/[2])",""); //for UL2016 endcap H
     //     leg->AddEntry(graphAlpha_,"[0]+([1]*x^[2]*exp(-x))","");//for 2016 endcap EH 
     //leg->AddEntry(histo,"[0]+[1]*x",""); //for UL 2016/2017 barrel
     if(tag=="EH_endcap" || tag=="H_endcap") leg->AddEntry(histo,"[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))","");
     if(tag=="EH_barrel") leg->AddEntry(histo,"[0]+[1]*exp(-x*[3]/[2])","");
     if(tag=="H_barrel") leg->AddEntry(histo,"[0]+[1]*x","");
     leg->SetTextSize(0.04);
     leg->Draw();

      saveString = "AlphaCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
      graphAlpha_->Write();
  }
   else if(graph == "beta" || graph == "Beta") 
   {
      histo->SetTitle("Beta parameter ");
      graphBeta_->SetTitle("Beta vs True Energy");
      histo->GetYaxis()->SetTitle( ("#beta coefficient ("+tag+")").c_str() );
      graphBeta_->SetMarkerStyle(22);
      graphBeta_->SetMarkerSize(1);
      graphBeta_->SetMarkerColor(2);
      graphBeta_->SetFillColor(0);
      
      graphBeta_->Draw("P");
      fbEtaBarrel->Draw("Lsame+");
      //    fbEtaBarrel52x->Draw("Lsame+");
     TLegend *leg=new TLegend(0.30,0.75,0.85,0.85);
     // leg->AddEntry(histo,"[0]+[1]*exp(-x/[2])",""); //for UL2017 endcap/barrel & for UL2016 barrel
     //     leg->AddEntry(histo,"[0]+[1]*x*exp(-x/[2])",""); //for UL 2016 endcap H 
     //     leg->AddEntry(histo,"[0]+[1]*(x^[3])*exp(-x/[2])",""); //for UL 2016 endcap EH
     if(tag=="EH_endcap") leg->AddEntry(histo,"[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))","");
     if(tag=="EH_endcap") leg->AddEntry(histo,"[0]+[1]*x*exp(-x/[2])","");
     if(tag=="EH_barrel") leg->AddEntry(histo,"[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))","");
     if(tag=="H_barrel") leg->AddEntry(histo,"[0]+[1]*exp(-x/[2])","");
     leg->SetTextSize(0.04);
     leg->Draw();
     saveString = "BetaCoefficient" + tag + ".gif";
     canvas->SaveAs(saveString.c_str());
     saveString = "BetaCoefficient" + tag + ".png";
     canvas->SaveAs(saveString.c_str());

      graphBeta_->Write();
   }  
 else if(graph == "gamma" || graph == "Gamma") 
   {
      histo->SetTitle("Gamma parameter");
      graphGamma_->SetTitle("Gamma vs True Energy");
      histo->GetYaxis()->SetTitle( ("#gamma coefficient ("+tag+")").c_str() );
      graphGamma_->SetMarkerStyle(22);
      graphGamma_->SetMarkerSize(1);
      graphGamma_->SetMarkerColor(2);
      graphGamma_->SetFillColor(0);
      
      graphGamma_->Draw("P");
      // fbEtaBarrel->Draw("Lsame+");
      // fbEtaBarrel52x->Draw("Lsame+");

      saveString = "GammaCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
      saveString = "GammaCoefficient" + tag + ".png";
      canvas->SaveAs(saveString.c_str());

   }  



   else cout << "No graph with that name" <<endl;
 
  TLine* line=new TLine(1,2,1,2);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);

 
  // leg->AddEntry(line,"coef 52X (new PF Ecal cluster calib)","l");
   file3->cd();
   file3->Write();
   file3->Close();

}
void Calibration::drawSigmaGraph(string graph)
{
  
   TCanvas* canvas2 = new TCanvas();
   TH2F* histo2 = new TH2F("histoS", "", 100, 0, ETrueMax_, 100,  0, .1); 
   
   canvas2->cd();
   histo2->SetStats(0);
   histo2->Draw();
   
   gPad->SetGridx();
   gPad->SetGridy();

   if(graph == "b" || graph =="B")
   {
      graphBError_->SetMarkerStyle(22);
      graphBError_->SetMarkerColor(2);
      graphBError_->SetMarkerSize(.5);
      graphBError_->Draw("P");

   }
   else if(graph == "c" || graph =="C")
   {
      graphCError_->SetMarkerStyle(22);
      graphCError_->SetMarkerColor(2);
      graphCError_->SetMarkerSize(.5);
      graphCError_->Draw("P");
   }
   else if(graph == "alpha" || graph =="Alpha")
   {
      graphAlphaError_->SetMarkerStyle(22);
      graphAlphaError_->SetMarkerColor(2);
      graphAlphaError_->SetMarkerSize(.5);
      graphAlphaError_->Draw("P");
   }
   else if(graph == "beta" || graph =="Beta")
   {
      graphBetaError_->SetMarkerStyle(22);
      graphBetaError_->SetMarkerColor(2);
      graphBetaError_->SetMarkerSize(.5);
      graphBetaError_->Draw("P");
   }
   else cout<<"No graph with that name"<<endl;
   
   delete histo2;

}
void Calibration::printBs()
{
   
   for(unsigned i = 0; i < as_.size(); i++)
   {
     cout<<bs_[i];
   }
}

void Calibration::printCs()
{
   
   for(unsigned i = 0; i < as_.size(); i++)
   {
      cout<<cs_[i]<<endl;
   }
}
void Calibration::printBetas()
{
   
   for(unsigned i = 0; i < betas_.size(); i++)
   {
      cout<<betas_[i];
   }
}
void Calibration::printAlphas()
{
   
   for(unsigned i = 0; i < alphas_.size(); i++)
   {
      cout<<"Alphas: "<<alphas_[i]<<endl;;
   }
}

void Calibration::printGammas()
{
   
   for(unsigned i = 0; i < gammas_.size(); i++)
   {
      cout<<"Gammas: "<<gammas_[i]<<endl;;
   }
}

///////////////////////////////////////////////////////////////////////////////
//All the needed variables for the main (calibChris) function.
///////////////////////////////////////////////////////////////////////////////

TFile* inputFile;
TTree* sTree;

vector<double> ETrueEnergies;  //The values that are taken from the root file
vector<double> ecalEnergies;
vector<double> hcalEnergies;
vector<double> etas;
vector<double> phis;


vector<ABC*> barrelABCEcalHcal; //Vectors of the ABC objects
vector<ABC*> barrelABCEcal;     //which hold all the calibration 
vector<ABC*> barrelABCHcal;     //constants for an individual bin.
vector<ABC*> endcapABCEcalHcal;
vector<ABC*> endcapABCEcal;
vector<ABC*> endcapABCHcal;

vector<AlphaBeta*> barrelAlphaBetaEcalHcal; //Vectors of the AlphaBeta objects which 
vector<AlphaBeta*> barrelAlphaBetaHcal;
vector<AlphaBeta*> endcapAlphaBetaEcalHcal; //hold all the calibration constants for 
vector<AlphaBeta*> endcapAlphaBetaHcal; //an individual bin.

TF1* functionBarrelEcalHcalA;     //Functions that the calibration equations
TF1* functionBarrelEcalHcalB;     //are fit to
TF1* functionBarrelEcalHcalC;  
TF1* functionEndcapEcalHcalA;
TF1* functionEndcapEcalHcalB;
TF1* functionEndcapEcalHcalC;  

TF1* functionBarrelHcalA;
TF1* functionBarrelHcalB;
TF1* functionBarrelHcalC;  
TF1* functionEndcapHcalA;
TF1* functionEndcapHcalB;
TF1* functionEndcapHcalC;  

TF1* functionBarrelAlphaEcalHcal;
TF1* functionBarrelBetaEcalHcal;
TF1* functionBarrelAlphaHcal;
TF1* functionBarrelBetaHcal;

TF1* functionEndcapAlphaEcalHcal;
TF1* functionEndcapBetaEcalHcal;
TF1* functionEndcapAlphaHcal;
TF1* functionEndcapBetaHcal;
TF1* functionEndcapGammaHcal;

//Calibration objects which hold the all the calibration costants as functions
//of ETrue. 
Calibration* barrelWithEcalHcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* barrelWithEcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* barrelWithHcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* endcapWithEcalHcalCalib = new Calibration(sampleRangeHigh, false);
Calibration* endcapWithEcalCalib = new Calibration(sampleRangeHigh, false);
Calibration* endcapWithHcalCalib = new Calibration(sampleRangeHigh, false);

//Temporary varibles that will be used in for loops just to cut down on the 
//length of the lines of code.
double etrue;
double ecal;
double hcal;
double eta;
double bpar;
double cpar;
double etrueMax;
double barrelEcalHcalB;
double barrelEcalHcalC;
double barrelHcalC;
double barrelAlpha;
double barrelBeta;
double endcapEcalHcalB;
double endcapEcalHcalC;
double endcapHcalC;
double endcapAlpha;
double endcapBeta;
double correctedE;
double correctedEta;

const char* functionEndcapEcalHcalB_e;
const char* functionEndcapEcalHcalC_e;
const char* functionEndcapHcalC_e;
const char* functionEndcapAlphaEH_e;
const char* functionEndcapBetaEH_e;
const char* functionEndcapAlphaH_e;
const char* functionEndcapBetaH_e;
const char* functionBarrelEcalHcalB_e;
const char* functionBarrelEcalHcalC_e;
const char* functionBarrelHcalC_e;
const char* functionBarrelAlphaEH_e;
const char* functionBarrelBetaEH_e;
const char* functionBarrelAlphaH_e;
const char* functionBarrelBetaH_e;

//All the differenct types of TH2's that will be filled in order to make 
//resolution and response plots

TH2F* h_trueE_vs_mod_eta_response = new TH2F("h_trueE_vs_mod_eta_response","True Energy vs eta distribution",6,0.0,3.0,33,2,sampleRangeHigh);
TH2F* h_trueE_vs_mod_eta_response_normalized = new TH2F("h_trueE_vs_mod_eta_response_normalized","True Energy vs eta distribution with response on z axis (EH-hadrons)",6,0.0,3.0,33,2,sampleRangeHigh);

TH2F * h_response_vs_phi_barrel_EH = new TH2F("h_response_vs_phi_barrel_EH", "EH hadrons Response vs #phi in Barrel", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham
TH2F * h_response_vs_phi_EndCap_EH_posZ = new TH2F("h_response_vs_phi_EndCap_EH_posZ", "EH hadrons Response vs #phi for +z in EndCap", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham
TH2F * h_response_vs_phi_EndCap_EH_negZ = new TH2F("h_response_vs_phi_EndCap_EH_negZ", "EH hadrons Response vs #phi for -z in EndCap", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham

TH2F * h_response_vs_phi_barrel_H = new TH2F("h_response_vs_phi_barrel_H", "H hadrons Response vs #phi in Barrel", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham
TH2F * h_response_vs_phi_EndCap_H_posZ = new TH2F("h_response_vs_phi_EndCap_H_posZ", "H hadrons Response vs #phi for +z in EndCap", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham
TH2F * h_response_vs_phi_EndCap_H_negZ = new TH2F("h_response_vs_phi_EndCap_H_negZ", "H hadrons Response vs #phi for -z in EndCap", 70, -3.5, 3.5, 30, -1.5, 1.5); //shubham


TH2F* raw = new TH2F("raw","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEta = new TH2F("corrEta", "", 1000, 0, 1000, 150, -1.5, 1.5);

TH2F* rawBarrel = new TH2F("rawBarrel","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrBarrel = new TH2F("corrBarrel", "", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEtaBarrel = new TH2F("corrEtaBarrel", "", 1000, 0, 1000, 150, -1.5, 
                               1.5);
TH2F* rawBarrelEcalHcal = new TH2F("rawBarrelEcalHcal","", 1000, 0, 1000, 150, 
                                   -1.5, 1.5);
TH2F* corrBarrelEcalHcal = new TH2F("corrBarrelEcalHcal", "", 1000, 0, 1000, 
                                    150, -1.5, 1.5);
TH2F* corrEtaBarrelEcalHcal = new TH2F("corrEtaBarrelEcalHcal", "", 1000, 0, 
                                       1000, 150, -1.5, 1.5);
TH2F* rawBarrelHcal = new TH2F("rawBarrelHcal","", 1000, 0, 1000, 150, -1.5, 1.5 );
			       //1000, 0.0, 5.0);
TH2F* corrBarrelHcal = new TH2F("corrBarrelHcal", "", 1000, 0, 1000, 150, -1.5,
                                1.5);
TH2F* corrEtaBarrelHcal = new TH2F("corrEtaBarrelHcal", "", 1000, 0, 1000, 150,
                                   -1.5, 1.5);

TH2F* rawEndcap = new TH2F("rawEndcap","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEndcap = new TH2F("corrEndcap", "", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEtaEndcap = new TH2F("corrEtaEndcap", "", 1000, 0, 1000, 150, -1.5, 
                               1.5);
//TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 150, -1.5, 1.5);
//TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 500, -1.5, 10.0);
TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 575, -1.5, 10.0);

TH2F* corrEndcapEcalHcal = new TH2F("corrEndcapEcalHcal", "", 1000, 0, 1000, 
                                    150, -1.5, 1.5);
// TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 
//                                        1000, 150, -1.5, 1.5);
//TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 1000, 500, -1.5, 10.0);
TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 1000, 575, -1.5, 10.0);

TH2F* rawEndcapHcal = new TH2F("rawEndcapHcal","", 1000, 0, 1000, 150, -1.5, 
                               1.5);
TH2F* corrEndcapHcal = new TH2F("corrEndcapHcal", "", 1000, 0, 1000, 150, -1.5,
                                1.5);
TH2F* corrEtaEndcapHcal = new TH2F("corrEtaEndcapHcal", "", 1000, 0, 1000, 150,-1.5, 1.5);


TH2F * rawEtaDependenceEH = new TH2F("rawEtaDependenceEH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * corrEtaDependenceEH = new TH2F("corrEtaDependenceEH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * hcorrEtaDependenceEH = new TH2F("hcorrEtaDependenceEH","Response vs. Eta", 75, 0., 3.0, 150, -1.0,1.0 );

TH2F * rawEtaDependenceH = new TH2F("rawEtaDependenceH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * corrEtaDependenceH = new TH2F("corrEtaDependenceH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * hcorrEtaDependenceH = new TH2F("hcorrEtaDependenceH","Response vs. Eta", 75, 0., 3.0, 150, -1.0,1.0 );

TH1F * trueTempHisto = new TH1F("trueTempHisto", "true", sampleRangeHigh, 0, sampleRangeHigh);
TH1F * ecalTempHisto = new TH1F("ecalTempHisto", "ecal", sampleRangeHigh, 0, sampleRangeHigh);

TProfile2D* corrEtaDependenceProfEH=new TProfile2D("EH","", 40, 0,100, 75, 0, 3.0);
TProfile2D* corrEtaDependenceProfH=new TProfile2D("H","", 40, 0,100, 75, 0, 3.0);

TH2F* bcplot= new TH2F("bcplot","bcplot",1000,-1.,1.,1000,-1.,1.);

//Temporary TGraphs to passed drawGausFit
TGraph response;
TGraph resolution;
TGraph responseRaw;
TGraph resolutionRaw;
TGraph responseCor;
TGraph resolutionCor;
TGraph responseEta;
TGraph resolutionEta;
TGraph responseEtaEtaEH;
TGraph responseEtaHCorrEtaEH;
TGraph responseEtaEtaH;
TGraph responseEtaHCorrEtaH;

/// Bin Manager =========================================
vector<double> BinsETrue;
vector<double> BinsETrueEta;

unsigned int GetETrueBin(double etrue) {

  bool find=false;

  int bm=0;
  int bM=BinsETrue.size()-1;
  if(etrue< BinsETrue[ bm ]) return -1;
  if(etrue> BinsETrue[ bM ]) return -1;
  int n=0;
  while(!find) {
    
    if(etrue <= BinsETrue[ bm+(bM-bm)/2 ] ) {
      bM = bm+(bM-bm)/2;
    }
    else {
      bm = bm+(bM-bm)/2;
    }
    if( fabs(bm-bM)==1 )  {
      return bm;
    }
   
    if(n>(int)BinsETrue.size()) return -1;
    n++;
  }
  return -1;
}


unsigned int GetETrueBinEta(double etrue) {

  bool find=false;

  int bm=0;
  int bM=BinsETrueEta.size()-1;
 
  if(etrue< BinsETrueEta[ bm ]) return -1;
  if(etrue> BinsETrueEta[ bM ]) return -1;
  int n=0;
  while(!find) {
    
    if(etrue <= BinsETrueEta[ bm+(bM-bm)/2 ] ) {
      bM = bm+(bM-bm)/2;
    }
    else {
      bm = bm+(bM-bm)/2;
    }
    if( fabs(bm-bM)==1 )  {
      return bm;
    }
   
    if(n>(int)BinsETrueEta.size()) return -1;
    n++;
  }
  return -1;
}

///======================================================

#endif
