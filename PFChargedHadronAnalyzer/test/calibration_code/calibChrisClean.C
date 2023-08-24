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
#include "TGraphErrors.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <string>
#include <iostream>
#include <math.h>

//#include "CrystalBall.C"

using namespace std;

double sigC_ = 5.;

/////////////MM
double aEH = 3.5;
double aE = 3.5;
double aH = 2.5;//3.0;

double aEHe = 3.5;
double aEe = 3.5;
double aHe = 2.5;
/////////////MM

double lBs = 0.25; //0.25 B //0.5 E
double mBs=1.; //1 B
double hBs=5.0; //10 40 B 50 E

double RBE = 4; //rebinning for alpha-beta curves

//Eta factor : ecal*factor + hcal
double factorB = 1.3;//Optimal //A factor put in by hand to make the eta dependence agree
double factorE = 1.3; //better. Should fit for the value later.

vector<double> sigmas;

void LoadOldThresholds() {
  aEH = 3.5; aE = 3.5; aH = 2.5; //3.5 2.5
  aEHe = 3.5; aEe = 3.5; aHe = 2.5;
}


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
      etaMaxFit_ = 1.0;
      etaMinEtaFit_ = 0.0;
      etaMaxEtaFit_ = 1.3;
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
   cout<<totalEMinusATimesHcal<<"   "<<totalHcalSquared<<"   "<<c_<<endl;
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
      etaMaxFit_ = 1.2;
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

     for(unsigned i = 0; i < etas_.size(); i++) 
       //etaPow.push_back((fabs(etas_[i]))*(fabs(etas_[i]) ) *  (fabs(etas_[i]))*(fabs(etas_[i]) ) );
       etaPow.push_back((fabs(etas_[i]) - 1.5)*(fabs(etas_[i]) - 1.5)  
			*(fabs(etas_[i]) - 1.5)*(fabs(etas_[i]) - 1.5));
     
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
      coeffs.Print(cout); cout<<endl;
      consts.Print(cout);
      
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
double Calibration::getCalibratedEnergy(double ETrue, double ecalEnergy, 
                                        double hcalEnergy)
{
  double a = functionA_->Eval(ETrue);
  double b = functionB_->Eval(ETrue);
  double c = functionC_->Eval(ETrue);
  
   return a+ b*ecalEnergy + c*hcalEnergy;
}
double Calibration::getCalibratedEnergy(double ETrue, double ecalEnergy, 
                                        double hcalEnergy, double eta)
{
   double etaPow;
   double factor_;
   double a = functionA_->Eval(ETrue);
   double b = functionB_->Eval(ETrue);
   double c = functionC_->Eval(ETrue);
   double alpha = functionAlpha_->Eval(ETrue);
   double beta = functionBeta_->Eval(ETrue);
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
   
   }
   else 
   {
     //  etaPow = fabs(eta)*fabs(eta)*fabs(eta)*fabs(eta);
     etaPow = (fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);
      factor_ = factorE;
   }

   // if(fabs(eta)>1.5 && ecalEnergy==0)
   //   cout<<alpha<<"  "<<beta<<"   "<<ETrue<<endl;
   
   return a + (1.0 + alpha + factor_*beta*etaPow)*b*ecalEnergy + 
      (1.0 + alpha + beta*etaPow - counterAlpha - counterBeta*etaPow)*
      c*hcalEnergy;
}

double Calibration::getCalibratedEnergyWithGamma(double ETrue, double ecalEnergy, 
						 double hcalEnergy, double eta)
{
   double etaPow;
   double etaPow2;
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
   graphA_->Fit(functionA_->GetName(), "Q", "", 0, ETrueMax_);
   return true;
}

bool Calibration::fitBsToFunction(TF1 *functionB)
{
   functionB_ = functionB;
   graphB_->Fit(functionB_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}
bool Calibration::fitBsToFunction()
{
   graphB_->Fit(functionB_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}
bool Calibration::fitCsToFunction(TF1 *functionC)
{
   functionC_ = functionC;
   graphC_->Fit(functionC_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}
bool Calibration::fitCsToFunction()
{
   graphC_->Fit(functionC_->GetName(), "Q", "", 1.5, ETrueMax_);
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
   graphGamma_->Fit(functionGamma_->GetName(), "Q", "", 2.0, ETrueMax_);
   return true;
}
bool Calibration::fitGammasToFunction()
{
   graphGamma_->Fit(functionGamma_->GetName(), "Q", "", 2.0, ETrueMax_);
   return true;
}



void Calibration::drawCoeffGraph(string graph, string tag)
{
  

  //  cout<<" Tag = "<<tag<<endl;

   string saveString;
   TCanvas* canvas = new TCanvas( (graph+tag).c_str() , graph.c_str(), 1200,600 );
   TH2F* histo = new TH2F("histoCG", "", 200, 0, 200, 200,  -2.0, 2.0); 

   canvas->cd();
   histo->SetStats(0);
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
       histo->SetTitle("A Parameter");
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

       saveString = "ACoefficient" + tag + ".gif";
       canvas->SaveAs(saveString.c_str());
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

      saveString = "BCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
   }
   else if(graph == "c" || graph == "C") 
   {
      histo->SetTitle("C parameter");
      graphC_->SetTitle("C vs True Energy");
      histo->GetYaxis()->SetTitle( ("C coefficient ("+tag+")").c_str() );
      graphC_->SetMarkerStyle(22);
      graphC_->SetMarkerSize(1);
      //graphC_->SetMarkerColor(1);
      graphC_->SetFillColor(0);
      
      graphC_->Draw("P");
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
   }
   else if(graph == "alpha" || graph == "Alpha") 
   {
      histo->SetTitle("Alpha parameter");
      graphAlpha_->SetTitle("Alpha vs True Energy");
      histo->GetYaxis()->SetTitle( ("#alpha coefficient ("+tag+")").c_str() );
      graphAlpha_->SetMarkerStyle(22);
      graphAlpha_->SetMarkerSize(1);
      graphAlpha_->SetMarkerColor(2);
      graphAlpha_->SetFillColor(0);

      graphAlpha_->Draw("P");
      faEtaBarrel->Draw("Lsame+");
      //   faEtaBarrel52x->Draw("Lsame+");

      saveString = "AlphaCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
   }
   else if(graph == "beta" || graph == "Beta") 
   {
      histo->SetTitle("Beta parameter");
      graphBeta_->SetTitle("Beta vs True Energy");
      histo->GetYaxis()->SetTitle( ("#beta coefficient ("+tag+")").c_str() );
      graphBeta_->SetMarkerStyle(22);
      graphBeta_->SetMarkerSize(1);
      graphBeta_->SetMarkerColor(2);
      graphBeta_->SetFillColor(0);
      
      graphBeta_->Draw("P");
      fbEtaBarrel->Draw("Lsame+");
      //    fbEtaBarrel52x->Draw("Lsame+");

      saveString = "BetaCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
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
   }  



   else cout << "No graph with that name" <<endl;
 
  TLine* line=new TLine(1,2,1,2);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);

 
  // leg->AddEntry(line,"coef 52X (new PF Ecal cluster calib)","l");

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
//Global functions used in the main of the code.
///////////////////////////////////////////////////////////////////////////////


//Takes apart a TH2 and creates response and resolution plots from it. Note: 
//the TGraphs will not draw correctly without passing the TGraphs as references
//to the function.

void drawGausFit(TH2F* inHisto, TGraph& response, TGraph& resolution)
{
  cout<<"check point 1"<<endl;
   if(inHisto->GetEntries() == 0) return;
   
   vector<TH1F*> ETrueBin;
   TF1* gaus; 
   string name;
   char num[4];
   float rebin = 1.;
   TCanvas* canvas;
   TCanvas* temp = new TCanvas();
   TLine *line = new TLine(0.0,0.0,200,0.0);

   TH2F* respHisto = new TH2F("respHisto", "", 200, 0, 200, 100, -0.5, 0.5);
   TH2F* resoHisto = new TH2F("resoHisto", "", 200, 0, 200, 100, 0.0, 0.5);

   TGraph averages;
   TGraph rmss;

   vector<double> ETrue;
   vector<double> gausMean; 
   vector<double> gausSigma;
   vector<double> average;
   vector<double> rms;
   
   
   temp->cd();//This TCanvas is only used since when we do the fit down below 
              //it creates an unwanted TCanvas. We will get rid of it later on 
              //in the function.

   // TCanvas* cccc=new TCanvas("balda","bacla");

   for(unsigned bin = 2; bin < 200; )
   {
      name = "histcorhybrid";
      sprintf(num, "%i", bin);
      name += num;
      //Split up the TH2 into many TH1's for each ETrue bin.
     
      ETrueBin.push_back((TH1F*)inHisto->ProjectionY(name.c_str(),bin, 
                                                     bin + 4*rebin));
      //cout<<" youpali "<<bin<<"   "<<ETrueBin.back()->GetEntries()<<endl;
      if(ETrueBin.back()->GetEntries() > 0)
	{
	  //Fit each ETrue bin to a gaus (iteratively done to get better fit)

	  if(bin > 0) {
	  
	    gaus =new TF1("gaus","gaus(0)",-3,3);
	    gaus->SetParameters(500.,0.,0.2);
	    ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);
	    gaus = ETrueBin.back()->GetFunction("gaus");
            
	    ETrueBin.back()->Fit("gaus", "Q", "",
				 gaus->GetParameter(1) - 2*gaus->
				 GetParameter(2), 1.0);  
	    gaus = ETrueBin.back()->GetFunction("gaus");

            ETrueBin.back()->Fit("gaus", "Q", "",
                                 gaus->GetParameter(1) - 2*gaus->
                                 GetParameter(2), 1.0);
            gaus = ETrueBin.back()->GetFunction("gaus");

	    // if(bin<=16)
	    // gausMean.push_back(ETrueBin.back()->GetMean());
	    // else

	    // cout<<ETrueBin.back()->GetMean()<<"   "<<gaus->GetParameter(1)<<endl;

	    gausMean.push_back(gaus->GetParameter(1));
            gausSigma.push_back(gaus->GetParameter(2)/
                                (1.0 + min(0.0, gaus->GetParameter(1))));

	   }
	   else {
	   
	     gaus =new TF1("gaus","gaus",-3,3);
	     gaus->SetParameters( 500, 10, 5, 0, 0.20 );
	     gaus->FixParameter(2,5);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);

	     gausMean.push_back(gaus->GetParameter(3));
	     gausSigma.push_back(gaus->GetParameter(4)/
				 (1.0 + min(0.0, gaus->GetParameter(3))));
	   }


	  // cout<<bin<<"   "<<median1(ETrueBin.back())<<endl;

	  // TFile oFile( ("tmp/"+name+".root").c_str() ,"RECREATE");
	  // ETrueBin.back()->Write();
	  // gaus->Write();
	  // oFile.Close();

	  // cccc->cd();
	  // ETrueBin.back()->Draw();
	  // // //   gaus->Draw("same");
	  // cccc->SaveAs( ("tmp/"+name+".png").c_str() );
	  // cccc->SaveAs( ("tmp/"+name+".C").c_str() );

            ETrue.push_back(bin + 2.0*rebin);

	    //  if(bin<=16)
	    //  gausMean.push_back(ETrueBin.back()->GetMean());
	    // else
	    //   gausMean.push_back(gaus->GetParameter(1));

	    
	   
            average.push_back(ETrueBin.back()->GetMean());
            rms.push_back(ETrueBin.back()->GetRMS());

	    //   cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gaus->GetParameter(1)<<"   "<<ETrueBin.back()->GetMeanError()<<"   "<<gaus->GetParError(1)<<endl;
         }
      
      bin += 2*rebin;
      
      //Increase bin size with increasing ETrue since there are fewer high 
      //energy events than low energy ones.
      if(bin > 10) rebin = 2.0;
      if(bin > 100) rebin = 5.0; //20
      if(bin > 1000) rebin = 20.0; //50
   }

   // delete cccc;

   response = TGraph(ETrue.size(), &ETrue[0], &gausMean[0]); //Fill the graphs
   //response = TGraph(ETrue.size(), &ETrue[0], &average[0]); //Fill the graphs
   resolution = TGraph(ETrue.size(),&ETrue[0], &gausSigma[0]);
   averages =  TGraph(ETrue.size(), &ETrue[0], &average[0]);
   rmss = TGraph(ETrue.size(), &ETrue[0], &rms[0]);

   //Set up the graphs to look how you want them to.
   response.SetMarkerStyle(22);
   response.SetMarkerSize(0.8);
   response.SetMarkerColor(4);

   resolution.SetMarkerStyle(22);
   resolution.SetMarkerSize(0.8);
   resolution.SetMarkerColor(4);

   averages.SetMarkerStyle(22);
   averages.SetMarkerSize(0.8);
   averages.SetMarkerColor(4);

   rmss.SetMarkerStyle(22);
   rmss.SetMarkerSize(0.8);
   rmss.SetMarkerColor(4);

   line->SetLineStyle(1);
   line->SetLineWidth(2);
   line->SetLineColor(2);


   //  gStyle->SetOptStat(0); 
   //gStyle->SetOptFit(0);
   canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 
			1000, 500);


 
   canvas->Divide(2, 1);
   temp->~TCanvas();  //destroy the TCanvas 

   canvas->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle("Response");
   respHisto->Draw();
   response.Draw("P");
   line->Draw();

   canvas->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   resoHisto->SetStats(0);
   resoHisto->SetTitle("Resolution");
   resoHisto->Draw();
   resolution.Draw("P");


   respHisto->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");
   respHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

   resoHisto->GetYaxis()->SetTitle("#sigma(E)/E_{true}");
   resoHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

   //MM Fit Resolution
   TF1* f=new TF1( ("ResoFit"+ (string)(inHisto->GetName())).c_str(),"sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",20,1000);// 3.5*4
   f->SetParameters(0.06,1.20,0.);
   f->SetParLimits(0,0,10);
   f->SetParLimits(1,0,10);
   f->SetParLimits(2,0,10);
   resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
   resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
   resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"R");

   

   string legend;
   int fres0 = (int)(f->GetParameter(0)*100.);
   int fres1 = (int)(10.*(f->GetParameter(0)*100.-fres0));
   int fres2 = (int)(f->GetParameter(1)*100.);
   // char text[100];
   // sprintf(text,"#sigma/E = %i%/#sqrt{E} + %i.%i%",fres2,fres0,fres1);
   TString text = "#sigma/E = ";
   text+=(int)fres2;
   text+="%/#sqrt{E} + ";
   text+=(int)fres0;
   text+=".";
   text+=(int)fres1;

  legend += text;
  TLegend *leg=new TLegend(0.30,0.75,0.85,0.85);
  leg->AddEntry((&resolution),legend.c_str(),"lp");
  leg->SetTextSize(0.04);
  leg->Draw();



}

void drawEtaDependence(TH2F* inHisto, TGraph& responseEta)
{
   if(inHisto->GetEntries() == 0) return;

   vector<TH1F*> etaBin;
   TF1* gaus; 
   TString name;
   char num[4];

   TCanvas* canvas;
   TCanvas* temp = new TCanvas();
   TLine* line = new TLine(0, 0, 3, 0);

   TH2F* respHisto = new TH2F("respHisto", "", 30, 0.0, 3.00, 100, -1.0, 1.0);

   TGraph averages;
   TGraph rmss;

   vector<double> etaAverage;
   vector<double> gausMean; 
   vector<double> gausSigma;
   vector<double> average;
   vector<double> etaRms;
   
   
   temp->cd();//This TCanvas is only used since when we do the fit down below 
              //it creates an unwanted TCanvas. We will get rid of it later on 
              //in the function.

   float mR=0;
   float mR2=0;
   int N=0;

   //  TCanvas* cccc=new TCanvas("bala","bala");

   for(unsigned bin = 1; bin < (unsigned)inHisto->GetNbinsX(); bin = bin + 1)
   {
      name = "histEta";
      //      sprintf(num, "%i", bin);
      TString name2 = name;
      name2 += bin;
      //Split up the TH2 into many TH1's for each eta bin.
   
      etaBin.push_back((TH1F*)inHisto->ProjectionY(name,bin, bin + 1));
      
      name += inHisto->GetXaxis()->GetBinCenter(bin);
      

      if(etaBin.back()->GetEntries() > 0)
      {
         //Fit each eta bin to a gaus (iteratively done to get better fit)

	gaus =new TF1("gaus","gaus(0)",-3,3);
	gaus->SetParameters(500.,0.,0.2);

         etaBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);
	 // gaus = etaBin.back()->GetFunction("gaus");
         
         etaBin.back()->Fit("gaus", "Q", "",
                            gaus->GetParameter(1) - 2*gaus->
                            GetParameter(2), 1.0);  
	 // gaus = etaBin.back()->GetFunction("gaus");
         
         etaBin.back()->Fit("gaus", "Q", "",
                            gaus->GetParameter(1) - 2*gaus->
                            GetParameter(2), 1.0);
	 // gaus = etaBin.back()->GetFunction("gaus");
         
         etaAverage.push_back(inHisto->GetXaxis()->GetBinCenter(bin));
         etaRms.push_back(0.1);
	 
	 gausMean.push_back( gaus->GetParameter(1) );
         //gausMean.push_back(etaBin.back()->GetMean());
         gausSigma.push_back(etaBin.back()->GetRMS());

	 if(etaAverage.back()>1.6) {
	   mR += gausMean.back();
	   mR2 += gausMean.back()*gausMean.back();
	   N++;
	 }

	 // cccc->cd();
	 // etaBin.back()->Draw();
	 // gaus->Draw("same");
	 // cccc->SaveAs( ("tmp/Eta"+name+".png") );
	 // cccc->SaveAs( ("tmp/Eta"+name+".root") );
	 cout<<gaus->GetParameter(1)<<endl;

      }
            
      //Increase bin size with increasing eta since there are fewer high 
      //energy events than low energy ones.
   }

   //  delete cccc;

   responseEta = TGraph(etaAverage.size(), &etaAverage[0], &gausMean[0]); 
//&etaRms[0], &gausSigma[0]); 

   responseEta.SetMarkerStyle(22);
   responseEta.SetMarkerSize(1);
   responseEta.SetMarkerColor(4);   

   line->SetLineStyle(1);
   line->SetLineWidth(2);
   line->SetLineColor(2);

 canvas = new TCanvas( ("canvas"+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName()) ).c_str(), 
                        1000, 500);

   
   temp->~TCanvas();  //destroy the TCanvas 
   
   canvas->cd();
   canvas->SetFillColor(0);

   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle("Response");
   respHisto->Draw();
   responseEta.Draw("P");
   line->Draw();   
   respHisto->GetXaxis()->SetRangeUser(0,3);
   respHisto->GetYaxis()->SetRangeUser(-0.2,0.2);
   respHisto->GetXaxis()->SetTitle("|#eta|");
   respHisto->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");


   //Spread
   cout<<" Endcap spread and mean for "<<inHisto->GetName()<<endl;

   mR/=N;
   mR2/=N;
   cout<<" mean = "<<mR<<endl;
   cout<<" rms = "<<sqrt(mR2- pow(mR,2) )<<endl;
   
}

void drawCompare(TGraph& response1, TGraph& response2, TGraph& resolution1, TGraph& resolution2)
{
   

   TCanvas* Compare = new TCanvas("Compare" ,"", 1000, 500);
   TH2F * respHisto = new TH2F("respHisto", "", 100, 0, 1000, 100, -1, 1);
   TH2F * resoHisto = new TH2F("resoHisto", "", 100, 0, 1000, 100, 0, 1);
   TLegend * legend1 = new TLegend(0.75, 0.75, 0.95, 0.9);
   TLegend * legend2 = new TLegend(0.75, 0.75, 0.95, 0.9);

   response1.SetMarkerColor(4);
   response1.SetMarkerStyle(22);
   response1.SetMarkerSize(0.8);

   resolution1.SetMarkerColor(4);
   resolution1.SetMarkerStyle(22);
   resolution1.SetMarkerSize(0.8);

   response2.SetMarkerColor(2);
   response2.SetMarkerStyle(22);
   response2.SetMarkerSize(0.8);

   resolution2.SetMarkerColor(2);
   resolution2.SetMarkerStyle(22);
   resolution2.SetMarkerSize(0.8);
   
   legend1->AddEntry(&response1, "Raw");
   legend1->AddEntry(&response2, "Corrected");
   legend2->AddEntry(&resolution1, "Raw");
   legend2->AddEntry(&resolution2, "Corrected");

   Compare->Divide(2,1);
  
   Compare->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle("Response");
   respHisto->Draw();
   response1.Draw("P");
   response2.Draw("P");
   legend1->Draw();

   Compare->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   resoHisto->SetStats(0);
   resoHisto->SetTitle("Resolution");
   resoHisto->Draw();
   resolution1.Draw("P");
   resolution2.Draw("P");
   legend2->Draw();

}

//Takes apart a TTree from a root file and puts the wanted information into 
//vectors. 
void getValuesFromTree(TTree* tree, vector<double>& ETrueEnergies, 
                       vector<double>& ecalEnergies, 
                       vector<double>& hcalEnergies, vector<double>& etas, 
                       vector<double>& phis)
{
   Float_t         true_;
   Float_t         p_;
   Float_t         ecal_;
   Float_t         hcal_;
   Float_t         eta_;
   Float_t         phi_;
   TBranch        *b_true; 
   TBranch        *b_p;   
   TBranch        *b_ecal;   
   TBranch        *b_hcal;   
   TBranch        *b_eta;    
   TBranch        *b_phi;    

   tree->SetMakeClass(1);
   
   if(tree->GetBranchStatus("true"))
      tree->SetBranchAddress("true", &true_, &b_true);
   tree->SetBranchAddress("p", &p_, &b_p);
   tree->SetBranchAddress("ecal", &ecal_, &b_ecal);
   tree->SetBranchAddress("hcal", &hcal_, &b_hcal);
   tree->SetBranchAddress("eta", &eta_, &b_eta);
   tree->SetBranchAddress("phi", &phi_, &b_phi);

   double sigmaEcalHcal=1;

   for( unsigned entry = 0; entry < std::min((int)5000000,(int)(tree->GetEntriesFast()) ); entry++)
     {
       tree->GetEntry(entry);


       // if(ecal_<0.4) continue; //FIXME MM

       if(tree->GetBranchStatus("true"))
	 ETrueEnergies.push_back(true_);
       else
	 ETrueEnergies.push_back(p_);
       ecalEnergies.push_back(ecal_ );
       hcalEnergies.push_back(hcal_);
       etas.push_back(eta_);
       phis.push_back(phi_);


       if(fabs(eta_)<1.5) 
	 sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal_ + hcal_), 1.0)));
       else
	 sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal_ + hcal_), 1.0)));

       sigmas.push_back(sigmaEcalHcal);

     }

   cout<<" Entries "<<ecalEnergies.size()<<endl;

}



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
Calibration* barrelWithEcalHcalCalib = new Calibration(200, true);
Calibration* barrelWithEcalCalib = new Calibration(200, true);
Calibration* barrelWithHcalCalib = new Calibration(200, true);
Calibration* endcapWithEcalHcalCalib = new Calibration(200, false);
Calibration* endcapWithEcalCalib = new Calibration(200, false);
Calibration* endcapWithHcalCalib = new Calibration(200, false);

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
TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 150, 
                                   -1.5, 1.5);
TH2F* corrEndcapEcalHcal = new TH2F("corrEndcapEcalHcal", "", 1000, 0, 1000, 
                                    150, -1.5, 1.5);
TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 
                                       1000, 150, -1.5, 1.5);
TH2F* rawEndcapHcal = new TH2F("rawEndcapHcal","", 1000, 0, 1000, 150, -1.5, 
                               1.5);
TH2F* corrEndcapHcal = new TH2F("corrEndcapHcal", "", 1000, 0, 1000, 150, -1.5,
                                1.5);
TH2F* corrEtaEndcapHcal = new TH2F("corrEtaEndcapHcal", "", 1000, 0, 1000, 150,
                                   -1.5, 1.5);

TH2F * rawEtaDependenceEH = new TH2F("rawEtaDependenceEH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * corrEtaDependenceEH = new TH2F("corrEtaDependenceEH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * hcorrEtaDependenceEH = new TH2F("hcorrEtaDependenceEH","Response vs. Eta", 75, 0., 3.0, 150, -1.0,1.0 );

TH2F * rawEtaDependenceH = new TH2F("rawEtaDependenceH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * corrEtaDependenceH = new TH2F("corrEtaDependenceH","Response vs. Eta", 75, 0.0, 3.0, 150, -1.0,1.0 );
TH2F * hcorrEtaDependenceH = new TH2F("hcorrEtaDependenceH","Response vs. Eta", 75, 0., 3.0, 150, -1.0,1.0 );

TH1F * trueTempHisto = new TH1F("trueTempHisto", "true", 200, 0, 200);
TH1F * ecalTempHisto = new TH1F("ecalTempHisto", "ecal", 200, 0, 200);

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

///////////////////////////////////////////////////////////////////////////////
//This is the main of the macro. Everything that you want output must be added 
//in here. I have it so that all the variables that I use were defined above 
//since it looks neater.
///////////////////////////////////////////////////////////////////////////////
void calibChris() 
{
   gROOT->Reset();
   gStyle->SetCanvasColor(0);

   InitBarrelAlpha();
   LoadOldThresholds();

   gStyle->SetOptFit(0);

   //Open the file, get the tree and fill of the vectors of values you need.
   //inputFile = TFile::Open("IsolatedChargedHadronsFromQCD.root");
   //inputFile = TFile::Open("pfcalibTestTag_all.root");
   // inputFile = TFile::Open("IsolatedChargedHadronsFromMinBias.root");
   // sTree = (TTree*)inputFile->Get("s;1");
   
   TChain* chain= new TChain("s");
   // chain->Add("../PFCalibration/ParticleGun.root");
   //chain->Add("../data/73XHcalCor.root");//PGun73XNoPu.root");
   //chain->Add("PGun73X_oct_18.root");
   //chain->Add("../PGun73X_oct_18.root");
   //chain->Add("../oct_21_PGun73X/PGun73X_97.root");
   //chain->Add("../final_v1.root");
   //chain->Add("/home/work/spandey/public/PF_Cal/samples/oct_31/PGun73X_oct_31.root");
   chain->Add("/home/work/spandey/public/PF_Cal/samples/nov_22/PGun_autoMC_200GeV_500GeV_nov_22.root");
   sTree = (TTree*)chain;

   getValuesFromTree(sTree, ETrueEnergies, ecalEnergies, 
                     hcalEnergies, etas, phis);

   
   //Create all the ABC objects you need with increasing bin size
   //since there are fewer events at higher energies. 
   cout<<"Creating abc and alphabeta objects..."<<endl;
  
   BinsETrue.clear();
   BinsETrueEta.clear();

   for(double bin = 0.0; bin < 10.0; bin = bin + lBs)
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + lBs, true));
       barrelABCEcal.push_back(new ABC(bin, bin + lBs, true));
       barrelABCHcal.push_back(new ABC(bin, bin + lBs, true));
       endcapABCEcalHcal.push_back(new ABC(bin, bin + lBs, false));
       endcapABCEcal.push_back(new ABC(bin, bin + lBs,false));
       endcapABCHcal.push_back(new ABC(bin, bin + lBs, false));
       BinsETrue.push_back(bin);
     }
   
   
   
   for(double bin = 10.0; bin < 100.0 ; bin = bin + mBs) //2
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + mBs, true));
       barrelABCEcal.push_back(new ABC(bin, bin + mBs, true));
      barrelABCHcal.push_back(new ABC(bin, bin + mBs, true));
      endcapABCEcalHcal.push_back(new ABC(bin, bin + mBs, false));
      endcapABCEcal.push_back(new ABC(bin, bin + mBs,false));
      endcapABCHcal.push_back(new ABC(bin, bin + mBs, false));
      BinsETrue.push_back(bin);
   }
   
   
   for(double bin = 100.0; bin < 1000.0 ; bin = bin + hBs) //10
   {
     barrelABCEcalHcal.push_back(new ABC(bin, bin + hBs, true));
     barrelABCEcal.push_back(new ABC(bin, bin + hBs, true));
     barrelABCHcal.push_back(new ABC(bin, bin + hBs, true));
     endcapABCEcalHcal.push_back(new ABC(bin, bin + hBs, false));
     endcapABCEcal.push_back(new ABC(bin, bin + hBs,false));
     endcapABCHcal.push_back(new ABC(bin, bin + hBs, false));  
     BinsETrue.push_back(bin);
   }
   BinsETrue.push_back( BinsETrue.back() + hBs );



   cout<<"barrelABCEcalHcal size: "<<barrelABCEcalHcal.size()<<endl;
   cout<<"barrelABCEcal size: "<<barrelABCEcal.size()<<endl;
   cout<<"barrelABCHcal size: "<<barrelABCHcal.size()<<endl;
   /*for(int i = 0; i < barrelABCEcalHcal.size(); i++) {
     if((barrelABCEcalHcal.at(i)->getSize()) != 0 )
       cout<<"EcalHcal: found one!! at "<<i<<endl;
     if((barrelABCEcal.at(i)->getSize()) != 0 )
       cout<<"Ecal: found one!! at "<<i<<endl;
     if((barrelABCHcal.at(i)->getSize()) != 0 )
       cout<<"Hcal: found one!! at "<<i<<endl;

       }*/
   //cout<<"(barrelABCEcalHcal.at(300)->getBinLowEdge()) : "<<(barrelABCEcalHcal.at(300)->getBinLowEdge())<<endl;
   //cout<<"(barrelABCEcalHcal.at(300)->getBinHighEdge()) : "<<(barrelABCEcalHcal.at(300)->getBinHighEdge())<<endl;
   //cout<<"(barrelABCEcalHcal.at(0)->getETrue(0)) : "<<(barrelABCEcalHcal.at(0)->getETrue(0))<<endl;




   
   for(double bin = 0.0; bin < 10.0; bin = bin + lBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, false));
       BinsETrueEta.push_back(bin);
       // barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs, true));
       // endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs, false));

     }
   for(double bin = 10.0; bin < 100.0 ; bin = bin + mBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       BinsETrueEta.push_back(bin);
     //   barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs, true));
     //   endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs, false));
     }
   
   for(double bin = 100.0; bin < 1000.0 ; bin = bin + hBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       BinsETrueEta.push_back(bin);
       // barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs, true));
       // endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs, false));
     }
   BinsETrueEta.push_back( BinsETrueEta.back() + hBs*RBE );
   
   //Fill all the ABC Objects with their respective events. They are all 
   //divided up into the six possible case ( (endcap or barrel)x(ecalhcal or 
   //ecal or hcal))
   

   TH1F* EcalSpectrum=new TH1F("EcalSpectrum","EcalSpectrum",1000,0,100);

   cout<<"Filling abc objects..."<<endl;
   for( unsigned bin = 0; bin < barrelABCEcal.size(); ++bin)
     {
       barrelABCEcalHcal[bin]->computeA(aEH);
       barrelABCEcal[bin]->computeA(aE);
       barrelABCHcal[bin]->computeA(aH);

       endcapABCEcalHcal[bin]->computeA(aEHe);
       endcapABCEcal[bin]->computeA(aEe);
       endcapABCHcal[bin]->computeA(aHe);
     }
   

   cout<<"ETrueEnergies size: "<<ETrueEnergies.size()<<endl;
   //cout<<"GetETrueBinEta(100): "<<GetETrueBinEta(99.9858)<<endl;
   cout<<"GetETrueBinEta(100): "<<GetETrueBinEta(100)<<endl;

   //Filling ==============================
   {

     unsigned bin = 0;
     for(unsigned entry = 0; entry < ETrueEnergies.size(); entry++)
       {
	 etrue = ETrueEnergies[entry];
	 ecal = ecalEnergies[entry];
	 hcal = hcalEnergies[entry];
	 eta = etas[entry];

	 if(hcal == 0.0) continue;
	 if( etrue <1 ) continue;
	 if( etrue >200 ) continue;
	 // if( etrue <10 ) continue;
	 // if( etrue >12 ) continue;

	 bin = GetETrueBin( etrue );

	 if( ecal > 0.0 && hcal > 0.0)
	   {
	     barrelABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
            
	     if(eta<1.3)
	       EcalSpectrum->Fill(ecal);

	   }
	 else if(ecal > 0.0)
	   {
	     barrelABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	     endcapABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	   }
	 else if(hcal > 0.0)
	   {
	     barrelABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	   }
         
	 bin = GetETrueBinEta( etrue );

	 if(bin < barrelAlphaBetaEcalHcal.size())
	   {

	    
	     if( ecal > 0.0 && hcal >= 0.0 )
	       {
		 endcapAlphaBetaEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
		 barrelAlphaBetaEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	       }
	     else {
	       endcapAlphaBetaHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	       barrelAlphaBetaHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     }
	   }
       }
   }


   //   cout<<"#####################################################################"<<endl;
   //for(unsigned bin = 2; bin < barrelABCEcalHcal.size() - 1; ++bin) {
   // cout<<"barrelABCEcalHcal["<<bin<<"]->isEmptyInFitRange(): "<<barrelABCEcalHcal[bin]->isEmptyInFitRange()<<endl;
   //}
   //cout<<"#####################################################################"<<endl;

   //Filling ==============================  
   
   TFile* file=new TFile("output.root","recreate");
       EcalSpectrum->Write();
   file->Close();

   
   //Compute the calibration constants along with their uncertainties for each
   //ETrue bin, then add their values to a Calibration object.
   cout<<"Computing a, b, c coefficients..."<<endl;

   for(unsigned bin = 2; bin < barrelABCEcalHcal.size() - 1; ++bin)
   {
      
      if(!barrelABCEcalHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcalHcal[bin]->computeETrueAverage();
         barrelABCEcalHcal[bin]->computeETrueRMS();
         barrelABCEcalHcal[bin]->computeA(aEH);
         barrelABCEcalHcal[bin]->computeBC();
	 //exit(0);
      }
      
      if(!barrelABCEcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcal[bin]->computeETrueAverage();
         barrelABCEcal[bin]->computeETrueRMS();
         barrelABCEcal[bin]->computeA(aEH);
         barrelABCEcal[bin]->computeB();
      }
      if(!barrelABCHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCHcal[bin]->computeETrueAverage();
         barrelABCHcal[bin]->computeETrueRMS();
         barrelABCHcal[bin]->computeA(aH);
         barrelABCHcal[bin]->computeC();
      }
      if(!endcapABCEcalHcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcalHcal[bin]->computeETrueAverage();
         endcapABCEcalHcal[bin]->computeETrueRMS();
         endcapABCEcalHcal[bin]->computeA(aEHe);
         endcapABCEcalHcal[bin]->computeBC();
      }
      if(!endcapABCEcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcal[bin]->computeETrueAverage();
         endcapABCEcal[bin]->computeETrueRMS();
         endcapABCEcal[bin]->computeA(aEe);
         endcapABCEcal[bin]->computeB();
      }
      if(!endcapABCHcal[bin]->isEmptyInFitRange())
      {
         endcapABCHcal[bin]->computeETrueAverage();
         endcapABCHcal[bin]->computeETrueRMS();
         endcapABCHcal[bin]->computeA(aHe);
         endcapABCHcal[bin]->computeC();
      }
      

      if(!barrelABCEcalHcal[bin]->isEmpty() && 
         barrelABCEcalHcal[bin]->getBinHighEdge() >
         barrelWithEcalHcalCalib->getETrueMax())
      {
         barrelWithEcalHcalCalib->setETrueMax(
            barrelABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!barrelABCEcal[bin]->isEmpty() && 
         barrelABCEcal[bin]->getBinHighEdge() >
         barrelWithEcalCalib->getETrueMax())
      {
         barrelWithEcalCalib->setETrueMax(
            barrelABCEcal[bin]->getBinHighEdge());
      }
      if(!barrelABCHcal[bin]->isEmpty() && 
         barrelABCHcal[bin]->getBinHighEdge() >
         barrelWithHcalCalib->getETrueMax())
      {
         barrelWithHcalCalib->setETrueMax(
            barrelABCHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcalHcal[bin]->isEmpty() && 
         endcapABCEcalHcal[bin]->getBinHighEdge() >
         endcapWithEcalHcalCalib->getETrueMax())
      {
         endcapWithEcalHcalCalib->setETrueMax(
            endcapABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcal[bin]->isEmpty() && 
         endcapABCEcal[bin]->getBinHighEdge() >
         endcapWithEcalCalib->getETrueMax())
      {
         endcapWithEcalCalib->setETrueMax(
            endcapABCEcal[bin]->getBinHighEdge());
      }
      if(!endcapABCHcal[bin]->isEmpty() && 
         endcapABCHcal[bin]->getBinHighEdge() >
         endcapWithHcalCalib->getETrueMax())
      {
         endcapWithHcalCalib->setETrueMax(
            endcapABCHcal[bin]->getBinHighEdge());
      }
 

      barrelWithEcalHcalCalib->addGraphPoints(barrelABCEcalHcal[bin]); 
      barrelWithEcalCalib->addGraphPoints(barrelABCEcal[bin]); 
      barrelWithHcalCalib->addGraphPoints(barrelABCHcal[bin]); 
      endcapWithEcalHcalCalib->addGraphPoints(endcapABCEcalHcal[bin]); 
      endcapWithEcalCalib->addGraphPoints(endcapABCEcal[bin]); 
      endcapWithHcalCalib->addGraphPoints(endcapABCHcal[bin]);                 


   }
   
   cout<<"Fitting a, b, c coefficients..."<<endl;
   //Initialize all the ABC graphs in the calibration objects.
   barrelWithEcalHcalCalib->initializeGraphs("abc");

   barrelWithEcalCalib->initializeGraphs("abc");

   barrelWithHcalCalib->initializeGraphs("abc");   

   endcapWithEcalHcalCalib->initializeGraphs("abc");

   endcapWithEcalCalib->initializeGraphs("abc");

   endcapWithHcalCalib->initializeGraphs("abc");   


   //Define the functions that you will fit your ABC calibration constants to.
   functionBarrelEcalHcalA = new TF1("functionBarrelEcalHcalA","[0]", 0, 1000);
   // functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
   //functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","[0]+(([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3])))",0,1000); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",0,1000);
  
   functionEndcapEcalHcalA = new TF1("functionEndcapEcalHcalA","[0]", 0, 1000);
   //functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
   functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);

   functionBarrelHcalA = new TF1("functionBarrelHcalA","[0]", 0, 1000);
   functionBarrelHcalB = new TF1("functionBarrelHcalB","[0]", 0, 1000);
   // functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
  
   functionEndcapHcalA = new TF1("functionEndcapHcalA","[0]", 0, 1000);
   functionEndcapHcalB = new TF1("functionEndcapHcalB","[0]", 0, 1000);
   functionEndcapHcalC = new TF1("functionEndcapHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])


   //Set the parameters of the function you just defined.
   functionBarrelEcalHcalA->FixParameter(0, aEH);
   // functionBarrelEcalHcalB->SetParameters(1.32991e+00,1.52538e+02,-7.89866e+02,2.56441e-01,
   // 					  -1.13519e+00,2.71646e+00,1.86857e-01,4.68538e-01);
   functionBarrelEcalHcalB->SetParameters(-1.38681e+01,1.49678e+01,3.45153e+00,1.04212e+00,-2.00910e-02,
					  9.41444e-16,-1.31641e+00,-7.07963e+00);
   
   functionBarrelEcalHcalC->SetParameters(1.70114,0.404676,-3.88962,1.2109e+06,
					  0.970741,0.0527482,2.60552,-0.8956);
   functionBarrelEcalHcalC->FixParameter(0, 1.70114 );
   functionBarrelEcalHcalC->FixParameter(1, 0.404676 );
   functionBarrelEcalHcalC->FixParameter(2, -3.88962);
   functionBarrelEcalHcalC->FixParameter(3, 1.2109e+06);
   functionBarrelEcalHcalC->FixParameter(4, 0.970741 );
   functionBarrelEcalHcalC->FixParameter(5, 0.0527482);
   functionBarrelEcalHcalC->FixParameter(6, 2.60552);
   functionBarrelEcalHcalC->FixParameter(7, -0.8956);	


   functionEndcapEcalHcalA->FixParameter(0, aEHe);
  
   functionEndcapEcalHcalB->SetParameters(0.930193,11.9536,-30.0337,0.76133,
					  0.0776373,7.3809e-10,0.158734,-6.92163);
   functionEndcapEcalHcalB->FixParameter(0, 0.930193);
   functionEndcapEcalHcalB->FixParameter(1, 11.9536);
   functionEndcapEcalHcalB->FixParameter(2, -30.0337);
   functionEndcapEcalHcalB->FixParameter(3, 0.76133);
   functionEndcapEcalHcalB->FixParameter(4, 0.0776373);
   functionEndcapEcalHcalB->FixParameter(5, 7.3809e-10);
   functionEndcapEcalHcalB->FixParameter(6, 0.158734);
   functionEndcapEcalHcalB->FixParameter(7, -6.92163);

   functionEndcapEcalHcalC->SetParameters(-0.436687,2.73698,-3.1509,1.20536,
					  -1.39685,0.0180331,0.270058,-2.30372);
   functionEndcapEcalHcalC->FixParameter(0, -0.436687 );
   functionEndcapEcalHcalC->FixParameter(1, 2.73698);
   functionEndcapEcalHcalC->FixParameter(2, -3.1509);
   functionEndcapEcalHcalC->FixParameter(3, 1.20536);
   functionEndcapEcalHcalC->FixParameter(4, -1.39685);
   functionEndcapEcalHcalC->FixParameter(5, 0.0180331);
   functionEndcapEcalHcalC->FixParameter(6, 0.270058);
   functionEndcapEcalHcalC->FixParameter(7, -2.30372);			  
  

   functionBarrelHcalA->FixParameter(0, aH);
   functionBarrelHcalB->FixParameter(0, 0.0);
   functionBarrelHcalC->SetParameters(1.58827e+00,4.06865e-01,-3.69939e+00,1.28926e+03,
					7.13400e-01,2.21925e-02,1.47842e+00,-1.22041e+00);
   
   functionEndcapHcalA->FixParameter(0, aHe);
   functionEndcapHcalB->FixParameter(0, 0.0);
   
   functionEndcapHcalC->SetParameters(1.13795,1.21698,-3.81192,115.409,
				      0.673456,0.217077,1.95596,-0.252215);
   functionEndcapHcalC->FixParameter(0, 1.13795);
   functionEndcapHcalC->FixParameter(1, 1.21698);
   functionEndcapHcalC->FixParameter(2, -3.81192);
   functionEndcapHcalC->FixParameter(3, 115.409);
   functionEndcapHcalC->FixParameter(4, 0.673456);
   functionEndcapHcalC->FixParameter(5, 0.217077);
   functionEndcapHcalC->FixParameter(6, 1.95596);
   functionEndcapHcalC->FixParameter(7, -0.252215);
   
   //Fit the functions that were defined above to the calibration graphs that
   //are stored in the calibration objects.

   barrelWithEcalHcalCalib->fitAsToFunction(functionBarrelEcalHcalA);
   barrelWithEcalHcalCalib->fitBsToFunction(functionBarrelEcalHcalB);
   barrelWithEcalHcalCalib->fitBsToFunction();
   barrelWithEcalHcalCalib->fitCsToFunction(functionBarrelEcalHcalC);
   barrelWithEcalHcalCalib->fitCsToFunction();
   barrelWithEcalHcalCalib->fitCsToFunction();
   endcapWithEcalHcalCalib->fitAsToFunction(functionEndcapEcalHcalA);
   endcapWithEcalHcalCalib->fitBsToFunction(functionEndcapEcalHcalB);
   endcapWithEcalHcalCalib->fitBsToFunction();
   endcapWithEcalHcalCalib->fitCsToFunction(functionEndcapEcalHcalC);
   endcapWithEcalHcalCalib->fitCsToFunction();
   endcapWithEcalHcalCalib->fitCsToFunction();
   cout<<"Fit check 1\n";
   barrelWithHcalCalib->fitAsToFunction(functionBarrelHcalA);
   barrelWithHcalCalib->fitBsToFunction(functionBarrelHcalB);
   cout<<"Fit check 1.a\n";
   barrelWithHcalCalib->fitBsToFunction();
   cout<<"Fit check 1.b\n";
   barrelWithHcalCalib->fitCsToFunction(functionBarrelHcalC);
   barrelWithHcalCalib->fitCsToFunction();
   
   barrelWithHcalCalib->fitCsToFunction();
   cout<<"Fit check 2\n";
   endcapWithHcalCalib->fitAsToFunction(functionEndcapHcalA);
   endcapWithHcalCalib->fitBsToFunction(functionEndcapHcalB);
   endcapWithHcalCalib->fitBsToFunction();
   endcapWithHcalCalib->fitCsToFunction(functionEndcapHcalC);
   endcapWithHcalCalib->fitCsToFunction();
   endcapWithHcalCalib->fitCsToFunction();
   cout<<"Fit check 3\n";
   //exit(0);
   //Here we fill up the AlphaBeta objects, compute alpha and beta, then add 
   //them to the Calibration objects. 
   cout<<"Computing alpha and beta coefficients..."<<endl;
   for(unsigned bin = 2; bin < barrelAlphaBetaEcalHcal.size() - 1; bin++)
   {
     for(unsigned entry = 0; entry < barrelAlphaBetaEcalHcal[bin]->getSize(); entry++)
      {
         
         etrue = barrelAlphaBetaEcalHcal[bin]->getETrue(entry);
         ecal = barrelAlphaBetaEcalHcal[bin]->getEcal(entry);
         hcal = barrelAlphaBetaEcalHcal[bin]->getHcal(entry);
         bpar = 1.0;
         cpar = 1.0;
         

         if(ecal > 0 && hcal > 0)
         {
            bpar = barrelWithEcalHcalCalib->getFunctionB()->Eval(etrue);
            cpar = barrelWithEcalHcalCalib->getFunctionC()->Eval(etrue);
         }
         else if(ecal > 0)
            bpar = barrelWithEcalHcalCalib->getFunctionB()->Eval(etrue);
       
         
         barrelAlphaBetaEcalHcal[bin]->correctEcal(entry, bpar);
         barrelAlphaBetaEcalHcal[bin]->correctHcal(entry, cpar);
      }

   
     for(unsigned entry = 0; entry < barrelAlphaBetaHcal[bin]->getSize(); entry++)
       {
         
	 etrue = barrelAlphaBetaHcal[bin]->getETrue(entry);
	 ecal = barrelAlphaBetaHcal[bin]->getEcal(entry);
	 hcal = barrelAlphaBetaHcal[bin]->getHcal(entry);
	 bpar = 1.0;
	 cpar = 1.0;
	 
	 if(hcal > 0 && ecal==0)
	   cpar = barrelWithHcalCalib->getFunctionC()->Eval(etrue);
      
	 // if(etrue<10)
	 //   cout<<etrue<<"   "<<cpar<<endl;

	 barrelAlphaBetaHcal[bin]->correctEcal(entry, bpar);
         barrelAlphaBetaHcal[bin]->correctHcal(entry, cpar);

       }     

     for(unsigned entry = 0; entry < endcapAlphaBetaEcalHcal[bin]->getSize(); entry++)
	{

	  etrue = endcapAlphaBetaEcalHcal[bin]->getETrue(entry);
	  ecal = endcapAlphaBetaEcalHcal[bin]->getEcal(entry);
	  hcal = endcapAlphaBetaEcalHcal[bin]->getHcal(entry);
	  bpar = 1.0;
	  cpar = 1.0;
         
	  if(ecal > 0 && hcal > 0)
	    {
	      bpar = endcapWithEcalHcalCalib->getFunctionB()->Eval(etrue);
	      cpar = endcapWithEcalHcalCalib->getFunctionC()->Eval(etrue);
	    }
	  else if(ecal > 0)
            bpar = endcapWithEcalHcalCalib->getFunctionB()->Eval(etrue);
	  
	  endcapAlphaBetaEcalHcal[bin]->correctEcal(entry, bpar);
	  endcapAlphaBetaEcalHcal[bin]->correctHcal(entry, cpar);
	}

      for(unsigned entry = 0; entry < endcapAlphaBetaHcal[bin]->getSize(); entry++)
	{

	  etrue = endcapAlphaBetaHcal[bin]->getETrue(entry);
	  ecal = endcapAlphaBetaHcal[bin]->getEcal(entry);
	  hcal = endcapAlphaBetaHcal[bin]->getHcal(entry);
	  bpar = 1.0;
	  cpar = 1.0;
         
	  if(ecal == 0 && hcal > 0)
	    {
	      cpar = endcapWithHcalCalib->getFunctionC()->Eval(etrue);
	    }
	  endcapAlphaBetaHcal[bin]->correctEcal(entry, bpar);
	  endcapAlphaBetaHcal[bin]->correctHcal(entry, cpar);
	}
      
      
      barrelAlphaBetaEcalHcal[bin]->computeSigmaEcalHcal();
      barrelAlphaBetaEcalHcal[bin]->computeETrueAverage();
      barrelAlphaBetaEcalHcal[bin]->computeETrueRMS();

      barrelAlphaBetaHcal[bin]->computeSigmaEcalHcal();
      barrelAlphaBetaHcal[bin]->computeETrueAverage();
      barrelAlphaBetaHcal[bin]->computeETrueRMS();
      
      endcapAlphaBetaEcalHcal[bin]->computeSigmaEcalHcal();
      endcapAlphaBetaEcalHcal[bin]->computeETrueAverage();
      endcapAlphaBetaEcalHcal[bin]->computeETrueRMS();

      endcapAlphaBetaHcal[bin]->computeSigmaEcalHcal();
      endcapAlphaBetaHcal[bin]->computeETrueAverage();
      endcapAlphaBetaHcal[bin]->computeETrueRMS();

      if(barrelAlphaBetaEcalHcal[bin]->computeAlphaBeta())
      {
	barrelWithEcalHcalCalib->addGraphPoints(barrelAlphaBetaEcalHcal[bin]);
	barrelWithEcalCalib->addGraphPoints(barrelAlphaBetaEcalHcal[bin]);
      }
      if(barrelAlphaBetaHcal[bin]->computeAlphaBeta()) {
	barrelWithHcalCalib->addGraphPoints(barrelAlphaBetaHcal[bin]);
      }

      if(endcapAlphaBetaEcalHcal[bin]->computeAlphaBeta())
	{
	  endcapWithEcalHcalCalib->addGraphPoints(endcapAlphaBetaEcalHcal[bin]);
	  endcapWithEcalCalib->addGraphPoints(endcapAlphaBetaEcalHcal[bin]);
	}

     
      if(endcapAlphaBetaHcal[bin]->computeAlphaBeta()) //FIXME
	{
	  endcapWithHcalCalib->addGraphPoints(endcapAlphaBetaHcal[bin]);
	}
   }
   
   cout<<"Fitting alpha, beta coefficients..."<<endl;
   barrelWithEcalHcalCalib->initializeGraphs("alphabeta");
   barrelWithEcalCalib->initializeGraphs("alphabeta");
   barrelWithHcalCalib->initializeGraphs("alphabeta");   
   endcapWithEcalHcalCalib->initializeGraphs("alphabeta");
   endcapWithEcalCalib->initializeGraphs("alphabeta");
   endcapWithHcalCalib->initializeGraphs("alphabeta");   

   functionBarrelAlphaEcalHcal = new TF1("functionBarrelAlphaEcalHcal","[0]+[1]*exp(-x/[2])", 0, 1000);
   //functionBarrelBeta = new TF1("functionBarrelBeta","[0]+[1]*exp(-x/[2])", 0, 1000);
   //functionBarrelAlpha = new TF1("functionBarrelAlpha","[0]+(([1]+x*[3])*exp(-(x/[2])))", 0, 1000);
   functionBarrelBetaEcalHcal = new TF1("functionBarrelBetaEcalHcal","[0]+[1]*exp(-x/[2])", 0, 1000);

   functionBarrelAlphaHcal = new TF1("functionBarrelAlphaHcal","[0]+[1]*x", 0, 1000);// +[1]*exp(-x/[2])
   functionBarrelBetaHcal = new TF1("functionBarrelBetaHcal","[0]+[1]*exp(-x/[2])", 0, 1000);

   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+[1]*exp(-x/[2])", 0, 1000);
   functionEndcapBetaEcalHcal = new TF1("functionEndcapBetaEcalHcal","[0]+[1]*exp(-x/[2])",0,1000); //+[3]*[3]*exp(-x*x/([4]*[4]))

   functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+[1]*x", 0, 1000);// +[1]*exp(-x/[2])
   functionEndcapBetaHcal = new TF1("functionEndcapBetaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",0,1000);
   functionEndcapGammaHcal = new TF1("functionEndcapGammaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",0,1000);


   functionBarrelAlphaEcalHcal->SetParameters(0.02, -0.1, 200);
   // functionBarrelBeta->SetParameters(-0.02, 0.4, 200);
   //  functionBarrelAlpha->SetParameters(-2.33313e-02,-7.56070e+00,1.39070e+00,1.13114e+00);
   functionBarrelBetaEcalHcal->SetParameters(1.17842e-01,1.71167e-01,5.88921e+00);

   functionBarrelAlphaHcal->SetParameters(0.02, -0.1, 200);
   functionBarrelBetaHcal->SetParameters(1.17842e-01,1.71167e-01,5.88921e+00);

   functionEndcapAlphaEcalHcal->SetParameters(0.02, -0.1, 200);
   functionEndcapBetaEcalHcal->SetParameters(0.07, -2.5, 6.0, 0.3, 175.0);

   functionEndcapAlphaHcal->SetParameters(0.02, -0.1, 200);
   functionEndcapBetaHcal->SetParameters(0.07, -2.5, 6.0, 0.3, 175.0);
   functionEndcapGammaHcal->SetParameters(-1.71458e+00,9.61337e+00,2.94747e+00);

   barrelWithEcalHcalCalib->fitAlphasToFunction(functionBarrelAlphaEcalHcal);
   barrelWithEcalHcalCalib->fitAlphasToFunction();
   barrelWithEcalHcalCalib->fitBetasToFunction(functionBarrelBetaEcalHcal);
   barrelWithEcalHcalCalib->fitBetasToFunction();
   endcapWithEcalHcalCalib->fitAlphasToFunction(functionEndcapAlphaEcalHcal);
   endcapWithEcalHcalCalib->fitAlphasToFunction();
   endcapWithEcalHcalCalib->fitBetasToFunction(functionEndcapBetaEcalHcal);
   endcapWithEcalHcalCalib->fitBetasToFunction();

   barrelWithHcalCalib->fitAlphasToFunction(functionBarrelAlphaHcal);
   barrelWithHcalCalib->fitAlphasToFunction();
   barrelWithHcalCalib->fitBetasToFunction(functionBarrelBetaHcal);
   barrelWithHcalCalib->fitBetasToFunction();
   endcapWithHcalCalib->fitAlphasToFunction(functionEndcapAlphaHcal);
   endcapWithHcalCalib->fitAlphasToFunction();
   endcapWithHcalCalib->fitBetasToFunction(functionEndcapBetaHcal);
   endcapWithHcalCalib->fitBetasToFunction();
   endcapWithHcalCalib->fitGammasToFunction(functionEndcapGammaHcal);
   endcapWithHcalCalib->fitGammasToFunction();


   
   //Fill all the TH2's that can be put into drawGausFit in order to produce 
   //response and resolution plots.
   cout<<"Making response and resolution plots..."<<endl;
   for(unsigned entry = 0; entry < ETrueEnergies.size(); ++entry)
   {

      etrue = ETrueEnergies[entry];
      ecal = ecalEnergies[entry];
      hcal = hcalEnergies[entry];
      eta = abs(etas[entry]);

      if((ecal + hcal) < 0.5) continue;
      if( etrue < 1.0) continue;
      if( hcal == 0) continue;
      // if( ecal > 0) continue;
      
      if(fabs(eta) < 1.5) //alpha beta fit range for barrel
      {     
         raw->Fill(etrue, (ecal + hcal - etrue)/etrue);
         if(ecal > 0)
         {
	   correctedEta = barrelWithEcalHcalCalib->
	     getCalibratedEnergy(etrue, ecal, hcal, eta);
            
	   correctedE = barrelWithEcalHcalCalib->
	     getCalibratedEnergy(etrue, ecal, hcal);

	   corrEta->Fill(etrue, (correctedEta - etrue)/etrue);
	   corrEtaBarrel->Fill(etrue, (correctedEta - etrue)/etrue);
	   corrEtaBarrelEcalHcal->Fill(etrue, (correctedEta - etrue)/etrue);

	   rawEtaDependenceEH->Fill(eta, (ecal + hcal - etrue)/etrue);
	   corrEtaDependenceEH->Fill(eta, (correctedEta - etrue)/etrue);
	   hcorrEtaDependenceEH->Fill(eta, (correctedE - etrue)/etrue);
	   corrEtaDependenceProfEH->Fill(etrue, eta, (correctedEta - etrue)/etrue);
         }
         else 
	   {

	     correctedEta = barrelWithHcalCalib->
	       getCalibratedEnergy(etrue, ecal, hcal, eta);
   
	     correctedE = barrelWithHcalCalib->
	       getCalibratedEnergy(etrue, ecal, hcal);

	     if(etrue>7 && etrue<9) {
	       for(int k=0;k<1000;k++) {
		 float step=k/500.-0.99995;
		 float b= step;
		 float a = (etrue - 3.5 )/hcal - 1 -b*eta*eta;
		 bcplot->Fill(b,a);
	       }
	     }



	     corrEta->Fill(etrue, (correctedEta - etrue)/etrue);
	     corrEtaBarrel->Fill(etrue, (correctedEta - etrue)/etrue);
	     corrEtaBarrelHcal->Fill(etrue, (correctedEta - etrue)/etrue);

	     rawEtaDependenceH->Fill(eta, (ecal + hcal - etrue)/etrue);
	     corrEtaDependenceH->Fill(eta, (correctedEta - etrue)/etrue);
	     hcorrEtaDependenceH->Fill(eta, (correctedE - etrue)/etrue);
	     corrEtaDependenceProfH->Fill(etrue, eta, (correctedEta - etrue)/etrue);

         }
         if(fabs(eta) < 1.0) //b, c fit range
         {
            rawBarrel->Fill(etrue, (ecal + hcal - etrue)/etrue);
            
            if(ecal > 0)
            {
               correctedE = barrelWithEcalHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);

	       rawBarrelEcalHcal->Fill(etrue, (ecal + hcal  - etrue)/etrue );
               corrBarrel->Fill(etrue, (correctedE - etrue)/etrue);
               corrBarrelEcalHcal->Fill(etrue, (correctedE - etrue)/etrue);
               
	       // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);

               // rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
               // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);

	       // if(entry<5000) 
	       // 	 cout<<entry<<"   "<<eta<<"   "<<etrue<<"   "<<ecal+hcal<<"   "<<correctedE<<"   "<<correctedEta<<endl;

	     

            }
            else
	      {
		correctedE = barrelWithHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);

		rawBarrelHcal->Fill(etrue, ( ecal + hcal - etrue)/etrue );// (etrue-3.0)/(ecal+hcal) );//, 11936/(3917*sigmas[entry]*sigmas[entry]) );// ( ecal + hcal - etrue)/etrue );

		corrBarrel->Fill(etrue, (correctedE- etrue )/etrue);// 
		corrBarrelHcal->Fill(etrue, (correctedE- etrue )/etrue);

	    
	      }
         }
      }
      
      if(fabs(eta) < 2.4 && fabs(eta) > 1.55) //alpha beta fit range for endcap
      {
         raw->Fill(etrue, (ecal + hcal - etrue)/etrue);
         if(ecal > 0)
         {
            correctedEta = endcapWithEcalHcalCalib->
               getCalibratedEnergy(etrue, ecal, hcal, eta);

	    correctedE = endcapWithEcalHcalCalib->
	      getCalibratedEnergy(etrue, ecal, hcal);

            corrEta->Fill(etrue, (correctedEta - etrue)/etrue);
            corrEtaEndcap->Fill(etrue, (correctedEta - etrue)/etrue);
            corrEtaEndcapEcalHcal->Fill(etrue, (correctedEta - etrue)/etrue);

	    rawEtaDependenceEH->Fill(eta, (ecal + hcal - etrue)/etrue);
	    corrEtaDependenceEH->Fill(eta, (correctedEta - etrue)/etrue);
	    hcorrEtaDependenceEH->Fill(eta, (correctedE - etrue)/etrue); //FIXME
	    corrEtaDependenceProfEH->Fill(etrue, eta, (correctedEta - etrue)/etrue);

         }
         else
         {
            correctedEta = endcapWithHcalCalib->
               getCalibratedEnergy(etrue, ecal, hcal, eta);

	    correctedE = endcapWithHcalCalib->
	      getCalibratedEnergy(etrue, ecal, hcal);

            corrEta->Fill(etrue, (correctedEta - etrue)/etrue);            
            corrEtaEndcap->Fill(etrue, (correctedEta - etrue)/etrue);
            corrEtaEndcapHcal->Fill(etrue, (correctedEta - etrue)/etrue);

	    rawEtaDependenceH->Fill(eta, (ecal + hcal - etrue)/etrue);
	    corrEtaDependenceH->Fill(eta, (correctedEta - etrue)/etrue);
	    hcorrEtaDependenceH->Fill(eta, (correctedE - etrue)/etrue);
	    corrEtaDependenceProfH->Fill(etrue, eta, (correctedEta - etrue)/etrue);
	      
         }
         if(fabs(eta) < 2.2) //b, c fi trange
         {
            rawEndcap->Fill(etrue, (ecal + hcal - etrue)/etrue);
            
            if(ecal > 0)
            {

	      correctedEta = endcapWithEcalHcalCalib->
		getCalibratedEnergy(etrue, ecal, hcal, eta);

               correctedE = endcapWithEcalHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);

	       rawEndcapEcalHcal->Fill(etrue, (ecal + hcal - etrue)/etrue);
	       corrEndcap->Fill(etrue, (correctedE - etrue)/etrue);
	       corrEndcapEcalHcal->Fill(etrue, (correctedE - etrue)/etrue);

	   

	       // rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
	       // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);
	       // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);


            }
            else 
            {
               correctedE = endcapWithHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);
               
               rawEndcapHcal->Fill(etrue, (ecal + hcal - etrue)/etrue);
               corrEndcap->Fill(etrue, (correctedE - etrue)/etrue);
               corrEndcapHcal->Fill(etrue, (correctedE - etrue)/etrue);

             
            }
         }
      }  
   }

   ////////////////////////////////////////////////////////////////////////////
   //Add all the draw functions that you would like here, as well as any 
   //additional output you would like.
   ////////////////////////////////////////////////////////////////////////////

   cout<<" Now Summary "<<endl;
   
   //   exit(0);
   //drawGausFit(rawBarrelHcal,responseRaw,resolutionRaw);
   //drawGausFit(corrBarrelEcalHcal,responseCor,resolutionCor);
   //drawGausFit(corrBarrelHcal,responseCor,resolutionCor);
   //drawGausFit(corrEtaBarrelEcalHcal, responseEta, resolutionEta);

   //cout<<"Check pt 1"<<endl;
    
   //drawEtaDependence(hcorrEtaDependenceEH, responseEtaHCorrEtaEH);
   //drawEtaDependence(corrEtaDependenceEH, responseEtaEtaEH);
   drawEtaDependence(hcorrEtaDependenceH, responseEtaHCorrEtaH);
   drawEtaDependence(corrEtaDependenceH, responseEtaEtaH);
   //   drawGausFit(corrEta,response, resolution);
   //   drawCompare(responseRaw, response, resolutionRaw, resolution);

   // drawGausFit(rawEndcapHcal,responseRaw,resolutionRaw);
   //drawGausFit(corrEndcapHcal,responseCor,resolutionCor);
   //drawGausFit(corrEndcapEcalHcal,responseCor,resolutionCor);
   // drawGausFit(corrEtaEndcapHcal, responseEta, resolutionEta);

   barrelWithEcalHcalCalib->drawCoeffGraph("A","EH_barrel");
   barrelWithEcalHcalCalib->drawCoeffGraph("B", "EH_barrel");
   // barrelWithHcalCalib->drawCoeffGraph("C", "H_barrel");
   // barrelWithHcalCalib->drawCoeffGraph("Alpha","H_barrel");
   // barrelWithHcalCalib->drawCoeffGraph("Beta", "H_barrel");
   // barrelWithEcalHcalCalib->drawCoeffGraph("Alpha","EH_barrel");
   // barrelWithEcalHcalCalib->drawCoeffGraph("Beta", "EH_barrel");
   // barrelWithEcalHcalCalib->drawCoeffGraph("Alpha","EH");
   // barrelWithEcalHcalCalib->drawCoeffGraph("Beta", "EH");
   //endcapWithHcalCalib->drawCoeffGraph("Gamma", "H");
   //endcapWithEcalHcalCalib->drawCoeffGraph("A","EH_endcap");
   // endcapWithEcalHcalCalib->drawCoeffGraph("B", "EH_endcap");
   //endcapWithHcalCalib->drawCoeffGraph("C", "H_endcap");
   // endcapWithHcalCalib->drawCoeffGraph("Alpha","H_endcap");
   // endcapWithHcalCalib->drawCoeffGraph("Beta", "H_endcap");
   // endcapWithEcalHcalCalib->drawCoeffGraph("Alpha","EH_endcap");
   // endcapWithEcalHcalCalib->drawCoeffGraph("Beta", "EH_endcap");


   // TCanvas *cded = new TCanvas("cefzf","ceced");
   // bcplot->Draw("colz");
   // corrEtaDependenceProf->Draw("colz");

   
   functionBarrelEcalHcalB_e = functionBarrelEcalHcalB->GetTitle();
   functionBarrelEcalHcalC_e = functionBarrelEcalHcalC->GetTitle();
   functionBarrelHcalC_e = functionBarrelHcalC->GetTitle();
   functionBarrelAlphaEH_e = functionBarrelAlphaEcalHcal->GetTitle();
   functionBarrelBetaEH_e = functionBarrelBetaEcalHcal->GetTitle();

   functionBarrelAlphaH_e = functionBarrelAlphaHcal->GetTitle();
   functionBarrelBetaH_e = functionBarrelBetaHcal->GetTitle();


   functionEndcapEcalHcalB_e = functionEndcapEcalHcalB->GetTitle();
   functionEndcapEcalHcalC_e = functionEndcapEcalHcalC->GetTitle();
   functionEndcapHcalC_e = functionEndcapHcalC->GetTitle();
   functionEndcapAlphaEH_e = functionEndcapAlphaEcalHcal->GetTitle();
   functionEndcapBetaEH_e = functionEndcapBetaEcalHcal->GetTitle();
 
   functionEndcapAlphaH_e = functionEndcapAlphaHcal->GetTitle();
   functionEndcapBetaH_e = functionEndcapBetaHcal->GetTitle();


   //FUnction printing =============================================
   //takes more place, but easier to read

   //Thresholds first
   cout<<"  threshE = "<<barrelABCEcalHcal[2]->getA()<<";"<<endl;
   cout<<"  threshH = "<<barrelABCHcal[2]->getA()<<";"<<endl; 

   //Now functions : Ecal coef barrel
   cout<<"  faBarrel = new TF1(\"faBarrel\",\""<<
     functionBarrelEcalHcalB_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelEcalHcalB = functionBarrelEcalHcalB->GetParameter(i);
       if ( barrelEcalHcalB != 0. ) 
         cout<<"  faBarrel->SetParameter("<<i<<","<<barrelEcalHcalB<<");"<<endl;
     }
   // Hcal coef barrel for ecalHcal
   cout<<"  fbBarrel = new TF1(\"fbBarrel\",\""<<
     functionBarrelEcalHcalC_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
       if ( barrelEcalHcalC != 0. ) 
         cout<<"  fbBarrel->SetParameter("<<i<<","<<barrelEcalHcalC<<");"<<endl;
     }
   // Hcal coef barrel for Hcal
   cout<<"  fcBarrel = new TF1(\"fcBarrel\",\""<<
     functionBarrelHcalC_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelHcalC = functionBarrelHcalC->GetParameter(i);
       if ( barrelHcalC != 0. ) 
         cout<<"  fcBarrel->SetParameter("<<i<<","<<barrelHcalC<<");"<<endl;
     }
   //alpha function EcalHcal
   cout<<"  faEtaBarrelEH = new TF1(\"faEtaBarrelEH\",\""<<
     functionBarrelAlphaEH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelAlpha = functionBarrelAlphaEcalHcal->GetParameter(i);
       if ( barrelAlpha != 0. ) 
         cout<<"  faEtaBarrelEH->SetParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   //beta function EcalHcal
   cout<<"  fbEtaBarrelEH = new TF1(\"fbEtaBarrelEH\",\""<<
     functionBarrelBetaEH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelBeta = functionBarrelBetaEcalHcal->GetParameter(i);
       if ( barrelBeta != 0. ) 
         cout<<"  fbEtaBarrelEH->SetParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }
  //alpha function Hcal
   cout<<"  faEtaBarrelH = new TF1(\"faEtaBarrelH\",\""<<
     functionBarrelAlphaH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelAlpha = functionBarrelAlphaHcal->GetParameter(i);
       if ( barrelAlpha != 0. ) 
         cout<<"  faEtaBarrelH->SetParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   //beta function Hcal
   cout<<"  fbEtaBarrelH = new TF1(\"fbEtaBarrelH\",\""<<
     functionBarrelBetaH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       barrelBeta = functionBarrelBetaHcal->GetParameter(i);
       if ( barrelBeta != 0. ) 
         cout<<"  fbEtaBarrelH->SetParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }


   //Now endcaps (just a copy)
   //Now functions : Ecal coef Endcap
   cout<<"  faEndcap = new TF1(\"faEndcap\",\""<<
     functionEndcapEcalHcalB_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapEcalHcalB = functionEndcapEcalHcalB->GetParameter(i);
       if ( endcapEcalHcalB != 0. ) 
         cout<<"  faEndcap->SetParameter("<<i<<","<<endcapEcalHcalB<<");"<<endl;
     }
   // Hcal coef endcap for ecalHcal
   cout<<"  fbEndcap = new TF1(\"fbEndcap\",\""<<
     functionEndcapEcalHcalC_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapEcalHcalC = functionEndcapEcalHcalC->GetParameter(i);
       if ( endcapEcalHcalC != 0. ) 
         cout<<"  fbEndcap->SetParameter("<<i<<","<<endcapEcalHcalC<<");"<<endl;
     }
   // Hcal coef endcap for Hcal
   cout<<"  fcEndcap = new TF1(\"fcEndcap\",\""<<
     functionEndcapHcalC_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapHcalC = functionEndcapHcalC->GetParameter(i);
       if ( endcapHcalC != 0. ) 
         cout<<"  fcEndcap->SetParameter("<<i<<","<<endcapHcalC<<");"<<endl;
     }
   //alpha function for EcalHcal
   cout<<"  faEtaEndcapEH = new TF1(\"faEtaEndcapEH\",\""<<
     functionEndcapAlphaEH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapAlpha = functionEndcapAlphaEcalHcal->GetParameter(i);
       if ( endcapAlpha != 0. ) 
         cout<<"  faEtaEndcapEH->SetParameter("<<i<<","<<endcapAlpha<<");"<<endl;
     }
   //beta function for EcalHcal
   cout<<"  fbEtaEndcapEH = new TF1(\"fbEtaEndcapEH\",\""<<
     functionEndcapBetaEH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapBeta = functionEndcapBetaEcalHcal->GetParameter(i);
       if ( endcapBeta != 0. ) 
         cout<<"  fbEtaEndcapEH->SetParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }
 //alpha function for Hcal
   cout<<"  faEtaEndcapH = new TF1(\"faEtaEndcapH\",\""<<
     functionEndcapAlphaH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapAlpha = functionEndcapAlphaHcal->GetParameter(i);
       if ( endcapAlpha != 0. ) 
         cout<<"  faEtaEndcapH->SetParameter("<<i<<","<<endcapAlpha<<");"<<endl;
     }
   //beta function for Hcal
   cout<<"  fbEtaEndcapH = new TF1(\"fbEtaEndcapH\",\""<<
     functionEndcapBetaH_e<<"\",1.,1000.);"<<endl;
   for ( unsigned i = 0; i < 10; ++i ) 
     { 
       endcapBeta = functionEndcapBetaHcal->GetParameter(i);
       if ( endcapBeta != 0. ) 
         cout<<"  fbEtaEndcapH->SetParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }


  //  cout<<"  threshE = "<<barrelABCEcalHcal[2]->getA()<<";"<<endl;
  //  cout<<"  threshH = "<<barrelABCHcal[2]->getA()<<";"<<endl; 
   
  //  cout<<"  faBarrel = new TF1(\"faBarrel\",\""<<
  //     functionBarrelEcalHcalB_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fbBarrel = new TF1(\"fbBarrel\",\""<<
  //     functionBarrelEcalHcalC_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fcBarrel = new TF1(\"fcBarrel\",\""<<
  //     functionBarrelHcalC_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  faEtaBarrel = new TF1(\"faEtaBarrel\",\""<<
  //     functionBarrelAlpha_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fbEtaBarrel = new TF1(\"fbEtaBarrel\",\""<<
  //     functionBarrelBeta_e<<"\",1.,1000.);"<<endl;

  // for ( unsigned i = 0; i < 10; ++i ) 
  // { 
  //    barrelEcalHcalB = functionBarrelEcalHcalB->GetParameter(i);
  //    barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
  //    barrelHcalC = functionBarrelHcalC->GetParameter(i);
  //    barrelAlpha = functionBarrelAlpha->GetParameter(i);
  //    barrelBeta = functionBarrelBeta->GetParameter(i);
     
  //    if ( barrelEcalHcalB != 0. ) 
  //       cout<<"  faBarrel->SetParameter("<<i<<","<<barrelEcalHcalB<<");"<<endl;
  //    if ( barrelEcalHcalC != 0. ) 
  //       cout<<"  fbBarrel->SetParameter("<<i<<","<<barrelEcalHcalC<<");"<<endl;
  //    if ( barrelHcalC != 0. ) 
  //       cout<<"  fcBarrel->SetParameter("<<i<<","<<barrelHcalC<<");"<<endl;
  //    if ( barrelAlpha != 0. ) 
  //       cout<<"  faEtaBarrel->SetParameter("<<i<<","<<barrelAlpha<<");"<<endl;
  //    if ( barrelBeta != 0. ) 
  //       cout<<"  fbEtaBarrel->SetParameter("<<i<<","<<barrelBeta<<");"<<endl;
  // }

  //  cout<<"  faEndcap = new TF1(\"faEndcap\",\""<<
  //     functionEndcapEcalHcalB_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fbEndcap = new TF1(\"fbEndcap\",\""<<
  //     functionEndcapEcalHcalC_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fcEndcap = new TF1(\"fcEndcap\",\""<<
  //     functionEndcapHcalC_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  faEtaEndcap = new TF1(\"faEtaEndcap\",\""<<
  //     functionEndcapAlpha_e<<"\",1.,1000.);"<<endl;
  //  cout<<"  fbEtaEndcap = new TF1(\"fbEtaEndcap\",\""<<
  //     functionEndcapBeta_e<<"\",1.,1000.);"<<endl;

  // for ( unsigned i = 0; i < 10; ++i ) 
  // { 
  //    endcapEcalHcalB = functionEndcapEcalHcalB->GetParameter(i);
  //    endcapEcalHcalC = functionEndcapEcalHcalC->GetParameter(i);
  //    endcapHcalC = functionEndcapHcalC->GetParameter(i);
  //    endcapAlpha = functionEndcapAlphaEcalHcal->GetParameter(i);
  //    endcapBeta = functionEndcapBetaEcalHcal->GetParameter(i);
     
  //    if ( endcapEcalHcalB != 0. ) 
  //       cout<<"  faEndcap->SetParameter("<<i<<","<<endcapEcalHcalB<<");"<<endl;
  //    if ( endcapEcalHcalC != 0. ) 
  //       cout<<"  fbEndcap->SetParameter("<<i<<","<<endcapEcalHcalC<<");"<<endl;
  //    if ( endcapHcalC != 0. ) 
  //       cout<<"  fcEndcap->SetParameter("<<i<<","<<endcapHcalC<<");"<<endl;
  //    if ( endcapAlpha != 0. ) 
  //       cout<<"  faEtaEndcap->SetParameter("<<i<<","<<endcapAlpha<<");"<<endl;
  //    if ( endcapBeta != 0. ) 
  //       cout<<"  fbEtaEndcap->SetParameter("<<i<<","<<endcapBeta<<");"<<endl;
  // }

}

