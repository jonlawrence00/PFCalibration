#ifndef RecoParticleFlow_PFPatProducer_PFChargedHadronAnalyzer_
#define RecoParticleFlow_PFPatProducer_PFChargedHadronAnalyzer_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <math.h>

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"


/**\class PFChargedHadronAnalyzer 
\brief selects isolated charged hadrons from PF Charged Hadrons

\author Patrick Janot
\date   September 13, 2010
*/




class PFChargedHadronAnalyzer : public edm::EDAnalyzer {
 public:

  typedef reco::PFCandidateCollection::const_iterator CI;

  explicit PFChargedHadronAnalyzer(const edm::ParameterSet&);

  ~PFChargedHadronAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

 private:
  

  /// PFCandidates in which we'll look for pile up particles 
  edm::InputTag   inputTagPFCandidates_;
  edm::InputTag   inputTagPFSimParticles_;
  
  edm::InputTag   inputTagEcalPFClusters_;
  edm::EDGetTokenT<reco::PFCandidateCollection> tokenPFCandidates_;
  edm::EDGetTokenT<reco::PFSimParticleCollection>   tokenPFSimParticles_;
  edm::EDGetTokenT<reco::PFClusterCollection>   tokenEcalPFClusters_;

  /// Min pt for charged hadrons
  double ptMin_;
  
  /// Min p for charged hadrons
  double pMin_;

  /// Min hcal raw energy for charged hadrons
  double hcalMin_;
  
  /// Max ecal raw energy to define a MIP
  double ecalMax_;
  
  /// Min number of pixel hits for charged hadrons
  int nPixMin_;
  
  // isMInbias simulation
  bool isMBMC_;

  /// Min number of track hits for charged hadrons
  std::vector<int> nHitMin_;
  std::vector<double> nEtaMin_;
  
  // Number of tracks after cuts
  std::vector<unsigned int> nCh;
  std::vector<unsigned int> nEv;
  
  std::string outputfile_;
  TFile *tf1;
  TTree* s;
  
  float true_,p_,ecal_,hcal_,eta_,phi_,ho_;
  float etaEcal_,phiEcal_;
  int charge_;
  std::vector<float> dr_,Eecal_,Ehcal_,pfcID_;
  //bhumika Nov 2018
  std::vector<float> correcal_,corrhcal_;
  float Ccorrecal_,Ccorrhcal_;
  size_t orun,oevt,olumiBlock,otime;

  /* TH1F* h_phi = new TH1F("h_phi","phi w/o cut",90,-4.0,4.0); */


  /* TH2F* h_pix_phi = new TH2F("h_pix_phi","Npix vs phi",90,-4.0,4.0,10,0,10); */


  /* TH2F* h_pix_phi_valid_hits = new TH2F("h_pix_phi_valid_hits","Npix vs phi before iteration",90,-4.0,4.0,10,0,10); */

  /* TH2F* h_pix_phi_Barrel = new TH2F("h_pix_phi_barrel","Npix vs phi in barrel",90,-4.0,4.0,10,0,10); */
  /* TH2F* h_pix_phi_inTrack_EC = new TH2F("h_pix_phi_inTrack_EC","Npix vs phi EC,in tracker",90,-4.0,4.0,10,0,10); */
  /* TH2F* h_pix_phi_outTrack_EC = new TH2F("h_pix_phi_outTrack_EC","Npix vs phi EC,outside tracker",90,-4.0,4.0,10,0,10); */


  /* TH2F* h_hit_phi_Barrel = new TH2F("h_hit_phi_barrel","Nhit vs phi in barrel",90,-4.0,4.0,50,0,50); */
  /* TH2F* h_hit_phi_inTrack_EC = new TH2F("h_hit_phi_inTrack_EC","Nhit vs phi EC,in tracker",90,-4.0,4.0,50,0,50); */
  /* TH2F* h_hit_phi_outTrack_EC = new TH2F("h_hit_phi_outTrack_EC","Nhit vs phi EC,outside tracker",90,-4.0,4.0,10,0,10); */


  /* TH1F* h_phi_1 = new TH1F("h_phi_1","phi w/o cut",90,-4.0,4.0); */
  /* TH1F* h_phi_2 = new TH1F("h_phi_2","phi w/o cut",90,-4.0,4.0); */
  /* TH1F* h_phi_3 = new TH1F("h_phi_3","phi w/o cut",90,-4.0,4.0); */
  /* TH1F* h_phi_4 = new TH1F("h_phi_4","phi w/o cut",90,-4.0,4.0); */
  /* TH1F* h_phi_5 = new TH1F("h_phi_5","phi w/o cut",90,-4.0,4.0); */

  
  edm::RunNumber_t run;
  edm::EventNumber_t evt;
  edm::LuminosityBlockNumber_t lumiBlock;
  edm::Timestamp time;

  std::vector<float> addDr,addEmE,addHadE,addEta,addPhi;
  std::vector<int> addPdgId;
  std::vector<float> genDr,genE,genEta,genPhi;
  std::vector<int> genPdgId;

  std::vector<float> cluEcalE;
  std::vector<float> cluEcalEta;
  std::vector<float> cluEcalPhi;

  std::vector<float> distEcalTrk;


  std::vector<float> cluHcalE;
  std::vector<float> cluHcalEta;
  std::vector<float> cluHcalPhi;

  std::vector<float> distHcalTrk;
  std::vector<std::vector<float> > distHcalEcal;

  std::vector<int> pfcsID;

  std::vector<float> cluHadE;

  std::vector<std::vector<float> > emHitX; //eta for barrel
  std::vector<std::vector<float> > emHitY; //phi for barrel
  std::vector<std::vector<float> > emHitZ;
  std::vector<std::vector<float> > emHitE;
  std::vector<std::vector<float> > emHitF;

  std::vector<std::vector<float> > hadHitX;
  std::vector<std::vector<float> > hadHitY;
  std::vector<std::vector<float> > hadHitZ;
  std::vector<std::vector<float> > hadHitE;
  std::vector<std::vector<float> > hadHitF;

  //Basic clusters ECAL
  std::vector<float> bcEcalE;
  std::vector<float> bcEcalEta;
  std::vector<float> bcEcalPhi;

  //SimHits
  std::vector<float> EcalSimHits;
  std::vector<float> ESSimHits;
  std::vector<float> HcalSimHits;

  std::vector<float> EcalRecHits;
  std::vector<float> ESRecHits;
  std::vector<float> HcalRecHits;
  
  std::vector<float> EcalRecHitsDr;
  std::vector<float> ESRecHitsDr;
  std::vector<float> HcalRecHitsDr;
  

  const CaloGeometry*    theCaloGeom;


  //  void SaveRecHits(const edm::Event& iEvent, float eta_, float phi_);


  void SaveSimHit(const edm::Event& iEvent, float eta_, float phi_);
  float Eta( float theta);
  //  float dR( float eta1, float eta2, float phi1, float phi2);
  // float dPhi( float phi1, float phi2 );
  // float phi( float x, float y );



  /// verbose ?
  bool   verbose_;
  
  float dR(float eta1, float eta2, float phi1, float phi2 );


  float phi(float, float);
  float dPhi(float, float);


};

#endif
