//#include "RecoParticleFlow/Configuration/test/PFChargedHadronAnalyzer.h"
#include "PFCalibration/PFChargedHadronAnalyzer/plugins/PFChargedHadronAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
//#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include <TROOT.h>
#include <TVector3.h>

//#include "PFChargedHadronAnalyzer.h"
using namespace std;
using namespace edm;
using namespace reco;

PFChargedHadronAnalyzer::PFChargedHadronAnalyzer(const edm::ParameterSet& iConfig) {

  nCh = std::vector<unsigned int>(11,static_cast<unsigned int>(0));
  nEv = std::vector<unsigned int>(2,static_cast<unsigned int>(0));

  inputTagPFCandidates_
    = iConfig.getParameter<InputTag>("PFCandidates");
  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(inputTagPFCandidates_);

  //std::cout << "Check point 1 " << std::endl;
  inputTagEcalPFRechit_ = iConfig.getParameter<InputTag>("EcalPFrechit");
  tokenEcalPFRechit_ = consumes<reco::PFRecHitCollection>(inputTagEcalPFRechit_);

  inputTagHcalPFRechit_ = iConfig.getParameter<InputTag>("HcalPFrechit");
  tokenHcalPFRechit_ = consumes<reco::PFRecHitCollection>(inputTagHcalPFRechit_);

  inputTagPFSimParticles_
    = iConfig.getParameter<InputTag>("PFSimParticles");
  tokenPFSimParticles_ = consumes<reco::PFSimParticleCollection>(inputTagPFSimParticles_);

  inputTagEcalPFClusters_
    = iConfig.getParameter<InputTag>("EcalPFClusters");
  tokenEcalPFClusters_ = consumes<reco::PFClusterCollection>(inputTagEcalPFClusters_);

  // Smallest track pt
  ptMin_ = iConfig.getParameter<double>("ptMin");

  // Smallest track p
  pMin_ = iConfig.getParameter<double>("pMin");

  // Smallest raw HCAL energy linked to the track
  hcalMin_ = iConfig.getParameter<double>("hcalMin");

  // Largest ECAL energy linked to the track to define a MIP
  ecalMax_ = iConfig.getParameter<double>("ecalMax");

  // Smallest number of pixel hits
  nPixMin_ = iConfig.getParameter<int>("nPixMin");

  // Smallest number of track hits in different eta ranges
  nHitMin_ = iConfig.getParameter< std::vector<int> > ("nHitMin");
  nEtaMin_ = iConfig.getParameter< std::vector<double> > ("nEtaMin");

  // Is minbias from simulation
  isMBMC_ = iConfig.getUntrackedParameter<bool>("isMinBiasMC",false);

  // Consider only leading PFNeutralHadron (for neutral case.)
  leadingPFNHOnly_ = iConfig.getUntrackedParameter<bool>("leadingPFNHOnly",false);

  // Consider only up to one leading PFPhoton (for neutral case.)
  leadingPFPhotonOnly_ = iConfig.getUntrackedParameter<bool>("leadingPFPhotonOnly",false);

  verbose_ =
    iConfig.getUntrackedParameter<bool>("verbose",false);

  LogDebug("PFChargedHadronAnalyzer")
    <<" input collection : "<<inputTagPFCandidates_ ;


  // The root tuple
  outputfile_ = iConfig.getParameter<std::string>("rootOutputFile");
  tf1 = new TFile(outputfile_.c_str(), "RECREATE");
  s = new TTree("s"," PFCalibration");

  s->Branch("true",&true_,"true/F");
  s->Branch("p",&p_,"p/F");
  s->Branch("ecal",&ecal_,"ecal/F");
  s->Branch("hcal",&hcal_,"hcal/F");
  s->Branch("ho",&ho_,"ho/F");
  s->Branch("eta",&eta_,"eta/F");
  s->Branch("phi",&phi_,"phi/F");
  s->Branch("charge",&charge_,"charge/I");

  //spandey Apr_27 dR
  s->Branch("dr",&dr_);
  s->Branch("Eecal",&Eecal_);
  s->Branch("Ehcal",&Ehcal_);
  s->Branch("pfcID",&pfcID_);

  s->Branch("pfcs",&pfcsID);

  //by bhumika Nov 2018
  s->Branch("correcal",&correcal_);
  s->Branch("corrhcal",&corrhcal_);
  s->Branch("Ccorrecal",&Ccorrecal_,"Ccorrecal/F");
  s->Branch("Ccorrhcal",&Ccorrhcal_,"Ccorrhcal/F");

  //PF RecHits added (by dR cut)
  s->Branch("EcalRechit_posx", &EcalRechit_posx_);
  s->Branch("EcalRechit_posy", &EcalRechit_posy_);
  s->Branch("EcalRechit_posz", &EcalRechit_posz_);
  s->Branch("EcalRechit_E", &EcalRechit_E_);
  s->Branch("EcalRechit_dr", &EcalRechit_dr_);
  s->Branch("EcalRechit_eta", &EcalRechit_eta_);
  s->Branch("EcalRechit_phi", &EcalRechit_phi_);
  s->Branch("EcalRechit_totE", &EcalRechit_totE_,"EcalRechit_totE/F");
  s->Branch("EcalRechit_depth", &EcalRechit_depth_);
  s->Branch("EcalPFclustereta", &EcalPFclustereta_);

  s->Branch("HcalRechit_posx", &HcalRechit_posx_);
  s->Branch("HcalRechit_posy", &HcalRechit_posy_);
  s->Branch("HcalRechit_posz", &HcalRechit_posz_);
  s->Branch("HcalRechit_E", &HcalRechit_E_);
  s->Branch("HcalRechit_dr", &HcalRechit_dr_);
  s->Branch("HcalRechit_eta", &HcalRechit_eta_);
  s->Branch("HcalRechit_phi", &HcalRechit_phi_);
  s->Branch("HcalRechit_totE", &HcalRechit_totE_,"HcalRechit_totE/F");
  s->Branch("HcalRechit_depth", &HcalRechit_depth_);
  s->Branch("HcalPFclustereta", &HcalPFclustereta_);

  // PF clusters based on pfc.elementsInBlocks
  s->Branch("cluEcalE",&cluEcalE);
  s->Branch("cluEcalEta",&cluEcalEta);
  s->Branch("cluEcalPhi",&cluEcalPhi);
  s->Branch("cluEcalEtot",&cluEcalEtot);

  s->Branch("cluHcalE",&cluHcalE);
  s->Branch("cluHcalEta",&cluHcalEta);
  s->Branch("cluHcalPhi",&cluHcalPhi);
  s->Branch("cluHcalEtot",&cluHcalEtot);

  // PFRecHit based on pfc.elementsInBlocks & clusterRef()->recHitFractions()
  s->Branch("emHitX",&emHitX );
  s->Branch("emHitY",&emHitY );
  s->Branch("emHitZ",&emHitZ );
  s->Branch("emHitE",&emHitE );
  s->Branch("emHitF",&emHitF );
  s->Branch("emHitEtot",&emHitEtot );

  s->Branch("hadHitX",&hadHitX );
  s->Branch("hadHitY",&hadHitY );
  s->Branch("hadHitZ",&hadHitZ );
  s->Branch("hadHitE",&hadHitE );
  s->Branch("hadHitF",&hadHitF );
  s->Branch("hadHitEtot",&hadHitEtot );

  s->Branch("run",&orun,"orun/l");
  s->Branch("evt",&oevt,"orun/l");
  s->Branch("lumiBlock",&olumiBlock,"orun/l");
  s->Branch("time",&otime,"orun/l");

}



PFChargedHadronAnalyzer::~PFChargedHadronAnalyzer() {

  std::cout << "Total number of events .............. " << nEv[0] << std::endl;
  std::cout << "Number of events with 1 Sim Particle  " << nEv[1] << std::endl;


  std::cout << "Number of PF candidates ............. " << nCh[0] << std::endl;
  std::cout << "Number of PF Charged Hadrons......... " << nCh[1] << std::endl;
  std::cout << " - With pt > " << ptMin_ << " GeV/c ................ " << nCh[2] << std::endl;
  std::cout << " - With E_HCAL > " << hcalMin_ << " GeV .............. " << nCh[3] << std::endl;
  std::cout << " - With only 1 track in the block ... " << nCh[4] << std::endl;
  std::cout << " - With p > " << pMin_ << " GeV/c ................. " << nCh[5] << std::endl;
  std::cout << " - With at least " << nPixMin_ << " pixel hits ....... " << nCh[6] << std::endl;
  std::cout << " - With more than "<< nHitMin_[0] << " track hits ..... " << nCh[7] << std::endl;
  std::cout << " - With E_ECAL < " << ecalMax_ << " GeV ............ " << nCh[8] << std::endl;
  std::cout << " ============================================= "<< std::endl;
  std::cout << "Number of PF rechit associated with PF Charged Hadrons......... " << nCh[9] << std
::endl;
  std::cout << "Number of PF rechit associated with PF neutral Hadrons......... " << nCh[10] << std
::endl;

  tf1->cd();
  s->Write();

  h_true_barrel->Write();
  h_true_ec_in->Write();
  h_true_ec_out->Write();
  h_true_hf->Write();

  tf1->Write();
  tf1->Close();


}



void
PFChargedHadronAnalyzer::beginRun(const edm::Run& run,
				  const edm::EventSetup & es) { }


void
PFChargedHadronAnalyzer::analyze(const Event& iEvent,
				 const EventSetup& iSetup) {

  LogDebug("PFChargedHadronAnalyzer")<<"START event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;

  /*
   edm::ESHandle<CaloGeometry> pCalo;
   iSetup.get<CaloGeometryRecord>().get( pCalo );
   theCaloGeom = pCalo.product();
  */


  run  = iEvent.id().run();
  evt  = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  time = iEvent.time();

  orun = (size_t)run;
  oevt = (size_t)evt;
  olumiBlock = (size_t)lumiBlock;
  otime = (size_t)((iEvent.time().value())>>32);


  // get PFCandidates
  Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByToken(tokenPFCandidates_, pfCandidates);

  //get Ecal PFClusters
  Handle<reco::PFClusterCollection> pfClustersEcal;
  iEvent.getByToken(tokenEcalPFClusters_, pfClustersEcal);

  Handle<PFSimParticleCollection> trueParticles;
  bool isSimu = iEvent.getByToken(tokenPFSimParticles_, trueParticles);

  Handle<reco::PFRecHitCollection> pfRechitEcal;
  iEvent.getByToken(tokenEcalPFRechit_, pfRechitEcal);

  Handle<reco::PFRecHitCollection> pfRechitHcal;
  iEvent.getByToken(tokenHcalPFRechit_, pfRechitHcal);

  //simHits
  EcalSimHits.clear();
  ESSimHits.clear();
  HcalSimHits.clear();

  //recHits
  EcalRecHits.clear();
  ESRecHits.clear();
  HcalRecHits.clear();
  EcalRecHitsDr.clear();
  ESRecHitsDr.clear();
  HcalRecHitsDr.clear();

  pfcsID.clear();

  charge_=0;
  dr_.clear();
  Eecal_.clear();
  Ehcal_.clear();
  pfcID_.clear();

  //bhumika Nov 2018
  correcal_.clear();
  corrhcal_.clear();
  //

  if(isMBMC_)
    isSimu=false;

  if ( isSimu ) {
    nEv[0]++;//  cout<<" True part size "<<(*trueParticles).size()<<"    "
    // 		  <<(*trueParticles)[0].pdgCode()<<"   "<<(*trueParticles)[1].pdgCode()<<endl;
    if ( (*trueParticles).size() != 1 ) return;
    nEv[1]++;

    // Check if there is a reconstructed track

    bool isCharged = false;
    for( CI ci  = pfCandidates->begin();
	 ci!=pfCandidates->end(); ++ci)  {
      const reco::PFCandidate& pfc = *ci;
      //if ( pfc.particleId() == 5 )
	pfcsID.push_back( pfc.particleId() );
      // std::cout << "Id = " << pfc.particleId() << std::endl;
      if ( pfc.particleId() < 4 ) {
	isCharged = true;
	break;
      }
    }

    //
    // (1) Non-charged case (no track)
    //
    if ( !isCharged ) {

      reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
      const reco::PFTrajectoryPoint& tpatecal = ((*trueParticles)[0]).extrapolatedPoint( ecalEntrance );
      eta_ = tpatecal.positionREP().Eta();
      if ( fabs(eta_) < 1E-10 ) return;
      phi_ = tpatecal.positionREP().Phi();
      true_ = std::sqrt(tpatecal.momentum().Vect().Mag2());
      p_ = 0.;
      charge_=0;
      ecal_ = 0.;
      hcal_ = 0.;
      Ccorrecal_=0.;
      Ccorrhcal_=0.;
      dr_.clear();  //spandey Apr_27 dR
      Eecal_.clear();
      Ehcal_.clear();
      pfcID_.clear();

      //bhumika Nov 2018
      correcal_.clear();
      corrhcal_.clear();
      //
      //cout<<"***********************"<<endl;
      int a=0;
      int b=0;

      // return; // KenH: comment out when just focusing on the case with a reconstructed track
      // KenH
      double emHitEtotTmp=0, hadHitEtotTmp=0, cluEcalEtotTmp=0, cluHcalEtotTmp=0;
	//CG
	  emHitF.clear();
	  emHitE.clear();
	  emHitX.clear();
	  emHitY.clear();
	  emHitZ.clear();
	  emHitIndex.clear();
	  hadHitF.clear();
	  hadHitE.clear();
	  hadHitX.clear();
	  hadHitY.clear();
	  hadHitZ.clear();
	  hadHitIndex.clear();

      // KenH - Prepare to possibly select only minDR PFPhoton and PFNeutralHadron
      int nPFNH=0;
      int nPFPhoton=0;
      float pfPhotonMinDR=99.;
      float pfNHMinDR=99.;
      int pfPhotonMinDRIndex=-1;
      int pfNHMinDRIndex=-1;
      int iPFCand=0;
      for( CI ci  = pfCandidates->begin();
	   ci!=pfCandidates->end(); ++ci)  {
	const reco::PFCandidate& pfc = *ci;
	iPFCand++;
	double deta = eta_ - pfc.eta();
	double dphi = dPhi(phi_, pfc.phi() );
	double dR = std::sqrt(deta*deta+dphi*dphi);
	if ( pfc.particleId() == 4 && dR < 0.2 ) {
	  nPFPhoton++;
	  if (pfPhotonMinDR>dR){
	    pfPhotonMinDR=dR;
	    pfPhotonMinDRIndex=(unsigned int)iPFCand;
	  }
	}
	if ( pfc.particleId() == 5 && dR < 0.4 ) {
	  nPFNH++;
	  if (pfNHMinDR>dR) {
	    pfNHMinDR=dR;
	    pfNHMinDRIndex=(unsigned int)iPFCand;
	  }
	}
      }
      /*
      if (nPFNH>1 || nPFPhoton>0){
	for( CI ci  = pfCandidates->begin();
	     ci!=pfCandidates->end(); ++ci)  {
	  const reco::PFCandidate& pfc = *ci;
	  double deta = eta_ - pfc.eta();
	  double dphi = dPhi(phi_, pfc.phi() );
	  double dR = std::sqrt(deta*deta+dphi*dphi);
	  if (( pfc.particleId() == 4 && dR < 0.2 ) ||
	      ( pfc.particleId() == 5 && dR < 0.4 )){
	    std::cout << dR << " " << pfc.particleId() << " "
		      << pfc.p() << " " << pfc.energy() << " "
		      << pfc.eta() << " " << pfc.phi() << " "
		      << pfc.rawEcalEnergy() << " " << pfc.ecalEnergy() << " "
		      << pfc.rawHcalEnergy() << " " << pfc.hcalEnergy() << " "
		      << std::endl;

	    const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();
	    if( theElements.empty() ) continue;
	    for (unsigned int j=0; j<theElements.size(); j++){
	      const reco::PFBlockRef blockRef = theElements[j].first;
	      PFBlock::LinkData linkData =  blockRef->linkData();
	      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	      unsigned iEle = theElements[j].second;
	      PFBlockElement::Type type = elements[iEle].type();
	      std::cout << "elements size: " << elements.size() << " " << theElements[j].second << " " << theElements.size() << " " << type << std::endl;
	      if (type == PFBlockElement::TRACK){
		const reco::PFBlockElementTrack& et =
		  dynamic_cast<const reco::PFBlockElementTrack &>( elements[iEle] );
		double p = et.trackRef()->p();
		std::cout << "track p: " << p << std::endl;
	      }
	    }
	  }
	}
      }
      */
      // KenH ends

      for( CI ci  = pfCandidates->begin();
	   ci!=pfCandidates->end(); ++ci)  {
	iPFCand++;

	const reco::PFCandidate& pfc = *ci;

	if (leadingPFPhotonOnly_ &&
	    pfc.particleId() == 4 &&
	    (pfPhotonMinDRIndex>=0) &&
	    iPFCand != pfPhotonMinDRIndex) continue;
	if (leadingPFNHOnly_ &&
	    pfc.particleId() == 5 &&
	    (pfNHMinDRIndex>=0) &&
	    iPFCand != pfNHMinDRIndex) continue;

	double deta = eta_ - pfc.eta();
	double dphi = dPhi(phi_, pfc.phi() );
	double dR = std::sqrt(deta*deta+dphi*dphi);
	if ( dR < 1.2 ) {
	  //spandey Apr_27 dR
	  dr_.push_back(dR);
	  pfcID_.push_back(pfc.particleId());
	  Eecal_.push_back(pfc.rawEcalEnergy());
	  Ehcal_.push_back(pfc.rawHcalEnergy());
	  //bhumika Nov 2018
	  correcal_.push_back(pfc.ecalEnergy());
	  corrhcal_.push_back(pfc.hcalEnergy());
	  //
	}

	if ( pfc.particleId() == 4 && dR < 0.2 ) ecal_ += pfc.rawEcalEnergy();
	//KenH if ( pfc.particleId() == 5 && dR < 0.4 ) hcal_ += pfc.rawHcalEnergy();
	if ( pfc.particleId() == 5 && dR < 0.4 ){
	  ecal_ += pfc.rawEcalEnergy();
	  hcal_ += pfc.rawHcalEnergy();
	}

	if(pfc.particleId() == 5 && dR < 0.4 && a==0){
	  a++; // counter
	  double Ecalrechit_en=0, Hcalrechit_en=0;
	  for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitEcal->begin(); it != pfRechitEcal->end(); ++it) {
            double deta_ = eta_ - it->positionREP().eta();
            double dphi_ = dPhi(phi_, it->positionREP().phi() );
            double dR_ = std::sqrt(deta_*deta_+dphi_*dphi_);
            if(dR_<=0.1)
	      {
		Ecalrechit_en+=it->energy();
	      }
	  }

	  for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitHcal->begin(); it != pfRechitHcal->end(); ++it) {
	    double deta_ = eta_ - it->positionREP().eta();
	    double dphi_ = dPhi(phi_, it->positionREP().phi() );
	    double dR_ = std::sqrt(deta_*deta_+dphi_*dphi_);
	    if(dR_<=0.2)
	      {
		Hcalrechit_en+=it->energy();
	      }
	  }
	} // if(pfc.particleId() == 5 && dR < 0.4 && a==1)
	//if((Ecalrechit_en + Hcalrechit_en)==0) continue;

	// Consider only nearby PF photon and PF hadron
	if (( pfc.particleId() == 4 && dR < 0.2 ) ||
	    ( pfc.particleId() == 5 && dR < 0.4 )){

/**
	  emHitF.clear();
	  emHitE.clear();
	  emHitX.clear();
	  emHitY.clear();
	  emHitZ.clear();
	  emHitIndex.clear();
	  hadHitF.clear();
	  hadHitE.clear();
	  hadHitX.clear();
	  hadHitY.clear();
	  hadHitZ.clear();
	  hadHitIndex.clear();
**/

	  const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();
	  if( theElements.empty() ) continue;

	    const reco::PFBlockRef blockRef = theElements[0].first;
	    PFBlock::LinkData linkData =  blockRef->linkData();
	    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();

	    // Check that there is only one track in the block.
	    unsigned int nTracks = 0;
	    unsigned int nEcal = 0;
	    unsigned int nHcal = 0;
	    unsigned iTrack = 999;
	    vector<unsigned> iECAL;// =999;
	    vector<unsigned> iHCAL;// =999;

	    //KenH for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    for(unsigned jEle=0; jEle<theElements.size(); jEle++) {

	      const reco::PFBlockRef blockRef2 = theElements[jEle].first;
	      // Ensure different PFElements of a PFCandidate belongs to the same PFBlock
	      if (blockRef.id()!=blockRef2.id() || blockRef.key()!=blockRef2.key() || blockRef.index()!=blockRef2.index() )
		std::cout << blockRef.id() << " "<< blockRef2.id() << " "
			  << blockRef.key() << " " << blockRef2.key() << " "
			  << blockRef.index() << " " << blockRef2.index() << std::endl;
	      unsigned iEle = theElements[jEle].second;

	      // Find the tracks in the block
	      PFBlockElement::Type type = elements[iEle].type();

	      // for type definition
	      // https://cmssdt.cern.ch/lxr/source/DataFormats/ParticleFlowReco/interface/PFBlockElement.h

	      switch( type ) {
	      case PFBlockElement::TRACK:
		iTrack = iEle;
		nTracks++;
		break;
	      case PFBlockElement::ECAL:
		iECAL.push_back( iEle );
		//cout<<"iEle "<<iEle<<endl;
		nEcal++;
		break;
	      case PFBlockElement::HCAL:
		iHCAL.push_back( iEle );
		nHcal++;
		break;
	      default:
		continue;
	      }

	    }

	    //ECAL element
	    // Loop over ECAL clusters
	    for(unsigned int ii=0;ii<nEcal;ii++) {
	      const reco::PFBlockElementCluster& eecal =
		dynamic_cast<const reco::PFBlockElementCluster &>( elements[ iECAL[ii] ] );
	      double E_ECAL = eecal.clusterRef()->energy();
	      double eta_ECAL = eecal.clusterRef()->eta();
	      double phi_ECAL = eecal.clusterRef()->phi();

	      cluEcalE.push_back( E_ECAL );
	      cluEcalEta.push_back( eta_ECAL );
	      cluEcalPhi.push_back( phi_ECAL );
	      cluEcalEtotTmp += E_ECAL;

	      double d = blockRef->dist(iTrack, iECAL[ii], linkData);
	      distEcalTrk.push_back( d );

	      if(isMBMC_ || isSimu) {
	      const std::vector< reco::PFRecHitFraction > erh=eecal.clusterRef()->recHitFractions();
	      // cout<<"Number of Rechits: "<<erh.size()<<endl;
	      // Loop over associated PFRecHits
	      for(unsigned int ieh=0;ieh<erh.size();ieh++) {

		auto it = find(emHitIndex.begin(), emHitIndex.end(), erh[ieh].recHitRef().index());
		if (it == emHitIndex.end()) { // new PFRecHit
		  emHitIndex.push_back( erh[ieh].recHitRef().index() );
		  emHitF.push_back( erh[ieh].fraction() );
		  emHitE.push_back( erh[ieh].recHitRef()->energy() );
		  emHitX.push_back( erh[ieh].recHitRef()->position().x() );
		  emHitY.push_back( erh[ieh].recHitRef()->position().y() );
		  emHitZ.push_back( erh[ieh].recHitRef()->position().z() );
		} else { // already-filled PFRecHit
		  //std::cout << "PFRecHit ECAL already stored" << std::endl;
		  unsigned int index = it - emHitIndex.begin();
		  //std::cout << index << " " << emHitF[index] << std::endl;
		  emHitF[index] += erh[ieh].fraction();
		  //std::cout << index << " " << emHitF[index] << std::endl;
		}
		emHitEtotTmp += erh[ieh].recHitRef()->energy() * erh[ieh].fraction();

	      }
	      }
	    } // nEcal loop
	    cluEcalEtot = cluEcalEtotTmp;
	    emHitEtot = emHitEtotTmp;

	    //HCAL element
	    // Loop over HCAL clusters
	      for(unsigned int ii=0;ii<nHcal;ii++) {
		const reco::PFBlockElementCluster& ehcal =
		  dynamic_cast<const reco::PFBlockElementCluster &>( elements[iHCAL[ii] ] );
		double E_HCAL = ehcal.clusterRef()->energy();
		double eta_HCAL = ehcal.clusterRef()->eta();
		double phi_HCAL = ehcal.clusterRef()->phi();

		cluHcalE.push_back( E_HCAL );
		cluHcalEta.push_back( eta_HCAL );
		cluHcalPhi.push_back( phi_HCAL );
		cluHcalEtotTmp += E_HCAL;

		double d = blockRef->dist(iTrack, iHCAL[ii], linkData);
		distHcalTrk.push_back( d );

		//ECAL-HCAL distance
		vector<float> tmp;
		distHcalEcal.push_back(tmp);
		for(unsigned int ij=0;ij<nEcal;ij++) {
		  d = blockRef->dist(iECAL[ij], iHCAL[ii], linkData);
		  distHcalEcal[ii].push_back( d );
		}

		if(isMBMC_ || isSimu) {
		  const std::vector< reco::PFRecHitFraction > erh=ehcal.clusterRef()->recHitFractions();
		  // Loop over associated PFRecHits
		  for(unsigned int ieh=0;ieh<erh.size();ieh++) {

		    auto it = find(hadHitIndex.begin(), hadHitIndex.end(), erh[ieh].recHitRef().index());
		    if (it == hadHitIndex.end()) { // new PFRecHit
		      hadHitIndex.push_back( erh[ieh].recHitRef().index() );
		      hadHitF.push_back( erh[ieh].fraction() );
		      hadHitE.push_back( erh[ieh].recHitRef()->energy() );
		      hadHitX.push_back( erh[ieh].recHitRef()->position().x() );
		      hadHitY.push_back( erh[ieh].recHitRef()->position().y() );
		      hadHitZ.push_back( erh[ieh].recHitRef()->position().z() );
		    } else { // already-filled PFRecHit
		      unsigned int index = it - hadHitIndex.begin();
		      // std::cout << "PFRecHit HCAL already stored" << std::endl;
		      // std::cout << index << " " << hadHitF[index] << " " << erh[ieh].fraction() << " "
		      // 		<< hadHitX[index] << " " << erh[ieh].recHitRef()->position().x() << " "
		      // 		<< hadHitY[index] << " " << erh[ieh].recHitRef()->position().y() << " "
		      // 		<< std::endl;
		      hadHitF[index] += erh[ieh].fraction();
		      //std::cout << index << " " << hadHitF[index] << std::endl;
		    }
		    hadHitEtotTmp += erh[ieh].recHitRef()->energy() * erh[ieh].fraction();

		  }
		}
	      }
	    cluHcalEtot = cluHcalEtotTmp;
	    hadHitEtot = hadHitEtotTmp;
        }

	// Now loop over any PFRecHits and select based on dR
	if(pfc.particleId() == 5 && dR < 0.4 && b==0){
	  b++;

	  EcalRechit_posx_.clear();
	  EcalRechit_posy_.clear();
	  EcalRechit_posz_.clear();
	  EcalRechit_E_.clear();
	  HcalRechit_posx_.clear();
	  HcalRechit_posy_.clear();
	  HcalRechit_posz_.clear();
	  HcalRechit_E_.clear();
	  EcalRechit_dr_.clear();
	  EcalRechit_depth_.clear();
	  HcalRechit_dr_.clear();
	  HcalRechit_depth_.clear();
	  EcalRechit_eta_.clear();
	  EcalRechit_phi_.clear();
	  HcalRechit_eta_.clear();
	  HcalRechit_phi_.clear();
	  EcalPFclustereta_.clear();
	  HcalPFclustereta_.clear();

	  double Ecalrechit_energy=0, Hcalrechit_energy=0;

	  for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitEcal->begin(); it != pfRechitEcal->end(); ++it) {
	    double deta_ = eta_ - it->positionREP().eta();
	    double dphi_ = dPhi(phi_, it->positionREP().phi() );
	    double dR_ = std::sqrt(deta_*deta_+dphi_*dphi_);
	    if(dR_<=0.1){ //Bhumika
	      Ecalrechit_energy+=it->energy();
	      EcalRechit_posx_.push_back(it->position().x());
	      EcalRechit_posy_.push_back(it->position().y());
	      EcalRechit_posz_.push_back(it->position().z());
	      EcalRechit_E_.push_back(it->energy());
	      EcalRechit_dr_.push_back(dR_);
	      EcalRechit_eta_.push_back(it->positionREP().eta());
	      EcalRechit_phi_.push_back(it->positionREP().phi());
	      EcalRechit_depth_.push_back(it->depth());
	      EcalPFclustereta_.push_back(eta_);
	    }
	  }
	  EcalRechit_totE_=Ecalrechit_energy;

	  for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitHcal->begin(); it != pfRechitHcal->end(); ++it) {
	    double deta_ = eta_ - it->positionREP().eta();
	    double dphi_ = dPhi(phi_, it->positionREP().phi() );
	    double dR_ = std::sqrt(deta_*deta_+dphi_*dphi_);
	    if(dR_<=0.2){ //Bhumika
	      Hcalrechit_energy+=it->energy();
	      HcalRechit_posx_.push_back(it->position().x());
	      HcalRechit_posy_.push_back(it->position().y());
	      HcalRechit_posz_.push_back(it->position().z());
	      HcalRechit_E_.push_back(it->energy());
	      HcalRechit_dr_.push_back(dR_);
	      HcalRechit_eta_.push_back(it->positionREP().eta());
	      HcalRechit_phi_.push_back(it->positionREP().phi());
	      HcalRechit_depth_.push_back(it->depth());
	      HcalPFclustereta_.push_back(eta_);
	    }
	  }
	  HcalRechit_totE_=Hcalrechit_energy;

	  //cout<<"Size of EcalRechit = "<<EcalRechit_posx_.size()<<" , Size of HcalRechit = "<<HcalRechit_posx_.size()<<endl;
	  //cout<<"True energy = "<<true_<<" , EcalRechit energy = "<<Ecalrechit_energy<<" ,  HcalRechit energy = "<<Hcalrechit_energy<<endl;
	  //cout<<"Size of ecal supercluster rechits = "<<sc_nRechits<<endl;
	  //cout<<"Total Ecal energy = "<<sc_totalE<<endl;
	  nCh[10]++;
	} // if(pfc.particleId() == 5 && dR < 0.4 && a==0){
      } // loop over PFCandidates
      if(eta_<1.5) h_true_barrel->Fill(true_);
      else if(eta_>=1.5 && eta_<2.5) h_true_ec_in->Fill(true_);
      else if(eta_>=2.5 && eta_<2.75) h_true_ec_out->Fill(true_);
      else h_true_hf->Fill(true_);

      s->Fill();
      return;
    }

  }

  //
  // (2) Track case
  //

  // Case of a reconstructed track.
  // Loop on pfCandidates
  for( CI ci  = pfCandidates->begin();
       ci!=pfCandidates->end(); ++ci)  {

    // The pf candidate
    const reco::PFCandidate& pfc = *ci;
    nCh[0]++;

    // Only charged hadrons (no PF muons, no PF electrons) 1 / 5
    if ( pfc.particleId() != 1 ) continue;
    nCh[1]++;

    // Charged hadron minimum pt (the track pt, to an excellent approximation)
    if ( pfc.pt() < ptMin_ ) continue;
    nCh[2]++;

    // At least 1 GeV in HCAL
    double ecalRaw = pfc.rawEcalEnergy();
    double hcalRaw = pfc.rawHcalEnergy();
    double hoRaw = pfc.rawHoEnergy();
    double ecalcorr = pfc.ecalEnergy();
    double hcalcorr = pfc.hcalEnergy();

    if ( ecalRaw + hcalRaw < hcalMin_ ) continue;
    nCh[3]++;

    // Find the corresponding PF block elements
    const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();
    if( theElements.empty() ) continue;
    double Ecalrechit_energy=0, Hcalrechit_energy=0, Ecalrechit_en=0, Hcalrechit_en=0,ecalrechit=0,hcalrechit=0;
    // KenH
    double emHitEtotTmp=0, hadHitEtotTmp=0, cluEcalEtotTmp=0, cluHcalEtotTmp=0;
	    emHitF.clear();
	    emHitE.clear();
	    emHitX.clear();
	    emHitY.clear();
	    emHitZ.clear();
	    hadHitF.clear();
	    hadHitE.clear();
	    hadHitX.clear();
	    hadHitY.clear();
	    hadHitZ.clear();

	    const reco::PFBlockRef blockRef = theElements[0].first;
	    PFBlock::LinkData linkData =  blockRef->linkData();
	    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();

	    // Check that there is only one track in the block.
	    unsigned int nTracks = 0;
	    unsigned int nEcal = 0;
	    unsigned int nHcal = 0;
	    unsigned iTrack = 999;
	    vector<unsigned> iECAL;// =999;
	    vector<unsigned> iHCAL;// =999;

	    //KenH for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    for(unsigned jEle=0; jEle<theElements.size(); jEle++) {

	      const reco::PFBlockRef blockRef2 = theElements[jEle].first;
	      // Ensure different PFElements of a PFCandidate belongs to the same PFBlock
	      if (blockRef.id()!=blockRef2.id() || blockRef.key()!=blockRef2.key() || blockRef.index()!=blockRef2.index() )
		std::cout << blockRef.id() << " "<< blockRef2.id() << " "
			  << blockRef.key() << " " << blockRef2.key() << " "
			  << blockRef.index() << " " << blockRef2.index() << std::endl;
	      unsigned iEle = theElements[jEle].second;

	      // Find the tracks in the block
	      PFBlockElement::Type type = elements[iEle].type();

	      // for type definition
	      // https://cmssdt.cern.ch/lxr/source/DataFormats/ParticleFlowReco/interface/PFBlockElement.h

	      switch( type ) {
	      case PFBlockElement::TRACK:
		iTrack = iEle;
		nTracks++;
		break;
	      case PFBlockElement::ECAL:
		iECAL.push_back( iEle );
		nEcal++;
		break;
	      case PFBlockElement::HCAL:
		iHCAL.push_back( iEle );
		nHcal++;
		break;
	      default:
		continue;
	      }

	    }

	    //KenH
	    //std::cout << "bbb: " << nTracks << " " << nEcal << " " << nHcal << std::endl;

	    if ( nTracks != 1 ) continue;
	    nCh[4]++;

	    // Characteristics of the track
	    const reco::PFBlockElementTrack& et =
	      dynamic_cast<const reco::PFBlockElementTrack &>( elements[iTrack] );
	    double p =/*pfc.energy();//*/ et.trackRef()->p();
	    double pt =/*pfc.pt();//*/ et.trackRef()->pt();
	    double eta =/*pfc.eta();//*/ et.trackRef()->eta();
	    double phi =/*pfc.phi();//*/ et.trackRef()->phi();

	    // A minimum p and pt
	    if ( p < pMin_ || pt < ptMin_ ) continue;
	    nCh[5]++;

	    // Count the number of valid hits (first three iteration only)
	    //unsigned int nHits = et.trackRef()->found();
	    unsigned int tobN = 0;
	    unsigned int tecN = 0;
	    unsigned int tibN = 0;
	    unsigned int tidN = 0;
	    unsigned int pxbN = 0;
	    unsigned int pxdN = 0;
	    unsigned int validPix = 0;
	    unsigned int validTrackerHits = 0;
	    const reco::HitPattern& hp = et.trackRef()->hitPattern();

	    // h_pix_phi_valid_hits->Fill(phi,hp.numberOfValidPixelHits());
	    validPix = hp.numberOfValidPixelHits();
	    validTrackerHits = hp.numberOfValidTrackerHits();

	    switch ( et.trackRef()->algo() ) {
	    case TrackBase::initialStep:
	      tobN += hp.numberOfValidStripTOBHits();
	      tecN += hp.numberOfValidStripTECHits();
	      tibN += hp.numberOfValidStripTIBHits();
	      tidN += hp.numberOfValidStripTIDHits();
	      pxbN += hp.numberOfValidPixelBarrelHits();
	      pxdN += hp.numberOfValidPixelEndcapHits();
	      //validPix += hp.numberOfValidPixelHits();
	      break;
	    case TrackBase::lowPtQuadStep:
	    case TrackBase::highPtTripletStep:
	    case TrackBase::lowPtTripletStep:
	      // tobN += hp.numberOfValidStripTOBHits();
	      // tecN += hp.numberOfValidStripTECHits();
	      // tibN += hp.numberOfValidStripTIBHits();
	      // tidN += hp.numberOfValidStripTIDHits();
	      // pxbN += hp.numberOfValidPixelBarrelHits();
	      // pxdN += hp.numberOfValidPixelEndcapHits();
	      // validPix += hp.numberOfValidPixelHits();
	      // break;
	    case TrackBase::detachedQuadStep:
	    case TrackBase::detachedTripletStep:
	    case TrackBase::pixelPairStep:
	    case TrackBase::mixedTripletStep:
	    case TrackBase::pixelLessStep:
	    case TrackBase::tobTecStep:
	    case TrackBase::jetCoreRegionalStep:
	    default:
	      break;
	    }
	    //int inner = pxbN+pxdN;
	    int inner = validPix;

	    //int outer = tibN+tobN+tidN+tecN;
	    int outer = validTrackerHits;

	    // // Number of pixel hits
	      if ( inner < nPixMin_ ) continue;
	      nCh[6]++;

	    // Number of tracker hits (eta-dependent cut)
	    bool trackerHitOK = false;
	    double etaMin = 0.;
	    for ( unsigned int ieta=0; ieta<nEtaMin_.size(); ++ieta ) {
	      if ( fabs(eta) < etaMin ) break;
	      double etaMax = nEtaMin_[ieta];
	      trackerHitOK =
		fabs(eta)>etaMin && fabs(eta)<etaMax && inner+outer>nHitMin_[ieta];
	      if ( trackerHitOK ) break;
	      etaMin = etaMax;
	    }
	      if ( !trackerHitOK ) continue;
	      nCh[7]++;

	    // Selects only ECAL MIPs
	    if ( ecalRaw > ecalMax_ ) continue;
	    nCh[8]++;


	    //extrapolate track to ECAL --> impact position
	    etaEcal_ = et.positionAtECALEntrance().Eta();
	    phiEcal_ = et.positionAtECALEntrance().Phi();

	    // Fill the root-tuple
	    p_ = p;
	    ecal_ = ecalRaw;
	    hcal_ = hcalRaw;
	    ho_ = hoRaw;
	    Ccorrecal_ = ecalcorr;
	    Ccorrhcal_ = hcalcorr;
	    charge_ = pfc.charge();

	    if( isSimu ) {
	      reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
	      const reco::PFTrajectoryPoint& tpatecal = ((*trueParticles)[0]).extrapolatedPoint( ecalEntrance );
	      eta_ = tpatecal.positionREP().Eta();
	      phi_ = tpatecal.positionREP().Phi();
	      true_ = std::sqrt(tpatecal.momentum().Vect().Mag2());
	      if(eta<1.5) h_true_barrel->Fill(std::sqrt(tpatecal.momentum().Vect().Mag2()));
	      else if(eta>=1.5 && eta<2.5) h_true_ec_in->Fill(std::sqrt(tpatecal.momentum().Vect().Mag2()));
	      else if(eta>=2.5 && eta<2.75) h_true_ec_out->Fill(std::sqrt(tpatecal.momentum().Vect().Mag2()));
	      else h_true_hf->Fill(std::sqrt(tpatecal.momentum().Vect().Mag2()));
	    }
	    else {
	      eta_ = eta; // tpatecal.positionREP().Eta();
	      phi_ = phi; //tpatecal.positionREP().Phi();
	      true_ = p; //std::sqrt(tpatecal.momentum().Vect().Mag2());
	      if(eta<1.5) h_true_barrel->Fill(p);
	      else if(eta>=1.5 && eta<2.5) h_true_ec_in->Fill(p);
	      else if(eta>=2.5 && eta<2.75) h_true_ec_out->Fill(p);
	      else h_true_hf->Fill(p);
	    }

	    for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitEcal->begin(); it != pfRechitEcal->end(); ++it) {
	      double deta = eta_ - it->positionREP().eta();
	      double dphi = dPhi(phi_, it->positionREP().phi() );
	      double dR = std::sqrt(deta*deta+dphi*dphi);
	      if(dR<=0.1)
		{
		  Ecalrechit_en+=it->energy();
		  ecalrechit++;
		}
	    }

	    for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitHcal->begin(); it != pfRechitHcal->end(); ++it) {
	      double deta = eta_ - it->positionREP().eta();
	      double dphi = dPhi(phi_, it->positionREP().phi() );
	      double dR = std::sqrt(deta*deta+dphi*dphi);
	      if(dR<=0.2)
	      {
		Hcalrechit_en+=it->energy();
		hcalrechit++;
	      }
	    }

	    if((Ecalrechit_en + Hcalrechit_en)==0) {
	      //      cout<<" Total Raw PFrechit energy = "<<(Ecalrechit_en+Hcalrechit_en)<<endl;
	      continue;
	    }
	    if(ecalrechit==0 && hcalrechit==0) continue;

	    EcalRechit_posx_.clear();
	    EcalRechit_posy_.clear();
	    EcalRechit_posz_.clear();
	    EcalRechit_E_.clear();
	    HcalRechit_posx_.clear();
	    HcalRechit_posy_.clear();
	    HcalRechit_posz_.clear();
	    HcalRechit_E_.clear();
	    EcalRechit_dr_.clear();
	    EcalRechit_depth_.clear();
	    HcalRechit_dr_.clear();
	    HcalRechit_depth_.clear();
	    EcalRechit_eta_.clear();
	    EcalRechit_phi_.clear();
	    HcalRechit_eta_.clear();
	    HcalRechit_phi_.clear();
	    EcalPFclustereta_.clear();
	    HcalPFclustereta_.clear();

	    //ECAL elements
            // Loop over ECAL clusters
	    for(unsigned int ii=0;ii<nEcal;ii++) {
	      const reco::PFBlockElementCluster& eecal =
		dynamic_cast<const reco::PFBlockElementCluster &>( elements[ iECAL[ii] ] );
	      double E_ECAL = eecal.clusterRef()->energy();
	      double eta_ECAL = eecal.clusterRef()->eta();
	      double phi_ECAL = eecal.clusterRef()->phi();

	      cluEcalE.push_back( E_ECAL );
	      cluEcalEta.push_back( eta_ECAL );
	      cluEcalPhi.push_back( phi_ECAL );
	      cluEcalEtotTmp += E_ECAL;

	      double d = blockRef->dist(iTrack, iECAL[ii], linkData);
	      distEcalTrk.push_back( d );

	      if(isMBMC_ || isSimu) {
		const std::vector< reco::PFRecHitFraction > erh=eecal.clusterRef()->recHitFractions();
		//cout<<"Number of Rechits: "<<erh.size()<<endl;
		//
		// Loop over associated PFRecHits
		for(unsigned int ieh=0;ieh<erh.size();ieh++) {
		emHitF.push_back( erh[ieh].fraction() );
		emHitE.push_back( erh[ieh].recHitRef()->energy() );
		emHitX.push_back( erh[ieh].recHitRef()->position().x() );
		emHitY.push_back( erh[ieh].recHitRef()->position().y() );
		emHitZ.push_back( erh[ieh].recHitRef()->position().z() );
		emHitEtotTmp += erh[ieh].recHitRef()->energy() * erh[ieh].fraction();
	      }
	      }
	    }
	    cluEcalEtot = cluEcalEtotTmp;
	    emHitEtot = emHitEtotTmp;

	    //HCAL element
            // Loop over HCAL clusters
	    for(unsigned int ii=0;ii<nHcal;ii++) {
	      const reco::PFBlockElementCluster& ehcal =
		dynamic_cast<const reco::PFBlockElementCluster &>( elements[iHCAL[ii] ] );
	      double E_HCAL = ehcal.clusterRef()->energy();
	      double eta_HCAL = ehcal.clusterRef()->eta();
	      double phi_HCAL = ehcal.clusterRef()->phi();

	      cluHcalE.push_back( E_HCAL );
	      cluHcalEta.push_back( eta_HCAL );
	      cluHcalPhi.push_back( phi_HCAL );
	      cluHcalEtotTmp += E_HCAL;

	      double d = blockRef->dist(iTrack, iHCAL[ii], linkData);
	      distHcalTrk.push_back( d );

	      //ECAL-HCAL distance
	      vector<float> tmp;
	      distHcalEcal.push_back(tmp);
	      for(unsigned int ij=0;ij<nEcal;ij++) {
		d = blockRef->dist(iECAL[ij], iHCAL[ii], linkData);
		distHcalEcal[ii].push_back( d );
	      }

	      if(isMBMC_ || isSimu) {
	      const std::vector< reco::PFRecHitFraction > erh=ehcal.clusterRef()->recHitFractions();
	      //cout<<"Number of Rechits: "<<erh.size()<<endl;
	      //
	      // Loop over associated PFRecHits
	      for(unsigned int ieh=0;ieh<erh.size();ieh++) {

		hadHitF.push_back( erh[ieh].fraction() );
		hadHitE.push_back( erh[ieh].recHitRef()->energy() );
		hadHitX.push_back( erh[ieh].recHitRef()->position().x() );
		hadHitY.push_back( erh[ieh].recHitRef()->position().y() );
		hadHitZ.push_back( erh[ieh].recHitRef()->position().z() );
		hadHitEtotTmp += erh[ieh].recHitRef()->energy() * erh[ieh].fraction();

	      }
	      }

	    }
	    cluHcalEtot = cluHcalEtotTmp;
	    hadHitEtot = hadHitEtotTmp;

    //Basic Ecal PFRecHits by dR
    for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitEcal->begin(); it != pfRechitEcal->end(); ++it) {
      double deta = eta_ - it->positionREP().eta();
      double dphi = dPhi(phi_, it->positionREP().phi() );
      double dR = std::sqrt(deta*deta+dphi*dphi);
      if(dR<=0.1){
	Ecalrechit_energy+=it->energy();
	EcalRechit_posx_.push_back(it->position().x());
	EcalRechit_posy_.push_back(it->position().y());
	EcalRechit_posz_.push_back(it->position().z());
	EcalRechit_E_.push_back(it->energy());
	EcalRechit_dr_.push_back(dR);
	EcalRechit_eta_.push_back(it->positionREP().eta());
	EcalRechit_phi_.push_back(it->positionREP().phi());
	EcalRechit_depth_.push_back(it->depth());
	EcalPFclustereta_.push_back(eta_);
      }
    }
    EcalRechit_totE_=Ecalrechit_energy;

    for (std::vector<reco::PFRecHit>::const_iterator it = pfRechitHcal->begin(); it != pfRechitHcal->end(); ++it) {
      double deta = eta_ - it->positionREP().eta();
      double dphi = dPhi(phi_, it->positionREP().phi() );
      double dR = std::sqrt(deta*deta+dphi*dphi);
      if(dR<=0.2){
	Hcalrechit_energy+=it->energy();
	HcalRechit_posx_.push_back(it->position().x());
	HcalRechit_posy_.push_back(it->position().y());
	HcalRechit_posz_.push_back(it->position().z());
	HcalRechit_E_.push_back(it->energy());
	HcalRechit_dr_.push_back(dR);
	HcalRechit_eta_.push_back(it->positionREP().eta());
	HcalRechit_phi_.push_back(it->positionREP().phi());
	HcalRechit_depth_.push_back(it->depth());
	HcalPFclustereta_.push_back(eta_);
      }
    }
    HcalRechit_totE_=Hcalrechit_energy;

    nCh[9]++;

    // KenH
    /*
    std::cout << ecal_ << " "  // rawECAL. It has been used for PF hadron calibration
	      << cluEcalEtot << " "
	      << emHitEtot << " "
	      << EcalRechit_totE_ << " " // deltaR-matched
	      << hcal_ << " " // rawHCAL. It has been used for PF hadron calibration
	      << cluHcalEtot << " "
	      << hadHitEtot << " "
	      << HcalRechit_totE_ << " " // deltaR-matched
	      << std::endl;
    */

    s->Fill();

    addDr.clear();
    addPdgId.clear();
    addEmE.clear();
    addHadE.clear();
    addEta.clear();
    addPhi.clear();

    cluEcalE.clear();
    cluEcalEta.clear();
    cluEcalPhi.clear();

    distEcalTrk.clear();

    cluHcalE.clear();
    cluHcalEta.clear();
    cluHcalPhi.clear();

    distHcalTrk.clear();
    distHcalEcal.clear();

    genDr.clear();
    genPdgId.clear();
    genE.clear();
    genEta.clear();
    genPhi.clear();


    bcEcalE.clear();
    bcEcalEta.clear();
    bcEcalPhi.clear();

  } // Looping over PF candidates (and select a good charged hadron candidate)

}


float PFChargedHadronAnalyzer::dR(float eta1, float eta2, float phi1, float phi2 ) {

  TVector3 v1(0,0,0),v2(0,0,0);

  v1.SetPtEtaPhi(1, eta1, phi1);
  v2.SetPtEtaPhi(1, eta2, phi2);

  return v1.DrEtaPhi( v2 );

}


void PFChargedHadronAnalyzer::SaveSimHit(const edm::Event& iEvent,  float eta_, float phi_) {

  //Access to simHits informations
  Handle<PCaloHitContainer> h_PCaloHitsEB;
  iEvent.getByLabel("g4SimHits","EcalHitsEB", h_PCaloHitsEB);

  Handle<PCaloHitContainer> h_PCaloHitsEE;
  iEvent.getByLabel("g4SimHits","EcalHitsEE", h_PCaloHitsEE);

  Handle<PCaloHitContainer> h_PCaloHitsES;
  iEvent.getByLabel("g4SimHits","EcalHitsES", h_PCaloHitsES);

  Handle<PCaloHitContainer> h_PCaloHitsH;
  iEvent.getByLabel("g4SimHits","HcalHits", h_PCaloHitsH);

  //iterator
  PCaloHitContainer::const_iterator genSH;

  //match hits... dR 0.2, should contains all simHits

  //ECAL
  if( fabs(eta_) <1.5 ) { //barrel

    for(genSH = h_PCaloHitsEB->begin(); genSH != h_PCaloHitsEB->end(); genSH++) {
      // float theta = genSH->thetaAtEntry();
      // float phi = genSH->phiAtEntry();
      // float eta = Eta( theta );
      // float dr = dR( eta, eta_, phi, phi_ );

      // if(dr > 0.2 ) continue;
      //cout<<" ecal hit : "<<genSH->energy()<<endl;
      EcalSimHits.push_back( genSH->energy() );
    }
  }
  else {

    for(genSH = h_PCaloHitsEE->begin(); genSH != h_PCaloHitsEE->end(); genSH++) {
      // float theta = genSH->thetaAtEntry();
      // float phi = genSH->phiAtEntry();
      // float eta = Eta( theta );
      // float dr = dR( eta, eta_, phi, phi_ );

      // if(dr > 0.2 ) continue;
      EcalSimHits.push_back( genSH->energy() );
    }

    for(genSH = h_PCaloHitsES->begin(); genSH != h_PCaloHitsES->end(); genSH++) {
       // float theta = genSH->thetaAtEntry();
       // float phi = genSH->phiAtEntry();
       // float eta = Eta( theta );
       // float dr = dR( eta, eta_, phi, phi_ );

       // if(dr > 0.2 ) continue;
       ESSimHits.push_back( genSH->energy() );
    }
  }

  //Hcal
  float sH=0;
  for(genSH = h_PCaloHitsH->begin(); genSH != h_PCaloHitsH->end(); genSH++) {
    // float theta = genSH->thetaAtEntry();
    // float phi = genSH->phiAtEntry();
    // float eta = Eta( theta );
    // float dr = dR( eta, eta_, phi, phi_ );
       // if(dr > 0.2 ) continue;
    sH += genSH->energy();
    //cout<<" ecal hit : "<<genSH->energy()<<"    "<<genSH->energyEM()<<"   "<<genSH->energyHad()<<"   "<<sH<<endl;
       HcalSimHits.push_back( genSH->energy() );
    }

}


float PFChargedHadronAnalyzer::Eta( float theta_ ) {
  if( sin(theta_/2.)==0 ) return 10000.*cos(theta_/2.);
  return -log(tan(theta_/2.0));
}

/*
void PFChargedHadronAnalyzer::SaveRecHits(const edm::Event& iEvent, float eta_, float phi_) {

  //get rechits

  edm::Handle< EcalRecHitCollection > ebRecHits_h;
  edm::Handle< EcalRecHitCollection > eeRecHits_h;
  edm::Handle< EcalRecHitCollection > esRecHits_h;
 // Barrel
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsEB", ebRecHits_h );
  // Endcaps
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsEE", eeRecHits_h );
  // Preshower
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsES", esRecHits_h );
  // Hcal
  edm::Handle< HBHERecHitCollection > hbheRecHits_h;
  iEvent.getByLabel( "hbhereco","", hbheRecHits_h );


  for( size_t ii =0; ii < ebRecHits_h->size(); ++ii )
    {
      EcalRecHitRef recHitRef( ebRecHits_h, ii );
      EBDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.1 ) continue;
      //cout<<"EB : "<<dr<<"  "<<recHitRef->energy()<<endl;
      EcalRecHits.push_back( recHitRef->energy() );
      EcalRecHitsDr.push_back( dr );
    }

  for( size_t ii =0; ii < eeRecHits_h->size(); ++ii )
    {
      EcalRecHitRef recHitRef( eeRecHits_h, ii );
      EEDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.1 ) continue;
      EcalRecHits.push_back( recHitRef->energy() );
      EcalRecHitsDr.push_back( dr );
    }


  for( size_t ii =0; ii < hbheRecHits_h->size(); ++ii )
    {
      HBHERecHitRef recHitRef( hbheRecHits_h, ii );
      HcalDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.15 ) continue;
      //cout<<"Hcal : "<<dr<<"  "<<recHitRef->energy()<<endl;
      HcalRecHits.push_back( recHitRef->energy() );
      HcalRecHitsDr.push_back( dr );
    }



}
*/
float
PFChargedHadronAnalyzer::phi( float x, float y ) {
  float phi_ =atan2(y, x);
  return (phi_>=0) ?  phi_ : phi_ + 2*3.141592;
}

float
PFChargedHadronAnalyzer::dPhi( float phi1, float phi2 )
{
  float phi1_= phi( cos(phi1), sin(phi1) );
  float phi2_= phi( cos(phi2), sin(phi2) );
  float dphi_= phi1_-phi2_;
  if( dphi_> 3.141592 ) dphi_-=2*3.141592;
  if( dphi_<-3.141592 ) dphi_+=2*3.141592;
  return dphi_;
}






// float PFChargedHadronAnalyzer::dPhi( float phi1, float phi2 )
// {
//   float phi1_= phi( cos(phi1), sin(phi1) );
//   float phi2_= phi( cos(phi2), sin(phi2) );
//   float dphi_= phi1_-phi2_;
//   if( dphi_> 3.141592 ) dphi_-=2*3.141592;
//   if( dphi_<-3.141592 ) dphi_+=2*3.141592;
//   return dphi_;
// }

// float
// PFChargedHadronAnalyzer::phi( float x, float y )
// {
//   float phi_ =atan2(y, x);
//   return (phi_>=0) ?  phi_ : phi_ + 2*3.141592;
// }

// float PFChargedHadronAnalyzer::dR(float eta1, float eta2, float phi1, float phi2) {

//   float deta = eta1-eta2;
//   float dphi = dPhi( phi1, phi2 );

//   return sqrt( pow( deta, 2) + pow( dphi, 2) );

// }


DEFINE_FWK_MODULE(PFChargedHadronAnalyzer);
