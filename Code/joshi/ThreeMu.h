#ifndef ThreeMu_h
#define ThreeMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ThreeMu : public Selection {

 public:
  ThreeMu(TString Name_, TString id_);
  virtual ~ThreeMu();

  virtual double deltaR(double, double, double, double);
  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {L1SeedOk=0,HLTOk,PrimeVtx,isThreeMu,fitVtxChiSq,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

 // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  
  std::vector<TH1D> Muon1Muon2dR;
  std::vector<TH1D> Muon2Muon3dR;
  std::vector<TH1D> Muon1Muon3dR;
  std::vector<TH1D> TripleMass;
  std::vector<TH1D> Category;

  //3 Muon variables
		
	std::vector<TH1D> Muon1_isGlobal;
   std::vector<TH1D> Muon2_isGlobal;
   std::vector<TH1D> Muon3_isGlobal;
   std::vector<TH1D> Muon1_isStandAlone;
   std::vector<TH1D> Muon2_isStandAlone;
   std::vector<TH1D> Muon3_isStandAlone;
   std::vector<TH1D> Muon1_isTracker;
   std::vector<TH1D> Muon2_isTracker;
   std::vector<TH1D> Muon3_isTracker;
   std::vector<TH1D> Muon1_isCalo;
   std::vector<TH1D> Muon2_isCalo;
   std::vector<TH1D> Muon3_isCalo;
   std::vector<TH1D> Muon1_isIsolationValid;
   std::vector<TH1D> Muon2_isIsolationValid;
   std::vector<TH1D> Muon3_isIsolationValid;
   std::vector<TH1D> Muon1_isTimeValid;
   std::vector<TH1D> Muon2_isTimeValid;
   std::vector<TH1D> Muon3_isTimeValid;
   std::vector<TH1D> Muon1_emEt03;
   std::vector<TH1D> Muon2_emEt03;
   std::vector<TH1D> Muon3_emEt03;
   std::vector<TH1D> Muon1_emVetoEt03;
   std::vector<TH1D> Muon2_emVetoEt03;
   std::vector<TH1D> Muon3_emVetoEt03;
   std::vector<TH1D> Muon1_hadEt03;
   std::vector<TH1D> Muon2_hadEt03;
   std::vector<TH1D> Muon3_hadEt03;
   std::vector<TH1D> Muon1_hadVetoEt03;
   std::vector<TH1D> Muon2_hadVetoEt03;
   std::vector<TH1D> Muon3_hadVetoEt03;
   std::vector<TH1D> Muon1_nJets03;
   std::vector<TH1D> Muon2_nJets03;
   std::vector<TH1D> Muon3_nJets03;
   std::vector<TH1D> Muon1_nTracks03;
   std::vector<TH1D> Muon2_nTracks03;
   std::vector<TH1D> Muon3_nTracks03;
   std::vector<TH1D> Muon1_sumPt03;
   std::vector<TH1D> Muon2_sumPt03;
   std::vector<TH1D> Muon3_sumPt03;
   std::vector<TH1D> Muon1_trackerVetoPt03;
   std::vector<TH1D> Muon2_trackerVetoPt03;
   std::vector<TH1D> Muon3_trackerVetoPt03;
   std::vector<TH1D> Muon1_emEt05;
   std::vector<TH1D> Muon2_emEt05;
   std::vector<TH1D> Muon3_emEt05;
   std::vector<TH1D> Muon1_emVetoEt05;
   std::vector<TH1D> Muon2_emVetoEt05;
   std::vector<TH1D> Muon3_emVetoEt05;
   std::vector<TH1D> Muon1_hadEt05;
   std::vector<TH1D> Muon2_hadEt05;
   std::vector<TH1D> Muon3_hadEt05;
   std::vector<TH1D> Muon1_hadVetoEt05;
   std::vector<TH1D> Muon2_hadVetoEt05;
   std::vector<TH1D> Muon3_hadVetoEt05;
   std::vector<TH1D> Muon1_nJets05;
   std::vector<TH1D> Muon2_nJets05;
   std::vector<TH1D> Muon3_nJets05;
   std::vector<TH1D> Muon1_nTracks05;
   std::vector<TH1D> Muon2_nTracks05;
   std::vector<TH1D> Muon3_nTracks05;
   std::vector<TH1D> Muon1_sumPt05;
   std::vector<TH1D> Muon2_sumPt05;
   std::vector<TH1D> Muon3_sumPt05;
   std::vector<TH1D> Muon1_trackerVetoPt05;
   std::vector<TH1D> Muon2_trackerVetoPt05;
   std::vector<TH1D> Muon3_trackerVetoPt05;
   std::vector<TH1D> Muon1_sumChargedHadronPt03;
   std::vector<TH1D> Muon2_sumChargedHadronPt03;
   std::vector<TH1D> Muon3_sumChargedHadronPt03;
   std::vector<TH1D> Muon1_sumChargedParticlePt03;
   std::vector<TH1D> Muon2_sumChargedParticlePt03;
   std::vector<TH1D> Muon3_sumChargedParticlePt03;
   std::vector<TH1D> Muon1_sumNeutralHadronEt03;
   std::vector<TH1D> Muon2_sumNeutralHadronEt03;
   std::vector<TH1D> Muon3_sumNeutralHadronEt03;
   std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold03;
   std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold03;
   std::vector<TH1D> Muon3_sumNeutralHadronEtHighThreshold03;
   std::vector<TH1D> Muon1_sumPhotonEt03;
   std::vector<TH1D> Muon2_sumPhotonEt03;
   std::vector<TH1D> Muon3_sumPhotonEt03;
   std::vector<TH1D> Muon1_sumPhotonEtHighThreshold03;
   std::vector<TH1D> Muon2_sumPhotonEtHighThreshold03;
   std::vector<TH1D> Muon3_sumPhotonEtHighThreshold03;
   std::vector<TH1D> Muon1_sumPUPt03;
   std::vector<TH1D> Muon2_sumPUPt03;
   std::vector<TH1D> Muon3_sumPUPt03;
   std::vector<TH1D> Muon1_sumChargedHadronPt04;
   std::vector<TH1D> Muon2_sumChargedHadronPt04;
   std::vector<TH1D> Muon3_sumChargedHadronPt04;
   std::vector<TH1D> Muon1_sumChargedParticlePt04;
   std::vector<TH1D> Muon2_sumChargedParticlePt04;
   std::vector<TH1D> Muon3_sumChargedParticlePt04;
   std::vector<TH1D> Muon1_sumNeutralHadronEt04;
   std::vector<TH1D> Muon2_sumNeutralHadronEt04;
   std::vector<TH1D> Muon3_sumNeutralHadronEt04;
   std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold04;
   std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold04;
   std::vector<TH1D> Muon3_sumNeutralHadronEtHighThreshold04;
   std::vector<TH1D> Muon1_sumPhotonEt04;
   std::vector<TH1D> Muon2_sumPhotonEt04;
   std::vector<TH1D> Muon3_sumPhotonEt04;
   std::vector<TH1D> Muon1_sumPhotonEtHighThreshold04;
   std::vector<TH1D> Muon2_sumPhotonEtHighThreshold04;
   std::vector<TH1D> Muon3_sumPhotonEtHighThreshold04;
   std::vector<TH1D> Muon1_sumPUPt04;
   std::vector<TH1D> Muon2_sumPUPt04;
   std::vector<TH1D> Muon3_sumPUPt04;
   
	std::vector<TH1D> Isolation_NTracks;
   std::vector<TH1D> Isolation_RelPt;
   std::vector<TH1D> Isolation_MinDist;
   std::vector<TH1D> Isolation05_RelPt;
   std::vector<TH1D> Isolation05_NTracks;
   std::vector<TH1D> Isolation05_MinDist;
   std::vector<TH1D> Isolation_Ntrk1;
   std::vector<TH1D> Isolation_Ntrk2;
   std::vector<TH1D> Isolation_Ntrk3;
   std::vector<TH1D> Isolation_Ntrk0p1;
   std::vector<TH1D> Isolation_Ntrk0p2;
   std::vector<TH1D> Isolation_Ntrk0p5;
   std::vector<TH1D> Isolation_maxdxy;
	
};
#endif
