#ifndef DsPhiPeak_h
#define DsPhiPeak_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class DsPhiPeak : public Selection {

 public:
  DsPhiPeak(TString Name_, TString id_);
  virtual ~DsPhiPeak();

  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {TriggerOk,TwoMuTrkCandidate,OSMuons,Mu1PtCut,Mu2PtCut,MuonID,MuMuMassCut,TrackPtCut,NTrackHits,ChiSqCut,TriggerMatch,GenMatch,NCuts}; 

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  float PhiMassHigh;
  float PhiMassLow;
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  
  std::vector<TH1D> MuonsPtRatio;
  std::vector<TH1D> DimuondR;
  std::vector<TH1D> Muon1TrkdR;
  std::vector<TH1D> Muon2TrkdR;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> TripleMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D> Category;

  //Muon variables
  
  //Track candidate variables
  std::vector<TH1D> Track_P;
  std::vector<TH1D> Track_E;
  std::vector<TH1D> Track_Pt;
  std::vector<TH1D> Track_Eta;
  std::vector<TH1D> Track_Phi;
  std::vector<TH1D> Track_vx;
  std::vector<TH1D> Track_vy;
  std::vector<TH1D> Track_vz;
  std::vector<TH1D> Track_normalizedChi2;
  std::vector<TH1D> Track_numberOfValidHits;
  std::vector<TH1D> Track_charge;
  std::vector<TH1D> Track_dxy;
  std::vector<TH1D> Track_dz;
  std::vector<TH1D> Track_dxyError;
  std::vector<TH1D> Track_dzError;

//Dimuon variables

std::vector<TH1D> Muon1_Pt;
std::vector<TH1D> Muon1_E;
std::vector<TH1D> Muon1_P;
std::vector<TH1D> Muon1_Eta;
std::vector<TH1D> Muon1_Phi;
std::vector<TH1D> Muon1_vx;
std::vector<TH1D> Muon1_vy;
std::vector<TH1D> Muon1_vz;

std::vector<TH1D> Muon2_Pt;
std::vector<TH1D> Muon2_E;
std::vector<TH1D> Muon2_P;
std::vector<TH1D> Muon2_Eta;
std::vector<TH1D> Muon2_Phi;
std::vector<TH1D> Muon2_vx;
std::vector<TH1D> Muon2_vy;
std::vector<TH1D> Muon2_vz;
		
  std::vector<TH1D> Muon1_isGlobal;
  std::vector<TH1D> Muon2_isGlobal;
  std::vector<TH1D> Muon1_isStandAlone;
  std::vector<TH1D> Muon2_isStandAlone;
  std::vector<TH1D> Muon1_isTracker;
  std::vector<TH1D> Muon2_isTracker;
  std::vector<TH1D> Muon1_isCalo;
  std::vector<TH1D> Muon2_isCalo;
  std::vector<TH1D> Muon1_isIsolationValid;
  std::vector<TH1D> Muon2_isIsolationValid;
  std::vector<TH1D> Muon1_isTimeValid;
  std::vector<TH1D> Muon2_isTimeValid;
  std::vector<TH1D> Muon1_emEt03;
  std::vector<TH1D> Muon2_emEt03;
  std::vector<TH1D> Muon1_emVetoEt03;
  std::vector<TH1D> Muon2_emVetoEt03;
  std::vector<TH1D> Muon1_hadEt03;
  std::vector<TH1D> Muon2_hadEt03;
  std::vector<TH1D> Muon1_hadVetoEt03;
  std::vector<TH1D> Muon2_hadVetoEt03;
  std::vector<TH1D> Muon1_nJets03;
  std::vector<TH1D> Muon2_nJets03;
  std::vector<TH1D> Muon1_nTracks03;
  std::vector<TH1D> Muon2_nTracks03;
  std::vector<TH1D> Muon1_sumPt03;
  std::vector<TH1D> Muon2_sumPt03;
  std::vector<TH1D> Muon1_trackerVetoPt03;
  std::vector<TH1D> Muon2_trackerVetoPt03;
  std::vector<TH1D> Muon1_emEt05;
  std::vector<TH1D> Muon2_emEt05;
  std::vector<TH1D> Muon1_emVetoEt05;
  std::vector<TH1D> Muon2_emVetoEt05;
  std::vector<TH1D> Muon1_hadEt05;
  std::vector<TH1D> Muon2_hadEt05;
  std::vector<TH1D> Muon1_hadVetoEt05;
  std::vector<TH1D> Muon2_hadVetoEt05;
  std::vector<TH1D> Muon1_nJets05;
  std::vector<TH1D> Muon2_nJets05;
  std::vector<TH1D> Muon1_nTracks05;
  std::vector<TH1D> Muon2_nTracks05;
  std::vector<TH1D> Muon1_sumPt05;
  std::vector<TH1D> Muon2_sumPt05;
  std::vector<TH1D> Muon1_trackerVetoPt05;
  std::vector<TH1D> Muon2_trackerVetoPt05;
  std::vector<TH1D> Muon1_sumChargedHadronPt03;
  std::vector<TH1D> Muon2_sumChargedHadronPt03;
  std::vector<TH1D> Muon1_sumChargedParticlePt03;
  std::vector<TH1D> Muon2_sumChargedParticlePt03;
  std::vector<TH1D> Muon1_sumNeutralHadronEt03;
  std::vector<TH1D> Muon2_sumNeutralHadronEt03;
  std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Muon1_sumPhotonEt03;
  std::vector<TH1D> Muon2_sumPhotonEt03;
  std::vector<TH1D> Muon1_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Muon2_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Muon1_sumPUPt03;
  std::vector<TH1D> Muon2_sumPUPt03;
  std::vector<TH1D> Muon1_sumChargedHadronPt04;
  std::vector<TH1D> Muon2_sumChargedHadronPt04;
  std::vector<TH1D> Muon1_sumChargedParticlePt04;
  std::vector<TH1D> Muon2_sumChargedParticlePt04;
  std::vector<TH1D> Muon1_sumNeutralHadronEt04;
  std::vector<TH1D> Muon2_sumNeutralHadronEt04;
  std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Muon1_sumPhotonEt04;
  std::vector<TH1D> Muon2_sumPhotonEt04;
  std::vector<TH1D> Muon1_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Muon2_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Muon1_sumPUPt04;
  std::vector<TH1D> Muon2_sumPUPt04;
  std::vector<TH1D> Muon1_TriggerMatchdR;
  std::vector<TH1D> Muon2_TriggerMatchdR;
  std::vector<TH1D> Track_TriggerMatchdR;

};
#endif
