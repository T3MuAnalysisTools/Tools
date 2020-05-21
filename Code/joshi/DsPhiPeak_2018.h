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

  enum cuts {TriggerOk,TwoMuTrkCandidate,OSMuons,Mu1PtCut,Mu2PtCut,MuonID,MuMuMassCut,TrackPtCut,NTrackHits, /*ChiSqCut,*/ DsMassCut,TriggerMatchMu1,TriggerMatchMu2,TriggerMatchTrack,GenMatch,NCuts}; 

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  // PU Weights
  TFile* PUWeightFile;
  TH1D* puWeights;

  // Selection Variables
  float PhiMassHigh;
  float PhiMassLow;
  float sidebandDsMin;
  float sidebandDsMax;
  float peakDsMin;
  float peakDsMax;
  float dsMassMin;
  float dsMassMax;

  double nSidebands;
  double nPeak;

  unsigned int Muon_1_idx, Muon_2_idx, Track_idx;

  TRandom rndm;
  Float_t random_num;

  // Build collection of histograms
  std::vector<std::vector<TH1D>*> peakCollection;
  std::vector<std::vector<TH1D>*> sidebandCollection;
  std::vector<std::vector<TH1D>*> validationCollection;

  // Initializhere your analysis histograms

  std::vector<TH1D> TripleMass;
  std::vector<TH2D> PhiMassVsDsMass;

  std::vector<TH1D> ChargeMisId;
  std::vector<TH1D> TriggerEfficiency;
  


  ///////////////////
  // validation block
  ///////////////////
  
  // More plots
  std::vector<TH1D> Muon1_isGlobal_peak;
  std::vector<TH1D> Muon2_isGlobal_peak;
  std::vector<TH1D> Muon1_isStandAlone_peak;
  std::vector<TH1D> Muon2_isStandAlone_peak;
  std::vector<TH1D> Muon1_isTracker_peak;
  std::vector<TH1D> Muon2_isTracker_peak;
  std::vector<TH1D> Muon1_isCalo_peak;
  std::vector<TH1D> Muon2_isCalo_peak;
  std::vector<TH1D> Muon1_isIsolationValid_peak;
  std::vector<TH1D> Muon2_isIsolationValid_peak;
  std::vector<TH1D> Muon1_TriggerMatchdR_peak;
  std::vector<TH1D> Muon2_TriggerMatchdR_peak;
  std::vector<TH1D> Track_TriggerMatchdR_peak;
  std::vector<TH1D> NVtx_peak;
  std::vector<TH1D> NVtx_woPUWeights_peak;
  std::vector<TH1D> PhiMass_peak;

  std::vector<TH1D> Muon1_isGlobal_sideband;
  std::vector<TH1D> Muon2_isGlobal_sideband;
  std::vector<TH1D> Muon1_isStandAlone_sideband;
  std::vector<TH1D> Muon2_isStandAlone_sideband;
  std::vector<TH1D> Muon1_isTracker_sideband;
  std::vector<TH1D> Muon2_isTracker_sideband;
  std::vector<TH1D> Muon1_isCalo_sideband;
  std::vector<TH1D> Muon2_isCalo_sideband;
  std::vector<TH1D> Muon1_isIsolationValid_sideband;
  std::vector<TH1D> Muon2_isIsolationValid_sideband;
  std::vector<TH1D> Muon1_TriggerMatchdR_sideband;
  std::vector<TH1D> Muon2_TriggerMatchdR_sideband;
  std::vector<TH1D> Track_TriggerMatchdR_sideband;
  std::vector<TH1D> NVtx_sideband;
  std::vector<TH1D> NVtx_woPUWeights_sideband;
  std::vector<TH1D> PhiMass_sideband;

  std::vector<TH1D> Muon1_isGlobal_validation;
  std::vector<TH1D> Muon2_isGlobal_validation;
  std::vector<TH1D> Muon1_isStandAlone_validation;
  std::vector<TH1D> Muon2_isStandAlone_validation;
  std::vector<TH1D> Muon1_isTracker_validation;
  std::vector<TH1D> Muon2_isTracker_validation;
  std::vector<TH1D> Muon1_isCalo_validation;
  std::vector<TH1D> Muon2_isCalo_validation;
  std::vector<TH1D> Muon1_isIsolationValid_validation;
  std::vector<TH1D> Muon2_isIsolationValid_validation;
  std::vector<TH1D> Muon1_TriggerMatchdR_validation;
  std::vector<TH1D> Muon2_TriggerMatchdR_validation;
  std::vector<TH1D> Track_TriggerMatchdR_validation;
  std::vector<TH1D> NVtx_validation;
  std::vector<TH1D> NVtx_woPUWeights_validation;
  std::vector<TH1D> PhiMass_validation;

  std::vector<TH1D> DimuondR_peak;
  std::vector<TH1D> Muon1TrkdR_peak;
  std::vector<TH1D> Muon2TrkdR_peak;

  //Track candidate variables
  std::vector<TH1D> Track_Pt_peak;
  std::vector<TH1D> Track_Eta_peak;
  std::vector<TH1D> Track_Phi_peak;
  std::vector<TH1D> Track_vx_peak;
  std::vector<TH1D> Track_vy_peak;
  std::vector<TH1D> Track_vz_peak;
  std::vector<TH1D> Track_normalizedChi2_peak;
  std::vector<TH1D> Track_numberOfValidHits_peak;
  std::vector<TH1D> Track_charge_peak;
  std::vector<TH1D> Track_dxy_peak;
  std::vector<TH1D> Track_dz_peak;

  //Dimuon variables

  std::vector<TH1D> Muon1_Pt_peak;
  std::vector<TH1D> Muon1_Eta_peak;
  std::vector<TH1D> Muon1_Phi_peak;
  std::vector<TH1D> Muon1_vx_peak;
  std::vector<TH1D> Muon1_vy_peak;
  std::vector<TH1D> Muon1_vz_peak;

  std::vector<TH1D> Muon2_Pt_peak;
  std::vector<TH1D> Muon2_Eta_peak;
  std::vector<TH1D> Muon2_Phi_peak;
  std::vector<TH1D> Muon2_vx_peak;
  std::vector<TH1D> Muon2_vy_peak;
  std::vector<TH1D> Muon2_vz_peak;
  
  // sidebands 
  std::vector<TH1D> DimuondR_sideband;
  std::vector<TH1D> Muon1TrkdR_sideband;
  std::vector<TH1D> Muon2TrkdR_sideband;

  //Track candidate variables
  std::vector<TH1D> Track_Pt_sideband;
  std::vector<TH1D> Track_Eta_sideband;
  std::vector<TH1D> Track_Phi_sideband;
  std::vector<TH1D> Track_vx_sideband;
  std::vector<TH1D> Track_vy_sideband;
  std::vector<TH1D> Track_vz_sideband;
  std::vector<TH1D> Track_normalizedChi2_sideband;
  std::vector<TH1D> Track_numberOfValidHits_sideband;
  std::vector<TH1D> Track_charge_sideband;
  std::vector<TH1D> Track_dxy_sideband;
  std::vector<TH1D> Track_dz_sideband;

  //Dimuon variables

  std::vector<TH1D> Muon1_Pt_sideband;
  std::vector<TH1D> Muon1_Eta_sideband;
  std::vector<TH1D> Muon1_Phi_sideband;
  std::vector<TH1D> Muon1_vx_sideband;
  std::vector<TH1D> Muon1_vy_sideband;
  std::vector<TH1D> Muon1_vz_sideband;

  std::vector<TH1D> Muon2_Pt_sideband;
  std::vector<TH1D> Muon2_Eta_sideband;
  std::vector<TH1D> Muon2_Phi_sideband;
  std::vector<TH1D> Muon2_vx_sideband;
  std::vector<TH1D> Muon2_vy_sideband;
  std::vector<TH1D> Muon2_vz_sideband;
  
  // validation
  std::vector<TH1D> DimuondR_validation;
  std::vector<TH1D> Muon1TrkdR_validation;
  std::vector<TH1D> Muon2TrkdR_validation;

  //Track candidate variables
  std::vector<TH1D> Track_Pt_validation;
  std::vector<TH1D> Track_Eta_validation;
  std::vector<TH1D> Track_Phi_validation;
  std::vector<TH1D> Track_vx_validation;
  std::vector<TH1D> Track_vy_validation;
  std::vector<TH1D> Track_vz_validation;
  std::vector<TH1D> Track_normalizedChi2_validation;
  std::vector<TH1D> Track_numberOfValidHits_validation;
  std::vector<TH1D> Track_charge_validation;
  std::vector<TH1D> Track_dxy_validation;
  std::vector<TH1D> Track_dz_validation;

  //Dimuon variables

  std::vector<TH1D> Muon1_Pt_validation;
  std::vector<TH1D> Muon1_Eta_validation;
  std::vector<TH1D> Muon1_Phi_validation;
  std::vector<TH1D> Muon1_vx_validation;
  std::vector<TH1D> Muon1_vy_validation;
  std::vector<TH1D> Muon1_vz_validation;

  std::vector<TH1D> Muon2_Pt_validation;
  std::vector<TH1D> Muon2_Eta_validation;
  std::vector<TH1D> Muon2_Phi_validation;
  std::vector<TH1D> Muon2_vx_validation;
  std::vector<TH1D> Muon2_vy_validation;
  std::vector<TH1D> Muon2_vz_validation;


  
  std::vector<TH1D> VertexKFChi2_peak;
  std::vector<TH1D> VertexKFChi2_sideband;
  std::vector<TH1D> VertexKFChi2_validation;

  std::vector<TH1D> SVPVDsDirAngle_peak;
  std::vector<TH1D> SVPVDsDirAngle_sideband;
  std::vector<TH1D> SVPVDsDirAngle_validation;

  std::vector<TH1D> NtracksClose_peak;
  std::vector<TH1D> NtracksClose_sideband;
  std::vector<TH1D> NtracksClose_validation;
  
  std::vector<TH1D> NSV_peak;
  std::vector<TH1D> NSV_sideband;
  std::vector<TH1D> NSV_validation;

  std::vector<TH1D> MinMuon_chi2LocalPosition_peak;
  std::vector<TH1D> MinMuon_chi2LocalPosition_sideband;
  std::vector<TH1D> MinMuon_chi2LocalPosition_validation;
  
  std::vector<TH1D> MindcaTrackSV_peak;
  std::vector<TH1D> MindcaTrackSV_sideband;
  std::vector<TH1D> MindcaTrackSV_validation;

  std::vector<TH1D> MinDca_peak;
  std::vector<TH1D> MinDca_sideband;
  std::vector<TH1D> MinDca_validation;

  std::vector<TH1D> MinD0SigSV_peak;
  std::vector<TH1D> MinD0SigSV_sideband;
  std::vector<TH1D> MinD0SigSV_validation;

  std::vector<TH1D> MinD0SigPV_peak;
  std::vector<TH1D> MinD0SigPV_sideband;
  std::vector<TH1D> MinD0SigPV_validation;

  std::vector<TH1D> MaxVertexPairQuality_peak;
  std::vector<TH1D> MaxVertexPairQuality_sideband;
  std::vector<TH1D> MaxVertexPairQuality_validation;

  std::vector<TH1D> MaxdeltaMuZ_peak;
  std::vector<TH1D> MaxdeltaMuZ_sideband;
  std::vector<TH1D> MaxdeltaMuZ_validation;

  std::vector<TH1D> MaxDca_peak;
  std::vector<TH1D> MaxDca_sideband;
  std::vector<TH1D> MaxDca_validation;

  std::vector<TH1D> MaxD0SigSV_peak;
  std::vector<TH1D> MaxD0SigSV_sideband;
  std::vector<TH1D> MaxD0SigSV_validation;

  std::vector<TH1D> MaxD0SigPV_peak;
  std::vector<TH1D> MaxD0SigPV_sideband;
  std::vector<TH1D> MaxD0SigPV_validation;

  std::vector<TH1D> Iso1_peak;
  std::vector<TH1D> Iso1_sideband;
  std::vector<TH1D> Iso1_validation;

  std::vector<TH1D> FLSignificance_peak;
  std::vector<TH1D> FLSignificance_sideband;
  std::vector<TH1D> FLSignificance_validation;
  
  std::vector<TH1D> DecayLength_peak;
  std::vector<TH1D> DecayLength_prompt_peak;
  std::vector<TH1D> DecayLength_non_prompt_peak;
  
  std::vector<TH1D> DecayLength_sideband;
  std::vector<TH1D> DecayLength_prompt_sideband;
  std::vector<TH1D> DecayLength_non_prompt_sideband;

  std::vector<TH1D> DecayLength_validation;
  std::vector<TH1D> DecayLength_prompt_validation;
  std::vector<TH1D> DecayLength_non_prompt_validation;

  std::vector<TH1D> TripleMass_sideband;
  std::vector<TH1D> TripleMass_peak;
  std::vector<TH1D> TripleMass_validation;

/*
  std::vector<TH1D> Muon1_isTimeValid;
  std::vector<TH1D> Muon1_emEt03;
  std::vector<TH1D> Muon2_emEt03;
  std::vector<TH1D> Muon1_emVetoEt03;
  std::vector<TH1D> Muon1_hadEt03;
  std::vector<TH1D> Muon2_hadEt03;
  std::vector<TH1D> Muon1_hadVetoEt03;
  std::vector<TH1D> Muon1_nJets03;
  std::vector<TH1D> Muon2_nJets03;
  std::vector<TH1D> Muon1_nTracks03;
  std::vector<TH1D> Muon1_sumPt03;
  std::vector<TH1D> Muon2_sumPt03;
  std::vector<TH1D> Muon1_trackerVetoPt03;
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
*/

};
#endif
