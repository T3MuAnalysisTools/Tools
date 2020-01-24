#ifndef DsToPhiPi_h
#define DsToPhiPi_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class DsToPhiPi : public Selection {

 public:
  DsToPhiPi(TString Name_, TString id_);
  virtual ~DsToPhiPi();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1TOk=0,HLTOk,is2MuTrk,GlobalMu,MuCharge,Mass2Mu,Mu1dR,Mu2dR,TrkdR,Mu1pt,Mu2pt,Trkpt,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

  TFile * file;
  TTree * Sync_tree;

 private:
 // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;

  std::vector<TH1D> DimuondR;
  std::vector<TH1D> Muon1TrkdR;
  std::vector<TH1D> Muon2TrkdR;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> PhiPlusTrackMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D> DsMass;
  std::vector<TH1D> Category;
  std::vector<TH1D> DsGenMatch;

  std::vector<TH1D> DsPt_peak;
  std::vector<TH1D> DsPt_sideband;
  std::vector<TH1D> Ds_Pt;
  std::vector<TH1D> DsP_peak;
  std::vector<TH1D> DsP_sideband;
  std::vector<TH1D> Ds_P;
  std::vector<TH1D> DsM_peak;
  std::vector<TH1D> DsM_sideband;
  std::vector<TH1D> Ds_M;
  std::vector<TH1D> DsL_peak;
  std::vector<TH1D> DsL_sideband;
  std::vector<TH1D> Ds_L;
  std::vector<TH1D> DsEta_peak;
  std::vector<TH1D> DsEta_sideband;
  std::vector<TH1D> Ds_eta;

  //Muon variables
  
  //Track candidate variables
  std::vector<TH1D> Track_P;
  std::vector<TH1D> Track_E;
  std::vector<TH1D> Track_Pt;
  std::vector<TH1D> Track_Eta;
  std::vector<TH1D> Track_Phi;
  std::vector<TH1D> Track_normalizedChi2;
  std::vector<TH1D> Track_numberOfValidHits;
  std::vector<TH1D> Track_charge;

  std::vector<TH1D> Muon1_Pt;
  std::vector<TH1D> Muon1_E;
  std::vector<TH1D> Muon1_P;
  std::vector<TH1D> Muon1_Eta;
  std::vector<TH1D> Muon1_Phi;

  std::vector<TH1D> Muon2_Pt;
  std::vector<TH1D> Muon2_E;
  std::vector<TH1D> Muon2_P;
  std::vector<TH1D> Muon2_Eta;
  std::vector<TH1D> Muon2_Phi;
		
  std::vector<TH1D> Muon1_isGlobal;
  std::vector<TH1D> Muon2_isGlobal;
  std::vector<TH1D> Muon1_isStandAlone;
  std::vector<TH1D> Muon2_isStandAlone;
  std::vector<TH1D> Muon1_isTracker;
  std::vector<TH1D> Muon2_isTracker;
  std::vector<TH1D> Muon1_isCalo;
  std::vector<TH1D> Muon2_isCalo;
  std::vector<TH1D> Muon1_TriggerMatchdR;
  std::vector<TH1D> Muon2_TriggerMatchdR;
  std::vector<TH1D> Track_TriggerMatchdR;

  std::vector<TH1D> DecayLength_peak;
  std::vector<TH1D> DecayLength_sideband;
  std::vector<TH1D> DecayLength_prompt;
  std::vector<TH1D> DecayLength_non_prompt;
  std::vector<TH1D> DecayLength;

  std::vector<TH1D> Muon1_Pt_peak;
  std::vector<TH1D> Muon1_Pt_sideband;
  std::vector<TH1D> Muon1_Eta_peak;
  std::vector<TH1D> Muon1_Eta_sideband;
  std::vector<TH1D> Muon1_Phi_peak;
  std::vector<TH1D> Muon1_Phi_sideband;
  std::vector<TH1D> control_Muon1_Pt;
  std::vector<TH1D> control_Muon1_Eta;
  std::vector<TH1D> control_Muon1_Phi;

  std::vector<TH1D> Muon2_Pt_peak;
  std::vector<TH1D> Muon2_Pt_sideband;
  std::vector<TH1D> Muon2_Eta_peak;
  std::vector<TH1D> Muon2_Eta_sideband;
  std::vector<TH1D> Muon2_Phi_peak;
  std::vector<TH1D> Muon2_Phi_sideband;
  std::vector<TH1D> control_Muon2_Pt;
  std::vector<TH1D> control_Muon2_Eta;
  std::vector<TH1D> control_Muon2_Phi;

  std::vector<TH1D> Track_Pt_peak;
  std::vector<TH1D> Track_Pt_sideband;
  std::vector<TH1D> Track_Eta_peak;
  std::vector<TH1D> Track_Eta_sideband;
  std::vector<TH1D> Track_Phi_peak;
  std::vector<TH1D> Track_Phi_sideband;
  std::vector<TH1D> control_Track_Pt;
  std::vector<TH1D> control_Track_Eta;
  std::vector<TH1D> control_Track_Phi;

  std::vector<TH1D> VertexKFChi2_peak;
  std::vector<TH1D> VertexKFChi2_sideband;
  std::vector<TH1D> VertexKFChi2;

  std::vector<TH1D> SVPVDsDirAngle_peak;
  std::vector<TH1D> SVPVDsDirAngle_sideband;
  std::vector<TH1D> SVPVDsDirAngle;

  std::vector<TH1D> NtracksClose_peak;
  std::vector<TH1D> NtracksClose_sideband;
  std::vector<TH1D> NtracksClose;
  
  std::vector<TH1D> NSV_peak;
  std::vector<TH1D> NSV_sideband;
  std::vector<TH1D> NSV;

  std::vector<TH1D> MinMuon_chi2LocalPosition_peak;
  std::vector<TH1D> MinMuon_chi2LocalPosition_sideband;
  std::vector<TH1D> MinMuon_chi2LocalPosition;
  
  std::vector<TH1D> MindcaTrackSV_peak;
  std::vector<TH1D> MindcaTrackSV_sideband;
  std::vector<TH1D> MindcaTrackSV;

  std::vector<TH1D> MinDca_peak;
  std::vector<TH1D> MinDca_sideband;
  std::vector<TH1D> MinDca;

  std::vector<TH1D> MinD0SigSV_peak;
  std::vector<TH1D> MinD0SigSV_sideband;
  std::vector<TH1D> MinD0SigSV;

  std::vector<TH1D> MinD0SigPV_peak;
  std::vector<TH1D> MinD0SigPV_sideband;
  std::vector<TH1D> MinD0SigPV;

  std::vector<TH1D> MaxVertexPairQuality_peak;
  std::vector<TH1D> MaxVertexPairQuality_sideband;
  std::vector<TH1D> MaxVertexPairQuality;

  std::vector<TH1D> MaxdeltaMuZ_peak;
  std::vector<TH1D> MaxdeltaMuZ_sideband;
  std::vector<TH1D> MaxdeltaMuZ;

  std::vector<TH1D> MaxDca_peak;
  std::vector<TH1D> MaxDca_sideband;
  std::vector<TH1D> MaxDca;

  std::vector<TH1D> MaxD0SigSV_peak;
  std::vector<TH1D> MaxD0SigSV_sideband;
  std::vector<TH1D> MaxD0SigSV;

  std::vector<TH1D> Iso1_peak;
  std::vector<TH1D> Iso1_sideband;
  std::vector<TH1D> Iso1;

  std::vector<TH1D> FLSignificance_peak;
  std::vector<TH1D> FLSignificance_sideband;
  std::vector<TH1D> FLSignificance;


};
#endif
