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

  enum cuts {HLTOk=0,is2MuTrk,GlobalMu,Chi2Cut,MuCharge,Mass2Mu,Mu1dR,Mu2dR,TrkdR,Mu1pt,Mu2pt,Trkpt,NCuts}; 


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

  //Sync variables
  double sync_pt_1;
  double sync_pt_2;
  double sync_pt_3;

  double sync_eta_1;
  double sync_eta_2;
  double sync_eta_3;

  double muon_1_isGlob;
  double muon_2_isGlob;

  double muon_1_isTrack;
  double muon_2_isTrack;

  double phi_mass;
  double ds_mass;

  double evt,run,lumi;

  double sync_DsPhiPiVtx_x;
  double sync_DsPhiPiVtx_y;
  double sync_DsPhiPiVtx_z;
  double sync_DsPhiPiVtx_Chi2;

};
#endif
