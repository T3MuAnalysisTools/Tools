#ifndef ZTau3MuTaumu_h
#define ZTau3MuTaumu_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"
#include "SimpleFits/FitSoftware/interface/GEFObject.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"
#include "SimpleFits/FitSoftware/interface/TPTRObject.h"


class ZTau3MuTaumu : public Selection {

 public:
  ZTau3MuTaumu(TString Name_, TString id_);
  virtual ~ZTau3MuTaumu();

  virtual void  Configure();
  virtual void  Finish();
  
  enum cuts {PassedFiducialCuts=0,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
	     TripletPT,
	     nMuons_PF_GL,
             nMuons_dR,
             nMuons_pT,
             nMuons_eta,
	     OSCharge,
	     Tau3MuIsolation,
	     MuonIsolation,
	     VisMass,
	     NCuts};


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  TString AnalysisName;
  
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;
  std::vector<TH1D>   PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PairMass_OppositeSign_dR13;
  
  std::vector<TH1D>   dR_betweenTruth_NeutrinoGuess;
  std::vector<TH1D>   dR_betweenTruth_Tau;
  std::vector<TH1D>   Z_Pt;
  std::vector<TH2D>   OS_vs_3mu_trigger;
  
  std::vector<TH1D>   Selection_Cut_3mu_Pt;
  std::vector<TH1D>   Selection_Cut_3mu_Rel_Iso;
  std::vector<TH1D>   Selection_Cut_muon_Pt;
  std::vector<TH1D>   Selection_Cut_muon_Eta;
  std::vector<TH1D>   Selection_Cut_muon_DeltaR_3mu;
  std::vector<TH1D>   Selection_Cut_muon_Rel_Iso;
  std::vector<TH1D>   Selection_Cut_Vis_InvM;
  
  std::vector<TH1D>   Selection_Cut_Mu1_dR;
  std::vector<TH1D>   Selection_Cut_Mu1_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_Mu2_dR;
  std::vector<TH1D>   Selection_Cut_Mu2_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_Mu3_dR;
  std::vector<TH1D>   Selection_Cut_Mu3_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_mu_dR;
  std::vector<TH1D>   Selection_Cut_mu_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_before;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_after;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_after_reco;
  
  std::vector<TH1D>   Selection_Cut_RecoMu_P;
  std::vector<TH1D>   Selection_Cut_RecoMu_Eta;
  
  std::vector<TH1D>   PostSelection_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostSelection_OppositeMuRelativeIsolation;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass;
  std::vector<TH1D>   PostSelection_MTT;
  std::vector<TH1D>   PostSelection_TripletMass;
  
  std::vector<TH1D>   PostSelection_TripletPt;
  std::vector<TH1D>   PostSelection_OppositeMuonPt;
  std::vector<TH1D>   PostSelection_TripletEta;
  std::vector<TH1D>   PostSelection_OppositeMuonEta;
  
  std::vector<TH1D>   PostSelection_MET_Et;
  std::vector<TH1D>   PostSelection_MET_Phi;
  std::vector<TH2D>   PostSelection_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostSelection_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostSelection_Mu1_Pt;
  std::vector<TH1D>   PostSelection_Mu1_Eta;
  std::vector<TH1D>   PostSelection_Mu2_Pt;
  std::vector<TH1D>   PostSelection_Mu2_Eta;
  std::vector<TH1D>   PostSelection_Mu3_Pt;
  std::vector<TH1D>   PostSelection_Mu3_Eta;
  std::vector<TH1D>   PostSelection_mu_Pt;
  std::vector<TH1D>   PostSelection_mu_Eta;
  
  std::vector<TH1D>   PostSelection_FLSignificance;
  std::vector<TH1D>   PostSelection_SVPVTauDirAngle;
  std::vector<TH1D>   PostSelection_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostSelection_VertexChi2KF;
  std::vector<TH1D>   PostSelection_MinDistToIsoTrack;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostSelection_prod_size;
  
  std::vector<TH1D>   PostSelection_Vertex_Dist;
  std::vector<TH1D>   PostSelection_Vertex_Chi2;
  
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_peakspectra1;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_peakspectra2;
  
  Double_t m3m;
  Double_t dataMCtype;
  Double_t event_weight;
  Double_t m12;
  Double_t m13;
  
  Double_t var_TripletPT;
  Double_t var_TripletEta;
  Double_t var_Tau3MuIsolation;
  Double_t var_mu1_pT;
  Double_t var_mu2_pT;
  Double_t var_mu3_pT;
  
  Double_t var_MuonIsolation;
  Double_t var_Muon_pT;
  
  Double_t var_FLSignificance;
  Double_t var_SVPVTauDirAngle;
  Double_t var_ThreeMuVertexChi2KF;
  Double_t var_MinDistToIsoTrack;
  Double_t var_DeltaPhi;
  
  Double_t var_MET_Et;
  Double_t var_MET_Phi;
  
  Double_t var_VisMass;
  Double_t var_DiTauMass_Collinear;
  Double_t var_4Mu_Chi2;
  
  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;
  
  TRandom rndm;
  double random_num;
  
  std::vector<std::vector<TH1D>*> InputFeatureCollection;

};
#endif
