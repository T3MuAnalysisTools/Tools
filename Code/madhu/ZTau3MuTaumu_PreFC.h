#ifndef ZTau3MuTaumu_PreFC_h
#define ZTau3MuTaumu_PreFC_h

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


class ZTau3MuTaumu_PreFC : public Selection {

 public:
  ZTau3MuTaumu_PreFC(TString Name_, TString id_);
  virtual ~ZTau3MuTaumu_PreFC();

  virtual void  Configure();
  virtual void  Finish();
  
  enum cuts {WhetherDecayFound=0,
	     Mu1_Candidate_p,
             Mu1_Candidate_eta,
             Mu2_Candidate_p,
             Mu2_Candidate_eta,
	     Mu3_Candidate_p,
             Mu3_Candidate_eta,
	     Tau_mu_Candidate_p,
             Tau_mu_Candidate_eta,
             Mu1_Candidate_recod,
	     Mu2_Candidate_recod,
	     Mu3_Candidate_recod,
	     Tau_mu_Candidate_recod,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
	     TripletPT,
	     nMuons_PF_GL,
             nMuons_pT,
             nMuons_eta,
             nMuons_dR,
	     OSCharge,
	     Tau3MuIsolation,
	     MuonIsolation,
	     VisMass,
	     NCuts};


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:


  std::vector<TH1D>   Tau3MuRelativeIsolation;
  std::vector<TH1D>   OppositeMuRelativeIsolation;
  std::vector<TH1D>   VisibleDiTauMass;
  std::vector<TH1D>   MTT;
  std::vector<TH1D>   TripletMass;
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;
  std::vector<TH1D>   PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PairMass_OppositeSign_dR13;

  std::vector<TH1D>   TripletPt;
  std::vector<TH1D>   OppositeMuonPt;

  std::vector<TH1D>   TripletEta;
  std::vector<TH1D>   OppositeMuonEta;
  
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
  
  std::vector<TH1D>   Selection_Cut_Mu1_P;
  std::vector<TH1D>   Selection_Cut_Mu1_Eta;
  std::vector<TH1D>   Selection_Cut_Mu1_dR;
  std::vector<TH1D>   Selection_Cut_Mu1_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_Mu2_P;
  std::vector<TH1D>   Selection_Cut_Mu2_Eta;
  std::vector<TH1D>   Selection_Cut_Mu2_dR;
  std::vector<TH1D>   Selection_Cut_Mu2_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_Mu3_P;
  std::vector<TH1D>   Selection_Cut_Mu3_Eta;
  std::vector<TH1D>   Selection_Cut_Mu3_dR;
  std::vector<TH1D>   Selection_Cut_Mu3_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after_reco;
  std::vector<TH1D>   Selection_Cut_mu_P;
  std::vector<TH1D>   Selection_Cut_mu_Eta;
  std::vector<TH1D>   Selection_Cut_mu_dR;
  std::vector<TH1D>   Selection_Cut_mu_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_before;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_after;
  std::vector<TH2D>   Selection_Cut_mu_p_eta_after_reco;
  
  std::vector<TH1D>   Selection_Cut_RecoMu_P;
  std::vector<TH1D>   Selection_Cut_RecoMu_Eta;
  
  TRandom rndm;
  double random_num;

};
#endif
