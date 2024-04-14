#ifndef ZTau3MuTaue_PreFC_h
#define ZTau3MuTaue_PreFC_h

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


class ZTau3MuTaue_PreFC : public Selection {

 public:
  ZTau3MuTaue_PreFC(TString Name_, TString id_);
  virtual ~ZTau3MuTaue_PreFC();

  virtual void  Configure();
  virtual void  Finish();
  
  /*
  enum cuts {WhetherDecayFound=0,
	     Mu1_Candidate_p,
             Mu1_Candidate_eta,
             Mu2_Candidate_p,
             Mu2_Candidate_eta,
	     Mu3_Candidate_p,
             Mu3_Candidate_eta,
	     Tau_e_Candidate_p,
             Tau_e_Candidate_eta,
             Mu1_Candidate_recod,
	     Mu2_Candidate_recod,
	     Mu3_Candidate_recod,
	     Tau_e_Candidate_recod,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     NCuts}; 
  */
           
           
  enum cuts {WhetherDecayFound=0,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     Mu1_Candidate_p,
             Mu1_Candidate_eta,
             Mu2_Candidate_p,
             Mu2_Candidate_eta,
	     Mu3_Candidate_p,
             Mu3_Candidate_eta,
	     Tau_e_Candidate_p,
             Tau_e_Candidate_eta,
             Mu1_Candidate_recod,
	     Mu2_Candidate_recod,
	     Mu3_Candidate_recod,
	     Tau_e_Candidate_recod,
	     NCuts}; 

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:


  std::vector<TH1D>   Tau3MuRelativeIsolation;
  std::vector<TH1D>   ElectronSumIsolation;
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
  
  std::vector<TH1D>   dR_betweenTruth_NeutrinoGuess;
  std::vector<TH1D>   dR_betweenTruth_Tau;
  std::vector<TH1D>   Z_Pt;
  std::vector<TH2D>   OS_vs_3mu_trigger;
  
  std::vector<TH1D>   Whether_4object_reconstructed;
  std::vector<TH2D>   Whether_triplet_fiducial_or_reco;
  std::vector<TH2D>   Whether_reco_triplet_trigger_L1_and_HLT;
  std::vector<TH2D>   Whether_Mu1_fiducial_and_reco;
  std::vector<TH2D>   Whether_Mu2_fiducial_and_reco;
  std::vector<TH2D>   Whether_Mu3_fiducial_and_reco;
  std::vector<TH2D>   Whether_Tau_e_fiducial_and_reco;
  
  std::vector<TH1D>   Selection_Cut_3mu_Pt;
  std::vector<TH1D>   Selection_Cut_3mu_Rel_Iso;
  std::vector<TH1D>   Selection_Cut_elect_Pt;
  std::vector<TH1D>   Selection_Cut_elect_Eta;
  std::vector<TH1D>   Selection_Cut_elect_DeltaR_3mu;
  std::vector<TH1D>   Selection_Cut_Vis_InvM;
  
  std::vector<TH1D>   Selection_Cut_Mu1_P;
  std::vector<TH1D>   Selection_Cut_Mu1_Eta;
  std::vector<TH1D>   Selection_Cut_Mu1_dR;
  std::vector<TH1D>   Selection_Cut_Mu1_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after_reco;
  std::vector<TH2D>   Selection_Cut_Mu1_p_eta_after_trigger;
  std::vector<TH1D>   Selection_Cut_Mu2_P;
  std::vector<TH1D>   Selection_Cut_Mu2_Eta;
  std::vector<TH1D>   Selection_Cut_Mu2_dR;
  std::vector<TH1D>   Selection_Cut_Mu2_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after_reco;
  std::vector<TH2D>   Selection_Cut_Mu2_p_eta_after_trigger;
  std::vector<TH1D>   Selection_Cut_Mu3_P;
  std::vector<TH1D>   Selection_Cut_Mu3_Eta;
  std::vector<TH1D>   Selection_Cut_Mu3_dR;
  std::vector<TH1D>   Selection_Cut_Mu3_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_before;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after_reco;
  std::vector<TH2D>   Selection_Cut_Mu3_p_eta_after_trigger;
  std::vector<TH1D>   Selection_Cut_El_Pt;
  std::vector<TH1D>   Selection_Cut_El_Eta;
  std::vector<TH1D>   Selection_Cut_El_dR;
  std::vector<TH1D>   Selection_Cut_El_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_El_pt_eta_before;
  std::vector<TH2D>   Selection_Cut_El_pt_eta_after;
  std::vector<TH2D>   Selection_Cut_El_pt_eta_after_reco;
  std::vector<TH2D>   Selection_Cut_El_pt_eta_after_trigger;
  
  std::vector<TH1D>   Selection_Cut_RecoMu_P;
  std::vector<TH1D>   Selection_Cut_RecoMu_Eta;
  std::vector<TH1D>   Selection_Cut_RecoEl_Pt;
  std::vector<TH1D>   Selection_Cut_RecoEl_Eta;
  
  TRandom rndm;
  double random_num;
  
};
#endif
