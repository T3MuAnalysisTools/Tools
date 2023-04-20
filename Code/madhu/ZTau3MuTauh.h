#ifndef ZTau3MuTauh_h
#define ZTau3MuTauh_h

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

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class ZTau3MuTauh : public Selection {

 public:
  ZTau3MuTauh(TString Name_, TString id_);
  virtual ~ZTau3MuTauh();

  virtual void  Configure();
  virtual void  Finish();
  
  enum cuts {PassedFiducialCuts=0,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
             TripletPT,
             Tau3MuIsolation,
             nTaus_dR,
             nTaus_pT,
             nTaus_eta,
	     OSCharge,
	     DeepTauJets,
	     DeepTauMuons,
	     DeepTauElectrons,
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
  std::vector<TH1D>   Selection_Cut_tauh_Pt;
  std::vector<TH1D>   Selection_Cut_tauh_Eta;
  std::vector<TH1D>   Selection_Cut_tauh_DeltaR_3mu;
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

  std::vector<TH1D>   Selection_Cut_h_dR;
  std::vector<TH1D>   Selection_Cut_h_dR_large_scale;
  std::vector<TH2D>   Selection_Cut_h_pt_eta_before;
  std::vector<TH2D>   Selection_Cut_h_pt_eta_after;
  std::vector<TH2D>   Selection_Cut_h_pt_eta_after_reco;
  
  std::vector<TH1D>   Selection_Cut_RecoMu_P;
  std::vector<TH1D>   Selection_Cut_RecoMu_Eta;
  std::vector<TH1D>   Selection_Cut_RecoH_Pt;
  std::vector<TH1D>   Selection_Cut_RecoH_Eta;
  
  //After Selection
  
  std::vector<TH1D>   PostSelection_NumberOfTaus;
  std::vector<TH1D>   PostSelection_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostSelection_TauHDecayMode;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass;
  std::vector<TH1D>   PostSelection_MTT;
  std::vector<TH1D>   PostSelection_TripletMass;
  
  std::vector<TH1D>   PostSelection_TripletPt;
  std::vector<TH1D>   PostSelection_OppositeTauPt;
  std::vector<TH1D>   PostSelection_TripletEta;
  std::vector<TH1D>   PostSelection_OppositeTauEta;
  
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
  std::vector<TH1D>   PostSelection_h_Pt;
  std::vector<TH1D>   PostSelection_h_Eta;
  
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
  
  std::vector<TH1D>   PostSelection_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostSelection_prod_size;
  
  std::vector<TH1D>   PostSelection_BDT_Output;
  
  std::vector<TH1D>   PostSelection_BDT_Output_MC_Bkg;
  
  //After BDT
  
  std::vector<TH1D>   PostBDT_NumberOfTaus;
  std::vector<TH1D>   PostBDT_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostBDT_TauHDecayMode;
  std::vector<TH1D>   PostBDT_VisibleDiTauMass;
  std::vector<TH1D>   PostBDT_MTT;
  std::vector<TH1D>   PostBDT_TripletMass;
  
  std::vector<TH1D>   PostBDT_TripletPt;
  std::vector<TH1D>   PostBDT_OppositeTauPt;
  std::vector<TH1D>   PostBDT_TripletEta;
  std::vector<TH1D>   PostBDT_OppositeTauEta;
  
  std::vector<TH1D>   PostBDT_MET_Et;
  std::vector<TH1D>   PostBDT_MET_Phi;
  std::vector<TH2D>   PostBDT_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostBDT_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostBDT_Mu1_Pt;
  std::vector<TH1D>   PostBDT_Mu1_Eta;
  std::vector<TH1D>   PostBDT_Mu2_Pt;
  std::vector<TH1D>   PostBDT_Mu2_Eta;
  std::vector<TH1D>   PostBDT_Mu3_Pt;
  std::vector<TH1D>   PostBDT_Mu3_Eta;
  std::vector<TH1D>   PostBDT_h_Pt;
  std::vector<TH1D>   PostBDT_h_Eta;
  
  std::vector<TH1D>   PostBDT_FLSignificance;
  std::vector<TH1D>   PostBDT_SVPVTauDirAngle;
  std::vector<TH1D>   PostBDT_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostBDT_VertexChi2KF;
  std::vector<TH1D>   PostBDT_MinDistToIsoTrack;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostBDT_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostBDT_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostBDT_prod_size;
  
  
  //Pre BDT 2D scan
  
  std::vector<TH2D>   BDT_2Dscan_NumberOfTaus;
  std::vector<TH2D>   BDT_2Dscan_Tau3MuRelativeIsolation;
  std::vector<TH2D>   BDT_2Dscan_TauHDecayMode;
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass;
  std::vector<TH2D>   BDT_2Dscan_MTT;
  std::vector<TH2D>   BDT_2Dscan_TripletMass;
  
  std::vector<TH2D>   BDT_2Dscan_TripletPt;
  std::vector<TH2D>   BDT_2Dscan_OppositeTauPt;
  std::vector<TH2D>   BDT_2Dscan_TripletEta;
  std::vector<TH2D>   BDT_2Dscan_OppositeTauEta;
  
  std::vector<TH2D>   BDT_2Dscan_MET_Et;
  std::vector<TH2D>   BDT_2Dscan_MET_Phi;
  std::vector<TH2D>   BDT_2Dscan_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   BDT_2Dscan_MET_vs_NeutrinoPt;
  
  std::vector<TH2D>   BDT_2Dscan_Mu1_Pt;
  std::vector<TH2D>   BDT_2Dscan_Mu1_Eta;
  std::vector<TH2D>   BDT_2Dscan_Mu2_Pt;
  std::vector<TH2D>   BDT_2Dscan_Mu2_Eta;
  std::vector<TH2D>   BDT_2Dscan_Mu3_Pt;
  std::vector<TH2D>   BDT_2Dscan_Mu3_Eta;
  std::vector<TH2D>   BDT_2Dscan_h_Pt;
  std::vector<TH2D>   BDT_2Dscan_h_Eta;
  
  std::vector<TH2D>   BDT_2Dscan_FLSignificance;
  std::vector<TH2D>   BDT_2Dscan_SVPVTauDirAngle;
  std::vector<TH2D>   BDT_2Dscan_SVPVTauDirAngle_largescale;
  std::vector<TH2D>   BDT_2Dscan_VertexChi2KF;
  std::vector<TH2D>   BDT_2Dscan_MinDistToIsoTrack;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_cos;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_pT;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_MET;
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass_Collinear;
  
  std::vector<TH2D>   BDT_2Dscan_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH2D>   BDT_2Dscan_prod_size;
  
  
  Float_t m3m;
  Float_t dataMCtype;
  Float_t event_weight;
  Float_t m12;
  Float_t m13;
  
  Float_t var_TripletPT;
  Float_t var_TripletEta;
  Float_t var_Tau3MuIsolation;
  Float_t var_mu1_pT;
  Float_t var_mu2_pT;
  Float_t var_mu3_pT;
  Float_t var_mu1_eta;
  Float_t var_mu2_eta;
  Float_t var_mu3_eta;
  
  Float_t var_Tau_pT;
  Float_t var_Tau_eta;
  
  Float_t var_FLSignificance;
  Float_t var_SVPVTauDirAngle;
  Float_t var_ThreeMuVertexChi2KF;
  Float_t var_MinDistToIsoTrack;
  Float_t var_DeltaPhi;
  Float_t var_Phi_To_Opposite_Side;
  
  
  Float_t var_MET_Et;
  Float_t var_MET_Phi;
  
  Float_t var_VisMass;
  Float_t var_DiTauMass_Collinear;
  
  Float_t BDT_Evaluated;
  
  Float_t BDT_Evaluated_MC_Bkg;
  
  //for combine
  Float_t combine_var_m3m;
  Float_t combine_var_ID;
  Float_t combine_var_BDTOutput;
  
  TMVA::Reader *reader_Tauh;
  TMVA::Reader *reader_Tauh_MC_Bkg;
  
  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;
  
  TTree *T3MCombineTree;
  TFile *T3MFCombineTree;
  
  TRandom rndm;
  double random_num;
  
  std::vector<std::vector<TH1D>*> InputFeatureCollection;
  std::vector<std::vector<TH2D>*> InputFeatureCollection_2D;

};
#endif
