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

#include "SimpleFits/FitSoftware/interface/Track_Functor.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "EventClassifier.h"


class ZTau3MuTauh : public Selection {

 public:
  ZTau3MuTauh(TString Name_, TString id_);
  virtual ~ZTau3MuTauh();

  virtual void  Configure();
  virtual void  Finish();
  
  enum cuts {WhetherZTTDecayFound=0,
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

  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double tauMassResCutLow, tauMassResCutHigh;
  double phiVetoSigma, omegaVetoSigma;
  
  float bdt_cut_1_;  float bdt_cut_2_;
  
  TString AnalysisName;
  
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;


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
  
  //Cat A
  std::vector<TH1D>   PostSelection_A_NumberOfTaus;
  std::vector<TH1D>   PostSelection_A_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostSelection_A_TauHDecayMode;
  std::vector<TH1D>   PostSelection_A_VisibleDiTauMass;
  std::vector<TH1D>   PostSelection_A_MTT;
  std::vector<TH1D>   PostSelection_A_TripletMass;
  
  std::vector<TH1D>   PostSelection_A_TripletPt;
  std::vector<TH1D>   PostSelection_A_OppositeTauPt;
  std::vector<TH1D>   PostSelection_A_TripletEta;
  std::vector<TH1D>   PostSelection_A_OppositeTauEta;
  
  std::vector<TH1D>   PostSelection_A_MET_Et;
  std::vector<TH1D>   PostSelection_A_MET_Phi;
  std::vector<TH2D>   PostSelection_A_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostSelection_A_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostSelection_A_Mu1_Pt;
  std::vector<TH1D>   PostSelection_A_Mu1_Eta;
  std::vector<TH1D>   PostSelection_A_Mu2_Pt;
  std::vector<TH1D>   PostSelection_A_Mu2_Eta;
  std::vector<TH1D>   PostSelection_A_Mu3_Pt;
  std::vector<TH1D>   PostSelection_A_Mu3_Eta;
  std::vector<TH1D>   PostSelection_A_h_Pt;
  std::vector<TH1D>   PostSelection_A_h_Eta;
  
  std::vector<TH1D>   PostSelection_A_FLSignificance;
  std::vector<TH1D>   PostSelection_A_SVPVTauDirAngle;
  std::vector<TH1D>   PostSelection_A_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostSelection_A_VertexChi2KF;
  std::vector<TH1D>   PostSelection_A_MinDistToIsoTrack;
  std::vector<TH1D>   PostSelection_A_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostSelection_A_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostSelection_A_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostSelection_A_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostSelection_A_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostSelection_A_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostSelection_A_prod_size;
  
  std::vector<TH1D>   PostSelection_A_BDT_Output;
  std::vector<TH1D>   PostSelection_A_BDT_Output_MC_Bkg;
  std::vector<TH2D>   PostSelection_A_BDT_Output_Data_vs_MC_Bkg;
  
  std::vector<TH1D>   PostSelection_A_HPS_Vtx_DiTauMass;
  std::vector<TH1D>   PostSelection_A_HPS_FL_Sigfig;
  std::vector<TH1D>   PostSelection_A_HPS_Vtx_Tau_to_Tauh;
  
  
  std::vector<TH1D>   PostSelection_A_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostSelection_A_PairMass_OppositeSign_dR13;
  
  //Cat B
  std::vector<TH1D>   PostSelection_B_NumberOfTaus;
  std::vector<TH1D>   PostSelection_B_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostSelection_B_TauHDecayMode;
  std::vector<TH1D>   PostSelection_B_VisibleDiTauMass;
  std::vector<TH1D>   PostSelection_B_MTT;
  std::vector<TH1D>   PostSelection_B_TripletMass;
  
  std::vector<TH1D>   PostSelection_B_TripletPt;
  std::vector<TH1D>   PostSelection_B_OppositeTauPt;
  std::vector<TH1D>   PostSelection_B_TripletEta;
  std::vector<TH1D>   PostSelection_B_OppositeTauEta;
  
  std::vector<TH1D>   PostSelection_B_MET_Et;
  std::vector<TH1D>   PostSelection_B_MET_Phi;
  std::vector<TH2D>   PostSelection_B_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostSelection_B_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostSelection_B_Mu1_Pt;
  std::vector<TH1D>   PostSelection_B_Mu1_Eta;
  std::vector<TH1D>   PostSelection_B_Mu2_Pt;
  std::vector<TH1D>   PostSelection_B_Mu2_Eta;
  std::vector<TH1D>   PostSelection_B_Mu3_Pt;
  std::vector<TH1D>   PostSelection_B_Mu3_Eta;
  std::vector<TH1D>   PostSelection_B_h_Pt;
  std::vector<TH1D>   PostSelection_B_h_Eta;
  
  std::vector<TH1D>   PostSelection_B_FLSignificance;
  std::vector<TH1D>   PostSelection_B_SVPVTauDirAngle;
  std::vector<TH1D>   PostSelection_B_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostSelection_B_VertexChi2KF;
  std::vector<TH1D>   PostSelection_B_MinDistToIsoTrack;
  std::vector<TH1D>   PostSelection_B_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostSelection_B_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostSelection_B_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostSelection_B_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostSelection_B_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostSelection_B_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostSelection_B_prod_size;
  
  std::vector<TH1D>   PostSelection_B_BDT_Output;
  std::vector<TH1D>   PostSelection_B_BDT_Output_MC_Bkg;
  std::vector<TH2D>   PostSelection_B_BDT_Output_Data_vs_MC_Bkg;
  
  std::vector<TH1D>   PostSelection_B_HPS_Vtx_DiTauMass;
  std::vector<TH1D>   PostSelection_B_HPS_FL_Sigfig;
  std::vector<TH1D>   PostSelection_B_HPS_Vtx_Tau_to_Tauh;
  
  
  std::vector<TH1D>   PostSelection_B_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostSelection_B_PairMass_OppositeSign_dR13;


  
  //After BDT
  
  //Cat A
  std::vector<TH1D>   PostBDT_A_TripletMass_VeryLooseCut;
  
  std::vector<TH1D>   PostBDT_A_NumberOfTaus;
  std::vector<TH1D>   PostBDT_A_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostBDT_A_TauHDecayMode;
  std::vector<TH1D>   PostBDT_A_VisibleDiTauMass;
  std::vector<TH1D>   PostBDT_A_MTT;
  std::vector<TH1D>   PostBDT_A_TripletMass;
  std::vector<TH1D>   PostBDT_A_TripletMass_a;
  std::vector<TH1D>   PostBDT_A_TripletMass_b;
  
  std::vector<TH1D>   PostBDT_A_TripletPt;
  std::vector<TH1D>   PostBDT_A_OppositeTauPt;
  std::vector<TH1D>   PostBDT_A_TripletEta;
  std::vector<TH1D>   PostBDT_A_OppositeTauEta;
  
  std::vector<TH1D>   PostBDT_A_MET_Et;
  std::vector<TH1D>   PostBDT_A_MET_Phi;
  std::vector<TH2D>   PostBDT_A_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostBDT_A_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostBDT_A_Mu1_Pt;
  std::vector<TH1D>   PostBDT_A_Mu1_Eta;
  std::vector<TH1D>   PostBDT_A_Mu2_Pt;
  std::vector<TH1D>   PostBDT_A_Mu2_Eta;
  std::vector<TH1D>   PostBDT_A_Mu3_Pt;
  std::vector<TH1D>   PostBDT_A_Mu3_Eta;
  std::vector<TH1D>   PostBDT_A_h_Pt;
  std::vector<TH1D>   PostBDT_A_h_Eta;
  
  std::vector<TH1D>   PostBDT_A_FLSignificance;
  std::vector<TH1D>   PostBDT_A_SVPVTauDirAngle;
  std::vector<TH1D>   PostBDT_A_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostBDT_A_VertexChi2KF;
  std::vector<TH1D>   PostBDT_A_MinDistToIsoTrack;
  std::vector<TH1D>   PostBDT_A_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostBDT_A_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostBDT_A_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostBDT_A_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostBDT_A_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostBDT_A_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostBDT_A_prod_size;
  
  std::vector<TH1D>   PostBDT_A_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostBDT_A_PairMass_OppositeSign_dR13;
  
  //Cat B
  std::vector<TH1D>   PostBDT_B_TripletMass_VeryLooseCut;
  
  std::vector<TH1D>   PostBDT_B_NumberOfTaus;
  std::vector<TH1D>   PostBDT_B_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostBDT_B_TauHDecayMode;
  std::vector<TH1D>   PostBDT_B_VisibleDiTauMass;
  std::vector<TH1D>   PostBDT_B_MTT;
  std::vector<TH1D>   PostBDT_B_TripletMass;
  std::vector<TH1D>   PostBDT_B_TripletMass_a;
  std::vector<TH1D>   PostBDT_B_TripletMass_b;
  
  std::vector<TH1D>   PostBDT_B_TripletPt;
  std::vector<TH1D>   PostBDT_B_OppositeTauPt;
  std::vector<TH1D>   PostBDT_B_TripletEta;
  std::vector<TH1D>   PostBDT_B_OppositeTauEta;
  
  std::vector<TH1D>   PostBDT_B_MET_Et;
  std::vector<TH1D>   PostBDT_B_MET_Phi;
  std::vector<TH2D>   PostBDT_B_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   PostBDT_B_MET_vs_NeutrinoPt;
  
  std::vector<TH1D>   PostBDT_B_Mu1_Pt;
  std::vector<TH1D>   PostBDT_B_Mu1_Eta;
  std::vector<TH1D>   PostBDT_B_Mu2_Pt;
  std::vector<TH1D>   PostBDT_B_Mu2_Eta;
  std::vector<TH1D>   PostBDT_B_Mu3_Pt;
  std::vector<TH1D>   PostBDT_B_Mu3_Eta;
  std::vector<TH1D>   PostBDT_B_h_Pt;
  std::vector<TH1D>   PostBDT_B_h_Eta;
  
  std::vector<TH1D>   PostBDT_B_FLSignificance;
  std::vector<TH1D>   PostBDT_B_SVPVTauDirAngle;
  std::vector<TH1D>   PostBDT_B_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostBDT_B_VertexChi2KF;
  std::vector<TH1D>   PostBDT_B_MinDistToIsoTrack;
  std::vector<TH1D>   PostBDT_B_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostBDT_B_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostBDT_B_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostBDT_B_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostBDT_B_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostBDT_B_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostBDT_B_prod_size;
  
  std::vector<TH1D>   PostBDT_B_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostBDT_B_PairMass_OppositeSign_dR13;

  
  //Pre BDT 2D scan
  
  
  //Cat A
  std::vector<TH2D>   BDT_2Dscan_A_NumberOfTaus;
  std::vector<TH2D>   BDT_2Dscan_A_Tau3MuRelativeIsolation;
  std::vector<TH2D>   BDT_2Dscan_A_TauHDecayMode;
  std::vector<TH2D>   BDT_2Dscan_A_VisibleDiTauMass;
  std::vector<TH2D>   BDT_2Dscan_A_MTT;
  std::vector<TH2D>   BDT_2Dscan_A_TripletMass;
  
  std::vector<TH2D>   BDT_2Dscan_A_TripletPt;
  std::vector<TH2D>   BDT_2Dscan_A_OppositeTauPt;
  std::vector<TH2D>   BDT_2Dscan_A_TripletEta;
  std::vector<TH2D>   BDT_2Dscan_A_OppositeTauEta;
  
  std::vector<TH2D>   BDT_2Dscan_A_MET_Et;
  std::vector<TH2D>   BDT_2Dscan_A_MET_Phi;
  std::vector<TH2D>   BDT_2Dscan_A_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   BDT_2Dscan_A_MET_vs_NeutrinoPt;
  
  std::vector<TH2D>   BDT_2Dscan_A_Mu1_Pt;
  std::vector<TH2D>   BDT_2Dscan_A_Mu1_Eta;
  std::vector<TH2D>   BDT_2Dscan_A_Mu2_Pt;
  std::vector<TH2D>   BDT_2Dscan_A_Mu2_Eta;
  std::vector<TH2D>   BDT_2Dscan_A_Mu3_Pt;
  std::vector<TH2D>   BDT_2Dscan_A_Mu3_Eta;
  std::vector<TH2D>   BDT_2Dscan_A_h_Pt;
  std::vector<TH2D>   BDT_2Dscan_A_h_Eta;
  
  std::vector<TH2D>   BDT_2Dscan_A_FLSignificance;
  std::vector<TH2D>   BDT_2Dscan_A_SVPVTauDirAngle;
  std::vector<TH2D>   BDT_2Dscan_A_SVPVTauDirAngle_largescale;
  std::vector<TH2D>   BDT_2Dscan_A_VertexChi2KF;
  std::vector<TH2D>   BDT_2Dscan_A_MinDistToIsoTrack;
  std::vector<TH2D>   BDT_2Dscan_A_Kinematics_MissingTrMass;
  std::vector<TH2D>   BDT_2Dscan_A_Kinematics_MissingTrMass_cos;
  std::vector<TH2D>   BDT_2Dscan_A_Kinematics_MissingTrMass_pT;
  std::vector<TH2D>   BDT_2Dscan_A_Kinematics_MissingTrMass_MET;
  std::vector<TH2D>   BDT_2Dscan_A_VisibleDiTauMass_Collinear;
  
  std::vector<TH2D>   BDT_2Dscan_A_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH2D>   BDT_2Dscan_A_prod_size;
  
  std::vector<TH2D>   BDT_2Dscan_A_PairMass_OppositeSign_dR12;
  std::vector<TH2D>   BDT_2Dscan_A_PairMass_OppositeSign_dR13;
  
  //Cat B
  std::vector<TH2D>   BDT_2Dscan_B_NumberOfTaus;
  std::vector<TH2D>   BDT_2Dscan_B_Tau3MuRelativeIsolation;
  std::vector<TH2D>   BDT_2Dscan_B_TauHDecayMode;
  std::vector<TH2D>   BDT_2Dscan_B_VisibleDiTauMass;
  std::vector<TH2D>   BDT_2Dscan_B_MTT;
  std::vector<TH2D>   BDT_2Dscan_B_TripletMass;
  
  std::vector<TH2D>   BDT_2Dscan_B_TripletPt;
  std::vector<TH2D>   BDT_2Dscan_B_OppositeTauPt;
  std::vector<TH2D>   BDT_2Dscan_B_TripletEta;
  std::vector<TH2D>   BDT_2Dscan_B_OppositeTauEta;
  
  std::vector<TH2D>   BDT_2Dscan_B_MET_Et;
  std::vector<TH2D>   BDT_2Dscan_B_MET_Phi;
  std::vector<TH2D>   BDT_2Dscan_B_MET_Phi_vs_NeutrinoPhi;
  std::vector<TH2D>   BDT_2Dscan_B_MET_vs_NeutrinoPt;
  
  std::vector<TH2D>   BDT_2Dscan_B_Mu1_Pt;
  std::vector<TH2D>   BDT_2Dscan_B_Mu1_Eta;
  std::vector<TH2D>   BDT_2Dscan_B_Mu2_Pt;
  std::vector<TH2D>   BDT_2Dscan_B_Mu2_Eta;
  std::vector<TH2D>   BDT_2Dscan_B_Mu3_Pt;
  std::vector<TH2D>   BDT_2Dscan_B_Mu3_Eta;
  std::vector<TH2D>   BDT_2Dscan_B_h_Pt;
  std::vector<TH2D>   BDT_2Dscan_B_h_Eta;
  
  std::vector<TH2D>   BDT_2Dscan_B_FLSignificance;
  std::vector<TH2D>   BDT_2Dscan_B_SVPVTauDirAngle;
  std::vector<TH2D>   BDT_2Dscan_B_SVPVTauDirAngle_largescale;
  std::vector<TH2D>   BDT_2Dscan_B_VertexChi2KF;
  std::vector<TH2D>   BDT_2Dscan_B_MinDistToIsoTrack;
  std::vector<TH2D>   BDT_2Dscan_B_Kinematics_MissingTrMass;
  std::vector<TH2D>   BDT_2Dscan_B_Kinematics_MissingTrMass_cos;
  std::vector<TH2D>   BDT_2Dscan_B_Kinematics_MissingTrMass_pT;
  std::vector<TH2D>   BDT_2Dscan_B_Kinematics_MissingTrMass_MET;
  std::vector<TH2D>   BDT_2Dscan_B_VisibleDiTauMass_Collinear;
  
  std::vector<TH2D>   BDT_2Dscan_B_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH2D>   BDT_2Dscan_B_prod_size;
  
  std::vector<TH2D>   BDT_2Dscan_B_PairMass_OppositeSign_dR12;
  std::vector<TH2D>   BDT_2Dscan_B_PairMass_OppositeSign_dR13;

  
  Float_t m3m;
  Float_t dataMCtype;
  Float_t ifCommonCV;
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
  Float_t var_MinDrToIsoTrack;
  Float_t var_MinDistToIsoTrack;
  Float_t var_DeltaPhi;
  Float_t var_Phi_To_Opposite_Side;
  
  
  Float_t var_MET_Et;
  Float_t var_MET_Phi;
  
  Float_t var_VisMass;
  Float_t var_DiTauMass_Collinear;
  
  Float_t var_Whether_HPS_Tau_Vtx_Exists;
  Float_t var_HPS_FL_Sig;
  Float_t var_HPS_Inv_Mass_Z_Tau3mu_SpecTau;
  Float_t var_HPS_GJ_Angle_Ratio;
  
  Float_t BDT_Evaluated_A;
  Float_t BDT_Evaluated_B;
  
  Float_t BDT_Evaluated_MC_Bkg;
  
  //For combine
  Float_t tripletMass;
  Float_t bdt_cv;
  Float_t category;
  Float_t isMC;
  Float_t weight;
  Float_t dimu_OS1;
  Float_t dimu_OS2;
  
  TMVA::Reader *reader_Tauh;
  TMVA::Reader *reader_Tauh_CV;
  TMVA::Reader *reader_Tauh_NoCV;
  TMVA::Reader *reader_Tauh_MC_Bkg;
  
  TTree *T3MMiniTree_A;
  TTree *T3MMiniTree_B;
  TFile *T3MFMiniTree;
  
  TTree *T3MCombineTree_A;
  TTree *T3MCombineTree_B;
  TFile *T3MFCombineTree;
  
  TRandom rndm;
  double random_num;
  
  std::vector<std::vector<TH1D>*> InputFeatureCollection;
  std::vector<std::vector<TH2D>*> InputFeatureCollection_2D;

};
#endif
