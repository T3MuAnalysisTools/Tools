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

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


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
  
  //After Selection
  
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
  std::vector<TH1D>   PostSelection_MinDistToIsoTrack_cm;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostSelection_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostSelection_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostSelection_Neutrino_Test_Phi;
  std::vector<TH1D>   PostSelection_Neutrino_Test_Et;
  std::vector<TH1D>   PostSelection_Neutrino_Test_Mass;
  std::vector<TH1D>   PostSelection_Neutrino_Test_Phi_Opposing;
  
  std::vector<TH1D>   PostSelection_prod_size;
  
  std::vector<TH1D>   PostSelection_Vertex_Dist;
  std::vector<TH1D>   PostSelection_Vertex_Chi2;
  
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_peakspectra1;
  std::vector<TH1D>   PostSelection_VisibleDiTauMass_peakspectra2;
  
  std::vector<TH1D>   PostSelection_BDT_Output;
  std::vector<TH1D>   PostSelection_BDT_Output_Without_Vis_Mass;
  std::vector<TH1D>   PostSelection_BDT_Output_MC_Bkg;
  std::vector<TH2D>   PostSelection_BDT_Output_Data_vs_MC_Bkg;
  
  std::vector<TH1D>   PostSelection_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostSelection_PairMass_OppositeSign_dR13;
  
  //After BDT
  
  std::vector<TH1D>   PostBDT_TripletMass_VeryLooseCut;
  
  std::vector<TH1D>   PostBDT_Tau3MuRelativeIsolation;
  std::vector<TH1D>   PostBDT_OppositeMuRelativeIsolation;
  std::vector<TH1D>   PostBDT_VisibleDiTauMass;
  std::vector<TH1D>   PostBDT_MTT;
  std::vector<TH1D>   PostBDT_TripletMass;
  std::vector<TH1D>   PostBDT_TripletMass_a;
  std::vector<TH1D>   PostBDT_TripletMass_b;
  
  std::vector<TH1D>   PostBDT_TripletPt;
  std::vector<TH1D>   PostBDT_OppositeMuonPt;
  std::vector<TH1D>   PostBDT_TripletEta;
  std::vector<TH1D>   PostBDT_OppositeMuonEta;
  
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
  std::vector<TH1D>   PostBDT_mu_Pt;
  std::vector<TH1D>   PostBDT_mu_Eta;
  
  std::vector<TH1D>   PostBDT_FLSignificance;
  std::vector<TH1D>   PostBDT_SVPVTauDirAngle;
  std::vector<TH1D>   PostBDT_SVPVTauDirAngle_largescale;
  std::vector<TH1D>   PostBDT_VertexChi2KF;
  std::vector<TH1D>   PostBDT_MinDistToIsoTrack;
  std::vector<TH1D>   PostBDT_MinDistToIsoTrack_cm;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_cos;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_pT;
  std::vector<TH1D>   PostBDT_Kinematics_MissingTrMass_MET;
  std::vector<TH1D>   PostBDT_VisibleDiTauMass_Collinear;
  
  std::vector<TH1D>   PostBDT_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH1D>   PostBDT_Neutrino_Test_Phi;
  std::vector<TH1D>   PostBDT_Neutrino_Test_Et;
  std::vector<TH1D>   PostBDT_Neutrino_Test_Mass;
  std::vector<TH1D>   PostBDT_Neutrino_Test_Phi_Opposing;
  
  std::vector<TH1D>   PostBDT_prod_size;
  
  std::vector<TH1D>   PostBDT_Vertex_Dist;
  std::vector<TH1D>   PostBDT_Vertex_Chi2;
  
  std::vector<TH1D>   PostBDT_VisibleDiTauMass_peakspectra1;
  std::vector<TH1D>   PostBDT_VisibleDiTauMass_peakspectra2;
  
  std::vector<TH1D>   PostBDT_PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PostBDT_PairMass_OppositeSign_dR13;
  
  //Pre BDT 2D scan
  
  std::vector<TH2D>   BDT_2Dscan_Tau3MuRelativeIsolation;
  std::vector<TH2D>   BDT_2Dscan_OppositeMuRelativeIsolation;
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass;
  std::vector<TH2D>   BDT_2Dscan_MTT;
  std::vector<TH2D>   BDT_2Dscan_TripletMass;
  
  std::vector<TH2D>   BDT_2Dscan_TripletPt;
  std::vector<TH2D>   BDT_2Dscan_OppositeMuonPt;
  std::vector<TH2D>   BDT_2Dscan_TripletEta;
  std::vector<TH2D>   BDT_2Dscan_OppositeMuonEta;
  
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
  std::vector<TH2D>   BDT_2Dscan_mu_Pt;
  std::vector<TH2D>   BDT_2Dscan_mu_Eta;
  
  std::vector<TH2D>   BDT_2Dscan_FLSignificance;
  std::vector<TH2D>   BDT_2Dscan_SVPVTauDirAngle;
  std::vector<TH2D>   BDT_2Dscan_SVPVTauDirAngle_largescale;
  std::vector<TH2D>   BDT_2Dscan_VertexChi2KF;
  std::vector<TH2D>   BDT_2Dscan_MinDistToIsoTrack;
  std::vector<TH2D>   BDT_2Dscan_MinDistToIsoTrack_cm;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_cos;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_pT;
  std::vector<TH2D>   BDT_2Dscan_Kinematics_MissingTrMass_MET;
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass_Collinear;
  
  std::vector<TH2D>   BDT_2Dscan_Phi_Triplet_to_Spectator_Tau;
  
  std::vector<TH2D>   BDT_2Dscan_Neutrino_Test_Phi;
  std::vector<TH2D>   BDT_2Dscan_Neutrino_Test_Et;
  std::vector<TH2D>   BDT_2Dscan_Neutrino_Test_Mass;
  std::vector<TH2D>   BDT_2Dscan_Neutrino_Test_Phi_Opposing;
  
  std::vector<TH2D>   BDT_2Dscan_prod_size;
  
  std::vector<TH2D>   BDT_2Dscan_Vertex_Dist;
  std::vector<TH2D>   BDT_2Dscan_Vertex_Chi2;
  
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass_peakspectra1;
  std::vector<TH2D>   BDT_2Dscan_VisibleDiTauMass_peakspectra2;
  
  std::vector<TH2D>   BDT_2Dscan_PairMass_OppositeSign_dR12;
  std::vector<TH2D>   BDT_2Dscan_PairMass_OppositeSign_dR13;
  
  
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
  
  Float_t var_MuonIsolation;
  Float_t var_Muon_pT;
  Float_t var_Muon_eta;
  
  Float_t var_FLSignificance;
  Float_t var_SVPVTauDirAngle;
  Float_t var_ThreeMuVertexChi2KF;
  Float_t var_MinDistToIsoTrack;
  Float_t var_DeltaPhi;//MET and Visible tau
  Float_t var_Phi_To_Opposite_Side;
  
  Float_t var_MET_Et;
  Float_t var_MET_Phi;
  
  Float_t var_VisMass;
  Float_t var_DiTauMass_Collinear;
  Float_t var_4Mu_Chi2;
  
  Float_t BDT_Evaluated;
  Float_t BDT_Evaluated_Without_Vis_Mass;
  Float_t BDT_Evaluated_MC_Bkg;
  
  //For combine
  Float_t tripletMass;
  Float_t bdt_cv;
  Float_t category;
  Float_t isMC;
  Float_t weight;
  Float_t dimu_OS1;
  Float_t dimu_OS2;
  
  TMVA::Reader *reader_Taumu;
  TMVA::Reader *reader_Taumu_Without_Vis_Mass;
  TMVA::Reader *reader_Taumu_MC_Bkg;
  
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
