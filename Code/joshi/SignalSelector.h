#ifndef SignalSelector_h
#define SignalSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class SignalSelector : public Selection {

 public:
  SignalSelector(TString Name_, TString id_);
  virtual ~SignalSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PVRefit, PhiVeto, OmegaVeto, TriggerMatchMu1, TriggerMatchMu2,TriggerMatchMu3,TauMassCut,GenMatch,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double tauMassResCutLow, tauMassResCutHigh;
  double phiVetoSigma, omegaVetoSigma;
  
  float bdt_cutA1_3glb;  float bdt_cutA2_3glb;  float bdt_cutB1_3glb;  float bdt_cutB2_3glb;  float bdt_cutC1_3glb;  float bdt_cutC2_3glb;  float bdt_cutA1_2glbTrk;  float bdt_cutA2_2glbTrk;  float bdt_cutB1_2glbTrk;  float bdt_cutB2_2glbTrk;  float bdt_cutC1_2glbTrk;  float bdt_cutC2_2glbTrk;

 // categorization variables
 bool MC;
 float category;
 Int_t threeGlobal;

 //commmon variables (2016 + 2017)
 float var_vertexKFChi2; // <= should be changed to normalized KF chi2
 float var_svpvTauAngle; 
 float var_flightLenSig;
 float var_segCompMuMin;
	 
 // 2016 variables
 float var_pmin; // Minimum p of the three muons
 float var_max_cLP; // Maximum chi square of the STA-TRK matching
 float var_max_tKink; // Maximum of the track kink of the 3 muons
 float var_MinD0Significance; // Minimum of the transverse IP significance of the 3 muons
 float var_mindca_iso; // Minimum DCA of tracks to muons with pT > 1 GeV (which muon?)
 float var_trk_relPt; // Ratio of sum of Pt of the tracks in muon isolation to muon (max value) [trk_pt>1 GeV, dR<0.03, dca<1 mm]
 float var_MinMIPLikelihood;


 // 2017 variables
 float var_MuMu_minKFChi2;
 float var_MuTau_maxdR;
 float var_sumMuTrkKinkChi2; // sum of chi square of STA-TRK matching of 3 muons
 float var_MaxD0Significance; // Minimum of the transverse IP significance of the 3 muons
 
 float var_Eta_Tau; 
 float var_tauMassRes;
 float var_tauMass;
 float var_MuMu_mindR;
 float var_RelPt_Mu1Ta;
 float var_maxdca;
 float var_RelPt_Mu1Tau;
 float var_tauMassRefit;

 // MiniTree variables
 Float_t m3m;
  Float_t dataMCtype;
  Float_t event_weight;
  Float_t bdt;
  Float_t rapidity;
  Float_t LumiScale;

  // file and ttree pointer
  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  // Initializhere your analysis histograms

  std::vector<TH1D> Muon1Pt;
  std::vector<TH1D> Muon2Pt;
  std::vector<TH1D> Muon3Pt;

  std::vector<TH1D> Muon1Eta;
  std::vector<TH1D> Muon2Eta;
  std::vector<TH1D> Muon3Eta;


  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauP;

  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefitA1_3glb;
  std::vector<TH1D> TauMassA1_3glb;
  std::vector<TH1D> TauMassRefitB1_3glb;
  std::vector<TH1D> TauMassB1_3glb;
  std::vector<TH1D> TauMassRefitC1_3glb;
  std::vector<TH1D> TauMassC1_3glb;

  std::vector<TH1D> TauMassRefitA2_3glb;
  std::vector<TH1D> TauMassA2_3glb;
  std::vector<TH1D> TauMassRefitB2_3glb;
  std::vector<TH1D> TauMassB2_3glb;
  std::vector<TH1D> TauMassRefitC2_3glb;
  std::vector<TH1D> TauMassC2_3glb;

  std::vector<TH1D> TauMassRefitA1_2glbtrk;
  std::vector<TH1D> TauMassA1_2glbtrk;
  std::vector<TH1D> TauMassRefitB1_2glbtrk;
  std::vector<TH1D> TauMassB1_2glbtrk;
  std::vector<TH1D> TauMassRefitC1_2glbtrk;
  std::vector<TH1D> TauMassC1_2glbtrk;

  std::vector<TH1D> TauMassRefitA2_2glbtrk;
  std::vector<TH1D> TauMassA2_2glbtrk;
  std::vector<TH1D> TauMassRefitB2_2glbtrk;
  std::vector<TH1D> TauMassB2_2glbtrk;
  std::vector<TH1D> TauMassRefitC2_2glbtrk;
  std::vector<TH1D> TauMassC2_2glbtrk;
  
  std::vector<TH1D> BDTOutputA_3glb;
  std::vector<TH1D> BDTOutputB_3glb;
  std::vector<TH1D> BDTOutputC_3glb;
  
  std::vector<TH1D> BDTOutputA_2glbtrk;
  std::vector<TH1D> BDTOutputB_2glbtrk;
  std::vector<TH1D> BDTOutputC_2glbtrk;
  
  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;
  std::vector<TH1D> TauMassRefitBarrel1;
  std::vector<TH1D> TauMassBarrel1;

  std::vector<TH1D> TauMassRefitEndcap1;
  std::vector<TH1D> TauMassEndcap1;

  std::vector<TH1D> TauMassRefitBarrel2;
  std::vector<TH1D> TauMassBarrel2;

  std::vector<TH1D> TauMassRefitEndcap2;
  std::vector<TH1D> TauMassEndcap2;
  std::vector<TH1D> TauMassResolution_3glb;
  std::vector<TH1D> TauMassResolution_2glbtrk;
  std::vector<TH1D> TauMassResolutionRefit_3glb;
  std::vector<TH1D> TauMassResolutionRefit_2glbtrk;

  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> TauMass_all;
  
  std::vector<TH2D> TauMass_allVsBDTA_3glb;
  std::vector<TH2D> TauMass_allVsBDTB_3glb;
  std::vector<TH2D> TauMass_allVsBDTC_3glb;
  
  std::vector<TH2D> TauMass_allVsBDTA_2glbtrk;
  std::vector<TH2D> TauMass_allVsBDTB_2glbtrk;
  std::vector<TH2D> TauMass_allVsBDTC_2glbtrk;

  std::vector<TH2D> TauMass_allVsBDTBarrel;
  std::vector<TH2D> TauMass_allVsBDTEndcap;
  std::vector<TH2D> EMR_tau_eta;
  std::vector<TH2D> L1Triggers;

  std::vector<TH2D> L1Triggers2017B;
  std::vector<TH2D> L1Triggers2017C;
  std::vector<TH2D> L1Triggers2017D;
  std::vector<TH2D> L1Triggers2017E;
  std::vector<TH2D> L1Triggers2017F;

  std::vector<TH2D> L1Triggers2018A;
  std::vector<TH2D> L1Triggers2018B;
  std::vector<TH2D> L1Triggers2018C;
  std::vector<TH2D> L1Triggers2018D;

  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> Muon1isStand;
  std::vector<TH1D> Muon2isStand;
  std::vector<TH1D> Muon3isStand;

  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;

  std::vector<TH1D> Muon1PtResolution;
  std::vector<TH1D> Muon2PtResolution;
  std::vector<TH1D> Muon3PtResolution;

  std::vector<TH1D> Muon1EtaResolution;
  std::vector<TH1D> Muon2EtaResolution;
  std::vector<TH1D> Muon3EtaResolution;

  std::vector<TH1D> Muon1DRToTruth_3glb;
  std::vector<TH1D> Muon2DRToTruth_3glb;
  std::vector<TH1D> Muon3DRToTruth_3glb;

  std::vector<TH1D> Muon1DRToTruth_2glbtrk;
  std::vector<TH1D> Muon2DRToTruth_2glbtrk;
  std::vector<TH1D> Muon3DRToTruth_2glbtrk;
  
  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;

  std::vector<TH1D> dR12;
  std::vector<TH1D> dR23;
  std::vector<TH1D> dR31;
  std::vector<TH1D> dR1Tau;
  std::vector<TH1D> dR2Tau;
  std::vector<TH1D> dR3Tau;

  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> FLSignificance;

  TMVA::Reader *readerA_3glb;
  TMVA::Reader *readerB_3glb;
  TMVA::Reader *readerC_3glb;
  TMVA::Reader *readerA_2glbTrk;
  TMVA::Reader *readerB_2glbTrk;
  TMVA::Reader *readerC_2glbTrk;

};
#endif
