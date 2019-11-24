#ifndef MinBiasSelector_h
#define MinBiasSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class MinBiasSelector : public Selection {

 public:
  MinBiasSelector(TString Name_, TString id_);
  virtual ~MinBiasSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  // Selection Variables
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

  std::vector<TH1D> TauMassRefitA1;
  std::vector<TH1D> TauMassA1;
  std::vector<TH1D> TauMassRefitB1;
  std::vector<TH1D> TauMassB1;
  std::vector<TH1D> TauMassRefitC1;
  std::vector<TH1D> TauMassC1;

  std::vector<TH1D> TauMassRefitA2;
  std::vector<TH1D> TauMassA2;
  std::vector<TH1D> TauMassRefitB2;
  std::vector<TH1D> TauMassB2;
  std::vector<TH1D> TauMassRefitC2;
  std::vector<TH1D> TauMassC2;




  std::vector<TH1D> TauMassRefitBarrel1;
  std::vector<TH1D> TauMassBarrel1;

  std::vector<TH1D> TauMassRefitEndcap1;
  std::vector<TH1D> TauMassEndcap1;

  std::vector<TH1D> TauMassRefitBarrel2;
  std::vector<TH1D> TauMassBarrel2;

  std::vector<TH1D> TauMassRefitEndcap2;
  std::vector<TH1D> TauMassEndcap2;

  std::vector<TH1D> TauMassResolutionRefit;

  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> TauMass_all;
  std::vector<TH2D> TauMass_allVsBDTA;
  std::vector<TH2D> TauMass_allVsBDTB;
  std::vector<TH2D> TauMass_allVsBDTC;
  std::vector<TH2D> TauMass_allVsBDTBarrel;
  std::vector<TH2D> TauMass_allVsBDTEndcap;
  std::vector<TH2D> EMR_tau_eta;
  std::vector<TH2D> L1Triggers;


  std::vector<TH2D> L1TriggersB;
  std::vector<TH2D> L1TriggersC;
  std::vector<TH2D> L1TriggersD;
  std::vector<TH2D> L1TriggersE;
  std::vector<TH2D> L1TriggersF;


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

  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;


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
  std::vector<TH1D> BDTOutputA;
  std::vector<TH1D> BDTOutputB;
  std::vector<TH1D> BDTOutputC;
  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;

  TMVA::Reader *readerA;
  TMVA::Reader *readerB;
  TMVA::Reader *readerC;
  TMVA::Reader *readerBarrel;
  TMVA::Reader *readerEndcap;

  Float_t var_vertexKFChi2;// (chi sq of the fit of the secondary vertex)
  Float_t var_svpvTauAngle;// (The angle between PV-SV vector and the tau vector)
  Float_t var_flightLenSig;// (Flight length significance of the tau candidate)
  Float_t var_sumMuTrkKinkChi2;// (sum of chi sq of the kink of all three muons)
  Float_t var_segCompMuMin;// (Minimum of the segment compatibility of the three muons)
  Float_t var_MinMIPLikelihood;// (Minimum of the calorimeter compatibility of the three muons)
  Float_t var_tauMass;
  Float_t var_MuMu_mindR;
  Float_t var_RelPt_Mu1Tau;
  Float_t var_Eta_au;
  Float_t var_MuMu_minKFChi2;
  Float_t var_maxdca;
  Float_t var_MuTau_maxdR;
  Float_t var_MaxD0Significance;
  Float_t var_IsolationMinDist;


  Float_t m3m;
  Float_t dataMCtype;
  Float_t event_weight;
  Float_t bdt;
  Float_t category;
  Float_t rapidity;
  Float_t LumiScale;




};
#endif
