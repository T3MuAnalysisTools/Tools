#ifndef BParkingSelector_h
#define BParkingSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class BParkingSelector : public Selection {

 public:
  BParkingSelector(TString Name_, TString id_);
  virtual ~BParkingSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch1, TriggerMatch2,TriggerMatch3,TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  // Trigger counters
  Int_t doubleMu_counter = 0;
  Int_t singleMu_counter = 0;
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
  std::vector<TH2D> TriggersDecision;
  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefit_3glb_A1;
  std::vector<TH1D> TauMass_3glb_A1;
  std::vector<TH1D> TauMassRefit_3glb_B1;
  std::vector<TH1D> TauMass_3glb_B1;
  std::vector<TH1D> TauMassRefit_3glb_C1;
  std::vector<TH1D> TauMass_3glb_C1;

  std::vector<TH1D> TauMassRefit_3glb_A2;
  std::vector<TH1D> TauMass_3glb_A2;
  std::vector<TH1D> TauMassRefit_3glb_B2;
  std::vector<TH1D> TauMass_3glb_B2;
  std::vector<TH1D> TauMassRefit_3glb_C2;
  std::vector<TH1D> TauMass_3glb_C2;

  std::vector<TH1D> TauMassRefit_2glbTrk_A1;
  std::vector<TH1D> TauMass_2glbTrk_A1;
  std::vector<TH1D> TauMassRefit_2glbTrk_B1;
  std::vector<TH1D> TauMass_2glbTrk_B1;
  std::vector<TH1D> TauMassRefit_2glbTrk_C1;
  std::vector<TH1D> TauMass_2glbTrk_C1;

  std::vector<TH1D> TauMassRefit_2glbTrk_A2;
  std::vector<TH1D> TauMass_2glbTrk_A2;
  std::vector<TH1D> TauMassRefit_2glbTrk_B2;
  std::vector<TH1D> TauMass_2glbTrk_B2;
  std::vector<TH1D> TauMassRefit_2glbTrk_C2;
  std::vector<TH1D> TauMass_2glbTrk_C2;

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
  std::vector<TH2D> TauMass_allVsBDT_3glb_A;
  std::vector<TH2D> TauMass_allVsBDT_3glb_B;
  std::vector<TH2D> TauMass_allVsBDT_3glb_C;
  std::vector<TH2D> TauMass_allVsBDT_2glbTrk_A;
  std::vector<TH2D> TauMass_allVsBDT_2glbTrk_B;
  std::vector<TH2D> TauMass_allVsBDT_2glbTrk_C;
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
  std::vector<TH1D> BDTOutput_3glb_A;
  std::vector<TH1D> BDTOutput_3glb_B;
  std::vector<TH1D> BDTOutput_3glb_C;
  std::vector<TH1D> BDTOutput_2glbTrk_A;
  std::vector<TH1D> BDTOutput_2glbTrk_B;
  std::vector<TH1D> BDTOutput_2glbTrk_C;
  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;

  TMVA::Reader *readerA_3glb;
  TMVA::Reader *readerB_3glb;
  TMVA::Reader *readerC_3glb;
  TMVA::Reader *readerA_2glbTrk;
  TMVA::Reader *readerB_2glbTrk;
  TMVA::Reader *readerC_2glbTrk;

  Float_t var_vertexKFChi2;
  Float_t var_svpvTauAngle;
  Float_t var_flightLenSig;
  Float_t var_segCompMuMin;
  Float_t var_MinD0Significance;
  Float_t var_mindca_iso;
  
  Float_t var_pmin;
  Float_t var_tauMass;
  Float_t var_tauMassRes;
  Float_t var_Eta_Tau;
  Float_t var_max_cLP;
  Float_t var_max_tKink;
  Float_t var_trk_relPt;

  Float_t m3m;
  Float_t dataMCtype;
  Float_t event_weight;
  Float_t bdt;
  Float_t category;
  Float_t rapidity;
  Float_t LumiScale;


};
#endif
