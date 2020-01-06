#ifndef Isolation_h
#define Isolation_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class Isolation : public Selection {

 public:
  Isolation(TString Name_, TString id_);
  virtual ~Isolation();

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
  std::vector<TH1D> SVSecondPVTauDirAngle;

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
  std::vector<TH2D> IsoVsDr;

  std::vector<TH1D> Iso02;
  std::vector<TH1D> Iso04;
  std::vector<TH1D> Iso06;
  std::vector<TH1D> Iso08;
  std::vector<TH1D> Iso1;
  std::vector<TH1D> Iso12;
  std::vector<TH1D> Iso14;
  std::vector<TH1D> Iso16;
  std::vector<TH1D> Iso18;
  std::vector<TH1D> Iso2;

  std::vector<TH2D> Iso18VSPU;
  std::vector<TH2D> Iso12VSPU;
  std::vector<TH2D> Iso02VSPU;



  std::vector<TH1D> Iso02Mu1;
  std::vector<TH1D> Iso04Mu1;
  std::vector<TH1D> Iso06Mu1;
  std::vector<TH1D> Iso08Mu1;
  std::vector<TH1D> Iso1Mu1;
  std::vector<TH1D> Iso12Mu1;
  std::vector<TH1D> Iso14Mu1;
  std::vector<TH1D> Iso16Mu1;
  std::vector<TH1D> Iso18Mu1;
  std::vector<TH1D> Iso2Mu1;
  std::vector<TH2D> Iso12Mu1VSPU;
  std::vector<TH2D> Iso02Mu1VSPU;
  std::vector<TH2D> Iso18Mu1VSPU;

  std::vector<TH1D> Iso02Mu2;
  std::vector<TH1D> Iso04Mu2;
  std::vector<TH1D> Iso06Mu2;
  std::vector<TH1D> Iso08Mu2;
  std::vector<TH1D> Iso1Mu2;
  std::vector<TH1D> Iso12Mu2;
  std::vector<TH1D> Iso14Mu2;
  std::vector<TH1D> Iso16Mu2;
  std::vector<TH1D> Iso18Mu2;
  std::vector<TH1D> Iso2Mu2;

  std::vector<TH1D> Iso02Mu3;
  std::vector<TH1D> Iso04Mu3;
  std::vector<TH1D> Iso06Mu3;
  std::vector<TH1D> Iso08Mu3;
  std::vector<TH1D> Iso1Mu3;
  std::vector<TH1D> Iso12Mu3;
  std::vector<TH1D> Iso14Mu3;
  std::vector<TH1D> Iso16Mu3;
  std::vector<TH1D> Iso18Mu3;
  std::vector<TH1D> Iso2Mu3;

  std::vector<TH1D> Iso08MuMax;
  std::vector<TH1D> Iso08MuMin;

  std::vector<TH1D> MindcaTrackSV;
  std::vector<TH2D> MindcaTrackSVNPU;
  std::vector<TH1D> MindcaTrackSV2;
  std::vector<TH1D> dcaTrackPV;
  std::vector<TH1D> PVVerticesDeltaZ;
  std::vector<TH1D> PVZResolution;
  std::vector<TH1D> PV2ZResolution;
  std::vector<TH1D> SVDeltaR;
  std::vector<TH1D> SVDistance;
  std::vector<TH1D> NSV;
  std::vector<TH2D> MatchedSV_Mass;
  std::vector<TH2D> SVMassVsDistanceToSV;
  std::vector<TH1D> MatchedSV_MassOS;

  std::vector<TH1D> MaxD0SigPV;
  std::vector<TH1D> MinD0SigPV;

  std::vector<TH1D> MaxD0SigBS;
  std::vector<TH1D> MinD0SigBS;

  std::vector<TH1D> MaxD0SigSV;
  std::vector<TH1D> MinD0SigSV;


  std::vector<TH1D> deltaMuZ12;
  std::vector<TH1D> deltaMuZ13;
  std::vector<TH1D> deltaMuZ23;
  std::vector<TH1D> MaxdeltaMuZ;
  std::vector<TH1D> MindeltaMuZ;



  std::vector<TH1D> PVTrackDz;
  std::vector<TH1D> NtracksClose;
  std::vector<TH2D> NtracksCloseVSPU;
  std::vector<TH1D> MuMatchedTrackMass;

  std::vector<TH1D> MuMuMindR;
  std::vector<TH1D> MuMuMaxdR;

  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> OS1_muon_Mass;
  std::vector<TH1D> OS2_muon_Mass;
  std::vector<TH1D> SS2_muon_Mass;
  std::vector<TH1D> OS_muon_Mass_max;


  std::vector<TH1D> Muon1isStand;
  std::vector<TH1D> Muon2isStand;
  std::vector<TH1D> Muon3isStand;


  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;

  std::vector<TH1D> MaxVertexPairQuality;
  std::vector<TH1D> MinVertexPairQuality;


  std::vector<TH1D> Muon1PtResolution;
  std::vector<TH1D> Muon2PtResolution;
  std::vector<TH1D> Muon3PtResolution;

  std::vector<TH1D> Muon1EtaResolution;
  std::vector<TH1D> Muon2EtaResolution;
  std::vector<TH1D> Muon3EtaResolution;

  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;


  std::vector<TH1D> Mu1TrackMass;
  std::vector<TH1D> Mu2TrackMass;
  std::vector<TH1D> Mu3TrackMass;
  std::vector<TH1D> SV_Mass;
  std::vector<TH1D> SV_Mass2OS;
  std::vector<TH1D> MinSegmentCompatibility;
  std::vector<TH1D> MinMuon_chi2LocalPosition;
  std::vector<TH1D> MaxMuon_chi2LocalPosition;

  std::vector<TH1D> MinMuon_chi2LocalMomentum;
  std::vector<TH1D> MaxMuon_chi2LocalMomentum;

  std::vector<TH1D> MaxDca;
  std::vector<TH1D> MinDca;
  std::vector<TH1D> MinMatchedStations;
  std::vector<TH1D> MinTrackerLayers;
 
  std::vector<TH1D> MintrkKink;
  std::vector<TH1D> MaxtrkKink;
  std::vector<TH1D> MinglbKink;
  std::vector<TH1D> MaxglbKink;




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
