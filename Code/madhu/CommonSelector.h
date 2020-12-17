#ifndef CommonSelector_h
#define CommonSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class CommonSelector : public Selection {

 public:
  CommonSelector(TString Name_, TString id_);
  virtual ~CommonSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1T=0,HLT,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID,PhiVeto1, OmegaVeto1, PhiVeto2, OmegaVeto2, TriggerMatch, TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double phiVetoCut1, phiVetoCut2, rmgCutVeto1, rmgCutVeto2;
  double PEMassResolutionCut1_,PEMassResolutionCut2_;
  double mvaA1_,mvaA2_,mvaB1_,mvaB2_,mvaC1_,mvaC2_;





  bool RunB, RunC, RunD, RunE, RunF = 0;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  TRandom rndm;
  double random_num;
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
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVA;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVA;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVASV;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVASV;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVASV;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVASVAngle0;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVASVAngle0;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVASVAngle0;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVASVAngle1;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVASVAngle1;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVASVAngle1;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVASVAngle2;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVASVAngle2;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVASVAngle2;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVASVAngle3;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVASVAngle3;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVASVAngle3;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarSV;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarSV;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarSV;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarSVAngle0;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarSVAngle0;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarSVAngle0;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarSVAngle1;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarSVAngle1;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarSVAngle1;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarSVAngle2;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarSVAngle2;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarSVAngle2;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarSVAngle3;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarSVAngle3;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarSVAngle3;  
  
  std::vector<TH1D> TauAngleTest;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAFiner;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAFiner;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAFiner;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVACoarser;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVACoarser;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVACoarser;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStar;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStar;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStar;
  
  std::vector<TH1D> MuOSSS1InvariantMassBeforeMVA;
  std::vector<TH1D> MuOSSS2InvariantMassBeforeMVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdR;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdR;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdR;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdR;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdR;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdR;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRIncrease;
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRDecrease;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd1;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd1;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd1;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd2;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd2;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd2;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd3;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd3;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd3;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd4;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd4;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd4;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd5;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd5;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd5;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd6;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd6;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd6;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestdRLtd7;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestdRLtd7;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestdRLtd7;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7;
  
  std::vector<TH1D> IsolationTrackCount;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVABestMass;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVABestMass;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVABestMass;
  
  std::vector<TH1D> Mu1TrackInvariantMassBeforeMVAKStarBestMass;
  std::vector<TH1D> Mu2TrackInvariantMassBeforeMVAKStarBestMass;
  std::vector<TH1D> Mu3TrackInvariantMassBeforeMVAKStarBestMass;
  
  
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterA1MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterA1MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterA1MVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterB1MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterB1MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterB1MVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterC1MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterC1MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterC1MVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterA2MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterA2MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterA2MVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterB2MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterB2MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterB2MVA;
  
  std::vector<TH1D> Mu1TrackInvariantMassAfterC2MVA;
  std::vector<TH1D> Mu2TrackInvariantMassAfterC2MVA;
  std::vector<TH1D> Mu3TrackInvariantMassAfterC2MVA;




  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefitA1;
  std::vector<TH1D> TauMassRefitA1MassCut;
  std::vector<TH1D> TauMassRefitA2MassCut;
  std::vector<TH1D> TauMassRefitA1HalfMassCut;
  std::vector<TH1D> TauMassRefitA2HalfMassCut;
  std::vector<TH1D> TauMassRefitA1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitA2FullEtaVetoCut;
  std::vector<TH1D> TauMassA1;
  std::vector<TH1D> TauMassRefitB1;
  std::vector<TH1D> TauMassRefitB1MassCut;
  std::vector<TH1D> TauMassRefitB2MassCut;
  std::vector<TH1D> TauMassRefitB1HalfMassCut;
  std::vector<TH1D> TauMassRefitB2HalfMassCut;

  std::vector<TH1D> TauMassRefitB1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitB2FullEtaVetoCut;
  std::vector<TH1D> TauMassB1;
  std::vector<TH1D> TauMassRefitC1;
  std::vector<TH1D> TauMassRefitC1MassCut;
  std::vector<TH1D> TauMassRefitC2MassCut;

  std::vector<TH1D> TauMassRefitC1HalfMassCut;
  std::vector<TH1D> TauMassRefitC2HalfMassCut;

  std::vector<TH1D> TauMassRefitC1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitC2FullEtaVetoCut;
  std::vector<TH1D> TauMassC1;

  std::vector<TH1D> TauMassRefitA2;
  std::vector<TH1D> TauMassA2;
  std::vector<TH1D> TauMassRefitB2;
  std::vector<TH1D> TauMassB2;
  std::vector<TH1D> TauMassRefitC2;
  std::vector<TH1D> TauMassC2;



  std::vector<TH1D> TauMassRefitABC1;
  std::vector<TH1D> TauMassRefitABC2;


  std::vector<TH1D> TauMassResolutionRefit;

  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> TauMass_all;
  std::vector<TH2D> TauMass_allVsBDTA;
  std::vector<TH2D> TauMass_allVsBDTB;
  std::vector<TH2D> TauMass_allVsBDTC;
  std::vector<TH2D> TauMass_allVsBDTBarrel;
  std::vector<TH2D> TauMass_allVsBDTEndcap;
  std::vector<TH2D> EMR_tau_eta;

  std::vector<TH2D> PairMass;
  std::vector<TH2D> PairMassFinalSel;
  std::vector<TH1D> PairMass1;
  std::vector<TH1D> PairMass2;

  std::vector<TH2D> PairMassEta;
  std::vector<TH2D> PairMassEtaPrime;

  std::vector<TH2D> PairMassWithCut;
  std::vector<TH2D> CategoryOverlap;
  std::vector<TH2D> Mu3IdOverlap;
  std::vector<TH1D> IDOriginOfOSMuon;

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


  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> BDTOutputA;
  std::vector<TH1D> BDTOutputB;
  std::vector<TH1D> BDTOutputC;
  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;


  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;
  std::vector<TH1D>  Muon3MVAID;



  TMVA::Reader *readerA;
  TMVA::Reader *readerB;
  TMVA::Reader *readerC;


  TMVA::Reader *readerMuIDBarrel;
  TMVA::Reader *readerMuIDEndcap;



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
  Float_t var_Muon1DetID;
  Float_t var_Muon2DetID;
  Float_t var_Muon3DetID;
  Float_t var_mass12_dRsorting;
  Float_t var_mass13_drSorting;


  

  float var_MaxtrkKink;
  float var_MaxD0SigSV;
  float var_MindcaTrackSV;
  float var_maxMuonsDca;
  float var_nsv;

  float var_MaxMuon_chi2LocalPosition;
  float var_MaxVertexPairQuality;
  float var_MuonglbkinkSum;
  float var_MaxMuon_chi2LocalMomentum;



  Double_t m3m;
  Double_t dataMCtype;
  Double_t event_weight;
  Double_t bdt;
  Double_t category;
  Double_t m12;
  Double_t m13;
  Double_t LumiScale;
  Double_t mvaA1;
  Double_t mvaA2;
  Double_t mvaB1;
  Double_t mvaB2;
  Double_t mvaC1;
  Double_t mvaC2;

  Float_t mu_combinedQuality_chi2LocalMomentum;
  Float_t mu_combinedQuality_chi2LocalPosition;
  Float_t mu_combinedQuality_staRelChi2;
  Float_t mu_combinedQuality_trkRelChi2;
  Float_t mu_combinedQuality_globalDeltaEtaPhi;
  Float_t mu_combinedQuality_trkKink;
  Float_t mu_combinedQuality_glbKink;
  Float_t mu_combinedQuality_glbTrackProbability;
  Float_t mu_Numberofvalidtrackerhits;
  Float_t mu_Numberofvalidpixelhits;
  Float_t mu_validMuonHitComb;
  Float_t mu_numberOfMatchedStations;
  Float_t mu_segmentCompatibility;
  Float_t mu_timeAtIpInOutErr;
  Float_t mu_GLnormChi2;
  Float_t mu_innerTrack_normalizedChi2;
  Float_t mu_outerTrack_normalizedChi2;
  Float_t mu_innerTrack_validFraction;

  Float_t mu_eta;
  Float_t mu_pt;
  Float_t mu_phi;
  Float_t mu_SoftMVA;


  Float_t Muon1DetID;
  Float_t Muon2DetID;
  Float_t Muon3DetID;



};
#endif
