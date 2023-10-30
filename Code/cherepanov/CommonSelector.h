#ifndef CommonSelector_h
#define CommonSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include <algorithm>
#include <set>
#include "EventClassifier.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "PDGInfo.h"



class CommonSelector : public Selection {

 public:
  CommonSelector(TString Name_, TString id_);
  virtual ~CommonSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1T=0,HLT,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID,PhiVeto1, OmegaVeto1, PhiVeto2, OmegaVeto2, TriggerMatch, TauMassCut,NCuts}; 
  enum DecayCategory{ThreeParticlesFrom1Origin=0, TwoParticlesFrom1Origin,  TwoParticlesFrom1OriginQED};

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double phiVetoCut1, phiVetoCut2, rmgCutVeto1, rmgCutVeto2;
  double PEMassResolutionCut1_,PEMassResolutionCut2_;
  double mvaA1_,mvaA2_,mvaB1_,mvaB2_,mvaC1_,mvaC2_;
  double mvaA1train11_,mvaA2train11_,mvaB1train11_,mvaB2train11_,mvaC1train11_,mvaC2train11_;
  double mvaA1train8_,mvaA2train8_,mvaB1train8_,mvaB2train8_,mvaC1train8_,mvaC2train8_;

  double mvaBTrainA1_,mvaBTrainA2_,mvaBTrainB1_,mvaBTrainB2_,mvaBTrainC1_,mvaBTrainC2_;
  double mvaDTrainA1_,mvaDTrainA2_,mvaDTrainB1_,mvaDTrainB2_,mvaDTrainC1_,mvaDTrainC2_;





  bool RunB, RunC, RunD, RunE, RunF = 0;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  TRandom rndm;
  double random_num;
  // Selection Variables
  // Initializhere your analysis histograms


  std::vector<TH1D> pTMu1OverMass_TRF;
  std::vector<TH1D> cTheta_MuonOS_TauPol_TRF;

  std::vector<TH1D> OSSS1Angle_TRF;
  std::vector<TH1D> OSSS2Angle_TRF;

  std::vector<TH1D> OSSS1Angle_RRF;
  std::vector<TH1D> OSSS2Angle_RRF;


  std::vector<TH1D> cTheta_TRF_SSSS;
  std::vector<TH1D> cTheta_TRF_OSSS;



  std::vector<TH1D> EtaGenSource;

  std::vector<TH1D> PairMassDRSorted1A;
  std::vector<TH1D> PairMassDRSorted2A;


  std::vector<TH1D> PairMassDRSorted1B;
  std::vector<TH1D> PairMassDRSorted2B;

  std::vector<TH1D> PairMassDRSorted1C;
  std::vector<TH1D> PairMassDRSorted2C;



  std::vector<TH1D> BetterMuMuVertex;
  std::vector<TH1D> WorseMuMuVertex;




  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;
  std::vector<TH1D> EventMassResolution_PtEtaPhi_TauEta1p2;

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefitA1;



  std::vector<TH1D> TauMassA1;
  std::vector<TH1D> TauMassRefitB1;



  std::vector<TH1D> TauMassB1;
  std::vector<TH1D> TauMassRefitC1;


  std::vector<TH1D> TauMassRefitABC1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitABC2FullEtaVetoCut;

  std::vector<TH1D> TauMassC1;

  std::vector<TH1D> TauMassRefitA2;
  std::vector<TH1D> TauMassA2;
  std::vector<TH1D> TauMassRefitB2;
  std::vector<TH1D> TauMassB2;
  std::vector<TH1D> TauMassRefitC2;
  std::vector<TH1D> TauMassC2;



  std::vector<TH1D> TauMassRefitABC1;
  std::vector<TH1D> TauMassRefitABC2;


  std::vector<TH1D> TauMassRefitABC1_BDSeparateTrain;
  std::vector<TH1D> TauMassRefitABC2_BDSeparateTrain;


  std::vector<TH2D> TauMassRefitABC1_eta;
  std::vector<TH2D> TauMassRefitABC2_eta;



  std::vector<TH1D> TauMassRefitABC1_train11;
  std::vector<TH1D> TauMassRefitABC2_train11;
  std::vector<TH1D> TauMassRefitA1_train11;
  std::vector<TH1D> TauMassRefitB1_train11;
  std::vector<TH1D> TauMassRefitC1_train11;
  std::vector<TH1D> TauMassRefitA2_train11;
  std::vector<TH1D> TauMassRefitB2_train11;
  std::vector<TH1D> TauMassRefitC2_train11;
  std::vector<TH1D> AllignSortMass1_train11;
  std::vector<TH1D> AllignSortMass2_train11;
  std::vector<TH1D> BDTOutputA_train11;
  std::vector<TH1D> BDTOutputB_train11;
  std::vector<TH1D> BDTOutputC_train11;

  std::vector<TH1D> TauMassRefitABC1_train8;
  std::vector<TH1D> TauMassRefitABC2_train8;
  std::vector<TH1D> TauMassRefitA1_train8;
  std::vector<TH1D> TauMassRefitB1_train8;
  std::vector<TH1D> TauMassRefitC1_train8;
  std::vector<TH1D> TauMassRefitA2_train8;
  std::vector<TH1D> TauMassRefitB2_train8;
  std::vector<TH1D> TauMassRefitC2_train8;
  std::vector<TH1D> AllignSortMass1_train8;
  std::vector<TH1D> AllignSortMass2_train8;
  std::vector<TH1D> BDTOutputA_train8;
  std::vector<TH1D> BDTOutputB_train8;
  std::vector<TH1D> BDTOutputC_train8;

  std::vector<TH1D> BDTOutputABC_train2;
  std::vector<TH1D> BDTOutputABC_train8;
  std::vector<TH1D> BDTOutputABC_train11;



  std::vector<TH1D> TauMassResolutionRefit;
  std::vector<TH1D> TauMassResolutionHelixRefit;

  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> TauMass_all;
  std::vector<TH2D> TauMass_allVsBDTA;
  std::vector<TH2D> TauMass_allVsBDTB;
  std::vector<TH2D> TauMass_allVsBDTC;
  std::vector<TH2D> TauMass_allVsBDTBarrel;
  std::vector<TH2D> TauMass_allVsBDTEndcap;
  std::vector<TH2D> EMR_tau_eta;

  std::vector<TH2D> PairMass;

  std::vector<TH1D> KpiIsolationMass_OS;
  std::vector<TH1D> KpiIsolationMass_SS1;
  std::vector<TH1D> KpiIsolationMass_SS2;





  std::vector<TH2D> PairMassdRSorted;
  std::vector<TH2D> PairMassVertexSorted;
  std::vector<TH1D> PairMass1VertexSorting;
  std::vector<TH1D> PairMass2VertexSorting;

  std::vector<TH2D> PairMassdRSortedXVeto;
  std::vector<TH2D> PairMassPhiMassSorting;
  std::vector<TH1D> PairMass1PhiMassSorting;
  std::vector<TH1D> PairMass2PhiMassSorting;

  std::vector<TH2D> PairMass1TauPhiMassSorting;
  std::vector<TH2D> PairMass2TauPhiMassSorting;


  std::vector<TH2D> PairMassFinalSel;
  std::vector<TH1D> PairMass1;
  std::vector<TH1D> PairMass2;

  std::vector<TH2D> PairMassEta;
  std::vector<TH2D> PairMassEtaPrime;

  std::vector<TH1D> AllignSortMass1;
  std::vector<TH1D> AllignSortMass2;


  std::vector<TH1D> PairMass1NoSorting;
  std::vector<TH1D> PairMass2NoSorting;
  std::vector<TH2D> MuMuMassNoSorting;

  std::vector<TH1D> PairMass1PTSorting;
  std::vector<TH1D> PairMass2PTSorting;
  std::vector<TH2D> MuMuMassPTSorting;

  std::vector<TH1D> PairMass1AllignedSorting;
  std::vector<TH1D> PairMass2AllignedSorting;
  std::vector<TH2D> MuMuMassAllignedSorting;

  std::vector<TH2D> PairMassWithCut;




  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;


  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;


  std::vector<TH1D> VertexChi2KF;
  std::vector<TH2D> VertexChi2KF_vs_HelixFit;

  std::vector<TH1D>   KF_Helix_deltaX;
  std::vector<TH1D>   KF_Helix_deltaY;
  std::vector<TH1D>   KF_Helix_deltaZ;

  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> BDTOutputA;
  std::vector<TH1D> BDTOutputB;
  std::vector<TH1D> BDTOutputC;
  std::vector<TH1D> BvsDBDTG;

  std::vector<TH1D> BvsDBDTG_ABC1;
  std::vector<TH1D> BvsDBDTG_ABC2;

  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;


  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;
  std::vector<TH1D>  Muon3MVAID;

  std::vector<TH1D> AncestorsSize;

  std::vector<TH1D> TypeI;
  std::vector<TH1D> TypeI_N_Fakes;
  std::vector<TH1D> TypeI_N_DecInFlight;



  TMVA::Reader *readerA_train11;
  TMVA::Reader *readerB_train11;
  TMVA::Reader *readerC_train11;

  TMVA::Reader *readerA_train8;
  TMVA::Reader *readerB_train8;
  TMVA::Reader *readerC_train8;

  TMVA::Reader *readerA_train12;
  TMVA::Reader *readerB_train12;
  TMVA::Reader *readerC_train12;

  TMVA::Reader *readerA_train13;
  TMVA::Reader *readerB_train13;
  TMVA::Reader *readerC_train13;



  TMVA::Reader *readerA;
  TMVA::Reader *readerB;
  TMVA::Reader *readerC;


  TMVA::Reader *readerMuIDBarrel;
  TMVA::Reader *readerMuIDEndcap;

  TMVA::Reader *readerBvsD;



  TMVA::Reader *readerBTrainA;
  TMVA::Reader *readerBTrainB;
  TMVA::Reader *readerBTrainC;


  TMVA::Reader *readerDTrainA;
  TMVA::Reader *readerDTrainB;
  TMVA::Reader *readerDTrainC;



  Float_t var_vertexKFChi2;// (chi sq of the fit of the secondary vertex)
  Float_t var_svpvTauAngle;// (The angle between PV-SV vector and the tau vector)
  Float_t var_flightLenSig;// (Flight length significance of the tau candidate)
  Float_t var_flightLenDist;
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

  Float_t var_NtracksClose;
  Float_t var_Iso08;
  Float_t var_MaxD0SigBS;
  Float_t var_MinD0SigBS;
  

  float var_MaxtrkKink;
  float var_MaxD0SigSV;
  float var_MindcaTrackSV;
  float var_maxMuonsDca;
  float var_nsv;
  float var_dcaTrackPV;

  float var_MaxMuon_chi2LocalPosition;
  float var_MaxVertexPairQuality;
  float var_MuonglbkinkSum;
  float var_MaxMuon_chi2LocalMomentum;
  float var_MinMuonImpactAngle;
  float var_MaxMuonImpactAngle;
  float var_Vertex2muTrkKF;



  Double_t m3m;
  Double_t dataMCtype;
  Double_t event_weight;
  Double_t bdt;
  Double_t category;
  Double_t m12;
  Double_t m13;
  Double_t mDr1;
  Double_t mDr2;
  Double_t xv;
  Double_t phiv;


  Double_t LumiScale;
  Double_t mvaA1;
  Double_t mvaA2;
  Double_t mvaB1;
  Double_t mvaB2;
  Double_t mvaC1;
  Double_t mvaC2;

  Double_t mvaBTrainA1;
  Double_t mvaBTrainA2;
  Double_t mvaBTrainB1;
  Double_t mvaBTrainB2;
  Double_t mvaBTrainC1;
  Double_t mvaBTrainC2;


  Double_t mvaDTrainA1;
  Double_t mvaDTrainA2;
  Double_t mvaDTrainB1;
  Double_t mvaDTrainB2;
  Double_t mvaDTrainC1;
  Double_t mvaDTrainC2;


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


  Float_t var_IsoPhiKKMass_Mu1;
  Float_t var_IsoKStarMass_Mu1;
  Float_t var_IsoMuMuMass_Mu1;

  Float_t var_IsoPhiKKMass_Mu2;
  Float_t var_IsoKStarMass_Mu2;
  Float_t var_IsoMuMuMass_Mu2;

  Float_t var_IsoPhiKKMass_Mu3;
  Float_t var_IsoKStarMass_Mu3;
  Float_t var_IsoMuMuMass_Mu3;

  float var_BvsDSeprator;


};
#endif
