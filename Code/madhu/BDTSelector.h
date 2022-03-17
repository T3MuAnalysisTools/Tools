#ifndef BDTSelector_h
#define BDTSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"




class BDTSelector : public Selection {

 public:
  BDTSelector(TString Name_, TString id_);
  virtual ~BDTSelector();

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

  double mvaBTrainA1_,mvaBTrainA2_,mvaBTrainB1_,mvaBTrainB2_,mvaBTrainC1_,mvaBTrainC2_;
  double mvaDTrainA1_,mvaDTrainA2_,mvaDTrainB1_,mvaDTrainB2_,mvaDTrainC1_,mvaDTrainC2_;





  bool RunB, RunC, RunD, RunE, RunF = 0;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  TRandom rndm;
  double random_num;
  // Selection Variables
  // Initializhere your analysis histograms

  std::vector<TH1D> PairMassDRSorted1A;
  std::vector<TH1D> PairMassDRSorted2A;


  std::vector<TH1D> PairMassDRSorted1B;
  std::vector<TH1D> PairMassDRSorted2B;

  std::vector<TH1D> PairMassDRSorted1C;
  std::vector<TH1D> PairMassDRSorted2C;


  std::vector<TH1D> Muon1Pt;
  std::vector<TH1D> Muon2Pt;
  std::vector<TH1D> Muon3Pt;

  std::vector<TH1D> Muon1Eta;
  std::vector<TH1D> Muon2Eta;
  std::vector<TH1D> Muon3Eta;


  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauP;

  std::vector<TH1D> BetterMuMuVertex;
  std::vector<TH1D> WorseMuMuVertex;
  
  
  
  std::vector<TH1D> TauAngleTest;
  
  
  
  std::vector<TH1D> MuOSSS1InvariantMassBeforeMVA;
  std::vector<TH1D> MuOSSS2InvariantMassBeforeMVA;
  
  
  
  std::vector<TH1D> IsolationTrackCount;
  
  
  
  
  
  
  
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
  std::vector<TH2D> KKMass_dR_sort;
  std::vector<TH1D> KKMass_dR_sort1;
  std::vector<TH1D> KKMass_dR_sort2;

  std::vector<TH2D> KKMass_pt_sort;
  std::vector<TH1D> KKMass_pt_sort1;
  std::vector<TH1D> KKMass_pt_sort2;


  std::vector<TH2D> KKMass_dR_sort_XVeto;
  std::vector<TH1D> KKMass_dR_sort1_XVeto;
  std::vector<TH1D> KKMass_dR_sort2_XVeto;

  std::vector<TH2D> KKMass_pt_sort_XVeto;
  std::vector<TH1D> KKMass_pt_sort1_XVeto;
  std::vector<TH1D> KKMass_pt_sort2_XVeto;

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
  std::vector<TH2D> VertexChi2KF_vs_HelixFit;

  std::vector<TH1D>   KF_Helix_deltaX;
  std::vector<TH1D>   KF_Helix_deltaY;
  std::vector<TH1D>   KF_Helix_deltaZ;

  std::vector<TH1D> FLSignificance;
  
  std::vector<TH1D> BDTOutputATrain0;
  std::vector<TH1D> BDTOutputBTrain0;
  std::vector<TH1D> BDTOutputCTrain0;
  
  std::vector<TH1D> BDTOutputATrain1;
  std::vector<TH1D> BDTOutputBTrain1;
  std::vector<TH1D> BDTOutputCTrain1;
  
  std::vector<TH1D> BDTOutputATrain2;
  std::vector<TH1D> BDTOutputBTrain2;
  std::vector<TH1D> BDTOutputCTrain2;
  
  std::vector<TH1D> BDTOutputATrain3;
  std::vector<TH1D> BDTOutputBTrain3;
  std::vector<TH1D> BDTOutputCTrain3;
  
  std::vector<TH1D> BDTOutputATrain4;
  std::vector<TH1D> BDTOutputBTrain4;
  std::vector<TH1D> BDTOutputCTrain4;
  
  std::vector<TH1D> BDTOutputATrain5;
  std::vector<TH1D> BDTOutputBTrain5;
  std::vector<TH1D> BDTOutputCTrain5;
  
  std::vector<TH1D> BDTOutputATrain6;
  std::vector<TH1D> BDTOutputBTrain6;
  std::vector<TH1D> BDTOutputCTrain6;
  
  std::vector<TH1D> BDTOutputATrain7;
  std::vector<TH1D> BDTOutputBTrain7;
  std::vector<TH1D> BDTOutputCTrain7;
  
  std::vector<TH1D> BDTOutputATrain8;
  std::vector<TH1D> BDTOutputBTrain8;
  std::vector<TH1D> BDTOutputCTrain8;
  
  std::vector<TH1D> BDTOutputATrain9;
  std::vector<TH1D> BDTOutputBTrain9;
  std::vector<TH1D> BDTOutputCTrain9;
  
  std::vector<TH1D> BDTOutputATrain10;
  std::vector<TH1D> BDTOutputBTrain10;
  std::vector<TH1D> BDTOutputCTrain10;
  
  std::vector<TH1D> BDTOutputATrain11;
  std::vector<TH1D> BDTOutputBTrain11;
  std::vector<TH1D> BDTOutputCTrain11;
  
  std::vector<TH1D> BDTOutputATrain12;
  std::vector<TH1D> BDTOutputBTrain12;
  std::vector<TH1D> BDTOutputCTrain12;
  
  std::vector<TH1D> BDTOutputATrain13;
  std::vector<TH1D> BDTOutputBTrain13;
  std::vector<TH1D> BDTOutputCTrain13;
  
  std::vector<TH1D> BDTOutputATrain14;
  std::vector<TH1D> BDTOutputBTrain14;
  std::vector<TH1D> BDTOutputCTrain14;
  
  std::vector<TH1D> BDTOutputATrain15;
  std::vector<TH1D> BDTOutputBTrain15;
  std::vector<TH1D> BDTOutputCTrain15;
  
  std::vector<TH1D> BDTOutputATrain16;
  std::vector<TH1D> BDTOutputBTrain16;
  std::vector<TH1D> BDTOutputCTrain16;
  
  std::vector<TH1D> BDTOutputATrain17;
  std::vector<TH1D> BDTOutputBTrain17;
  std::vector<TH1D> BDTOutputCTrain17;
  
  std::vector<TH1D> BDTOutputATrain18;
  std::vector<TH1D> BDTOutputBTrain18;
  std::vector<TH1D> BDTOutputCTrain18;
  
  std::vector<TH1D> BDTOutputATrain19;
  std::vector<TH1D> BDTOutputBTrain19;
  std::vector<TH1D> BDTOutputCTrain19;
  
  std::vector<TH1D> BvsDBDTG;

  std::vector<TH1D> BvsDBDTG_ABC1;
  std::vector<TH1D> BvsDBDTG_ABC2;

  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;


  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;
  std::vector<TH1D>  Muon3MVAID;
  
  std::vector<TH2D> BDTOutputCTrain0_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain1_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain2_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain3_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain4_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain5_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain6_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain7_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain8_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain9_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain10_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain11_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain12_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain13_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain14_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain15_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain16_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain17_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain18_Vs_TauMass;
  std::vector<TH2D> BDTOutputCTrain19_Vs_TauMass;



  TMVA::Reader *readerATrain0;
  TMVA::Reader *readerBTrain0;
  TMVA::Reader *readerCTrain0;
  
  TMVA::Reader *readerATrain1;
  TMVA::Reader *readerBTrain1;
  TMVA::Reader *readerCTrain1;
  
  TMVA::Reader *readerATrain2;
  TMVA::Reader *readerBTrain2;
  TMVA::Reader *readerCTrain2;
  
  TMVA::Reader *readerATrain3;
  TMVA::Reader *readerBTrain3;
  TMVA::Reader *readerCTrain3;
  
  TMVA::Reader *readerATrain4;
  TMVA::Reader *readerBTrain4;
  TMVA::Reader *readerCTrain4;
  
  TMVA::Reader *readerATrain5;
  TMVA::Reader *readerBTrain5;
  TMVA::Reader *readerCTrain5;
  
  TMVA::Reader *readerATrain6;
  TMVA::Reader *readerBTrain6;
  TMVA::Reader *readerCTrain6;
  
  TMVA::Reader *readerATrain7;
  TMVA::Reader *readerBTrain7;
  TMVA::Reader *readerCTrain7;
  
  TMVA::Reader *readerATrain8;
  TMVA::Reader *readerBTrain8;
  TMVA::Reader *readerCTrain8;
  
  TMVA::Reader *readerATrain9;
  TMVA::Reader *readerBTrain9;
  TMVA::Reader *readerCTrain9;
  
  TMVA::Reader *readerATrain10;
  TMVA::Reader *readerBTrain10;
  TMVA::Reader *readerCTrain10;
  
  TMVA::Reader *readerATrain11;
  TMVA::Reader *readerBTrain11;
  TMVA::Reader *readerCTrain11;
  
  TMVA::Reader *readerATrain12;
  TMVA::Reader *readerBTrain12;
  TMVA::Reader *readerCTrain12;
  
  TMVA::Reader *readerATrain13;
  TMVA::Reader *readerBTrain13;
  TMVA::Reader *readerCTrain13;
  
  TMVA::Reader *readerATrain14;
  TMVA::Reader *readerBTrain14;
  TMVA::Reader *readerCTrain14;
  
  TMVA::Reader *readerATrain15;
  TMVA::Reader *readerBTrain15;
  TMVA::Reader *readerCTrain15;
  
  TMVA::Reader *readerATrain16;
  TMVA::Reader *readerBTrain16;
  TMVA::Reader *readerCTrain16;
  
  TMVA::Reader *readerATrain17;
  TMVA::Reader *readerBTrain17;
  TMVA::Reader *readerCTrain17;
  
  TMVA::Reader *readerATrain18;
  TMVA::Reader *readerBTrain18;
  TMVA::Reader *readerCTrain18;
  
  TMVA::Reader *readerATrain19;
  TMVA::Reader *readerBTrain19;
  TMVA::Reader *readerCTrain19;


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
  
  float var_BvsDSeprator;
  float var_Vertex2muTrkKF;

  float var_MaxMuon_chi2LocalPosition;
  float var_MaxVertexPairQuality;
  float var_MuonglbkinkSum;
  float var_MaxMuon_chi2LocalMomentum;
  float var_MinMuonImpactAngle;
  float var_MaxMuonImpactAngle;
  
  float var_IsoPhiKKMass_Mu1;
  float var_IsoKStarMass_Mu1;
  float var_IsoMuMuMass_Mu1;

  float var_IsoPhiKKMass_Mu2;
  float var_IsoKStarMass_Mu2;
  float var_IsoMuMuMass_Mu2;

  float var_IsoPhiKKMass_Mu3;
  float var_IsoKStarMass_Mu3;
  float var_IsoMuMuMass_Mu3;

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


};
#endif
