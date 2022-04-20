#ifndef SignalVertexSelector_h
#define SignalVertexSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TError.h"
#include "PDGInfo.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "TDecompBK.h"
#include <TMatrixT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
#include "TMatrixDSymEigen.h"




class SignalVertexSelector : public Selection {

 public:
  SignalVertexSelector(TString Name_, TString id_);
  virtual ~SignalVertexSelector();

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
  
  std::vector<TH1D>  WhetherdRMatch;
  std::vector<TH1D>  WhetherTau3Mu;
  
  std::vector<TH1D>  IsoTrackToMCdR01;
  std::vector<TH1D>  IsoTrackToMCdR08;
  std::vector<TH1D>  IsoTrackToMCAngle01;
  
  std::vector<TH1D>  NumberOfFS_ChargedParticles;
  std::vector<TH1D>  NumberOfFS_ChargedParticles_RecoMatch;
  
  std::vector<TH1D>  NumberOfRecoChargedParticlesIfMC1;
  std::vector<TH1D>  NumberOfRecoChargedParticlesIfMC2;
  std::vector<TH1D>  NumberOfRecoChargedParticlesIfMC3;
  
  std::vector<TH1D>  TwoProngInvariantMassReco;
  std::vector<TH1D>  TwoProngInvariantMassMC;
  
  std::vector<TH1D>  ThreeProngInvariantMassReco;
  std::vector<TH1D>  ThreeProngInvariantMassMC;
  
  std::vector<TH1D>  TwoProngInvariantMassReco005;
  std::vector<TH1D>  TwoProngInvariantMassMC005;
  
  std::vector<TH1D>  ThreeProngInvariantMassReco005;
  std::vector<TH1D>  ThreeProngInvariantMassMC005;
  
  std::vector<TH1D>  TwoProngTrackPt;
  std::vector<TH1D>  TwoProngTrack2Pt;
  
  std::vector<TH1D>  ThreeProngTrackPt;
  std::vector<TH1D>  ThreeProngTrack2Pt;
  std::vector<TH1D>  ThreeProngTrack3Pt;
  
  std::vector<TH1D>  TrackToTauDr2Prong;
  std::vector<TH1D>  TrackToTauDr3Prong;
  std::vector<TH1D>  TrackToTauDrAll;
  
  std::vector<TH1D>  TwoProngTrackEta;
  std::vector<TH1D>  TwoProngTrack2Eta;
  
  std::vector<TH2D>  dR_vs_dP;
  
  std::vector<TH2D>  InvMass2_vs_pdgid;
  std::vector<TH2D>  InvMass3_vs_pdgid;
  
  std::vector<TH2D>  dRmin_sum_vs_InvariantMass_2prong;
  std::vector<TH2D>  dRmin_sum_vs_InvariantMass_3prong;
  
  std::vector<TH1D>  NoOfIsoTracks2Prong;
  std::vector<TH1D>  NoOfIsoTracks3Prong;
  
  std::vector<TH1D>  RankMatchedTrackpT;
  std::vector<TH1D>  RankMatchedTrackdR;
  std::vector<TH1D>  RankMatchedTrackdR_trim;
  
  std::vector<TH1D>  RankMatchedTrackAvgDiff;
  std::vector<TH1D>  RankMatchedTrackCombn;
  
  std::vector<TH1D>  var_All_7_Iso_dR;
  std::vector<TH1D>  var_Correct_Iso_dR;
  
  std::vector<TH1D>  var_All_7_Iso_AvgDiff;
  std::vector<TH1D>  var_Correct_Iso_AvgDiff;
  
  std::vector<TH1D>  var_All_7_Iso_Avg;
  std::vector<TH1D>  var_Correct_Iso_Avg;
  
  std::vector<TH1D>  var_All_7_Iso_Combn;
  std::vector<TH1D>  var_Correct_Iso_Combn;
  
  std::vector<TH1D>  var_All_7_Iso_Trial;
  std::vector<TH1D>  var_Correct_Iso_Trial;
  
  std::vector<TH1D>  RankMatchedTrackdR_cut;
  
  std::vector<TH1D>  RankMatchedTrackPairdR;
  std::vector<TH1D>  TrackPairdR_Crt;
  std::vector<TH1D>  TrackPairdR_Bkg;
  
  std::vector<TH1D>  IsoTrackMatchedToSV_1;
  std::vector<TH1D>  IsoTrackMatchedToSV_TwoMatched;
  std::vector<TH1D>  IsoTrackMatchedToSV_TwoCrtIso;
  std::vector<TH1D>  IsoTrackMatchedToSV_MassMatch;
  std::vector<TH1D>  IsoTrackMatchedToSV_MassMatch1;
  std::vector<TH1D>  IsoTrackMatchedToSV_CombMatch;
  
  std::vector<TH1D>  IsoTrackMatchedToSV_Count;
  
  std::vector<TH1D>  CombMatch_Avg1;
  std::vector<TH1D>  CombMatch_Avg2;
  
  std::vector<TH1D>  Angle_SVPV_iSVSV;
  std::vector<TH1D>  Angle_SVPV_isvSV;
  
  std::vector<TH1D>  IsoTrackMatchedToSV_ThreeMassMatch;
  std::vector<TH1D>  IsoTrackMatchedToSV_ThreeMassMatch1;
  
  std::vector<TH1D>  iSVSV_Distance;
  std::vector<TH1D>  iSVSV_Distance_Sig;
  
  std::vector<TH1D>  isvSV_Distance;
  std::vector<TH1D>  isvSV_Distance_Sig;
  
  std::vector<TH1D>  InvMass2ProngMatched;
  std::vector<TH1D>  InvMass2ProngMatchedSV;
  std::vector<TH1D>  InvMass2ProngNotMatched;
  
  std::vector<TH1D>  SVSize;
  std::vector<TH1D>  SVNoOfTracksMatched;
  std::vector<TH1D>  SVNoOfTracksMatchedThree;
  std::vector<TH1D>  SVNoOfTracksUnmatched;
  
  std::vector<TH1D>  InvMass3ProngMatched;
  std::vector<TH1D>  InvMass3ProngMatchedSV;
  std::vector<TH1D>  InvMass3ProngNotMatched;
  
  std::vector<TH1D>  InvMassTotal;
  std::vector<TH1D>  InvMassTotal1;
  std::vector<TH1D>  InvMassTotal2;
  
  std::vector<TH1D>  SVCollectionNoOfSignalMu;
  std::vector<TH1D>  SVCollectionNoOfNeither;
  
  std::vector<TH1D>  SVCollectionNoOfSignalMu_if1;
  std::vector<TH1D>  SVCollectionNoOfNeither_if1;
  
  std::vector<TH1D>  SVCollectionNoOfSignalMu_ifmore1;
  std::vector<TH1D>  SVCollectionNoOfNeither_ifmore1;
  
  std::vector<TH1D>  SVCollectionNoOfCrt;
  std::vector<TH1D>  SVCollectionNoOfCrt_if3;
  
  std::vector<TH1D>  Whether_Lowest_Chi2_is_Correct_2iso;
  
  std::vector<TH1D>  Rank_Correct_1iso3mu_Chi2;
  std::vector<TH1D>  Rank_Correct_2iso_Chi2;
  std::vector<TH1D>  Rank_Correct_2iso_PairAngle;
  std::vector<TH1D>  Rank_Correct_2iso3mu_Chi2;
  std::vector<TH1D>  Rank_Correct_3iso_Chi2;
  std::vector<TH1D>  Rank_Correct_3iso3mu_Chi2;
  
  std::vector<TH2D> TauEnergyVsPairdR;
  std::vector<TH2D> EnergyVsPairdR;
  std::vector<TH2D> EnergyVsPairdR_1;
  std::vector<TH2D> EnergyVsPairdR_2;
  
  std::vector<TH2D> EnergyVsPairdR_Circular;
  std::vector<TH2D> EnergyVsPairdR_Tau;
  std::vector<TH2D> EnergyVsPairdR_Circular_reco;
  std::vector<TH2D> EnergyVsPairdR_Tau_reco;
  std::vector<TH2D> EnergyVsPairdR_Circular_dR;
  std::vector<TH2D> EnergyVsPairdR_Circular_Incorrect;
  std::vector<TH2D> EnergyVsPairdR_Circular_Incorrect_reco;
  
  
  std::vector<TH1D>  PairAngle_2Cor;
  std::vector<TH1D>  PairAngle_1Cor;
  std::vector<TH1D>  PairAngle_0Cor;
  
  std::vector<TH1D>  PairMass_2Cor;
  std::vector<TH1D>  PairMass_1Cor;
  std::vector<TH1D>  PairMass_0Cor;
  
  std::vector<TH2D> GammaVs_Tau_Energy;
  std::vector<TH2D> GammaVs_Pair_Energy;
  
  std::vector<TH1D>  Angle_SVPV_BMeson;
  std::vector<TH1D>  Angle_Reco_Comparison;
  std::vector<TH1D>  Angle_Tau_Vtx;
  
  std::vector<TH1D>  Angle_dVtx_Pair_Crt;
  std::vector<TH1D>  Angle_dVtx_Pair_InCrt;
  
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
