//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TSystem.h"

// Include files (C & C++ libraries)
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <utility>      // std::pair
#include <tuple>
#include <functional>

#include "NtupleReader.h"

#include "HistoConfig.h"


// small struct needed to allow sorting indices by some value
struct sortIdxByValue {
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
};

///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Vladimir Cherepanov
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  int currentEvent;

  bool cannotObtainHiggsMass; // avoid repeated printing of warning when running locally

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // helper functions for internal calculations
  void printMCDecayChain(unsigned int par, unsigned int level = 0, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false);


  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;

  // Interfaces
  HistoConfig HistoC;

  // Fit Variables


  // muon correction related objects
  //  rochcor2012*   rmcor;
  std::vector<TLorentzVector> Muon_corrected_p4;
  void           CorrectMuonP4();
  bool           Muon_isCorrected;


 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  virtual ~Ntuple_Controller() ;

  // Event initializer
  void InitEvent();

  //TauSpiner function
  double TauSpinerGet(int SpinType);
  void TauSpinerSetSignal(int signalcharge){
#ifdef USE_TauSpinner
TauSpinerInt.SetTauSignalCharge(signalcharge);
#endif
}
   enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
   enum TrackQuality {
     undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
     confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
     qualitySize = 7
   };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};
  enum MuonQualityBitMask{Bit_MuonLoose=0,Bit_MuonSoft,Bit_MuonMedium,Bit_MuonTight, Bit_MuonHighPt, Bit_MuonTight_noVtx};
  enum MuonStandardSelectors{CutBasedIdLoose=0,
    CutBasedIdMedium,
    CutBasedIdMediumPrompt,
    CutBasedIdTight,
    CutBasedIdGlobalHighPt,
    CutBasedIdTrkHighPt,
    PFIsoVeryLoose,
    PFIsoLoose,
    PFIsoMedium,
    PFIsoTight,
    PFIsoVeryTight,
    TkIsoLoose,
    TkIsoTight,
    SoftCutBasedId,
    SoftMvaId,
    MvaLoose,
    MvaMedium,
    MvaTight,
    MiniIsoLoose,
    MiniIsoMedium,
    MiniIsoTight,
    MiniIsoVeryTight};

  // Ntuple Access Functions 
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline  void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}
 
  // Systematic controls
  enum    Systematic {Default=0,NSystematics};

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


  // Data/MC switch and thin
  bool isData()  {return (bool)Ntp->Event_isRealData;}
  void ThinTree();


  bool                                isInit;

  // Information from input Ntuple path name
  TString GetInputNtuplePath();
  TString GetInputDatasetName();
  TString GetInputPublishDataName();


  // Physics Variable Get Functions
  // Event Variables
  Long64_t GetMCID();
  int GetStrippedMCID();
  /* // Vertex Information */

   ULong64_t EventNumber(){return Ntp->Event_EventNumber;}
   Int_t     RunNumber(){return Ntp->Event_RunNumber;}
   Int_t     DataMC_Type(){return Ntp->Event_DataMC_Type;}
   Int_t     LuminosityBlock(){return Ntp->Event_luminosityBlock;}
   float     NVtx(){return Ntp->Vertex_N_primary;}
   double    DeltaPhi(double, double);

   unsigned int   NTracks(){return Ntp->Track_p4->size();}
   TLorentzVector Track_P4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2), Ntp->Track_p4->at(i).at(3), Ntp->Track_p4->at(i).at(0));}
   TVector3       Track_Poca(unsigned int i){return TVector3(Ntp->Track_poca->at(i).at(0),Ntp->Track_poca->at(i).at(1),Ntp->Track_poca->at(i).at(2));}
   double         Track_normalizedChi2(unsigned int i){return Ntp->Track_normalizedChi2->at(i);}
   double         Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);}
   double         Track_charge(unsigned int i){return Ntp->Track_charge->at(i);}
   double         Track_dxy(unsigned int i){return  Ntp->Track_dxy->at(i);}
   double         Track_dz(unsigned int i){return  Ntp->Track_dz->at(i);}
   double         Track_dxyError(unsigned int i){return  Ntp->Track_dxyError->at(i);}
   double         Track_dzError(unsigned int i){return  Ntp->Track_dzError->at(i);}

   unsigned int   NMuons(){return Ntp->Muon_p4->size();}
   TLorentzVector Muon_P4(unsigned int i){return TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2), Ntp->Muon_p4->at(i).at(3), Ntp->Muon_p4->at(i).at(0));}
   TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));}
   bool    Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);} 
   bool    Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);}
   bool    Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
   bool    Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);}
   bool    Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);}
   bool    Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);}
   bool    Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);}
   float   Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);}
   float   Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);}
   float   Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);}
   float   Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);}
   int     Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);}
   int     Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);}
   float   Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);}
   float   Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);}
   float   Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);}
   float   Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);}
   float   Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);}
   float   Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);}
   int     Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);}
   int     Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);}
   float   Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);}
   float   Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);}
   float   Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}
   float   Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}
   float   Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}
   float   Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}
   float   Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}
   float   Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}
   float   Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}
   float   Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}
   float   Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}
   float   Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}
   float   Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}
   float   Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}
   float   Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}
   float   Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}

   TLorentzVector Muon_outerTrack_p4(unsigned int i){return TLorentzVector(Ntp->Muon_outerTrack_p4->at(i).at(1), Ntp->Muon_outerTrack_p4->at(i).at(2), Ntp->Muon_outerTrack_p4->at(i).at(3),Ntp->Muon_outerTrack_p4->at(i).at(0));}

   TLorentzVector Muon_innerTrack_p4(unsigned int i){return TLorentzVector(Ntp->Muon_innerTrack_p4->at(i).at(1), Ntp->Muon_innerTrack_p4->at(i).at(2), Ntp->Muon_innerTrack_p4->at(i).at(3),Ntp->Muon_innerTrack_p4->at(i).at(0));}
   unsigned int Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);}
   int     Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);}
   int     Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);}
   float   Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);}
   int     Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);}
   int     Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);}
   int     Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);}
   int     Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);}
   bool    Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}
   bool    Muon_isRPCMuon(unsigned int i){return Ntp->Muon_isRPCMuon->at(i);}
   int     Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);} 
   int     Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);}
   bool    Muon_combinedQuality_updatedSta(unsigned int i){return Ntp->Muon_combinedQuality_updatedSta->at(i);}
   double  Muon_combinedQuality_trkKink(unsigned int i){return Ntp->Muon_combinedQuality_trkKink->at(i);}
   double  Muon_combinedQuality_glbKink(unsigned int i){return Ntp->Muon_combinedQuality_glbKink->at(i);}
   double  Muon_combinedQuality_trkRelChi2(unsigned int i){return Ntp->Muon_combinedQuality_trkRelChi2->at(i);}
   double  Muon_combinedQuality_staRelChi2(unsigned int i){return Ntp->Muon_combinedQuality_staRelChi2->at(i);}
   double  Muon_combinedQuality_chi2LocalPosition(unsigned int i){return Ntp->Muon_combinedQuality_chi2LocalPosition->at(i);}
   double  Muon_combinedQuality_chi2LocalMomentum(unsigned int i){return Ntp->Muon_combinedQuality_chi2LocalMomentum->at(i);}
   double  Muon_combinedQuality_localDistance(unsigned int i){return Ntp->Muon_combinedQuality_localDistance->at(i);}
   double  Muon_combinedQuality_globalDeltaEtaPhi(unsigned int i){return Ntp->Muon_combinedQuality_globalDeltaEtaPhi->at(i);}
   bool    Muon_combinedQuality_tightMatch(unsigned int i){return Ntp->Muon_combinedQuality_tightMatch->at(i);}
   double  Muon_combinedQuality_glbTrackProbability(unsigned int i){return Ntp->Muon_combinedQuality_glbTrackProbability->at(i);}
   double  Muon_prod_inner_outer_charge(unsigned int i){return Ntp->Muon_prod_inner_outer_charge->at(i);}

   double  Muon_innerTrack_quality(unsigned int i){return Ntp->Muon_innerTrack_quality->at(i);}
   double  Muon_ptErrOverPt(unsigned int i){return Ntp->Muon_ptErrOverPt->at(i);}
   double  Muon_calEnergy_hadS9(unsigned int i){return Ntp->Muon_calEnergy_hadS9->at(i);}
   double  Muon_calEnergy_had(unsigned int i){return Ntp->Muon_calEnergy_had->at(i);}
   double  Muon_calEnergy_emS25(unsigned int i){return Ntp->Muon_calEnergy_emS25->at(i);}
   double  Muon_calEnergy_emS9(unsigned int i){return Ntp->Muon_calEnergy_emS9->at(i);}
   double  Muon_calEnergy_em(unsigned int i){return Ntp->Muon_calEnergy_em->at(i);}
   bool    Muon_segmentCompatibility(unsigned int i){return Ntp->Muon_segmentCompatibility->at(i);}
   bool    Muon_caloCompatibility(unsigned int i){return Ntp->Muon_caloCompatibility->at(i);}
   bool    Muon_isGoodMuon_TM2DCompatibility(unsigned int i){return Ntp->Muon_isGoodMuon_TM2DCompatibility->at(i);}
   double  Muon_innerTrack_validFraction(unsigned int i){return Ntp->Muon_innerTrack_validFraction->at(i);}
   double  Muon_innerTrack_pixelLayersWithMeasurement(unsigned int i){return Ntp->Muon_innerTrack_pixelLayersWithMeasurement->at(i);}
   double  Muon_innerTrack_numberOfValidTrackerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfValidTrackerHits->at(i);}
   double  Muon_innerTrack_numberOfLostTrackerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerHits->at(i);}
   double  Muon_innerTrack_numberOfLostTrackerInnerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits->at(i);}
   double  Muon_innerTrack_numberOfLostTrackerOuterHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits->at(i);}
   double  Muon_innerTrack_normalizedChi2(unsigned int i){return Ntp->Muon_innerTrack_normalizedChi2->at(i);}
   double  Muon_outerTrack_normalizedChi2(unsigned int i){return Ntp->Muon_outerTrack_normalizedChi2->at(i);}
   double  Muon_outerTrack_muonStationsWithValidHits(unsigned int i){return Ntp->Muon_outerTrack_muonStationsWithValidHits->at(i);}
   bool    Muon_isGoodMuon_TrackerMuonArbitrated(unsigned int i){return Ntp->Muon_isGoodMuon_TrackerMuonArbitrated->at(i);}
   bool    Muon_isGoodMuon_TMOneStationTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMOneStationTight->at(i);}
   bool    Muon_isGoodMuon_TMOneStationAngTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMOneStationAngTight->at(i);}
   bool    Muon_isGoodMuon_TMLastStationTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationTight->at(i);}
   bool    Muon_isGoodMuon_TMLastStationAngTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationAngTight->at(i);}
   bool    Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight->at(i);}
   bool    Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight->at(i);}
   double  Muon_vmuonhitcomb_reco(unsigned int i){return Ntp->Muon_vmuonhitcomb_reco->at(i);}
   double  Muon_rpchits_reco(unsigned int i){return Ntp->Muon_rpchits_reco->at(i);}
   int     Muon_charge(unsigned int i){return Ntp->Muon_charge->at(i);}
   int     Muon_trackCharge(unsigned int i){return Ntp->Muon_trackCharge->at(i);}
   int     Muon_pdgid(unsigned int i){return Ntp->Muon_pdgid->at(i);}
   double  Muon_B(unsigned int i){return Ntp->Muon_B->at(i);}
   double  Muon_M(unsigned int i){return Ntp->Muon_M->at(i);}
   std::vector<unsigned int> SortedPtMuons(std::vector<unsigned int> indixes);
   /*    will be fixed later
   bool Muon_TrackParticleHasMomentum(unsigned int i){if(Ntp->Muon_par->at(i).size()!=0)return true; return false;} 
   TrackParticle Muon_TrackParticle(unsigned int i){ 
     TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1); 
     TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar); 
     unsigned int l=0;
     for(int k=0; k<TrackParticle::NHelixPar; k++){ 
       mu_par(k,0)=Ntp->Muon_par->at(i).at(k); 
       for(int j=k; j<TrackParticle::NHelixPar; j++){ 
         mu_cov(k,j)=Ntp->Muon_cov->at(i).at(l); 
         l++; 
       } 
     } 
     return TrackParticle(mu_par,mu_cov,Ntp->Muon_pdgid->at(i),Ntp->Muon_M->at(i),Ntp->Muon_trackCharge->at(i),Ntp->Muon_B->at(i));
   }  
   */

   unsigned int   NThreeMuons(){return Ntp->ThreeMuons_index->size();}
   std::vector<unsigned int> ThreeMuonIndices(unsigned int i){return Ntp->ThreeMuons_index->at(i);}

   unsigned int   NTwoMuonsTrack(){return Ntp->TwoMuonsTrack_Muonsindex->size();}
   std::vector<unsigned int> TwoMuonsTrackMuonIndices(unsigned int i){return Ntp->TwoMuonsTrack_Muonsindex->at(i);}
   std::vector<unsigned int> TwoMuonsTrackTrackIndex(unsigned int i){return Ntp->TwoMuonsTrack_Trackindex->at(i);}

   int        NumberOfSVertices(){return Ntp->Vertex_signal_KF_Chi2->size();}// should coincide with number of candidate
   double     Vertex_Signal_KF_Chi2(unsigned int i){return Ntp->Vertex_signal_KF_Chi2->at(i);}
   TVector3   Vertex_Signal_KF_pos(unsigned int i){return TVector3(Ntp->Vertex_signal_KF_pos->at(i).at(0), Ntp->Vertex_signal_KF_pos->at(i).at(1),Ntp->Vertex_signal_KF_pos->at(i).at(2));}

   ///// closest distance between the tracks of a candidate
   double     Vertex_DCA12(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(0);}
   double     Vertex_DCA23(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(1);}
   double     Vertex_DCA31(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(2);}
   double     Vertex_DCAMax(unsigned int i){return Ntp->Vertex_signal_dca_reco->at(i).at(3);}

   TLorentzVector Vertex_signal_KF_refittedTracksP4(unsigned int i, unsigned int j){
     return TLorentzVector(Ntp->Vertex_signal_KF_refittedTracksP4->at(i).at(j).at(1),Ntp->Vertex_signal_KF_refittedTracksP4->at(i).at(j).at(2),Ntp->Vertex_signal_KF_refittedTracksP4->at(i).at(j).at(3),Ntp->Vertex_signal_KF_refittedTracksP4->at(i).at(j).at(0));}


   double     Vertex_d0_reco(unsigned int i, unsigned int j){return Ntp->Vertex_d0_reco->at(i).at(j);}//  j - is a track number; j=0,1,2
   double     Vertex_d0sig_reco(unsigned int i, unsigned int j){return Ntp->Vertex_d0sig_reco->at(i).at(j);}//  j - is a track number; j=0,1,2
   double     Vertex_2Ddisplacement(unsigned int i, unsigned int j){return Ntp->Vertex_2Ddisplacement->at(i).at(j);}// j=0 - value, j=1 - significane, j=3 - distance
   double     Vertex_3Ddisplacement(unsigned int i, unsigned int j){return Ntp->Vertex_3Ddisplacement->at(i).at(j);}// j=0 - value, j=1 - significane, j=3 - distance

   double     Vertex_pair_quality(unsigned int i, unsigned int j){return Ntp->Vertex_pair_quality->at(i).at(j);}//the second index should not exceed 3
   bool       Vertex_pairfit_status(unsigned int i, unsigned int j){return Ntp->Vertex_pairfit_status->at(i).at(j);}//the second index should not exceed 3
   double     Vertex_signal_KF_Chi2(unsigned int i){return Ntp->Vertex_signal_KF_Chi2->at(i); }
   TVector3   Vertex_signal_AF_pos(unsigned int i){return TVector3(Ntp->Vertex_signal_AF_pos->at(i).at(0),Ntp->Vertex_signal_AF_pos->at(i).at(1),Ntp->Vertex_signal_AF_pos->at(i).at(2));}
   double     Vertex_signal_AF_Chi2(unsigned int i){return Ntp->Vertex_signal_AF_Chi2->at(i); }
   double     Vertex_signal_AF_Ndf(unsigned int i){return Ntp->Vertex_signal_AF_Ndf->at(i); }
   TVector3   Vertex_MatchedPrimaryVertex(unsigned int i){return TVector3(Ntp->Vertex_MatchedPrimaryVertex->at(i).at(0), Ntp->Vertex_MatchedPrimaryVertex->at(i).at(1),Ntp->Vertex_MatchedPrimaryVertex->at(i).at(2));}
   bool       Vertex_RefitPVisValid(unsigned int i){return Ntp->Vertex_RefitPVisValid->at(i);}
   TVector3   Vertex_MatchedRefitPrimaryVertex(unsigned int i){return TVector3(Ntp->Vertex_MatchedRefitPrimaryVertex->at(i).at(0), Ntp->Vertex_MatchedRefitPrimaryVertex->at(i).at(1),Ntp->Vertex_MatchedRefitPrimaryVertex->at(i).at(2));}

   float      Isolation_RelPt(unsigned int i){return Ntp->Vertex_Isolation1->at(i).at(0);}
   float      Isolation_NTracks(unsigned int i){return Ntp->Vertex_Isolation1->at(i).at(1);}
   float      Isolation_MinDist(unsigned int i){return Ntp->Vertex_Isolation1->at(i).at(2);}

   float      Isolation05_RelPt(unsigned int i){return Ntp->Vertex_Isolation2->at(i).at(0);}
   float      Isolation05_NTracks(unsigned int i){return Ntp->Vertex_Isolation2->at(i).at(1);}
   float      Isolation05_MinDist(unsigned int i){return Ntp->Vertex_Isolation2->at(i).at(2);}

   float      Isolation_Ntrk1(unsigned int i){return Ntp->Vertex_Isolation3->at(i).at(0);}
   float      Isolation_Ntrk2(unsigned int i){return Ntp->Vertex_Isolation3->at(i).at(1);}
   float      Isolation_Ntrk3(unsigned int i){return Ntp->Vertex_Isolation3->at(i).at(2);}

   float      Isolation_Ntrk0p1(unsigned int i){return Ntp->Vertex_Isolation4->at(i).at(0);}
   float      Isolation_Ntrk0p2(unsigned int i){return Ntp->Vertex_Isolation4->at(i).at(1);}
   float      Isolation_Ntrk0p5(unsigned int i){return Ntp->Vertex_Isolation4->at(i).at(2);}
   float      Isolation_maxdy(unsigned int i){return Ntp->Vertex_Isolation4->at(i).at(3);}


   int        NL1Seeds(){return Ntp->Trigger_l1name->size();}
   string     L1Name(unsigned int i){return Ntp->Trigger_l1name->at(i);}
   int        L1Decision(unsigned int i){return Ntp->Trigger_l1decision->at(i);}
   int        Trigger_l1prescale(unsigned int i){return Ntp->Trigger_l1prescale->at(i);}
   int        NHLT(){return Ntp->Trigger_hltname->size();}
   string     HLTName(unsigned int i){return Ntp->Trigger_hltname->at(i);}
   int        HLTDecision(unsigned int i){return Ntp->Trigger_hltdecision->at(i);}

   bool CHECK_BIT(int var, int pos){  return ((var & (1 << pos)) == (1 << pos)); }

};
#endif
