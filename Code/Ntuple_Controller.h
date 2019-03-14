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

   //   std::vector<std::vector<double> > *Muon_outerTrack_p4;
   //   std::vector<std::vector<double> > *Muon_innerTrack_p4;
   /*   unsigned int> *Muon_Track_idx;
   int>          *Muon_hitPattern_pixelLayerwithMeas;
   int>          *Muon_numberOfMatchedStations;
   float>   *Muon_normChi2;
   int>     *Muon_hitPattern_numberOfValidMuonHits;
   int>     *Muon_innerTrack_numberofValidHits;
   int>     *Muon_numberOfMatches;
   int>     *Muon_numberOfChambers;
   bool>    *Muon_isPFMuon;
   bool>    *Muon_isRPCMuon;
   int>     *Muon_numberofValidPixelHits;
   int>     *Muon_trackerLayersWithMeasurement;
   bool>    *Muon_combinedQuality_updatedSta;
   double>  *Muon_combinedQuality_trkKink;
   double>  *Muon_combinedQuality_glbKink;
   double>  *Muon_combinedQuality_trkRelChi2;
   double>  *Muon_combinedQuality_staRelChi2;
   double>  *Muon_combinedQuality_chi2LocalPosition;
   double>  *Muon_combinedQuality_chi2LocalMomentum;
   double>  *Muon_combinedQuality_localDistance;
   double>  *Muon_combinedQuality_globalDeltaEtaPhi;
   bool>    *Muon_combinedQuality_tightMatch;
   double>  *Muon_combinedQuality_glbTrackProbability;
   double>  *Muon_prod_inner_outer_charge;

   double>  *Muon_innerTrack_quality;
   double>  *Muon_ptErrOverPt;
   double>  *Muon_calEnergy_hadS9;
   double>  *Muon_calEnergy_had;
   double>  *Muon_calEnergy_emS25;
   double>  *Muon_calEnergy_emS9;
   double>  *Muon_calEnergy_em;
   bool>    *Muon_segmentCompatibility;
   bool>    *Muon_caloCompatibility;
   bool>    *Muon_isGoodMuon_TM2DCompatibility;
   double>  *Muon_innerTrack_validFraction;
   double>  *Muon_innerTrack_pixelLayersWithMeasurement;
   double>  *Muon_innerTrack_numberOfValidTrackerHits;
   double>  *Muon_innerTrack_numberOfLostTrackerHits;
   double>  *Muon_innerTrack_numberOfLostTrackerInnerHits;
   double>  *Muon_innerTrack_numberOfLostTrackerOuterHits;
   double>  *Muon_innerTrack_normalizedChi2;
   double>  *Muon_outerTrack_normalizedChi2;
   double>  *Muon_outerTrack_muonStationsWithValidHits;
   bool>    *Muon_isGoodMuon_TrackerMuonArbitrated;
   bool>    *Muon_isGoodMuon_TMOneStationTight;
   bool>    *Muon_isGoodMuon_TMOneStationAngTight;
   bool>    *Muon_isGoodMuon_TMLastStationTight;
   bool>    *Muon_isGoodMuon_TMLastStationAngTight;
   bool>    *Muon_isGoodMuon_TMLastStationOptimizedLowPtTight;
   bool>    *Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;
   double>  *Muon_vmuonhitcomb_reco;
   double>  *Muon_rpchits_reco; */
   int      Muon_charge(unsigned int i){return Ntp->Muon_charge->at(i);}
   /*   int>     *Muon_trackCharge;
   int>     *Muon_pdgid;
   double>  *Muon_B;
   double>  *Muon_M;
   std::vector<std::vector<double> > *Muon_par;
   std::vector<std::vector<double> > *Muon_cov;*/
   





   unsigned int   NThreeMuons(){return Ntp->ThreeMuons_index->size();}
   std::vector<unsigned int> ThreeMuonIndices(unsigned int i){return Ntp->ThreeMuons_index->at(i);}


   unsigned int   NTwoMuonsTrack(){return Ntp->TwoMuonsTrack_Muonsindex->size();}
   std::vector<unsigned int> TwoMuonsTrackMuonIndices(unsigned int i){return Ntp->TwoMuonsTrack_Muonsindex->at(i);}
   std::vector<unsigned int> TwoMuonsTrackTrackIndex(unsigned int i){return Ntp->TwoMuonsTrack_Trackindex->at(i);}



   int      NumberOfSVertices(){return Ntp->Vertex_signal_KF_Chi2->size();}
   double   Vertex_Signal_KF_Chi2(unsigned int i){return Ntp->Vertex_signal_KF_Chi2->at(i);}
   TVector3 Vertex_Signal_KF_pos(unsigned int i){return TVector3(Ntp->Vertex_signal_KF_pos->at(i).at(0), Ntp->Vertex_signal_KF_pos->at(i).at(1),Ntp->Vertex_signal_KF_pos->at(i).at(2))}

   ///// closest distance between the tracks of a candidate
   double Vertex_DCA12(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(0);}
   double Vertex_DCA23(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(1);}
   double Vertex_DCA31(unsigned int i) {return Ntp->Vertex_signal_dca_reco->at(i).at(2);}
   double Vertex_DCAMax(unsigned int i){return Ntp->Vertex_signal_dca_reco->at(i).at(3);}


   /*
   std::vector<std::vector<std::vector<double> > > *Vertex_signal_KF_refittedTracksP4;
   std::vector<double>  *Vertex_signal_KF_Chi2;
   std::vector<std::vector<double> > *Vertex_signal_AF_pos;
   std::vector<double>  *Vertex_signal_AF_Chi2;
   std::vector<double>  *Vertex_signal_AF_Ndf;
   std::vector<std::vector<double> > *Vertex_pair_quality;
   std::vector<std::vector<double> > *Vertex_pairfit_status;
   std::vector<std::vector<double> > *Vertex_MatchedPrimaryVertex;
   std::vector<bool>    *Vertex_RefitPVisValid;
   std::vector<std::vector<double> > *Vertex_MatchedRefitPrimaryVertex;
   std::vector<std::vector<double> > *Vertex_d0_reco;
   std::vector<std::vector<double> > *Vertex_d0sig_reco;
   std::vector<std::vector<double> > *Vertex_2Ddisplacement;
   std::vector<std::vector<double> > *Vertex_3Ddisplacement;
   std::vector<std::vector<float> > *Vertex_Isolation1;
   std::vector<std::vector<float> > *Vertex_Isolation2;
   std::vector<std::vector<float> > *Vertex_Isolation3;
   std::vector<std::vector<float> > *Vertex_Isolation4;
   */


   int NL1Seeds(){return Ntp->Trigger_l1name->size();}
   string L1Name(unsigned int i){return Ntp->Trigger_l1name->at(i);}
   int L1Decision(unsigned int i){return Ntp->Trigger_l1decision->at(i);}
   int Trigger_l1prescale(unsigned int i){return Ntp->Trigger_l1prescale->at(i);}
   int NHLT(){return Ntp->Trigger_hltname->size();}
   string HLTName(unsigned int i){return Ntp->Trigger_hltname->at(i);}
   int HLTDecision(unsigned int i){return Ntp->Trigger_hltdecision->at(i);}


};
#endif
