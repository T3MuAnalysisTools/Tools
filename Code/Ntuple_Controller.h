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

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // TString flags for object corrections
  TString tauCorrection;
  TString muonCorrection;
  TString elecCorrection;
  TString jetCorrection;

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
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

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
   Int_t RunNumber(){return Ntp->Event_RunNumber;}
   Int_t           LuminosityBlock(){return Ntp->Event_luminosityBlock;}
   //   Long64_t        triggerbit(){return Ntp->triggerbit;}
   //   Long64_t           metfilterbit(){return Ntp->metfilterbit;}

   //   Float_t         npu(){return Ntp->npu;}
   //   Int_t           npv(){return Ntp->npv;}
   //   Float_t         PUReweight(){return Ntp->PUReweight;}
   //   Int_t           PUNumInteractions(){return Ntp->PUNumInteractions;}
   //   Float_t         rho(){return Ntp->rho;}



   float NVtx(){return Ntp->Vertex_N_primary;}


   double DeltaPhi(double, double);
   /*
   int NTriggers(){return Ntp->trigger_accept->size();}
   bool TriggerAccept(unsigned int i){return Ntp->trigger_accept->at(i);}
   TString TriggerName(unsigned int i){return Ntp->trigger_name->at(i);}
   bool         GetTriggerIndex(TString n,  int &i);
   std::vector<int> GetVectorTriggers(TString n); 
   std::vector<int> GetVectorTriggers(std::vector<TString> v);
   std::vector<int> GetVectorTriggersFullMatch(std::vector<TString> v);
   std::vector<int> GetVectorCrossTriggers(TString n1,TString n2,TString f1,TString f2);
   bool  CheckIfAnyPassed(  std::vector<int> list);
 
   Int_t            DataMC_Type(){return Ntp->DataMC_Type_idx;}

   unsigned int    NGenParts(){return Ntp->genpart_px->size();}
   TLorentzVector Genpart_P4(unsigned int i){return TLorentzVector(Ntp->genpart_px->at(i), Ntp->genpart_py->at(i), Ntp->genpart_pz->at(i),Ntp->genpart_e->at(i));}
   int Genpart_pdg(unsigned int i){return Ntp->genpart_pdg->at(i);}
   int Genpart_status(unsigned int i){return Ntp->genpart_status->at(i);}
   int Genpart_HMothInd(unsigned int i){return Ntp->genpart_HMothInd->at(i);}
   int Genpart_MSSMHMothInd(unsigned int i){return Ntp->genpart_MSSMHMothInd->at(i);}
   int Genpart_TopMothInd(unsigned int i){return Ntp->genpart_TopMothInd->at(i);}
   int Genpart_TauMothInd(unsigned int i){return Ntp->genpart_TauMothInd->at(i);}
   int Genpart_ZMothInd(unsigned int i){return Ntp->genpart_ZMothInd->at(i);}
   int Genpart_WMothInd(unsigned int i){return Ntp->genpart_WMothInd->at(i);}
   int Genpart_bMothInd(unsigned int i){return Ntp->genpart_bMothInd->at(i);}
   int Genpart_HZDecayMode(unsigned int i){return Ntp->genpart_HZDecayMode->at(i);}
   int Genpart_TopDecayMode(unsigned int i){return Ntp->genpart_TopDecayMode->at(i);}
   int Genpart_WDecayMode(unsigned int i){return Ntp->genpart_WDecayMode->at(i);}
   int Genpart_TauGenDecayMode(unsigned int i){return Ntp->genpart_TauGenDecayMode->at(i);}
   int Genpart_TauGenDetailedDecayMode(unsigned int i){return Ntp->genpart_TauGenDetailedDecayMode->at(i);}
   int Genpart_flags(unsigned int i){return Ntp->genpart_flags->at(i);}



 float  dxy(unsigned int i){return Ntp->dxy->at(i);}
 float  dz(unsigned int i){return Ntp->dz->at(i);}
 float  dxy_innerTrack(unsigned int i){return Ntp->dxy_innerTrack->at(i);}
 float  dz_innerTrack(unsigned int i){return Ntp->dz_innerTrack->at(i);}
 float  Daughters_rel_error_trackpt(unsigned int i){return Ntp->daughters_rel_error_trackpt->at(i);}
 TLorentzVector             MCParticle_p4(unsigned int i){return TLorentzVector(Ntp->MC_p4->at(i).at(1),Ntp->MC_p4->at(i).at(2),Ntp->MC_p4->at(i).at(3),Ntp->MC_p4->at(i).at(0));}
 int                        MCParticle_pdgid(unsigned int i){return Ntp->MC_pdgid->at(i);}
 int                        MCParticle_charge(unsigned int i){return Ntp->MC_charge->at(i);}
 //int              		  MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
 int              		  MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
 std::vector<int>           MCParticle_childpdgid(unsigned int i){return Ntp->MC_childpdgid->at(i);}
 std::vector<int>           MCParticle_childidx(unsigned int i){return Ntp->MC_childidx->at(i);}
 int						  MCParticle_status(unsigned int i){return Ntp->MC_status->at(i);}
 int 						  getMatchTruthIndex(TLorentzVector tvector, int pid, double dr);
 int						  matchTruth(TLorentzVector tvector);
 bool						  matchTruth(TLorentzVector tvector, int pid, double dr);
 // decay tree functionality
 bool						  MCParticle_hasMother(unsigned int i){return Ntp->MC_midx->at(i) >= 0;}
 void						  printMCDecayChainOfMother(unsigned int i, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // decay chain of object i
 void						  printMCDecayChainOfEvent(bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // full event decay chain
 std::string				  MCParticleToString(unsigned int par, bool printStatus = false, bool printPt = false, bool printEtaPhi = false);
 bool CheckDecayID( unsigned int jak1,unsigned  int jak2);
 TLorentzVector GetTruthTauLV(unsigned int jak, unsigned  int number);
 TLorentzVector GetTruthTauProductLV(unsigned int jak, int pdgID, unsigned  int number);
 std::vector<TLorentzVector> GetTruthPionsFromA1(unsigned int number=0);
 std::vector<TLorentzVector> GetTruthPionsFromRho(unsigned int number=0);
 



 TLorentzVector MCTau_invisiblePart(unsigned int i);
 TLorentzVector MCTau_visiblePart(unsigned int i);
 
 //Tau and decay products
 int NMCTauDecayProducts(unsigned int i){if(0<=i && i<(unsigned int)NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
 TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
 int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
 unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
 int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
 TVector3 MCTauandProd_Vertex(unsigned int i, unsigned int j){
   return TVector3(Ntp->MCTauandProd_Vertex->at(i).at(j).at(0),Ntp->MCTauandProd_Vertex->at(i).at(j).at(1),Ntp->MCTauandProd_Vertex->at(i).at(j).at(2));
 }
 bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx);
 bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx);
   */
};
#endif
