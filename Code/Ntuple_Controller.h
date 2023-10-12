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
#include "TMatrixDSymEigen.h"
#include "TDecompBK.h"
// Include files (C & C++ libraries)
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <utility>      // std::pair
#include <tuple>
#include <functional>
#include <algorithm>
#include <set>
#include "NtupleReader.h"
#include "HistoConfig.h"
#include "PDG_Var.h"

#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"


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
//*   Purpose: The purpose of this class is to provide an interface to the
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

      // Ntuple Access Functions
      virtual void Branch_Setup(TString B_Name, int type);
      virtual void Branch_Setup(){}

      // Functions to configure objects
      virtual void ConfigureObjects(); 
      unsigned int ObjEvent;



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
      virtual ~Ntuple_Controller();

      // Event initializer
      void InitEvent();

      enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
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


      bool isInit;

      // Information from input Ntuple path name
      TString GetInputNtuplePath();
      TString GetInputDatasetName();
      TString GetInputPublishDataName();


      // Physics Variable Get Functions
      // Event Variables
      Long64_t GetMCID();
      int GetStrippedMCID();

      ULong64_t EventNumber(){return Ntp->Event_EventNumber;}
      Int_t     RunNumber(){return Ntp->Event_RunNumber;}
      Int_t     DataMC_Type(){return Ntp->Event_DataMC_Type;}
      Int_t     LuminosityBlock(){return Ntp->Event_luminosityBlock;}
      float     METEt(){return Ntp->Event_METEt ;}
      float     METPhi(){return Ntp->Event_METPhi ;}
      float     METXX(){return Ntp->Event_METXX ;}
      float     METXY(){return Ntp->Event_METXY ;}
      float     METYY(){return Ntp->Event_METYY ;}


      float     NVtx(){return Ntp->Vertex_N_primary;}
      double    DeltaPhi(double, double);
      double    TruthNumberOfInteraction(){return Ntp->puN;}
      TString   WhichEra(int year);



      unsigned int   NTracks(){return Ntp->Track_p4->size();}
      TLorentzVector Track_P4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2), Ntp->Track_p4->at(i).at(3), Ntp->Track_p4->at(i).at(0));}
      TLorentzVector KaonTrack_P4(unsigned int index);
      TVector3       Track_Poca(unsigned int i){return TVector3(Ntp->Track_poca->at(i).at(0),Ntp->Track_poca->at(i).at(1),Ntp->Track_poca->at(i).at(2));}
      double         Track_normalizedChi2(unsigned int i){return Ntp->Track_normalizedChi2->at(i);}
      double         Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);}
      double         Track_charge(unsigned int i){return Ntp->Track_charge->at(i);}
      double         Track_dxy(unsigned int i){return  Ntp->Track_dxy->at(i);}
      double         Track_dz(unsigned int i){return  Ntp->Track_dz->at(i);}
      double         Track_dxyError(unsigned int i){return  Ntp->Track_dxyError->at(i);}
      double         Track_dzError(unsigned int i){return  Ntp->Track_dzError->at(i);}






      bool                     Vertex_HighestPtVertexIsValid(){return Ntp->Vertex_HighestPt_PrimaryVertex->size();}
      TVector3                 Vertex_HighestPt_PrimaryVertex(){

	if(Ntp->Vertex_HighestPt_PrimaryVertex->size() == 0) return TVector3(0,0,0);
	return TVector3(Ntp->Vertex_HighestPt_PrimaryVertex->at(0).at(0),
			Ntp->Vertex_HighestPt_PrimaryVertex->at(0).at(1),
			Ntp->Vertex_HighestPt_PrimaryVertex->at(0).at(2));
      } 
      TMatrixTSym<double>  Vertex_HighestPt_PrimaryVertex_Covariance();
	

      unsigned int             NTaus(){return Ntp->Tau_p4->size();}
      int                      Tau_charge(unsigned int i){return Ntp->Tau_charge->at(i);}
      TLorentzVector           Tau_P4(unsigned int i){return TLorentzVector(Ntp->Tau_p4->at(i).at(1),Ntp->Tau_p4->at(i).at(2),Ntp->Tau_p4->at(i).at(3),Ntp->Tau_p4->at(i).at(0));}
      int                      Tau_DecayMode(unsigned int i){return Ntp->Tau_DecayMode->at(i);}
      int                      Tau_DecayModeFinding(unsigned int i){return Ntp->Tau_DecayModeFinding->at(i);}
      int                      Tau_NewDecayModeFinding(unsigned int i){return Ntp->Tau_NewDecayModeFinding->at(i);}
      unsigned int             Tau_byLooseDeepTau2017v2p1VSe(unsigned int i){return Ntp->Tau_byLooseDeepTau2017v2p1VSe->at(i);}
      unsigned int             Tau_byMediumDeepTau2017v2p1VSe(unsigned int i){return Ntp->Tau_byMediumDeepTau2017v2p1VSe->at(i);}
      unsigned int             Tau_byTightDeepTau2017v2p1VSe(unsigned int i){return Ntp->Tau_byTightDeepTau2017v2p1VSe->at(i);}
      unsigned int             Tau_byLooseDeepTau2017v2p1VSmu(unsigned int i){return Ntp->Tau_byLooseDeepTau2017v2p1VSmu->at(i);}
      unsigned int             Tau_byMediumDeepTau2017v2p1VSmu(unsigned int i){return Ntp->Tau_byMediumDeepTau2017v2p1VSmu->at(i);}
      unsigned int             Tau_byTightDeepTau2017v2p1VSmu(unsigned int i){return Ntp->Tau_byTightDeepTau2017v2p1VSmu->at(i);}
      unsigned int             Tau_byLooseDeepTau2017v2p1VSjet(unsigned int i){return Ntp->Tau_byLooseDeepTau2017v2p1VSjet->at(i);}
      unsigned int             Tau_byMediumDeepTau2017v2p1VSjet(unsigned int i){return Ntp->Tau_byMediumDeepTau2017v2p1VSjet->at(i);}
      unsigned int             Tau_byTightDeepTau2017v2p1VSjet(unsigned int i){return Ntp->Tau_byTightDeepTau2017v2p1VSjet->at(i);}
      unsigned int             Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits(unsigned int i){return Ntp->Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(i);}
      unsigned int             Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits(unsigned int i){return Ntp->Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(i);}
      unsigned int             Tau_byTightCombinedIsolationDeltaBetaCorr3Hits(unsigned int i){return Ntp->Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(i);}
      float                    Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits(unsigned int i){return Ntp->Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(i);}
      float                    Tau_RawDiscriminator(unsigned int i, unsigned int j){return Ntp->Tau_FloatDiscriminants->at(i).at(j);}
      int                      Tau_Discriminator(unsigned int i, unsigned int j){return Ntp->Tau_IntDiscriminants->at(i).at(j);}



      int                      Tau_NTracks(){return Ntp->Tau_PFTauTrack_p4->size() ; }  //   This is either 0 or 1  
      TLorentzVector           Tau_PFTauTrack_p4(unsigned int i){return TLorentzVector(Ntp->Tau_PFTauTrack_p4->at(i).at(1),
										       Ntp->Tau_PFTauTrack_p4->at(i).at(2),
										       Ntp->Tau_PFTauTrack_p4->at(i).at(3),
										       Ntp->Tau_PFTauTrack_p4->at(i).at(0));}

      bool                     Tau_SV_isValid(unsigned int i){return Ntp->Tau_SVPos->at(i).size();}
      TVector3                 Tau_SVPos(unsigned int i){
	if(Ntp->Tau_SVPos->at(i).size() == 0) return TVector3(0,0,0);
	return TVector3(Ntp->Tau_SVPos->at(i).at(0),
			Ntp->Tau_SVPos->at(i).at(1),
			Ntp->Tau_SVPos->at(i).at(2));
      } 
      TMatrixTSym<float>  Tau_SVCov(unsigned int i);


      LorentzVectorParticle    Tau_a1_LVP(unsigned int i);
      LorentzVectorParticle    Tau_Track_LVP(unsigned int i);


      unsigned int             NElectrons(){return Ntp->Electron_p4->size();}
      int                      Electron_charge(unsigned int i){return Ntp->Electron_Charge->at(i);}
      TLorentzVector           Electron_P4(unsigned int i){return TLorentzVector(Ntp->Electron_p4->at(i).at(1),
										 Ntp->Electron_p4->at(i).at(2),
										 Ntp->Electron_p4->at(i).at(3),
										 Ntp->Electron_p4->at(i).at(0));}

      float                    Electron_puppiNeutralHadronIso(unsigned int i){return Ntp->Electron_puppiNeutralHadronIso->at(i);}
      float                    Electron_puppiChargedHadronIso(unsigned int i){return Ntp->Electron_puppiChargedHadronIso->at(i);}
      float                    Electron_puppiPhotonIso(unsigned int i){return Ntp->Electron_puppiPhotonIso->at(i);}
      float                    Electron_trackIso(unsigned int i){return Ntp->Electron_trackIso->at(i);}
      float                    Electron_relativeIsolation(unsigned int i){return Ntp->Electron_relativeIsolation->at(i);}
      int                      Electron_isPF(unsigned int i){return Ntp->Electron_isPF->at(i);}
      int                      Electron_cutBasedElectronID_Fall17_94X_V2_veto(unsigned int i){return Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_veto->at(i);}
      int                      Electron_cutBasedElectronID_Fall17_94X_V2_loose(unsigned int i){return Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_loose->at(i);}
      int                      Electron_cutBasedElectronID_Fall17_94X_V2_medium(unsigned int i){return Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_medium->at(i);}
      int                      Electron_cutBasedElectronID_Fall17_94X_V2_tight(unsigned int i){return Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_tight->at(i);}


      unsigned int   NMuons(){return Ntp->Muon_p4->size();}
      TLorentzVector Muon_P4(unsigned int i){return TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2), Ntp->Muon_p4->at(i).at(3), Ntp->Muon_p4->at(i).at(0));}
      TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));}
      bool           Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);} 
      bool           Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);}
      bool           Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
      bool           Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);}
      bool           Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);}
      bool           Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);}
      bool           Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);}
      bool           isLooseMuon(unsigned int i);
      float          Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);}
      float          Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);}
      float          Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);}
      float          Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);}
      int            Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);}
      int            Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);}
      float          Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);}
      float          Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);}
      float          Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);}
      float          Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);}
      float          Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);}
      float          Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);}
      int            Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);}
      int            Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);}
      float          Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);}
      float          Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);}
      float          Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}
      float          Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}
      float          Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}
      float          Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}
      float          Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}
      float          Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}
      float          Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}
      float          Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}
      float          Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}
      float          Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}
      float          Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}
      float          Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}
      float          Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}
      float          Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}
      float          Muon_RelIso(unsigned int i){return Muon_sumChargedHadronPt04(i) + std::max(0., Muon_sumNeutralHadronEt04(i) + Muon_sumPhotonEt04(i) -0.5 * Muon_sumPUPt04(i)) ;}
      double         Muon_ptError(unsigned int i){return Ntp->Muon_ptError->at(i);}
      double         Muon_phiError(unsigned int i){return Ntp->Muon_phiError->at(i);}
      double         Muon_etaError(unsigned int i){return Ntp->Muon_etaError->at(i);}

      TLorentzVector Muon_outerTrack_p4(unsigned int i){return TLorentzVector(Ntp->Muon_outerTrack_p4->at(i).at(1), 
									      Ntp->Muon_outerTrack_p4->at(i).at(2), 
									      Ntp->Muon_outerTrack_p4->at(i).at(3),
									      Ntp->Muon_outerTrack_p4->at(i).at(0));}

      TLorentzVector Muon_innerTrack_p4(unsigned int i){return TLorentzVector(Ntp->Muon_innerTrack_p4->at(i).at(1), 
									      Ntp->Muon_innerTrack_p4->at(i).at(2), 
									      Ntp->Muon_innerTrack_p4->at(i).at(3),
									      Ntp->Muon_innerTrack_p4->at(i).at(0));}
      unsigned int Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);}
      int          Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);}
      int          Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);}
      float        Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);}
      int          Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);}
      int          Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);}
      int          Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);}
      int          Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);}
      bool         Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}
      bool         Muon_isRPCMuon(unsigned int i){return Ntp->Muon_isRPCMuon->at(i);}
      int          Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);} 
      int          Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);}
      bool         Muon_combinedQuality_updatedSta(unsigned int i){return Ntp->Muon_combinedQuality_updatedSta->at(i);}
      double       Muon_combinedQuality_trkKink(unsigned int i){return Ntp->Muon_combinedQuality_trkKink->at(i);}
      double       Muon_combinedQuality_glbKink(unsigned int i){return Ntp->Muon_combinedQuality_glbKink->at(i);}
      double       Muon_combinedQuality_trkRelChi2(unsigned int i){return Ntp->Muon_combinedQuality_trkRelChi2->at(i);}
      double       Muon_combinedQuality_staRelChi2(unsigned int i){return Ntp->Muon_combinedQuality_staRelChi2->at(i);}
      double       Muon_combinedQuality_chi2LocalPosition(unsigned int i){return Ntp->Muon_combinedQuality_chi2LocalPosition->at(i);}
      double       Muon_combinedQuality_chi2LocalMomentum(unsigned int i){return Ntp->Muon_combinedQuality_chi2LocalMomentum->at(i);}
      double       Muon_combinedQuality_localDistance(unsigned int i){return Ntp->Muon_combinedQuality_localDistance->at(i);}
      double       Muon_combinedQuality_globalDeltaEtaPhi(unsigned int i){return Ntp->Muon_combinedQuality_globalDeltaEtaPhi->at(i);}
      bool         Muon_combinedQuality_tightMatch(unsigned int i){return Ntp->Muon_combinedQuality_tightMatch->at(i);}
      double       Muon_combinedQuality_glbTrackProbability(unsigned int i){return Ntp->Muon_combinedQuality_glbTrackProbability->at(i);}
      double       Muon_prod_inner_outer_charge(unsigned int i){return Ntp->Muon_prod_inner_outer_charge->at(i);}

      double       Muon_innerTrack_quality(unsigned int i){return Ntp->Muon_innerTrack_quality->at(i);}
      double       Muon_ptErrOverPt(unsigned int i){return Ntp->Muon_ptErrOverPt->at(i);}
      double       Muon_calEnergy_hadS9(unsigned int i){return Ntp->Muon_calEnergy_hadS9->at(i);}
      double       Muon_calEnergy_had(unsigned int i){return Ntp->Muon_calEnergy_had->at(i);}
      double       Muon_calEnergy_emS25(unsigned int i){return Ntp->Muon_calEnergy_emS25->at(i);}
      double       Muon_calEnergy_emS9(unsigned int i){return Ntp->Muon_calEnergy_emS9->at(i);}
      double       Muon_calEnergy_em(unsigned int i){return Ntp->Muon_calEnergy_em->at(i);}
      double       Muon_segmentCompatibility(unsigned int i){return Ntp->Muon_segmentCompatibility->at(i);}
      double       Muon_caloCompatibility(unsigned int i){return Ntp->Muon_caloCompatibility->at(i);}
      bool         Muon_isGoodMuon_TM2DCompatibility(unsigned int i){return Ntp->Muon_isGoodMuon_TM2DCompatibility->at(i);}
      double       Muon_innerTrack_validFraction(unsigned int i){return Ntp->Muon_innerTrack_validFraction->at(i);}
      double       Muon_innerTrack_pixelLayersWithMeasurement(unsigned int i){return Ntp->Muon_innerTrack_pixelLayersWithMeasurement->at(i);}
      double       Muon_innerTrack_numberOfValidTrackerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfValidTrackerHits->at(i);}
      double       Muon_innerTrack_numberOfLostTrackerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerHits->at(i);}
      double       Muon_innerTrack_numberOfLostTrackerInnerHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits->at(i);}
      double       Muon_innerTrack_numberOfLostTrackerOuterHits(unsigned int i){return Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits->at(i);}
      double       Muon_innerTrack_normalizedChi2(unsigned int i){return Ntp->Muon_innerTrack_normalizedChi2->at(i);}
      double       Muon_outerTrack_normalizedChi2(unsigned int i){return Ntp->Muon_outerTrack_normalizedChi2->at(i);}
      double       Muon_outerTrack_muonStationsWithValidHits(unsigned int i){return Ntp->Muon_outerTrack_muonStationsWithValidHits->at(i);}
      bool         Muon_isGoodMuon_TrackerMuonArbitrated(unsigned int i){return Ntp->Muon_isGoodMuon_TrackerMuonArbitrated->at(i);}
      bool         Muon_isGoodMuon_TMOneStationTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMOneStationTight->at(i);}
      bool         Muon_isGoodMuon_TMOneStationAngTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMOneStationAngTight->at(i);}
      bool         Muon_isGoodMuon_TMLastStationTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationTight->at(i);}
      bool         Muon_isGoodMuon_TMLastStationAngTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationAngTight->at(i);}
      bool         Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight->at(i);}
      bool         Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(unsigned int i){return Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight->at(i);}
      double       Muon_vmuonhitcomb_reco(unsigned int i){return Ntp->Muon_vmuonhitcomb_reco->at(i);}
      double       Muon_rpchits_reco(unsigned int i){return Ntp->Muon_rpchits_reco->at(i);}
      int          Muon_charge(unsigned int i){return Ntp->Muon_charge->at(i);}
      int          Muon_trackCharge(unsigned int i){return Ntp->Muon_trackCharge->at(i);}
      int          Muon_pdgid(unsigned int i){return Ntp->Muon_pdgid->at(i);}
      float        Muon_B(unsigned int i){return Ntp->Muon_B->at(i);}
      float        Muon_M(unsigned int i){return Ntp->Muon_M->at(i);}

      int          Muon_ID(unsigned int i){return Ntp->Muon_ID->at(i);}
      int          Muon_StandardSelection(unsigned int i){return Ntp->Muon_StandardSelection->at(i);}

      float Muon_TrackX(unsigned int i, unsigned int det){return Ntp->Muon_TrackX->at(i).at(det);}
      float Muon_TrackY(unsigned int i, unsigned int det){return Ntp->Muon_TrackY->at(i).at(det);}
      float Muon_dDxDz(unsigned int i, unsigned int det){return Ntp->Muon_dDxDz->at(i).at(det);} //  det = 0,1,2,3 - DT, det = 4,5,6,7 - CSC
      float Muon_dDyDz(unsigned int i, unsigned int det){return Ntp->Muon_dDyDz->at(i).at(det);}
      float Muon_dX(unsigned int i, unsigned int det){return Ntp->Muon_dX->at(i).at(det);}
      float Muon_dY(unsigned int i, unsigned int det){return Ntp->Muon_dY->at(i).at(det);}
      float Muon_pullX(unsigned int i, unsigned int det){return Ntp->Muon_pullX->at(i).at(det);}
      float Muon_pullY(unsigned int i, unsigned int det){return Ntp->Muon_pullY->at(i).at(det);}
      float Muon_pullDxDz(unsigned int i, unsigned int det){return Ntp->Muon_pullDxDz->at(i).at(det);}
      float Muon_pullDyDz(unsigned int i, unsigned int det){return Ntp->Muon_pullDyDz->at(i).at(det);}
      float numberOfSegments(unsigned int i, unsigned int det){return Ntp->numberOfSegments->at(i).at(det);}

      float Muon_timeAtIpInOut(unsigned int i){return Ntp->Muon_timeAtIpInOut->at(i);}
      float Muon_timeAtIpInOutErr(unsigned int i){return Ntp->Muon_timeAtIpInOutErr->at(i);}

      TLorentzVector MuonToKaon(unsigned int index);
      TLorentzVector MuonToPion(unsigned int index);
      //Vertex_NMuonsAssocWithPV

      int NPhotons(){return Ntp->Gamma_P4->size();}
      TLorentzVector Photon_p4(unsigned int i){return TLorentzVector(Ntp->Gamma_P4->at(i).at(1),Ntp->Gamma_P4->at(i).at(2), Ntp->Gamma_P4->at(i).at(3), Ntp->Gamma_P4->at(i).at(0));}
      float Photon_e1x5(unsigned int i){return Ntp->Gamma_e1x5->at(i);}
      float Photon_e2x5(unsigned int i){return Ntp->Gamma_e2x5->at(i);}
      float Photon_e3x3(unsigned int i){return Ntp->Gamma_e3x3->at(i);}
      float Photon_e5x5(unsigned int i){return Ntp->Gamma_e5x5->at(i);}

      int Photon_hasPixelSeed(unsigned int i){return Ntp->Gamma_hasPixelSeed->at(i);}
      int Photon_hasConversionTracks(unsigned int i){return Ntp->Gamma_hasConversionTracks->at(i);}
      int Photon_isPF(unsigned int i){return Ntp->Gamma_isPFPhoton->at(i);}

      LorentzVectorParticle Tau3mu_LVP(unsigned int i);






      bool Muon_TrackParticleHasMomentum(unsigned int i){if(Ntp->Muon_par->at(i).size()!=0)return true; return false;} 
      TrackParticle Muon_TrackParticle(unsigned int i){
	TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1);
	TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar);
	unsigned int l=0;
	for(int k=0; k<TrackParticle::NHelixPar; k++){
	  mu_par(k,0)=Ntp->Muon_par->at(i).at(k);
	  for(int j=k; j<TrackParticle::NHelixPar; j++){
	    mu_cov(k,j)=Ntp->Muon_cov->at(i).at(l);
	    mu_cov(j,k)=Ntp->Muon_cov->at(i).at(l);
	    l++;
	  }
	}
	return TrackParticle(mu_par,mu_cov,Ntp->Muon_pdgid->at(i),Ntp->Muon_M->at(i),Ntp->Muon_charge->at(i),Ntp->Muon_B->at(i));
	//	return TrackParticle(mu_par,mu_cov,0,0,0,0);
      }


      bool MCEventIsReconstructed(){return Ntp->MC_isReco;}

      unsigned int   NThreeMuons(){return Ntp->ThreeMuons_index->size();}
      std::vector<unsigned int> ThreeMuonIndices(unsigned int i){return Ntp->ThreeMuons_index->at(i);}

      double ThreeMuons_SV_Chi2(unsigned int i){return Ntp->ThreeMuons_SV_Chi2->at(i);}
      std::vector<float> ThreeMuons_TriggerMatch_dR(unsigned int i){return Ntp->ThreeMuons_TriggerMatch_dR->at(i);}

      double TwoMuonsTrack_SV_Chi2(unsigned int i){return Ntp->TwoMuonsTrack_SV_Chi2->at(i);}
      std::vector<float> TwoMuonsTrack_TriggerMatch_dR(unsigned int i){return Ntp->TwoMuonsTrack_TriggerMatch_dR->at(i);}

      double FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov );
      double TransverseFlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov );


      unsigned int   NTwoMuonsTrack(){return Ntp->TwoMuonsTrack_Muonsindex->size();}
      std::vector<unsigned int> TwoMuonsTrackMuonIndices(unsigned int i){return Ntp->TwoMuonsTrack_Muonsindex->at(i);}
      std::vector<unsigned int> TwoMuonsTrackTrackIndex(unsigned int i){return Ntp->TwoMuonsTrack_Trackindex->at(i);}

      // Vertex variables take an additional argument corresponding to the channel (Singal(false), Normalization(true))
      // default channel is set to false
      int        NumberOfSVertices(){return Ntp->Vertex_signal_KF_Chi2->size();}// should coincide with number of candidate
      double     Vertex_Signal_KF_Chi2(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_KF_Chi2->at(index);
      }
      TVector3   Vertex_Signal_KF_pos(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TVector3(Ntp->Vertex_signal_KF_pos->at(index).at(0), Ntp->Vertex_signal_KF_pos->at(index).at(1),Ntp->Vertex_signal_KF_pos->at(index).at(2));
      }
      
      TMatrixTSym<double>   Vertex_Signal_KF_Covariance(unsigned int i, bool channel=false);
      TMatrixTSym<double>   Vertex_PrimaryVertex_Covariance(unsigned int i, bool channel=false);



      double     Vertex_2MuonsIsoTrack_KF_Chi2(unsigned int i, bool channel=false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->Vertex_2MuonsIsoTrack_KF_Chi2->at(index);
      }
      TVector3   Vertex_2MuonsIsoTrack_KF_pos(unsigned int i, bool channel=false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return TVector3(Ntp->Vertex_2MuonsIsoTrack_KF_pos->at(index).at(0), Ntp->Vertex_2MuonsIsoTrack_KF_pos->at(index).at(1),Ntp->Vertex_2MuonsIsoTrack_KF_pos->at(index).at(2));
      }

      TMatrixTSym<double>   Vertex_2MuonsIsoTrack_KF_Covariance(unsigned int i, bool channel=false);

      //   int NTracksInThePV(unsigned int i){return Ntp->IsolationBranch_Trackp4->at(i).size();}
      //   TLorentzVector PrimaryVertexTrack_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->IsolationBranch_Trackp4->at(i).at(j).at(1), 
      //												      Ntp->IsolationBranch_Trackp4->at(i).at(j).at(2), 
      //												      Ntp->IsolationBranch_Trackp4->at(i).at(j).at(3), 
      //												      Ntp->IsolationBranch_Trackp4->at(i).at(j).at(0));}

      // ----------  Tracks in the signal tau cone
      int NIsolationTrack(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_p4->at(index).size();
      }
      TLorentzVector IsolationTrack_p4(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TLorentzVector(Ntp->IsolationTrack_p4->at(index).at(j).at(1), 
            Ntp->IsolationTrack_p4->at(index).at(j).at(2),
            Ntp->IsolationTrack_p4->at(index).at(j).at(3),
            Ntp->IsolationTrack_p4->at(index).at(j).at(0));
      }



      //------------------  functions in the box are valid ONLY for a signal category, not control channel
      int IsolationTrack_VertexWithSignalMuon1IsValid(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonIsValid->at(index).at(j).at(0);
      }
      
      int IsolationTrack_VertexWithSignalMuon2IsValid(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonIsValid->at(index).at(j).at(1);
      }
      
      int IsolationTrack_VertexWithSignalMuon3IsValid(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonIsValid->at(index).at(j).at(2);
      }

      float IsolationTrack_VertexWithSignalMuon1Chi2(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonChi2->at(index).at(j).at(0);
      }

      float IsolationTrack_VertexWithSignalMuon2Chi2(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonChi2->at(index).at(j).at(1);
      }

      float IsolationTrack_VertexWithSignalMuon3Chi2(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return Ntp->IsolationTrack_VertexWithSignalMuonChi2->at(index).at(j).at(2);
      }


      TVector3 IsolationTrack_VertexWithSignalMuon1Position(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return TVector3(Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(0), 
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(1),
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(2));

      }


      TVector3 IsolationTrack_VertexWithSignalMuon2Position(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return TVector3(Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(3), 
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(4),
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(5));

      }

      TVector3 IsolationTrack_VertexWithSignalMuon3Position(unsigned int i, unsigned int j, bool channel = false){
	unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	return TVector3(Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(6), 
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(7),
			Ntp->IsolationTrack_VertexWithSignalMuonPosition->at(index).at(j).at(8));

      }

      //------------------  functions in the box are valid ONLY for a signal category, not control channel


      int IsolationTrack_charge(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_charge->at(index).at(j);
      }
      double IsolationTrack_dxySV(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_dxySV->at(index).at(j);
      }
      double IsolationTrack_dzSV(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_dzSV->at(index).at(j);
      }
      double IsolationTrack_dxyPV(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_dxyPV->at(index).at(j);
      }
      double IsolationTrack_dzPV(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_dzPV->at(index).at(j);
      }
      double IsolationTrack_DocaMu1(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_DocaMu1->at(index).at(j);
      }
      double IsolationTrack_DocaMu2(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_DocaMu2->at(index).at(j);
      }
      double IsolationTrack_DocaMu3(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->IsolationTrack_DocaMu3->at(index).at(j);
      }
      

      int      IsolationTrack_NTracks(){return Ntp->IsolationTrack_Helcharge->size();}
      int      IsolationTrack_Helcharge(unsigned int i){return Ntp->IsolationTrack_Helcharge->at(i);}
      int      IsolationTrack_pdgid(unsigned int i){return Ntp->IsolationTrack_pdgid->at(i);}
      float      IsolationTrack_B(unsigned int i){return Ntp->IsolationTrack_B->at(i);}
      float      IsolationTrack_M(unsigned int i){return Ntp->IsolationTrack_M->at(i);}

      bool IsolationTrack_TrackParticleHasMomentum(unsigned int i){if(Ntp->IsolationTrack_par->at(i).size()!=0)return true; return false;}
      TrackParticle IsolationTrack_TrackParticle(unsigned int i){
        TMatrixT<float>    IsolationTrack_par(TrackParticle::NHelixPar,1);
        TMatrixTSym<float> IsolationTrack_cov(TrackParticle::NHelixPar);
        unsigned int l=0;
        for(int k=0; k<TrackParticle::NHelixPar; k++){
          IsolationTrack_par(k,0)=Ntp->IsolationTrack_par->at(i).at(k);
          for(int j=k; j<TrackParticle::NHelixPar; j++){
            IsolationTrack_cov(k,j)=Ntp->IsolationTrack_cov->at(i).at(l);
	    IsolationTrack_cov(j,k)=Ntp->IsolationTrack_cov->at(i).at(l);
            l++;
          }
        }
	TMatrixTSym<float> IsolationTrack_cov_regularised = RegulariseCovariance(IsolationTrack_cov, 1.0);
        return TrackParticle(IsolationTrack_par,IsolationTrack_cov_regularised,Ntp->IsolationTrack_pdgid->at(i),Ntp->IsolationTrack_M->at(i),Ntp->IsolationTrack_Helcharge->at(i),Ntp->IsolationTrack_B->at(i));
	//        return TrackParticle(IsolationTrack_par,IsolationTrack_cov,Ntp->IsolationTrack_pdgid->at(i),Ntp->IsolationTrack_M->at(i),Ntp->IsolationTrack_Helcharge->at(i),Ntp->IsolationTrack_B->at(i));
      }


      TVectorD EigenValues(TMatrixTSym<double> M);
      TMatrixD EigenVectors(TMatrixTSym<double> M);
      TMatrixTSym<double> RegulariseCovariance(TMatrixTSym<double>  M, double coef);
      TMatrixTSym<double> InvertMatrix(TMatrixTSym<double>   M);


      ///// closest distance between the tracks of a candidate
      double     Vertex_DCA12(unsigned int i, bool channel=false) {
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_dca_reco->at(index).at(0);
      }
      double     Vertex_DCA23(unsigned int i, bool channel=false) {
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_dca_reco->at(index).at(1);
      }
      double     Vertex_DCA31(unsigned int i, bool channel=false) {
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_dca_reco->at(index).at(2);
      }
      double     Vertex_DCAMax(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_dca_reco->at(index).at(3);
      }

      TLorentzVector Vertex_signal_KF_refittedTracksP4(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); 
         return TLorentzVector(Ntp->Vertex_signal_KF_refittedTracksP4->at(index).at(j).at(1),Ntp->Vertex_signal_KF_refittedTracksP4->at(index).at(j).at(2),Ntp->Vertex_signal_KF_refittedTracksP4->at(index).at(j).at(3),Ntp->Vertex_signal_KF_refittedTracksP4->at(index).at(j).at(0));
      }

      double     Vertex_d0_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2
      double     Vertex_d0sig_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0sig_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_dz_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_dz_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_d0SV_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0SV_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_dzSV_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_dzSV_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_d0sigSV_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0sigSV_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_d0BeamSpot_reco(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0BeamSpot_reco->at(index).at(j);
      }//  j - is a track number; j=0,1,2
      double     Vertex_d0BeamSpot_reco_sig(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_d0BeamSpot_reco_sig->at(index).at(j);
      }//  j - is a track number; j=0,1,2

      double     Vertex_2Ddisplacement(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_2Ddisplacement->at(index).at(j);
      }// j=0 - value, j=1 - significane, j=3 - distance

      double     Vertex_3Ddisplacement(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_3Ddisplacement->at(index).at(j);
      }// j=0 - value, j=1 - significane, j=3 - distance

      double     Vertex_pair_quality(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_pair_quality->at(index).at(j);
      }//the second index should not exceed 3

      bool       Vertex_pairfit_status(unsigned int i, unsigned int j, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_pairfit_status->at(index).at(j);
      }//the second index should not exceed 3


      TVector3  Vertex_Pair12_Pos(unsigned int i, bool channel = false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	 return TVector3(Ntp->Vertex_Pair12_Pos->at(index).at(0),Ntp->Vertex_Pair12_Pos->at(index).at(1),Ntp->Vertex_Pair12_Pos->at(index).at(2));
      }

      TVector3  Vertex_Pair23_Pos(unsigned int i, bool channel = false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	 return TVector3(Ntp->Vertex_Pair23_Pos->at(index).at(0),Ntp->Vertex_Pair23_Pos->at(index).at(1),Ntp->Vertex_Pair23_Pos->at(index).at(2));
      }

      TVector3  Vertex_Pair31_Pos(unsigned int i, bool channel = false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
	 return TVector3(Ntp->Vertex_Pair31_Pos->at(index).at(0),Ntp->Vertex_Pair31_Pos->at(index).at(1),Ntp->Vertex_Pair31_Pos->at(index).at(2));
      }

      double     Vertex_signal_KF_Chi2(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_KF_Chi2->at(index); 
      }

      double     Vertex_signal_KF_BS_2Ddistance(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_KF_BS_2Ddistance->at(index); 
      }

      double     Vertex_signal_KF_BS_error(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_KF_BS_error->at(index); 
      }

      double     Vertex_signal_KF_BS_significance(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_KF_BS_significance->at(index); 
      }

      TVector3   Vertex_signal_AF_pos(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TVector3(Ntp->Vertex_signal_AF_pos->at(index).at(0),Ntp->Vertex_signal_AF_pos->at(index).at(1),Ntp->Vertex_signal_AF_pos->at(index).at(2));
      }

      double     Vertex_signal_AF_Chi2(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_AF_Chi2->at(index);
      }

      double     Vertex_signal_AF_Ndf(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_signal_AF_Ndf->at(index); }

      TVector3   Vertex_MatchedPrimaryVertex(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TVector3(Ntp->Vertex_MatchedPrimaryVertex->at(index).at(0), Ntp->Vertex_MatchedPrimaryVertex->at(index).at(1),Ntp->Vertex_MatchedPrimaryVertex->at(index).at(2));
      }
      bool       Vertex_RefitPVisValid(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return Ntp->Vertex_RefitPVisValid->at(index);
      }

      TVector3   Vertex_MatchedRefitPrimaryVertex(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TVector3(Ntp->Vertex_MatchedRefitPrimaryVertex->at(index).at(0), 
            Ntp->Vertex_MatchedRefitPrimaryVertex->at(index).at(1),
            Ntp->Vertex_MatchedRefitPrimaryVertex->at(index).at(2));}

      TVector3   Vertex_SecondBestPrimaryVertex(unsigned int i, bool channel=false){
         unsigned int index = i + channel*Ntuple_Controller::NThreeMuons();
         return TVector3(Ntp->Vertex_SecondBestPrimaryVertex->at(index).at(0),
            Ntp->Vertex_SecondBestPrimaryVertex->at(index).at(1),
            Ntp->Vertex_SecondBestPrimaryVertex->at(index).at(2));}

      TVector3   SVPVDirection(TVector3 SV, TVector3 PV);

      float      Isolation_RelPt(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation1->at(index).at(0);}
      float      Isolation_NTracks(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation1->at(index).at(1);}
      float      Isolation_MinDist(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation1->at(index).at(2);}

      float      Isolation05_RelPt(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation2->at(index).at(0);}
      float      Isolation05_NTracks(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation2->at(index).at(1);}
      float      Isolation05_MinDist(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation2->at(index).at(3);}

      float      Isolation_Mu1RelIso(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation3->at(index).at(0);}
      float      Isolation_Mu2RelIso(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation3->at(index).at(1);}
      float      Isolation_Mu3RelIso(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation3->at(index).at(2);}
      float      Isolation_MuMaxRelIso(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation3->at(index).at(3);}

      float      Isolation_Ntrk1(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(0);}
      float      Isolation_Ntrk2(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(1);}
      float      Isolation_Ntrk3(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(2);}
      float      Isolation_Ntrk0p1(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(3);}
      float      Isolation_Ntrk0p2(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(4);}
      float      Isolation_Ntrk0p5(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(5);}
      float      Isolation_maxdy(unsigned int i, bool channel=false){unsigned int index = i + channel*Ntuple_Controller::NThreeMuons(); return Ntp->Vertex_Isolation4->at(index).at(6);}
      
      // ==================== SV from secondary vertex collection =====================
      // --------------------  secondary vertices found by the secondray vertex finder; Note: this is not the signal SV, but the signal SV is included.
      unsigned int   NSecondaryVertices(){return Ntp->SV_pos->size();}
      unsigned int   NTracksAtSecondaryVertex(unsigned int i){return Ntp->SV_Track_P4->at(i).size();}
      float          SecondaryVertexMass(unsigned int i){return Ntp->SV_Mass->at(i);}

      TVector3       SecondaryVertexPosition(unsigned int i){return TVector3(Ntp->SV_pos->at(i).at(0), Ntp->SV_pos->at(i).at(1), Ntp->SV_pos->at(i).at(2));}
      TLorentzVector SecondaryVertexTrack_P4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->SV_Track_P4->at(i).at(j).at(1),
            Ntp->SV_Track_P4->at(i).at(j).at(2),
            Ntp->SV_Track_P4->at(i).at(j).at(3),
            Ntp->SV_Track_P4->at(i).at(j).at(0));}
      
      int            SecondaryVertexTrackCharge(unsigned int i, unsigned int j){return Ntp->SV_TrackCharge->at(i).at(j);}
      TMatrixTSym<float>   SecondaryVertexCovariance(unsigned int i);
      // ==============================================================================


      int NJets(){return Ntp->Jet_p4->size();}
      TLorentzVector Jet_P4(unsigned int i){return TLorentzVector(Ntp->Jet_p4->at(i).at(1), Ntp->Jet_p4->at(i).at(2), Ntp->Jet_p4->at(i).at(3), Ntp->Jet_p4->at(i).at(0));}
      double JetBTagCVSB(unsigned int i){return Ntp->Jet_BTagCVSB->at(i);}
      double JetBTagMVA(unsigned int i){return Ntp->Jet_BTagMVA->at(i);}
      double JetBTagCSV(unsigned int i){return Ntp->Jet_BTagCSV->at(i);}



      int        NL1Seeds(){return Ntp->Trigger_l1name->size();}
      std::string     L1Name(unsigned int i){return Ntp->Trigger_l1name->at(i);}
      int        L1Decision(unsigned int i){return Ntp->Trigger_l1decision->at(i);}
      int        Trigger_l1prescale(unsigned int i){return Ntp->Trigger_l1prescale->at(i);}
      int        NHLT(){return Ntp->Trigger_hltname->size();}
      std::string     HLTName(unsigned int i){return Ntp->Trigger_hltname->at(i);}
      int        HLTDecision(unsigned int i){return Ntp->Trigger_hltdecision->at(i);}
      std::string TriggerObject_name(unsigned int i){return Ntp->TriggerObject_name->at(i);}
      double TriggerObject_pt(unsigned int i){return Ntp->TriggerObject_pt->at(i);}
      double TriggerObject_eta(unsigned int i){return Ntp->TriggerObject_eta->at(i);}
      double TriggerObject_phi(unsigned int i){return Ntp->TriggerObject_phi->at(i);}
      int NTriggerObjects(){return Ntp->TriggerObject_name->size();}
      std::pair<bool, std::vector<float>> triggerMatchTriplet(std::vector<TLorentzVector>, std::vector<TLorentzVector>);


      unsigned int               NMCParticles(){return Ntp->MC_p4->size();}
      TLorentzVector             MCParticle_p4(unsigned int i){return TLorentzVector(Ntp->MC_p4->at(i).at(1),Ntp->MC_p4->at(i).at(2),Ntp->MC_p4->at(i).at(3),Ntp->MC_p4->at(i).at(0));}
      int                        MCParticle_pdgid(unsigned int i){return Ntp->MC_pdgid->at(i);}
      int                        MCParticle_charge(unsigned int i){return Ntp->MC_charge->at(i);}
      int                        MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
      std::vector<int>           MCParticle_childpdgid(unsigned int i){return Ntp->MC_childpdgid->at(i);}
      std::vector<int>           MCParticle_childidx(unsigned int i){return Ntp->MC_childidx->at(i);}
      int  MCParticle_status(unsigned int i){return Ntp->MC_status->at(i);}

      int   getMatchTruthIndex(TLorentzVector tvector);
      int    matchTruth(TLorentzVector tvector);




      unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
      TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
      TVector3                      MCSignalParticle_SourceVertex(unsigned int i){return TVector3(Ntp->MCSignalParticle_SourceVertex->at(i).at(0),
            Ntp->MCSignalParticle_SourceVertex->at(i).at(1),
            Ntp->MCSignalParticle_SourceVertex->at(i).at(2));}
      TVector3                      MCSignalParticle_Vertex(unsigned int i){return TVector3(Ntp->MCSignalParticle_Vertex->at(i).at(0),
            Ntp->MCSignalParticle_Vertex->at(i).at(1),
            Ntp->MCSignalParticle_Vertex->at(i).at(2));}

      int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
      int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
      std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}


      int                        NMCSignalParticleSources(unsigned int i){return  Ntp->MCSignalParticle_Sourcepdgid->at(i).size();}
      int                        MCSignalParticle_Sourcepdgid(unsigned int i, unsigned int j){return Ntp->MCSignalParticle_Sourcepdgid->at(i).at(j);}
      TLorentzVector             MCSignalParticle_Sourcep4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCSignalParticle_Sourcep4->at(i).at(j).at(1), Ntp->MCSignalParticle_Sourcep4->at(i).at(j).at(2), Ntp->MCSignalParticle_Sourcep4->at(i).at(j).at(3), Ntp->MCSignalParticle_Sourcep4->at(i).at(j).at(0));}



      int                        MCSignalParticle_Nchilds(unsigned int i){return Ntp->MCSignalParticle_childp4->at(i).size();}
      int                        MCSignalParticle_childpdgid(unsigned int i,unsigned int j){return Ntp->MCSignalParticle_childpdgid->at(i).at(j);}
      TLorentzVector             MCSignalParticle_child_p4(unsigned int i, unsigned int j){
         return TLorentzVector(Ntp->MCSignalParticle_childp4->at(i).at(j).at(1),
               Ntp->MCSignalParticle_childp4->at(i).at(j).at(2),
               Ntp->MCSignalParticle_childp4->at(i).at(j).at(3),
               Ntp->MCSignalParticle_childp4->at(i).at(j).at(0));}



      // Tau decays (Tau is first element of vector)
      int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
      TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
      int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
      int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}
      int MCTau_midx(unsigned int i){return Ntp->MCTauandProd_midx->at(i);}
      bool  MCParticle_hasMother(unsigned int i){return Ntp->MC_midx->at(i) >= 0;}
      //Tau and decay products
      int NMCTauDecayProducts(unsigned int i){if(0<=i && i<(unsigned int)NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
      TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}

      TVector3 MCTauandProd_Vertex(unsigned int i, unsigned int j){return TVector3(Ntp->MCTauandProd_Vertex->at(i).at(j).at(0), 
            Ntp->MCTauandProd_Vertex->at(i).at(j).at(1),
            Ntp->MCTauandProd_Vertex->at(i).at(j).at(2));}
      int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
      int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}


      //Tool functions
      std::vector<unsigned int> SortedPtMuons(std::vector<unsigned int> indixes);
      std::vector<int> SortedPtMuons_MC(std::vector<int> indices);
      std::vector<unsigned int> SortedEtaMuons(std::vector<unsigned int> indixes);
      std::vector<unsigned int> SortedChargeMuons(std::vector<unsigned int> indices);


      double TauMassResolution(std::vector<unsigned int>  indices, int type, bool UseRefited); // type = 0 - only pt propagation; type = 1 - indluding direction
      TLorentzVector MatchedLV(std::vector<TLorentzVector> list, unsigned int index);

      float DsGenMatch(unsigned int tmp_idx);
      float ThreeMuonDsGenMatch(unsigned int tmp_idx);
      float TauGenMatch(unsigned int tmp_idx);
      bool isPromptDs();

      int DsMotherPdgId(unsigned int index);
      std::pair<int, int> GENMatchedPdgId(TLorentzVector vec);
      TLorentzVector GENMatchedLV(TLorentzVector vec);
      double deltaR(double eta1, double phi1, double eta2, double phi2); 
      TLorentzVector matchToTruthTauDecay(TLorentzVector vector);



      std::vector<int> MuonStandardSelectorBitMask(unsigned int MuonIndex);
      std::vector<int> MuonCustomID(unsigned int MuonIndex);
      TMatrixT<double> convertToMatrix(TVectorT<double> V);

      TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
      TVector3 Rotate(TVector3 LVec, TVector3 Rot);


      template<typename T>
         void printVec(int size, T& vec);


      bool CHECK_BIT(int var, int pos){  return ((var & (1 << pos)) == (1 << pos)); }
      void  printMCDecayChainOfParticle(unsigned int index, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // full event decay chain

      void printMCDecayChain(unsigned int par, unsigned int level = 0, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false);
      void printMCDecayChainOfMother(unsigned int i, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // decay chain of object i
      void printMCDecayChainOfEvent(bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // full event decay chain
      std::string  MCParticleToString(unsigned int par, bool printStatus = false, bool printPt = false, bool printEtaPhi = false);



};
#endif
