#ifndef ZTau3MuTauh_h
#define ZTau3MuTauh_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"
#include "SimpleFits/FitSoftware/interface/GEFObject.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"
#include "SimpleFits/FitSoftware/interface/TPTRObject.h"


class ZTau3MuTauh : public Selection {

 public:
  ZTau3MuTauh(TString Name_, TString id_);
  virtual ~ZTau3MuTauh();

  virtual void  Configure();
  virtual void  Finish();
  /*
  enum cuts {L1_TriggerOk=0,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
	     TripletPT,
	     Tau3MuIsolation,
	     nTaus,
	     OSCharge,
	     DeepTauJets,
	     DeepTauMuons,
	     DeepTauElectrons,
	     VisMass,
	     NCuts}; 
  */
  
  enum cuts {WhetherDecayFound=0,
	     Mu1_Candidate,
             Mu2_Candidate,
	     Mu3_Candidate,
	     Tau_h_Candidate,
             Mu1_Candidate_recod,
	     Mu2_Candidate_recod,
	     Mu3_Candidate_recod,
	     Tau_h_Candidate_recod,
             L1_TriggerOk,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
             TripletPT,
             Tau3MuIsolation,
	     nTaus_pT,
             nTaus_eta,
             nTaus_dR,
	     OSCharge,
	     DeepTauJets,
	     DeepTauMuons,
	     DeepTauElectrons,
	     VisMass,
	     NCuts}; 

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  std::vector<TH1D>   NumberOfTaus;
  std::vector<TH1D>   Tau3MuRelativeIsolation;
  std::vector<TH1D>   TauHDecayMode;
  std::vector<TH1D>   VisibleDiTauMass;
  std::vector<TH1D>   MTT;
  std::vector<TH1D>   TripletMass;
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;

  std::vector<TH1D>   PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PairMass_OppositeSign_dR13;

  std::vector<TH1D>   TripletPt;
  std::vector<TH1D>   OppositeTauPt;
  
  std::vector<TH1D>   dR_betweenTruth_NeutrinoGuess;
  std::vector<TH1D>   dR_betweenTruth_Tau;
  std::vector<TH1D>   Z_Pt;
  std::vector<TH2D>   OS_vs_3mu_trigger;
  
  std::vector<TH1D>   Selection_Cut_3mu_Pt;
  std::vector<TH1D>   Selection_Cut_3mu_Rel_Iso;
  std::vector<TH1D>   Selection_Cut_tauh_Pt;
  std::vector<TH1D>   Selection_Cut_tauh_Eta;
  std::vector<TH1D>   Selection_Cut_tauh_DeltaR_3mu;
  std::vector<TH1D>   Selection_Cut_Vis_InvM;
  
  std::vector<TH1D>   Selection_Cut_Mu1_P;
  std::vector<TH1D>   Selection_Cut_Mu1_Eta;
  std::vector<TH1D>   Selection_Cut_Mu1_dR;
  std::vector<TH1D>   Selection_Cut_Mu2_P;
  std::vector<TH1D>   Selection_Cut_Mu2_Eta;
  std::vector<TH1D>   Selection_Cut_Mu2_dR;
  std::vector<TH1D>   Selection_Cut_Mu3_P;
  std::vector<TH1D>   Selection_Cut_Mu3_Eta;
  std::vector<TH1D>   Selection_Cut_Mu3_dR;
  std::vector<TH1D>   Selection_Cut_h_Pt;
  std::vector<TH1D>   Selection_Cut_h_Eta;
  std::vector<TH1D>   Selection_Cut_h_dR;
  
  std::vector<TH1D>   Selection_Cut_RecoMu_P;
  std::vector<TH1D>   Selection_Cut_RecoMu_Eta;
  std::vector<TH1D>   Selection_Cut_RecoH_Pt;
  std::vector<TH1D>   Selection_Cut_RecoH_Eta;
  
  Double_t m3m;
  Double_t dataMCtype;
  Double_t event_weight;
  Double_t m12;
  Double_t m13;
  Double_t mDr1;
  Double_t mDr2;
  Double_t LumiScale;
  
  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;
  
  TRandom rndm;
  double random_num;

};
#endif
