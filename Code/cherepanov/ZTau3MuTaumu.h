#ifndef ZTau3MuTaumu_h
#define ZTau3MuTaumu_h

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


class ZTau3MuTaumu : public Selection {

 public:
  ZTau3MuTaumu(TString Name_, TString id_);
  virtual ~ZTau3MuTaumu();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1_TriggerOk=0,
	     HLT_TriggerOk,
	     SignalCandidate,
	     TriggerMatch,
	     TripletPT,
	     nMuons,
	     OSCharge,
	     Tau3MuIsolation,
	     MuonIsolation,
	     VisMass,
	     NCuts}; 

  TString AnalysisName;
 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  std::vector<TH1D>   Tau3MuRelativeIsolation;
  std::vector<TH1D>   OppositeMuRelativeIsolation;
  std::vector<TH1D>   VisibleDiTauMass;
  std::vector<TH1D>   MTT;
  std::vector<TH1D>   TripletMass;
  std::vector<TH1D>   TripletMassWR;
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;
  std::vector<TH1D>   PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PairMass_OppositeSign_dR13;

  std::vector<TH1D>   TripletPt;
  std::vector<TH1D>   OppositeMuonPt;

  std::vector<TH1D>   TripletEta;
  std::vector<TH1D>   OppositeMuonEta;

  std::vector<TH2D>   SingleVsThreeMuTrigger;

  
  TRandom rndm;
  double random_num;


  float  mini_m3m;
  float  mini_dataMCtype;
  float  mini_event_weight;
  float  mini_TripletPt;
  float  mini_MuonPt;
  float  mini_Tau3MuRelativeIsolation;
  float  mini_OppositeMuRelativeIsolation;
  float  mini_ditaumass;


};
#endif
