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

  TString AnalysisName;
 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;




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
  std::vector<TH2D>   SingleVsThreeMuTrigger;

  TRandom rndm;
  double random_num;


  float  mini_m3m;
  float  mini_dataMCtype;
  float  mini_event_weight;
  float  mini_TripletPt;
  float  mini_Tau3MuRelativeIsolation;
  float  mini_TauPt;
  float  mini_ditaumass;





};
#endif
