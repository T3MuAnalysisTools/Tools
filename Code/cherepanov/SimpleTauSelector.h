#ifndef SimpleTauSelector_h
#define SimpleTauSelector_h

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



class SimpleTauSelector : public Selection {

 public:
  SimpleTauSelector(TString Name_, TString id_);
  virtual ~SimpleTauSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx, 
	     SignalCandidate,
	     MuonCandidate,
	     OSCharge,
	     nTaus,
	     NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  std::vector<TH1D>   NumberOfTaus;
  std::vector<TH1D>   Mu3MuVisibleMass;
  std::vector<TH1D>   Mu3MudPhi;
  std::vector<TH2D>   NumTausvsNumMuons;
  std::vector<TH2D> VertexChi2KF_vs_HelixFit;

};
#endif
