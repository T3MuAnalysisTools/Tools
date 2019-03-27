#ifndef MCStudy_h
#define MCStudy_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MCStudy : public Selection {

 public:
  MCStudy(TString Name_, TString id_);
  virtual ~MCStudy();

  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {L1SeedOk=0,HLTOk,PrimeVtx,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> TripleMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D>  MuonsPtRatio;
  std::vector<TH1D> Category;

  std::vector<TH1D> FirstMuonsPt;
  std::vector<TH1D> SecondMuonsPt;
  std::vector<TH1D> ThirdMuonsPt;

  std::vector<TH1D> FirstMuonsEta;
  std::vector<TH1D> SecondMuonsEta;
  std::vector<TH1D> ThirdMuonsEta;


};
#endif