#ifndef MyTest_h
#define MyTest_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MyTest : public Selection {

 public:
  MyTest(TString Name_, TString id_);
  virtual ~MyTest();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts}; // this is a cut enumerator, for other event cuts please put a new enumarator with an understandanle name, for exmaple
  //   enum cuts {TriggerOk=0, PrimeVts, EventCut1, EventCut2, ..., NCuts};  // Do not remove/rename  the last enumerator   NCuts;

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  //  Initialize  rout histrogram

  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> TripleMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D>  MuonsPtRatio;

};
#endif
