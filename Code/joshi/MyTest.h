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

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;


};
#endif
