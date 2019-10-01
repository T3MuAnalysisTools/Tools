#ifndef PUReweight_h
#define PUReweight_h

#include <TH1F.h>
#include <string>
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
class PUReweight {
public:
  enum Type {NONE=0,  RUN2ANALYSIS=1};

  PUReweight(Type type=RUN2ANALYSIS); //default RUN2Analysis 

  ~PUReweight();


  float PUweightHTT(float npu);
public:
  Type theType;
  TH1D *PUWeightHiso;
  
};
#endif
