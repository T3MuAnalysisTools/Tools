#include "PUReweight.h"

#include <iostream>
#include <cstdlib>


float PUReweight::PU_weight(float npu) {

  int ibin = PUWeightHiso->FindFixBin(npu);
  return PUWeightHiso->GetBinContent(ibin);

}


PUReweight::PUReweight(Type type) : 
theType(type) {

  TString rfilename = (std::string)std::getenv("workdir")+"/Code/CommonFiles/" + "weights/PU/Data_Pileup_2016_271036-284044_13TeVMoriond17_23Sep2016ReReco_69p2mbMinBiasXS.root";


  TFile *file = TFile::Open(rfilename);
  PUWeightHiso = (TH1D*)file->Get("pileup");
}

PUReweight::~PUReweight() {}
