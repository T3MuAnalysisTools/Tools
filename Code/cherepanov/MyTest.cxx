#include "MyTest.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

MyTest::MyTest(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

MyTest::~MyTest(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  MyTest::Configure(){
  ////////////////////////////////////////////////////////////////////////
  // Here you can defined your cuts. There are three vector for cuts:
  // std::vector<double> cut, std::vector<double> value and std::vecto<bool> pass.  
  // For exmaple if you want to aplly a  selection to variables Var1,Var2,Var3
  // with selection values Val1, Val2, Val3, then you have: vector cut = (Val1,Val2,Val3),
  // vector value  = (actual_value1, actual_value2, actual_value3), where the actuala_value
  // is an actual value of Var1,2,3 in a given event (this vector will be filled later)
  // vector pass contains boolean of the cut status.

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  MyTest::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  MyTest::doEvent(){ 
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Here the index t belongs to sample type, this value originally filled at Ntuple filling
  // level: https://github.com/T3MuAnalysisTools/DsTau23Mu/blob/master/T3MNtuple/interface/DataMCType.h
  // but you can flexibly redefined this and make combinations like, Data, MC1, MC2, MC3+MC4, etc ...
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection


  value.at(PrimeVtx)=Ntp->NVtx(); // Here the actual_value of a cut is set
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); // Here we check that the actuall value of PrimeVrtices is above 5.
  
  value.at(TriggerOk)=(Ntp->EventNumber()%1000)==1;
  pass.at(TriggerOk)=true; // always true
  
  double wobs=1;
  double w;  //  This is an event weights, one may intorduce any weights to the event, for exmaple PU. 
             //  there can be several weights, e.g. w = w1*w2*w3 ...
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}



  bool status=AnalysisCuts(t,w,wobs);
  ///////////////////////////////////////////////////////////
  // Add plots
  // The status boolean is true if all elements in pass are true
  // and false if at least one is false: status = true if 
  // pass = (true, true, true ..., true)  and status = false
  // if pass = (true, true, false, ..., true)

  if(status){ // Only selected events pass this if statement
    // Lets fill below some plots ...
    // All available get functions can be found in https://github.com/T3MuAnalysisTools/Tools/blob/master/Code/Ntuple_Controller.h

  }
}


void  MyTest::Finish(){
  Selection::Finish();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





