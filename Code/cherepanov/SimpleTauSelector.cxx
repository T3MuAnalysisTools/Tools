#include "SimpleTauSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

SimpleTauSelector::SimpleTauSelector(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

SimpleTauSelector::~SimpleTauSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SimpleTauSelector::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)        cut.at(TriggerOk)=1;
    if(i==PrimeVtx)         cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
    if(i==nTaus)            cut.at(nTaus)=1;

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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nTaus){
      title.at(i)=" nTaus ";
      hlabel="number of taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,3,-0.5,2.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,3,-0.5,2.5,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams

  NumberOfTaus=HConfig.GetTH1D(Name+"_NumberOfTaus","NumberOfTaus",5,-0.5,4.5,"Number of #tau ","Events");
  NumTausvsNumMuons=HConfig.GetTH2D(Name+"_NumTausvsNumMuons","NumTausvsNumMuons",5,-0.5,4.5,5,-0.5,4.5,"N 3#mu ","N #tau");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SimpleTauSelector::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfTaus);
  Extradist2d.push_back(&NumTausvsNumMuons);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  SimpleTauSelector::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection


  value.at(PrimeVtx)=Ntp->NVtx(); 
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  
  value.at(TriggerOk)=(Ntp->EventNumber()%1000)==1;
  pass.at(TriggerOk)=true; 


  value.at(nTaus) = Ntp->NTaus();
  pass.at(nTaus)  = ( value.at(nTaus) >= cut.at(nTaus) );

  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  
  
  bool status=AnalysisCuts(t,w,wobs);
  NumTausvsNumMuons.at(t).Fill(Ntp->NTaus(), Ntp->NThreeMuons());
  if(status){ 
    NumberOfTaus.at(t).Fill(Ntp->NTaus());

  }
}


void  SimpleTauSelector::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





