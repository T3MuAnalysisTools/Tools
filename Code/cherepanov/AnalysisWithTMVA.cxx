#include "AnalysisWithTMVA.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


AnalysisWithTMVA::AnalysisWithTMVA(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

AnalysisWithTMVA::~AnalysisWithTMVA(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  AnalysisWithTMVA::Configure(){

  reader = new TMVA::Reader( "!Color:!Silent" );

  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  reader->AddVariable( "var1", &var1 );
  reader->AddVariable( "var2", &var2 );
  //  reader->AddSpectator("var1*2", &spec1 );  // add spectattor if needed
  reader->BookMVA( "BDT", basedir+"TMVAClassification_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles

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
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  MuonsPt=HConfig.GetTH1D(Name+"_MuonsPt","All Muons Pt",25,0,50,"p_{T} of all stored Muons, GeV","Events");
  MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of 2 muons pt",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
  MuonsEta=HConfig.GetTH1D(Name+"_MuonsEta","All Muons #eta",25,-2.5,2.5,"Rapidity of all stored Muons","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.2,1.5,"Mass of the #mu#mu pair, GeV","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track, GeV","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  AnalysisWithTMVA::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&MuonsPt);
  Extradist1d.push_back(&MuonsEta);
  Extradist1d.push_back(&MuonsPtRatio);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&TripleMass);
  Extradist2d.push_back(&PhiMassVsDsMass);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  AnalysisWithTMVA::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection


  value.at(PrimeVtx)=Ntp->NVtx(); // Here the actual_value of a cut is set
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); // Here we check that the actuall value of PrimeVrtices is above 5.
  
  value.at(TriggerOk)=(Ntp->EventNumber()%1000)==1;
  pass.at(TriggerOk)=true; // always true
  
  double wobs=1;
  double w;
           
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}



  bool status=AnalysisCuts(t,w,wobs);

  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);



    for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++){ // loop over muons
      MuonsPt.at(t).Fill(Ntp->Muon_P4(iMuon).Pt(),w);
      MuonsEta.at(t).Fill(Ntp->Muon_P4(iMuon).Eta(),w);
    }


    double deltaMass(999.);
    unsigned int pair_index(0);
    if(Ntp->NTwoMuonsTrack()!=0){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

      unsigned int muon_1 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
      unsigned int muon_2 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);

      if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
	deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
	pair_index = i2M; // this is an index of the candidate with best mumu mass
      }
    }
    

    MuonsPtRatio.at(t).Fill(Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0)).Pt()/Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1)).Pt(),w );
    PhiMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M(), w);

    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M();


    var1=phimass;
    var2=dsmass;
    spec1 = var1*2;

    //    Float_t var1, var2, spec1;

    std::cout<<"   MVA classifier of the event :  " << reader->EvaluateMVA("BDT") << std::endl;
      


    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    }


  }
}



void  AnalysisWithTMVA::Finish(){


  Selection::Finish();

}





