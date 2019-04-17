#include "TwoMuTrack.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

TwoMuTrack::TwoMuTrack(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

TwoMuTrack::~TwoMuTrack(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  TwoMuTrack::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==HLTOk)               cut.at(HLTOk)=1;
    if(i==TwoMuTrackCandidate)    cut.at(TwoMuTrackCandidate)=1;
  }
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut

    if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==TwoMuTrackCandidate){
      title.at(i)="TwoMuTrackCandidate is found ";
      hlabel="is 2 Mu + 1 Track  Candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TwoMuTrackCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TwoMuTrackCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  MuonsPt=HConfig.GetTH1D(Name+"_MuonsPt","All Muons Pt",25,0,50,"p_{T} of all stored Muons, GeV","Events");
  Category=HConfig.GetTH1D(Name+"_Category","Category",2,0.5,2.5,"1 - #mu#mu#mu; 2 - #mu#mu+track ","Events");
  MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of 2 muons pt",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
  MuonsEta=HConfig.GetTH1D(Name+"_MuonsEta","All Muons #eta",25,-2.5,2.5,"Rapidity of all stored Muons","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",75,0.2,1.4,"Mass of the #mu#mu pair, GeV","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track, GeV","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",75,0.2,1.4,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");

  FirstMuonsPt=HConfig.GetTH1D(Name+"_FirstMuonsPt","FirstMuonPt",25,0,50,"First  #mu p_{T}, GeV","Events");
  SecondMuonsPt=HConfig.GetTH1D(Name+"_SecondMuonsPt","SecondMuonPt",25,0,50,"Second #mu p_{T}, GeV","Events");
  ThirdMuonsPt=HConfig.GetTH1D(Name+"_ThirdMuonsPt","ThirdMuonPt",25,0,50,"Third  #mu p_{T}, GeV","Events");


  FirstMuonsEta=HConfig.GetTH1D(Name+"_FirstMuonsEta","FirstMuonEta",25,-2.5,2.5,"First  #mu  rapidity","Events");
  SecondMuonsEta=HConfig.GetTH1D(Name+"_SecondMuonsEta","SecondMuonEta",25,-2.5,2.5,"Second #mu  rapidity","Events");
  ThirdMuonsEta=HConfig.GetTH1D(Name+"_ThirdMuonsEta","ThirdMuonEta",25,-2.5,2.5,"Third  #mu rapidity","Events");





  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  TwoMuTrack::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&MuonsPt);
  Extradist1d.push_back(&MuonsEta);
  Extradist1d.push_back(&MuonsPtRatio);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&TripleMass);
  Extradist1d.push_back(&Category);
  Extradist2d.push_back(&PhiMassVsDsMass);

  Extradist1d.push_back(&FirstMuonsPt);
  Extradist1d.push_back(&SecondMuonsPt);
  Extradist1d.push_back(&ThirdMuonsPt);

  Extradist1d.push_back(&FirstMuonsEta);
  Extradist1d.push_back(&SecondMuonsEta);
  Extradist1d.push_back(&ThirdMuonsEta);

}


void  TwoMuTrack::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  value.at(HLTOk) = 0;
  value.at(TwoMuTrackCandidate) = 0;
  if(Ntp->NTwoMuonsTrack()>0)  value.at(TwoMuTrackCandidate)=1;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  pass.at(TwoMuTrackCandidate)= (value.at(TwoMuTrackCandidate)==cut.at(TwoMuTrackCandidate)); 




  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);

    //    std::cout<<"N 2 mu + track  "<< Ntp->NTwoMuonsTrack() << " N three mu  "<< Ntp->NTwoMuTrackons() << "  Number of Vertices  " << Ntp->NumberOfSVertices() <<std::endl;
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

    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    }
  }
}



void  TwoMuTrack::Finish(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 

  if(mode == RECONSTRUCT){
    for(unsigned int i=0; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = 1/Nminus0.at(0).at(i).Integral();
      ScaleAllHistOfType(i,scale);
    }
  }


  Selection::Finish();

}





