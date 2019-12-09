#include "TriggerStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

TriggerStudy::TriggerStudy(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.731),
  tauMaxMass_(1.823),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

TriggerStudy::~TriggerStudy(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  TriggerStudy::Configure(){
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
	 if(i==TriggerOk){
      title.at(i)="Pass HLT";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
	 else if(i==SignalCandidate){
      title.at(i)="is signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  TriggerStudy::Store_ExtraDist(){ 
   cout<<"Nothing to store"<<endl;
}


void  TriggerStudy::doEvent(){ 
  
  unsigned int t;
  bool doubleMu_flag = false;
  bool singleMu_flag = false;
  bool muon_flag = false;

  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if((HLT.Contains("DoubleMu3_Trk_Tau3mu_v") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v"))) {
		if (Ntp->HLTDecision(iTrigger)) doubleMu_flag = true;
		}
	if(HLT.Contains("HLT_Mu8_v")) {
		if (Ntp->HLTDecision(iTrigger)) singleMu_flag = true;
		}
	 }
  
  for (unsigned int iMuon=0; iMuon<Ntp->NMuons()-1; iMuon++){
	for (unsigned int jMuon=iMuon+1; jMuon<Ntp->NMuons(); jMuon++){
  	// Mu_pt8_eta1p5_dxySig5
	/*
	if ( abs(Ntp->Muon_simPdgId(iMuon))!=13 ) continue;
	if ( Ntp->Muon_P4(iMuon).Pt()<=8.0 ) continue;
	if ( abs(Ntp->Muon_P4(iMuon).Eta()>=1.5) ) continue;
	if ( Ntp->Muon_innerTrack_numberofValidHits(iMuon)<1 ) continue;
   if ( Ntp->Muon_dxyError(iMuon)==0) continue;
	if ( (Ntp->Muon_dxy_beamSpot(iMuon))/(Ntp->Muon_dxyError(iMuon))<=5) continue;
	*/
   if (Ntp->Muon_simPdgId(iMuon)!=13 || Ntp->Muon_simPdgId(jMuon)==13) continue;
   if (!((Ntp->Muon_P4(iMuon).Pt()>3.0 && Ntp->Muon_P4(jMuon).Pt()>2.5) || (Ntp->Muon_P4(iMuon).Pt()>2.5 && Ntp->Muon_P4(jMuon).Pt()>3.0))) continue; // find a muon pair (3,2.5) GeV
   muon_flag = true;
	}
  }
  
  value.at(TriggerOk) = (doubleMu_flag || singleMu_flag);
  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));

  value.at(SignalCandidate) = Ntp->NThreeMuons();
  pass.at(SignalCandidate) = (value.at(SignalCandidate) == cut.at(SignalCandidate));
  
  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);

  if(status && muon_flag){
		if (doubleMu_flag) n_doubleMu++;
  		if (singleMu_flag) n_singleMu++;
  		if (doubleMu_flag && singleMu_flag) n_overlap++;
	 }   
}

void  TriggerStudy::Finish(){

  cout<<"[TriggerStudy]: Number of events passing the double mu trigger = "<<n_doubleMu<<endl;
  cout<<"[TriggerStudy]: Number of events passing the single mu trigger = "<<n_singleMu<<endl;
  cout<<"[TriggerStudy]: Number of events passing both triggers = "<<n_overlap<<endl;

  if(mode == RECONSTRUCT){
    
    int id(Ntp->GetMCID());
    double scale(1.);
    double scaleDsTau(0.637);
    double scaleBpTau(0.262);
    double scaleB0Tau(0.099);

    if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    ScaleAllHistOfType(2,scale*scaleDsTau);

    if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    ScaleAllHistOfType(3,scale*scaleB0Tau);

    if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    ScaleAllHistOfType(4,scale*scaleBpTau);

  }

    Selection::Finish();
}





