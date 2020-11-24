#include "MyNewAnalysis.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

MyNewAnalysis::MyNewAnalysis(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

using namespace std;
TLorentzVector Tau_P4;
TLorentzVector TauG_P4;
TLorentzVector TauT_P4;
TLorentzVector TauS_P4;

MyNewAnalysis::~MyNewAnalysis(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  MyNewAnalysis::Configure(){
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
    if(i==TriggerOk)        cut.at(TriggerOk)=1;
    if(i==PrimeVtx)         cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
    if(i==SignalCandidate)  cut.at(SignalCandidate)=1;
    if(i==LeadingMuonPt)    cut.at(LeadingMuonPt)=2.5;
    if(i==LeadingMuonPt1)   cut.at(LeadingMuonPt1)=2.5;
    if(i==LeadingMuonPt2)   cut.at(LeadingMuonPt2)=2.5;

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
    else if(i==SignalCandidate){
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==LeadingMuonPt){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==LeadingMuonPt1){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt1_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt1_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==LeadingMuonPt2){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt2_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt2_",htitle,80,0,20,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams

  LeadMuonPt=HConfig.GetTH1D(Name+"_LeadMuonPt","LeadMuonPt",40,0,20,"p_{T}(#mu_{1}), GeV","Events");
  LeadMuonEta=HConfig.GetTH1D(Name+"_LeadMuonEta","LeadMuonEta",40,-2.6,2.6,"#eta(#mu_{1})","Events");
  LeadMuonPhi=HConfig.GetTH1D(Name+"_LeadMuonPhi","LeadMuonPhi",40,-3.15,3.15,"#phi(#mu_{1})","Events");
  
  LeadMuonPt1=HConfig.GetTH1D(Name+"_LeadMuonPt1","LeadMuonPt1",40,0,20,"p_{T}(#mu_{2}), GeV","Events");
  LeadMuonEta1=HConfig.GetTH1D(Name+"_LeadMuonEta1","LeadMuonEta1",40,-2.6,2.6,"#eta(#mu_{2})","Events");
  LeadMuonPhi1=HConfig.GetTH1D(Name+"_LeadMuonPhi1","LeadMuonPhi1",40,-3.15,3.15,"#phi(#mu_{2})","Events");
  
  LeadMuonPt2=HConfig.GetTH1D(Name+"_LeadMuonPt2","LeadMuonPt2",40,0,20,"p_{T}(#mu_{3}), GeV","Events");
  LeadMuonEta2=HConfig.GetTH1D(Name+"_LeadMuonEta2","LeadMuonEta2",40,-2.6,2.6,"#eta(#mu_{3})","Events");
  LeadMuonPhi2=HConfig.GetTH1D(Name+"_LeadMuonPhi2","LeadMuonPhi2",40,-3.15,3.15,"#phi(#mu_{3})","Events");
  
  TauPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",40,0,20,"p_{T}(#tau_{1}), GeV","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",40,-2.6,2.6,"#eta(#tau_{1})","Events");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","TauPhi",40,-3.15,3.15,"#phi(#tau_{1})","Events");
  
  InvMu1G=HConfig.GetTH1D(Name+"_InvMu1G","First muon is global",40,0,0.5,"LEading Muon is Global","Events");
  InvMu1T=HConfig.GetTH1D(Name+"_InvMu1T","First muon is tracker",40,0,0.5,"LEading Muon is Tracker","Events");
  InvMu1S=HConfig.GetTH1D(Name+"_InvMu1S","First muon is standalone",40,0,0.5,"LEading Muon is Standalone","Events");
  
  InvMu2G=HConfig.GetTH1D(Name+"_InvMu2G","Second muon is global",40,0,0.5,"LEading Muon 1 is Global","Events");
  InvMu2T=HConfig.GetTH1D(Name+"_InvMu2T","Second muon is tracker",40,0,0.5,"LEading Muon 1 is Tracker","Events");
  InvMu2S=HConfig.GetTH1D(Name+"_InvMu2S","Second muon is standalone",40,0,0.5,"LEading Muon 1 is Standalone","Events");
  
  InvMu3G=HConfig.GetTH1D(Name+"_InvMu3G","Third muon is global",40,0,0.5,"LEading Muon 2 is Global","Events");
  InvMu3T=HConfig.GetTH1D(Name+"_InvMu3T","Third muon is tracker",40,0,0.5,"LEading Muon 3 is Tracker","Events");
  InvMu3S=HConfig.GetTH1D(Name+"_InvMu3S","Third muon is standalone",40,0,0.5,"LEading Muon 3 is Standalone","Events");
  
  InvTG=HConfig.GetTH1D(Name+"_InvTG","Tau is global",40,0,3,"Tau is Global","Events");
  InvTT=HConfig.GetTH1D(Name+"_InvTT","Tau is tracker",40,0,3,"Tau is Tracker","Events");
  InvTS=HConfig.GetTH1D(Name+"_InvTS","Tau is standalone",40,0,3,"Tau is Standalone","Events");  


  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  MyNewAnalysis::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&LeadMuonPt);
  Extradist1d.push_back(&LeadMuonEta);
  Extradist1d.push_back(&LeadMuonPhi);
  
  Extradist1d.push_back(&LeadMuonPt1);
  Extradist1d.push_back(&LeadMuonEta1);
  Extradist1d.push_back(&LeadMuonPhi1);
  
  Extradist1d.push_back(&LeadMuonPt2);
  Extradist1d.push_back(&LeadMuonEta2);
  Extradist1d.push_back(&LeadMuonPhi2);
  
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPhi);
  
  Extradist1d.push_back(&InvMu1G);
  Extradist1d.push_back(&InvMu1T);
  Extradist1d.push_back(&InvMu1S);
  
  Extradist1d.push_back(&InvMu2G);
  Extradist1d.push_back(&InvMu2T);
  Extradist1d.push_back(&InvMu2S);
  
  Extradist1d.push_back(&InvMu3G);
  Extradist1d.push_back(&InvMu3T);
  Extradist1d.push_back(&InvMu3S);
  
  Extradist1d.push_back(&InvTG);
  Extradist1d.push_back(&InvTT);
  Extradist1d.push_back(&InvTS);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  MyNewAnalysis::doEvent(){ 
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


  unsigned int  final_idx=0;


  value.at(SignalCandidate) = Ntp->NThreeMuons();
  if(Ntp->NThreeMuons()>0){  // Check if this is a signal category (take the first triplet only in this example)

    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon 0
    value.at(LeadingMuonPt) = Ntp->Muon_P4(mu1_pt_idx).Pt();
    
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);  // leading pT muon 1
    value.at(LeadingMuonPt1) = Ntp->Muon_P4(mu2_pt_idx).Pt();
    
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);  // leading pT muon 2
    value.at(LeadingMuonPt2) = Ntp->Muon_P4(mu3_pt_idx).Pt();
    
    Tau_P4=Ntp->Muon_P4(mu1_pt_idx)+Ntp->Muon_P4(mu2_pt_idx)+Ntp->Muon_P4(mu3_pt_idx);

  }


  pass.at(SignalCandidate) = (value.at(SignalCandidate)  > 0 );
  pass.at(LeadingMuonPt)   = (value.at(LeadingMuonPt)    > cut.at(LeadingMuonPt));
  pass.at(LeadingMuonPt1)   = (value.at(LeadingMuonPt1)    > cut.at(LeadingMuonPt1));
  pass.at(LeadingMuonPt2)   = (value.at(LeadingMuonPt2)    > cut.at(LeadingMuonPt2));


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
    // Lets plot the pT, phi and eta of the leading muon
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon 0
    LeadMuonPt.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Pt(),1);
    LeadMuonEta.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Eta(),1);
    LeadMuonPhi.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Phi(),1);
    
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);  // leading pT muon 1
    LeadMuonPt1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Pt(),1);
    LeadMuonEta1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Eta(),1);
    LeadMuonPhi1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Phi(),1);
    
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);  // leading pT muon 2
    LeadMuonPt2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Pt(),1);
    LeadMuonEta2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Eta(),1);
    LeadMuonPhi2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Phi(),1); 
    
    TauPt.at(t).Fill(Tau_P4.Pt(),1);
    TauEta.at(t).Fill(Tau_P4.Eta(),1);
    TauPhi.at(t).Fill(Tau_P4.Phi(),1);
    
    if(Ntp->Muon_isGlobalMuon(mu1_pt_idx)){
    InvMu1G.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu1_pt_idx)){
    InvMu1T.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu1_pt_idx)){
    InvMu1S.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu2_pt_idx)){
    InvMu2G.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu2_pt_idx)){
    InvMu2T.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu2_pt_idx)){
    InvMu2S.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu3_pt_idx)){
    InvMu3G.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu3_pt_idx)){
    InvMu3T.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu3_pt_idx)){
    InvMu3S.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu1_pt_idx)&&Ntp->Muon_isGlobalMuon(mu2_pt_idx)&&Ntp->Muon_isGlobalMuon(mu3_pt_idx)){
    InvTG.at(t).Fill(Tau_P4.M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu1_pt_idx)&&Ntp->Muon_isTrackerMuon(mu2_pt_idx)&&Ntp->Muon_isTrackerMuon(mu3_pt_idx)){
    InvTT.at(t).Fill(Tau_P4.M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu1_pt_idx)&&Ntp->Muon_isStandAloneMuon(mu2_pt_idx)&&Ntp->Muon_isStandAloneMuon(mu3_pt_idx)){
    InvTS.at(t).Fill(Tau_P4.M(),1);
    }


  }
}


void  MyNewAnalysis::Finish(){
  Selection::Finish();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





