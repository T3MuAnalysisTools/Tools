#include "SignalSelectorArun.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

SignalSelectorArun::SignalSelectorArun(TString Name_, TString id_):
   Selection(Name_,id_),
   tauMinMass_(1.73),
   tauMaxMass_(1.82),
   tauMinSideBand_(1.65),
   tauMaxSideBand_(2.02),
   tauMassResCutLow(0.007),
   tauMassResCutHigh(0.01),
   phiVetoSigma(0.03),
   omegaVetoSigma(0.03),
   bdt_cutA2_3glb(0.0178036087),
   bdt_cutA1_3glb(0.13855843498),
   bdt_cutA2_2glbTrk(-0.00931573728943),
   bdt_cutA1_2glbTrk(0.0900170166134),
   bdt_cutB2_3glb(0.0364416637001),
   bdt_cutB1_3glb(0.1669982497),
   bdt_cutB2_2glbTrk(0.00258006234741),
   bdt_cutB1_2glbTrk(0.117157105615),
   bdt_cutC2_3glb(0.032225350338),
   bdt_cutC1_3glb(0.121164828298),
   bdt_cutC2_2glbTrk(-0.0483397164297),
   bdt_cutC1_2glbTrk(0.0750508898802)
{
   // This is a class constructor;
   T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
}

SignalSelectorArun::~SignalSelectorArun(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SignalSelectorArun::Configure(){
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
      if(i==TriggerOk)          cut.at(TriggerOk)=1;
      if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
      if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
      if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
      if(i==MuonID)             cut.at(MuonID)=1;
      if(i==PVRefit) 			  cut.at(PVRefit)=1;
      if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
      if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
      if(i==TriggerMatchMu1)    cut.at(TriggerMatchMu1)=0.03;
      if(i==TriggerMatchMu2)    cut.at(TriggerMatchMu2)=0.03;
      if(i==TriggerMatchMu3)    cut.at(TriggerMatchMu3)=0.03;
      if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
      if(i==GenMatch)	  		  cut.at(GenMatch)=0.03;

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
      else if(i==Mu1PtCut){
         title.at(i)="$p_{T}(\\mu_{1}) >$";
         title.at(i)+=cut.at(Mu1PtCut);
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="$p_{T}(\\mu_{2}) >$";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      }
      else if(i==Mu3PtCut){
         title.at(i)="$p_{T}(\\mu_{3}) >$";
         title.at(i)+=cut.at(Mu3PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon3 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="gl,gl,(gl/trk)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==PVRefit){
         title.at(i)="PV refit valid";
         hlabel="Primary Vertex refit valid";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==PhiVeto){
         title.at(i)="$\\phi$ mass veto";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Phi mass Veto, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
      }
      else if(i==OmegaVeto){
         title.at(i)="$\\omega$ mass veto";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Omega mass veto, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu1){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{1}} <$ 0.03";
         hlabel="dR muon 1";
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");

         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu1_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu1_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu2){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{2}} <$ 0.03";
         hlabel="dR muon 2";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu2_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu2_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu3){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{3}} <$ 0.03";
         hlabel="dR muon 3";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu3_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu3_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TauMassCut){
         title.at(i)="$\\tau$ mass (sideband in data)";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");

         hlabel="three mu mass, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,80,1.4,2.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,80,1.4,2.2,hlabel,"Events"));
      }
      else if(i==GenMatch){
         title.at(i)="GEN matching";
         hlabel="GEN match";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
  
  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");
  
  TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SignalSelectorArun::Store_ExtraDist(){ 
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

void  SignalSelectorArun::doEvent(){ 

   TLorentzVector Tau_P4;
   TLorentzVector TauG_P4;
   TLorentzVector TauT_P4;
   TLorentzVector TauS_P4;

   value.at(TriggerOk)=0;
   value.at(SignalCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(Mu3PtCut)=0;
   value.at(MuonID)=0;
   value.at(TriggerMatchMu1)=1.0;
   value.at(TriggerMatchMu2)=1.0;
   value.at(TriggerMatchMu3)=1.0;
   value.at(PVRefit)=0;
   value.at(PhiVeto)=99.0;
   value.at(OmegaVeto)=99.0;
   value.at(TauMassCut)=0;
   value.at(GenMatch)=1;
   
   
   bool threeGlobal = false;
   unsigned int t;
   int id(Ntp->GetMCID());
   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
   bool HLTOk(false);
   bool L1Ok(false);
   for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);

      if(id==1){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }
      if( (id==1 && Ntp->WhichEra(2017).Contains("RunF")) || (id==1 && Ntp->WhichEra(2018).Contains("Run")) ){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }
      if(id!=1){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }

   }

   // ------------ HLT selection -----------------
   bool DoubleMuFired(0);
   bool TripleMuFired(0);
   for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
      TString L1TriggerName = Ntp->L1Name(il1);

      if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17"))      TripleMuFired = Ntp-> L1Decision(il1);
      }

      if(id!=1){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
      }
      if (id==1 && Ntp->WhichEra(2018).Contains("Run")){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))							DoubleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))	TripleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2"))								DoubleMuFired = Ntp->L1Decision(il1);
      }
      if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
      }
   }



   if(DoubleMuFired  or TripleMuFired) L1Ok = true;
   value.at(TriggerOk)=(HLTOk && L1Ok);
   pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));


   // ------------ Define new variables -----------------
   unsigned int  final_idx=0;

   double min_chi2(100.);
   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_tau = 99.0;


   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_tau = 0.;
   double sumPtTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_tau = 0;
   // ---------------------------------------------------

   // ------------ Chose candidate with minimum chi square -----------------
   for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
      if(Ntp->Vertex_Signal_KF_Chi2(i_idx,false) < min_chi2){
         min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx,false);
         final_idx = i_idx;
      }
   }
   // ---------------------------------------------------


   NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);
   value.at(SignalCandidate) = Ntp->NThreeMuons();

   if(Ntp->NThreeMuons()>0){
      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
      //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
      //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
      //
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) && 
            Ntp->Muon_isGlobalMuon(Muon_index_2) &&
            ( Ntp->Muon_isGlobalMuon(Muon_index_3) || Ntp->Muon_isTrackerMuon(Muon_index_3)));	 

      if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal = true;
      //------------------------------------------------------------------------------------------------------

      value.at(PVRefit) = (Ntp->Vertex_RefitPVisValid(final_idx)==1);
      value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();

      vector<unsigned int> idx_vec;

      idx_vec.push_back(mu1_idx);
      idx_vec.push_back(mu2_idx);
      idx_vec.push_back(mu3_idx);

      unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
      unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
      unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

      TLorentzVector TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);
      TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);


      double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
      double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

      value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
      value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

      value.at(TriggerMatchMu1) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(0);
      value.at(TriggerMatchMu2) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(1);
      value.at(TriggerMatchMu3) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(2);
      //value.at(TauMassCut) = TauRefittedLV.M();
      value.at(TauMassCut) = TauLV.M();
   }
   if (id!=1) value.at(GenMatch) = Ntp->TauGenMatch(final_idx);
   else value.at(GenMatch) = 0;

   pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
   pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
   pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
   pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
   pass.at(PVRefit) = (value.at(PVRefit) == cut.at(PVRefit));
   pass.at(TriggerMatchMu1) = (value.at(TriggerMatchMu1)  <  cut.at(TriggerMatchMu1));
   pass.at(TriggerMatchMu2) = (value.at(TriggerMatchMu2)  <  cut.at(TriggerMatchMu2));
   pass.at(TriggerMatchMu3) = (value.at(TriggerMatchMu3)  <  cut.at(TriggerMatchMu3));
   //pass.at(TriggerMatchMu1) = true;
   //pass.at(TriggerMatchMu2) = true;
   //pass.at(TriggerMatchMu3) = true;
   pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 8*PDG_Var::Phi_width());
   pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 3*PDG_Var::Omega_width()); 

   if(id!=1) pass.at(TauMassCut) = true;
   else  pass.at(TauMassCut) = (value.at(TauMassCut)>tauMinSideBand_  && value.at(TauMassCut) < tauMaxSideBand_);


   //pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
   pass.at(GenMatch) = true;
   std::vector<unsigned int> exclude_cuts;
   exclude_cuts.push_back(PhiVeto);
   exclude_cuts.push_back(OmegaVeto);

   if(passAllBut(exclude_cuts)){


      TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);

      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
      vector<unsigned int> idx_vec;

      idx_vec.push_back(mu1_idx);
      idx_vec.push_back(mu2_idx);
      idx_vec.push_back(mu3_idx);

      unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
      unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
      unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

      double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
      double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

      double pmass  = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 

      TauMass_all_nophiVeto.at(t).Fill(TauRefittedLV.M(),pmass,1);
   }
   
   
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Here the index t belongs to sample type, this value originally filled at Ntuple filling
  // level: https://github.com/T3MuAnalysisTools/DsTau23Mu/blob/master/T3MNtuple/interface/DataMCType.h
  // but you can flexibly redefined this and make combinations like, Data, MC1, MC2, MC3+MC4, etc ...
  
  // Apply Selection


  value.at(PrimeVtx)=Ntp->NVtx(); // Here the actual_value of a cut is set
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); // Here we check that the actuall value of PrimeVrtices is above 5.
  
  if(Ntp->NThreeMuons()>0){  // Check if this is a signal category (take the first triplet only in this example)

    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon 0
    value.at(LeadingMuonPt) = Ntp->Muon_P4(mu1_pt_idx).Pt();
    
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);  // leading pT muon 1
    value.at(LeadingMuonPt1) = Ntp->Muon_P4(mu2_pt_idx).Pt();
    
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);  // leading pT muon 2
    value.at(LeadingMuonPt2) = Ntp->Muon_P4(mu3_pt_idx).Pt();
    
    Tau_P4=Ntp->Muon_P4(mu1_pt_idx)+Ntp->Muon_P4(mu2_pt_idx)+Ntp->Muon_P4(mu3_pt_idx);

  }

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


void  SignalSelectorArun::Finish(){
  Selection::Finish();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





