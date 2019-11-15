#include "MuIDStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "SkimConfig.h"

using namespace std;

MuIDStudy::MuIDStudy(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

MuIDStudy::~MuIDStudy(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  MuIDStudy::Configure(){


  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    //    if(i==VertChi2)           cut.at(VertChi2)=20.0;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=0.03;
    if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==TriggerOk){
      title.at(i)="Pass HLT and L1";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    /*  else if(i==VertChi2){
      title.at(i)="Triple Vertex Chi Squared $<$ 15";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Chi squared of triple vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VertChi2_",htitle,50,0,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VertChi2_",htitle,50,0,25,hlabel,"Events"));
      }*/


    else if(i==Mu1PtCut){
      title.at(i)="$p_{T}(\\mu_{1}) >$ 3.0 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Muon1 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
    }

    else if(i==Mu2PtCut){
      title.at(i)="$p_{T}(\\mu_{2}) >$ 3.0 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Muon2 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
    }

    else if(i==Mu3PtCut){
      title.at(i)="$p_{T}(\\mu_{3}) >$ 2.0 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Muon3 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
    }

    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="gl,gl,tr";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
      hlabel="Rho mass veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
    }

    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
    }

    else if(i==TauMassCut){
      title.at(i)="$\\tau$ mass (sideband in data)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,60,2.1,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,60,1.4,2.2,hlabel,"Events"));
    }
  } 


  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms

  Muon1isGlob =HConfig.GetTH1D(Name+"_Muon1isGlob","Muon1isGlob",2,-0.5,1.5,"  #mu_{1} is global muon","Events");
  Muon2isGlob =HConfig.GetTH1D(Name+"_Muon2isGlob","Muon2isGlob",2,-0.5,1.5,"  #mu_{2} is global muon","Events");
  Muon3isGlob =HConfig.GetTH1D(Name+"_Muon3isGlob","Muon3isGlob",2,-0.5,1.5,"  #mu_{3} is global muon","Events");

  Muon1isStand =HConfig.GetTH1D(Name+"_Muon1isStand","Muon1isStand",2,-0.5,1.5,"  #mu_{1} is a standalone muon","Events");
  Muon2isStand =HConfig.GetTH1D(Name+"_Muon2isStand","Muon2isStand",2,-0.5,1.5,"  #mu_{2} is a standalone muon","Events");
  Muon3isStand =HConfig.GetTH1D(Name+"_Muon3isStand","Muon3isStand",2,-0.5,1.5,"  #mu_{3} is a standalone muon","Events");

  Muon1isTrack =HConfig.GetTH1D(Name+"_Muon1isTrack","Muon1isTrack",2,-0.5,1.5,"  #mu_{1} is tracker muon","Events");
  Muon2isTrack =HConfig.GetTH1D(Name+"_Muon2isTrack","Muon2isTrack",2,-0.5,1.5,"  #mu_{2} is tracker muon","Events");
  Muon3isTrack =HConfig.GetTH1D(Name+"_Muon3isTrack","Muon3isTrack",2,-0.5,1.5,"  #mu_{3} is tracker muon","Events");

  Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
  Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
  Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");

  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#eta(#mu_{1})","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#eta(#mu_{2})","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#eta(#mu_{3})","Events");

  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
  TauP  =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");

  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

  EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");

  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");


  CustMuonIdMC1=HConfig.GetTH1D(Name+"_CustMuonIdMC1","MCCustMuonIdMC1",6,0,6,"#mu_{1} MuonId","Events");
  CustMuonIdData1=HConfig.GetTH1D(Name+"_CustMuonIdData1","CustMuonIdData1",6,0,6,"#mu_{1} MuonId","Events");

  CustMuonIdMC2=HConfig.GetTH1D(Name+"_CustMuonIdMC2","MCCustMuonIdMC2",6,0,6,"#mu_{2} MuonId","Events");
  CustMuonIdData2=HConfig.GetTH1D(Name+"_CustMuonIdData2","CustMuonIdData2",6,0,6,"#mu_{2} MuonId","Events");

  CustMuonIdMC3=HConfig.GetTH1D(Name+"_CustMuonIdMC3","MCCustMuonIdMC3",6,0,6,"#mu_{3} MuonId","Events");
  CustMuonIdData3=HConfig.GetTH1D(Name+"_CustMuonIdData3","CustMuonIdData3",6,0,6,"#mu_{3} MuonId","Events");




  CustMuonIdMC1A=HConfig.GetTH1D(Name+"_CustMuonIdMC1A","MCCustMuonIdMC1A",6,0,6,"#mu_{1} MuonId","Events");
  CustMuonIdData1A=HConfig.GetTH1D(Name+"_CustMuonIdData1A","CustMuonIdData1A",6,0,6,"#mu_{1} MuonId","Events");

  CustMuonIdMC2A=HConfig.GetTH1D(Name+"_CustMuonIdMC2A","MCCustMuonIdMC2A",6,0,6,"#mu_{2} MuonId","Events");
  CustMuonIdData2A=HConfig.GetTH1D(Name+"_CustMuonIdData2A","CustMuonIdData2A",6,0,6,"#mu_{2} MuonId","Events");

  CustMuonIdMC3A=HConfig.GetTH1D(Name+"_CustMuonIdMC3A","MCCustMuonIdMC3A",6,0,6,"#mu_{3} MuonId","Events");
  CustMuonIdData3A=HConfig.GetTH1D(Name+"_CustMuonIdData3A","CustMuonIdData3A",6,0,6,"#mu_{3} MuonId","Events");






  CustMuonIdMC1B=HConfig.GetTH1D(Name+"_CustMuonIdMC1B","MCCustMuonIdMC1B",6,0,6,"#mu_{1} MuonId","Events");
  CustMuonIdData1B=HConfig.GetTH1D(Name+"_CustMuonIdData1B","CustMuonIdData1B",6,0,6,"#mu_{1} MuonId","Events");

  CustMuonIdMC2B=HConfig.GetTH1D(Name+"_CustMuonIdMC2B","MCCustMuonIdMC2B",6,0,6,"#mu_{2} MuonId","Events");
  CustMuonIdData2B=HConfig.GetTH1D(Name+"_CustMuonIdData2B","CustMuonIdData2B",6,0,6,"#mu_{2} MuonId","Events");

  CustMuonIdMC3B=HConfig.GetTH1D(Name+"_CustMuonIdMC3B","MCCustMuonIdMC3B",6,0,6,"#mu_{3} MuonId","Events");
  CustMuonIdData3B=HConfig.GetTH1D(Name+"_CustMuonIdData3B","CustMuonIdData3B",6,0,6,"#mu_{3} MuonId","Events");





  CustMuonIdMC1C=HConfig.GetTH1D(Name+"_CustMuonIdMC1C","MCCustMuonIdMC1C",6,0,6,"#mu_{1} MuonId","Events");
  CustMuonIdData1C=HConfig.GetTH1D(Name+"_CustMuonIdData1C","CustMuonIdData1C",6,0,6,"#mu_{1} MuonId","Events");

  CustMuonIdMC2C=HConfig.GetTH1D(Name+"_CustMuonIdMC2C","MCCustMuonIdMC2C",6,0,6,"#mu_{2} MuonId","Events");
  CustMuonIdData2C=HConfig.GetTH1D(Name+"_CustMuonIdData2C","CustMuonIdData2C",6,0,6,"#mu_{2} MuonId","Events");

  CustMuonIdMC3C=HConfig.GetTH1D(Name+"_CustMuonIdMC3C","MCCustMuonIdMC3C",6,0,6,"#mu_{3} MuonId","Events");
  CustMuonIdData3C=HConfig.GetTH1D(Name+"_CustMuonIdData3C","CustMuonIdData3C",6,0,6,"#mu_{3} MuonId","Events");




  Significance3C=HConfig.GetTH1D(Name+"_Significance3C","Significance3C",6,0,6,"Category C #mu_{3} MuonId","S/#sqrt{B}");
  Significance2C=HConfig.GetTH1D(Name+"_Significance2C","Significance2C",6,0,6,"Category C #mu_{2} MuonId","S/#sqrt{B}");
  Significance1C=HConfig.GetTH1D(Name+"_Significance1C","Significance1C",6,0,6,"Category C #mu_{1} MuonId","S/#sqrt{B}");


  Significance3B=HConfig.GetTH1D(Name+"_Significance3B","Significance3B",6,0,6,"Category B #mu_{3} MuonId","S/#sqrt{B}");
  Significance2B=HConfig.GetTH1D(Name+"_Significance2B","Significance2B",6,0,6,"Category B #mu_{2} MuonId","S/#sqrt{B}");
  Significance1B=HConfig.GetTH1D(Name+"_Significance1B","Significance1B",6,0,6,"Category B #mu_{1} MuonId","S/#sqrt{B}");


  Significance3A=HConfig.GetTH1D(Name+"_Significance3A","Significance3A",6,0,6,"Category A #mu_{3} MuonId","S/#sqrt{B}");
  Significance2A=HConfig.GetTH1D(Name+"_Significance2A","Significance2A",6,0,6,"Category A #mu_{2} MuonId","S/#sqrt{B}");
  Significance1A=HConfig.GetTH1D(Name+"_Significance1A","Significance1A",6,0,6,"Category A #mu_{1} MuonId","S/#sqrt{B}");


  Significance3=HConfig.GetTH1D(Name+"_Significance3","Significance3",6,0,6,"#mu_{3} MuonId","S/#sqrt{B}");
  Significance2=HConfig.GetTH1D(Name+"_Significance2","Significance2",6,0,6,"#mu_{2} MuonId","S/#sqrt{B}");
  Significance1=HConfig.GetTH1D(Name+"_Significance1","Significance1",6,0,6,"#mu_{1} MuonId","S/#sqrt{B}");





  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");

  L1Triggers=HConfig.GetTH2D(Name+"_L1Triggers","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  L1TriggersB=HConfig.GetTH2D(Name+"_L1TriggersB","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersC=HConfig.GetTH2D(Name+"_L1TriggersC","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersD=HConfig.GetTH2D(Name+"_L1TriggersD","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersE=HConfig.GetTH2D(Name+"_L1TriggersE","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersF=HConfig.GetTH2D(Name+"_L1TriggersF","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  MuIDStudy::Store_ExtraDist(){ 

  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);
  Extradist1d.push_back(&Significance3C);
  Extradist1d.push_back(&Significance2C);
  Extradist1d.push_back(&Significance1C);


  Extradist1d.push_back(&Significance3B);
  Extradist1d.push_back(&Significance2B);
  Extradist1d.push_back(&Significance1B);



  Extradist1d.push_back(&Significance3A);
  Extradist1d.push_back(&Significance2A);
  Extradist1d.push_back(&Significance1A);


  Extradist1d.push_back(&Significance3);
  Extradist1d.push_back(&Significance2);
  Extradist1d.push_back(&Significance1);
  //Extradist1d.push_back(&Muon1isGlob);
  //Extradist1d.push_back(&Muon2isGlob);
  //Extradist1d.push_back(&Muon3isGlob);


  //Extradist1d.push_back(&Muon1isStand);
  //Extradist1d.push_back(&Muon2isStand);
  //Extradist1d.push_back(&Muon3isStand);
    

  //Extradist1d.push_back(&Muon1isTrack);
  //Extradist1d.push_back(&Muon2isTrack);
  //Extradist1d.push_back(&Muon3isTrack);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  // Extradist1d.push_back(&TauP);

  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist2d.push_back(&EMR_tau_eta);
  Extradist1d.push_back(&SVPVTauDirAngle);

  //  Extradist1d.push_back(&Muon1DRToTruth);
  //  Extradist1d.push_back(&Muon2DRToTruth);
  //  Extradist1d.push_back(&Muon3DRToTruth);


  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);


  // Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&EventMassResolution_PtEtaPhi);

  //Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&NSignalCandidates);
  Extradist2d.push_back(&L1Triggers);

  Extradist2d.push_back(&L1TriggersB);
  Extradist2d.push_back(&L1TriggersC);
  Extradist2d.push_back(&L1TriggersD);
  Extradist2d.push_back(&L1TriggersE);
  Extradist2d.push_back(&L1TriggersF);

  Extradist1d.push_back(&CustMuonIdMC1);
  Extradist1d.push_back(&CustMuonIdData1);

  Extradist1d.push_back(&CustMuonIdMC2);
  Extradist1d.push_back(&CustMuonIdData2);

  Extradist1d.push_back(&CustMuonIdMC3);
  Extradist1d.push_back(&CustMuonIdData3);




  Extradist1d.push_back(&CustMuonIdMC1A);
  Extradist1d.push_back(&CustMuonIdData1A);

  Extradist1d.push_back(&CustMuonIdMC2A);
  Extradist1d.push_back(&CustMuonIdData2A);

  Extradist1d.push_back(&CustMuonIdMC3A);
  Extradist1d.push_back(&CustMuonIdData3A);



  Extradist1d.push_back(&CustMuonIdMC1B);
  Extradist1d.push_back(&CustMuonIdData1B);

  Extradist1d.push_back(&CustMuonIdMC2B);
  Extradist1d.push_back(&CustMuonIdData2B);

  Extradist1d.push_back(&CustMuonIdMC3B);
  Extradist1d.push_back(&CustMuonIdData3B);



  Extradist1d.push_back(&CustMuonIdMC1C);
  Extradist1d.push_back(&CustMuonIdData1C);

  Extradist1d.push_back(&CustMuonIdMC2C);
  Extradist1d.push_back(&CustMuonIdData2C);

  Extradist1d.push_back(&CustMuonIdMC3C);
  Extradist1d.push_back(&CustMuonIdData3C);

}


void  MuIDStudy::doEvent(){ 

  
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

    //    if(id==1 && Ntp->WhichEra(2017).Contains("RunF")){
    //      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    //    }	 

    if(id!=1){
      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    }
    
  }


  bool DoubleMuFired(0);
  bool TripleMuFired(0);
  for(unsigned int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17"))      TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }
  }

  if(DoubleMuFired  or TripleMuFired) L1Ok = true;
  value.at(TriggerOk)=(HLTOk and L1Ok);
  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));

  value.at(SignalCandidate)=0;
  unsigned int  signal_idx=0;
  value.at(TriggerMatch)=0;

  double min_chi2(99.);
  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }


  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();
    //    value.at(VertChi2) = Ntp->Vertex_Signal_KF_Chi2(signal_idx);
    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
    //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
    //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);
    //
    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
    			(Ntp->Muon_isGlobalMuon(mu3_pt_idx) or Ntp->Muon_isTrackerMuon(mu3_pt_idx)));
    //------------------------------------------------------------------------------------------------------
  
    value.at(Mu1PtCut) = Ntp->Muon_P4(mu1_pt_idx).Pt();
    value.at(Mu2PtCut) = Ntp->Muon_P4(mu2_pt_idx).Pt();
    value.at(Mu3PtCut) = Ntp->Muon_P4(mu3_pt_idx).Pt();
  
    vector<unsigned int> idx_vec;
    
    idx_vec.push_back(mu1_idx);
    idx_vec.push_back(mu2_idx);
    idx_vec.push_back(mu3_idx);

    unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    TLorentzVector TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);
    TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);


    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

    value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx)){
      value.at(TriggerMatch)+=i; 
    }
    
    //    value.at(TauMassCut) = TauLV.M();
    value.at(TauMassCut) = TauRefittedLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  //  pass.at(VertChi2) = (value.at(VertChi2) <= cut.at(VertChi2));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   = true;//(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = (value.at(TriggerMatch)  <  cut.at(TriggerMatch));
  pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 8*PDG_Var::Phi_width());
  pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 3*PDG_Var::Omega_width());

  if(id!=1) pass.at(TauMassCut) = true;
  else  pass.at(TauMassCut) =((value.at(TauMassCut) < tauMinMass_ && value.at(TauMassCut) > tauMinSideBand_) || (value.at(TauMassCut)> tauMaxMass_ && value.at(TauMassCut) < tauMaxSideBand_) );


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);

  if(status){
    NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);


    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);


    unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    std::vector<unsigned int> EtaSortedIndices;
    
    EtaSortedIndices.push_back(Muon_Eta_index_1);
    EtaSortedIndices.push_back(Muon_Eta_index_2);
    EtaSortedIndices.push_back(Muon_Eta_index_3);


    if(id==1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdData1.at(t).Fill(iMuId,1);

      }
    }
    if(id !=1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdMC1.at(t).Fill(iMuId,1);
      }
    }

    if(id==1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdData2.at(t).Fill(iMuId,1);
      }
    }
    if(id !=1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdMC2.at(t).Fill(iMuId,1);
      }
    }

    if(id==1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdData3.at(t).Fill(iMuId,1);
      }
    }
    if(id !=1){
      for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdMC3.at(t).Fill(iMuId,1);
      }
    }




    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007){

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdData1A.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdMC1A.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdData2A.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdMC2A.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdData3A.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdMC3A.at(t).Fill(iMuId,1);
	}
      }
      

    }
    

    
    
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.01){

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdData1B.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdMC1B.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdData2B.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdMC2B.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdData3B.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdMC3B.at(t).Fill(iMuId,1);
	}
      }
      
    }
    

    
    
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01){

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdData1C.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_1).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_1).at(iMuId) == 1) CustMuonIdMC1C.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdData2C.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_2).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_2).at(iMuId) == 1) CustMuonIdMC2C.at(t).Fill(iMuId,1);
	}
      }

      if(id==1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdData3C.at(t).Fill(iMuId,1);
	}
      }
      if(id !=1){
	for(unsigned int iMuId =0 ; iMuId < Ntp->MuonCustomID(Muon_index_3).size(); iMuId++){
	  if( Ntp-> MuonCustomID(Muon_index_3).at(iMuId) == 1) CustMuonIdMC3C.at(t).Fill(iMuId,1);
	}
      }
    }


    /*    std::cout<<"is Loose    "<< isLoose  <<"    BitMask   "  << Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_1) , 
	  Ntp->MuonQualityBitMask::Bit_MuonTight) 
	  <<"   SS   "  << Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdLoose) 
	  <<"   cut based tight   " <<  Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdTight)
	  <<"   cut based medium  " <<  Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdMedium) 
	  <<"   cut based soft    " <<  Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::SoftCutBasedId) 
	  <<"   cut based loose   " <<  Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdLoose) << std::endl;*/
    

    
    
    EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);
    

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    vector<unsigned int> idx_vec;

    idx_vec.push_back(Muon_index_1);
    idx_vec.push_back(Muon_index_2);
    idx_vec.push_back(Muon_index_3);

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);


    TLorentzVector MuonOS = Ntp->Muon_P4(os_mu_idx);
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);


    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);



    if(id==1){
      if((TauRefitLV.M() > 1.65 && TauRefitLV.M() < 1.73 )  ||   (TauRefitLV.M() > 1.82 && TauRefitLV.M() < 1.9)){
	L1Triggers.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunB"))L1TriggersB.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunC"))L1TriggersC.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunD"))L1TriggersD.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunE"))L1TriggersE.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunF"))L1TriggersF.at(t).Fill(DoubleMuFired,TripleMuFired);
      }
    }

    if(id!=1)L1Triggers.at(t).Fill(DoubleMuFired,TripleMuFired);




    Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
    Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
    Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);

    Muon1Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Eta(),1);
    Muon2Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Eta(),1);
    Muon3Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Eta(),1);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    TauEta.at(t).Fill(TauLV.Eta(),1);
    TauPt.at(t).Fill(TauLV.Pt(),1);
    TauP.at(t).Fill(TauLV.P(),1);

    EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());
    //    TauMass.at(t).Fill(TauLV.M(),1);
    //    TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    

    Muon1isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_1),1);
    Muon2isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_2),1);
    Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);


    Muon1isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_1),w);
    Muon2isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_2),w);
    Muon3isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_3),w);


    Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
    Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
    Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);


    TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(0),1);
    TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(1),1);
    TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(2),1);

    VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(signal_idx),w);
    FLSignificance.at(t).Fill(sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(signal_idx),Ntp->Vertex_PrimaryVertex_Covariance(signal_idx),
								   Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx))),w);
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    SVPVTauDirAngle.at(t).Fill(SVPV.Angle(TauLV.Vect()),w);


    var_vertexKFChi2 =Ntp->Vertex_signal_KF_Chi2(signal_idx);
    var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
    var_flightLenSig = sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(signal_idx),Ntp->Vertex_PrimaryVertex_Covariance(signal_idx),
							    Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx)));
    var_sumMuTrkKinkChi2= (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));
    var_segCompMuMin  = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
    var_MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});


    var_MuMu_mindR = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
    var_RelPt_Mu1Tau = Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt();
    var_Eta_au = TauLV.Eta();
    var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(signal_idx,0), Ntp->Vertex_pair_quality(signal_idx,1), Ntp->Vertex_pair_quality(signal_idx,2)});
    var_maxdca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});
    var_MuTau_maxdR = std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});

    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(signal_idx,0),
	  Ntp->Vertex_d0sig_reco(signal_idx,1),
	  Ntp->Vertex_d0sig_reco(signal_idx,2)});

    var_MaxD0Significance = MaxD0Significance;
    var_IsolationMinDist = Ntp->Isolation_MinDist(signal_idx);



    var_tauMass=TauRefitLV.M();
    //    TauMass_all.at(t).Fill(TauRefitLV.M(),1);

    //---------------- define per event resolution categroies 


    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){
	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;

	//	TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	//	TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }
  }
}


void  MuIDStudy::Finish(){


  //  Extradist1d.push_back(&CustMuonIdMC3C);
  //  Extradist1d.push_back(&CustMuonIdData3C);
  //  Significance3C=HConfig.GetTH1D(Name+"_Significance3C","Significance3C",6,0,6,"#mu_{3} MuonId","Events");




  if(mode == RECONSTRUCT){


    SkimConfig SC;
    SC.ApplySkimEfficiency(types, Npassed, Npassed_noweight);

    //    std::cout<<"  CrossSectionandAcceptance.at(2)*Lumi/Npassed.at(i).GetBinContent(0) " << CrossSectionandAcceptance.at(2)*Lumi/Npassed.at(2).GetBinContent(0) <<std::endl;


    Significance3C.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData3C.at(0).GetNbinsX()+1; j++){
      CustMuonIdData3C.at(0).SetBinContent(j, sqrt(CustMuonIdData3C.at(0).GetBinContent(j)));
      CustMuonIdData3C.at(0).SetBinError(j, sqrt(CustMuonIdData3C.at(0).GetBinContent(j)));
    }
    Significance3C.at(2).Divide(&CustMuonIdMC3C.at(2),&CustMuonIdData3C.at(0),1.0,1.0,"B");


    Significance2C.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData2C.at(0).GetNbinsX()+1; j++){
      CustMuonIdData2C.at(0).SetBinContent(j, sqrt(CustMuonIdData2C.at(0).GetBinContent(j)));
      CustMuonIdData2C.at(0).SetBinError(j, sqrt(CustMuonIdData2C.at(0).GetBinContent(j)));
    }
    Significance2C.at(2).Divide(&CustMuonIdMC2C.at(2),&CustMuonIdData2C.at(0),1.0,1.0,"B");


    Significance1C.at(1).Reset();
    for(unsigned int j=0; j< CustMuonIdData1C.at(0).GetNbinsX()+1; j++){
      CustMuonIdData1C.at(0).SetBinContent(j, sqrt(CustMuonIdData1C.at(0).GetBinContent(j)));
      CustMuonIdData1C.at(0).SetBinError(j, sqrt(CustMuonIdData1C.at(0).GetBinContent(j)));
    }
    Significance1C.at(2).Divide(&CustMuonIdMC1C.at(2),&CustMuonIdData1C.at(0),1.0,1.0,"B");






    Significance3B.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData3B.at(0).GetNbinsX()+1; j++){
      CustMuonIdData3B.at(0).SetBinContent(j, sqrt(CustMuonIdData3B.at(0).GetBinContent(j)));
      CustMuonIdData3B.at(0).SetBinError(j, sqrt(CustMuonIdData3B.at(0).GetBinContent(j)));
    }
    Significance3B.at(2).Divide(&CustMuonIdMC3B.at(2),&CustMuonIdData3B.at(0),1.0,1.0,"B");


    Significance2B.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData2B.at(0).GetNbinsX()+1; j++){
      CustMuonIdData2B.at(0).SetBinContent(j, sqrt(CustMuonIdData2B.at(0).GetBinContent(j)));
      CustMuonIdData2B.at(0).SetBinError(j, sqrt(CustMuonIdData2B.at(0).GetBinContent(j)));
    }
    Significance2B.at(2).Divide(&CustMuonIdMC2B.at(2),&CustMuonIdData2B.at(0),1.0,1.0,"B");



    Significance1B.at(1).Reset();
    for(unsigned int j=0; j< CustMuonIdData1B.at(0).GetNbinsX()+1; j++){
      CustMuonIdData1B.at(0).SetBinContent(j, sqrt(CustMuonIdData1B.at(0).GetBinContent(j)));
      CustMuonIdData1B.at(0).SetBinError(j, sqrt(CustMuonIdData1B.at(0).GetBinContent(j)));
    }
    Significance1B.at(2).Divide(&CustMuonIdMC1B.at(2),&CustMuonIdData1B.at(0),1.0,1.0,"B");



    Significance3A.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData3A.at(0).GetNbinsX()+1; j++){
      CustMuonIdData3A.at(0).SetBinContent(j, sqrt(CustMuonIdData3A.at(0).GetBinContent(j)));
      CustMuonIdData3A.at(0).SetBinError(j, sqrt(CustMuonIdData3A.at(0).GetBinContent(j)));
    }
    Significance3A.at(2).Divide(&CustMuonIdMC3A.at(2),&CustMuonIdData3A.at(0),1.0,1.0,"B");


    Significance2A.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData2A.at(0).GetNbinsX()+1; j++){
      CustMuonIdData2A.at(0).SetBinContent(j, sqrt(CustMuonIdData2A.at(0).GetBinContent(j)));
      CustMuonIdData2A.at(0).SetBinError(j, sqrt(CustMuonIdData2A.at(0).GetBinContent(j)));
    }
    Significance2A.at(2).Divide(&CustMuonIdMC2A.at(2),&CustMuonIdData2A.at(0),1.0,1.0,"B");



    Significance1A.at(1).Reset();
    for(unsigned int j=0; j< CustMuonIdData1A.at(0).GetNbinsX()+1; j++){
      CustMuonIdData1A.at(0).SetBinContent(j, sqrt(CustMuonIdData1A.at(0).GetBinContent(j)));
      CustMuonIdData1A.at(0).SetBinError(j, sqrt(CustMuonIdData1A.at(0).GetBinContent(j)));
    }
    Significance1A.at(2).Divide(&CustMuonIdMC1A.at(2),&CustMuonIdData1A.at(0),1.0,1.0,"B");



    Significance3.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData3.at(0).GetNbinsX()+1; j++){
      CustMuonIdData3.at(0).SetBinContent(j, sqrt(CustMuonIdData3.at(0).GetBinContent(j)));
      CustMuonIdData3.at(0).SetBinError(j, sqrt(CustMuonIdData3.at(0).GetBinContent(j)));
    }
    Significance3.at(2).Divide(&CustMuonIdMC3.at(2),&CustMuonIdData3.at(0),1.0,1.0,"B");


    Significance2.at(2).Reset();
    for(unsigned int j=0; j< CustMuonIdData2.at(0).GetNbinsX()+1; j++){
      CustMuonIdData2.at(0).SetBinContent(j, sqrt(CustMuonIdData2.at(0).GetBinContent(j)));
      CustMuonIdData2.at(0).SetBinError(j, sqrt(CustMuonIdData2.at(0).GetBinContent(j)));
    }
    Significance2.at(2).Divide(&CustMuonIdMC2.at(2),&CustMuonIdData2.at(0),1.0,1.0,"B");



    Significance1.at(1).Reset();
    for(unsigned int j=0; j< CustMuonIdData1.at(0).GetNbinsX()+1; j++){
      CustMuonIdData1.at(0).SetBinContent(j, sqrt(CustMuonIdData1.at(0).GetBinContent(j)));
      CustMuonIdData1.at(0).SetBinError(j, sqrt(CustMuonIdData1.at(0).GetBinContent(j)));
    }
    Significance1.at(2).Divide(&CustMuonIdMC1.at(2),&CustMuonIdData1.at(0),1.0,1.0,"B");







    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    //    int id(Ntp->GetMCID());
    //    double scale(1.);
    //    double scaleDsTau(0.637);
    //    double scaleBpTau(0.262);
    //    double scaleB0Tau(0.099);

    //total xsection of producing taus is 12.848 ub 


    //    if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    // ScaleAllHistOfType(2,scale*scaleDsTau);
    
    //if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    // ScaleAllHistOfType(3,scale*scaleB0Tau);

    //if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    //ScaleAllHistOfType(4,scale*scaleBpTau);

    //    }
  }
  Selection::Finish();
}





