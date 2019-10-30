#include "DsToPhiPi.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"



DsToPhiPi::DsToPhiPi(TString Name_, TString id_):
  Selection(Name_,id_),
  RelLumiB(0.2087),
  RelLumiC(0.1038),
  RelLumiD(0.2358),
  RelLumiE(0.1075),
  RelLumiF(0.0996)
{
  // This is a class constructor;
}

DsToPhiPi::~DsToPhiPi(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  DsToPhiPi::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)           cut.at(TriggerOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==GlobalMu)        cut.at(GlobalMu)=1;
    if(i==Chi2Cut)         cut.at(Chi2Cut)=1;
    //    if(i==MuCharge)        cut.at(MuCharge)=1;
    if(i==Mass2Mu)         cut.at(Mass2Mu)=1;
    if(i==Mu1dR)           cut.at(Mu1dR)=1;
    if(i==Mu2dR)           cut.at(Mu2dR)=1;
    if(i==TrkdR)           cut.at(TrkdR)=1;
    if(i==Mu1pt)           cut.at(Mu1pt)=1;
    if(i==Mu2pt)           cut.at(Mu2pt)=1;
    if(i==Trkpt)           cut.at(Trkpt)=1;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==TriggerOk){
      title.at(i)="Trigger";
      hlabel="Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==GlobalMu){
      title.at(i)="Muons are Global";
      hlabel="Muons are Global";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
    else if(i==Chi2Cut){
      title.at(i)="Triple Vertex Chi Squared $<$ 15";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Chi squared of triple vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Chi2Cut_",htitle,50,0,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Chi2Cut_",htitle,50,0,25,hlabel,"Events"));
    }
    // else if(i==MuCharge){
    //   title.at(i)="Muons opposite charge";
    //   hlabel="Muons have opposite charge";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==Mass2Mu){
      title.at(i)="1.00 $<$ $M_{2\\mu}$ $<$ 1.04 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Invariant mass of 2 muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1dR){
      title.at(i)="Mu01 dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Muon 1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1dR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1dR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==Mu2dR){
      title.at(i)="Mu02 dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Muon 2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2dR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2dR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==TrkdR){
      title.at(i)="Tr dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Track";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrkdR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrkdR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==Mu1pt){
      title.at(i)="$\\mu_{1}$ Pt $>$ 3 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Muon 1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1pt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1pt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==Mu2pt){
      title.at(i)="$\\mu_{2}$ Pt $>$ 3 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Muon 2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2pt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2pt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==Trkpt){
      title.at(i)="Track Pt $>$ 2 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Track";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trkpt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trkpt_",htitle,80,0,20,hlabel,"Events"));
    }
  }

  // Track Candidate Information
  Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",36,-0.5,35.5,"p (track)","Events");
  Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",36,-0.5,35.5,"E (track)","Events");
  Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",30,-3.5,3.5,"#phi","Events");
  Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",20,-0.5,4.5,"#chi^{2} (track fit)","Events");
  Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",36,-0.5,35.5,"n valid track hits","Events");
  Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",3,-1.5,1.5,"Track charge","Events");
  
  Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); 
  Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
  Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  
  
  Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
  Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{2} #phi","Events"); 
  Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
  Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
  
  Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","Trigger Matching mu1",50,0,0.05,"#Delta R Trigger Match #mu_{1}","Events");
  Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","Trigger Matching mu2",50,0,0.05,"#Delta R Trigger Match #mu_{2}","Events");
  Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","Trigger Matching track",50,0,0.05,"#Delta R Trigger Match track","Events");
  
  Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  
  
  Track_PtB=HConfig.GetTH1D(Name+"_Track_PtB","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_EtaB=HConfig.GetTH1D(Name+"_Track_EtaB","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Muon1_PtB=HConfig.GetTH1D(Name+"_Muon1_PtB","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_EtaB=HConfig.GetTH1D(Name+"_Muon1_EtaB","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon2_PtB=HConfig.GetTH1D(Name+"_Muon2_PtB","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_EtaB=HConfig.GetTH1D(Name+"_Muon2_EtaB","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");

  Track_PtC=HConfig.GetTH1D(Name+"_Track_PtC","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_EtaC=HConfig.GetTH1D(Name+"_Track_EtaC","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Muon1_PtC=HConfig.GetTH1D(Name+"_Muon1_PtC","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_EtaC=HConfig.GetTH1D(Name+"_Muon1_EtaC","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon2_PtC=HConfig.GetTH1D(Name+"_Muon2_PtC","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_EtaC=HConfig.GetTH1D(Name+"_Muon2_EtaC","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");

  Track_PtD=HConfig.GetTH1D(Name+"_Track_PtD","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_EtaD=HConfig.GetTH1D(Name+"_Track_EtaD","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Muon1_PtD=HConfig.GetTH1D(Name+"_Muon1_PtD","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_EtaD=HConfig.GetTH1D(Name+"_Muon1_EtaD","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon2_PtD=HConfig.GetTH1D(Name+"_Muon2_PtD","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_EtaD=HConfig.GetTH1D(Name+"_Muon2_EtaD","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");


  Track_PtE=HConfig.GetTH1D(Name+"_Track_PtE","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_EtaE=HConfig.GetTH1D(Name+"_Track_EtaE","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Muon1_PtE=HConfig.GetTH1D(Name+"_Muon1_PtE","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_EtaE=HConfig.GetTH1D(Name+"_Muon1_EtaE","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon2_PtE=HConfig.GetTH1D(Name+"_Muon2_PtE","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_EtaE=HConfig.GetTH1D(Name+"_Muon2_EtaE","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");


  Track_PtF=HConfig.GetTH1D(Name+"_Track_PtF","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_EtaF=HConfig.GetTH1D(Name+"_Track_EtaF","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Muon1_PtF=HConfig.GetTH1D(Name+"_Muon1_PtF","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_EtaF=HConfig.GetTH1D(Name+"_Muon1_EtaF","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon2_PtF=HConfig.GetTH1D(Name+"_Muon2_PtF","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_EtaF=HConfig.GetTH1D(Name+"_Muon2_EtaF","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");


  DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,1,"dR","Events");
  Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,1,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu invariant mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  PhiPlusTrackMass=HConfig.GetTH1D(Name+"_PhiPlusTrackMass","#mu#mu + track invariant mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu invariant Mass vs. #mu#mu + track invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
 
  DsMass=HConfig.GetTH1D(Name+"_DsMass","Ds invariant mass 2017",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 

  DsMassB=HConfig.GetTH1D(Name+"_DsMassB","Ds invariant mass Run B",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 
  DsMassC=HConfig.GetTH1D(Name+"_DsMassC","Ds invariant mass Run C",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 
  DsMassD=HConfig.GetTH1D(Name+"_DsMassD","Ds invariant mass Run D",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 
  DsMassE=HConfig.GetTH1D(Name+"_DsMassE","Ds invariant mass Run E",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 
  DsMassF=HConfig.GetTH1D(Name+"_DsMassF","Ds invariant mass Run F",80,1.64,2.3,"M_{Ds} (GeV)", "Events"); 


  Muon1_PtF_peak=HConfig.GetTH1D(Name+"_Muon1_PtF_peak","Transverse Pt (muon 1) in ds peak",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_PtF_sideband=HConfig.GetTH1D(Name+"_Muon1_PtF_sideband","Transverse Pt (muon 1) in sideband",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_PtF_substracted=HConfig.GetTH1D(Name+"_Muon1_PtF_substracted","Transverse Pt (muon 1) in subs",25,0,30,"#mu_{1} p_{T} (GeV)","Events");


  DsDecayL=HConfig.GetTH1D(Name+"_DsDecayL","Decay Length",25,0,0.1,"Decay length, cm","Events");





  DsGenMatch=HConfig.GetTH1D(Name+"_DsGenMatch","dR between Gen Ds to Track",50,0,.1,"dR","Events");

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsToPhiPi::Store_ExtraDist(){ 
  
  //  Extradist1d.push_back(&Track_P);
  //  Extradist1d.push_back(&Track_E);


  Extradist1d.push_back(&Track_PtB);
  Extradist1d.push_back(&Track_EtaB);
  Extradist1d.push_back(&Muon1_PtB);
  Extradist1d.push_back(&Muon1_EtaB);
  Extradist1d.push_back(&Muon2_PtB);
  Extradist1d.push_back(&Muon2_EtaB);


  Extradist1d.push_back(&Track_PtC);
  Extradist1d.push_back(&Track_EtaC);
  Extradist1d.push_back(&Muon1_PtC);
  Extradist1d.push_back(&Muon1_EtaC);
  Extradist1d.push_back(&Muon2_PtC);
  Extradist1d.push_back(&Muon2_EtaC);


  Extradist1d.push_back(&Track_PtD);
  Extradist1d.push_back(&Track_EtaD);
  Extradist1d.push_back(&Muon1_PtD);
  Extradist1d.push_back(&Muon1_EtaD);
  Extradist1d.push_back(&Muon2_PtD);
  Extradist1d.push_back(&Muon2_EtaD);


  Extradist1d.push_back(&Track_PtE);
  Extradist1d.push_back(&Track_EtaE);
  Extradist1d.push_back(&Muon1_PtE);
  Extradist1d.push_back(&Muon1_EtaE);
  Extradist1d.push_back(&Muon2_PtE);
  Extradist1d.push_back(&Muon2_EtaE);


  Extradist1d.push_back(&Track_PtF);
  Extradist1d.push_back(&Track_EtaF);
  Extradist1d.push_back(&Muon1_PtF);
  Extradist1d.push_back(&Muon1_EtaF);
  Extradist1d.push_back(&Muon2_PtF);
  Extradist1d.push_back(&Muon2_EtaF);

  Extradist1d.push_back(&Muon1_PtF_substracted);
  Extradist1d.push_back(&Muon1_PtF_peak);
  Extradist1d.push_back(&Muon1_PtF_sideband);


  //  Extradist1d.push_back(&Track_Phi);
  //  Extradist1d.push_back(&Track_normalizedChi2);
  //  Extradist1d.push_back(&Track_numberOfValidHits);
  //  Extradist1d.push_back(&Track_charge);
  //  Extradist1d.push_back(&Muon1_E);
  //  Extradist1d.push_back(&Muon1_Phi);
  //  Extradist1d.push_back(&Muon2_E);
  //  Extradist1d.push_back(&Muon2_Phi);
  //  Extradist1d.push_back(&DimuondR);
  //  Extradist1d.push_back(&Muon1TrkdR);
  //  Extradist1d.push_back(&Muon2TrkdR);
  Extradist1d.push_back(&PhiMass);

  Extradist1d.push_back(&DsDecayL);

  //  Extradist1d.push_back(&PhiPlusTrackMass);
  //  Extradist2d.push_back(&PhiMassVsDsMass);
  //  Extradist1d.push_back(&Muon1_isGlobal);
  //  Extradist1d.push_back(&Muon2_isGlobal);
  //  Extradist1d.push_back(&Muon1_isStandAlone);
  //  Extradist1d.push_back(&Muon2_isStandAlone);
  //  Extradist1d.push_back(&Muon1_isTracker);
  //  Extradist1d.push_back(&Muon2_isTracker);
  //  Extradist1d.push_back(&Muon1_isCalo);
  //  Extradist1d.push_back(&Muon2_isCalo);
  Extradist1d.push_back(&Track_TriggerMatchdR);
  Extradist1d.push_back(&Muon1_TriggerMatchdR);
  Extradist1d.push_back(&Muon2_TriggerMatchdR);
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&DsMass);
  Extradist1d.push_back(&DsMassB);
  Extradist1d.push_back(&DsMassC);
  Extradist1d.push_back(&DsMassD);
  Extradist1d.push_back(&DsMassE);
  Extradist1d.push_back(&DsMassF);
  //  Extradist1d.push_back(&DsGenMatch);
	 
}


void  DsToPhiPi::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  value.at(TriggerOk) = 0;

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
  for(unsigned int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }
  }

  if(DoubleMuFired) L1Ok = true;
  value.at(TriggerOk)=(HLTOk and L1Ok);


    
  int mu1=-1, mu2=-1, track=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999.0;
  double check_PhiMass(999.0);
  value.at(is2MuTrk) = 0; 
  value.at(Chi2Cut) = 0;
  value.at(Mass2Mu) = 0;
  //  value.at(MuCharge) = 0;
  value.at(Mu1dR) = 0;
  value.at(Mu2dR) = 0;
  value.at(TrkdR) = 0;
  if(Ntp->NTwoMuonsTrack()!=0) value.at(is2MuTrk) = 1;


  //  std::cout<<"N two muons and tracks:   "<< Ntp->NTwoMuonsTrack()<< " N signal    "<< Ntp-> NThreeMuons() << "  Num vert   "<< Ntp->NumberOfSVertices() << std::endl;
  if (value.at(is2MuTrk)==1 && Ntp-> NThreeMuons()==0){

    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);
      double tmp_PhiMass = (Ntp->Muon_P4(tmp_mu1)+Ntp->Muon_P4(tmp_mu2)).M();

      if (abs(tmp_PhiMass-1.02)<=check_PhiMass || (tmp_PhiMass > .95 && tmp_PhiMass < 1.1)) {
	if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M)){
	  tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
	  check_PhiMass = abs(tmp_PhiMass-1.01);
	  mu1 = tmp_mu1;
	  mu2 = tmp_mu2;
	  track = tmp_track;
	  tmp_idx = i2M;
	}
      } 
    }

    value.at(GlobalMu) = Ntp->Muon_isGlobalMuon(mu1)==1 && Ntp->Muon_isGlobalMuon(mu2)==1;
    value.at(Chi2Cut)  = tmp_chisq;
    value.at(Mass2Mu)  = (Ntp->Muon_P4(mu1) + Ntp->Muon_P4(mu2)).M();
    //    value.at(MuCharge) = Ntp->Muon_charge(mu1)!=Ntp->Muon_charge(mu2);
    value.at(Mu1dR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0);
    value.at(Mu2dR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1);
    value.at(TrkdR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2);
    value.at(Mu1pt) = Ntp->Muon_P4(mu1).Pt();
    value.at(Mu2pt) = Ntp->Muon_P4(mu2).Pt();
    value.at(Trkpt) = Ntp->Track_P4(track).Pt();
    
}
  
  pass.at(is2MuTrk) = (value.at(is2MuTrk) >0 );
  pass.at(TriggerOk)= (value.at(TriggerOk)==cut.at(TriggerOk));
  pass.at(GlobalMu) = value.at(GlobalMu)==cut.at(GlobalMu);
  pass.at(Chi2Cut)  = value.at(Chi2Cut) >= 0 && value.at(Chi2Cut) < 15.01;
  pass.at(Mass2Mu)  = value.at(Mass2Mu) > 0.99 && value.at(Mass2Mu) < 1.041;
  //  pass.at(MuCharge) = value.at(MuCharge)==cut.at(MuCharge);
  pass.at(Mu1dR)    = value.at(Mu1dR) < .03;
  pass.at(Mu2dR)    = value.at(Mu2dR) < .03;
  pass.at(TrkdR)    = value.at(TrkdR) < .03;
  pass.at(Mu1pt)    = value.at(Mu1pt) > 3;
  pass.at(Mu2pt)    = value.at(Mu2pt) > 3;
  pass.at(Trkpt)    = value.at(Trkpt) > 1.5;
  
  


  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){



    NVtx.at(t).Fill(Ntp->NVtx(),w);

    // Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    // Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);

    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    //    Track_charge.at(t).Fill(Ntp->Track_charge(track),w);

    // Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    // Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    // Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    // Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    // Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    // Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
 

    Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0),w);
    Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1),w);
    Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2),w);
    //    std::cout<<"   "<< Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(tmp_idx),Ntp->Vertex_MatchedPrimaryVertex(tmp_idx)) << std::endl;
    TLorentzVector DsP4 = Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1)) + Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0));
    DsDecayL.at(t).Fill(Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(tmp_idx),Ntp->Vertex_MatchedPrimaryVertex(tmp_idx)).Mag()  *  DsP4.M()/DsP4.P(),1);

    bool isPromt(true);
    //    std::cout<<"reco signal particle P  "<< DsP4.P() << std::endl;
    if(id != 1){
      //std::cout<<" Loop Over Signal Particles  "<< std::endl;
            for(unsigned int isigp = 0; isigp < Ntp-> NMCSignalParticles() ; isigp++){
	      //std::cout<<" Id  "<< Ntp->MCSignalParticle_pdgid(isigp)<< "   MC P4  "<< Ntp->MCSignalParticle_p4(isigp).P() <<  "  N sources   "<< Ntp->NMCSignalParticleSources(isigp) << std::endl;

	//	std::cout<<"loop over sources  "<< std::endl;
	for(unsigned int is =0; is< Ntp->NMCSignalParticleSources(isigp); is++){

	  //	  std::cout<<" source ID   "<< Ntp->MCSignalParticle_Sourcepdgid(isigp,is)<< std::endl;
	  if( abs(  Ntp->MCSignalParticle_Sourcepdgid(isigp,is)  )  > 400 ){ // then it means that Ds is coming from B decay
	    isPromt = false;
	  }

	}
      }
    }

    //    if()



    // Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
    // Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
    // Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
    // Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);


    // Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
    // Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
    // Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
    // Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);


   
    //  DimuondR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Muon_P4(mu2)),w);
    //   Muon1TrkdR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Track_P4(track)),w);
    //  Muon2TrkdR.at(t).Fill(Ntp->Muon_P4(mu2).DeltaR(Ntp->Track_P4(track)),w);
    
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M(), w);
    PhiPlusTrackMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M(), w);
    
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M();
    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);



    DsMass.at(t).Fill(dsmass,w);



    if(id==1){

      if(dsmass > 1.93 && dsmass < 2.01){
	Muon1_PtF_peak.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      }
      if(dsmass > 1.70 && dsmass < 1.80){
	Muon1_PtF_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      }
      

      if(Ntp->WhichEra(2017).Contains("RunB") ){
	DsMassB.at(t).Fill(dsmass,w);
	Track_PtB.at(t).Fill(Ntp->Track_P4(track).Pt(),w*RelLumiB);
	Track_EtaB.at(t).Fill(Ntp->Track_P4(track).Eta(),w*RelLumiB);
	Muon1_PtB.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*RelLumiB);
	Muon1_EtaB.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*RelLumiB);
	Muon2_PtB.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*RelLumiB);
	Muon2_EtaB.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*RelLumiB);
      }

      if(Ntp->WhichEra(2017).Contains("RunC") ) {
	DsMassC.at(t).Fill(dsmass,w);
	Track_PtC.at(t).Fill(Ntp->Track_P4(track).Pt(),w*RelLumiC);
	Track_EtaC.at(t).Fill(Ntp->Track_P4(track).Eta(),w*RelLumiC);
	Muon1_PtC.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*RelLumiC);
	Muon1_EtaC.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*RelLumiC);
	Muon2_PtC.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*RelLumiC);
	Muon2_EtaC.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*RelLumiC);
      }

      if(Ntp->WhichEra(2017).Contains("RunD") ){
	DsMassD.at(t).Fill(dsmass,w);
	Track_PtD.at(t).Fill(Ntp->Track_P4(track).Pt(),w*RelLumiD);
	Track_EtaD.at(t).Fill(Ntp->Track_P4(track).Eta(),w*RelLumiD);
	Muon1_PtD.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*RelLumiD);
	Muon1_EtaD.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*RelLumiD);
	Muon2_PtD.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*RelLumiD);
	Muon2_EtaD.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*RelLumiD);
      }

      if(Ntp->WhichEra(2017).Contains("RunE") ){
	DsMassE.at(t).Fill(dsmass,w);
	Track_PtE.at(t).Fill(Ntp->Track_P4(track).Pt(),w*RelLumiE);
	Track_EtaE.at(t).Fill(Ntp->Track_P4(track).Eta(),w*RelLumiE);
	Muon1_PtE.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*RelLumiE);
	Muon1_EtaE.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*RelLumiE);
	Muon2_PtE.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*RelLumiE);
	Muon2_EtaE.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*RelLumiE);
      }

      if(Ntp->WhichEra(2017).Contains("RunF") ){
	DsMassF.at(t).Fill(dsmass,w);    
	//	std::cout<<"DsMassF  "<< dsmass << std::endl;
	Track_PtF.at(t).Fill(Ntp->Track_P4(track).Pt(),w*RelLumiF);
	Track_EtaF.at(t).Fill(Ntp->Track_P4(track).Eta(),w*RelLumiF);
	Muon1_PtF.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*RelLumiF);
	Muon1_EtaF.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*RelLumiF);
	Muon2_PtF.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*RelLumiF);
	Muon2_EtaF.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*RelLumiF);


      }
    }

    if(id==30){

      DsGenMatch.at(t).Fill(Ntp->DsGenMatch(tmp_idx),w);

    }

  }
}

void  DsToPhiPi::Finish(){

  double scaleRunF(15.);
  Muon1_PtF_sideband.at(0).Scale(scaleRunF/Muon1_PtF_sideband.at(0).Integral());

  Muon1_PtF_substracted.at(0).Add(&Muon1_PtF_peak.at(0));
  Muon1_PtF_substracted.at(0).Add(&Muon1_PtF_sideband.at(0),-1);

  //  if(mode == RECONSTRUCT){
    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    //      double scale(1.);
    //      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/1;
    //      ScaleAllHistOfType(i,scale);
    //    }
    //}
  Selection::Finish();

}



