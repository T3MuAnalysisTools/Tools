#include "Validation.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

double Validation::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}

Validation::Validation(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

Validation::~Validation(){
  for(unsigned int j=0; j<Npassed.size(); j++){
  Logger(Logger::Info) << "Selection Summary before: "
    << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
    << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  Validation::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
  cut.push_back(0);
  value.push_back(0);
  pass.push_back(false);
  if(i==TriggerOk) cut.at(TriggerOk)=1;
  if(i==PrimeVtx)  cut.at(PrimeVtx)=1;
  if(i==VertexChi2) cut.at(VertexChi2)=5.;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
  title.push_back("");
  distindx.push_back(false);
  dist.push_back(std::vector<float>());
  TString c="_Cut_";c+=i;

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
  
  else if(i==VertexChi2){
    title.at(i)="VertexChi2 ";
    hlabel="Vertex #chi^{2} ";
    Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VertexChi2_",htitle,50,0,50,hlabel,"Events"));
    Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VertexChi2__",htitle,50,0,50,hlabel,"Events"));
  }

  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  MuonsPt=HConfig.GetTH1D(Name+"_MuonsPt","All Muons Pt",25,0,50,"Transverse Momenta of all stored Muons","Events");
  MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of 2 muons pt",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
  MuonsEta=HConfig.GetTH1D(Name+"_MuonsEta","All Muons #eta",25,-2.5,2.5,"Rapidity of all stored Muons","Events");
  MuonsPhi=HConfig.GetTH1D(Name+"_MuonsPhi","All Muons #phi",25,-3.15,3.15,"Azimuthal angle of all stored Muons","Events");
  DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR","Events");
  Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muons status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"","Events");
  Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"","Events");
  Muon1_isTimeValid=HConfig.GetTH1D(Name+"_Muon1_isTimeValid","#mu_{1} isTimevalid",2,-0.5,1.5,"","Events");
  Muon2_isTimeValid=HConfig.GetTH1D(Name+"_Muon2_isTimeValid","#mu_{2} isTimeValid",2,-0.5,1.5,"","Events");
  Muon1_emEt03=HConfig.GetTH1D(Name+"_Muon1_emEt03","",10,0,10,"","Events");
  Muon2_emEt03=HConfig.GetTH1D(Name+"_Muon2_emEt03","",10,0,10,"","Events");
  Muon1_emVetoEt03=HConfig.GetTH1D(Name+"_Muon1_emVetoEt03","",10,0,10,"","Events");
  Muon2_emVetoEt03=HConfig.GetTH1D(Name+"_Muon2_emVetoEt03","",10,0,10,"","Events");
  Muon1_hadEt03=HConfig.GetTH1D(Name+"_Muon1_hadEt03","",10,0,10,"","Events");
  Muon2_hadEt03=HConfig.GetTH1D(Name+"_Muon2_hadEt03","",10,0,10,"","Events");
  Muon1_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt03","",10,0,10,"","Events");
  Muon2_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt03","",10,0,10,"","Events");
  Muon1_nJets03=HConfig.GetTH1D(Name+"_Muon1_nJets03","",10,0,10,"","Events");
  Muon2_nJets03=HConfig.GetTH1D(Name+"_Muon2_nJets03","",10,0,10,"","Events");
  Muon1_nTracks03=HConfig.GetTH1D(Name+"_Muon1_nTracks03","",10,0,10,"","Events");
  Muon2_nTracks03=HConfig.GetTH1D(Name+"_Muon2_nTracks03","",10,0,10,"","Events");
  Muon1_sumPt03=HConfig.GetTH1D(Name+"_Muon1_sumPt03","",10,0,10,"","Events");
  Muon2_sumPt03=HConfig.GetTH1D(Name+"_Muon2_sumPt03","",10,0,10,"","Events");
  Muon1_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt03","",10,0,10,"","Events");
  Muon2_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt03","",10,0,10,"","Events");
  Muon1_emEt05=HConfig.GetTH1D(Name+"_Muon1_emEt05","",10,0,10,"","Events");
  Muon2_emEt05=HConfig.GetTH1D(Name+"_Muon2_emEt05","",10,0,10,"","Events");
  Muon1_emVetoEt05=HConfig.GetTH1D(Name+"_Muon1_emVetoEt05","",10,0,10,"","Events");
  Muon2_emVetoEt05=HConfig.GetTH1D(Name+"_Muon2_emVetoEt05","",10,0,10,"","Events");
  Muon1_hadEt05=HConfig.GetTH1D(Name+"_Muon1_hadEt05","",10,0,10,"","Events");
  Muon2_hadEt05=HConfig.GetTH1D(Name+"_Muon2_hadEt05","",10,0,10,"","Events");
  Muon1_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt05","",10,0,10,"","Events");
  Muon2_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt05","",10,0,10,"","Events");
  Muon1_nJets05=HConfig.GetTH1D(Name+"_Muon1_nJets05","",10,0,10,"","Events");
  Muon2_nJets05=HConfig.GetTH1D(Name+"_Muon2_nJets05","",10,0,10,"","Events");
  Muon1_nTracks05=HConfig.GetTH1D(Name+"_Muon1_nTracks05","",10,0,10,"","Events");
  Muon2_nTracks05=HConfig.GetTH1D(Name+"_Muon2_nTracks05","",10,0,10,"","Events");
  Muon1_sumPt05=HConfig.GetTH1D(Name+"_Muon1_sumPt05","",10,0,10,"","Events");
  Muon2_sumPt05=HConfig.GetTH1D(Name+"_Muon2_sumPt05","",10,0,10,"","Events");
  Muon1_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt05","",10,0,10,"","Events");
  Muon2_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt05","",10,0,10,"","Events");
  Muon1_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt03","",10,0,10,"","Events");
  Muon2_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt03","",10,0,10,"","Events");
  Muon1_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt03","",10,0,10,"","Events");
  Muon2_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt03","",10,0,10,"","Events");
  Muon1_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt03","",10,0,10,"","Events");
  Muon2_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt03","",10,0,10,"","Events");
  Muon1_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold03","",10,0,10,"","Events");
  Muon2_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold03","",10,0,10,"","Events");
  Muon1_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt03","",10,0,10,"","Events");
  Muon2_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt03","",10,0,10,"","Events");
  Muon1_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold03","",10,0,10,"","Events");
  Muon2_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold03","",10,0,10,"","Events");
  Muon1_sumPUPt03=HConfig.GetTH1D(Name+"_Muon1_sumPUPt03","",10,0,10,"","Events");
  Muon2_sumPUPt03=HConfig.GetTH1D(Name+"_Muon2_sumPUPt03","",10,0,10,"","Events");
  Muon1_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt04","",10,0,10,"","Events");
  Muon2_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt04","",10,0,10,"","Events");
  Muon1_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt04","",10,0,10,"","Events");
  Muon2_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt04","",10,0,10,"","Events");
  Muon1_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt04","",10,0,10,"","Events");
  Muon2_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt04","",10,0,10,"","Events");
  Muon1_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold04","",10,0,10,"","Events");
  Muon2_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold04","",10,0,10,"","Events");
  Muon1_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt04","",10,0,10,"","Events");
  Muon2_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt04","",10,0,10,"","Events");
  Muon1_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold04","",10,0,10,"","Events");
  Muon2_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold04","",10,0,10,"","Events");
  Muon1_sumPUPt04=HConfig.GetTH1D(Name+"_Muon1_sumPUPt04","",10,0,10,"","Events");
  Muon2_sumPUPt04=HConfig.GetTH1D(Name+"_Muon2_sumPUPt04","",10,0,10,"","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}



void  Validation::Store_ExtraDist(){
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&MuonsPt);
  Extradist1d.push_back(&MuonsEta);
  Extradist1d.push_back(&MuonsPhi);
  Extradist1d.push_back(&MuonsPtRatio);
  Extradist1d.push_back(&DimuondR);
  Extradist1d.push_back(&Muon1TrkdR);
  Extradist1d.push_back(&Muon2TrkdR);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&TripleMass);
  Extradist1d.push_back(&Muon1_isGlobal);
  Extradist2d.push_back(&PhiMassVsDsMass);
  Extradist1d.push_back(&Muon2_isGlobal);
  Extradist1d.push_back(&Muon1_isStandAlone);
  Extradist1d.push_back(&Muon2_isStandAlone);
  Extradist1d.push_back(&Muon1_isTracker);
  Extradist1d.push_back(&Muon2_isTracker);
  Extradist1d.push_back(&Muon1_isCalo);
  Extradist1d.push_back(&Muon2_isCalo);
  Extradist1d.push_back(&Muon1_isIsolationValid);
  Extradist1d.push_back(&Muon2_isIsolationValid);
  Extradist1d.push_back(&Muon1_isTimeValid);
  Extradist1d.push_back(&Muon2_isTimeValid);
  Extradist1d.push_back(&Muon1_emEt03);
  Extradist1d.push_back(&Muon2_emEt03);
  Extradist1d.push_back(&Muon1_emVetoEt03);
  Extradist1d.push_back(&Muon2_emVetoEt03);
  Extradist1d.push_back(&Muon1_hadEt03);
  Extradist1d.push_back(&Muon2_hadEt03);
  Extradist1d.push_back(&Muon1_hadVetoEt03);
  Extradist1d.push_back(&Muon2_hadVetoEt03);
  Extradist1d.push_back(&Muon1_nJets03);
  Extradist1d.push_back(&Muon2_nJets03);
  Extradist1d.push_back(&Muon1_nTracks03);
  Extradist1d.push_back(&Muon2_nTracks03);
  Extradist1d.push_back(&Muon1_sumPt03);
  Extradist1d.push_back(&Muon2_sumPt03);
  Extradist1d.push_back(&Muon1_trackerVetoPt03);
  Extradist1d.push_back(&Muon2_trackerVetoPt03);
  Extradist1d.push_back(&Muon1_emEt05);
  Extradist1d.push_back(&Muon2_emEt05);
  Extradist1d.push_back(&Muon1_emVetoEt05);
  Extradist1d.push_back(&Muon2_emVetoEt05);
  Extradist1d.push_back(&Muon1_hadEt05);
  Extradist1d.push_back(&Muon2_hadEt05);
  Extradist1d.push_back(&Muon1_hadVetoEt05);
  Extradist1d.push_back(&Muon2_hadVetoEt05);
  Extradist1d.push_back(&Muon1_nJets05);
  Extradist1d.push_back(&Muon2_nJets05);
  Extradist1d.push_back(&Muon1_nTracks05);
  Extradist1d.push_back(&Muon2_nTracks05);
  Extradist1d.push_back(&Muon1_sumPt05);
  Extradist1d.push_back(&Muon2_sumPt05);
  Extradist1d.push_back(&Muon1_trackerVetoPt05);
  Extradist1d.push_back(&Muon2_trackerVetoPt05);
  Extradist1d.push_back(&Muon1_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon2_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon1_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon2_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon1_sumPhotonEt03);
  Extradist1d.push_back(&Muon2_sumPhotonEt03);
  Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon1_sumPUPt03);
  Extradist1d.push_back(&Muon2_sumPUPt03);
  Extradist1d.push_back(&Muon1_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon2_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon1_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon2_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon1_sumPhotonEt04);
  Extradist1d.push_back(&Muon2_sumPhotonEt04);
  Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon1_sumPUPt04);
  Extradist1d.push_back(&Muon2_sumPUPt04);

}

void  Validation::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  value.at(TriggerOk)=(Ntp->EventNumber()%1000)==1;
  pass.at(TriggerOk)=true;


  double wobs=1;
  double w;
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */}
  else{w=1;}



  bool status=AnalysisCuts(t,w,wobs);
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){

  NVtx.at(t).Fill(Ntp->NVtx(),w);

  for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++){
    MuonsPt.at(t).Fill(Ntp->Muon_P4(iMuon).Pt(),w);
    MuonsEta.at(t).Fill(Ntp->Muon_P4(iMuon).Eta(),w);
    MuonsPhi.at(t).Fill(Ntp->Muon_P4(iMuon).Phi(),w);

  }
  double deltaMass(999.);
  unsigned int pair_index(0);
  if(Ntp->NTwoMuonsTrack()!=0){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
    unsigned int muon_1 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
    unsigned int muon_2 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);
	 Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(muon_1),w);
    Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(muon_2),w);
    Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(muon_1),w);
    Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(muon_2),w);
    Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(muon_1),w);
    Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(muon_2),w);
    Muon1_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(muon_1),w);
    Muon2_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(muon_2),w);
    Muon1_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(muon_1),w);
    Muon2_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(muon_2),w);
    Muon1_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(muon_1),w);
    Muon2_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(muon_2),w);
    Muon1_emEt03.at(t).Fill(Ntp->Muon_emEt03(muon_1),w);
    Muon2_emEt03.at(t).Fill(Ntp->Muon_emEt03(muon_2),w);
    Muon1_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(muon_1),w);
    Muon2_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(muon_2),w);
    Muon1_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(muon_1),w);
    Muon2_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(muon_2),w);
    Muon1_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(muon_1),w);
    Muon2_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(muon_2),w);
    Muon1_nJets03.at(t).Fill(Ntp->Muon_nJets03(muon_1),w);
    Muon2_nJets03.at(t).Fill(Ntp->Muon_nJets03(muon_2),w);
    Muon1_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(muon_1),w);
    Muon2_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(muon_2),w);
    Muon1_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(muon_1),w);
    Muon2_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(muon_2),w);
    Muon1_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(muon_1),w);
    Muon2_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(muon_2),w);
    Muon1_emEt05.at(t).Fill(Ntp->Muon_emEt05(muon_1),w);
    Muon2_emEt05.at(t).Fill(Ntp->Muon_emEt05(muon_2),w);
    Muon1_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(muon_1),w);
    Muon2_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(muon_2),w);
    Muon1_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(muon_1),w);
    Muon2_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(muon_2),w);
    Muon1_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(muon_1),w);
    Muon2_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(muon_2),w);
    Muon1_nJets05.at(t).Fill(Ntp->Muon_nJets05(muon_1),w);
    Muon2_nJets05.at(t).Fill(Ntp->Muon_nJets05(muon_2),w);
    Muon1_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(muon_1),w);
    Muon2_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(muon_2),w);
    Muon1_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(muon_1),w);
    Muon2_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(muon_2),w);
    Muon1_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(muon_1),w);
    Muon2_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(muon_2),w);
    Muon1_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(muon_1),w);
    Muon2_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(muon_2),w);
    Muon1_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(muon_1),w);
    Muon2_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(muon_2),w);
    Muon1_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(muon_1),w);
    Muon2_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(muon_2),w);
    Muon1_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(muon_1),w);
    Muon2_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(muon_2),w);
    Muon1_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(muon_1),w);
    Muon2_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(muon_2),w);
    Muon1_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(muon_1),w);
    Muon2_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(muon_2),w);
    Muon1_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(muon_1),w);
    Muon2_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(muon_2),w);
    Muon1_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(muon_1),w);
    Muon2_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(muon_2),w);
    Muon1_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(muon_1),w);
    Muon2_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(muon_2),w);
    Muon1_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(muon_1),w);
    Muon2_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(muon_2),w);
    Muon1_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(muon_1),w);
    Muon2_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(muon_2),w);
    Muon1_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(muon_1),w);
    Muon2_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(muon_2),w);
    Muon1_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(muon_1),w);
    Muon2_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(muon_2),w);
    Muon1_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(muon_1),w);
    Muon2_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(muon_2),w);
    
	 if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
	 	deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
      pair_index = i2M;
    }
    }

	 DimuondR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Phi(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Phi()));
	 Muon1TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Phi()));
	 Muon2TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Phi()));
    MuonsPtRatio.at(t).Fill(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Pt()/Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1)).Pt(),w );
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+ 
      Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M(), w);
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+
      Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M();
    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);
  }
  }
}



void  Validation::Finish(){
  Selection::Finish();
}





