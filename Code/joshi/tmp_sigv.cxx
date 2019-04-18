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
  if(i==L1TOk) cut.at(L1TOk)=1;
  if(i==HLTOk) cut.at(HLTOk)=1;
  if(i==PrimeVtx)  cut.at(PrimeVtx)=1;
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

    else if(i==L1TOk){
      title.at(i)="L1 Trigger";
      hlabel="L1 Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     
	  else if(i==HLTOk){
      title.at(i)="HLT";
      hlabel="HLT ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }   

  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms

  //Event Information
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  NTracks=HConfig.GetTH1D(Name+"_NTracks","NTracks",31,-0.5,30.5,"Number of Tracks","Events");
  HLT_Tau3Mu=HConfig.GetTH1D(Name+"_HLT_Tau3Mu","_HLT_Tau3Mu",3,-0.5,2.5);
 
  //Dimuon Information (Muons from dimuon + track candidates)
  Muon1Muon2dR=HConfig.GetTH1D(Name+"_Muon1Muon2","dR between the muon pair",20,0,1,"dR","Events");
  Muon2Muon3dR=HConfig.GetTH1D(Name+"_Muon2Muon3","dR between the highest p muon and the track",100,0,5,"dR (#mu_{1}#mu_{2})","Events");
  Muon2Muon3dR=HConfig.GetTH1D(Name+"_Muon2Muon3dR","dR between the lowest p muon and the track",100,0,5,"dR (#mu_{2}#mu_{3})","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
  
  // Three muon variables
  Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon3_isGlobal=HConfig.GetTH1D(Name+"_Muon3_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon3_isStandAlone=HConfig.GetTH1D(Name+"_Muon3_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon3_isTracker=HConfig.GetTH1D(Name+"_Muon3_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  Muon3_isCalo=HConfig.GetTH1D(Name+"_Muon3_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"#mu_{1} is Iso valid","Events");
  Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"#mu_{2} is Iso valid","Events");
  Muon3_isIsolationValid=HConfig.GetTH1D(Name+"_Muon3_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"#mu_{3} is Iso valid","Events");
  Muon1_isTimeValid=HConfig.GetTH1D(Name+"_Muon1_isTimeValid","#mu_{1} isTimevalid",2,-0.5,1.5,"#mu_{1} is Time valid","Events");
  Muon2_isTimeValid=HConfig.GetTH1D(Name+"_Muon2_isTimeValid","#mu_{2} isTimeValid",2,-0.5,1.5,"#mu_{2} is Time valid","Events");
  Muon3_isTimeValid=HConfig.GetTH1D(Name+"_Muon3_isTimeValid","#mu_{2} isTimeValid",2,-0.5,1.5,"#mu_{3} is Time valid","Events");
  Muon1_emEt03=HConfig.GetTH1D(Name+"_Muon1_emEt03","",10,0,10,"#mu_{1} EM E_{T}03","Events");
  Muon2_emEt03=HConfig.GetTH1D(Name+"_Muon2_emEt03","",10,0,10,"#mu_{2} EM E_{T}03","Events");
  Muon3_emEt03=HConfig.GetTH1D(Name+"_Muon3_emEt03","",10,0,10,"#mu_{3} EM E_{T}03","Events");
  Muon1_emVetoEt03=HConfig.GetTH1D(Name+"_Muon1_emVetoEt03","",10,0,10,"#mu_{1} EM veto E_{T}03","Events");
  Muon2_emVetoEt03=HConfig.GetTH1D(Name+"_Muon2_emVetoEt03","",10,0,10,"#mu_{2} EM veto E_{T}03","Events");
  Muon3_emVetoEt03=HConfig.GetTH1D(Name+"_Muon3_emVetoEt03","",10,0,10,"#mu_{3} EM veto E_{T}03","Events");
  Muon1_hadEt03=HConfig.GetTH1D(Name+"_Muon1_hadEt03","",10,0,10,"#mu_{1} hadron E_{T}03","Events");
  Muon2_hadEt03=HConfig.GetTH1D(Name+"_Muon2_hadEt03","",10,0,10,"#mu_{2} hadron E_{T}03","Events");
  Muon3_hadEt03=HConfig.GetTH1D(Name+"_Muon3_hadEt03","",10,0,10,"#mu_{3} hadron E_{T}03","Events");
  Muon1_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt03","",10,0,10,"#mu_{1} hadron veto E_{T}03","Events");
  Muon2_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt03","",10,0,10,"#mu_{2} hadron veto E_{T}03","Events");
  Muon3_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon3_hadVetoEt03","",10,0,10,"#mu_{3} hadron veto E_{T}03","Events");
  Muon1_nJets03=HConfig.GetTH1D(Name+"_Muon1_nJets03","",10,0,10,"#mu_{1} njets03","Events");
  Muon2_nJets03=HConfig.GetTH1D(Name+"_Muon2_nJets03","",10,0,10,"#mu_{2} njets03","Events");
  Muon3_nJets03=HConfig.GetTH1D(Name+"_Muon3_nJets03","",10,0,10,"#mu_{3} njets03","Events");
  Muon1_nTracks03=HConfig.GetTH1D(Name+"_Muon1_nTracks03","",10,0,10,"#mu_{1} ntracks03","Events");
  Muon2_nTracks03=HConfig.GetTH1D(Name+"_Muon2_nTracks03","",10,0,10,"#mu_{2} ntracks03","Events");
  Muon3_nTracks03=HConfig.GetTH1D(Name+"_Muon3_nTracks03","",10,0,10,"#mu_{3} ntracks03","Events");
  Muon1_sumPt03=HConfig.GetTH1D(Name+"_Muon1_sumPt03","",10,0,10,"","Events");
  Muon2_sumPt03=HConfig.GetTH1D(Name+"_Muon2_sumPt03","",10,0,10,"","Events");
  Muon3_sumPt03=HConfig.GetTH1D(Name+"_Muon3_sumPt03","",10,0,10,"","Events");
  Muon1_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt03","",10,0,10,"","Events");
  Muon2_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt03","",10,0,10,"","Events");
  Muon3_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon3_trackerVetoPt03","",10,0,10,"","Events");
  Muon1_emEt05=HConfig.GetTH1D(Name+"_Muon1_emEt05","",10,0,10,"#mu_{1} EM E_{T}05","Events");
  Muon2_emEt05=HConfig.GetTH1D(Name+"_Muon2_emEt05","",10,0,10,"#mu_{2} EM E_{T}05","Events");
  Muon3_emEt05=HConfig.GetTH1D(Name+"_Muon3_emEt05","",10,0,10,"#mu_{3} EM E_{T}05","Events");
  Muon1_emVetoEt05=HConfig.GetTH1D(Name+"_Muon1_emVetoEt05","#mu_{1} EM veto E_{T}05",10,0,10,"","Events");
  Muon2_emVetoEt05=HConfig.GetTH1D(Name+"_Muon2_emVetoEt05","#mu_{2} EM veto E_{T}05",10,0,10,"","Events");
  Muon3_emVetoEt05=HConfig.GetTH1D(Name+"_Muon3_emVetoEt05","#mu_{3} EM veto E_{T}05",10,0,10,"","Events");
  Muon1_hadEt05=HConfig.GetTH1D(Name+"_Muon1_hadEt05","",10,0,10,"#mu_{1} hadron E_{T}05","Events");
  Muon2_hadEt05=HConfig.GetTH1D(Name+"_Muon2_hadEt05","",10,0,10,"#mu_{2} hadron E_{T}05","Events");
  Muon3_hadEt05=HConfig.GetTH1D(Name+"_Muon3_hadEt05","",10,0,10,"#mu_{3} hadron E_{T}05","Events");
  Muon1_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt05","",10,0,10,"#mu_{1} hadron veto E_{T}05","Events");
  Muon2_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt05","",10,0,10,"#mu_{2} hadron veto E_{T}05","Events");
  Muon3_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon3_hadVetoEt05","",10,0,10,"#mu_{3} hadron veto E_{T}05","Events");
  Muon1_nJets05=HConfig.GetTH1D(Name+"_Muon1_nJets05","",10,0,10,"#mu_{1} njets 05","Events");
  Muon2_nJets05=HConfig.GetTH1D(Name+"_Muon2_nJets05","",10,0,10,"#mu_{2} njets 05","Events");
  Muon3_nJets05=HConfig.GetTH1D(Name+"_Muon3_nJets05","",10,0,10,"#mu_{3} njets 05","Events");
  Muon1_nTracks05=HConfig.GetTH1D(Name+"_Muon1_nTracks05","",10,0,10,"#mu_{1} ntracks05","Events");
  Muon2_nTracks05=HConfig.GetTH1D(Name+"_Muon2_nTracks05","",10,0,10,"#mu_{2} ntracks05","Events");
  Muon3_nTracks05=HConfig.GetTH1D(Name+"_Muon3_nTracks05","",10,0,10,"#mu_{3} ntracks05","Events");
  Muon1_sumPt05=HConfig.GetTH1D(Name+"_Muon1_sumPt05","",10,0,10,"#mu_{1} #Sigmap_{T}05","Events");
  Muon2_sumPt05=HConfig.GetTH1D(Name+"_Muon2_sumPt05","",10,0,10,"#mu_{2} #Sigmap_{T}05","Events");
  Muon3_sumPt05=HConfig.GetTH1D(Name+"_Muon3_sumPt05","",10,0,10,"#mu_{3} #Sigmap_{T}05","Events");
  Muon1_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt05","",10,0,10,"#mu_{1} tracker veto p_{T}05","Events");
  Muon2_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt05","",10,0,10,"#mu_{2} tracker veto p_{T}05","Events");
  Muon3_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon3_trackerVetoPt05","",10,0,10,"#mu_{3} tracker veto p_{T}05","Events");
  Muon1_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt03","#mu_{1} #Sigma charged had p_{T}03",10,0,10,"","Events");
  Muon2_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt03","#mu_{2} #Sigma charged had p_{T}03",10,0,10,"","Events");
  Muon3_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon3_sumChargedHadronPt03","#mu_{3} #Sigma charged had p_{T}03",10,0,10,"","Events");
  Muon1_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt03","",10,0,10,"#mu_{1} #Sigma charged particle p_{T}03","Events");
  Muon2_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt03","",10,0,10,"#mu_{2} #Sigma charged particle p_{T}03","Events");
  Muon3_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon3_sumChargedParticlePt03","",10,0,10,"#mu_{3} #Sigma charged particle p_{T}03","Events");
  Muon1_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt03","",10,0,10,"#mu_{1} #Sigma neutral had E_{T}03","Events");
  Muon2_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt03","",10,0,10,"#mu_{2} #Sigma neutral had E_{T}03","Events");
  Muon3_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon3_sumNeutralHadronEt03","",10,0,10,"#mu_{3} #Sigma neutral had E_{T}03","Events");
  Muon1_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold03","",10,0,10,"#mu_{1} #Sigma neutral had E_{T} HT03","Events");
  Muon2_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold03","",10,0,10,"#mu_{2} #Sigma neutral had E_{T} HT03","Events");
  Muon3_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon3_sumNeutralHadronEtHighThreshold03","",10,0,10,"#mu_{3} #Sigma neutral had E_{T} HT03","Events");
  Muon1_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt03","",10,0,10,"#mu_{1} #Sigma#gamma E_{T}03","Events");
  Muon2_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt03","",10,0,10,"#mu_{2} #Sigma#gamma E_{T}03","Events");
  Muon3_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon3_sumPhotonEt03","",10,0,10,"#mu_{3} #Sigma#gamma E_{T}03","Events");
  Muon1_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold03","",10,0,10,"#mu_{1} #Sigma#gamma E_{T}03","Events");
  Muon2_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold03","",10,0,10,"#mu_{2} #Sigma#gamma E_{T}03","Events");
  Muon3_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon3_sumPhotonEtHighThreshold03","",10,0,10,"#mu_{3} #Sigma#gamma E_{T}03","Events");
  Muon1_sumPUPt03=HConfig.GetTH1D(Name+"_Muon1_sumPUPt03","",10,0,10,"#mu_{1} #Sigma PU p_{T}","Events");
  Muon2_sumPUPt03=HConfig.GetTH1D(Name+"_Muon2_sumPUPt03","",10,0,10,"#mu_{2} #Sigma PU p_{T}","Events");
  Muon3_sumPUPt03=HConfig.GetTH1D(Name+"_Muon3_sumPUPt03","",10,0,10,"#mu_{3} #Sigma PU p_{T}","Events");
  Muon1_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt04","#mu_{1} #Sigma charged had p_{T}04",10,0,10,"","Events");
  Muon2_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt04","#mu_{2} #Sigma charged had p_{T}04",10,0,10,"","Events");
  Muon3_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon3_sumChargedHadronPt04","#mu_{3} #Sigma charged had p_{T}04",10,0,10,"","Events");
  Muon1_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt04","",10,0,10,"#mu_{1} #Sigma charged particle p_{T}04","Events");
  Muon2_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt04","",10,0,10,"#mu_{2} #Sigma charged particle p_{T}04","Events");
  Muon3_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon3_sumChargedParticlePt04","",10,0,10,"#mu_{3} #Sigma charged particle p_{T}04","Events");
  Muon1_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt04","",10,0,10,"#mu_{1} #Sigma neutral had E_{T}04","Events");
  Muon2_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt04","",10,0,10,"#mu_{2} #Sigma neutral had E_{T}04","Events");
  Muon3_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon3_sumNeutralHadronEt04","",10,0,10,"#mu_{3} #Sigma neutral had E_{T}04","Events");
  Muon1_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold04","",10,0,10,"#mu_{1} #Sigma neutral had E_{T} HT04","Events");
  Muon2_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold04","",10,0,10,"#mu_{2} #Sigma neutral had E_{T} HT04","Events");
  Muon3_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon3_sumNeutralHadronEtHighThreshold04","",10,0,10,"#mu_{3} #Sigma neutral had E_{T} HT04","Events");
  Muon1_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt04","",10,0,10,"#mu_{1} #Sigma#gammaE_{T}04","Events");
  Muon2_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt04","",10,0,10,"#mu_{2} #Sigma#gammaE_{T}04","Events");
  Muon3_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon3_sumPhotonEt04","",10,0,10,"#mu_{3} #Sigma#gammaE_{T}04","Events");
  Muon1_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold04","",10,0,10,"#mu_{1} #Sigma#gammaE_{T} HT04","Events");
  Muon2_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold04","",10,0,10,"#mu_{2} #Sigma#gammaE_{T} HT04","Events");
  Muon3_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon3_sumPhotonEtHighThreshold04","",10,0,10,"#mu_{3} #Sigma#gammaE_{T} HT04","Events");
  Muon1_sumPUPt04=HConfig.GetTH1D(Name+"_Muon1_sumPUPt04","",10,0,10,"#mu_{1} #SigmaPUPt04","Events");
  Muon2_sumPUPt04=HConfig.GetTH1D(Name+"_Muon2_sumPUPt04","",10,0,10,"#mu_{2} #SigmaPUPt04","Events");
  Muon3_sumPUPt04=HConfig.GetTH1D(Name+"_Muon3_sumPUPt04","",10,0,10,"#mu_{3} #sigmaPUPt04","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  Validation::Store_ExtraDist(){

  //
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&NTracks);
  Extradist1d.push_back(&HLT_Tau3Mu);  

  // ThreeMuon variables
  Extradist1d.push_back(&Muon1_isGlobal);
  Extradist1d.push_back(&Muon2_isGlobal);
  Extradist1d.push_back(&Muon3_isGlobal);
  Extradist1d.push_back(&Muon1_isStandAlone);
  Extradist1d.push_back(&Muon2_isStandAlone);
  Extradist1d.push_back(&Muon3_isStandAlone);
  Extradist1d.push_back(&Muon1_isTracker);
  Extradist1d.push_back(&Muon2_isTracker);
  Extradist1d.push_back(&Muon3_isTracker);
  Extradist1d.push_back(&Muon1_isCalo);
  Extradist1d.push_back(&Muon2_isCalo);
  Extradist1d.push_back(&Muon3_isCalo);
  Extradist1d.push_back(&Muon1_isIsolationValid);
  Extradist1d.push_back(&Muon2_isIsolationValid);
  Extradist1d.push_back(&Muon3_isIsolationValid);
  Extradist1d.push_back(&Muon1_isTimeValid);
  Extradist1d.push_back(&Muon2_isTimeValid);
  Extradist1d.push_back(&Muon3_isTimeValid);
  Extradist1d.push_back(&Muon1_emEt03);
  Extradist1d.push_back(&Muon2_emEt03);
  Extradist1d.push_back(&Muon3_emEt03);
  Extradist1d.push_back(&Muon1_emVetoEt03);
  Extradist1d.push_back(&Muon2_emVetoEt03);
  Extradist1d.push_back(&Muon3_emVetoEt03);
  Extradist1d.push_back(&Muon1_hadEt03);
  Extradist1d.push_back(&Muon2_hadEt03);
  Extradist1d.push_back(&Muon3_hadEt03);
  Extradist1d.push_back(&Muon1_hadVetoEt03);
  Extradist1d.push_back(&Muon2_hadVetoEt03);
  Extradist1d.push_back(&Muon3_hadVetoEt03);
  Extradist1d.push_back(&Muon1_nJets03);
  Extradist1d.push_back(&Muon2_nJets03);
  Extradist1d.push_back(&Muon3_nJets03);
  Extradist1d.push_back(&Muon1_nTracks03);
  Extradist1d.push_back(&Muon2_nTracks03);
  Extradist1d.push_back(&Muon3_nTracks03);
  Extradist1d.push_back(&Muon1_sumPt03);
  Extradist1d.push_back(&Muon2_sumPt03);
  Extradist1d.push_back(&Muon3_sumPt03);
  Extradist1d.push_back(&Muon1_trackerVetoPt03);
  Extradist1d.push_back(&Muon2_trackerVetoPt03);
  Extradist1d.push_back(&Muon3_trackerVetoPt03);
  Extradist1d.push_back(&Muon1_emEt05);
  Extradist1d.push_back(&Muon2_emEt05);
  Extradist1d.push_back(&Muon3_emEt05);
  Extradist1d.push_back(&Muon1_emVetoEt05);
  Extradist1d.push_back(&Muon2_emVetoEt05);
  Extradist1d.push_back(&Muon3_emVetoEt05);
  Extradist1d.push_back(&Muon1_hadEt05);
  Extradist1d.push_back(&Muon2_hadEt05);
  Extradist1d.push_back(&Muon3_hadEt05);
  Extradist1d.push_back(&Muon1_hadVetoEt05);
  Extradist1d.push_back(&Muon2_hadVetoEt05);
  Extradist1d.push_back(&Muon3_hadVetoEt05);
  Extradist1d.push_back(&Muon1_nJets05);
  Extradist1d.push_back(&Muon2_nJets05);
  Extradist1d.push_back(&Muon3_nJets05);
  Extradist1d.push_back(&Muon1_nTracks05);
  Extradist1d.push_back(&Muon2_nTracks05);
  Extradist1d.push_back(&Muon3_nTracks05);
  Extradist1d.push_back(&Muon1_sumPt05);
  Extradist1d.push_back(&Muon2_sumPt05);
  Extradist1d.push_back(&Muon3_sumPt05);
  Extradist1d.push_back(&Muon1_trackerVetoPt05);
  Extradist1d.push_back(&Muon2_trackerVetoPt05);
  Extradist1d.push_back(&Muon3_trackerVetoPt05);
  Extradist1d.push_back(&Muon1_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon2_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon3_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon1_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon2_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon3_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon3_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon3_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon1_sumPhotonEt03);
  Extradist1d.push_back(&Muon2_sumPhotonEt03);
  Extradist1d.push_back(&Muon3_sumPhotonEt03);
  Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon3_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon1_sumPUPt03);
  Extradist1d.push_back(&Muon2_sumPUPt03);
  Extradist1d.push_back(&Muon3_sumPUPt03);
  Extradist1d.push_back(&Muon1_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon2_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon3_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon1_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon2_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon3_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon3_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon3_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon1_sumPhotonEt04);
  Extradist1d.push_back(&Muon2_sumPhotonEt04);
  Extradist1d.push_back(&Muon3_sumPhotonEt04);
  Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon3_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon1_sumPUPt04);
  Extradist1d.push_back(&Muon2_sumPUPt04);
  Extradist1d.push_back(&Muon3_sumPUPt04);
  
  Extradist1d.push_back(&TripleMass);

}

void  Validation::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  bool l1_pass = false;
  bool hlt_pass = false;

  for (int j=0; j<Ntp->NL1Seeds(); ++j){
  	auto l1_name = Ntp->L1Name(j);
	if (l1_name.find("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4")!=std::string::npos || l1_name.find("L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17")!=std::string::npos) {l1_pass = true;}
  }
  
  for (int j=0; j<Ntp->NHLT(); ++j){
	auto hlt_name = Ntp->HLTName(j);
	if (hlt_name.find("HLT_DoubleMu3_Trk")!=std::string::npos) { hlt_pass = true; }
  }
  
  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(L1TOk)=l1_pass;
  pass.at(L1TOk)=value.at(L1TOk);
  
  value.at(HLTOk)=hlt_pass;
  pass.at(HLTOk)=value.at(HLTOk);

  double wobs=1;
  double w;
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */}
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
  NVtx.at(t).Fill(Ntp->NVtx(),w);
  NTracks.at(t).Fill(Ntp->NTracks(),w);
  
   if (Ntp->NThreeMuons()>0)
	for(unsigned int i3M=0; i3M < Ntp->NThreeMuons(); i3M++){
       unsigned int mu1 = Ntp-> ThreeMuonIndices(i3M).at(0);
       unsigned int mu2 = Ntp-> ThreeMuonIndices(i3M).at(1);
       unsigned int mu3 = Ntp-> ThreeMuonIndices(i3M).at(2);
   	 
   	 //Three muon variables
   	 Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
       Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
       Muon3_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu3),w);
       Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
       Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
       Muon3_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu3),w);
       Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
       Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
       Muon3_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu3),w);
       Muon1_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu1),w);
       Muon2_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu2),w);
       Muon3_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu3),w);
       Muon1_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu1),w);
       Muon2_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu2),w);
       Muon3_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu3),w);
       Muon1_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu1),w);
       Muon2_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu2),w);
       Muon3_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu3),w);
       Muon1_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu1),w);
       Muon2_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu2),w);
       Muon3_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu3),w);
       Muon1_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu1),w);
       Muon2_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu2),w);
       Muon3_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu3),w);
       Muon1_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu1),w);
       Muon2_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu2),w);
       Muon3_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu3),w);
       Muon1_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu1),w);
       Muon2_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu2),w);
       Muon3_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu3),w);
       Muon1_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu1),w);
       Muon2_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu2),w);
       Muon3_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu3),w);
       Muon1_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu1),w);
       Muon2_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu2),w);
       Muon3_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu3),w);
       Muon1_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu1),w);
       Muon2_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu2),w);
       Muon3_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu3),w);
       Muon1_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu1),w);
       Muon2_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu2),w);
       Muon3_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu3),w);
       Muon1_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu1),w);
       Muon2_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu2),w);
       Muon3_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu3),w);
       Muon1_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu1),w);
       Muon2_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu2),w);
       Muon3_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu3),w);
       Muon1_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu1),w);
       Muon2_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu2),w);
       Muon3_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu3),w);
       Muon1_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu1),w);
       Muon2_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu2),w);
       Muon3_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu3),w);
       Muon1_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu1),w);
       Muon2_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu2),w);
       Muon3_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu3),w);
       Muon1_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu1),w);
       Muon2_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu2),w);
       Muon3_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu3),w);
       Muon1_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu1),w);
       Muon2_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu2),w);
       Muon3_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu3),w);
       Muon1_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu1),w);
       Muon2_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu2),w);
       Muon3_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu3),w);
       Muon1_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu1),w);
       Muon2_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu2),w);
       Muon3_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu3),w);
       Muon1_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu1),w);
       Muon2_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu2),w);
       Muon3_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu3),w);
       Muon1_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu1),w);
       Muon2_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu2),w);
       Muon3_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu3),w);
       Muon1_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu1),w);
       Muon2_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu2),w);
       Muon3_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu3),w);
       Muon1_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu1),w);
       Muon2_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu2),w);
       Muon3_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu3),w);
       Muon1_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu1),w);
       Muon2_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu2),w);
       Muon3_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu3),w);
       Muon1_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu1),w);
       Muon2_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu2),w);
       Muon3_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu3),w);
       Muon1_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu1),w);
       Muon2_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu2),w);
       Muon3_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu3),w);
       Muon1_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu1),w);
       Muon2_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu2),w);
       Muon3_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu3),w);
       Muon1_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu1),w);
       Muon2_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu2),w);
       Muon3_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu3),w);
       Muon1_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu1),w);
       Muon2_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu2),w);
       Muon3_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu3),w);
       Muon1_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu1),w);
       Muon2_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu2),w);
       Muon3_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu3),w);
       Muon1_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu1),w);
       Muon2_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu2),w);
       Muon3_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu3),w);
       Muon1_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu1),w);
       Muon2_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu2),w);
       Muon3_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu3),w);
   	 

	 Muon1Muon2dR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(0)).Eta(),Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(0)).Phi(),Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(1)).Eta(),Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(1)).Phi()));
	 Muon2Muon3dR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(0)).Eta(),Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(0)).Phi(),Ntp->Track_P4(Ntp->ThreeMuonIndices(i3M).at(2)).Eta(),Ntp->Track_P4(Ntp->ThreeMuonIndices(i3M).at(2)).Phi()));
	 Muon1Muon3dR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(1)).Eta(),Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(1)).Phi(),Ntp->Track_P4(Ntp->ThreeMuonIndices(i3M).at(2)).Eta(),Ntp->Track_P4(Ntp->ThreeMuonIndices(i3M).at(2)).Phi()));
    TripleMass.at(t).Fill((Ntp->Muon_P4(Ntp->ThreeMuonIndices(i3M).at(0))  + Ntp->Muon_P4(Ntp-> ThreeMuonIndices(i3M).at(1))+Ntp->Track_P4(Ntp->ThreeMuonIndices(i3M).at(2))).M(), w);
	 }
  }
  }


void  Validation::Finish(){
  Selection::Finish();
}
