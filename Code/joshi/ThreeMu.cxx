#include "ThreeMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

double ThreeMu::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}
ThreeMu::ThreeMu(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

ThreeMu::~ThreeMu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }

  Logger(Logger::Info) << "complete." << std::endl;
}

void  ThreeMu::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1SeedOk)     cut.at(L1SeedOk)=1;
    if(i==HLTOk)        cut.at(HLTOk)=1;
    if(i==isThreeMu)        cut.at(isThreeMu)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
	 if(i==fitVtxChiSq)	cut.at(fitVtxChiSq)=5.0;
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
    }
    else if(i==L1SeedOk){
      title.at(i)="L1 seed ";
      hlabel="L1 triggers";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==isThreeMu){
      title.at(i)="isThreeMu ";
      hlabel="3mu category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isThreeMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isThreeMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
	 
	 else if(i==fitVtxChiSq){
      title.at(i)="Normalized chi sq fit vertex $(>5)$ ";
      hlabel="Normalized chi sq fit vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
    }
  } 

// Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  //Muon variables (Muons from dimuon + track candidates)
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
  
  Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","",10,0,10,"#mu_{1} Iso ntrks","Events");
  Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","",10,0,10,"#mu_{1} Iso rel p_{T}","Events");
  Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","",10,0,10,"#mu_{1} Iso MinDist","Events");
  Isolation05_RelPt=HConfig.GetTH1D(Name+"_Isolation05_RelPt","",10,0,10,"#mu_{1} Iso05 rel p_{T}","Events");
  Isolation05_NTracks=HConfig.GetTH1D(Name+"_Isolation05_NTracks","",10,0,10,"#mu_{1} Iso05 ntrks","Events");
  Isolation05_MinDist=HConfig.GetTH1D(Name+"_Isolation05_MinDist","",10,0,10,"#mu_{1} Iso05 MinDist","Events");
  Isolation_Ntrk1=HConfig.GetTH1D(Name+"_Isolation_Ntrk1","",10,0,10,"#mu_{1} Iso ntrk 1","Events");
  Isolation_Ntrk2=HConfig.GetTH1D(Name+"_Isolation_Ntrk2","",10,0,10,"#mu_{1} Iso ntrk 2","Events");
  Isolation_Ntrk3=HConfig.GetTH1D(Name+"_Isolation_Ntrk3","",10,0,10,"#mu_{1} Iso ntrk 3","Events");
  Isolation_Ntrk0p1=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p1","",10,0,10,"#mu_{1} Iso ntrk0p1","Events");
  Isolation_Ntrk0p2=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p2","",10,0,10,"#mu_{1} Iso ntrk0p2","Events");
  Isolation_Ntrk0p5=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p5","",10,0,10,"#mu_{1} Iso ntrk0p5","Events");
  Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","",10,0,10,"#mu_{1} Iso max(dxy)","Events");
  
  //Dimuon Information (Muons from dimuon + track candidates)
  Muon1Muon3dR=HConfig.GetTH1D(Name+"_Muon1Muon3dR","dR between the highest p muon and the lowest pt muon",100,0,5,"dR","Events");
  Muon2Muon3dR=HConfig.GetTH1D(Name+"_Muon2Muon3dR","dR between the second highest p muon and the lowest pt muon",100,0,5,"dR","Events");
  Muon1Muon2dR=HConfig.GetTH1D(Name+"_Muon1Muon1dR","dR between the highest p muon and the second highest p muon",100,0,5,"dR","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu#mu mass",50,1.6,2.1,"Mass of the #mu#mu#mu","Events");

Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  ThreeMu::Store_ExtraDist(){ 
  Extradist1d.push_back(&Muon1Muon3dR);
  Extradist1d.push_back(&Muon2Muon3dR);
  Extradist1d.push_back(&Muon1Muon2dR);
  Extradist1d.push_back(&TripleMass);
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
  
  Extradist1d.push_back(&Isolation_NTracks);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_MinDist);
  Extradist1d.push_back(&Isolation05_RelPt);
  Extradist1d.push_back(&Isolation05_NTracks);
  Extradist1d.push_back(&Isolation05_MinDist);
  Extradist1d.push_back(&Isolation_Ntrk1);
  Extradist1d.push_back(&Isolation_Ntrk2);
  Extradist1d.push_back(&Isolation_Ntrk3);
  Extradist1d.push_back(&Isolation_Ntrk0p1);
  Extradist1d.push_back(&Isolation_Ntrk0p2);
  Extradist1d.push_back(&Isolation_Ntrk0p5);
  Extradist1d.push_back(&Isolation_maxdxy);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  Extradist1d.push_back(&NVtx);
}

void  ThreeMu::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  value.at(L1SeedOk) = 0;
  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  for(int l1iTrigger=0; l1iTrigger < Ntp->NL1Seeds(); l1iTrigger++){
    TString L1 = Ntp->L1Name(l1iTrigger);
    if(L1.Contains("L1_DoubleMu0er1p4_dEta_Max1p8_OS") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_DoubleMu_10_0_dEta_Max1p8") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_TripleMu0") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
  }

  value.at(PrimeVtx)=Ntp->NVtx();
  value.at(fitVtxChiSq)=0;

  value.at(isThreeMu) = 0;
  if(Ntp->NTwoMuonsTrack()==0 && Ntp->NThreeMuons() != 0) value.at(isThreeMu) = 1;

  pass.at(isThreeMu) = (value.at(isThreeMu) == cut.at(isThreeMu));
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk));

  int mu1=-1, mu2=-1, mu3=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999;

  if (value.at(isThreeMu)==1){
  	for (unsigned int i3M=0; i3M<Ntp->NThreeMuons(); i3M++){
		int tmp_mu1 =  Ntp->ThreeMuonIndices(i3M).at(0);
      int tmp_mu2 =  Ntp->ThreeMuonIndices(i3M).at(1);
      int tmp_mu3 =  Ntp->ThreeMuonIndices(i3M).at(2);
		if (tmp_chisq>Ntp->ThreeMuons_SV_Chi2(i3M)){
			tmp_idx = i3M;
			tmp_chisq = Ntp->ThreeMuons_SV_Chi2(i3M);
			mu1 = tmp_mu1;
			mu2 = tmp_mu2;
			mu3 = tmp_mu3;
		}
	}
  }
  if (tmp_idx==-1) value.at(fitVtxChiSq)=999.0;
  else value.at(fitVtxChiSq) = tmp_chisq;
  pass.at(fitVtxChiSq) = (value.at(fitVtxChiSq)<cut.at(fitVtxChiSq));
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  if(status){
  if (DEBUG) cout<<mu1<<" "<<mu2<<" "<<mu3<<" "<<Ntp->NMuons()<<endl;
	 NVtx.at(t).Fill(Ntp->NVtx(),w);
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
	 
	 Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(tmp_idx),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(tmp_idx),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(tmp_idx),w);
    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(tmp_idx),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(tmp_idx),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(tmp_idx),w);
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(tmp_idx),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(tmp_idx),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(tmp_idx),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(tmp_idx),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(tmp_idx),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(tmp_idx),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(tmp_idx),w);
	 
	 TripleMass.at(t).Fill((Ntp->Muon_P4(mu1)+Ntp->Muon_P4(mu2)+Ntp->Muon_P4(mu3)).M(), w);
/* 
//	 Muon1Muon2dR.at(t).Fill(deltaR(Ntp->Muon_P4(mu1).Eta(),Ntp->Muon_P4(mu1).Phi(),Ntp->Track_P4(mu2).Eta(),Ntp->Track_P4(mu2).Phi()));
//	 Muon2Muon3dR.at(t).Fill(deltaR(Ntp->Muon_P4(mu2).Eta(),Ntp->Muon_P4(mu2).Phi(),Ntp->Track_P4(mu3).Eta(),Ntp->Track_P4(mu3).Phi()));
//	 Muon1Muon3dR.at(t).Fill(deltaR(Ntp->Muon_P4(mu1).Eta(),Ntp->Muon_P4(mu1).Phi(),Ntp->Track_P4(mu3).Eta(),Ntp->Track_P4(mu3).Phi()));
cout<<mu3<<endl;

    if(Ntp->NThreeMuons()!=0){
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);

      double pt1 = Ntp->Muon_P4(Muon_index_1).Pt();
      double pt2 = Ntp->Muon_P4(Muon_index_2).Pt();
      double pt3 = Ntp->Muon_P4(Muon_index_3).Pt();
      FirstMuonsPt.at(t).Fill(pt1,1);
      SecondMuonsPt.at(t).Fill(pt2,1);
      ThirdMuonsPt.at(t).Fill(pt3,1);

      FirstMuonsEta.at(t).Fill( Ntp->Muon_P4(Muon_index_1).Eta(),1);
      SecondMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_2).Eta(),1);
      ThirdMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_3).Eta(),1);
    }

    double deltaMass(999.);
    unsigned int pair_index(0);
    if(Ntp->NTwoMuonsTrack()!=0){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

      unsigned int muon_1 =  Ntp-> ThreeMuonsMuonIndices(i2M).at(0);
      unsigned int muon_2 =  Ntp-> ThreeMuonsMuonIndices(i2M).at(1);

      if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
	deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
	pair_index = i2M; // this is an index of the candidate with best mumu mass
      }
    }
    
    PhiMass.at(t).Fill((Ntp->Muon_P4( Ntp-> ThreeMuonsMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> ThreeMuonsMuonIndices(pair_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4( Ntp-> ThreeMuonsMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> ThreeMuonsMuonIndices(pair_index).at(1))+ 
			   Ntp->Track_P4(Ntp->ThreeMuonsMuonIndices(pair_index).at(0))).M(), w);

    double phimass = (Ntp->Muon_P4( Ntp-> ThreeMuonsMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> ThreeMuonsMuonIndices(pair_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> ThreeMuonsMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> ThreeMuonsMuonIndices(pair_index).at(1))+
		     Ntp->Track_P4(Ntp->ThreeMuonsMuonIndices(pair_index).at(0))).M();

    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

}*/
  }
}


void  ThreeMu::Finish(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(mode == RECONSTRUCT){
    for(unsigned int i=0; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = 1/Nminus0.at(0).at(i).Integral();
      ScaleAllHistOfType(HConfig.GetType(i),scale);
    }
  }


  Selection::Finish();

}
