#include "TripleMuIsolation.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

TripleMuIsolation::TripleMuIsolation(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.731),
  tauMaxMass_(1.823),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

TripleMuIsolation::~TripleMuIsolation(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  TripleMuIsolation::Configure(){
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=2.5;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.5;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.5;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=0.03;
    if(i==ThreeMuMass)        cut.at(ThreeMuMass)=1;// true for MC and mass side band for data
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
    else if(i==Mu1PtCut){
      title.at(i)="Mu1 Pt";
      hlabel="Muon1 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
    }
    else if(i==Mu2PtCut){
      title.at(i)="Mu2 Pt";
      hlabel="Muon2 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
    }
    else if(i==Mu3PtCut){
      title.at(i)="Mu3 Pt";
      hlabel="Muon3 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
    }
    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="pass MuonID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PhiVeto){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
    }
    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==ThreeMuMass){
      title.at(i)="Tau Mass";
      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
    }


  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
 
  // New isolation based histograms
Muon1_Isolation_emEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_emEt03","Muon1_Isolation_emEt03",50,0,25,"#mu_{1} ECAL deposits in Iso03 (GeV)","Events");
  Muon1_Isolation_emVetoEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_emVetoEt03","Muon1_Isolation_emVetoEt03",50,0,10,"#mu_{1} ECAL deposits in Iso03 Veto (GeV)","Events");
  Muon1_Isolation_hadEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_hadEt03","Muon1_Isolation_hadEt03",30,0,15,"#mu_{1} HCAL deposits in Iso03 (GeV)","Events");
  Muon1_Isolation_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_hadVetoEt03","Muon1_Isolation_hadVetoEt03",50,0,10,"#mu_{1} HCAL deposits in Iso03 Veto (GeV)","Events");
  Muon1_Isolation_nJets03=HConfig.GetTH1D(Name+"_Muon1_Isolation_nJets03","Muon1_Isolation_nJets03",5,0,5,"#mu_{1} number of Jets in Iso03","Events");
  Muon1_Isolation_nTracks03=HConfig.GetTH1D(Name+"_Muon1_Isolation_nTracks03","Muon1_Isolation_nTracks03",20,0,20,"#mu_{1} number of tracks in ECAL Iso03","Events");
  Muon1_Isolation_sumPt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPt03","Muon1_Isolation_sumPt03",50,0,50,"#mu_{1} #Sigmap_{T} in Iso03 (GeV)","Events");
  Muon1_Isolation_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_trackerVetoPt03","Muon1_Isolation_trackerVetoPt03",10,0,10,"#mu_{1} p_{T} in tracker Iso03 Veto (GeV)","Events");
  Muon1_Isolation_emEt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_emEt05","Muon1_Isolation_emEt05",50,0,25,"#mu_{1} ECAL deposits in Iso05 (GeV)","Events");
  Muon1_Isolation_emVetoEt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_emVetoEt05","Muon1_Isolation_emVetoEt05",50,0,10,"#mu_{1} ECAL deposits in Iso05 Veto (GeV)","Events");
  Muon1_Isolation_hadEt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_hadEt05","Muon1_Isolation_hadEt05",50,0,25,"#mu_{1} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon1_Isolation_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_hadVetoEt05","Muon1_Isolation_hadVetoEt05",50,0,10,"#mu_{1} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon1_Isolation_nJets05=HConfig.GetTH1D(Name+"_Muon1_Isolation_nJets05","Muon1_Isolation_nJets05",5,0,5,"#mu_{1} number of Jets in Iso05","Events");
  Muon1_Isolation_nTracks05=HConfig.GetTH1D(Name+"_Muon1_Isolation_nTracks05","Muon1_Isolation_nTracks05",20,0,20,"#mu_{1} number of tracks in Iso05","Events");
  Muon1_Isolation_sumPt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPt05","Muon1_Isolation_sumPt05",50,0,50,"#mu_{1} #Sigmap_{T} in Iso05 (GeV)","Events");
  Muon1_Isolation_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon1_Isolation_trackerVetoPt05","Muon1_Isolation_trackerVetoPt05",50,0,25,"#mu_{1} p_{T} in the tracker Iso05 Veto (GeV)","Events");
  Muon1_Isolation_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumChargedHadronPt03","Muon1_Isolation_sumChargedHadronPt03",30,0,15,"#mu_{1} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon1_Isolation_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumChargedParticlePt03","Muon1_Isolation_sumChargedParticlePt03",30,0,15,"#mu_{1} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon1_Isolation_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumNeutralHadronEt03","Muon1_Isolation_sumNeutralHadronEt03",30,0,15,"#mu_{1} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon1_Isolation_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumNeutralHadronEtHighThreshold03","Muon1_Isolation_sumNeutralHadronEtHighThreshold03",30,0,15,"#mu_{1} #SigmaE_{T} (Neutral Hadron) (HT) in Iso03 (GeV)","Events");
  Muon1_Isolation_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPhotonEt03","Muon1_Isolation_sumPhotonEt03",30,0,15,"#mu_{1} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon1_Isolation_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPhotonEtHighThreshold03","Muon1_Isolation_sumPhotonEtHighThreshold03",30,0,15,"#mu_{1} #SigmaE_{T}#gamma (high threshold) in Iso03 (GeV)","Events");
  Muon1_Isolation_sumPUPt03=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPUPt03","Muon1_Isolation_sumPUPt03",30,0,15,"#mu_{1} #SigmaPUp_{T} in Iso03 (GeV)","Events");
  Muon1_Isolation_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumChargedHadronPt04","Muon1_Isolation_sumChargedHadronPt04",30,0,15,"#mu_{1} #Sigmap_{T} (charged hadrons) in Iso04 (GeV)","Events");
  Muon1_Isolation_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumChargedParticlePt04","Muon1_Isolation_sumChargedParticlePt04",30,0,15,"#mu_{1} #Sigmap_{T} (charged particles) (GeV)","Events");
  Muon1_Isolation_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumNeutralHadronEt04","Muon1_Isolation_sumNeutralHadronEt04",30,0,15,"#mu_{1} #SigmaE_{T} (neutral hadrons) in Iso04 (GeV)","Events");
  Muon1_Isolation_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumNeutralHadronEtHighThreshold04","Muon1_Isolation_sumNeutralHadronEtHighThreshold04",30,0,15,"#mu_{1} #SigmaE_{T} (neutral hadrons) (HT) in Iso04 (GeV)","Events");
  Muon1_Isolation_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPhotonEt04","Muon1_Isolation_sumPhotonEt04",100,0,50,"#mu_{1} #SigmaE_{T} #gamma in Iso04 (GeV)","Events");
  Muon1_Isolation_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPhotonEtHighThreshold04","Muon1_Isolation_sumPhotonEtHighThreshold04",30,0,15,"#mu_{1} #SigmaE_{T}#gamma (high threshold) in Iso04 (GeV)","Events");
  Muon1_Isolation_sumPUPt04=HConfig.GetTH1D(Name+"_Muon1_Isolation_sumPUPt04","Muon1_Isolation_sumPUPt04",30,0,15,"#mu_{1} #SigmaPUp_{T} in Iso04 (GeV)","Events");
  
  // ------------------------------------------------------------------------------------------------------------------------

Muon2_Isolation_emEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_emEt03","Muon2_Isolation_emEt03",50,0,25,"#mu_{2} ECAL deposits in Iso03 (GeV)","Events");
  Muon2_Isolation_emVetoEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_emVetoEt03","Muon2_Isolation_emVetoEt03",50,0,10,"#mu_{2} ECAL deposits in Iso03 Veto (GeV)","Events");
  Muon2_Isolation_hadEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_hadEt03","Muon2_Isolation_hadEt03",30,0,15,"#mu_{2} HCAL deposits in Iso03 (GeV)","Events");
  Muon2_Isolation_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_hadVetoEt03","Muon2_Isolation_hadVetoEt03",50,0,10,"#mu_{2} HCAL deposits in Iso03 Veto (GeV)","Events");
  Muon2_Isolation_nJets03=HConfig.GetTH1D(Name+"_Muon2_Isolation_nJets03","Muon2_Isolation_nJets03",5,0,5,"#mu_{2} number of Jets in Iso03","Events");
  Muon2_Isolation_nTracks03=HConfig.GetTH1D(Name+"_Muon2_Isolation_nTracks03","Muon2_Isolation_nTracks03",20,0,20,"#mu_{2} number of tracks in ECAL Iso03","Events");
  Muon2_Isolation_sumPt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPt03","Muon2_Isolation_sumPt03",50,0,50,"#mu_{2} #Sigmap_{T} in Iso03 (GeV)","Events");
  Muon2_Isolation_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_trackerVetoPt03","Muon2_Isolation_trackerVetoPt03",10,0,10,"#mu_{2} p_{T} in tracker Iso03 Veto (GeV)","Events");
  Muon2_Isolation_emEt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_emEt05","Muon2_Isolation_emEt05",50,0,25,"#mu_{2} ECAL deposits in Iso05 (GeV)","Events");
  Muon2_Isolation_emVetoEt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_emVetoEt05","Muon2_Isolation_emVetoEt05",50,0,10,"#mu_{2} ECAL deposits in Iso05 Veto (GeV)","Events");
  Muon2_Isolation_hadEt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_hadEt05","Muon2_Isolation_hadEt05",50,0,25,"#mu_{2} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon2_Isolation_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_hadVetoEt05","Muon2_Isolation_hadVetoEt05",50,0,10,"#mu_{2} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon2_Isolation_nJets05=HConfig.GetTH1D(Name+"_Muon2_Isolation_nJets05","Muon2_Isolation_nJets05",5,0,5,"#mu_{2} number of Jets in Iso05","Events");
  Muon2_Isolation_nTracks05=HConfig.GetTH1D(Name+"_Muon2_Isolation_nTracks05","Muon2_Isolation_nTracks05",20,0,20,"#mu_{2} number of tracks in Iso05","Events");
  Muon2_Isolation_sumPt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPt05","Muon2_Isolation_sumPt05",50,0,50,"#mu_{2} #Sigmap_{T} in Iso05 (GeV)","Events");
  Muon2_Isolation_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon2_Isolation_trackerVetoPt05","Muon2_Isolation_trackerVetoPt05",50,0,25,"#mu_{2} p_{T} in the tracker Iso05 Veto (GeV)","Events");
  Muon2_Isolation_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumChargedHadronPt03","Muon2_Isolation_sumChargedHadronPt03",30,0,15,"#mu_{2} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon2_Isolation_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumChargedParticlePt03","Muon2_Isolation_sumChargedParticlePt03",30,0,15,"#mu_{2} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon2_Isolation_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumNeutralHadronEt03","Muon2_Isolation_sumNeutralHadronEt03",30,0,15,"#mu_{2} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon2_Isolation_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumNeutralHadronEtHighThreshold03","Muon2_Isolation_sumNeutralHadronEtHighThreshold03",30,0,15,"#mu_{2} #SigmaE_{T} (Neutral Hadron) (HT) in Iso03 (GeV)","Events");
  Muon2_Isolation_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPhotonEt03","Muon2_Isolation_sumPhotonEt03",30,0,15,"#mu_{2} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon2_Isolation_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPhotonEtHighThreshold03","Muon2_Isolation_sumPhotonEtHighThreshold03",30,0,15,"#mu_{2} #SigmaE_{T}#gamma (high threshold) in Iso03 (GeV)","Events");
  Muon2_Isolation_sumPUPt03=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPUPt03","Muon2_Isolation_sumPUPt03",30,0,15,"#mu_{2} #SigmaPUp_{T} in Iso03 (GeV)","Events");
  Muon2_Isolation_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumChargedHadronPt04","Muon2_Isolation_sumChargedHadronPt04",30,0,15,"#mu_{2} #Sigmap_{T} (charged hadrons) in Iso04 (GeV)","Events");
  Muon2_Isolation_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumChargedParticlePt04","Muon2_Isolation_sumChargedParticlePt04",30,0,15,"#mu_{2} #Sigmap_{T} (charged particles) (GeV)","Events");
  Muon2_Isolation_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumNeutralHadronEt04","Muon2_Isolation_sumNeutralHadronEt04",30,0,15,"#mu_{2} #SigmaE_{T} (neutral hadrons) in Iso04 (GeV)","Events");
  Muon2_Isolation_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumNeutralHadronEtHighThreshold04","Muon2_Isolation_sumNeutralHadronEtHighThreshold04",30,0,15,"#mu_{2} #SigmaE_{T} (neutral hadrons) (HT) in Iso04 (GeV)","Events");
  Muon2_Isolation_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPhotonEt04","Muon2_Isolation_sumPhotonEt04",100,0,50,"#mu_{2} #SigmaE_{T} #gamma in Iso04 (GeV)","Events");
  Muon2_Isolation_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPhotonEtHighThreshold04","Muon2_Isolation_sumPhotonEtHighThreshold04",30,0,15,"#mu_{2} #SigmaE_{T}#gamma (high threshold) in Iso04 (GeV)","Events");
  Muon2_Isolation_sumPUPt04=HConfig.GetTH1D(Name+"_Muon2_Isolation_sumPUPt04","Muon2_Isolation_sumPUPt04",30,0,15,"#mu_{2} #SigmaPUp_{T} in Iso04 (GeV)","Events");

  // ------------------------------------------------------------------------------------------------------------------------

Muon3_Isolation_emEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_emEt03","Muon3_Isolation_emEt03",50,0,25,"#mu_{3} ECAL deposits in Iso03 (GeV)","Events");
  Muon3_Isolation_emVetoEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_emVetoEt03","Muon3_Isolation_emVetoEt03",50,0,10,"#mu_{3} ECAL deposits in Iso03 Veto (GeV)","Events");
  Muon3_Isolation_hadEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_hadEt03","Muon3_Isolation_hadEt03",30,0,15,"#mu_{3} HCAL deposits in Iso03 (GeV)","Events");
  Muon3_Isolation_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_hadVetoEt03","Muon3_Isolation_hadVetoEt03",50,0,10,"#mu_{3} HCAL deposits in Iso03 Veto (GeV)","Events");
  Muon3_Isolation_nJets03=HConfig.GetTH1D(Name+"_Muon3_Isolation_nJets03","Muon3_Isolation_nJets03",5,0,5,"#mu_{3} number of Jets in Iso03","Events");
  Muon3_Isolation_nTracks03=HConfig.GetTH1D(Name+"_Muon3_Isolation_nTracks03","Muon3_Isolation_nTracks03",20,0,20,"#mu_{3} number of tracks in ECAL Iso03","Events");
  Muon3_Isolation_sumPt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPt03","Muon3_Isolation_sumPt03",50,0,50,"#mu_{3} #Sigmap_{T} in Iso03 (GeV)","Events");
  Muon3_Isolation_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_trackerVetoPt03","Muon3_Isolation_trackerVetoPt03",10,0,10,"#mu_{3} p_{T} in tracker Iso03 Veto (GeV)","Events");
  Muon3_Isolation_emEt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_emEt05","Muon3_Isolation_emEt05",50,0,25,"#mu_{3} ECAL deposits in Iso05 (GeV)","Events");
  Muon3_Isolation_emVetoEt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_emVetoEt05","Muon3_Isolation_emVetoEt05",50,0,10,"#mu_{3} ECAL deposits in Iso05 Veto (GeV)","Events");
  Muon3_Isolation_hadEt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_hadEt05","Muon3_Isolation_hadEt05",50,0,25,"#mu_{3} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon3_Isolation_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_hadVetoEt05","Muon3_Isolation_hadVetoEt05",50,0,10,"#mu_{3} HCAL deposits in Iso05 Veto (GeV)","Events");
  Muon3_Isolation_nJets05=HConfig.GetTH1D(Name+"_Muon3_Isolation_nJets05","Muon3_Isolation_nJets05",5,0,5,"#mu_{3} number of Jets in Iso05","Events");
  Muon3_Isolation_nTracks05=HConfig.GetTH1D(Name+"_Muon3_Isolation_nTracks05","Muon3_Isolation_nTracks05",20,0,20,"#mu_{3} number of tracks in Iso05","Events");
  Muon3_Isolation_sumPt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPt05","Muon3_Isolation_sumPt05",50,0,50,"#mu_{3} #Sigmap_{T} in Iso05 (GeV)","Events");
  Muon3_Isolation_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon3_Isolation_trackerVetoPt05","Muon3_Isolation_trackerVetoPt05",50,0,25,"#mu_{3} p_{T} in the tracker Iso05 Veto (GeV)","Events");
  Muon3_Isolation_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumChargedHadronPt03","Muon3_Isolation_sumChargedHadronPt03",30,0,15,"#mu_{3} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon3_Isolation_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumChargedParticlePt03","Muon3_Isolation_sumChargedParticlePt03",30,0,15,"#mu_{3} #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Muon3_Isolation_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumNeutralHadronEt03","Muon3_Isolation_sumNeutralHadronEt03",30,0,15,"#mu_{3} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon3_Isolation_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumNeutralHadronEtHighThreshold03","Muon3_Isolation_sumNeutralHadronEtHighThreshold03",30,0,15,"#mu_{3} #SigmaE_{T} (Neutral Hadron) (HT) in Iso03 (GeV)","Events");
  Muon3_Isolation_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPhotonEt03","Muon3_Isolation_sumPhotonEt03",30,0,15,"#mu_{3} #SigmaE_{T} in Iso03 (GeV)","Events");
  Muon3_Isolation_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPhotonEtHighThreshold03","Muon3_Isolation_sumPhotonEtHighThreshold03",30,0,15,"#mu_{3} #SigmaE_{T}#gamma (high threshold) in Iso03 (GeV)","Events");
  Muon3_Isolation_sumPUPt03=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPUPt03","Muon3_Isolation_sumPUPt03",30,0,15,"#mu_{3} #SigmaPUp_{T} in Iso03 (GeV)","Events");
  Muon3_Isolation_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumChargedHadronPt04","Muon3_Isolation_sumChargedHadronPt04",30,0,15,"#mu_{3} #Sigmap_{T} (charged hadrons) in Iso04 (GeV)","Events");
  Muon3_Isolation_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumChargedParticlePt04","Muon3_Isolation_sumChargedParticlePt04",30,0,15,"#mu_{3} #Sigmap_{T} (charged particles) (GeV)","Events");
  Muon3_Isolation_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumNeutralHadronEt04","Muon3_Isolation_sumNeutralHadronEt04",30,0,15,"#mu_{3} #SigmaE_{T} (neutral hadrons) in Iso04 (GeV)","Events");
  Muon3_Isolation_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumNeutralHadronEtHighThreshold04","Muon3_Isolation_sumNeutralHadronEtHighThreshold04",30,0,15,"#mu_{3} #SigmaE_{T} (neutral hadrons) (HT) in Iso04 (GeV)","Events");
  Muon3_Isolation_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPhotonEt04","Muon3_Isolation_sumPhotonEt04",100,0,50,"#mu_{3} #SigmaE_{T} #gamma in Iso04 (GeV)","Events");
  Muon3_Isolation_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPhotonEtHighThreshold04","Muon3_Isolation_sumPhotonEtHighThreshold04",30,0,15,"#mu_{3} #SigmaE_{T}#gamma (high threshold) in Iso04 (GeV)","Events");
  Muon3_Isolation_sumPUPt04=HConfig.GetTH1D(Name+"_Muon3_Isolation_sumPUPt04","Muon3_Isolation_sumPUPt04",30,0,15,"#mu_{3} #SigmaPUp_{T} in Iso04 (GeV)","Events");

  // ------------------------------------------------------------------------------------------------------------------------

  Min_Isolation_emEt03=HConfig.GetTH1D(Name+"_Min_Isolation_emEt03","Min_Isolation_emEt03",50,0,25,"Min ECAL deposits in Iso03 (GeV)","Events");
  Min_Isolation_emVetoEt03=HConfig.GetTH1D(Name+"_Min_Isolation_emVetoEt03","Min_Isolation_emVetoEt03",50,0,10,"Min ECAL deposits in Iso03 Veto (GeV)","Events");
  Min_Isolation_hadEt03=HConfig.GetTH1D(Name+"_Min_Isolation_hadEt03","Min_Isolation_hadEt03",30,0,15,"Min HCAL deposits in Iso03 (GeV)","Events");
  Min_Isolation_hadVetoEt03=HConfig.GetTH1D(Name+"_Min_Isolation_hadVetoEt03","Min_Isolation_hadVetoEt03",50,0,10,"Min HCAL deposits in Iso03 Veto (GeV)","Events");
  Min_Isolation_nJets03=HConfig.GetTH1D(Name+"_Min_Isolation_nJets03","Min_Isolation_nJets03",5,0,5,"Min number of Jets in Iso03","Events");
  Min_Isolation_nTracks03=HConfig.GetTH1D(Name+"_Min_Isolation_nTracks03","Min_Isolation_nTracks03",20,0,20,"Min number of tracks in ECAL Iso03","Events");
  Min_Isolation_sumPt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumPt03","Min_Isolation_sumPt03",50,0,50,"Min #Sigmap_{T} in Iso03 (GeV)","Events");
  Min_Isolation_trackerVetoPt03=HConfig.GetTH1D(Name+"_Min_Isolation_trackerVetoPt03","Min_Isolation_trackerVetoPt03",10,0,10,"Min p_{T} in tracker Iso03 Veto (GeV)","Events");
  Min_Isolation_emEt05=HConfig.GetTH1D(Name+"_Min_Isolation_emEt05","Min_Isolation_emEt05",50,0,25,"Min ECAL deposits in Iso05 (GeV)","Events");
  Min_Isolation_emVetoEt05=HConfig.GetTH1D(Name+"_Min_Isolation_emVetoEt05","Min_Isolation_emVetoEt05",50,0,10,"Min ECAL deposits in Iso05 Veto (GeV)","Events");
  Min_Isolation_hadEt05=HConfig.GetTH1D(Name+"_Min_Isolation_hadEt05","Min_Isolation_hadEt05",50,0,25,"Min HCAL deposits in Iso05 Veto (GeV)","Events");
  Min_Isolation_hadVetoEt05=HConfig.GetTH1D(Name+"_Min_Isolation_hadVetoEt05","Min_Isolation_hadVetoEt05",50,0,10,"Min HCAL deposits in Iso05 Veto (GeV)","Events");
  Min_Isolation_nJets05=HConfig.GetTH1D(Name+"_Min_Isolation_nJets05","Min_Isolation_nJets05",5,0,5,"Min number of Jets in Iso05","Events");
  Min_Isolation_nTracks05=HConfig.GetTH1D(Name+"_Min_Isolation_nTracks05","Min_Isolation_nTracks05",20,0,20,"Min number of tracks in Iso05","Events");
  Min_Isolation_sumPt05=HConfig.GetTH1D(Name+"_Min_Isolation_sumPt05","Min_Isolation_sumPt05",50,0,50,"Min #Sigmap_{T} in Iso05 (GeV)","Events");
  Min_Isolation_trackerVetoPt05=HConfig.GetTH1D(Name+"_Min_Isolation_trackerVetoPt05","Min_Isolation_trackerVetoPt05",50,0,25,"Min p_{T} in the tracker Iso05 Veto (GeV)","Events");
  Min_Isolation_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumChargedHadronPt03","Min_Isolation_sumChargedHadronPt03",30,0,15,"Min #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Min_Isolation_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumChargedParticlePt03","Min_Isolation_sumChargedParticlePt03",30,0,15,"Min #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Min_Isolation_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumNeutralHadronEt03","Min_Isolation_sumNeutralHadronEt03",30,0,15,"Min #SigmaE_{T} in Iso03 (GeV)","Events");
  Min_Isolation_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Min_Isolation_sumNeutralHadronEtHighThreshold03","Min_Isolation_sumNeutralHadronEtHighThreshold03",30,0,15,"Min #SigmaE_{T} (Neutral Hadron) (HT) in Iso03 (GeV)","Events");
  Min_Isolation_sumPhotonEt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumPhotonEt03","Min_Isolation_sumPhotonEt03",30,0,15,"Min #SigmaE_{T} in Iso03 (GeV)","Events");
  Min_Isolation_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Min_Isolation_sumPhotonEtHighThreshold03","Min_Isolation_sumPhotonEtHighThreshold03",30,0,15,"Min #SigmaE_{T}#gamma (high threshold) in Iso03 (GeV)","Events");
  Min_Isolation_sumPUPt03=HConfig.GetTH1D(Name+"_Min_Isolation_sumPUPt03","Min_Isolation_sumPUPt03",30,0,15,"Min #SigmaPUp_{T} in Iso03 (GeV)","Events");
  Min_Isolation_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Min_Isolation_sumChargedHadronPt04","Min_Isolation_sumChargedHadronPt04",30,0,15,"Min #Sigmap_{T} (charged hadrons) in Iso04 (GeV)","Events");
  Min_Isolation_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Min_Isolation_sumChargedParticlePt04","Min_Isolation_sumChargedParticlePt04",30,0,15,"Min #Sigmap_{T} (charged particles) (GeV)","Events");
  Min_Isolation_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Min_Isolation_sumNeutralHadronEt04","Min_Isolation_sumNeutralHadronEt04",30,0,15,"Min #SigmaE_{T} (neutral hadrons) in Iso04 (GeV)","Events");
  Min_Isolation_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Min_Isolation_sumNeutralHadronEtHighThreshold04","Min_Isolation_sumNeutralHadronEtHighThreshold04",30,0,15,"Min #SigmaE_{T} (neutral hadrons) (HT) in Iso04 (GeV)","Events");
  Min_Isolation_sumPhotonEt04=HConfig.GetTH1D(Name+"_Min_Isolation_sumPhotonEt04","Min_Isolation_sumPhotonEt04",100,0,50,"Min #SigmaE_{T} #gamma in Iso04 (GeV)","Events");
  Min_Isolation_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Min_Isolation_sumPhotonEtHighThreshold04","Min_Isolation_sumPhotonEtHighThreshold04",30,0,15,"Min #SigmaE_{T}#gamma (high threshold) in Iso04 (GeV)","Events");
  Min_Isolation_sumPUPt04=HConfig.GetTH1D(Name+"_Min_Isolation_sumPUPt04","Min_Isolation_sumPUPt04",30,0,15,"Min #SigmaPUp_{T} in Iso04 (GeV)","Events");

  // ------------------------------------------------------------------------------------------------------------------------

  Max_Isolation_emEt03=HConfig.GetTH1D(Name+"_Max_Isolation_emEt03","Max_Isolation_emEt03",50,0,25,"Max ECAL deposits in Iso03 (GeV)","Events");
  Max_Isolation_emVetoEt03=HConfig.GetTH1D(Name+"_Max_Isolation_emVetoEt03","Max_Isolation_emVetoEt03",50,0,10,"Max ECAL deposits in Iso03 Veto (GeV)","Events");
  Max_Isolation_hadEt03=HConfig.GetTH1D(Name+"_Max_Isolation_hadEt03","Max_Isolation_hadEt03",30,0,15,"Max HCAL deposits in Iso03 (GeV)","Events");
  Max_Isolation_hadVetoEt03=HConfig.GetTH1D(Name+"_Max_Isolation_hadVetoEt03","Max_Isolation_hadVetoEt03",50,0,10,"Max HCAL deposits in Iso03 Veto (GeV)","Events");
  Max_Isolation_nJets03=HConfig.GetTH1D(Name+"_Max_Isolation_nJets03","Max_Isolation_nJets03",5,0,5,"Max number of Jets in Iso03","Events");
  Max_Isolation_nTracks03=HConfig.GetTH1D(Name+"_Max_Isolation_nTracks03","Max_Isolation_nTracks03",20,0,20,"Max number of tracks in ECAL Iso03","Events");
  Max_Isolation_sumPt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumPt03","Max_Isolation_sumPt03",50,0,50,"Max #Sigmap_{T} in Iso03 (GeV)","Events");
  Max_Isolation_trackerVetoPt03=HConfig.GetTH1D(Name+"_Max_Isolation_trackerVetoPt03","Max_Isolation_trackerVetoPt03",10,0,10,"Max p_{T} in tracker Iso03 Veto (GeV)","Events");
  Max_Isolation_emEt05=HConfig.GetTH1D(Name+"_Max_Isolation_emEt05","Max_Isolation_emEt05",50,0,25,"Max ECAL deposits in Iso05 (GeV)","Events");
  Max_Isolation_emVetoEt05=HConfig.GetTH1D(Name+"_Max_Isolation_emVetoEt05","Max_Isolation_emVetoEt05",50,0,10,"Max ECAL deposits in Iso05 Veto (GeV)","Events");
  Max_Isolation_hadEt05=HConfig.GetTH1D(Name+"_Max_Isolation_hadEt05","Max_Isolation_hadEt05",50,0,25,"Max HCAL deposits in Iso05 Veto (GeV)","Events");
  Max_Isolation_hadVetoEt05=HConfig.GetTH1D(Name+"_Max_Isolation_hadVetoEt05","Max_Isolation_hadVetoEt05",50,0,10,"Max HCAL deposits in Iso05 Veto (GeV)","Events");
  Max_Isolation_nJets05=HConfig.GetTH1D(Name+"_Max_Isolation_nJets05","Max_Isolation_nJets05",5,0,5,"Max number of Jets in Iso05","Events");
  Max_Isolation_nTracks05=HConfig.GetTH1D(Name+"_Max_Isolation_nTracks05","Max_Isolation_nTracks05",20,0,20,"Max number of tracks in Iso05","Events");
  Max_Isolation_sumPt05=HConfig.GetTH1D(Name+"_Max_Isolation_sumPt05","Max_Isolation_sumPt05",50,0,50,"Max #Sigmap_{T} in Iso05 (GeV)","Events");
  Max_Isolation_trackerVetoPt05=HConfig.GetTH1D(Name+"_Max_Isolation_trackerVetoPt05","Max_Isolation_trackerVetoPt05",50,0,25,"Max p_{T} in the tracker Iso05 Veto (GeV)","Events");
  Max_Isolation_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumChargedHadronPt03","Max_Isolation_sumChargedHadronPt03",30,0,15,"Max #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Max_Isolation_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumChargedParticlePt03","Max_Isolation_sumChargedParticlePt03",30,0,15,"Max #Sigmap_{T} (charged hadrons) in Iso03 (GeV)","Events");
  Max_Isolation_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumNeutralHadronEt03","Max_Isolation_sumNeutralHadronEt03",30,0,15,"Max #SigmaE_{T} in Iso03 (GeV)","Events");
  Max_Isolation_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Max_Isolation_sumNeutralHadronEtHighThreshold03","Max_Isolation_sumNeutralHadronEtHighThreshold03",30,0,15,"Max #SigmaE_{T} (Neutral Hadron) (HT) in Iso03 (GeV)","Events");
  Max_Isolation_sumPhotonEt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumPhotonEt03","Max_Isolation_sumPhotonEt03",30,0,15,"Max #SigmaE_{T} in Iso03 (GeV)","Events");
  Max_Isolation_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Max_Isolation_sumPhotonEtHighThreshold03","Max_Isolation_sumPhotonEtHighThreshold03",30,0,15,"Max #SigmaE_{T}#gamma (high threshold) in Iso03 (GeV)","Events");
  Max_Isolation_sumPUPt03=HConfig.GetTH1D(Name+"_Max_Isolation_sumPUPt03","Max_Isolation_sumPUPt03",30,0,15,"Max #SigmaPUp_{T} in Iso03 (GeV)","Events");
  Max_Isolation_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Max_Isolation_sumChargedHadronPt04","Max_Isolation_sumChargedHadronPt04",30,0,15,"Max #Sigmap_{T} (charged hadrons) in Iso04 (GeV)","Events");
  Max_Isolation_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Max_Isolation_sumChargedParticlePt04","Max_Isolation_sumChargedParticlePt04",30,0,15,"Max #Sigmap_{T} (charged particles) (GeV)","Events");
  Max_Isolation_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Max_Isolation_sumNeutralHadronEt04","Max_Isolation_sumNeutralHadronEt04",30,0,15,"Max #SigmaE_{T} (neutral hadrons) in Iso04 (GeV)","Events");
  Max_Isolation_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Max_Isolation_sumNeutralHadronEtHighThreshold04","Max_Isolation_sumNeutralHadronEtHighThreshold04",30,0,15,"Max #SigmaE_{T} (neutral hadrons) (HT) in Iso04 (GeV)","Events");
  Max_Isolation_sumPhotonEt04=HConfig.GetTH1D(Name+"_Max_Isolation_sumPhotonEt04","Max_Isolation_sumPhotonEt04",100,0,50,"Max #SigmaE_{T} #gamma in Iso04 (GeV)","Events");
  Max_Isolation_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Max_Isolation_sumPhotonEtHighThreshold04","Max_Isolation_sumPhotonEtHighThreshold04",30,0,15,"Max #SigmaE_{T}#gamma (high threshold) in Iso04 (GeV)","Events");
  Max_Isolation_sumPUPt04=HConfig.GetTH1D(Name+"_Max_Isolation_sumPUPt04","Max_Isolation_sumPUPt04",30,0,15,"Max #SigmaPUp_{T} in Iso04 (GeV)","Events");
   
	Muon1_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Muon1_calEnergy_hadS9","Muon1_calEnergy_hadS9",30,0,30,"#mu_{1} Energy deposit in crossed 3x3 HCAL crystal","Events");
   Muon1_calEnergy_had=HConfig.GetTH1D(Name+"_Muon1_calEnergy_had","Muon1_calEnergy_had",30,0,30,"#mu_{1} Energy deposit in crossed HCAL tower","Events");
   Muon1_calEnergy_emS25=HConfig.GetTH1D(Name+"_Muon1_calEnergy_emS25","Muon1_calEnergy_emS25",30,0,30,"#mu_{1} Energy deposit in 5x5 ECAL crystal","Events");
   Muon1_calEnergy_emS9=HConfig.GetTH1D(Name+"_Muon1_calEnergy_emS9","Muon1_calEnergy_emS9",30,0,30,"#mu_{1} Energy deposit in 3x3 ECAL crystal","Events");
   Muon1_calEnergy_em=HConfig.GetTH1D(Name+"_Muon1_calEnergy_em","Muon1_calEnergy_em",30,0,30,"Enenrgy deposit in crossed ECAL crystals","Events");

	Muon2_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Muon2_calEnergy_hadS9","Muon2_calEnergy_hadS9",30,0,30,"#mu_{2} Energy deposit in crossed 3x3 HCAL crystal","Events");
   Muon2_calEnergy_had=HConfig.GetTH1D(Name+"_Muon2_calEnergy_had","Muon2_calEnergy_had",30,0,30,"#mu_{2} Energy deposit in crossed HCAL tower","Events");
   Muon2_calEnergy_emS25=HConfig.GetTH1D(Name+"_Muon2_calEnergy_emS25","Muon2_calEnergy_emS25",30,0,30,"#mu_{2} Energy deposit in 5x5 ECAL crystal","Events");
   Muon2_calEnergy_emS9=HConfig.GetTH1D(Name+"_Muon2_calEnergy_emS9","Muon2_calEnergy_emS9",30,0,30,"#mu_{2} Energy deposit in 3x3 ECAL crystal","Events");
   Muon2_calEnergy_em=HConfig.GetTH1D(Name+"_Muon2_calEnergy_em","Muon2_calEnergy_em",30,0,30,"Enenrgy deposit in crossed ECAL crystals","Events");

	Muon3_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Muon3_calEnergy_hadS9","Muon3_calEnergy_hadS9",30,0,30,"#mu_{3} Energy deposit in crossed 3x3 HCAL crystal","Events");
   Muon3_calEnergy_had=HConfig.GetTH1D(Name+"_Muon3_calEnergy_had","Muon3_calEnergy_had",30,0,30,"#mu_{3} Energy deposit in crossed HCAL tower","Events");
   Muon3_calEnergy_emS25=HConfig.GetTH1D(Name+"_Muon3_calEnergy_emS25","Muon3_calEnergy_emS25",30,0,30,"#mu_{3} Energy deposit in 5x5 ECAL crystal","Events");
   Muon3_calEnergy_emS9=HConfig.GetTH1D(Name+"_Muon3_calEnergy_emS9","Muon3_calEnergy_emS9",30,0,30,"#mu_{3} Energy deposit in 3x3 ECAL crystal","Events");
   Muon3_calEnergy_em=HConfig.GetTH1D(Name+"_Muon3_calEnergy_em","Muon3_calEnergy_em",30,0,30,"Enenrgy deposit in crossed ECAL crystals","Events");

	Min_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Min_calEnergy_hadS9","Min_calEnergy_hadS9",30,0,30,"Min Energy deposit in crossed 3x3 HCAL crystal","Events");
   Min_calEnergy_had=HConfig.GetTH1D(Name+"_Min_calEnergy_had","Min_calEnergy_had",30,0,30,"Min Energy deposit in crossed HCAL tower","Events");
   Min_calEnergy_emS25=HConfig.GetTH1D(Name+"_Min_calEnergy_emS25","Min_calEnergy_emS25",30,0,30,"Min Energy deposit in 5x5 ECAL crystal","Events");
   Min_calEnergy_emS9=HConfig.GetTH1D(Name+"_Min_calEnergy_emS9","Min_calEnergy_emS9",30,0,30,"Min Energy deposit in 3x3 ECAL crystal","Events");
   Min_calEnergy_em=HConfig.GetTH1D(Name+"_Min_calEnergy_em","Min_calEnergy_em",30,0,30,"Enenrgy deposit in crossed ECAL crystals","Events");

	Max_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Max_calEnergy_hadS9","Max_calEnergy_hadS9",30,0,30,"Max Energy deposit in crossed 3x3 HCAL crystal","Events");
   Max_calEnergy_had=HConfig.GetTH1D(Name+"_Max_calEnergy_had","Max_calEnergy_had",30,0,30,"Max Energy deposit in crossed HCAL tower","Events");
   Max_calEnergy_emS25=HConfig.GetTH1D(Name+"_Max_calEnergy_emS25","Max_calEnergy_emS25",30,0,30,"Max Energy deposit in 5x5 ECAL crystal","Events");
   Max_calEnergy_emS9=HConfig.GetTH1D(Name+"_Max_calEnergy_emS9","Max_calEnergy_emS9",30,0,30,"Max Energy deposit in 3x3 ECAL crystal","Events");
   Max_calEnergy_em=HConfig.GetTH1D(Name+"_Max_calEnergy_em","Max_calEnergy_em",30,0,30,"Enenrgy deposit in crossed ECAL crystals","Events");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  TripleMuIsolation::Store_ExtraDist(){ 

  // Muon isolation variables
  Extradist1d.push_back(&Muon1_Isolation_emEt03);
  Extradist1d.push_back(&Muon1_Isolation_emVetoEt03);
  Extradist1d.push_back(&Muon1_Isolation_hadEt03);
  Extradist1d.push_back(&Muon1_Isolation_hadVetoEt03);
  Extradist1d.push_back(&Muon1_Isolation_nJets03);
  Extradist1d.push_back(&Muon1_Isolation_nTracks03);
  Extradist1d.push_back(&Muon1_Isolation_sumPt03);
  Extradist1d.push_back(&Muon1_Isolation_trackerVetoPt03);
  Extradist1d.push_back(&Muon1_Isolation_emEt05);
  Extradist1d.push_back(&Muon1_Isolation_emVetoEt05);
  Extradist1d.push_back(&Muon1_Isolation_hadEt05);
  Extradist1d.push_back(&Muon1_Isolation_hadVetoEt05);
  Extradist1d.push_back(&Muon1_Isolation_nJets05);
  Extradist1d.push_back(&Muon1_Isolation_nTracks05);
  Extradist1d.push_back(&Muon1_Isolation_sumPt05);
  Extradist1d.push_back(&Muon1_Isolation_trackerVetoPt05);
  Extradist1d.push_back(&Muon1_Isolation_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon1_Isolation_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon1_Isolation_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon1_Isolation_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon1_Isolation_sumPhotonEt03);
  Extradist1d.push_back(&Muon1_Isolation_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon1_Isolation_sumPUPt03);
  Extradist1d.push_back(&Muon1_Isolation_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon1_Isolation_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon1_Isolation_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon1_Isolation_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon1_Isolation_sumPhotonEt04);
  Extradist1d.push_back(&Muon1_Isolation_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon1_Isolation_sumPUPt04);

  Extradist1d.push_back(&Muon2_Isolation_emEt03);
  Extradist1d.push_back(&Muon2_Isolation_emVetoEt03);
  Extradist1d.push_back(&Muon2_Isolation_hadEt03);
  Extradist1d.push_back(&Muon2_Isolation_hadVetoEt03);
  Extradist1d.push_back(&Muon2_Isolation_nJets03);
  Extradist1d.push_back(&Muon2_Isolation_nTracks03);
  Extradist1d.push_back(&Muon2_Isolation_sumPt03);
  Extradist1d.push_back(&Muon2_Isolation_trackerVetoPt03);
  Extradist1d.push_back(&Muon2_Isolation_emEt05);
  Extradist1d.push_back(&Muon2_Isolation_emVetoEt05);
  Extradist1d.push_back(&Muon2_Isolation_hadEt05);
  Extradist1d.push_back(&Muon2_Isolation_hadVetoEt05);
  Extradist1d.push_back(&Muon2_Isolation_nJets05);
  Extradist1d.push_back(&Muon2_Isolation_nTracks05);
  Extradist1d.push_back(&Muon2_Isolation_sumPt05);
  Extradist1d.push_back(&Muon2_Isolation_trackerVetoPt05);
  Extradist1d.push_back(&Muon2_Isolation_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon2_Isolation_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon2_Isolation_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon2_Isolation_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon2_Isolation_sumPhotonEt03);
  Extradist1d.push_back(&Muon2_Isolation_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon2_Isolation_sumPUPt03);
  Extradist1d.push_back(&Muon2_Isolation_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon2_Isolation_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon2_Isolation_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon2_Isolation_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon2_Isolation_sumPhotonEt04);
  Extradist1d.push_back(&Muon2_Isolation_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon2_Isolation_sumPUPt04);

  Extradist1d.push_back(&Muon3_Isolation_emEt03);
  Extradist1d.push_back(&Muon3_Isolation_emVetoEt03);
  Extradist1d.push_back(&Muon3_Isolation_hadEt03);
  Extradist1d.push_back(&Muon3_Isolation_hadVetoEt03);
  Extradist1d.push_back(&Muon3_Isolation_nJets03);
  Extradist1d.push_back(&Muon3_Isolation_nTracks03);
  Extradist1d.push_back(&Muon3_Isolation_sumPt03);
  Extradist1d.push_back(&Muon3_Isolation_trackerVetoPt03);
  Extradist1d.push_back(&Muon3_Isolation_emEt05);
  Extradist1d.push_back(&Muon3_Isolation_emVetoEt05);
  Extradist1d.push_back(&Muon3_Isolation_hadEt05);
  Extradist1d.push_back(&Muon3_Isolation_hadVetoEt05);
  Extradist1d.push_back(&Muon3_Isolation_nJets05);
  Extradist1d.push_back(&Muon3_Isolation_nTracks05);
  Extradist1d.push_back(&Muon3_Isolation_sumPt05);
  Extradist1d.push_back(&Muon3_Isolation_trackerVetoPt05);
  Extradist1d.push_back(&Muon3_Isolation_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon3_Isolation_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon3_Isolation_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon3_Isolation_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon3_Isolation_sumPhotonEt03);
  Extradist1d.push_back(&Muon3_Isolation_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon3_Isolation_sumPUPt03);
  Extradist1d.push_back(&Muon3_Isolation_sumChargedHadronPt04);
  Extradist1d.push_back(&Muon3_Isolation_sumChargedParticlePt04);
  Extradist1d.push_back(&Muon3_Isolation_sumNeutralHadronEt04);
  Extradist1d.push_back(&Muon3_Isolation_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Muon3_Isolation_sumPhotonEt04);
  Extradist1d.push_back(&Muon3_Isolation_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Muon3_Isolation_sumPUPt04);

  Extradist1d.push_back(&Min_Isolation_emEt03);
  Extradist1d.push_back(&Min_Isolation_emVetoEt03);
  Extradist1d.push_back(&Min_Isolation_hadEt03);
  Extradist1d.push_back(&Min_Isolation_hadVetoEt03);
  Extradist1d.push_back(&Min_Isolation_nJets03);
  Extradist1d.push_back(&Min_Isolation_nTracks03);
  Extradist1d.push_back(&Min_Isolation_sumPt03);
  Extradist1d.push_back(&Min_Isolation_trackerVetoPt03);
  Extradist1d.push_back(&Min_Isolation_emEt05);
  Extradist1d.push_back(&Min_Isolation_emVetoEt05);
  Extradist1d.push_back(&Min_Isolation_hadEt05);
  Extradist1d.push_back(&Min_Isolation_hadVetoEt05);
  Extradist1d.push_back(&Min_Isolation_nJets05);
  Extradist1d.push_back(&Min_Isolation_nTracks05);
  Extradist1d.push_back(&Min_Isolation_sumPt05);
  Extradist1d.push_back(&Min_Isolation_trackerVetoPt05);
  Extradist1d.push_back(&Min_Isolation_sumChargedHadronPt03);
  Extradist1d.push_back(&Min_Isolation_sumChargedParticlePt03);
  Extradist1d.push_back(&Min_Isolation_sumNeutralHadronEt03);
  Extradist1d.push_back(&Min_Isolation_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Min_Isolation_sumPhotonEt03);
  Extradist1d.push_back(&Min_Isolation_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Min_Isolation_sumPUPt03);
  Extradist1d.push_back(&Min_Isolation_sumChargedHadronPt04);
  Extradist1d.push_back(&Min_Isolation_sumChargedParticlePt04);
  Extradist1d.push_back(&Min_Isolation_sumNeutralHadronEt04);
  Extradist1d.push_back(&Min_Isolation_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Min_Isolation_sumPhotonEt04);
  Extradist1d.push_back(&Min_Isolation_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Min_Isolation_sumPUPt04);

  Extradist1d.push_back(&Max_Isolation_emEt03);
  Extradist1d.push_back(&Max_Isolation_emVetoEt03);
  Extradist1d.push_back(&Max_Isolation_hadEt03);
  Extradist1d.push_back(&Max_Isolation_hadVetoEt03);
  Extradist1d.push_back(&Max_Isolation_nJets03);
  Extradist1d.push_back(&Max_Isolation_nTracks03);
  Extradist1d.push_back(&Max_Isolation_sumPt03);
  Extradist1d.push_back(&Max_Isolation_trackerVetoPt03);
  Extradist1d.push_back(&Max_Isolation_emEt05);
  Extradist1d.push_back(&Max_Isolation_emVetoEt05);
  Extradist1d.push_back(&Max_Isolation_hadEt05);
  Extradist1d.push_back(&Max_Isolation_hadVetoEt05);
  Extradist1d.push_back(&Max_Isolation_nJets05);
  Extradist1d.push_back(&Max_Isolation_nTracks05);
  Extradist1d.push_back(&Max_Isolation_sumPt05);
  Extradist1d.push_back(&Max_Isolation_trackerVetoPt05);
  Extradist1d.push_back(&Max_Isolation_sumChargedHadronPt03);
  Extradist1d.push_back(&Max_Isolation_sumChargedParticlePt03);
  Extradist1d.push_back(&Max_Isolation_sumNeutralHadronEt03);
  Extradist1d.push_back(&Max_Isolation_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Max_Isolation_sumPhotonEt03);
  Extradist1d.push_back(&Max_Isolation_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Max_Isolation_sumPUPt03);
  Extradist1d.push_back(&Max_Isolation_sumChargedHadronPt04);
  Extradist1d.push_back(&Max_Isolation_sumChargedParticlePt04);
  Extradist1d.push_back(&Max_Isolation_sumNeutralHadronEt04);
  Extradist1d.push_back(&Max_Isolation_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Max_Isolation_sumPhotonEt04);
  Extradist1d.push_back(&Max_Isolation_sumPhotonEtHighThreshold04);
  Extradist1d.push_back(&Max_Isolation_sumPUPt04);
  // ------------------------------------------------------------------------------------------------------------------------

  Extradist1d.push_back(&Muon1_calEnergy_hadS9);
  Extradist1d.push_back(&Muon1_calEnergy_had);
  Extradist1d.push_back(&Muon1_calEnergy_emS25);
  Extradist1d.push_back(&Muon1_calEnergy_emS9);
  Extradist1d.push_back(&Muon1_calEnergy_em);

  Extradist1d.push_back(&Muon2_calEnergy_hadS9);
  Extradist1d.push_back(&Muon2_calEnergy_had);
  Extradist1d.push_back(&Muon2_calEnergy_emS25);
  Extradist1d.push_back(&Muon2_calEnergy_emS9);
  Extradist1d.push_back(&Muon2_calEnergy_em);

  Extradist1d.push_back(&Muon3_calEnergy_hadS9);
  Extradist1d.push_back(&Muon3_calEnergy_had);
  Extradist1d.push_back(&Muon3_calEnergy_emS25);
  Extradist1d.push_back(&Muon3_calEnergy_emS9);
  Extradist1d.push_back(&Muon3_calEnergy_em);
 
  Extradist1d.push_back(&Min_calEnergy_hadS9);
  Extradist1d.push_back(&Min_calEnergy_had);
  Extradist1d.push_back(&Min_calEnergy_emS25);
  Extradist1d.push_back(&Min_calEnergy_emS9);
  Extradist1d.push_back(&Min_calEnergy_em);

  Extradist1d.push_back(&Max_calEnergy_hadS9);
  Extradist1d.push_back(&Max_calEnergy_had);
  Extradist1d.push_back(&Max_calEnergy_emS25);
  Extradist1d.push_back(&Max_calEnergy_emS9);
  Extradist1d.push_back(&Max_calEnergy_em);
}


void  TripleMuIsolation::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu"))) value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);
  }

  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));
  unsigned int  final_idx=0;
  
  value.at(TriggerMatch)=0;
  value.at(SignalCandidate)=0;
  value.at(Mu1PtCut)=0;
  value.at(Mu2PtCut)=0;
  value.at(Mu3PtCut)=0;
  value.at(PhiVeto)=0;
  value.at(OmegaVeto)=0;
  value.at(TriggerMatch)=0;
  value.at(MuonID)=0;
  value.at(ThreeMuMass)=0;

  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();
    unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
    //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
    //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
    //
    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
    			Ntp->Muon_isTrackerMuon(mu3_pt_idx));
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

    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

    value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(final_idx)){
      value.at(TriggerMatch)+=i; 
    }
    
    value.at(ThreeMuMass) = TauLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) == cut.at(SignalCandidate));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
  pass.at(MuonID) =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = (value.at(TriggerMatch)  <  cut.at(TriggerMatch));
  pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*PDG_Var::Phi_width());
  pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*PDG_Var::Omega_width());

  if(id!=1) pass.at(ThreeMuMass) = true;
  else  pass.at(ThreeMuMass) = ( (value.at(ThreeMuMass) > tauMinSideBand_ && value.at(ThreeMuMass) < tauMinMass_)  ||   (value.at(ThreeMuMass)> tauMaxMass_ && value.at(ThreeMuMass) < tauMaxSideBand_));


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);

  if(status){

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
    
    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);


    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    
 	 // Muon isolation variables
    Muon1_Isolation_emEt03.at(t).Fill(Ntp->Muon_emEt03(Muon_index_1),w);
    Muon1_Isolation_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(Muon_index_1),w);
    Muon1_Isolation_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(Muon_index_1),w);
    Muon1_Isolation_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(Muon_index_1),w);
    Muon1_Isolation_nJets03.at(t).Fill(Ntp->Muon_nJets03(Muon_index_1),w);
    Muon1_Isolation_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(Muon_index_1),w);
    Muon1_Isolation_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(Muon_index_1),w);
    Muon1_Isolation_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(Muon_index_1),w);
    Muon1_Isolation_emEt05.at(t).Fill(Ntp->Muon_emEt05(Muon_index_1),w);
    Muon1_Isolation_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(Muon_index_1),w);
    Muon1_Isolation_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(Muon_index_1),w);
    Muon1_Isolation_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(Muon_index_1),w);
    Muon1_Isolation_nJets05.at(t).Fill(Ntp->Muon_nJets05(Muon_index_1),w);
    Muon1_Isolation_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(Muon_index_1),w);
    Muon1_Isolation_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(Muon_index_1),w);
    Muon1_Isolation_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(Muon_index_1),w);
    Muon1_Isolation_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(Muon_index_1),w);
    Muon1_Isolation_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(Muon_index_1),w);
    Muon1_Isolation_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(Muon_index_1),w);
    Muon1_Isolation_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_1),w);
    Muon1_Isolation_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(Muon_index_1),w);
    Muon1_Isolation_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_1),w);
    Muon1_Isolation_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(Muon_index_1),w);
    Muon1_Isolation_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(Muon_index_1),w);
    Muon1_Isolation_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(Muon_index_1),w);
    Muon1_Isolation_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(Muon_index_1),w);
    Muon1_Isolation_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_1),w);
    Muon1_Isolation_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(Muon_index_1),w);
    Muon1_Isolation_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_1),w);
    Muon1_Isolation_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(Muon_index_1),w);

    Muon2_Isolation_emEt03.at(t).Fill(Ntp->Muon_emEt03(Muon_index_2),w);
    Muon2_Isolation_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(Muon_index_2),w);
    Muon2_Isolation_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(Muon_index_2),w);
    Muon2_Isolation_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(Muon_index_2),w);
    Muon2_Isolation_nJets03.at(t).Fill(Ntp->Muon_nJets03(Muon_index_2),w);
    Muon2_Isolation_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(Muon_index_2),w);
    Muon2_Isolation_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(Muon_index_2),w);
    Muon2_Isolation_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(Muon_index_2),w);
    Muon2_Isolation_emEt05.at(t).Fill(Ntp->Muon_emEt05(Muon_index_2),w);
    Muon2_Isolation_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(Muon_index_2),w);
    Muon2_Isolation_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(Muon_index_2),w);
    Muon2_Isolation_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(Muon_index_2),w);
    Muon2_Isolation_nJets05.at(t).Fill(Ntp->Muon_nJets05(Muon_index_2),w);
    Muon2_Isolation_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(Muon_index_2),w);
    Muon2_Isolation_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(Muon_index_2),w);
    Muon2_Isolation_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(Muon_index_2),w);
    Muon2_Isolation_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(Muon_index_2),w);
    Muon2_Isolation_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(Muon_index_2),w);
    Muon2_Isolation_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(Muon_index_2),w);
    Muon2_Isolation_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_2),w);
    Muon2_Isolation_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(Muon_index_2),w);
    Muon2_Isolation_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_2),w);
    Muon2_Isolation_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(Muon_index_2),w);
    Muon2_Isolation_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(Muon_index_2),w);
    Muon2_Isolation_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(Muon_index_2),w);
    Muon2_Isolation_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(Muon_index_2),w);
    Muon2_Isolation_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_2),w);
    Muon2_Isolation_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(Muon_index_2),w);
    Muon2_Isolation_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_2),w);
    Muon2_Isolation_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(Muon_index_2),w);

    Muon3_Isolation_emEt03.at(t).Fill(Ntp->Muon_emEt03(Muon_index_3),w);
    Muon3_Isolation_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(Muon_index_3),w);
    Muon3_Isolation_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(Muon_index_3),w);
    Muon3_Isolation_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(Muon_index_3),w);
    Muon3_Isolation_nJets03.at(t).Fill(Ntp->Muon_nJets03(Muon_index_3),w);
    Muon3_Isolation_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(Muon_index_3),w);
    Muon3_Isolation_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(Muon_index_3),w);
    Muon3_Isolation_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(Muon_index_3),w);
    Muon3_Isolation_emEt05.at(t).Fill(Ntp->Muon_emEt05(Muon_index_3),w);
    Muon3_Isolation_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(Muon_index_3),w);
    Muon3_Isolation_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(Muon_index_3),w);
    Muon3_Isolation_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(Muon_index_3),w);
    Muon3_Isolation_nJets05.at(t).Fill(Ntp->Muon_nJets05(Muon_index_3),w);
    Muon3_Isolation_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(Muon_index_3),w);
    Muon3_Isolation_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(Muon_index_3),w);
    Muon3_Isolation_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(Muon_index_3),w);
    Muon3_Isolation_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(Muon_index_3),w);
    Muon3_Isolation_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(Muon_index_3),w);
    Muon3_Isolation_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(Muon_index_3),w);
    Muon3_Isolation_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_3),w);
    Muon3_Isolation_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(Muon_index_3),w);
    Muon3_Isolation_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_3),w);
    Muon3_Isolation_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(Muon_index_3),w);
    Muon3_Isolation_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(Muon_index_3),w);
    Muon3_Isolation_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(Muon_index_3),w);
    Muon3_Isolation_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(Muon_index_3),w);
    Muon3_Isolation_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_3),w);
    Muon3_Isolation_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(Muon_index_3),w);
    Muon3_Isolation_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_3),w);
    Muon3_Isolation_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(Muon_index_3),w);

    Max_Isolation_emEt03.at(t).Fill(std::max(Ntp->Muon_emEt03(Muon_index_1), std::max(Ntp->Muon_emEt03(Muon_index_2),Ntp->Muon_emEt03(Muon_index_3))));
    Max_Isolation_emVetoEt03.at(t).Fill(std::max(Ntp->Muon_emVetoEt03(Muon_index_1), std::max(Ntp->Muon_emVetoEt03(Muon_index_2),Ntp->Muon_emVetoEt03(Muon_index_3))));
    Max_Isolation_hadEt03.at(t).Fill(std::max(Ntp->Muon_hadEt03(Muon_index_1), std::max(Ntp->Muon_hadEt03(Muon_index_2),Ntp->Muon_hadEt03(Muon_index_3))));
    Max_Isolation_hadVetoEt03.at(t).Fill(std::max(Ntp->Muon_hadVetoEt03(Muon_index_1), std::max(Ntp->Muon_hadVetoEt03(Muon_index_2),Ntp->Muon_hadVetoEt03(Muon_index_3))));
    Max_Isolation_nJets03.at(t).Fill(std::max(Ntp->Muon_nJets03(Muon_index_1), std::max(Ntp->Muon_nJets03(Muon_index_2),Ntp->Muon_nJets03(Muon_index_3))));
    Max_Isolation_nTracks03.at(t).Fill(std::max(Ntp->Muon_nTracks03(Muon_index_1), std::max(Ntp->Muon_nTracks03(Muon_index_2),Ntp->Muon_nTracks03(Muon_index_3))));
    Max_Isolation_sumPt03.at(t).Fill(std::max(Ntp->Muon_sumPt03(Muon_index_1), std::max(Ntp->Muon_sumPt03(Muon_index_2),Ntp->Muon_sumPt03(Muon_index_3))));
    Max_Isolation_trackerVetoPt03.at(t).Fill(std::max(Ntp->Muon_trackerVetoPt03(Muon_index_1), std::max(Ntp->Muon_trackerVetoPt03(Muon_index_2),Ntp->Muon_trackerVetoPt03(Muon_index_3))));
    Max_Isolation_emEt05.at(t).Fill(std::max(Ntp->Muon_emEt05(Muon_index_1), std::max(Ntp->Muon_emEt05(Muon_index_2),Ntp->Muon_emEt05(Muon_index_3))));
    Max_Isolation_emVetoEt05.at(t).Fill(std::max(Ntp->Muon_emVetoEt05(Muon_index_1), std::max(Ntp->Muon_emVetoEt05(Muon_index_2),Ntp->Muon_emVetoEt05(Muon_index_3))));
    Max_Isolation_hadEt05.at(t).Fill(std::max(Ntp->Muon_hadEt05(Muon_index_1), std::max(Ntp->Muon_hadEt05(Muon_index_2),Ntp->Muon_hadEt05(Muon_index_3))));
    Max_Isolation_hadVetoEt05.at(t).Fill(std::max(Ntp->Muon_hadVetoEt05(Muon_index_1), std::max(Ntp->Muon_hadVetoEt05(Muon_index_2),Ntp->Muon_hadVetoEt05(Muon_index_3))));
    Max_Isolation_nJets05.at(t).Fill(std::max(Ntp->Muon_nJets05(Muon_index_1), std::max(Ntp->Muon_nJets05(Muon_index_2),Ntp->Muon_nJets05(Muon_index_3))));
    Max_Isolation_nTracks05.at(t).Fill(std::max(Ntp->Muon_nTracks05(Muon_index_1), std::max(Ntp->Muon_nTracks05(Muon_index_2),Ntp->Muon_nTracks05(Muon_index_3))));
    Max_Isolation_sumPt05.at(t).Fill(std::max(Ntp->Muon_sumPt05(Muon_index_1), std::max(Ntp->Muon_sumPt05(Muon_index_2),Ntp->Muon_sumPt05(Muon_index_3))));
    Max_Isolation_trackerVetoPt05.at(t).Fill(std::max(Ntp->Muon_trackerVetoPt05(Muon_index_1), std::max(Ntp->Muon_trackerVetoPt05(Muon_index_2),Ntp->Muon_trackerVetoPt05(Muon_index_3))));
    Max_Isolation_sumChargedHadronPt03.at(t).Fill(std::max(Ntp->Muon_sumChargedHadronPt03(Muon_index_1), std::max(Ntp->Muon_sumChargedHadronPt03(Muon_index_2),Ntp->Muon_sumChargedHadronPt03(Muon_index_3))));
    Max_Isolation_sumChargedParticlePt03.at(t).Fill(std::max(Ntp->Muon_sumChargedParticlePt03(Muon_index_1), std::max(Ntp->Muon_sumChargedParticlePt03(Muon_index_2),Ntp->Muon_sumChargedParticlePt03(Muon_index_3))));
    Max_Isolation_sumNeutralHadronEt03.at(t).Fill(std::max(Ntp->Muon_sumNeutralHadronEt03(Muon_index_1), std::max(Ntp->Muon_sumNeutralHadronEt03(Muon_index_2),Ntp->Muon_sumNeutralHadronEt03(Muon_index_3))));
    Max_Isolation_sumNeutralHadronEtHighThreshold03.at(t).Fill(std::max(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_1), std::max(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_2),Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_3))));
    Max_Isolation_sumPhotonEt03.at(t).Fill(std::max(Ntp->Muon_sumPhotonEt03(Muon_index_1), std::max(Ntp->Muon_sumPhotonEt03(Muon_index_2),Ntp->Muon_sumPhotonEt03(Muon_index_3))));
    Max_Isolation_sumPhotonEtHighThreshold03.at(t).Fill(std::max(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_1), std::max(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_2),Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_3))));
    Max_Isolation_sumPUPt03.at(t).Fill(std::max(Ntp->Muon_sumPUPt03(Muon_index_1), std::max(Ntp->Muon_sumPUPt03(Muon_index_2),Ntp->Muon_sumPUPt03(Muon_index_3))));
    Max_Isolation_sumChargedHadronPt04.at(t).Fill(std::max(Ntp->Muon_sumChargedHadronPt04(Muon_index_1), std::max(Ntp->Muon_sumChargedHadronPt04(Muon_index_2),Ntp->Muon_sumChargedHadronPt04(Muon_index_3))));
    Max_Isolation_sumChargedParticlePt04.at(t).Fill(std::max(Ntp->Muon_sumChargedParticlePt04(Muon_index_1), std::max(Ntp->Muon_sumChargedParticlePt04(Muon_index_2),Ntp->Muon_sumChargedParticlePt04(Muon_index_3))));
    Max_Isolation_sumNeutralHadronEt04.at(t).Fill(std::max(Ntp->Muon_sumNeutralHadronEt04(Muon_index_1), std::max(Ntp->Muon_sumNeutralHadronEt04(Muon_index_2),Ntp->Muon_sumNeutralHadronEt04(Muon_index_3))));
    Max_Isolation_sumNeutralHadronEtHighThreshold04.at(t).Fill(std::max(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_1), std::max(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_2),Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_3))));
    Max_Isolation_sumPhotonEt04.at(t).Fill(std::max(Ntp->Muon_sumPhotonEt04(Muon_index_1), std::max(Ntp->Muon_sumPhotonEt04(Muon_index_2),Ntp->Muon_sumPhotonEt04(Muon_index_3))));
    Max_Isolation_sumPhotonEtHighThreshold04.at(t).Fill(std::max(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_1), std::max(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_2),Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_3))));
    Max_Isolation_sumPUPt04.at(t).Fill(std::max(Ntp->Muon_sumPUPt04(Muon_index_1), std::max(Ntp->Muon_sumPUPt04(Muon_index_2),Ntp->Muon_sumPUPt04(Muon_index_3))));

    Min_Isolation_emEt03.at(t).Fill(std::min(Ntp->Muon_emEt03(Muon_index_1), std::min(Ntp->Muon_emEt03(Muon_index_2),Ntp->Muon_emEt03(Muon_index_3))));
 Min_Isolation_emVetoEt03.at(t).Fill(std::min(Ntp->Muon_emVetoEt03(Muon_index_1), std::min(Ntp->Muon_emVetoEt03(Muon_index_2),Ntp->Muon_emVetoEt03(Muon_index_3))));
    Min_Isolation_hadEt03.at(t).Fill(std::min(Ntp->Muon_hadEt03(Muon_index_1), std::min(Ntp->Muon_hadEt03(Muon_index_2),Ntp->Muon_hadEt03(Muon_index_3))));
    Min_Isolation_hadVetoEt03.at(t).Fill(std::min(Ntp->Muon_hadVetoEt03(Muon_index_1), std::min(Ntp->Muon_hadVetoEt03(Muon_index_2),Ntp->Muon_hadVetoEt03(Muon_index_3))));
    Min_Isolation_nJets03.at(t).Fill(std::min(Ntp->Muon_nJets03(Muon_index_1), std::min(Ntp->Muon_nJets03(Muon_index_2),Ntp->Muon_nJets03(Muon_index_3))));
    Min_Isolation_nTracks03.at(t).Fill(std::min(Ntp->Muon_nTracks03(Muon_index_1), std::min(Ntp->Muon_nTracks03(Muon_index_2),Ntp->Muon_nTracks03(Muon_index_3))));
    Min_Isolation_sumPt03.at(t).Fill(std::min(Ntp->Muon_sumPt03(Muon_index_1), std::min(Ntp->Muon_sumPt03(Muon_index_2),Ntp->Muon_sumPt03(Muon_index_3))));
    Min_Isolation_trackerVetoPt03.at(t).Fill(std::min(Ntp->Muon_trackerVetoPt03(Muon_index_1), std::min(Ntp->Muon_trackerVetoPt03(Muon_index_2),Ntp->Muon_trackerVetoPt03(Muon_index_3))));
    Min_Isolation_emEt05.at(t).Fill(std::min(Ntp->Muon_emEt05(Muon_index_1), std::min(Ntp->Muon_emEt05(Muon_index_2),Ntp->Muon_emEt05(Muon_index_3))));
    Min_Isolation_emVetoEt05.at(t).Fill(std::min(Ntp->Muon_emVetoEt05(Muon_index_1), std::min(Ntp->Muon_emVetoEt05(Muon_index_2),Ntp->Muon_emVetoEt05(Muon_index_3))));
    Min_Isolation_hadEt05.at(t).Fill(std::min(Ntp->Muon_hadEt05(Muon_index_1), std::min(Ntp->Muon_hadEt05(Muon_index_2),Ntp->Muon_hadEt05(Muon_index_3))));
    Min_Isolation_hadVetoEt05.at(t).Fill(std::min(Ntp->Muon_hadVetoEt05(Muon_index_1), std::min(Ntp->Muon_hadVetoEt05(Muon_index_2),Ntp->Muon_hadVetoEt05(Muon_index_3))));
    Min_Isolation_nJets05.at(t).Fill(std::min(Ntp->Muon_nJets05(Muon_index_1), std::min(Ntp->Muon_nJets05(Muon_index_2),Ntp->Muon_nJets05(Muon_index_3))));
    Min_Isolation_nTracks05.at(t).Fill(std::min(Ntp->Muon_nTracks05(Muon_index_1), std::min(Ntp->Muon_nTracks05(Muon_index_2),Ntp->Muon_nTracks05(Muon_index_3))));
    Min_Isolation_sumPt05.at(t).Fill(std::min(Ntp->Muon_sumPt05(Muon_index_1), std::min(Ntp->Muon_sumPt05(Muon_index_2),Ntp->Muon_sumPt05(Muon_index_3))));
    Min_Isolation_trackerVetoPt05.at(t).Fill(std::min(Ntp->Muon_trackerVetoPt05(Muon_index_1), std::min(Ntp->Muon_trackerVetoPt05(Muon_index_2),Ntp->Muon_trackerVetoPt05(Muon_index_3))));
    Min_Isolation_sumChargedHadronPt03.at(t).Fill(std::min(Ntp->Muon_sumChargedHadronPt03(Muon_index_1), std::min(Ntp->Muon_sumChargedHadronPt03(Muon_index_2),Ntp->Muon_sumChargedHadronPt03(Muon_index_3))));
    Min_Isolation_sumChargedParticlePt03.at(t).Fill(std::min(Ntp->Muon_sumChargedParticlePt03(Muon_index_1), std::min(Ntp->Muon_sumChargedParticlePt03(Muon_index_2),Ntp->Muon_sumChargedParticlePt03(Muon_index_3))));
    Min_Isolation_sumNeutralHadronEt03.at(t).Fill(std::min(Ntp->Muon_sumNeutralHadronEt03(Muon_index_1), std::min(Ntp->Muon_sumNeutralHadronEt03(Muon_index_2),Ntp->Muon_sumNeutralHadronEt03(Muon_index_3))));
    Min_Isolation_sumNeutralHadronEtHighThreshold03.at(t).Fill(std::min(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_1), std::min(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_2),Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_3))));
    Min_Isolation_sumPhotonEt03.at(t).Fill(std::min(Ntp->Muon_sumPhotonEt03(Muon_index_1), std::min(Ntp->Muon_sumPhotonEt03(Muon_index_2),Ntp->Muon_sumPhotonEt03(Muon_index_3))));
    Min_Isolation_sumPhotonEtHighThreshold03.at(t).Fill(std::min(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_1), std::min(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_2),Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_3))));
    Min_Isolation_sumPUPt03.at(t).Fill(std::min(Ntp->Muon_sumPUPt03(Muon_index_1), std::min(Ntp->Muon_sumPUPt03(Muon_index_2),Ntp->Muon_sumPUPt03(Muon_index_3))));
    Min_Isolation_sumChargedHadronPt04.at(t).Fill(std::min(Ntp->Muon_sumChargedHadronPt04(Muon_index_1), std::min(Ntp->Muon_sumChargedHadronPt04(Muon_index_2),Ntp->Muon_sumChargedHadronPt04(Muon_index_3))));
    Min_Isolation_sumChargedParticlePt04.at(t).Fill(std::min(Ntp->Muon_sumChargedParticlePt04(Muon_index_1), std::min(Ntp->Muon_sumChargedParticlePt04(Muon_index_2),Ntp->Muon_sumChargedParticlePt04(Muon_index_3))));
    Min_Isolation_sumNeutralHadronEt04.at(t).Fill(std::min(Ntp->Muon_sumNeutralHadronEt04(Muon_index_1), std::min(Ntp->Muon_sumNeutralHadronEt04(Muon_index_2),Ntp->Muon_sumNeutralHadronEt04(Muon_index_3))));
    Min_Isolation_sumNeutralHadronEtHighThreshold04.at(t).Fill(std::min(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_1), std::min(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_2),Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_3))));
    Min_Isolation_sumPhotonEt04.at(t).Fill(std::min(Ntp->Muon_sumPhotonEt04(Muon_index_1), std::min(Ntp->Muon_sumPhotonEt04(Muon_index_2),Ntp->Muon_sumPhotonEt04(Muon_index_3))));
    Min_Isolation_sumPhotonEtHighThreshold04.at(t).Fill(std::min(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_1), std::min(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_2),Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_3))));
    Min_Isolation_sumPUPt04.at(t).Fill(std::min(Ntp->Muon_sumPUPt04(Muon_index_1), std::min(Ntp->Muon_sumPUPt04(Muon_index_2),Ntp->Muon_sumPUPt04(Muon_index_3))));

	 // ------------------------
  Muon1_calEnergy_hadS9.at(t).Fill(Ntp->Muon_calEnergy_hadS9(Muon_index_1));
  Muon1_calEnergy_had.at(t).Fill(Ntp->Muon_calEnergy_had(Muon_index_1));
  Muon1_calEnergy_emS25.at(t).Fill(Ntp->Muon_calEnergy_emS25(Muon_index_1));
  Muon1_calEnergy_emS9.at(t).Fill(Ntp->Muon_calEnergy_emS9(Muon_index_1));
  Muon1_calEnergy_em.at(t).Fill(Ntp->Muon_calEnergy_em(Muon_index_1));

  Muon2_calEnergy_hadS9.at(t).Fill(Ntp->Muon_calEnergy_hadS9(Muon_index_2));
  Muon2_calEnergy_had.at(t).Fill(Ntp->Muon_calEnergy_had(Muon_index_2));
  Muon2_calEnergy_emS25.at(t).Fill(Ntp->Muon_calEnergy_emS25(Muon_index_2));
  Muon2_calEnergy_emS9.at(t).Fill(Ntp->Muon_calEnergy_emS9(Muon_index_2));
  Muon2_calEnergy_em.at(t).Fill(Ntp->Muon_calEnergy_em(Muon_index_2));

  Muon3_calEnergy_hadS9.at(t).Fill(Ntp->Muon_calEnergy_hadS9(Muon_index_3));
  Muon3_calEnergy_had.at(t).Fill(Ntp->Muon_calEnergy_had(Muon_index_3));
  Muon3_calEnergy_emS25.at(t).Fill(Ntp->Muon_calEnergy_emS25(Muon_index_3));
  Muon3_calEnergy_emS9.at(t).Fill(Ntp->Muon_calEnergy_emS9(Muon_index_3));
  Muon3_calEnergy_em.at(t).Fill(Ntp->Muon_calEnergy_em(Muon_index_3));

  Min_calEnergy_hadS9.at(t).Fill(std::min(Ntp->Muon_calEnergy_hadS9(Muon_index_1), std::min(Ntp->Muon_calEnergy_hadS9(Muon_index_2),Ntp->Muon_calEnergy_hadS9(Muon_index_3))));
  Min_calEnergy_had.at(t).Fill(std::min(Ntp->Muon_calEnergy_had(Muon_index_1), std::min(Ntp->Muon_calEnergy_had(Muon_index_2),Ntp->Muon_calEnergy_had(Muon_index_3))));
  Min_calEnergy_emS25.at(t).Fill(std::min(Ntp->Muon_calEnergy_emS25(Muon_index_1), std::min(Ntp->Muon_calEnergy_emS25(Muon_index_2),Ntp->Muon_calEnergy_emS25(Muon_index_3))));
  Min_calEnergy_emS9.at(t).Fill(std::min(Ntp->Muon_calEnergy_emS9(Muon_index_1), std::min(Ntp->Muon_calEnergy_emS9(Muon_index_2),Ntp->Muon_calEnergy_emS9(Muon_index_3))));
  Min_calEnergy_em.at(t).Fill(std::min(Ntp->Muon_calEnergy_em(Muon_index_1), std::min(Ntp->Muon_calEnergy_em(Muon_index_2),Ntp->Muon_calEnergy_em(Muon_index_3))));

  Max_calEnergy_hadS9.at(t).Fill(std::max(Ntp->Muon_calEnergy_hadS9(Muon_index_1), std::max(Ntp->Muon_calEnergy_hadS9(Muon_index_2),Ntp->Muon_calEnergy_hadS9(Muon_index_3))));
  Max_calEnergy_had.at(t).Fill(std::max(Ntp->Muon_calEnergy_had(Muon_index_1), std::max(Ntp->Muon_calEnergy_had(Muon_index_2),Ntp->Muon_calEnergy_had(Muon_index_3))));
  Max_calEnergy_emS25.at(t).Fill(std::max(Ntp->Muon_calEnergy_emS25(Muon_index_1), std::max(Ntp->Muon_calEnergy_emS25(Muon_index_2),Ntp->Muon_calEnergy_emS25(Muon_index_3))));
  Max_calEnergy_emS9.at(t).Fill(std::max(Ntp->Muon_calEnergy_emS9(Muon_index_1), std::max(Ntp->Muon_calEnergy_emS9(Muon_index_2),Ntp->Muon_calEnergy_emS9(Muon_index_3))));
  Max_calEnergy_em.at(t).Fill(std::max(Ntp->Muon_calEnergy_em(Muon_index_1), std::max(Ntp->Muon_calEnergy_em(Muon_index_2),Ntp->Muon_calEnergy_em(Muon_index_3))));

	 // ------------------------
  }

}


void  TripleMuIsolation::Finish(){
 if(mode == RECONSTRUCT){
    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
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

    //    }
  }
  Selection::Finish();
}





