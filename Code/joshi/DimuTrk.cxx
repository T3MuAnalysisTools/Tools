#include "DimuTrk.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

double DimuTrk::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}

DimuTrk::DimuTrk(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

DimuTrk::~DimuTrk(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  DimuTrk::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1SeedOk)     cut.at(L1SeedOk)=1;
    if(i==HLTOk)        cut.at(HLTOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
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

    else if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

      // Track Candidate Information
      Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",66,-0.5,65.5,"p (track)","Events");
      Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",66,-0.5,65.5,"E (track)","Events");
      Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",66,-0.5,65.5,"p_{T} (track)","Events");
      Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",66,-0.5,65.5,"#eta","Events");
      Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",66,-0.5,65.5,"#phi","Events");
      Track_vx=HConfig.GetTH1D(Name+"_Track_vx","X coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex x coordinate (cm)","Events");
      Track_vy=HConfig.GetTH1D(Name+"_Track_vy","Y coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex y coordinate (cm)","Events");
      Track_vz=HConfig.GetTH1D(Name+"_Track_vz","Z coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex z coordinate (cm)","Events");
      Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",66,-0.5,65.5,"#chi^{2}/ndf (track fit)","Events");
      Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",66,-0.5,65.5,"n valid track hits","Events");
      Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",66,-0.5,65.5,"Number of Vertices","Events");
      Track_dxy=HConfig.GetTH1D(Name+"_Track_dxy","Transverse displacement of the parent vertex from the bs",66,-0.5,65.5,"dxy (cm)","Events");
      Track_dz=HConfig.GetTH1D(Name+"_Track_dz","Longitudnal displacement of the parent vertex from the bs",66,-0.5,65.5,"dz (cm)","Events");
      Track_dxyError=HConfig.GetTH1D(Name+"_Track_dxyError","dxy Error",66,-0.5,65.5,"#Deltadxy (cm)","Events");
      Track_dzError=HConfig.GetTH1D(Name+"_Track_dzError","dz Error",66,-0.5,65.5,"#Deltadz (cm)","Events");
      // Muon variables (Muons from dimuon + track candidates)
		 Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,50,"#mu_{1} p_{T} (GeV)","Events");
		 Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
		 Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{1} #phi","Events"); 
		 Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
		 Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  
		 Muon1_vx=HConfig.GetTH1D(Name+"_Muon1_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vx","Events"); 
		 Muon1_vy=HConfig.GetTH1D(Name+"_Muon1_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vy","Events"); 
		 Muon1_vz=HConfig.GetTH1D(Name+"_Muon1_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vz","Events");
  
  		Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,50,"#mu_{2} p_{T} (GeV)","Events");
  		Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
  		Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{2} #phi","Events"); 
  		Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
 		Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
  		Muon2_vx=HConfig.GetTH1D(Name+"_Muon2_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vx","Events"); 
  		Muon2_vy=HConfig.GetTH1D(Name+"_Muon2_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vy","Events"); 
  		Muon2_vz=HConfig.GetTH1D(Name+"_Muon2_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vz","Events");

//****
//Muon variables
  Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"#mu_{1} is Iso valid","Events");
  Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"#mu_{2} is Iso valid","Events");
  Muon1_isTimeValid=HConfig.GetTH1D(Name+"_Muon1_isTimeValid","#mu_{1} isTimevalid",2,-0.5,1.5,"#mu_{1} is Time valid","Events");
  Muon2_isTimeValid=HConfig.GetTH1D(Name+"_Muon2_isTimeValid","#mu_{2} isTimeValid",2,-0.5,1.5,"#mu_{2} is Time valid","Events");
  Muon1_emEt03=HConfig.GetTH1D(Name+"_Muon1_emEt03","",10,0,10,"#mu_{1} EM E_{T}03","Events");
  Muon2_emEt03=HConfig.GetTH1D(Name+"_Muon2_emEt03","",10,0,10,"#mu_{2} EM E_{T}03","Events");
  Muon1_emVetoEt03=HConfig.GetTH1D(Name+"_Muon1_emVetoEt03","",10,0,10,"#mu_{1} EM veto E_{T}03","Events");
  Muon2_emVetoEt03=HConfig.GetTH1D(Name+"_Muon2_emVetoEt03","",10,0,10,"#mu_{2} EM veto E_{T}03","Events");
  Muon1_hadEt03=HConfig.GetTH1D(Name+"_Muon1_hadEt03","",10,0,10,"#mu_{1} hadron E_{T}03","Events");
  Muon2_hadEt03=HConfig.GetTH1D(Name+"_Muon2_hadEt03","",10,0,10,"#mu_{2} hadron E_{T}03","Events");
  Muon1_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt03","",10,0,10,"#mu_{1} hadron veto E_{T}03","Events");
  Muon2_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt03","",10,0,10,"#mu_{2} hadron veto E_{T}03","Events");
  Muon1_nJets03=HConfig.GetTH1D(Name+"_Muon1_nJets03","",10,0,10,"#mu_{1} njets03","Events");
  Muon2_nJets03=HConfig.GetTH1D(Name+"_Muon2_nJets03","",10,0,10,"#mu_{2} njets03","Events");
  Muon1_nTracks03=HConfig.GetTH1D(Name+"_Muon1_nTracks03","",10,0,10,"#mu_{1} ntracks03","Events");
  Muon2_nTracks03=HConfig.GetTH1D(Name+"_Muon2_nTracks03","",10,0,10,"#mu_{2} ntracks03","Events");
  Muon1_sumPt03=HConfig.GetTH1D(Name+"_Muon1_sumPt03","",10,0,10,"","Events");
  Muon2_sumPt03=HConfig.GetTH1D(Name+"_Muon2_sumPt03","",10,0,10,"","Events");
  Muon1_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt03","",10,0,10,"","Events");
  Muon2_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt03","",10,0,10,"","Events");
  Muon1_emEt05=HConfig.GetTH1D(Name+"_Muon1_emEt05","",10,0,10,"#mu_{1} EM E_{T}05","Events");
  Muon2_emEt05=HConfig.GetTH1D(Name+"_Muon2_emEt05","",10,0,10,"#mu_{2} EM E_{T}05","Events");
  Muon1_emVetoEt05=HConfig.GetTH1D(Name+"_Muon1_emVetoEt05","#mu_{1} EM veto E_{T}05",10,0,10,"","Events");
  Muon2_emVetoEt05=HConfig.GetTH1D(Name+"_Muon2_emVetoEt05","#mu_{2} EM veto E_{T}05",10,0,10,"","Events");
  Muon1_hadEt05=HConfig.GetTH1D(Name+"_Muon1_hadEt05","",10,0,10,"#mu_{1} hadron E_{T}05","Events");
  Muon2_hadEt05=HConfig.GetTH1D(Name+"_Muon2_hadEt05","",10,0,10,"#mu_{2} hadron E_{T}05","Events");
  Muon1_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt05","",10,0,10,"#mu_{1} hadron veto E_{T}05","Events");
  Muon2_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt05","",10,0,10,"#mu_{2} hadron veto E_{T}05","Events");
  Muon1_nJets05=HConfig.GetTH1D(Name+"_Muon1_nJets05","",10,0,10,"#mu_{1} njets 05","Events");
  Muon2_nJets05=HConfig.GetTH1D(Name+"_Muon2_nJets05","",10,0,10,"#mu_{2} njets 05","Events");
  Muon1_nTracks05=HConfig.GetTH1D(Name+"_Muon1_nTracks05","",10,0,10,"#mu_{1} ntracks05","Events");
  Muon2_nTracks05=HConfig.GetTH1D(Name+"_Muon2_nTracks05","",10,0,10,"#mu_{2} ntracks05","Events");
  Muon1_sumPt05=HConfig.GetTH1D(Name+"_Muon1_sumPt05","",10,0,10,"#mu_{1} #Sigmap_{T}05","Events");
  Muon2_sumPt05=HConfig.GetTH1D(Name+"_Muon2_sumPt05","",10,0,10,"#mu_{2} #Sigmap_{T}05","Events");
  Muon1_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt05","",10,0,10,"#mu_{1} tracker veto p_{T}05","Events");
  Muon2_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt05","",10,0,10,"#mu_{2} tracker veto p_{T}05","Events");
  Muon1_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt03","#mu_{1} #Sigma charged had p_{T}03",10,0,10,"","Events");
  Muon2_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt03","#mu_{2} #Sigma charged had p_{T}03",10,0,10,"","Events");
  Muon1_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt03","",10,0,10,"#mu_{1} #Sigma charged particle p_{T}03","Events");
  Muon2_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt03","",10,0,10,"#mu_{2} #Sigma charged particle p_{T}03","Events");
  Muon1_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt03","",10,0,10,"#mu_{1} #Sigma neutral had E_{T}03","Events");
  Muon2_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt03","",10,0,10,"#mu_{2} #Sigma neutral had E_{T}03","Events");
  Muon1_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold03","",10,0,10,"#mu_{1} #Sigma neutral had E_{T} HT03","Events");
  Muon2_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold03","",10,0,10,"#mu_{2} #Sigma neutral had E_{T} HT03","Events");
  Muon1_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt03","",10,0,10,"#mu_{1} #Sigma#gamma E_{T}03","Events");
  Muon2_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt03","",10,0,10,"#mu_{2} #Sigma#gamma E_{T}03","Events");
  Muon1_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold03","",10,0,10,"#mu_{1} #Sigma#gamma E_{T}03","Events");
  Muon2_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold03","",10,0,10,"#mu_{2} #Sigma#gamma E_{T}03","Events");
  Muon1_sumPUPt03=HConfig.GetTH1D(Name+"_Muon1_sumPUPt03","",10,0,10,"#mu_{1} #Sigma PU p_{T}","Events");
  Muon2_sumPUPt03=HConfig.GetTH1D(Name+"_Muon2_sumPUPt03","",10,0,10,"#mu_{2} #Sigma PU p_{T}","Events");
  Muon1_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt04","#mu_{1} #Sigma charged had p_{T}04",10,0,10,"","Events");
  Muon2_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt04","#mu_{2} #Sigma charged had p_{T}04",10,0,10,"","Events");
  Muon1_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt04","",10,0,10,"#mu_{1} #Sigma charged particle p_{T}04","Events");
  Muon2_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt04","",10,0,10,"#mu_{2} #Sigma charged particle p_{T}04","Events");
  Muon1_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt04","",10,0,10,"#mu_{1} #Sigma neutral had E_{T}04","Events");
  Muon2_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt04","",10,0,10,"#mu_{2} #Sigma neutral had E_{T}04","Events");
  Muon1_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold04","",10,0,10,"#mu_{1} #Sigma neutral had E_{T} HT04","Events");
  Muon2_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold04","",10,0,10,"#mu_{2} #Sigma neutral had E_{T} HT04","Events");
  Muon1_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt04","",10,0,10,"#mu_{1} #Sigma#gammaE_{T}04","Events");
  Muon2_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt04","",10,0,10,"#mu_{2} #Sigma#gammaE_{T}04","Events");
  Muon1_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold04","",10,0,10,"#mu_{1} #Sigma#gammaE_{T} HT04","Events");
  Muon2_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold04","",10,0,10,"#mu_{2} #Sigma#gammaE_{T} HT04","Events");
  Muon1_sumPUPt04=HConfig.GetTH1D(Name+"_Muon1_sumPUPt04","",10,0,10,"#mu_{1} #SigmaPUPt04","Events");
  Muon2_sumPUPt04=HConfig.GetTH1D(Name+"_Muon2_sumPUPt04","",10,0,10,"#mu_{2} #SigmaPUPt04","Events");
  
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
      MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of Pt of two muons",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
      DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
      Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR","Events");
      Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR","Events");
      PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu invariant mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
      TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track invariant mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
      PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu invariant Mass vs. #mu#mu + track invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
    } 
    // Setup NPassed Histogams
    Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
    // Setup Extra Histograms
    // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
    NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");

    Selection::ConfigureHistograms(); //do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
 }

void  DimuTrk::Store_ExtraDist(){ 

//Track candidate variables 
    Extradist1d.push_back(&Track_P);
    Extradist1d.push_back(&Track_E);
    Extradist1d.push_back(&Track_Pt);
    Extradist1d.push_back(&Track_Eta);
    Extradist1d.push_back(&Track_Phi);
    Extradist1d.push_back(&Track_vx);
    Extradist1d.push_back(&Track_vy);
    Extradist1d.push_back(&Track_vz);
    Extradist1d.push_back(&Track_normalizedChi2);
    Extradist1d.push_back(&Track_numberOfValidHits);
    Extradist1d.push_back(&Track_charge);
    Extradist1d.push_back(&Track_dxy);
    Extradist1d.push_back(&Track_dz);
    Extradist1d.push_back(&Track_dxyError);
    Extradist1d.push_back(&Track_dzError);

    //Dimuon variable
	 cout<<"finni tracks"<<endl;
    Extradist1d.push_back(&Muon1_P);
    Extradist1d.push_back(&Muon1_E);
    Extradist1d.push_back(&Muon1_Pt);
    Extradist1d.push_back(&Muon1_Phi);
    Extradist1d.push_back(&Muon1_Eta);
    Extradist1d.push_back(&Muon1_vx);
    Extradist1d.push_back(&Muon1_vy);
    Extradist1d.push_back(&Muon1_vz);

	 Extradist1d.push_back(&Muon2_P);
    Extradist1d.push_back(&Muon2_E);
    Extradist1d.push_back(&Muon2_Pt);
    Extradist1d.push_back(&Muon2_Phi);
    Extradist1d.push_back(&Muon2_Eta);
    Extradist1d.push_back(&Muon2_vx);
    Extradist1d.push_back(&Muon2_vy);
    Extradist1d.push_back(&Muon2_vz);
    Extradist1d.push_back(&MuonsPtRatio);
    Extradist1d.push_back(&DimuondR);
    Extradist1d.push_back(&Muon1TrkdR);
    Extradist1d.push_back(&Muon2TrkdR);
    Extradist1d.push_back(&PhiMass);
    Extradist1d.push_back(&TripleMass);
    Extradist2d.push_back(&PhiMassVsDsMass);
    Extradist1d.push_back(&Muon1_isGlobal);
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
    //Extradist1d.push_back(&Track_TriggerMatchdR);
	 //Extradist1d.push_back(&Muon1_TriggerMatchdR);
	 //Extradist1d.push_back(&Muon2_TriggerMatchdR);

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


void  DimuTrk::doEvent(){ 
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

  value.at(is2MuTrk) = 0;
  if(Ntp->NTwoMuonsTrack()!=0 && Ntp->NThreeMuons() == 0) value.at(is2MuTrk) = 1;

  pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    double deltaMass(999.);
    unsigned int tmp_idx(0);

	for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
		tmp_idx = i2M;
		int mu1 = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0);
      int mu2 = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1);
      int track = Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0);

Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
    Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
    Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    Track_E.at(t).Fill(Ntp->Track_P4(track).E(),w);
    Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);
    Track_vx.at(t).Fill(Ntp->Track_Poca(track).X(),w);
    Track_vy.at(t).Fill(Ntp->Track_Poca(track).Y(),w);
    Track_vz.at(t).Fill(Ntp->Track_Poca(track).Z(),w);
    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
    Track_dxy.at(t).Fill(Ntp->Track_dxy(track),w);
    Track_dz.at(t).Fill(Ntp->Track_dz(track),w);
    Track_dxyError.at(t).Fill(Ntp->Track_dxyError(track),w);
    Track_dzError.at(t).Fill(Ntp->Track_dzError(track),w);
	 
	 Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
    Muon1_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu1),w);
    Muon2_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu2),w);
    Muon1_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu1),w);
    Muon2_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu2),w);
    Muon1_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu1),w);
    Muon2_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu2),w);
    Muon1_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu1),w);
    Muon2_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu2),w);
    Muon1_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu1),w);
    Muon2_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu2),w);
    Muon1_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu1),w);
    Muon2_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu2),w);
    Muon1_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu1),w);
    Muon2_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu2),w);
    Muon1_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu1),w);
    Muon2_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu2),w);
    Muon1_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu1),w);
    Muon2_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu2),w);
    Muon1_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu1),w);
    Muon2_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu2),w);
    Muon1_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu1),w);
    Muon2_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu2),w);
    Muon1_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu1),w);
    Muon2_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu2),w);
    Muon1_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu1),w);
    Muon2_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu2),w);
    Muon1_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu1),w);
    Muon2_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu2),w);
    Muon1_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu1),w);
    Muon2_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu2),w);
    Muon1_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu1),w);
    Muon2_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu2),w);
    Muon1_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu1),w);
    Muon2_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu2),w);
    Muon1_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu1),w);
    Muon2_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu2),w);
    Muon1_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu1),w);
    Muon2_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu2),w);
    Muon1_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu1),w);
    Muon2_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu2),w);
    Muon1_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu1),w);
    Muon2_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu2),w);
    Muon1_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu1),w);
    Muon2_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu2),w);
    Muon1_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu1),w);
    Muon2_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu2),w);
    Muon1_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu1),w);
    Muon2_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu2),w);
    Muon1_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu1),w);
    Muon2_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu2),w);
    Muon1_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu1),w);
    Muon2_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu2),w);
    Muon1_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu1),w);
    Muon2_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu2),w);
    Muon1_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu1),w);
    Muon2_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu2),w);
    Muon1_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu1),w);
    Muon2_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu2),w);
    Muon1_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu1),w);
    Muon2_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu2),w);
    Muon1_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu1),w);
    Muon2_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu2),w);
    Muon1_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu1),w);
    Muon2_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu2),w);
    Muon1_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu1),w);
    Muon2_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu2),w);
	 
	 Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(i2M),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(i2M),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(i2M),w);
    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(i2M),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(i2M),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(i2M),w);
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(i2M),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(i2M),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(i2M),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(i2M),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(i2M),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(i2M),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(i2M),w);
 
 // cout<<(Ntp->TwoMuonsTrack_TriggerMatch_dR).size()<<endl;
 // Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0),w);
 // Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1),w);
 // Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2),w);
    
DimuondR.at(t).Fill(deltaR(Ntp->Muon_P4(mu1).Eta(),Ntp->Muon_P4(mu1).Phi(),Ntp->Muon_P4(mu2).Eta(),Ntp->Muon_P4(mu2).Phi()));
    Muon1TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(mu1).Eta(),Ntp->Muon_P4(mu1).Phi(),Ntp->Track_P4(track).Eta(),Ntp->Track_P4(track).Phi()));
    Muon2TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(mu2).Eta(),Ntp->Muon_P4(mu2).Phi(),Ntp->Track_P4(track).Eta(),Ntp->Track_P4(track).Phi()));
    MuonsPtRatio.at(t).Fill(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0)).Pt()/Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1)).Pt(),w );
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+ 
    Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M(), w);
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
    Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M();
      PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

}
    /*
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
    unsigned int tmp_index(0);
    if(Ntp->NTwoMuonsTrack()!=0){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

      unsigned int muon_1 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
      unsigned int muon_2 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);

      if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
	deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
	tmp_index = i2M; // this is an index of the candidate with best mumu mass
      }
    }
    
    PhiMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_index).at(0))).M(), w);

    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_index).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_index).at(0))).M();

    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);
    }*/
  }
}

void  DimuTrk::Finish(){
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


