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
      title.at(i)="L1 Trigger ";
      hlabel="L1 Trigger ";
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
 
  // Track Information
  Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (all tracks)",66,-0.5,65.5,"","Events");
  Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (all tracks)",66,-0.5,65.5,"","Events");
  Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_vx=HConfig.GetTH1D(Name+"_Track_vx","X coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_vy=HConfig.GetTH1D(Name+"_Track_vy","Y coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_vz=HConfig.GetTH1D(Name+"_Track_vz","Z coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
  Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",66,-0.5,65.5,"Number of Vertices","Events");
  Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",66,-0.5,65.5,"Number of Vertices","Events");
  Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",66,-0.5,65.5,"Number of Vertices","Events");
  Track_dxy=HConfig.GetTH1D(Name+"_Track_dxy","Transverse displacement of the parent vertex from the bs",66,-0.5,65.5,"Number of Vertices","Events");
  Track_dz=HConfig.GetTH1D(Name+"_Track_dz","Longitudnal displacement of the parent vertex from the bs",66,-0.5,65.5,"Number of Vertices","Events");
  Track_dxyError=HConfig.GetTH1D(Name+"_Track_dxyError","dxy Error",66,-0.5,65.5,"Number of Vertices","Events");
  Track_dzError=HConfig.GetTH1D(Name+"_Track_dzError","dz Error",66,-0.5,65.5,"Number of Vertices","Events");

  //Muon Information
  Muon_Pt=HConfig.GetTH1D(Name+"_Muon_Pt","Transverse Pt (all muons)",25,0,50,"p_{T} (GeV)","Events");
  Muon_Eta=HConfig.GetTH1D(Name+"_Muon_Eta","Psuedorapidity (all muons)",25,-2.5,2.5,"#eta","Events");
  Muon_Phi=HConfig.GetTH1D(Name+"_Muon_Phi","Azimuthal angle of all muons",25,-3.15,3.15,"#phi","Events"); 
  Muon_E=HConfig.GetTH1D(Name+"_Muon_E","Energy of all muons",20,0,40,"Energy stored Muon_","Events");
  Muon_P=HConfig.GetTH1D(Name+"_Muon_P","Magnitude of momentum of all muons",20,0,40,"Azimuthal angle of all stored Muon_","Events");
  
  Muon_vx=HConfig.GetTH1D(Name+"_Muon_Vx","X coordinate of the parent vertex all muons",100,0,5,"vx","Events"); 
  Muon_vy=HConfig.GetTH1D(Name+"_Muon_Vy","Y coordinate of the parent vertex all muons",100,0,5,"vx","Events"); 
  Muon_vz=HConfig.GetTH1D(Name+"_Muon_Vz","Z coordinate of the parent vertex all muons",100,0,5,"vx","Events"); 
  
  Muon_IsGlobalMuon=HConfig.GetTH1D(Name+"_Muon_IsGlobal","GlobalMuon status (all muons)",2,0,2,"Global Muon","Events");
  Muon_IsStandAloneMuon=HConfig.GetTH1D(Name+"_Muon_IsStandAlone","StandAloneMuon status (all muons)",2,0,2,"StandAlone Muon","Events");
  Muon_IsTrackerMuon=HConfig.GetTH1D(Name+"_Muon_IsTrackerMuon","TrackerMuon Status (all muons)",2,0,2,"Tracker Muon","Events");
  Muon_IsCaloMuon=HConfig.GetTH1D(Name+"_Muon_","CaloMuon status (all muons)",2,0,2,"Calo Muon","Events");
  Muon_IsIsolationValid=HConfig.GetTH1D(Name+"_Muon_Phi","Isolation valid status",2,0,2,"Isolation Valid","Events");
  Muon_IsQualityValid=HConfig.GetTH1D(Name+"_Muon_Phi","Quality valid status",2,0,2,"Quality Valid","Events");
  Muon_IsTimeValid=HConfig.GetTH1D(Name+"_Muon_Phi","Time valid status",2,0,2,"Time Valid","Events");
  Muon_IsPFMuon=HConfig.GetTH1D(Name+"_Muon_Phi","ParticleFlowMuon status",2,0,2,"PF Muon","Events");
  Muon_IsRPCMuon=HConfig.GetTH1D(Name+"_Muon_Phi","RPCMuon status",2,0,2,"RPC Muon","Events");
  
  Muon_emEt03=HConfig.GetTH1D(Name+"_Muon_emEt03","Transverse energy in ECAL all muons",10,0,10,"EM E_{T}(03) (GeV)","Events");
  Muon_emVetoEt03=HConfig.GetTH1D(Name+"_Muon_VetoEt03","?",10,0,10,"EM Veto E_{T}(03) (GeV)","Events");
  Muon_hadEt03=HConfig.GetTH1D(Name+"_Muon_hadEt03","Tranverse HCAL energy of all muons",25,-3.15,3.15,"","Events");
  Muon_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon_VetoEt03","",25,-3.15,3.15,"","Events");
  Muon_nJets03=HConfig.GetTH1D(Name+"_Muon_nJets03","",25,-3.15,3.15,"","Events");
  Muon_nTracks03=HConfig.GetTH1D(Name+"_Muon_Tracks03","",25,-3.15,3.15,"","Events");
  Muon_sumPt03=HConfig.GetTH1D(Name+"_Muon_sumPt03","",25,-3.15,3.15,"","Events");
  Muon_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon_trackerVetoPt03","",25,-3.15,3.15,"","Events");
  
  Muon_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon_sumChargedHadronPt03","",25,-3.15,3.15,"","Events");              // sum-pt of charged Hadron
  Muon_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon_sumChargedParticlePt03","",25,-3.15,3.15,"","Events");            // sum-pt of charged Particles(inludes e/mu)
  Muon_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon_sumNeutralHadronEt03","",25,-3.15,3.15,"","Events");              // sum pt of neutral hadrons
  Muon_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon_sumNeutralHadronEtHighThreshold03","",25,-3.15,3.15,"","Events"); // sum pt of neutral hadrons with a higher threshold
  Muon_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon_sumPhotonEt03","",25,-3.15,3.15,"","Events");                     // sum pt of PF photons
  Muon_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon_PhotonEtHighThreshold03","",25,-3.15,3.15,"","Events");        // sum pt of PF photons with a higher threshold
  Muon_sumPUPt03=HConfig.GetTH1D(Name+"_Muon_sumPUPt03","",25,-3.15,3.15,"","Events");                         // sum pt of charged Particles not from PV (for Pu corrections)
  
  Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Muon_Phi","",25,-3.15,3.15,"","Events");
  Muon_Track_idx=HConfig.GetTH1D(Name+"_Muon_Phi","",25,-3.15,3.15,"","Events");
  
  Muon_combinedQuality_updatedSta=HConfig.GetTH1D(Name+"_combinedQuality_updatedSta","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_trkKink=HConfig.GetTH1D(Name+"_combinedQuality_trkKink","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_glbKink=HConfig.GetTH1D(Name+"_combinedQuality_glbKink","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_trkRelChi2=HConfig.GetTH1D(Name+"_combinedQuality_trkRelChi2","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_staRelChi2=HConfig.GetTH1D(Name+"_combinedQuality_staRelChi2","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_chi2LocalPosition=HConfig.GetTH1D(Name+"_combinedQuality_chi2LocalPosition","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_chi2LocalMomentum=HConfig.GetTH1D(Name+"_combinedQuality_chi2LocalMomentum","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_localDistance=HConfig.GetTH1D(Name+"_combinedQuality_localDistance","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_combinedQuality_globalDeltaEtaPhi","",25,-3.15,3.15,"","Events");
  Muon_combinedQuality_tightMatch=HConfig.GetTH1D(Name+"_combinedQuality_tightMatch","",2,0,2,"","Events");
  Muon_combinedQuality_glbTrackProbability=HConfig.GetTH1D(Name+"_combinedQuality_glbTrackProbability","",10,0,1,"P(global track)","Events");
  
  Muon_prod_inner_outer_charge=HConfig.GetTH1D(Name+"_Muon_ProdInnerOuterCharge","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_quality=HConfig.GetTH1D(Name+"_Muon_InnerTrack","",25,-3.15,3.15,"","Events");
  Muon_ptErrOverPt=HConfig.GetTH1D(Name+"_Muon_PtErrOverPt","",25,-3.15,3.15,"","Events");
  Muon_calEnergy_hadS9=HConfig.GetTH1D(Name+"_Muon_calEnergyhadS9","",25,-3.15,3.15,"","Events");
  Muon_calEnergy_had=HConfig.GetTH1D(Name+"_Muon_calEnergyhad","",25,-3.15,3.15,"","Events");
  Muon_calEnergy_emS25=HConfig.GetTH1D(Name+"_Muon_calEnergy","",25,-3.15,3.15,"","Events");
  Muon_calEnergy_emS9=HConfig.GetTH1D(Name+"_Muon_calEnergy","",25,-3.15,3.15,"","Events");
  Muon_calEnergy_em=HConfig.GetTH1D(Name+"_Muon_calEnergy","",25,-3.15,3.15,"","Events");
  
  Muon_charge=HConfig.GetTH1D(Name+"_Muon_charge","",25,-3.15,3.15,"","Events");
  Muon_trackCharge=HConfig.GetTH1D(Name+"_Muon_TrackCharge","",25,-3.15,3.15,"","Events");
  Muon_pdgid=HConfig.GetTH1D(Name+"_Muon_pdgid","",25,-3.15,3.15,"","Events");
  Muon_B=HConfig.GetTH1D(Name+"_Muon_B","",25,-3.15,3.15,"","Events");
  Muon_M=HConfig.GetTH1D(Name+"_Muon_M","",25,-3.15,3.15,"","Events");
  
  Muon_hitPattern_pixelLayerwithMeas=HConfig.GetTH1D(Name+"_Muon_hitPattern_pixelLayerwithMeas","",25,-3.15,3.15,"","Events");
  Muon_numberOfMatchedStations=HConfig.GetTH1D(Name+"_Muon_numberOfMatchedStations","",25,-3.15,3.15,"","Events");
  Muon_normChi2=HConfig.GetTH1D(Name+"_Muon_normChi2","",25,-3.15,3.15,"","Events");
  Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Muon_hitPattern_numberOfValidMuonHits","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_numberofValidHits=HConfig.GetTH1D(Name+"_Muon_innerTrack_numberofValidHits","",25,-3.15,3.15,"","Events");
  Muon_numberofValidPixelHits=HConfig.GetTH1D(Name+"_Muon_numberofValidPixelHits","",25,-3.15,3.15,"","Events");
  Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Muon_numberOfMatches","",25,-3.15,3.15,"","Events");
  Muon_trackerLayersWithMeasurement=HConfig.GetTH1D(Name+"_Muon_trackerLayersWithMeasurement","",25,-3.15,3.15,"","Events");
  Muon_segmentCompatibility=HConfig.GetTH1D(Name+"_Muon_segmentCompatibility","",25,-3.15,3.15,"","Events");
  Muon_caloCompatibility=HConfig.GetTH1D(Name+"_Muon_caloCompatibility","",25,-3.15,3.15,"","Events");
  
  Muon_innerTrack_validFraction=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_pixelLayersWithMeasurement=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_numberOfValidTrackerHits=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_numberOfLostTrackerHits=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_numberOfLostTrackerInnerHits=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_numberOfLostTrackerOuterHits=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  Muon_innerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Muon_innerTrack","",25,-3.15,3.15,"","Events");
  
  Muon_vmuonhitcomb_reco=HConfig.GetTH1D(Name+"_Muon_vmuonhitcomb","",25,-3.15,3.15,"","Events");
  Muon_rpchits_reco=HConfig.GetTH1D(Name+"_Muon_rpchits","",25,-3.15,3.15,"","Events");
  
  Muon_outerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Muon_outerTrack","",25,-3.15,3.15,"","Events");
  Muon_outerTrack_muonStationsWithValidHits=HConfig.GetTH1D(Name+"_Muon_outerTrack","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TM2DCompatibility=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TM2DCompatibility","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TrackerMuonArbitrated=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TrackerMuonArbitrated","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMOneStationTight=HConfig.GetTH1D(Name+"_Muon_isGood_TMOneStationTight","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMOneStationAngTight=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TMOneStationAngTight","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMLastStationTight=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TMLastStationTight","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMLastStationAngTight=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TMLastStationOptimizedLowTight","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMLastStationOptimizedLowPtTight=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TMLastStationOptimizedLowPtTight","",25,-3.15,3.15,"","Events");
  Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight=HConfig.GetTH1D(Name+"_Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight","",25,-3.15,3.15,"","Events");
  
  //Dimuon Information (Muons from dimuon + track candidates)
  MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of Pt of two muons",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
  DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR","Events");
  Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  
  //Muon variables (Muons from dimuon + track candidates)
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

  //
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&NTracks);
  Extradist1d.push_back(&HLT_Tau3Mu);  

  //Track variables 
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
 
  //Muon variables
  Extradist1d.push_back(&Muon_Pt);
  Extradist1d.push_back(&Muon_Eta);
  Extradist1d.push_back(&Muon_Phi);
  Extradist1d.push_back(&Muon_E);
  Extradist1d.push_back(&Muon_P);
  Extradist1d.push_back(&Muon_vx);
  Extradist1d.push_back(&Muon_vy);
  Extradist1d.push_back(&Muon_vz);
  Extradist1d.push_back(&Muon_IsGlobalMuon);
  Extradist1d.push_back(&Muon_IsStandAloneMuon);
  Extradist1d.push_back(&Muon_IsTrackerMuon);
  Extradist1d.push_back(&Muon_IsCaloMuon);
  Extradist1d.push_back(&Muon_IsIsolationValid);
  Extradist1d.push_back(&Muon_IsQualityValid);
  Extradist1d.push_back(&Muon_IsTimeValid);
  Extradist1d.push_back(&Muon_IsPFMuon);
  Extradist1d.push_back(&Muon_IsRPCMuon);
  Extradist1d.push_back(&Muon_emEt03);
  Extradist1d.push_back(&Muon_emVetoEt03);
  Extradist1d.push_back(&Muon_hadEt03);
  Extradist1d.push_back(&Muon_hadVetoEt03);
  Extradist1d.push_back(&Muon_nJets03);
  Extradist1d.push_back(&Muon_nTracks03);
  Extradist1d.push_back(&Muon_sumPt03);
  Extradist1d.push_back(&Muon_trackerVetoPt03);
  Extradist1d.push_back(&Muon_sumChargedHadronPt03);
  Extradist1d.push_back(&Muon_sumChargedParticlePt03);
  Extradist1d.push_back(&Muon_sumNeutralHadronEt03);
  Extradist1d.push_back(&Muon_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Muon_sumPhotonEt03);
  Extradist1d.push_back(&Muon_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Muon_sumPUPt03);
  Extradist1d.push_back(&Muon_numberOfChambers);
  Extradist1d.push_back(&Muon_Track_idx);
  Extradist1d.push_back(&Muon_combinedQuality_updatedSta);
  Extradist1d.push_back(&Muon_combinedQuality_trkKink);
  Extradist1d.push_back(&Muon_combinedQuality_glbKink);
  Extradist1d.push_back(&Muon_combinedQuality_trkRelChi2);
  Extradist1d.push_back(&Muon_combinedQuality_staRelChi2);
  Extradist1d.push_back(&Muon_combinedQuality_chi2LocalPosition);
  Extradist1d.push_back(&Muon_combinedQuality_chi2LocalMomentum);
  Extradist1d.push_back(&Muon_combinedQuality_localDistance);
  Extradist1d.push_back(&Muon_combinedQuality_globalDeltaEtaPhi);
  Extradist1d.push_back(&Muon_combinedQuality_tightMatch);
  Extradist1d.push_back(&Muon_combinedQuality_glbTrackProbability);
  Extradist1d.push_back(&Muon_prod_inner_outer_charge);
  Extradist1d.push_back(&Muon_innerTrack_quality);
  Extradist1d.push_back(&Muon_ptErrOverPt);
  Extradist1d.push_back(&Muon_calEnergy_hadS9);
  Extradist1d.push_back(&Muon_calEnergy_had);
  Extradist1d.push_back(&Muon_calEnergy_emS25);
  Extradist1d.push_back(&Muon_calEnergy_emS9);
  Extradist1d.push_back(&Muon_calEnergy_em);
  Extradist1d.push_back(&Muon_charge);
  Extradist1d.push_back(&Muon_trackCharge);
  Extradist1d.push_back(&Muon_pdgid);
  Extradist1d.push_back(&Muon_B);
  Extradist1d.push_back(&Muon_M);
  Extradist1d.push_back(&Muon_hitPattern_pixelLayerwithMeas);
  Extradist1d.push_back(&Muon_numberOfMatchedStations);
  Extradist1d.push_back(&Muon_normChi2);
  Extradist1d.push_back(&Muon_hitPattern_numberOfValidMuonHits);
  Extradist1d.push_back(&Muon_innerTrack_numberofValidHits);
  Extradist1d.push_back(&Muon_numberofValidPixelHits);
  Extradist1d.push_back(&Muon_numberOfMatches);
  Extradist1d.push_back(&Muon_trackerLayersWithMeasurement);
  Extradist1d.push_back(&Muon_segmentCompatibility);
  Extradist1d.push_back(&Muon_caloCompatibility);
  Extradist1d.push_back(&Muon_innerTrack_validFraction);
  Extradist1d.push_back(&Muon_innerTrack_pixelLayersWithMeasurement);
  Extradist1d.push_back(&Muon_innerTrack_numberOfValidTrackerHits);
  Extradist1d.push_back(&Muon_innerTrack_numberOfLostTrackerHits);
  Extradist1d.push_back(&Muon_innerTrack_numberOfLostTrackerInnerHits);
  Extradist1d.push_back(&Muon_innerTrack_numberOfLostTrackerOuterHits);
  Extradist1d.push_back(&Muon_innerTrack_normalizedChi2);
  Extradist1d.push_back(&Muon_vmuonhitcomb_reco);
  Extradist1d.push_back(&Muon_rpchits_reco);
  Extradist1d.push_back(&Muon_outerTrack_normalizedChi2);
  Extradist1d.push_back(&Muon_outerTrack_muonStationsWithValidHits);
  Extradist1d.push_back(&Muon_isGoodMuon_TM2DCompatibility);
  Extradist1d.push_back(&Muon_isGoodMuon_TrackerMuonArbitrated);
  Extradist1d.push_back(&Muon_isGoodMuon_TMOneStationTight);
  Extradist1d.push_back(&Muon_isGoodMuon_TMOneStationAngTight);
  Extradist1d.push_back(&Muon_isGoodMuon_TMLastStationTight);
  Extradist1d.push_back(&Muon_isGoodMuon_TMLastStationAngTight);
  Extradist1d.push_back(&Muon_isGoodMuon_TMLastStationOptimizedLowPtTight);
  Extradist1d.push_back(&Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);

  //Dimuon variables
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
  for (unsigned int iTrk=0; iTrk < Ntp->NTracks(); ++iTrk){
	 Track_Pt.at(t).Fill(Ntp->Track_P4(iTrk).Pt(),w);
    Track_Eta.at(t).Fill(Ntp->Track_P4(iTrk).Eta(),w);
    Track_Phi.at(t).Fill(Ntp->Track_P4(iTrk).Phi(),w);
    Track_E.at(t).Fill(Ntp->Track_P4(iTrk).E(),w);
    Track_P.at(t).Fill(Ntp->Track_P4(iTrk).P(),w);
    Track_vx.at(t).Fill(Ntp->Track_Poca(iTrk).X(),w);
    Track_vy.at(t).Fill(Ntp->Track_Poca(iTrk).Y(),w);
    Track_vz.at(t).Fill(Ntp->Track_Poca(iTrk).Z(),w);
 	 Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(iTrk),w);
    Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(iTrk),w);
    Track_charge.at(t).Fill(Ntp->Track_charge(iTrk),w);
    Track_dxy.at(t).Fill(Ntp->Track_dxy(iTrk),w);
    Track_dz.at(t).Fill(Ntp->Track_dz(iTrk),w);
    Track_dxyError.at(t).Fill(Ntp->Track_dxyError(iTrk),w);
    Track_dzError.at(t).Fill(Ntp->Track_dzError(iTrk),w);
  }

  for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++){
    Muon_Pt.at(t).Fill(Ntp->Muon_P4(iMuon).Pt(),w);
    Muon_Eta.at(t).Fill(Ntp->Muon_P4(iMuon).Eta(),w);
    Muon_Phi.at(t).Fill(Ntp->Muon_P4(iMuon).Phi(),w);
    Muon_E.at(t).Fill(Ntp->Muon_P4(iMuon).E(),w);
    Muon_P.at(t).Fill(Ntp->Muon_P4(iMuon).P(),w);
    Muon_vx.at(t).Fill(Ntp->Muon_Poca(iMuon).X(),w);
    Muon_vy.at(t).Fill(Ntp->Muon_Poca(iMuon).Y(),w);
    Muon_vz.at(t).Fill(Ntp->Muon_Poca(iMuon).Z(),w);
    Muon_IsGlobalMuon.at(t).Fill(Ntp->Muon_isGlobalMuon(iMuon),w);
    Muon_IsStandAloneMuon.at(t).Fill(Ntp->Muon_isStandAloneMuon(iMuon),w);
    Muon_IsTrackerMuon.at(t).Fill(Ntp->Muon_isTrackerMuon(iMuon),w);
    Muon_IsCaloMuon.at(t).Fill(Ntp->Muon_isCaloMuon(iMuon),w);
    Muon_IsIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(iMuon),w);
    Muon_IsQualityValid.at(t).Fill(Ntp->Muon_isQualityValid(iMuon),w);
    Muon_IsTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(iMuon),w);
    Muon_IsPFMuon.at(t).Fill(Ntp->Muon_isPFMuon(iMuon),w);
    Muon_IsRPCMuon.at(t).Fill(Ntp->Muon_isRPCMuon(iMuon),w);
    Muon_emEt03.at(t).Fill(Ntp->Muon_emEt03(iMuon),w);
    Muon_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(iMuon),w);
    Muon_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(iMuon),w);
    Muon_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(iMuon),w);
    Muon_nJets03.at(t).Fill(Ntp->Muon_nJets03(iMuon),w);
    Muon_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(iMuon),w);
    Muon_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(iMuon),w);
    Muon_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(iMuon),w);
    Muon_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(iMuon),w);
    Muon_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(iMuon),w);
    Muon_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(iMuon),w);
    Muon_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(iMuon),w);
    Muon_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(iMuon),w);
    Muon_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(iMuon),w);
    Muon_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(iMuon),w);
    Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(iMuon),w);
    Muon_Track_idx.at(t).Fill(Ntp->Muon_Track_idx(iMuon),w);
    Muon_combinedQuality_updatedSta.at(t).Fill(Ntp->Muon_combinedQuality_updatedSta(iMuon),w);
    Muon_combinedQuality_trkKink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(iMuon),w);
    Muon_combinedQuality_glbKink.at(t).Fill(Ntp->Muon_combinedQuality_glbKink(iMuon),w);
    Muon_combinedQuality_trkRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(iMuon),w);
    Muon_combinedQuality_staRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(iMuon),w);
    Muon_combinedQuality_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(iMuon),w);
    Muon_combinedQuality_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(iMuon),w);
    Muon_combinedQuality_localDistance.at(t).Fill(Ntp->Muon_combinedQuality_localDistance(iMuon),w);
    Muon_combinedQuality_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(iMuon),w);
    Muon_combinedQuality_tightMatch.at(t).Fill(Ntp->Muon_combinedQuality_tightMatch(iMuon),w);
    Muon_combinedQuality_glbTrackProbability.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(iMuon),w);
    Muon_prod_inner_outer_charge.at(t).Fill(Ntp->Muon_prod_inner_outer_charge(iMuon),w);
    Muon_innerTrack_quality.at(t).Fill(Ntp->Muon_innerTrack_quality(iMuon),w);
    Muon_ptErrOverPt.at(t).Fill(Ntp->Muon_ptErrOverPt(iMuon),w);
    Muon_calEnergy_hadS9.at(t).Fill(Ntp->Muon_calEnergy_hadS9(iMuon),w);
    Muon_calEnergy_had.at(t).Fill(Ntp->Muon_calEnergy_had(iMuon),w);
    Muon_calEnergy_emS25.at(t).Fill(Ntp->Muon_calEnergy_emS25(iMuon),w);
    Muon_calEnergy_emS9.at(t).Fill(Ntp->Muon_calEnergy_emS9(iMuon),w);
    Muon_calEnergy_em.at(t).Fill(Ntp->Muon_calEnergy_em(iMuon),w);
    Muon_charge.at(t).Fill(Ntp->Muon_charge(iMuon),w);
    Muon_trackCharge.at(t).Fill(Ntp->Muon_trackCharge(iMuon),w);
    Muon_pdgid.at(t).Fill(Ntp->Muon_pdgid(iMuon),w);
    Muon_B.at(t).Fill(Ntp->Muon_B(iMuon),w);
    Muon_M.at(t).Fill(Ntp->Muon_M(iMuon),w);
    Muon_hitPattern_pixelLayerwithMeas.at(t).Fill(Ntp->Muon_hitPattern_pixelLayerwithMeas(iMuon),w);
    Muon_numberOfMatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(iMuon),w);
    Muon_normChi2.at(t).Fill(Ntp->Muon_normChi2(iMuon),w);
    Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(iMuon),w);
    Muon_innerTrack_numberofValidHits.at(t).Fill(Ntp->Muon_innerTrack_numberofValidHits(iMuon),w);
    Muon_numberofValidPixelHits.at(t).Fill(Ntp->Muon_numberofValidPixelHits(iMuon),w);
    Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(iMuon),w);
    Muon_trackerLayersWithMeasurement.at(t).Fill(Ntp->Muon_trackerLayersWithMeasurement(iMuon),w);
    Muon_segmentCompatibility.at(t).Fill(Ntp->Muon_segmentCompatibility(iMuon),w);
    Muon_caloCompatibility.at(t).Fill(Ntp->Muon_caloCompatibility(iMuon),w);
    Muon_innerTrack_validFraction.at(t).Fill(Ntp->Muon_innerTrack_validFraction(iMuon),w);
    Muon_innerTrack_pixelLayersWithMeasurement.at(t).Fill(Ntp->Muon_innerTrack_pixelLayersWithMeasurement(iMuon),w);
    Muon_innerTrack_numberOfValidTrackerHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(iMuon),w);
    Muon_innerTrack_numberOfLostTrackerHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerHits(iMuon),w);
    Muon_innerTrack_numberOfLostTrackerInnerHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(iMuon),w);
    Muon_innerTrack_numberOfLostTrackerOuterHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(iMuon),w);
    Muon_innerTrack_normalizedChi2.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(iMuon),w);
    Muon_vmuonhitcomb_reco.at(t).Fill(Ntp->Muon_vmuonhitcomb_reco(iMuon),w);
    Muon_rpchits_reco.at(t).Fill(Ntp->Muon_rpchits_reco(iMuon),w);
    Muon_outerTrack_normalizedChi2.at(t).Fill(Ntp->Muon_outerTrack_normalizedChi2(iMuon),w);
    Muon_outerTrack_muonStationsWithValidHits.at(t).Fill(Ntp->Muon_outerTrack_muonStationsWithValidHits(iMuon),w);
    Muon_isGoodMuon_TM2DCompatibility.at(t).Fill(Ntp->Muon_isGoodMuon_TM2DCompatibility(iMuon),w);
    Muon_isGoodMuon_TrackerMuonArbitrated.at(t).Fill(Ntp->Muon_isGoodMuon_TrackerMuonArbitrated(iMuon),w);
    Muon_isGoodMuon_TMOneStationTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMOneStationTight(iMuon),w);
    Muon_isGoodMuon_TMOneStationAngTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMOneStationAngTight(iMuon),w);
    Muon_isGoodMuon_TMLastStationTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMLastStationTight(iMuon),w);
    Muon_isGoodMuon_TMLastStationAngTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMLastStationAngTight(iMuon),w);
    Muon_isGoodMuon_TMLastStationOptimizedLowPtTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(iMuon),w);
    Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.at(t).Fill(Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(iMuon),w);
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
