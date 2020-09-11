#include "DsPhiPeak_2018.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


DsPhiPeak_2018::DsPhiPeak_2018(TString Name_, TString id_):
  Selection(Name_,id_),
  dsMassMin(1.68),
  dsMassMax(2.02)
{



  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2018";
  PUWeightFile = new TFile(basedir+"/PUWeights_Run2018.root");
  puWeights = (TH1D*)PUWeightFile->Get("h1_weights");

  // This is a class constructor;
}

DsPhiPeak_2018::~DsPhiPeak_2018(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  DsPhiPeak_2018::Configure(){

  readerMuIDBarrel= new TMVA::Reader( "!Color:!Silent" );
  TString basedir = "";



  basedir = (TString)std::getenv("workdir")+"Code/CommonFiles/";
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  //  readerMuIDBarrel->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDBarrel->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDBarrel->AddVariable("mu_trackerLayersWithMeasurement" ,&mu_trackerLayersWithMeasurement );
  readerMuIDBarrel->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDBarrel->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDBarrel->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDBarrel->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDBarrel->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDBarrel->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDBarrel->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDBarrel->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDBarrel->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDBarrel->BookMVA( "BDT", basedir+"MuonID/MuonMVA_04june_TLWM_barrel/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerMuIDEndcap= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  //  readerMuIDEndcap->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDEndcap->AddVariable("mu_trackerLayersWithMeasurement" ,&mu_trackerLayersWithMeasurement );
  readerMuIDEndcap->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDEndcap->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDEndcap->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDEndcap->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDEndcap->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDEndcap->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDEndcap->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDEndcap->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDEndcap->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDEndcap->BookMVA( "BDT", basedir+"MuonID/MuonMVA_04june_TLWM_endcap/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles






  reader_trackerMuonId = new TMVA::Reader("!Color:!Silent");
  reader_trackerMuonId->AddSpectator("fake",&fake);
  reader_trackerMuonId->AddSpectator("muonPt",&muonPt);
  reader_trackerMuonId->AddSpectator("muonEta",&muonEta);
  reader_trackerMuonId->AddSpectator("muonPhi",&muonPhi);
  reader_trackerMuonId->AddVariable("muonInnerNC2",&muonInnerNC2);
  reader_trackerMuonId->AddVariable("muonInnerNValidHits",&muonInnerNValidHits);
  reader_trackerMuonId->AddVariable("muonValidFraction", &muonValidFraction);
  reader_trackerMuonId->AddVariable("muonNLostTrackerHits",&muonNLostTrackerHits);
  reader_trackerMuonId->AddVariable("muonNLostTrackerHitsInner",&muonNLostTrackerHitsInner);
  reader_trackerMuonId->AddVariable("muonNLostTrackerHitsOuter",&muonNLostTrackerHitsOuter);
  reader_trackerMuonId->AddVariable("muonPixelLayers",&muonPixelLayers);
  reader_trackerMuonId->AddVariable("muonNMatchedStations",&muonNMatchedStations);
  reader_trackerMuonId->AddVariable("muonPtErrPt",&muonPtErrPt);
  reader_trackerMuonId->AddVariable("muonSegComp",&muonSegComp);
  reader_trackerMuonId->AddVariable("muonCaloComp",&muonCaloComp);
  reader_trackerMuonId->AddVariable("muonHad",&muonHad);
  reader_trackerMuonId->AddVariable("muonEM",&muonEM);

  reader_trackerMuonId->BookMVA( "BDT", basedir+"weights/TrackerMuonBDT_2018/weights/MuPiTMVA_2018_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles


  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1TOk)           cut.at(L1TOk)=1;
    if(i==HLTOk)           cut.at(HLTOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==GlobalMu)        cut.at(GlobalMu)=1;
    if(i==MuCharge)        cut.at(MuCharge)=1;
    if(i==Mass2Mu)         cut.at(Mass2Mu)=1;
    if(i==TriggerMatch)    cut.at(TriggerMatch)=1;
    if(i==Mu1pt)           cut.at(Mu1pt)=3;
    if(i==Mu2pt)           cut.at(Mu2pt)=3;
    if(i==Trkpt)           cut.at(Trkpt)=1.8;
    if(i==DsMassCut)       cut.at(DsMassCut)=1;

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==L1TOk){
      title.at(i)="L1 DoubleMu0er1p5";
      hlabel="Level 1 Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1TOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1TOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLTOk){
      title.at(i)="DoubleMu3TrkTau3mu";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==GlobalMu){
      title.at(i)="GL and PF muons";
      hlabel="Muons are Global";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
    else if(i==MuCharge){
      title.at(i)="Muons opposite charge";
      hlabel="Muons have opposite charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mass2Mu){
      title.at(i)="0.99 $<$ $M_{2\\mu}$ $<$ 1.048 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Invariant mass of 2 muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
    }
    else if(i==TriggerMatch){
      title.at(i)="TriggerMatch";
      hlabel="TriggerMatch";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
    /*    else if(i==Mu1dR){
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
      }*/
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
    else if(i==DsMassCut){
      title.at(i)="Ds Mass Cut [1.68, 2.02]";
      htitle=title.at(i);
      hlabel="Ds Mass Cut";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
    }


  }

  // Track Candidate Information
  //Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",36,-0.5,35.5,"p (track)","Events");
  //Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",36,-0.5,35.5,"E (track)","Events");
  //Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  //Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  //Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",30,-3.5,3.5,"#phi","Events");
  Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",20,-0.5,4.5,"#chi^{2} (track fit)","Events");
  //Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",36,-0.5,35.5,"n valid track hits","Events");
  //Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",3,-1.5,1.5,"Track charge","Events");
  
  //Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  //Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  //Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); 
  //Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
  //Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  
  //
  //Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  //Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
  //Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{2} #phi","Events"); 
  //Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
  //Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
  
  Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","Trigger Matching mu1",50,0,0.05,"#Delta R Trigger Match #mu_{1}","Events");
  Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","Trigger Matching mu2",50,0,0.05,"#Delta R Trigger Match #mu_{2}","Events");
  Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","Trigger Matching track",50,0,0.05,"#Delta R Trigger Match track","Events");
  
  //Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  //Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  //Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  //Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  //Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  //Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  //Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  //Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  
  
  //DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  //Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,1,"dR","Events");
  //Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,1,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu invariant mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  PhiPlusTrackMass=HConfig.GetTH1D(Name+"_PhiPlusTrackMass","#mu#mu + track invariant mass",100,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu invariant Mass vs. #mu#mu + track invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
 
  DsMass=HConfig.GetTH1D(Name+"_DsMass","Ds invariant mass",100,1.7,2.1,"M_{Ds} (GeV)", "Events"); 
  DsPt_peak=HConfig.GetTH1D(Name+"_DsPt_peak","Transverse Pt (Ds) in Ds Peak",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events"); peakCollection.push_back(&DsPt_peak);
  DsPt_sideband=HConfig.GetTH1D(Name+"_DsPt_sideband","Transverse Pt (Ds) in Ds Sideband",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events"); sidebandCollection.push_back(&DsPt_sideband);
  Ds_Pt=HConfig.GetTH1D(Name+"_DsPt","Transverse Pt (Ds)",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events"); validationCollection.push_back(&Ds_Pt);

  DsP_peak=HConfig.GetTH1D(Name+"_DsP_peak","Transverse P (Ds) in Ds Peak",51,-0.5,50.5,"Ds p (GeV)", "Events"); peakCollection.push_back(&DsP_peak);
  DsP_sideband=HConfig.GetTH1D(Name+"_DsP_sideband","Transverse P (Ds) in Ds Sideband",51,-0.5,50.5,"Ds p (GeV)", "Events"); sidebandCollection.push_back(&DsP_sideband);
  Ds_P=HConfig.GetTH1D(Name+"_Ds_P","Transverse P (Ds)",51,-0.5,50.5,"Ds p (GeV)", "Events"); validationCollection.push_back(&Ds_P);

  DsM_peak=HConfig.GetTH1D(Name+"_DsM_peak","Ds invariant mass in Ds Peak",100,1.7,2.1,"M_{Ds} (GeV)", "Events"); peakCollection.push_back(&DsM_peak);
  DsM_sideband=HConfig.GetTH1D(Name+"_DsM_sideband","Ds invariant mass in Ds Sideband",100,1.7,2.1,"M_{Ds} (GeV)", "Events"); sidebandCollection.push_back(&DsM_sideband);
  Ds_M=HConfig.GetTH1D(Name+"_Ds_M","Ds invariant mass",100,1.7,2.1,"M_{Ds} (GeV)", "Events"); validationCollection.push_back(&Ds_M);

  DsL_peak=HConfig.GetTH1D(Name+"_DsL_peak","Length of Ds Decay in Ds Peak",25,0,1,"L_{Ds} (cm)", "Events"); peakCollection.push_back(&DsL_peak);
  DsL_sideband=HConfig.GetTH1D(Name+"_DsL_sideband","Length of Ds Decay in Ds Sideband",25,0,1,"L_{Ds} (cm)", "Events"); sidebandCollection.push_back(&DsL_sideband);
  Ds_L=HConfig.GetTH1D(Name+"_Ds_L","Length of Ds Decay",25,0,1,"L_{Ds} (cm)", "Events"); validationCollection.push_back(&Ds_L);

  DsEta_peak=HConfig.GetTH1D(Name+"_DsEta_peak","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events"); peakCollection.push_back(&DsEta_peak);
  DsEta_sideband=HConfig.GetTH1D(Name+"_DsEta_sideband","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events"); sidebandCollection.push_back(&DsEta_sideband);
  Ds_eta=HConfig.GetTH1D(Name+"_Ds_eta","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events"); validationCollection.push_back(&Ds_eta);


  Muon1MVAID_peak=HConfig.GetTH1D(Name+"_Muon1MVAID_peak","Muon1MVAID_peak",50,0.0,1.0,"#mu_{1} MVA","Events");peakCollection.push_back(&Muon1MVAID_peak);
  Muon2MVAID_peak=HConfig.GetTH1D(Name+"_Muon2MVAID_peak","Muon2MVAID_peak",50,0.0,1.0,"#mu_{2} MVA","Events");peakCollection.push_back(&Muon2MVAID_peak);

  Muon1MVAID_sideband=HConfig.GetTH1D(Name+"_Muon1MVAID_sideband","Muon1MVAID_sideband",50,0.0,1.0,"#mu_{1} MVA","Events");sidebandCollection.push_back(&Muon1MVAID_sideband);
  Muon2MVAID_sideband=HConfig.GetTH1D(Name+"_Muon2MVAID_sideband","Muon2MVAID_sideband",50,0.0,1.0,"#mu_{2} MVA","Events");sidebandCollection.push_back(&Muon2MVAID_sideband);


  Muon1MVAID=HConfig.GetTH1D(Name+"_Muon1MVAID","Muon1MVAID",50,0.0,1.0,"#mu_{1} MVA","Events");validationCollection.push_back(&Muon1MVAID);
  Muon2MVAID=HConfig.GetTH1D(Name+"_Muon2MVAID","Muon2MVAID",50,0.0,1.0,"#mu_{2} MVA","Events");validationCollection.push_back(&Muon2MVAID);

  mu1segmentCompatibility= HConfig.GetTH1D(Name+"_mu1segmentCompatibility","mu1segmentCompatibility",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
  validationCollection.push_back(&mu1segmentCompatibility);
  mu2segmentCompatibility= HConfig.GetTH1D(Name+"_mu2segmentCompatibility","mu2segmentCompatibility",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
  validationCollection.push_back(&mu2segmentCompatibility);


  mu1segmentCompatibility_peak= HConfig.GetTH1D(Name+"_mu1segmentCompatibility_peak","mu1segmentCompatibility_peak",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
  peakCollection.push_back(&mu1segmentCompatibility_peak);
  mu2segmentCompatibility_peak= HConfig.GetTH1D(Name+"_mu2segmentCompatibility_peak","mu2segmentCompatibility_peak",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
  peakCollection.push_back(&mu2segmentCompatibility_peak);


  mu1segmentCompatibility_sideband= HConfig.GetTH1D(Name+"_mu1segmentCompatibility_sideband","mu1segmentCompatibility_sideband",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
  sidebandCollection.push_back(&mu1segmentCompatibility_sideband);
  mu2segmentCompatibility_sideband= HConfig.GetTH1D(Name+"_mu2segmentCompatibility_sideband","mu2segmentCompatibility_sideband",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
  sidebandCollection.push_back(&mu2segmentCompatibility_sideband);



  DsGenMatch=HConfig.GetTH1D(Name+"_DsGenMatch","dR between Gen Ds to Track",50,0,.1,"dR","Events");


  DecayLength_peak=HConfig.GetTH1D(Name+"_DecayLength_peak","Proper Decay Length of Ds in Ds Peak",20,0,.1,"Proper Decay Length (cm)","Events"); peakCollection.push_back(&DecayLength_peak);
  DecayLength_sideband=HConfig.GetTH1D(Name+"_DecayLength_sideband","Proper Decay Length of Ds in sideband",20,0,.1,"Proper Decay Length (cm)","Events"); sidebandCollection.push_back(&DecayLength_sideband);
  DecayLength_prompt=HConfig.GetTH1D(Name+"_DecayLength_prompt","Proper Decay Length of Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_non_prompt=HConfig.GetTH1D(Name+"_DecayLength_non_prompt","Proper Decay Length of Non-Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength=HConfig.GetTH1D(Name+"_DecayLength","Proper Decay Length of Ds",20,0,.1,"Proper Decay Length (cm)","Events"); validationCollection.push_back(&DecayLength);

  Muon1_Pt_peak=HConfig.GetTH1D(Name+"_Muon1_Pt_peak","Transverse Pt in Ds Peak (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events"); peakCollection.push_back(&Muon1_Pt_peak);
  Muon1_Pt_sideband=HConfig.GetTH1D(Name+"_Muon1_Pt_sideband","Transverse Pt in Ds Sideband (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events"); sidebandCollection.push_back(&Muon1_Pt_sideband);
  Muon1_Eta_peak=HConfig.GetTH1D(Name+"_Muon1_Eta_peak","Psuedorapidity in Ds Peak (muon 1)",30,-2,2,"#mu_{1} #eta","Events"); peakCollection.push_back(&Muon1_Eta_peak);
  Muon1_Eta_sideband=HConfig.GetTH1D(Name+"_Muon1_Eta_sideband","Psuedorapidity in Ds Sideband (muon 1)",30,-2,2,"#mu_{1} #eta","Events"); sidebandCollection.push_back(&Muon1_Eta_sideband);
  Muon1_Phi_peak=HConfig.GetTH1D(Name+"_Muon1_Phi_peak","Azimuthal angle in Ds Peak of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); peakCollection.push_back(&Muon1_Phi_peak);
  Muon1_Phi_sideband=HConfig.GetTH1D(Name+"_Muon1_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); sidebandCollection.push_back(&Muon1_Phi_sideband);
  control_Muon1_Pt=HConfig.GetTH1D(Name+"_control_Muon1_Pt","Transverse Pt (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events"); validationCollection.push_back(&control_Muon1_Pt);
  control_Muon1_Eta=HConfig.GetTH1D(Name+"_control_Muon1_Eta","Psuedorapidity (muon 1)",30,-2,2,"#mu_{1} #eta","Events"); validationCollection.push_back(&control_Muon1_Eta);
  control_Muon1_Phi=HConfig.GetTH1D(Name+"_control_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); validationCollection.push_back(&control_Muon1_Phi);

  Muon2_Pt_peak=HConfig.GetTH1D(Name+"_Muon2_Pt_peak","Transverse Pt in Ds Peak (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events"); peakCollection.push_back(&Muon2_Pt_peak);
  Muon2_Pt_sideband=HConfig.GetTH1D(Name+"_Muon2_Pt_sideband","Transverse Pt in Ds Sideband (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events"); sidebandCollection.push_back(&Muon2_Pt_sideband);
  Muon2_Eta_peak=HConfig.GetTH1D(Name+"_Muon2_Eta_peak","Psuedorapidity in Ds Peak (muon 2)",30,-2,2,"#mu_{2} #eta","Events"); peakCollection.push_back(&Muon2_Eta_peak);
  Muon2_Eta_sideband=HConfig.GetTH1D(Name+"_Muon2_Eta_sideband","Psuedorapidity in Ds Sideband (muon 2)",30,-2,2,"#mu_{2} #eta","Events"); sidebandCollection.push_back(&Muon2_Eta_sideband);
  Muon2_Phi_peak=HConfig.GetTH1D(Name+"_Muon2_Phi_peak","Azimuthal angle in Ds Peak of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events"); peakCollection.push_back(&Muon2_Phi_peak);
  Muon2_Phi_sideband=HConfig.GetTH1D(Name+"_Muon2_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events"); sidebandCollection.push_back(&Muon2_Phi_sideband);
  control_Muon2_Pt=HConfig.GetTH1D(Name+"_control_Muon2_Pt","Transverse Pt (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events"); validationCollection.push_back(&control_Muon2_Pt);
  control_Muon2_Eta=HConfig.GetTH1D(Name+"_control_Muon2_Eta","Psuedorapidity (muon 2)",30,-2,2,"#mu_{2} #eta","Events"); validationCollection.push_back(&control_Muon2_Eta);
  control_Muon2_Phi=HConfig.GetTH1D(Name+"_control_Muon2_Phi","Azimuthal angle of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events"); validationCollection.push_back(&control_Muon2_Phi);

  Track_Pt_peak=HConfig.GetTH1D(Name+"_Track_Pt_peak","Transverse Pt in Ds Peak (Track)",25,0,20,"Track p_{T} (GeV)","Events"); peakCollection.push_back(&Track_Pt_peak);
  Track_Pt_sideband=HConfig.GetTH1D(Name+"_Track_Pt_sideband","Transverse Pt in Ds Sideband (Track)",25,0,20,"Track p_{T} (GeV)","Events"); sidebandCollection.push_back(&Track_Pt_sideband);
  Track_Eta_peak=HConfig.GetTH1D(Name+"_Track_Eta_peak","Psuedorapidity in Ds Peak (Track)",30,-2,2,"Track #eta","Events"); peakCollection.push_back(&Track_Eta_peak);
  Track_Eta_sideband=HConfig.GetTH1D(Name+"_Track_Eta_sideband","Psuedorapidity in Ds Sideband (Track)",30,-2,2,"Track #eta","Events"); sidebandCollection.push_back(&Track_Eta_sideband);
  Track_Phi_peak=HConfig.GetTH1D(Name+"_Track_Phi_peak","Azimuthal angle in Ds Peak of (Track)",25,-3.4,3.4,"Track #phi","Events"); peakCollection.push_back(&Track_Phi_peak);
  Track_Phi_sideband=HConfig.GetTH1D(Name+"_Track_Phi_sideband","Azimuthal angle in Ds Sideband of (Track)",25,-3.4,3.4,"Track #phi","Events"); sidebandCollection.push_back(&Track_Phi_sideband);
  control_Track_Pt=HConfig.GetTH1D(Name+"_control_Track_Pt","Transverse Pt (Track)",25,0,20,"Track p_{T} (GeV)","Events"); validationCollection.push_back(&control_Track_Pt);
  control_Track_Eta=HConfig.GetTH1D(Name+"_control_Track_Eta","Psuedorapidity (Track)",30,-2,2,"Track #eta","Events"); validationCollection.push_back(&control_Track_Eta);
  control_Track_Phi=HConfig.GetTH1D(Name+"_control_Track_Phi","Azimuthal angle of (Track)",25,-3.4,3.4,"Track #phi","Events"); validationCollection.push_back(&control_Track_Phi);

  Mu1combinedQuality_chi2LocalMomentum_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalMomentum_peak","Mu1combinedQuality_chi2LocalMomentum_peak",25,0,20,"#chi^{2} local momentum",""); peakCollection.push_back(&Mu1combinedQuality_chi2LocalMomentum_peak);

  Mu1combinedQuality_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalPosition_peak","Mu1combinedQuality_chi2LocaPosition_peak",25,0,20,"#chi^{2} local position",""); peakCollection.push_back(&Mu1combinedQuality_chi2LocalPosition_peak);

  Mu1combinedQuality_staRelChi2_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_staRelChi2_peak","Mu1combinedQuality_staRelChi2_peak",25,0,2,"station relative #chi^{2}",""); peakCollection.push_back(&Mu1combinedQuality_staRelChi2_peak);
  Mu1combinedQuality_trkRelChi2_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_trkRelChi2_peak","Mu1combinedQuality_trkRelChi2_peak",25,0,2,"track relative #chi^{2}",""); peakCollection.push_back(&Mu1combinedQuality_trkRelChi2_peak);
  Mu1combinedQuality_globalDeltaEtaPhi_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_globalDeltaEtaPhi_peak","Mu1combinedQuality_globalDeltaEtaPhi_peak",25,0,0.3,"gl #Delta#eta#phi",""); peakCollection.push_back(&Mu1combinedQuality_globalDeltaEtaPhi_peak);
  Mu1mu_combinedQuality_trkKink_peak=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_trkKink_peak","Mu1mu_combinedQuality_trkKink_peak",25,0,20,"trk kink",""); peakCollection.push_back(&Mu1mu_combinedQuality_trkKink_peak);
  
  Mu1mu_combinedQuality_glbKink_peak=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_glbKink_peak","Mu1mu_combinedQuality_glbKink_peak",25,0,5,"glb kink",""); peakCollection.push_back(&Mu1mu_combinedQuality_glbKink_peak);
  Mu1combinedQuality_glbTrackProbability_peak=HConfig.GetTH1D(Name+"_Mu1combinedQuality_glbTrackProbability_peak","Mu1combinedQuality_glbTrackProbability_peak",25,0,10,"glb track probability",""); peakCollection.push_back(&Mu1combinedQuality_glbTrackProbability_peak);
  Mu1Numberofvalidtrackerhits_peak=HConfig.GetTH1D(Name+"_Mu1Numberofvalidtrackerhits_peak","Mu1Numberofvalidtrackerhits_peak",20,-0.5,19.5,"n valid tracker hits",""); peakCollection.push_back(&Mu1Numberofvalidtrackerhits_peak);
  Mu1Numberofvalidpixelhits_peak=HConfig.GetTH1D(Name+"_Mu1Numberofvalidpixelhits_peak","Mu1Numberofvalidpixelhits_peak",20,-0.5,19.5,"n valid pixel hits",""); peakCollection.push_back(&Mu1Numberofvalidpixelhits_peak);

  Mu1trackerLayersWithMeasurement_peak=HConfig.GetTH1D(Name+"_Mu1trackerLayersWithMeasurement_peak","Mu1trackerLayersWithMeasurement_peak",15,-0.5,14.5,"#mu_{1} tracker layer with measurements",""); peakCollection.push_back(&Mu1trackerLayersWithMeasurement_peak);

  Mu1validMuonHitComb_peak=HConfig.GetTH1D(Name+"_Mu1validMuonHitComb_peak","Mu1validMuonHitComb_peak",70,-0.5,69.5,"valid hit combined",""); peakCollection.push_back(&Mu1validMuonHitComb_peak);
  Mu1numberOfMatchedStations_peak=HConfig.GetTH1D(Name+"_Mu1numberOfMatchedStations_peak","Mu1numberOfMatchedStations_peak",10,-0.5,9.5,"N matched stations",""); peakCollection.push_back(&Mu1numberOfMatchedStations_peak);
  Mu1segmentCompatibility_peak=HConfig.GetTH1D(Name+"_Mu1segmentCompatibility_peak","Mu1segmentCompatibility_peak",25,0,1,"#mu_{1} segment compatibility",""); peakCollection.push_back(&Mu1segmentCompatibility_peak);
  Mu1timeAtIpInOutErr_peak=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOutErr_peak","Mu1timeAtIpInOutErr_peak",25,0,3,"timeAtIpInOutErr",""); peakCollection.push_back(&Mu1timeAtIpInOutErr_peak);
  Mu1timeAtIpInOut_peak=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOut_peak","Mu1timeAtIpInOut_peak",25,0,3,"timeAtIpInOut",""); peakCollection.push_back(&Mu1timeAtIpInOut_peak);
  Mu1GLnormChi2_peak=HConfig.GetTH1D(Name+"_Mu1GLnormChi2_peak","Mu1GLnormChi2_peak",25,0,2,"global normalised #chi^{2}",""); peakCollection.push_back(&Mu1GLnormChi2_peak);
  Mu1innerTrack_normalizedChi2_peak=HConfig.GetTH1D(Name+"_Mu1innerTrack_normalizedChi2_peak","Mu1innerTrack_normalizedChi2_peak",25,0,20,"inner track #chi^{2}",""); peakCollection.push_back(&Mu1innerTrack_normalizedChi2_peak);
  Mu1outerTrack_normalizedChi2_peak=HConfig.GetTH1D(Name+"_Mu1outerTrack_normalizedChi2_peak","Mu1outerTrack_normalizedChi2_peak",25,0,20,"outer track #chi^{2}",""); peakCollection.push_back(&Mu1outerTrack_normalizedChi2_peak);
  Mu1innerTrack_validFraction_peak=HConfig.GetTH1D(Name+"_Mu1innerTrack_validFraction_peak","Mu1innerTrack_validFraction_peak",25,0,2,"inner Track valid Fraction",""); peakCollection.push_back(&Mu1innerTrack_validFraction_peak);
  
  Mu1combinedQuality_chi2LocalMomentum_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalMomentum_sideband","Mu1combinedQuality_chi2LocalMomentum_sideband",25,0,20,"#chi^{2} local momentum",""); sidebandCollection.push_back(&Mu1combinedQuality_chi2LocalMomentum_sideband);
  Mu1combinedQuality_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalPosition_sideband","Mu1combinedQuality_chi2LocalPosition_sideband",25,0,20,"#chi^{2} local position",""); sidebandCollection.push_back(&Mu1combinedQuality_chi2LocalPosition_sideband);
  Mu1combinedQuality_staRelChi2_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_staRelChi2_sideband","Mu1combinedQuality_staRelChi2_sideband",25,0,2,"station relative #chi^{2}",""); sidebandCollection.push_back(&Mu1combinedQuality_staRelChi2_sideband);
  Mu1combinedQuality_trkRelChi2_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_trkRelChi2_sideband","Mu1combinedQuality_trkRelChi2_sideband",25,0,2,"track relative #chi^{2}",""); sidebandCollection.push_back(&Mu1combinedQuality_trkRelChi2_sideband);
  Mu1combinedQuality_globalDeltaEtaPhi_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_globalDeltaEtaPhi_sideband","Mu1combinedQuality_globalDeltaEtaPhi_sideband",25,0,0.3,"gl #Delta#eta#phi",""); sidebandCollection.push_back(&Mu1combinedQuality_globalDeltaEtaPhi_sideband);
  Mu1mu_combinedQuality_trkKink_sideband=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_trkKink_sideband","mMu1u_combinedQuality_trkKink_sideband",25,0,20,"trk kink",""); sidebandCollection.push_back(&Mu1mu_combinedQuality_trkKink_sideband);
  Mu1mu_combinedQuality_glbKink_sideband=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_glbKink_sideband","Mu1mu_combinedQuality_glbKink_sideband",25,0,5,"glb kink",""); sidebandCollection.push_back(&Mu1mu_combinedQuality_glbKink_sideband);
  Mu1combinedQuality_glbTrackProbability_sideband=HConfig.GetTH1D(Name+"_Mu1combinedQuality_glbTrackProbability_sideband","Mu1combinedQuality_glbTrackProbability_sideband",25,0,10,"glb track probability",""); sidebandCollection.push_back(&Mu1combinedQuality_glbTrackProbability_sideband);
  Mu1Numberofvalidtrackerhits_sideband=HConfig.GetTH1D(Name+"_Mu1Numberofvalidtrackerhits_sideband","Mu1Numberofvalidtrackerhits_sideband",20,-0.5,19.5,"n valid tracker hits",""); sidebandCollection.push_back(&Mu1Numberofvalidtrackerhits_sideband);
  Mu1Numberofvalidpixelhits_sideband=HConfig.GetTH1D(Name+"_Mu1Numberofvalidpixelhits_sideband","Mu1Numberofvalidpixelhits_sideband",20,-0.5,19.5,"n valid pixel hits",""); sidebandCollection.push_back(&Mu1Numberofvalidpixelhits_sideband);
  Mu1trackerLayersWithMeasurement_sideband=HConfig.GetTH1D(Name+"_Mu1trackerLayersWithMeasurement_sideband","Mu1trackerLayersWithMeasurement_sideband",15,-0.5,14.5,"#mu_{1} tracker layer with measurements",""); sidebandCollection.push_back(&Mu1trackerLayersWithMeasurement_sideband);

  Mu1validMuonHitComb_sideband=HConfig.GetTH1D(Name+"_Mu1validMuonHitComb_sideband","Mu1validMuonHitComb_sideband",70,-0.5,69.5,"valid hit combined",""); sidebandCollection.push_back(&Mu1validMuonHitComb_sideband);
  Mu1numberOfMatchedStations_sideband=HConfig.GetTH1D(Name+"_Mu1numberOfMatchedStations_sideband","Mu1numberOfMatchedStations_sideband",10,-0.5,9.5,"N matched stations",""); sidebandCollection.push_back(&Mu1numberOfMatchedStations_sideband);
  Mu1segmentCompatibility_sideband=HConfig.GetTH1D(Name+"_Mu1segmentCompatibility_sideband","Mu1segmentCompatibility_sideband",25,0,1,"#mu_{1} segment compatibility",""); sidebandCollection.push_back(&Mu1segmentCompatibility_sideband);
  Mu1timeAtIpInOutErr_sideband=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOutErr_sideband","Mu1timeAtIpInOutErr_sideband",25,0,3,"timeAtIpInOutErr",""); sidebandCollection.push_back(&Mu1timeAtIpInOutErr_sideband);
  Mu1timeAtIpInOut_sideband=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOut_sideband","Mu1timeAtIpInOut_sideband",25,0,3,"timeAtIpInOut",""); sidebandCollection.push_back(&Mu1timeAtIpInOut_sideband);
  Mu1GLnormChi2_sideband=HConfig.GetTH1D(Name+"_Mu1GLnormChi2_sideband","Mu1GLnormChi2_sideband",25,0,2,"global normalised #chi^{2}",""); sidebandCollection.push_back(&Mu1GLnormChi2_sideband);
  Mu1innerTrack_normalizedChi2_sideband=HConfig.GetTH1D(Name+"_Mu1innerTrack_normalizedChi2_sideband","Mu1innerTrack_normalizedChi2_sideband",25,0,20,"inner track #chi^{2}",""); sidebandCollection.push_back(&Mu1innerTrack_normalizedChi2_sideband);
  Mu1outerTrack_normalizedChi2_sideband=HConfig.GetTH1D(Name+"_Mu1outerTrack_normalizedChi2_sideband","Mu1outerTrack_normalizedChi2_sideband",25,0,20,"outer track #chi^{2}",""); sidebandCollection.push_back(&Mu1outerTrack_normalizedChi2_sideband);
  Mu1innerTrack_validFraction_sideband=HConfig.GetTH1D(Name+"_Mu1innerTrack_validFraction_sideband","Mu1innerTrack_validFraction_sideband",25,0,2,"inner Track valid Fraction",""); sidebandCollection.push_back(&Mu1innerTrack_validFraction_sideband);
  



  Mu1combinedQuality_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalMomentum","Mu1combinedQuality_chi2LocalMomentum",25,0,20,"#chi^{2} local momentum",""); validationCollection.push_back(&Mu1combinedQuality_chi2LocalMomentum);
  Mu1combinedQuality_chi2LocalPosition=HConfig.GetTH1D(Name+"_Mu1combinedQuality_chi2LocalPosition","Mu1combinedQuality_chi2LocalPosition",25,0,20,"#chi^{2} local position",""); validationCollection.push_back(&Mu1combinedQuality_chi2LocalPosition);
  Mu1combinedQuality_staRelChi2=HConfig.GetTH1D(Name+"_Mu1combinedQuality_staRelChi2","Mu1combinedQuality_staRelChi2",25,0,2,"station relative #chi^{2}",""); validationCollection.push_back(&Mu1combinedQuality_staRelChi2);
  Mu1combinedQuality_trkRelChi2=HConfig.GetTH1D(Name+"_Mu1combinedQuality_trkRelChi2","Mu1combinedQuality_trkRelChi2",25,0,2,"track relative #chi^{2}",""); validationCollection.push_back(&Mu1combinedQuality_trkRelChi2);
  Mu1combinedQuality_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_Mu1combinedQuality_globalDeltaEtaPhi","Mu1combinedQuality_globalDeltaEtaPhi",25,0,0.3,"gl #Delta#eta#phi",""); validationCollection.push_back(&Mu1combinedQuality_globalDeltaEtaPhi);
  Mu1mu_combinedQuality_trkKink=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_trkKink","Mu1mu_combinedQuality_trkKink",25,0,20,"trk kink",""); validationCollection.push_back(&Mu1mu_combinedQuality_trkKink);
  Mu1mu_combinedQuality_glbKink=HConfig.GetTH1D(Name+"_Mu1mu_combinedQuality_glbKink","Mu1mu_combinedQuality_glbKink",25,0,5,"glb kink",""); validationCollection.push_back(&Mu1mu_combinedQuality_glbKink);
  Mu1combinedQuality_glbTrackProbability=HConfig.GetTH1D(Name+"_Mu1combinedQuality_glbTrackProbability","Mu1combinedQuality_glbTrackProbability",25,0,10,"glb track probability",""); validationCollection.push_back(&Mu1combinedQuality_glbTrackProbability);
  Mu1Numberofvalidtrackerhits=HConfig.GetTH1D(Name+"_Mu1Numberofvalidtrackerhits","Mu1Numberofvalidtrackerhits",20,-0.5,19.5,"n valid tracker hits",""); validationCollection.push_back(&Mu1Numberofvalidtrackerhits);
  Mu1Numberofvalidpixelhits=HConfig.GetTH1D(Name+"_Mu1Numberofvalidpixelhits","Mu1Numberofvalidpixelhits",20,-0.5,19.5,"n valid pixel hits",""); validationCollection.push_back(&Mu1Numberofvalidpixelhits);
  Mu1trackerLayersWithMeasurement=HConfig.GetTH1D(Name+"_Mu1trackerLayersWithMeasurement","Mu1trackerLayersWithMeasurement",15,-0.5,14.5,"#mu_{1} tracker layer with measurements",""); validationCollection.push_back(&Mu1trackerLayersWithMeasurement);
  Mu1validMuonHitComb=HConfig.GetTH1D(Name+"_Mu1validMuonHitComb","Mu1validMuonHitComb",70,-0.5,69.5,"valid hit combined",""); validationCollection.push_back(&Mu1validMuonHitComb);
  Mu1numberOfMatchedStations=HConfig.GetTH1D(Name+"_Mu1numberOfMatchedStations","Mu1numberOfMatchedStations",10,-0.5,9.5,"N matched stations",""); validationCollection.push_back(&Mu1numberOfMatchedStations);
  Mu1segmentCompatibility=HConfig.GetTH1D(Name+"_Mu1segmentCompatibility","Mu1segmentCompatibility",25,0,1,"#mu_{1} segment compatibility",""); validationCollection.push_back(&Mu1segmentCompatibility);
  Mu1timeAtIpInOutErr=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOutErr","Mu1timeAtIpInOutErr",25,0,3,"timeAtIpInOutErr",""); validationCollection.push_back(&Mu1timeAtIpInOutErr);
  Mu1timeAtIpInOut=HConfig.GetTH1D(Name+"_Mu1timeAtIpInOut","Mu1timeAtIpInOut",25,0,3,"timeAtIpInOut",""); validationCollection.push_back(&Mu1timeAtIpInOut);
  Mu1GLnormChi2=HConfig.GetTH1D(Name+"_Mu1GLnormChi2","Mu1GLnormChi2",25,0,2,"global normalised #chi^{2}",""); validationCollection.push_back(&Mu1GLnormChi2);
  Mu1innerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Mu1innerTrack_normalizedChi2","Mu1innerTrack_normalizedChi2",25,0,20,"inner track #chi^{2}",""); validationCollection.push_back(&Mu1innerTrack_normalizedChi2);
  Mu1outerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Mu1outerTrack_normalizedChi2","Mu1outerTrack_normalizedChi2",25,0,20,"outer track #chi^{2}",""); validationCollection.push_back(&Mu1outerTrack_normalizedChi2);
  Mu1innerTrack_validFraction=HConfig.GetTH1D(Name+"_Mu1innerTrack_validFraction","Mu1innerTrack_validFraction",25,0,2,"inner track valid Fraction",""); validationCollection.push_back(&Mu1innerTrack_validFraction);





  Mu2combinedQuality_chi2LocalMomentum_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalMomentum_peak","Mu2combinedQuality_chi2LocalMomentum_peak",25,0,20,"#chi^{2} local momentum",""); peakCollection.push_back(&Mu2combinedQuality_chi2LocalMomentum_peak);
  Mu2combinedQuality_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalPosition_peak","Mu2combinedQuality_chi2LocalMomentum_peak",25,0,20,"#chi^{2} local position",""); peakCollection.push_back(&Mu2combinedQuality_chi2LocalPosition_peak);
  Mu2combinedQuality_staRelChi2_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_staRelChi2_peak","Mu2combinedQuality_staRelChi2_peak",25,0,2,"station relative #chi^{2}",""); peakCollection.push_back(&Mu2combinedQuality_staRelChi2_peak);
  Mu2combinedQuality_trkRelChi2_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_trkRelChi2_peak","Mu2combinedQuality_trkRelChi2_peak",25,0,2,"track relative #chi^{2}",""); peakCollection.push_back(&Mu2combinedQuality_trkRelChi2_peak);
  Mu2combinedQuality_globalDeltaEtaPhi_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_globalDeltaEtaPhi_peak","Mu2combinedQuality_globalDeltaEtaPhi_peak",25,0,0.3,"gl #Delta#eta#phi",""); peakCollection.push_back(&Mu2combinedQuality_globalDeltaEtaPhi_peak);
  Mu2mu_combinedQuality_trkKink_peak=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_trkKink_peak","Mu2mu_combinedQuality_trkKink_peak",25,0,20,"trk kink",""); peakCollection.push_back(&Mu2mu_combinedQuality_trkKink_peak);
  Mu2mu_combinedQuality_glbKink_peak=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_glbKink_peak","Mu2mu_combinedQuality_glbKink_peak",25,0,5,"glb kink",""); peakCollection.push_back(&Mu2mu_combinedQuality_glbKink_peak);
  Mu2combinedQuality_glbTrackProbability_peak=HConfig.GetTH1D(Name+"_Mu2combinedQuality_glbTrackProbability_peak","Mu2combinedQuality_glbTrackProbability_peak",25,0,10,"",""); peakCollection.push_back(&Mu2combinedQuality_glbTrackProbability_peak);
  Mu2Numberofvalidtrackerhits_peak=HConfig.GetTH1D(Name+"_Mu2Numberofvalidtrackerhits_peak","Mu2Numberofvalidtrackerhits_peak",20,-0.5,19.5,"n valid tracker hits",""); peakCollection.push_back(&Mu2Numberofvalidtrackerhits_peak);
  Mu2Numberofvalidpixelhits_peak=HConfig.GetTH1D(Name+"_Mu2Numberofvalidpixelhits_peak","Mu2Numberofvalidpixelhits_peak",20,-0.5,19.5,"n valid pixel hits",""); peakCollection.push_back(&Mu2Numberofvalidpixelhits_peak);
  Mu2validMuonHitComb_peak=HConfig.GetTH1D(Name+"_Mu2validMuonHitComb_peak","Mu2validMuonHitComb_peak",70,-0.5,69.5,"valid hit combined",""); peakCollection.push_back(&Mu2validMuonHitComb_peak);
  Mu2numberOfMatchedStations_peak=HConfig.GetTH1D(Name+"_Mu2numberOfMatchedStations_peak","Mu2numberOfMatchedStations_peak",10,-0.5,9.5,"N matched stations",""); peakCollection.push_back(&Mu2numberOfMatchedStations_peak);
  Mu2segmentCompatibility_peak=HConfig.GetTH1D(Name+"_Mu2segmentCompatibility_peak","Mu2segmentCompatibility_peak",25,0,1,"#mu_{1} segment compatibility",""); peakCollection.push_back(&Mu2segmentCompatibility_peak);
  Mu2timeAtIpInOutErr_peak=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOutErr_peak","Mu2timeAtIpInOutErr_peak",25,0,3,"timeAtIpInOutErr",""); peakCollection.push_back(&Mu2timeAtIpInOutErr_peak);
  Mu2timeAtIpInOut_peak=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOut_peak","Mu2timeAtIpInOut_peak",25,0,3,"timeAtIpInOut",""); peakCollection.push_back(&Mu2timeAtIpInOut_peak);
  Mu2GLnormChi2_peak=HConfig.GetTH1D(Name+"_Mu2GLnormChi2_peak","Mu2GLnormChi2_peak",25,0,2,"global normalised #chi^{2}",""); peakCollection.push_back(&Mu2GLnormChi2_peak);
  Mu2innerTrack_normalizedChi2_peak=HConfig.GetTH1D(Name+"_Mu2innerTrack_normalizedChi2_peak","Mu2innerTrack_normalizedChi2_peak",25,0,20,"inner track #chi^{2}",""); peakCollection.push_back(&Mu2innerTrack_normalizedChi2_peak);
  Mu2outerTrack_normalizedChi2_peak=HConfig.GetTH1D(Name+"_Mu2outerTrack_normalizedChi2_peak","Mu2outerTrack_normalizedChi2_peak",25,0,20,"outer track #chi^{2}",""); peakCollection.push_back(&Mu2outerTrack_normalizedChi2_peak);
  Mu2innerTrack_validFraction_peak=HConfig.GetTH1D(Name+"_Mu2innerTrack_validFraction_peak","Mu2innerTrack_validFraction_peak",25,0,2,"inner track valid Fraction",""); peakCollection.push_back(&Mu2innerTrack_validFraction_peak);
  
  Mu2combinedQuality_chi2LocalMomentum_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalMomentum_sideband","Mu2combinedQuality_chi2LocalMomentum_sideband",25,0,20,"#chi^{2} local momentum",""); sidebandCollection.push_back(&Mu2combinedQuality_chi2LocalMomentum_sideband);
  Mu2combinedQuality_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalPosition_sideband","Mu2combinedQuality_chi2LocalPosition_sideband",25,0,20,"#chi^{2} local position",""); sidebandCollection.push_back(&Mu2combinedQuality_chi2LocalPosition_sideband);
  Mu2combinedQuality_staRelChi2_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_staRelChi2_sideband","Mu2combinedQuality_staRelChi2_sideband",25,0,2,"station relative #chi^{2}",""); sidebandCollection.push_back(&Mu2combinedQuality_staRelChi2_sideband);
  Mu2combinedQuality_trkRelChi2_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_trkRelChi2_sideband","Mu2combinedQuality_trkRelChi2_sideband",25,0,2,"track relative #chi^{2}",""); sidebandCollection.push_back(&Mu2combinedQuality_trkRelChi2_sideband);
  Mu2combinedQuality_globalDeltaEtaPhi_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_globalDeltaEtaPhi_sideband","Mu2combinedQuality_globalDeltaEtaPhi_sideband",25,0,0.3,"gl #Delta#eta#phi",""); sidebandCollection.push_back(&Mu2combinedQuality_globalDeltaEtaPhi_sideband);
  Mu2mu_combinedQuality_trkKink_sideband=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_trkKink_sideband","mMu2u_combinedQuality_trkKink_sideband",25,0,20,"trk kink",""); sidebandCollection.push_back(&Mu2mu_combinedQuality_trkKink_sideband);
  Mu2mu_combinedQuality_glbKink_sideband=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_glbKink_sideband","Mu2mu_combinedQuality_glbKink_sideband",25,0,5,"glb kink",""); sidebandCollection.push_back(&Mu2mu_combinedQuality_glbKink_sideband);
  Mu2combinedQuality_glbTrackProbability_sideband=HConfig.GetTH1D(Name+"_Mu2combinedQuality_glbTrackProbability_sideband","Mu2combinedQuality_glbTrackProbability_sideband",25,0,10,"",""); sidebandCollection.push_back(&Mu2combinedQuality_glbTrackProbability_sideband);
  Mu2Numberofvalidtrackerhits_sideband=HConfig.GetTH1D(Name+"_Mu2Numberofvalidtrackerhits_sideband","Mu2Numberofvalidtrackerhits_sideband",20,-0.5,19.5,"n valid tracker hits",""); sidebandCollection.push_back(&Mu2Numberofvalidtrackerhits_sideband);
  Mu2Numberofvalidpixelhits_sideband=HConfig.GetTH1D(Name+"_Mu2Numberofvalidpixelhits_sideband","Mu2Numberofvalidpixelhits_sideband",20,-0.5,19.5,"n valid pixel hits",""); sidebandCollection.push_back(&Mu2Numberofvalidpixelhits_sideband);
  Mu2validMuonHitComb_sideband=HConfig.GetTH1D(Name+"_Mu2validMuonHitComb_sideband","Mu2validMuonHitComb_sideband",70,-0.5,69.5,"valid hit combined",""); sidebandCollection.push_back(&Mu2validMuonHitComb_sideband);
  Mu2numberOfMatchedStations_sideband=HConfig.GetTH1D(Name+"_Mu2numberOfMatchedStations_sideband","Mu2numberOfMatchedStations_sideband",10,-0.5,9.5,"N matched stations",""); sidebandCollection.push_back(&Mu2numberOfMatchedStations_sideband);
  Mu2segmentCompatibility_sideband=HConfig.GetTH1D(Name+"_Mu2segmentCompatibility_sideband","Mu2segmentCompatibility_sideband",25,0,1,"#mu_{1} segment compatibility",""); sidebandCollection.push_back(&Mu2segmentCompatibility_sideband);
  Mu2timeAtIpInOutErr_sideband=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOutErr_sideband","Mu2timeAtIpInOutErr_sideband",25,0,3,"timeAtIpInOutErr",""); sidebandCollection.push_back(&Mu2timeAtIpInOutErr_sideband);
  Mu2timeAtIpInOut_sideband=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOut_sideband","Mu2timeAtIpInOut_sideband",25,0,3,"timeAtIpInOut",""); sidebandCollection.push_back(&Mu2timeAtIpInOut_sideband);
  Mu2GLnormChi2_sideband=HConfig.GetTH1D(Name+"_Mu2GLnormChi2_sideband","Mu2GLnormChi2_sideband",25,0,2,"global normalised #chi^{2}",""); sidebandCollection.push_back(&Mu2GLnormChi2_sideband);
  Mu2innerTrack_normalizedChi2_sideband=HConfig.GetTH1D(Name+"_Mu2innerTrack_normalizedChi2_sideband","Mu2innerTrack_normalizedChi2_sideband",25,0,20,"inner track #chi^{2}",""); sidebandCollection.push_back(&Mu2innerTrack_normalizedChi2_sideband);
  Mu2outerTrack_normalizedChi2_sideband=HConfig.GetTH1D(Name+"_Mu2outerTrack_normalizedChi2_sideband","Mu2outerTrack_normalizedChi2_sideband",25,0,20,"outer track #chi^{2}",""); sidebandCollection.push_back(&Mu2outerTrack_normalizedChi2_sideband);
  Mu2innerTrack_validFraction_sideband=HConfig.GetTH1D(Name+"_Mu2innerTrack_validFraction_sideband","Mu2innerTrack_validFraction_sideband",25,0,2,"inner track valid Fraction",""); sidebandCollection.push_back(&Mu2innerTrack_validFraction_sideband);
  


  Mu2combinedQuality_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalMomentum","Mu2combinedQuality_chi2LocalMomentum",25,0,20,"#chi^{2} local momentum",""); validationCollection.push_back(&Mu2combinedQuality_chi2LocalMomentum);
  Mu2combinedQuality_chi2LocalPosition=HConfig.GetTH1D(Name+"_Mu2combinedQuality_chi2LocalPosition","Mu2combinedQuality_chi2LocalPosition",25,0,20,"#chi^{2} local position",""); validationCollection.push_back(&Mu2combinedQuality_chi2LocalPosition);
  Mu2combinedQuality_staRelChi2=HConfig.GetTH1D(Name+"_Mu2combinedQuality_staRelChi2","Mu2combinedQuality_staRelChi2",25,0,2,"station relative #chi^{2}",""); validationCollection.push_back(&Mu2combinedQuality_staRelChi2);
  Mu2combinedQuality_trkRelChi2=HConfig.GetTH1D(Name+"_Mu2combinedQuality_trkRelChi2","Mu2combinedQuality_trkRelChi2",25,0,2,"track relative #chi^{2}",""); validationCollection.push_back(&Mu2combinedQuality_trkRelChi2);
  Mu2combinedQuality_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_Mu2combinedQuality_globalDeltaEtaPhi","Mu2combinedQuality_globalDeltaEtaPhi",25,0,0.3,"gl #Delta#eta#phi",""); validationCollection.push_back(&Mu2combinedQuality_globalDeltaEtaPhi);
  Mu2mu_combinedQuality_trkKink=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_trkKink","Mu2mu_combinedQuality_trkKink",25,0,20,"trk kink",""); validationCollection.push_back(&Mu2mu_combinedQuality_trkKink);
  Mu2mu_combinedQuality_glbKink=HConfig.GetTH1D(Name+"_Mu2mu_combinedQuality_glbKink","Mu2mu_combinedQuality_glbKink",25,0,5,"glb kink",""); validationCollection.push_back(&Mu2mu_combinedQuality_glbKink);
  Mu2combinedQuality_glbTrackProbability=HConfig.GetTH1D(Name+"_Mu2combinedQuality_glbTrackProbability","Mu2combinedQuality_glbTrackProbability",25,0,10,"",""); validationCollection.push_back(&Mu2combinedQuality_glbTrackProbability);
  Mu2Numberofvalidtrackerhits=HConfig.GetTH1D(Name+"_Mu2Numberofvalidtrackerhits","Mu2Numberofvalidtrackerhits",20,-0.5,19.5,"n valid tracker hits",""); validationCollection.push_back(&Mu2Numberofvalidtrackerhits);
  Mu2Numberofvalidpixelhits=HConfig.GetTH1D(Name+"_Mu2Numberofvalidpixelhits","Mu2Numberofvalidpixelhits",20,-0.5,19.5,"n valid pixel hits",""); validationCollection.push_back(&Mu2Numberofvalidpixelhits);
  Mu2validMuonHitComb=HConfig.GetTH1D(Name+"_Mu2validMuonHitComb","Mu2validMuonHitComb",70,-0.5,69.5,"valid hit combined",""); validationCollection.push_back(&Mu2validMuonHitComb);
  Mu2numberOfMatchedStations=HConfig.GetTH1D(Name+"_Mu2numberOfMatchedStations","Mu2numberOfMatchedStations",10,-0.5,9.5,"N matched stations",""); validationCollection.push_back(&Mu2numberOfMatchedStations);
  Mu2segmentCompatibility=HConfig.GetTH1D(Name+"_Mu2segmentCompatibility","Mu2segmentCompatibility",25,0,1,"#mu_{1} segment compatibility",""); validationCollection.push_back(&Mu2segmentCompatibility);
  Mu2timeAtIpInOutErr=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOutErr","Mu2timeAtIpInOutErr",25,0,3,"timeAtIpInOutErr",""); validationCollection.push_back(&Mu2timeAtIpInOutErr);
  Mu2timeAtIpInOut=HConfig.GetTH1D(Name+"_Mu2timeAtIpInOut","Mu2timeAtIpInOut",25,0,3,"timeAtIpInOut",""); validationCollection.push_back(&Mu2timeAtIpInOut);
  Mu2GLnormChi2=HConfig.GetTH1D(Name+"_Mu2GLnormChi2","Mu2GLnormChi2",25,0,2,"global normalised #chi^{2}",""); validationCollection.push_back(&Mu2GLnormChi2);
  Mu2innerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Mu2innerTrack_normalizedChi2","Mu2innerTrack_normalizedChi2",25,0,20,"inner track #chi^{2}",""); validationCollection.push_back(&Mu2innerTrack_normalizedChi2);
  Mu2outerTrack_normalizedChi2=HConfig.GetTH1D(Name+"_Mu2outerTrack_normalizedChi2","Mu2outerTrack_normalizedChi2",25,0,20,"outer track #chi^{2}",""); validationCollection.push_back(&Mu2outerTrack_normalizedChi2);
  Mu2innerTrack_validFraction=HConfig.GetTH1D(Name+"_Mu2innerTrack_validFraction","Mu2innerTrack_validFraction",25,0,2,"inner Track valid Fraction",""); validationCollection.push_back(&Mu2innerTrack_validFraction);



  VertexKFChi2_peak=HConfig.GetTH1D(Name+"VertexKFChi2_peak","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events"); peakCollection.push_back(&VertexKFChi2_peak);
  VertexKFChi2_sideband=HConfig.GetTH1D(Name+"VertexKFChi2_sideband","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events"); sidebandCollection.push_back(&VertexKFChi2_sideband);
  VertexKFChi2=HConfig.GetTH1D(Name+"VertexKFChi2","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events"); validationCollection.push_back(&VertexKFChi2);

  SVPVDsDirAngle_peak=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_peak","SVPVDsDirAngle",25,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); peakCollection.push_back(&SVPVDsDirAngle_peak);
  SVPVDsDirAngle_sideband=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_sideband","SVPVDsDirAngle",25,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); sidebandCollection.push_back(&SVPVDsDirAngle_sideband);
  SVPVDsDirAngle=HConfig.GetTH1D(Name+"_SVPVDsDirAngle","SVPVDsDirAngle",25,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); validationCollection.push_back(&SVPVDsDirAngle);

  NtracksClose_peak=HConfig.GetTH1D(Name+"_NtracksClose_peak","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","Events"); peakCollection.push_back(&NtracksClose_peak);
  NtracksClose_sideband=HConfig.GetTH1D(Name+"_NtracksClose_sideband","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","Events"); sidebandCollection.push_back(&NtracksClose_sideband);
  NtracksClose=HConfig.GetTH1D(Name+"_NtracksClose","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","Events"); validationCollection.push_back(&NtracksClose);

  NSV_peak=HConfig.GetTH1D(Name+"_NSV_peak","NSV",8,-0.5,7.5,"N vertices in the Ds cone","Events"); peakCollection.push_back(&NSV_peak);
  NSV_sideband=HConfig.GetTH1D(Name+"_NSV_sideband","NSV",8,-0.5,7.5,"N vertices in the Ds cone","Events"); sidebandCollection.push_back(&NSV_sideband);
  NSV=HConfig.GetTH1D(Name+"_NSV","NSV",8,-0.5,7.5,"N vertices in the Ds cone","Events"); validationCollection.push_back(&NSV);

  MinMuon_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_peak","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","Events"); peakCollection.push_back(&MinMuon_chi2LocalPosition_peak);
  MinMuon_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_sideband","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","Events"); sidebandCollection.push_back(&MinMuon_chi2LocalPosition_sideband);
  MinMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","Events"); validationCollection.push_back(&MinMuon_chi2LocalPosition);


  Mu1Muon_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_Mu1Muon_chi2LocalPosition_peak","Mu1Muon_chi2LocalPosition",50,0,5,"Mu1 Inner/Outer track #chi^{2}","Events"); peakCollection.push_back(&Mu1Muon_chi2LocalPosition_peak);
  Mu1Muon_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_Mu1Muon_chi2LocalPosition_sideband","Mu1Muon_chi2LocalPosition",50,0,5,"Mu1 Inner/Outer track #chi^{2}","Events"); sidebandCollection.push_back(&Mu1Muon_chi2LocalPosition_sideband);
  Mu1Muon_chi2LocalPosition=HConfig.GetTH1D(Name+"_Mu1Muon_chi2LocalPosition","Mu1Muon_chi2LocalPosition",50,0,5,"Mu1 Inner/Outer track #chi^{2}","Events"); validationCollection.push_back(&Mu1Muon_chi2LocalPosition);


  Mu2Muon_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_Mu2Muon_chi2LocalPosition_peak","Mu2Muon_chi2LocalPosition",50,0,5,"Mu2 Inner/Outer track #chi^{2}","Events"); peakCollection.push_back(&Mu2Muon_chi2LocalPosition_peak);
  Mu2Muon_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_Mu2Muon_chi2LocalPosition_sideband","Mu2Muon_chi2LocalPosition",50,0,5,"Mu2 Inner/Outer track #chi^{2}","Events"); sidebandCollection.push_back(&Mu2Muon_chi2LocalPosition_sideband);
  Mu2Muon_chi2LocalPosition=HConfig.GetTH1D(Name+"_Mu2Muon_chi2LocalPosition","Mu2Muon_chi2LocalPosition",50,0,5,"Mu2 Inner/Outer track #chi^{2}","Events"); validationCollection.push_back(&Mu2Muon_chi2LocalPosition);
 
  MindcaTrackSV_peak=HConfig.GetTH1D(Name+"_MindcaTrackSV_peak","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","Events"); peakCollection.push_back(&MindcaTrackSV_peak);
  MindcaTrackSV_sideband=HConfig.GetTH1D(Name+"_MindcaTrackSV_sideband","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","Events"); sidebandCollection.push_back(&MindcaTrackSV_sideband);
  MindcaTrackSV=HConfig.GetTH1D(Name+"_MindcaTrackSV","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","Events"); validationCollection.push_back(&MindcaTrackSV); 

  MinDca_peak=HConfig.GetTH1D(Name+"_MinDca_peak","MinDca",50,0,0.02,"Min distance between muons and track","Events"); peakCollection.push_back(&MinDca_peak);
  MinDca_sideband=HConfig.GetTH1D(Name+"_MinDca_sideband","MinDca",50,0,0.02,"Min distance between muons and track","Events"); sidebandCollection.push_back(&MinDca_sideband);
  MinDca=HConfig.GetTH1D(Name+"_MinDca","MinDca",50,0,0.02,"Min distance between muons and track","Events"); validationCollection.push_back(&MinDca);

  MinD0SigSV_peak=HConfig.GetTH1D(Name+"_MinD0SigSV_peak","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","Events"); peakCollection.push_back(&MinD0SigSV_peak);
  MinD0SigSV_sideband=HConfig.GetTH1D(Name+"_MinD0SigSV_sideband","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","Events"); sidebandCollection.push_back(&MinD0SigSV_sideband);
  MinD0SigSV=HConfig.GetTH1D(Name+"_MinD0SigSV","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","Events"); validationCollection.push_back(&MinD0SigSV);

  MinD0SigPV_peak=HConfig.GetTH1D(Name+"_MinD0SigPV_peak","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","Events"); peakCollection.push_back(&MinD0SigPV_peak);
  MinD0SigPV_sideband=HConfig.GetTH1D(Name+"_MinD0SigPV_sideband","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","Events"); sidebandCollection.push_back(&MinD0SigPV_sideband);
  MinD0SigPV=HConfig.GetTH1D(Name+"_MinD0SigPV","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","Events"); validationCollection.push_back(&MinD0SigPV);

  MaxVertexPairQuality_peak=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_peak","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events"); peakCollection.push_back(&MaxVertexPairQuality_peak);
  MaxVertexPairQuality_sideband=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_sideband","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events"); sidebandCollection.push_back(&MaxVertexPairQuality_sideband);
  MaxVertexPairQuality=HConfig.GetTH1D(Name+"_MaxVertexPairQuality","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events"); validationCollection.push_back(&MaxVertexPairQuality);

  MaxMuDistanceToSV_peak=HConfig.GetTH1D(Name+"_MaxMuDistanceToSV_peak","MaxMuDistanceToSV_peak",50,0,0.03,"#mu - SV max  distance, cm","Events");peakCollection.push_back(&MaxMuDistanceToSV_peak);
  MinMuDistanceToSV_peak=HConfig.GetTH1D(Name+"_MinMuDistanceToSV_peak","MinMuDistanceToSV_peak",50,0,0.03,"#mu - SV min  distance, cm","Events");peakCollection.push_back(&MinMuDistanceToSV_peak);

  MaxMuDistanceToSV_sideband=HConfig.GetTH1D(Name+"_MaxMuDistanceToSV_sideband","MaxMuDistanceToSV_sideband",50,0,0.03,"#mu - SV max  distance, cm","Events");sidebandCollection.push_back(&MaxMuDistanceToSV_sideband);
  MinMuDistanceToSV_sideband=HConfig.GetTH1D(Name+"_MinMuDistanceToSV_sideband","MinMuDistanceToSV_sideband",50,0,0.03,"#mu - SV min  distance, cm","Events");sidebandCollection.push_back(&MinMuDistanceToSV_sideband);

  MaxMuDistanceToSV=HConfig.GetTH1D(Name+"_MaxMuDistanceToSV","MaxMuDistanceToSV",50,0,0.03,"#mu - SV max  distance, cm","Events");validationCollection.push_back(&MaxMuDistanceToSV);
  MinMuDistanceToSV=HConfig.GetTH1D(Name+"_MinMuDistanceToSV","MinMuDistanceToSV",50,0,0.03,"#mu - SV min  distance, cm","Events");validationCollection.push_back(&MinMuDistanceToSV);

  /*  MinMuonImpacAngle_peak =HConfig.GetTH1D(Name+"_MinMuonImpacAngle_peak","MinMuonImpacAngle_peak",30,-1,1,"Min #mu impact angle","");;peakCollection.push_back(&MinMuonImpacAngle_peak);
  MaxMuonImpacAngle_peak =HConfig.GetTH1D(Name+"_MaxMuonImpacAngle_peak","MaxMuonImpacAngle_peak",30,-1,1,"Min #mu impact angle","");;peakCollection.push_back(&MaxMuonImpacAngle_peak);

  MinMuonImpacAngle_sideband =HConfig.GetTH1D(Name+"_MinMuonImpacAngle_sideband","MinMuonImpacAngle_sideband",30,-1,1,"Min #mu impact angle","");sidebandCollection.push_back(&MinMuonImpacAngle_sideband);
  MaxMuonImpacAngle_sideband =HConfig.GetTH1D(Name+"_MaxMuonImpacAngle_sideband","MaxMuonImpacAngle_sideband",30,-1,1,"Min #mu impact angle","");sidebandCollection.push_back(&MaxMuonImpacAngle_sideband);

  MinMuonImpacAngle =HConfig.GetTH1D(Name+"_MinMuonImpacAngle","MinMuonImpacAngle",30,-1,1,"Min #mu impact angle","");validationCollection.push_back(&MinMuonImpacAngle);
  MaxMuonImpacAngle =HConfig.GetTH1D(Name+"_MaxMuonImpacAngle","MaxMuonImpacAngle",30,-1,1,"Min #mu impact angle","");validationCollection.push_back(&MaxMuonImpacAngle);
  */


  Muon1ImpactAngle_peak =HConfig.GetTH1D(Name+"_Muon1ImpactAngle_peak","Muon1ImpactAngle_peak",30,-1,1,"#mu_{1} impact angle",""); peakCollection.push_back(&Muon1ImpactAngle_peak);
  Muon2ImpactAngle_peak =HConfig.GetTH1D(Name+"_Muon2ImpactAngle_peak","Muon2ImpactAngle_peak",30,-1,1,"#mu_{2} impact angle",""); peakCollection.push_back(&Muon2ImpactAngle_peak);

  Muon1ImpactAngle_sideband =HConfig.GetTH1D(Name+"_Muon1ImpactAngle_sideband","Muon1ImpactAngle_sideband",30,-1,1,"#mu_{1} impact angle","");sidebandCollection.push_back(&Muon1ImpactAngle_sideband);
  Muon2ImpactAngle_sideband =HConfig.GetTH1D(Name+"_Muon2ImpactAngle_sideband","Muon2ImpactAngle_sideband",30,-1,1,"#mu_{2} impact angle","");sidebandCollection.push_back(&Muon2ImpactAngle_sideband);

  Muon1ImpactAngle =HConfig.GetTH1D(Name+"_Muon1ImpactAngle","Muon1ImpactAngle",30,-1,1,"#mu_{1} impact angle","");validationCollection.push_back(&Muon1ImpactAngle);
  Muon2ImpactAngle =HConfig.GetTH1D(Name+"_Muon2ImpactAngle","Muon2ImpactAngle",30,-1,1,"#mu_{2} impact angle","");validationCollection.push_back(&Muon2ImpactAngle);

  MaxdeltaMuZ_peak = HConfig.GetTH1D(Name+"_MaxdeltaMuZ_peak","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","Events"); peakCollection.push_back(&MaxdeltaMuZ_peak);
  MaxdeltaMuZ_sideband = HConfig.GetTH1D(Name+"_MaxdeltaMuZ_sideband","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","Events"); sidebandCollection.push_back(&MaxdeltaMuZ_sideband);
  MaxdeltaMuZ = HConfig.GetTH1D(Name+"_MaxdeltaMuZ","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","Events"); validationCollection.push_back(&MaxdeltaMuZ);

  MaxDca_peak=HConfig.GetTH1D(Name+"_MaxDca_peak","MaxDca",20,0,0.10,"Max distance between muons","Events"); peakCollection.push_back(&MaxDca_peak);
  MaxDca_sideband=HConfig.GetTH1D(Name+"_MaxDca_sideband","MaxDca",20,0,0.10,"Max distance between muons","Events"); sidebandCollection.push_back(&MaxDca_sideband);
  MaxDca=HConfig.GetTH1D(Name+"_MaxDca","MaxDca",20,0,0.10,"Max distance between muons","Events"); validationCollection.push_back(&MaxDca);

  MaxD0SigSV_peak=HConfig.GetTH1D(Name+"_MaxD0SigSV_peak","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","Events"); peakCollection.push_back(&MaxD0SigSV_peak);
  MaxD0SigSV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigSV_sideband","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","Events"); sidebandCollection.push_back(&MaxD0SigSV_sideband);
  MaxD0SigSV=HConfig.GetTH1D(Name+"_MaxD0SigSV","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","Events"); validationCollection.push_back(&MaxD0SigSV);

  MaxD0SigPV_peak=HConfig.GetTH1D(Name+"_MaxD0SigPV_peak","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","Events"); peakCollection.push_back(&MaxD0SigPV_peak);
  MaxD0SigPV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigPV_sideband","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","Events"); sidebandCollection.push_back(&MaxD0SigPV_sideband);
  MaxD0SigPV=HConfig.GetTH1D(Name+"_MaxD0SigPV","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","Events"); validationCollection.push_back(&MaxD0SigPV);

  MaxD0SigBS_peak=HConfig.GetTH1D(Name+"_MaxD0SigBS_peak","MaxD0SigBS",50,0,15,"Max Transverse Impact significance w.r.t BS",""); peakCollection.push_back(&MaxD0SigBS_peak);
  MaxD0SigBS_sideband=HConfig.GetTH1D(Name+"_MaxD0SigBS_sideband","MaxD0SigBS",50,0,15,"Max Transverse Impact significance w.r.t BS","");sidebandCollection.push_back(&MaxD0SigBS_sideband);
  MaxD0SigBS=HConfig.GetTH1D(Name+"_MaxD0SigBS","MaxD0SigBS",50,0,15,"Max Transverse Impact significance w.r.t BS",""); validationCollection.push_back(&MaxD0SigBS);


  Iso1_peak=HConfig.GetTH1D(Name+"_Iso1_peak","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0"); peakCollection.push_back(&Iso1_peak);
  Iso1_sideband=HConfig.GetTH1D(Name+"_Iso1_sideband","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0"); sidebandCollection.push_back(&Iso1_sideband);
  Iso1=HConfig.GetTH1D(Name+"_Iso1","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0"); validationCollection.push_back(&Iso1);

  Iso1Mu1_peak=HConfig.GetTH1D(Name+"_Iso1Mu1_peak","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0"); peakCollection.push_back(&Iso1Mu1_peak);
  Iso1Mu1_sideband=HConfig.GetTH1D(Name+"_Iso1Mu1_sideband","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0"); sidebandCollection.push_back(&Iso1Mu1_sideband);
  Iso1Mu1=HConfig.GetTH1D(Name+"_Iso1Mu1","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0"); validationCollection.push_back(&Iso1Mu1);

  Iso8Mu1_peak=HConfig.GetTH1D(Name+"_Iso8Mu1_peak","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8"); peakCollection.push_back(&Iso8Mu1_peak);
  Iso8Mu1_sideband=HConfig.GetTH1D(Name+"_Iso8Mu1_sideband","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8"); sidebandCollection.push_back(&Iso8Mu1_sideband);
  Iso8Mu1=HConfig.GetTH1D(Name+"_Iso8Mu1","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8"); validationCollection.push_back(&Iso8Mu1);

  FLSignificance_peak=HConfig.GetTH1D(Name+"_FLSignificance_peak","FLSignificance",60,0,60,"PV - SV distance  significance","Events"); peakCollection.push_back(&FLSignificance_peak);
  FLSignificance_sideband=HConfig.GetTH1D(Name+"_FLSignificance_sideband","FLSignificance",60,0,60,"PV - SV distance  significance","Events"); sidebandCollection.push_back(&FLSignificance_sideband);
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",60,0,60,"PV - SV distance  significance","Events"); validationCollection.push_back(&FLSignificance);

  MatchedMuonTrackdPT=HConfig.GetTH1D(Name+"_MatchedMuonTrackdPT","MatchedMuonTrackdPT",100,0,0.5,"|p_{T}(#mu) - p_{T}(track)|, GeV","Events");



  TrkIDTrackmuonInnerNC2_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNC2_peak","TrkIDTrackmuonInnerNC2_peak",50,0,10,"track normalized #chi^{2}",""); peakCollection.push_back(&TrkIDTrackmuonInnerNValidHits_peak);
  TrkIDTrackmuonInnerNValidHits_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNValidHits_peak","TrkIDTrackmuonInnerNValidHits_peak",30,-0.5,29.5,"valid hits in inner tracker",""); peakCollection.push_back(&TrkIDTrackmuonInnerNValidHits_peak);
  TrkIDTrackmuonValidFraction_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonValidFraction_peak","TrkIDTrackmuonValidFraction_peak",50,0.5,2,"TrkMu ID, valid fraction",""); peakCollection.push_back(&TrkIDTrackmuonValidFraction_peak);
  TrkIDTrackmuonNLostTrackerHits_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHits_peak","TrkIDTrackmuonNLostTrackerHits_peak",5,-0.5,4.5,"lost tracker hits",""); peakCollection.push_back(&TrkIDTrackmuonNLostTrackerHits_peak);
 
 TrkIDTrackmuonNLostTrackerHitsInner_peak   =HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsInner_peak","TrkIDTrackmuonNLostTrackerHitsInner_peak",15,-0.5,14.5,"n lost tracker inner tracker hits",""); peakCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsInner_peak);
  TrkIDTrackmuonNLostTrackerHitsOuter_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsOuter_peak","TrkIDTrackmuonNLostTrackerHitsOuter_peak",10,-0.5,9.5,"lost outer tracker hits",""); peakCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsOuter_peak);
  TrkIDTrackmuonPixelLayers_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPixelLayers_peak","TrkIDTrackmuonPixelLayers_peak",10,-0.5,9.5,"N pixel layers",""); peakCollection.push_back(&TrkIDTrackmuonPixelLayers_peak);
  TrkIDTrackmuonNMatchedStations_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNMatchedStations_peak","TrkIDTrackmuonNMatchedStations_peak",10,-0.5,9.5,"N matched stations",""); peakCollection.push_back(&TrkIDTrackmuonNMatchedStations_peak);
  TrkIDTrackmuonPtErrPt_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPtErrPt_peak","TrkIDTrackmuonPtErrPt_peak",20,0,0.05,"#delta p_{T}/p_{T}",""); peakCollection.push_back(&TrkIDTrackmuonPtErrPt_peak);
  TrkIDTrackmuonCaloComp_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonCaloComp_peak","TrkIDTrackmuonCaloComp_peak",50,0,1,"calorimeter compatibility",""); peakCollection.push_back(&TrkIDTrackmuonCaloComp_peak);
  TrkIDTrackmuonSegComp_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonSegComp_peak","TrkIDTrackmuonSegComp_peak",50,0,1,"segment compatibility",""); peakCollection.push_back(&TrkIDTrackmuonSegComp_peak);
  TrkIDTrackmuonHad_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonHad_peak","TrkIDTrackmuonHad_peak",50,0,15,"HCAL responce, GeV",""); peakCollection.push_back(&TrkIDTrackmuonHad_peak);
  TrkIDTrackmuonEM_peak=HConfig.GetTH1D(Name+"_TrkIDTrackmuonEM_peak","TrkIDTrackmuonEM_peak",50,0,4,"ECAL responce, GeV",""); peakCollection.push_back(&TrkIDTrackmuonEM_peak);




  TrkIDTrackmuonInnerNC2_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNC2_sideband","TrkIDTrackmuonInnerNC2_sideband",50,0,10,"track normalized #chi^{2}",""); sidebandCollection.push_back(&TrkIDTrackmuonInnerNC2_sideband);
  TrkIDTrackmuonInnerNValidHits_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNValidHits_sideband","TrkIDTrackmuonInnerNValidHits_sideband",30,-0.5,29.5,"valid hits in inner tracker",""); sidebandCollection.push_back(&TrkIDTrackmuonInnerNValidHits_sideband);
  TrkIDTrackmuonValidFraction_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonValidFraction_sideband","TrkIDTrackmuonValidFraction_sideband",50,0.5,2,"TrkMu ID, valid fraction",""); sidebandCollection.push_back(&TrkIDTrackmuonValidFraction_sideband);
  TrkIDTrackmuonNLostTrackerHits_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHits_sideband","TrkIDTrackmuonNLostTrackerHits_sideband",5,-0.5,4.5,"lost tracker hits",""); sidebandCollection.push_back(&TrkIDTrackmuonNLostTrackerHits_sideband);
  TrkIDTrackmuonNLostTrackerHitsInner_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsInner_sideband","TrkIDTrackmuonNLostTrackerHitsInner_sideband",15,-0.5,14.5,"n lost tracker inner tracker hits",""); sidebandCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsInner_sideband);
  TrkIDTrackmuonNLostTrackerHitsOuter_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsOuter_sideband","TrkIDTrackmuonNLostTrackerHitsOuter_sideband",10,-0.5,9.5,"lost outer tracker hits",""); sidebandCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsOuter_sideband);
  TrkIDTrackmuonPixelLayers_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPixelLayers_sideband","TrkIDTrackmuonPixelLayers_sideband",10,-0.5,9.5,"N pixel layers",""); sidebandCollection.push_back(&TrkIDTrackmuonPixelLayers_sideband);
  TrkIDTrackmuonNMatchedStations_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNMatchedStations_sideband","TrkIDTrackmuonNMatchedStations_sideband",10,-0.5,9.5,"N matched stations",""); sidebandCollection.push_back(&TrkIDTrackmuonNMatchedStations_sideband);
  TrkIDTrackmuonPtErrPt_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPtErrPt_sideband","TrkIDTrackmuonPtErrPt_sideband",20,0,0.05,"#delta p_{T}/p_{T}",""); sidebandCollection.push_back(&TrkIDTrackmuonPtErrPt_sideband);
  TrkIDTrackmuonCaloComp_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonCaloComp_sideband","TrkIDTrackmuonCaloComp_sideband",50,0,1,"calorimeter compatibility",""); sidebandCollection.push_back(&TrkIDTrackmuonCaloComp_sideband);
  TrkIDTrackmuonSegComp_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonSegComp_sideband","TrkIDTrackmuonSegComp_sideband",50,0,1,"segment compatibility",""); sidebandCollection.push_back(&TrkIDTrackmuonSegComp_peak);
  TrkIDTrackmuonHad_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonHad_sideband","TrkIDTrackmuonHad_sideband",50,0,15,"HCAL responce, GeV",""); sidebandCollection.push_back(&TrkIDTrackmuonHad_sideband);
  TrkIDTrackmuonEM_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackmuonEM_sideband","TrkIDTrackmuonEM_sideband",50,0,4,"ECAL responce, GeV",""); sidebandCollection.push_back(&TrkIDTrackmuonEM_sideband);
  




  TrkIDTrackmuonInnerNC2=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNC2","TrkIDTrackmuonInnerNC2",50,0,10,"track normalized #chi^{2}",""); validationCollection.push_back(&TrkIDTrackmuonInnerNC2);
  TrkIDTrackmuonInnerNValidHits=HConfig.GetTH1D(Name+"_TrkIDTrackmuonInnerNValidHits","TrkIDTrackmuonInnerNValidHits",30,-0.5,29.5,"valid hits in inner tracker",""); validationCollection.push_back(&TrkIDTrackmuonInnerNValidHits);
  TrkIDTrackmuonValidFraction=HConfig.GetTH1D(Name+"_TrkIDTrackmuonValidFraction","TrkIDTrackmuonValidFraction",50,0.5,2,"TrkMu ID, valid fraction",""); validationCollection.push_back(&TrkIDTrackmuonValidFraction);
  TrkIDTrackmuonNLostTrackerHits=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHits","TrkIDTrackmuonNLostTrackerHits",5,-0.5,4.5,"lost tracker hits",""); validationCollection.push_back(&TrkIDTrackmuonNLostTrackerHits);
  TrkIDTrackmuonNLostTrackerHitsInner=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsInner","TrkIDTrackmuonNLostTrackerHitsInner",15,-0.5,14.5,"n lost tracker inner tracker hits",""); validationCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsInner);
  TrkIDTrackmuonNLostTrackerHitsOuter=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNLostTrackerHitsOuter","TrkIDTrackmuonNLostTrackerHitsOuter",10,-0.5,9.5,"lost outer tracker hits",""); validationCollection.push_back(&TrkIDTrackmuonNLostTrackerHitsOuter);
  TrkIDTrackmuonPixelLayers=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPixelLayers","TrkIDTrackmuonPixelLayers",10,-0.5,9.5,"N pixel layers",""); validationCollection.push_back(&TrkIDTrackmuonPixelLayers);
  TrkIDTrackmuonNMatchedStations=HConfig.GetTH1D(Name+"_TrkIDTrackmuonNMatchedStations","TrkIDTrackmuonNMatchedStations",10,-0.5,9.5,"N matched stations",""); validationCollection.push_back(&TrkIDTrackmuonNMatchedStations);
  TrkIDTrackmuonPtErrPt=HConfig.GetTH1D(Name+"_TrkIDTrackmuonPtErrPt","TrkIDTrackmuonPtErrPt",20,0,0.05,"#delta p_{T}/p_{T}",""); validationCollection.push_back(&TrkIDTrackmuonPtErrPt);
  TrkIDTrackmuonCaloComp=HConfig.GetTH1D(Name+"_TrkIDTrackmuonCaloComp","TrkIDTrackmuonCaloComp",50,0,1,"calorimeter compatibility",""); validationCollection.push_back(&TrkIDTrackmuonCaloComp);
  TrkIDTrackmuonSegComp=HConfig.GetTH1D(Name+"_TrkIDTrackmuonSegComp","TrkIDTrackmuonSegComp",50,0,1,"segment compatibility",""); validationCollection.push_back(&TrkIDTrackmuonSegComp_peak);
  TrkIDTrackmuonHad=HConfig.GetTH1D(Name+"_TrkIDTrackmuonHad","TrkIDTrackmuonHad",50,0,15,"HCAL responce, GeV",""); validationCollection.push_back(&TrkIDTrackmuonHad);
  TrkIDTrackmuonEM=HConfig.GetTH1D(Name+"_TrkIDTrackmuonEM","TrkIDTrackmuonEM",50,0,4,"ECAL responce, GeV",""); validationCollection.push_back(&TrkIDTrackmuonEM);
  


  TrkIDTrackBDT_peak=HConfig.GetTH1D(Name+"_TrkIDTrackBDT_peak","TrkIDTrackBDT_peak",30,-0.3,0.3,"BDT output",""); peakCollection.push_back(&TrkIDTrackBDT_peak);
  TrkIDTrackBDT_sideband=HConfig.GetTH1D(Name+"_TrkIDTrackBDT_sideband","TrkIDTrackBDT_sideband",30,-0.3,0.3,"BDT output",""); sidebandCollection.push_back(&TrkIDTrackBDT_sideband);
  TrkIDTrackBDT=HConfig.GetTH1D(Name+"_TrkIDTrackBDT","TrkIDTrackBDT",30,-0.3,0.3,"BDT output",""); validationCollection.push_back(&TrkIDTrackBDT);


  SVPosX_peak=HConfig.GetTH1D(Name+"_SVPosX_peak","SVPosX_peak",60,-1,1,"SV pos X, cm",""); peakCollection.push_back(&SVPosX_peak);
  SVPosX_sideband=HConfig.GetTH1D(Name+"_SVPosX_sideband","SVPosX_sideband",60,-1,1,"SV pos X, cm",""); sidebandCollection.push_back(&SVPosX_sideband);
  SVPosX=HConfig.GetTH1D(Name+"_SVPosX","SVPosX",60,-1,1,"SV pos X, cm",""); validationCollection.push_back(&SVPosX);

  SVPosY_peak=HConfig.GetTH1D(Name+"_SVPosY_peak","SVPosY_peak",60,-1,1,"SV pos Y, cm",""); peakCollection.push_back(&SVPosY_peak);
  SVPosY_sideband=HConfig.GetTH1D(Name+"_SVPosY_sideband","SVPosY_sideband",60,-1,1,"SV pos Y, cm",""); sidebandCollection.push_back(&SVPosY_sideband);
  SVPosY=HConfig.GetTH1D(Name+"_SVPosY","SVPosY",60,-1,1,"SV pos Y, cm",""); validationCollection.push_back(&SVPosY);

  SVPosZ_peak=HConfig.GetTH1D(Name+"_SVPosZ_peak","SVPosZ_peak",50,-15,15,"SV pos Z, cm",""); peakCollection.push_back(&SVPosZ_peak);
  SVPosZ_sideband=HConfig.GetTH1D(Name+"_SVPosZ_sideband","SVPosZ_sideband",50,-15,15,"SV pos Z, cm",""); sidebandCollection.push_back(&SVPosZ_sideband);
  SVPosZ=HConfig.GetTH1D(Name+"_SVPosZ","SVPosZ",50,-15,15,"SV pos Z, cm",""); validationCollection.push_back(&SVPosZ);

  /*  BSPosX_peak=HConfig.GetTH1D(Name+"_BSPosX_peak","BSPosX_peak",30,-2,2,"BS pos X, cm",""); peakCollection.push_back(&BSPosX_peak);
  BSPosX_sideband=HConfig.GetTH1D(Name+"_BSPosX_sideband","BSPosX_sideband",30,-2,2,"BS pos X, cm",""); sidebandCollection.push_back(&BSPosX_sideband);
  BSPosX=HConfig.GetTH1D(Name+"_BSPosX","BSPosX",30,-2,2,"BS pos X, cm",""); validationCollection.push_back(&BSPosX);

  BSPosY_peak=HConfig.GetTH1D(Name+"_BSPosY_peak","BSPosY_peak",30,-2,2,"BS pos Y, cm",""); peakCollection.push_back(&BSPosY_peak);
  BSPosY_sideband=HConfig.GetTH1D(Name+"_BSPosY_sideband","BSPosY_sideband",30,-2,2,"BS pos Y, cm",""); sidebandCollection.push_back(&BSPosY_sideband);
  BSPosY=HConfig.GetTH1D(Name+"_BSPosY","BSPosY",30,-2,2,"BS pos Y, cm",""); validationCollection.push_back(&BSPosY);*/


  PVPosX_peak=HConfig.GetTH1D(Name+"_PVPosX_peak","PVPosX_peak",50,-0.5,0.5,"PV pos X, cm",""); peakCollection.push_back(&PVPosX_peak);
  PVPosX_sideband=HConfig.GetTH1D(Name+"_PVPosX_sideband","PVPosX_sideband",50,-0.5,0.5,"PV pos X, cm",""); sidebandCollection.push_back(&PVPosX_sideband);
  PVPosX=HConfig.GetTH1D(Name+"_PVPosX","PVPosX",50,-0.5,0.5,"PV pos X, cm",""); validationCollection.push_back(&PVPosX);

  PVPosY_peak=HConfig.GetTH1D(Name+"_PVPosY_peak","PVPosY_peak",50,-0.5,0.5,"PV pos Y, cm",""); peakCollection.push_back(&PVPosY_peak);
  PVPosY_sideband=HConfig.GetTH1D(Name+"_PVPosY_sideband","PVPosY_sideband",50,-0.5,0.5,"PV pos Y, cm",""); sidebandCollection.push_back(&PVPosY_sideband);
  PVPosY=HConfig.GetTH1D(Name+"_PVPosY","PVPosY",50,-0.5,0.5,"PV pos Y, cm",""); validationCollection.push_back(&PVPosY);

  PVPosZ_peak=HConfig.GetTH1D(Name+"_PVPosZ_peak","PVPosZ_peak",50,-15,15,"PV pos Z, cm",""); peakCollection.push_back(&PVPosZ_peak);
  PVPosZ_sideband=HConfig.GetTH1D(Name+"_PVPosZ_sideband","PVPosZ_sideband",50,-15,15,"PV pos Z, cm",""); sidebandCollection.push_back(&PVPosZ_sideband);
  PVPosZ=HConfig.GetTH1D(Name+"_PVPosZ","PVPosZ",50,-15,15,"PV pos Z, cm",""); validationCollection.push_back(&PVPosZ);

  SVPVEta_peak=HConfig.GetTH1D(Name+"_SVPVEta_peak","SVPVEta_peak",40,-2.5,2.5,"#eta (#vec{SV}-#vec{PV})",""); peakCollection.push_back(&SVPVEta_peak);
  SVPVPhi_peak=HConfig.GetTH1D(Name+"_SVPVPhi_peak","SVPVPhi_peak",40,-3.14,3.14,"#phi (#vec{SV}-#vec{PV})",""); peakCollection.push_back(&SVPVPhi_peak);

  SVPVEta_sideband=HConfig.GetTH1D(Name+"_SVPVEta_sideband","SVPVEta_sideband",40,-2.5,2.5,"#eta (#vec{SV}-#vec{PV})",""); sidebandCollection.push_back(&SVPVEta_sideband);
  SVPVPhi_sideband=HConfig.GetTH1D(Name+"_SVPVPhi_sideband","SVPVPhi_sideband",40,-3.14,3.14,"#phi (#vec{SV}-#vec{PV})",""); sidebandCollection.push_back(&SVPVPhi_sideband);
  
  SVPVEta=HConfig.GetTH1D(Name+"_SVPVEta","SVPVEta",40,-2.5,2.5,"#eta (#vec{SV}-#vec{PV})",""); validationCollection.push_back(&SVPVEta);
  SVPVPhi=HConfig.GetTH1D(Name+"_SVPVPhi","SVPVPhi",40,-3.14,3.14,"#phi (#vec{SV}-#vec{PV})",""); validationCollection.push_back(&SVPVPhi);




  SVPVTransverse_peak=HConfig.GetTH1D(Name+"_SVPVTransverse_peak","SVPVTransverse_peak",40,0,1,"d_{T} (SV-PV), cm",""); peakCollection.push_back(&SVPVTransverse_peak);
  SVPVTransverse_sideband=HConfig.GetTH1D(Name+"_SVPVTransverse_sideband","SVPVTransverse_sideband",40,0,1,"d_{T} (SV-PV), cm",""); sidebandCollection.push_back(&SVPVTransverse_sideband);
  SVPVTransverse=HConfig.GetTH1D(Name+"_SVPVTransverse","SVPVTransverse",40,0,1,"d_{T} (SV-PV), cm",""); validationCollection.push_back(&SVPVTransverse);



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsPhiPeak_2018::Store_ExtraDist(){ 
  
  Extradist1d.push_back(&Track_normalizedChi2);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&PhiPlusTrackMass);
  Extradist2d.push_back(&PhiMassVsDsMass);
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&MatchedMuonTrackdPT);

  // Validation plots

  for (unsigned int j=0; j<peakCollection.size(); ++j){
    //    Extradist1d.push_back(peakCollection.at(j));
    //    Extradist1d.push_back(sidebandCollection.at(j));
    Extradist1d.push_back(validationCollection.at(j));
  }

  Extradist1d.push_back(&DecayLength_prompt);
  Extradist1d.push_back(&DecayLength_non_prompt);

}


void  DsPhiPeak_2018::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  //  std::cout<<"--- "<< std::endl;
  //bool RunB,   RunC, RunD, RunE, RunF = 0;

  if(id==1 && Ntp->WhichEra(2017).Contains("RunB") ){ RunB=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunC") ){ RunC=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunD") ){ RunD=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunE") ){ RunE=1;}
  if(id==1 && Ntp->WhichEra(2017).Contains("RunF") ){ RunF=1;}

  random_num = rndm.Rndm();


  value.at(L1TOk) = 0;
  value.at(HLTOk) = 0;

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) {value.at(HLTOk)=Ntp->HLTDecision(iTrigger);}
  }


  bool DoubleMuFired(false);
  //  bool TripleMuFired(false);
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);


    //    if( Ntp->L1Decision(il1)) std::cout<<" Trigger  "<< Ntp->L1Name(il1) <<std::endl;
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    //    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    //    if( id!=1 && random_num>0.3516 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    //    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
  }

  // if (DoubleMuFired || TripleMuFired) value.at(L1TOk)=1;; // Remove all the events seeded by triple mu
  if (DoubleMuFired) value.at(L1TOk)=1;; // Remove all the events seeded by triple mu




  int mu1=-1, mu2=-1, track=-1;
  int final_idx = -1;
  double minChiSq = 999.0;
  //double check_PhiMass = 999.0; 
  std::vector<unsigned int> selectedIndices;
  std::vector<unsigned int> candidateRank;

  value.at(is2MuTrk) = 0; 
  value.at(Mass2Mu) = 0;
  value.at(MuCharge) = 0;
  value.at(TriggerMatch) = 0;

  if(Ntp->NTwoMuonsTrack()!=0/* && Ntp->NThreeMuons() == 0*/){
    value.at(is2MuTrk) = 1;
  }

  pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  pass.at(L1TOk)= (value.at(L1TOk)==cut.at(L1TOk));
  pass.at(HLTOk)= (value.at(HLTOk)/*==cut.at(HLTOk)*/);



  if (value.at(is2MuTrk)==1){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);
      double mumutrkmass = (Ntp->Muon_P4(tmp_mu1)+Ntp->Muon_P4(tmp_mu2)+Ntp->Track_P4(tmp_track)).M();    



      std::vector<TLorentzVector> trigobjTriplet;
      for (int i=0; i<Ntp->NTriggerObjects(); i++){
	TString name = Ntp->TriggerObject_name(i);
	if (!(name.Contains("hltTau3muTkVertexFilter"))) continue;

	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
	trigobjTriplet.push_back(tmp);
      }


      std::vector<TLorentzVector> targetTriplet;
      targetTriplet.push_back(Ntp->Muon_P4(tmp_mu1));
      targetTriplet.push_back(Ntp->Muon_P4(tmp_mu2));
      targetTriplet.push_back(Ntp->Track_P4(tmp_track));

      bool triggerCheck = false;
      if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(targetTriplet, trigobjTriplet);
      value.at(TriggerMatch) = triggerCheck;


      value.at(GlobalMu) = (Ntp->Muon_isGlobalMuon(tmp_mu1)==1 && Ntp->Muon_isGlobalMuon(tmp_mu2)==1 && Ntp->Muon_isPFMuon(tmp_mu1) &&  Ntp->Muon_isPFMuon(tmp_mu2)) ;
      value.at(Mass2Mu) = (Ntp->Muon_P4(tmp_mu1) + Ntp->Muon_P4(tmp_mu2)).M();
      value.at(MuCharge) = Ntp->Muon_charge(tmp_mu1)!=Ntp->Muon_charge(tmp_mu2);
      value.at(Mu1pt) = Ntp->Muon_P4(tmp_mu1).Pt();
      value.at(Mu2pt) = Ntp->Muon_P4(tmp_mu2).Pt();
      value.at(Trkpt) = Ntp->Track_P4(tmp_track).Pt();
      value.at(DsMassCut) = mumutrkmass;

      pass.at(GlobalMu) = value.at(GlobalMu)==cut.at(GlobalMu);
      pass.at(Mass2Mu) = value.at(Mass2Mu) >= 0.99 && value.at(Mass2Mu) <= 1.048;
      pass.at(MuCharge) = value.at(MuCharge)==cut.at(MuCharge);

      pass.at(TriggerMatch) = (value.at(TriggerMatch) == cut.at(TriggerMatch));
      pass.at(Mu1pt) = (value.at(Mu1pt) >= cut.at(Mu1pt));
      pass.at(Mu2pt) = (value.at(Mu2pt) >= cut.at(Mu2pt));
      pass.at(Trkpt) = (value.at(Trkpt) >= cut.at(Trkpt)); 
      pass.at(DsMassCut) = value.at(DsMassCut) < dsMassMax && value.at(DsMassCut) >= dsMassMin;
 
      unsigned int score = 0;
      for (unsigned int k=0; k<NCuts; ++k) if (pass.at(k)) score++;

      if (score==NCuts) selectedIndices.push_back(i2M);
      candidateRank.push_back(score);

      //if (abs(tmp_PhiMass-1.01)<=check_PhiMass || (tmp_PhiMass > .95 && tmp_PhiMass < 1.1)) {
      //if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M)){
      //  tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
      //  check_PhiMass = abs(tmp_PhiMass-1.01);
      //  
      //  if (Ntp->Muon_P4(tmp_mu1).Pt() > Ntp->Muon_P4(tmp_mu2).Pt()) {
      //    mu1 = tmp_mu1;
      //    mu2 = tmp_mu2;
      //  } else {
      //    mu1 = tmp_mu2;
      //    mu2 = tmp_mu1;
      //  }
      //  track = tmp_track;
      //  final_idx = i2M;
      //}
      //}

      if (selectedIndices.size()==0) final_idx=std::distance(candidateRank.begin(), std::max_element(candidateRank.begin(), candidateRank.end()));
      else{
        for (size_t j=0; j<selectedIndices.size(); ++j){
          int _ = selectedIndices.at(j);
            double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(_);
            if (tmpchi<minChiSq) {
               minChiSq = tmpchi;
               final_idx = _;
            }
        }
      }
    }




    mu1 = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
    mu2 = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
    track = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);
    //double PhiMass = (Ntp->Muon_P4(mu1)+Ntp->Muon_P4(mu2)).M();
    double mumutrkmass = (Ntp->Muon_P4(mu1)+Ntp->Muon_P4(mu2)+Ntp->Track_P4(track)).M();


    std::vector<TLorentzVector> trigobjTriplet;
    for (int i=0; i<Ntp->NTriggerObjects(); i++){
      TString name = Ntp->TriggerObject_name(i);
      if (!(name.Contains("hltTau3muTkVertexFilter"))) continue;  // two mu plus track trigger objects
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
      trigobjTriplet.push_back(tmp);
    }


    std::vector<TLorentzVector> targetTriplet;
    targetTriplet.push_back(Ntp->Muon_P4(mu1));
    targetTriplet.push_back(Ntp->Muon_P4(mu2));
    targetTriplet.push_back(Ntp->Track_P4(track));

    bool triggerCheck = false;
    if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(targetTriplet, trigobjTriplet);
    value.at(TriggerMatch) = triggerCheck;


    value.at(GlobalMu) = (Ntp->Muon_isGlobalMuon(mu1)==1 && Ntp->Muon_isGlobalMuon(mu2)==1 && Ntp->Muon_isPFMuon(mu1) &&  Ntp->Muon_isPFMuon(mu2));
    value.at(Mass2Mu) = (Ntp->Muon_P4(mu1) + Ntp->Muon_P4(mu2)).M();
    value.at(MuCharge) = Ntp->Muon_charge(mu1)!=Ntp->Muon_charge(mu2);
    value.at(Mu1pt) = Ntp->Muon_P4(mu1).Pt();
    value.at(Mu2pt) = Ntp->Muon_P4(mu2).Pt();
    value.at(Trkpt) = Ntp->Track_P4(track).Pt();
    value.at(DsMassCut) = mumutrkmass;

    pass.at(GlobalMu) = value.at(GlobalMu)==cut.at(GlobalMu);
    pass.at(Mass2Mu) = value.at(Mass2Mu) >= 1 && value.at(Mass2Mu) <= 1.04;
    pass.at(MuCharge) = value.at(MuCharge)==cut.at(MuCharge);
    pass.at(Mu1pt) = (value.at(Mu1pt) >= cut.at(Mu1pt));
    pass.at(Mu2pt) = (value.at(Mu2pt) >= cut.at(Mu2pt));
    pass.at(Trkpt) = (value.at(Trkpt) >= cut.at(Trkpt));
    pass.at(TriggerMatch) = (value.at(TriggerMatch) == cut.at(TriggerMatch));
    pass.at(DsMassCut) = value.at(DsMassCut) < dsMassMax && value.at(DsMassCut) >= dsMassMin;

  }

  //pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  //pass.at(L1TOk)= (value.at(L1TOk)/*==cut.at(L1TOk)*/);
  //pass.at(HLTOk)= (value.at(HLTOk)/*==cut.at(HLTOk)*/);
  //pass.at(GlobalMu) = value.at(GlobalMu)==cut.at(GlobalMu);
  //pass.at(Mass2Mu) = value.at(Mass2Mu) >= 1 && value.at(Mass2Mu) <= 1.04;
  //pass.at(MuCharge) = value.at(MuCharge)==cut.at(MuCharge);
  //pass.at(Mu1dR) = 1; //value.at(Mu1dR) < .03;
  //pass.at(Mu2dR) = 1; //value.at(Mu2dR) < .03;
  //pass.at(TrkdR) = 1; //value.at(TrkdR) < .03;
  //pass.at(Mu1pt) = value.at(Mu1pt) > 3;
  //pass.at(Mu2pt) = value.at(Mu2pt) > 3;
  //pass.at(Trkpt) = value.at(Trkpt) > 2;
  //pass.at(DsMassCut) = value.at(DsMassCut) < dsMassMax && value.at(DsMassCut) >= dsMassMin;

  double wobs=1;
  double w;  
  double w_peak;     




  if(!Ntp->isData()){w = puWeights->GetBinContent(Ntp->TruthNumberOfInteraction()); w_peak = 1;}//.083;} //.358;}//Ntp->PUReweight(); } //  No weights to data
  else{w=1; w_peak=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){

    int Nvertices(0);
    for(unsigned int l=0; l < Ntp->NSecondaryVertices(); l++){
      TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
      TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(l),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
      if(SVfakePV.DeltaR(SVsignalPV) < 1 && (Ntp->Vertex_Signal_KF_pos(final_idx,true) - Ntp->SecondaryVertexPosition(l)).Mag() > 0.05){
	Nvertices++;
      }
    }
    
    NVtx.at(t).Fill(Ntp->NVtx(),w);

    //Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
    //Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
    //Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    //Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);

    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    //Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    //Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
    //Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    //Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    //Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    //Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    //Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    //Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
 

    Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(0),w);
    Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(1),w);
    Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(2),w);


    //Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
    //Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
    //Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
    //Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);


    //Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
    //Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
    //Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
    //Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);

    //
    //DimuondR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Muon_P4(mu2)),w);
    //Muon1TrkdR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Track_P4(track)),w);
    //Muon2TrkdR.at(t).Fill(Ntp->Muon_P4(mu2).DeltaR(Ntp->Track_P4(track)),w);
    


    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))).M(), w);
    PhiPlusTrackMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))).M(), w);
    
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))).M();
    double dsPt = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))).Pt();
    double dsP = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))).P();
    double dsEta = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))).Eta();
    TLorentzVector DsLV = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(final_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)));
    TLorentzVector DsRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,true);
    TLorentzVector Muon1LV = Ntp->Muon_P4(mu1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(mu2);




    int vertex_idx = final_idx;// + Ntp->NThreeMuons();
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(vertex_idx,true),Ntp->Vertex_MatchedPrimaryVertex(vertex_idx,true));
    TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(mu1),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
    TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(mu2),Ntp->Vertex_MatchedPrimaryVertex(final_idx));



    // Find a tracker muon matched to the track
    // Track/Muon ID
    double dR(999.);
    int MuonClosestToTrack;
    bool MatchedMuonToTrackisFound(false);
    for(int im=0; im<Ntp->NMuons(); im++){
      if(Ntp->Muon_P4(im).DeltaR( Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0))) < dR ) {
	dR=Ntp->Muon_P4(im).DeltaR( Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)));
	MuonClosestToTrack = im;
      }
    }

    if(Ntp->Muon_isTrackerMuon(MuonClosestToTrack) && fabs(Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)).Pt() - Ntp->Muon_P4(MuonClosestToTrack).Pt()) < 0.05){MatchedMuonToTrackisFound=true;
      //        std::cout<<"  "<<fabs(Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)).Pt() - Ntp->Muon_P4(MuonClosestToTrack).Pt()) <<std::endl;
    }

    if(Ntp->Muon_isTrackerMuon(MuonClosestToTrack))MatchedMuonTrackdPT.at(t).Fill(fabs(Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)).Pt() - Ntp->Muon_P4(MuonClosestToTrack).Pt()),1);


    if(MatchedMuonToTrackisFound){
      muonInnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(MuonClosestToTrack);
      muonInnerNValidHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(MuonClosestToTrack);
      muonNLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(MuonClosestToTrack);
      muonValidFraction = Ntp->Muon_innerTrack_validFraction(MuonClosestToTrack);
      muonNLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(MuonClosestToTrack);
      muonNLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(MuonClosestToTrack);
      muonPixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(MuonClosestToTrack);
      muonNMatchedStations = Ntp->Muon_numberOfMatchedStations(MuonClosestToTrack);
      muonPtErrPt = Ntp->Muon_ptErrOverPt(MuonClosestToTrack);
      muonSegComp = Ntp->Muon_segmentCompatibility(MuonClosestToTrack);
      muonCaloComp = Ntp->Muon_caloCompatibility(MuonClosestToTrack);
      muonHad = Ntp->Muon_calEnergy_had(MuonClosestToTrack);
      muonEM = Ntp->Muon_calEnergy_em(MuonClosestToTrack);
      TrackMuBDT_Track = reader_trackerMuonId->EvaluateMVA("BDT");
    }




    for(unsigned int imu=0; imu<2;imu++){
      //      Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu)
      mu_combinedQuality_chi2LocalMomentum=Ntp->Muon_combinedQuality_chi2LocalMomentum(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_combinedQuality_chi2LocalPosition=Ntp->Muon_combinedQuality_chi2LocalPosition(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_combinedQuality_staRelChi2=Ntp->Muon_combinedQuality_staRelChi2(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_combinedQuality_trkRelChi2=Ntp->Muon_combinedQuality_trkRelChi2(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_combinedQuality_globalDeltaEtaPhi=Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_combinedQuality_trkKink=log(Ntp->Muon_combinedQuality_glbKink(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu)));
      mu_combinedQuality_glbKink=log(Ntp->Muon_combinedQuality_trkKink(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu)));
      mu_combinedQuality_glbTrackProbability=Ntp->Muon_combinedQuality_glbTrackProbability(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_Numberofvalidtrackerhits=Ntp->Muon_numberofValidPixelHits(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_Numberofvalidpixelhits=Ntp->Muon_innerTrack_numberOfValidTrackerHits(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_validMuonHitComb=Ntp->Muon_hitPattern_numberOfValidMuonHits(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_numberOfMatchedStations=Ntp->Muon_numberOfMatchedStations(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_segmentCompatibility=Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_timeAtIpInOutErr=Ntp->Muon_timeAtIpInOutErr(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_GLnormChi2=Ntp->Muon_normChi2(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      
      mu_innerTrack_normalizedChi2=Ntp->Muon_innerTrack_normalizedChi2(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_outerTrack_normalizedChi2= Ntp->Muon_outerTrack_normalizedChi2(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
      mu_innerTrack_validFraction=Ntp->Muon_innerTrack_validFraction(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu));
    
      if(fabs(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(imu)).Eta()) < 1.2   )
	{
	  if(imu==0)
	    {
	      Muon1DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }
	  if(imu==1)
	    {
	      Muon2DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }
	}
      else
	{
	  if(imu==0)
	    {
	      Muon1DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	  if(imu==1)
	    {
	      Muon2DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	  
	}
    }
    //    std::cout<<"  "<< Muon1DetID << std::endl;
  











    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    DsMass.at(t).Fill(dsmass,w);

    double SumPT1(0), SumPT8(0);
    int NcloseTracksCount(0);
    int TrackIndex(0);
    double dca_temp(999.);
    for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){
      if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.03){
	  NcloseTracksCount++;
      }

      if( sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2) ) <  dca_temp){
	dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2));
	TrackIndex = i;
      }
      if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.7 && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05 && 
	 sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.05){
        if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(DsRefitLV) < 1.0){
	  SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	}
        if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(DsRefitLV) < 0.8){
	  SumPT8 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	}
      }
    }



    /*    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));

    TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_1),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
    TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_2),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
    TVector3 Mu3ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_3),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
    */

    float var_MaxMuDistanceToSV=std::max({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0,true),2)) ,
	  sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1,true),2) ),
	  sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2,true),2))});
    float var_MinMuDistanceToSV=std::min({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0,true),2)) ,
	  sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1,true),2) ),
	  sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2,true),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2,true),2))});



    float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0,true),
	  Ntp->Vertex_d0sig_reco(final_idx,1,true),
	  Ntp->Vertex_d0sig_reco(final_idx,2,true)});

    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,true),
	  Ntp->Vertex_d0sig_reco(final_idx,1,true),
	  Ntp->Vertex_d0sig_reco(final_idx,2,true)});

    float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
          Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
          Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});
  
    float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
          Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
          Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});    


    float MaxD0BSSignificance = std::max({Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0,true),
	  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1,true),
	  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2,true)});


    if(id==30){

      DsGenMatch.at(t).Fill(Ntp->DsGenMatch(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0) + Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1) + Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0)));

    }
  

  
    double DecayL = SVPV.Mag()*dsmass/dsP;
    if(dsmass > 1.93 && dsmass < 2.01){

      DecayLength_peak.at(t).Fill(DecayL,w);
      Muon1_Pt_peak.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta_peak.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi_peak.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon2_Pt_peak.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta_peak.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi_peak.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Track_Pt_peak.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
      Track_Eta_peak.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
      Track_Phi_peak.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
      
      DsPt_peak.at(t).Fill(dsPt,w);
      DsP_peak.at(t).Fill(dsP,w);
      DsM_peak.at(t).Fill(dsmass,w);
      DsL_peak.at(t).Fill(SVPV.Mag(),w);
      DsEta_peak.at(t).Fill(dsEta,w);

      VertexKFChi2_peak.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,true));
      SVPVDsDirAngle_peak.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
      NtracksClose_peak.at(t).Fill(NcloseTracksCount,1);
      NSV_peak.at(t).Fill(Nvertices,1);
      MinMuon_chi2LocalPosition_peak.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w);

      Mu1Muon_chi2LocalPosition_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),w);
      Mu2Muon_chi2LocalPosition_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),w);



      if(Ntp->NIsolationTrack(final_idx)!=0) MindcaTrackSV_peak.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2)),1);
      MinDca_peak.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1);
      MinD0SigSV_peak.at(t).Fill(MinD0SVSignificance,1);
      MinD0SigPV_peak.at(t).Fill(MinD0Significance,1);

      MaxMuDistanceToSV_peak.at(t).Fill(var_MaxMuDistanceToSV,1);
      MinMuDistanceToSV_peak.at(t).Fill(var_MinMuDistanceToSV,1);


      mu1segmentCompatibility_peak.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0)),1);
      mu2segmentCompatibility_peak.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1)),1);


      MaxVertexPairQuality_peak.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,0,true),Ntp->Vertex_pair_quality(final_idx,1,true),Ntp->Vertex_pair_quality(final_idx,2,true)}),1);
      MaxdeltaMuZ_peak.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
	    fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
	    fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1);
      MaxDca_peak.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1);
      MaxD0SigSV_peak.at(t).Fill(MaxD0SVSignificance,1);
      MaxD0SigPV_peak.at(t).Fill(MaxD0Significance,1);


      MaxD0SigBS_peak.at(t).Fill(MaxD0BSSignificance,1);
      //      MinMuonImpacAngle_peak.at(t).Fill(,1);
      //      MaxMuonImpacAngle_peak.at(t).Fill(,1);

      Muon1ImpactAngle_peak.at(t).Fill(SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),1);
      Muon2ImpactAngle_peak.at(t).Fill(SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),1);


      Mu1trackerLayersWithMeasurement_peak.at(t).Fill(Ntp->Muon_trackerLayersWithMeasurement(mu1),1);

      Iso1_peak.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1);
      Iso1Mu1_peak.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1);
      Iso8Mu1_peak.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1);

      Muon1MVAID_peak.at(t).Fill(Muon1DetID,1);
      Muon2MVAID_peak.at(t).Fill(Muon2DetID,1);

      FLSignificance_peak.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
								      Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_Signal_KF_Covariance(final_idx,true))),w);





      Mu1combinedQuality_chi2LocalMomentum_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu1),1);
      Mu1combinedQuality_chi2LocalPosition_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),1);
      Mu1combinedQuality_staRelChi2_peak.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu1),1);
      Mu1combinedQuality_trkRelChi2_peak.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu1),1);
      Mu1combinedQuality_globalDeltaEtaPhi_peak.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu1),1);
      Mu1mu_combinedQuality_trkKink_peak.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu1)),1);
      Mu1mu_combinedQuality_glbKink_peak.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu1)),1);
      Mu1combinedQuality_glbTrackProbability_peak.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu1),1);
      Mu1Numberofvalidtrackerhits_peak.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu1),1);
      Mu1Numberofvalidpixelhits_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu1),1);
      Mu1validMuonHitComb_peak.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu1),1);
      Mu1numberOfMatchedStations_peak.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu1),1);
      Mu1segmentCompatibility_peak.at(t).Fill(Ntp->Muon_segmentCompatibility(mu1),1);
      Mu1timeAtIpInOutErr_peak.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu1),1);
      Mu1timeAtIpInOut_peak.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu1),1);
      Mu1GLnormChi2_peak.at(t).Fill(Ntp->Muon_normChi2(mu1),1);
      Mu1innerTrack_normalizedChi2_peak.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu1),1);
      Mu1outerTrack_normalizedChi2_peak.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu1),1);
      Mu1innerTrack_validFraction_peak.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu1),1);



      Mu2combinedQuality_chi2LocalMomentum_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu2),1);
      Mu2combinedQuality_chi2LocalPosition_peak.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),1);
      Mu2combinedQuality_staRelChi2_peak.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu2),1);
      Mu2combinedQuality_trkRelChi2_peak.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu2),1);
      Mu2combinedQuality_globalDeltaEtaPhi_peak.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu2),1);
      Mu2mu_combinedQuality_trkKink_peak.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu2)),1);
      Mu2mu_combinedQuality_glbKink_peak.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu2)),1);
      Mu2combinedQuality_glbTrackProbability_peak.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu2),1);
      Mu2Numberofvalidtrackerhits_peak.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu2),1);
      Mu2Numberofvalidpixelhits_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu2),1);
      Mu2validMuonHitComb_peak.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu2),1);
      Mu2numberOfMatchedStations_peak.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu2),1);
      Mu2segmentCompatibility_peak.at(t).Fill(Ntp->Muon_segmentCompatibility(mu2),1);
      Mu2timeAtIpInOutErr_peak.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu2),1);
      Mu2timeAtIpInOut_peak.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu2),1);
      Mu2GLnormChi2_peak.at(t).Fill(Ntp->Muon_normChi2(mu2),1);
      Mu2innerTrack_normalizedChi2_peak.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu2),1);
      Mu2outerTrack_normalizedChi2_peak.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu2),1);
      Mu2innerTrack_validFraction_peak.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu2),1);


      if(MatchedMuonToTrackisFound){
	TrkIDTrackmuonInnerNC2_peak.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(MuonClosestToTrack),1);
	TrkIDTrackmuonInnerNValidHits_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(MuonClosestToTrack),1);
	TrkIDTrackmuonNLostTrackerHits_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerHits(MuonClosestToTrack),1);
	TrkIDTrackmuonValidFraction_peak.at(t).Fill(Ntp->Muon_innerTrack_validFraction(MuonClosestToTrack),1);
	TrkIDTrackmuonNLostTrackerHitsInner_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(MuonClosestToTrack),1);
	TrkIDTrackmuonNLostTrackerHitsOuter_peak.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(MuonClosestToTrack),1);
	TrkIDTrackmuonPixelLayers_peak.at(t).Fill(Ntp->Muon_innerTrack_pixelLayersWithMeasurement(MuonClosestToTrack),1);
	TrkIDTrackmuonNMatchedStations_peak.at(t).Fill(Ntp->Muon_numberOfMatchedStations(MuonClosestToTrack),1);
	TrkIDTrackmuonPtErrPt_peak.at(t).Fill(Ntp->Muon_ptErrOverPt(MuonClosestToTrack),1);
	TrkIDTrackmuonCaloComp_peak.at(t).Fill(Ntp->Muon_caloCompatibility(MuonClosestToTrack),1);
	TrkIDTrackmuonSegComp_peak.at(t).Fill(Ntp->Muon_segmentCompatibility(MuonClosestToTrack),1);
	TrkIDTrackmuonHad_peak.at(t).Fill(Ntp->Muon_calEnergy_had(MuonClosestToTrack),1);
	TrkIDTrackmuonEM_peak.at(t).Fill(Ntp->Muon_calEnergy_em(MuonClosestToTrack),1);
	TrkIDTrackBDT_peak.at(t).Fill(TrackMuBDT_Track,1);
      }



      SVPosX_peak.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).X(),1);
      SVPosY_peak.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Y(),1);
      SVPosZ_peak.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Z(),1);
      PVPosX_peak.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).X(),1);
      PVPosY_peak.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Y(),1);
      PVPosZ_peak.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Z(),1);

      //      SVPosX_peak.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).X(),1);
      //      SVPosY_peak.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Y(),1);

      SVPVEta_peak.at(t).Fill(SVPV.Eta(),1);
      SVPVPhi_peak.at(t).Fill(SVPV.Phi(),1);
      SVPVTransverse_peak.at(t).Fill(SVPV.Perp(),1);


    }

    if(dsmass > 1.70 && dsmass < 1.80){

      DecayLength_sideband.at(t).Fill(DecayL,w);
      Muon1_Pt_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon2_Pt_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Track_Pt_sideband.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
      Track_Eta_sideband.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
      Track_Phi_sideband.at(t).Fill(Ntp->Track_P4(track).Phi(),w);

      DsPt_sideband.at(t).Fill(dsPt,w);
      DsP_sideband.at(t).Fill(dsP,w);
      DsM_sideband.at(t).Fill(dsmass,w);
      DsL_sideband.at(t).Fill(SVPV.Mag(),w);
      DsEta_sideband.at(t).Fill(dsEta,w);

      VertexKFChi2_sideband.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,true));
      SVPVDsDirAngle_sideband.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
      NtracksClose_sideband.at(t).Fill(NcloseTracksCount,1);
      NSV_sideband.at(t).Fill(Nvertices,1);
      MinMuon_chi2LocalPosition_sideband.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w);

      Mu1Muon_chi2LocalPosition_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),w);
      Mu2Muon_chi2LocalPosition_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),w);


      if(Ntp->NIsolationTrack(final_idx)!=0) MindcaTrackSV_sideband.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2)),1);
      MinDca_sideband.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1);
      MinD0SigSV_sideband.at(t).Fill(MinD0SVSignificance,1);
      MinD0SigPV_sideband.at(t).Fill(MinD0Significance,1);
      MaxVertexPairQuality_sideband.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,0,true),Ntp->Vertex_pair_quality(final_idx,1,true),Ntp->Vertex_pair_quality(final_idx,2,true)}),1);
      MaxdeltaMuZ_sideband.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1);
      MaxDca_sideband.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1);
      MaxD0SigSV_sideband.at(t).Fill(MaxD0SVSignificance,1);
      MaxD0SigPV_sideband.at(t).Fill(MaxD0Significance,1);
      MaxD0SigBS_sideband.at(t).Fill(MaxD0BSSignificance,1);
      mu1segmentCompatibility_sideband.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0)),1);
      mu2segmentCompatibility_sideband.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1)),1);

      MaxMuDistanceToSV_sideband.at(t).Fill(var_MaxMuDistanceToSV,1);
      MinMuDistanceToSV_sideband.at(t).Fill(var_MinMuDistanceToSV,1);
      Iso1_sideband.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1);
      Iso1Mu1_sideband.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1);
      Iso8Mu1_sideband.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1);
      FLSignificance_sideband.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
									  Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_Signal_KF_Covariance(final_idx,true))),w);

      Muon1MVAID_sideband.at(t).Fill(Muon1DetID,1);
      Muon2MVAID_sideband.at(t).Fill(Muon2DetID,1);



      Mu1combinedQuality_chi2LocalMomentum_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu1),1);
      Mu1combinedQuality_chi2LocalPosition_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),1);
      Mu1combinedQuality_staRelChi2_sideband.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu1),1);
      Mu1combinedQuality_trkRelChi2_sideband.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu1),1);
      Mu1combinedQuality_globalDeltaEtaPhi_sideband.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu1),1);
      Mu1mu_combinedQuality_trkKink_sideband.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu1)),1);
      Mu1mu_combinedQuality_glbKink_sideband.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu1)),1);
      Mu1combinedQuality_glbTrackProbability_sideband.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu1),1);
      Mu1Numberofvalidtrackerhits_sideband.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu1),1);
      Mu1trackerLayersWithMeasurement_sideband.at(t).Fill(Ntp->Muon_trackerLayersWithMeasurement(mu1),1);
      Mu1Numberofvalidpixelhits_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu1),1);
      Mu1validMuonHitComb_sideband.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu1),1);
      Mu1numberOfMatchedStations_sideband.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu1),1);
      Mu1segmentCompatibility_sideband.at(t).Fill(Ntp->Muon_segmentCompatibility(mu1),1);
      Mu1timeAtIpInOutErr_sideband.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu1),1);
      Mu1timeAtIpInOut_sideband.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu1),1);
      Mu1GLnormChi2_sideband.at(t).Fill(Ntp->Muon_normChi2(mu1),1);
      Mu1innerTrack_normalizedChi2_sideband.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu1),1);
      Mu1outerTrack_normalizedChi2_sideband.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu1),1);
      Mu1innerTrack_validFraction_sideband.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu1),1);





      Mu2combinedQuality_chi2LocalMomentum_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu2),1);
      Mu2combinedQuality_chi2LocalPosition_sideband.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),1);
      Mu2combinedQuality_staRelChi2_sideband.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu2),1);
      Mu2combinedQuality_trkRelChi2_sideband.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu2),1);
      Mu2combinedQuality_globalDeltaEtaPhi_sideband.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu2),1);
      Mu2mu_combinedQuality_trkKink_sideband.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu2)),1);
      Mu2mu_combinedQuality_glbKink_sideband.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu2)),1);
      Mu2combinedQuality_glbTrackProbability_sideband.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu2),1);
      Mu2Numberofvalidtrackerhits_sideband.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu2),1);
      Mu2Numberofvalidpixelhits_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu2),1);
      Mu2validMuonHitComb_sideband.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu2),1);
      Mu2numberOfMatchedStations_sideband.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu2),1);
      Mu2segmentCompatibility_sideband.at(t).Fill(Ntp->Muon_segmentCompatibility(mu2),1);
      Mu2timeAtIpInOutErr_sideband.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu2),1);
      Mu2timeAtIpInOut_sideband.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu2),1);
      Mu2GLnormChi2_sideband.at(t).Fill(Ntp->Muon_normChi2(mu2),1);
      Mu2innerTrack_normalizedChi2_sideband.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu2),1);
      Mu2outerTrack_normalizedChi2_sideband.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu2),1);
      Mu2innerTrack_validFraction_sideband.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu2),1);


      Muon1ImpactAngle_sideband.at(t).Fill(SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),1);
      Muon2ImpactAngle_sideband.at(t).Fill(SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),1);


      if(MatchedMuonToTrackisFound){
        TrkIDTrackmuonInnerNC2_sideband.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(MuonClosestToTrack),1);
        TrkIDTrackmuonInnerNValidHits_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHits_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonValidFraction_sideband.at(t).Fill(Ntp->Muon_innerTrack_validFraction(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHitsInner_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHitsOuter_sideband.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(MuonClosestToTrack),1);
        TrkIDTrackmuonPixelLayers_sideband.at(t).Fill(Ntp->Muon_innerTrack_pixelLayersWithMeasurement(MuonClosestToTrack),1);
        TrkIDTrackmuonNMatchedStations_sideband.at(t).Fill(Ntp->Muon_numberOfMatchedStations(MuonClosestToTrack),1);
        TrkIDTrackmuonPtErrPt_sideband.at(t).Fill(Ntp->Muon_ptErrOverPt(MuonClosestToTrack),1);
        TrkIDTrackmuonCaloComp_sideband.at(t).Fill(Ntp->Muon_caloCompatibility(MuonClosestToTrack),1);
	TrkIDTrackmuonSegComp_sideband.at(t).Fill(Ntp->Muon_segmentCompatibility(MuonClosestToTrack),1);
        TrkIDTrackmuonHad_sideband.at(t).Fill(Ntp->Muon_calEnergy_had(MuonClosestToTrack),1);
        TrkIDTrackmuonEM_sideband.at(t).Fill(Ntp->Muon_calEnergy_em(MuonClosestToTrack),1);
        TrkIDTrackBDT_sideband.at(t).Fill(TrackMuBDT_Track,1);
      }


      SVPosX_sideband.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).X(),1);
      SVPosY_sideband.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Y(),1);
      SVPosZ_sideband.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Z(),1);
      PVPosX_sideband.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).X(),1);
      PVPosY_sideband.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Y(),1);
      PVPosZ_sideband.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Z(),1);
      SVPVEta_sideband.at(t).Fill(SVPV.Eta(),1);
      SVPVPhi_sideband.at(t).Fill(SVPV.Phi(),1);
      SVPVTransverse_sideband.at(t).Fill(SVPV.Perp(),1);


    }
  
    bool isPrompt(true);
    if(id!=1){

      control_Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*w_peak);
      control_Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*w_peak);
      control_Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w*w_peak);
      control_Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*w_peak);
      control_Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*w_peak);
      control_Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w*w_peak);
      control_Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w*w_peak);
      control_Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w*w_peak);
      control_Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w*w_peak);

      Ds_Pt.at(t).Fill(dsPt,w*w_peak);
      Ds_P.at(t).Fill(dsP,w*w_peak);
      Ds_M.at(t).Fill(dsmass,w*w_peak);
      Ds_L.at(t).Fill(SVPV.Mag(),w*w_peak);
      Ds_eta.at(t).Fill(dsEta,w*w_peak);

      VertexKFChi2.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,true),w*w_peak);
      SVPVDsDirAngle.at(t).Fill(SVPV.Angle(DsLV.Vect()),w*w_peak);
      NtracksClose.at(t).Fill(NcloseTracksCount,1*w_peak);
      NSV.at(t).Fill(Nvertices,1*w_peak);
      MinMuon_chi2LocalPosition.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w*w_peak);
      Mu1Muon_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),1*w_peak);
      Mu2Muon_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),1*w_peak);




      if(Ntp->NIsolationTrack(final_idx)!=0) MindcaTrackSV.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2)),1*w_peak);
      MinDca.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1*w_peak);
      MinD0SigSV.at(t).Fill(MinD0SVSignificance,1*w_peak);
      MinD0SigPV.at(t).Fill(MinD0Significance,1*w_peak);
      MaxVertexPairQuality.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,0,true),Ntp->Vertex_pair_quality(final_idx,1,true),Ntp->Vertex_pair_quality(final_idx,2,true)}),1*w_peak);
      MaxdeltaMuZ.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1*w_peak);
      MaxDca.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),1*w_peak);
      MaxD0SigSV.at(t).Fill(MaxD0SVSignificance,1*w_peak);
      MaxD0SigPV.at(t).Fill(MaxD0Significance,1*w_peak);
      MaxD0SigBS.at(t).Fill(MaxD0BSSignificance,1*w_peak);

      mu1segmentCompatibility.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0)),1*w_peak);
      mu2segmentCompatibility.at(t).Fill(Ntp->Muon_segmentCompatibility(Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1)),1*w_peak);

      MaxMuDistanceToSV.at(t).Fill(var_MaxMuDistanceToSV,1*w_peak);
      MinMuDistanceToSV.at(t).Fill(var_MinMuDistanceToSV,1*w_peak);

      Iso1.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1*w_peak);
      Iso1Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1*w_peak);
      Iso8Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1*w_peak);
      FLSignificance.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
                                                                   Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_Signal_KF_Covariance(final_idx,true))),w*w_peak);



      Muon1MVAID.at(t).Fill(Muon1DetID,1*w_peak);
      Muon2MVAID.at(t).Fill(Muon2DetID,1*w_peak);


      Mu1combinedQuality_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu1),1*w_peak);
      Mu1combinedQuality_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),1*w_peak);
      Mu1combinedQuality_staRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu1),1*w_peak);
      Mu1combinedQuality_trkRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu1),1*w_peak);
      Mu1combinedQuality_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu1),1*w_peak);
      Mu1mu_combinedQuality_trkKink.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu1)),1*w_peak);
      Mu1mu_combinedQuality_glbKink.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu1)),1*w_peak);
      Mu1combinedQuality_glbTrackProbability.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu1),1*w_peak);
      Mu1Numberofvalidtrackerhits.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu1),1*w_peak);
      Mu1Numberofvalidpixelhits.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu1),1*w_peak);
      Mu1validMuonHitComb.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu1),1*w_peak);
      Mu1trackerLayersWithMeasurement.at(t).Fill(Ntp->Muon_trackerLayersWithMeasurement(mu1),1*w_peak);
      Mu1numberOfMatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu1),1*w_peak);
      Mu1segmentCompatibility.at(t).Fill(Ntp->Muon_segmentCompatibility(mu1),1*w_peak);
      Mu1timeAtIpInOutErr.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu1),1*w_peak);
      Mu1timeAtIpInOut.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu1),1*w_peak);
      Mu1GLnormChi2.at(t).Fill(Ntp->Muon_normChi2(mu1),1*w_peak);
      Mu1innerTrack_normalizedChi2.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu1),1*w_peak);
      Mu1outerTrack_normalizedChi2.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu1),1*w_peak);
      Mu1innerTrack_validFraction.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu1),1*w_peak);




      Mu2combinedQuality_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(mu2),1*w_peak);
      Mu2combinedQuality_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(mu2),1*w_peak);
      Mu2combinedQuality_staRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_staRelChi2(mu2),1*w_peak);
      Mu2combinedQuality_trkRelChi2.at(t).Fill(Ntp->Muon_combinedQuality_trkRelChi2(mu2),1*w_peak);
      Mu2combinedQuality_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(mu2),1*w_peak);
      Mu2mu_combinedQuality_trkKink.at(t).Fill(log(Ntp->Muon_combinedQuality_glbKink(mu2)),1*w_peak);
      Mu2mu_combinedQuality_glbKink.at(t).Fill(log(Ntp->Muon_combinedQuality_trkKink(mu2)),1*w_peak);
      Mu2combinedQuality_glbTrackProbability.at(t).Fill(Ntp->Muon_combinedQuality_glbTrackProbability(mu2),1*w_peak);
      Mu2Numberofvalidtrackerhits.at(t).Fill(Ntp->Muon_numberofValidPixelHits(mu2),1*w_peak);
      Mu2Numberofvalidpixelhits.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu2),1*w_peak);
      Mu2validMuonHitComb.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(mu2),1*w_peak);
      Mu2numberOfMatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(mu2),1*w_peak);
      Mu2segmentCompatibility.at(t).Fill(Ntp->Muon_segmentCompatibility(mu2),1*w_peak);
      Mu2timeAtIpInOutErr.at(t).Fill(Ntp->Muon_timeAtIpInOutErr(mu2),1*w_peak);
      Mu2timeAtIpInOut.at(t).Fill(Ntp->Muon_timeAtIpInOut(mu2),1*w_peak);
      Mu2GLnormChi2.at(t).Fill(Ntp->Muon_normChi2(mu2),1*w_peak);
      Mu2innerTrack_normalizedChi2.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(mu2),1*w_peak);
      Mu2outerTrack_normalizedChi2.at(t).Fill( Ntp->Muon_outerTrack_normalizedChi2(mu2),1*w_peak);
      Mu2innerTrack_validFraction.at(t).Fill(Ntp->Muon_innerTrack_validFraction(mu2),1*w_peak);



      Muon1ImpactAngle.at(t).Fill(SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),1);
      Muon2ImpactAngle.at(t).Fill(SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),1);
      
      if(MatchedMuonToTrackisFound){
        TrkIDTrackmuonInnerNC2.at(t).Fill(Ntp->Muon_innerTrack_normalizedChi2(MuonClosestToTrack),1);
        TrkIDTrackmuonInnerNValidHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfValidTrackerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHits.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonValidFraction.at(t).Fill(Ntp->Muon_innerTrack_validFraction(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHitsInner.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(MuonClosestToTrack),1);
        TrkIDTrackmuonNLostTrackerHitsOuter.at(t).Fill(Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(MuonClosestToTrack),1);
        TrkIDTrackmuonPixelLayers.at(t).Fill(Ntp->Muon_innerTrack_pixelLayersWithMeasurement(MuonClosestToTrack),1);
        TrkIDTrackmuonNMatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(MuonClosestToTrack),1);
        TrkIDTrackmuonPtErrPt.at(t).Fill(Ntp->Muon_ptErrOverPt(MuonClosestToTrack),1);
        TrkIDTrackmuonCaloComp.at(t).Fill(Ntp->Muon_caloCompatibility(MuonClosestToTrack),1);
	TrkIDTrackmuonSegComp.at(t).Fill(Ntp->Muon_segmentCompatibility(MuonClosestToTrack),1);
        TrkIDTrackmuonHad.at(t).Fill(Ntp->Muon_calEnergy_had(MuonClosestToTrack),1);
        TrkIDTrackmuonEM.at(t).Fill(Ntp->Muon_calEnergy_em(MuonClosestToTrack),1);
        TrkIDTrackBDT.at(t).Fill(TrackMuBDT_Track,1);
      }


     SVPosX.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).X(),1);
     SVPosY.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Y(),1);
     SVPosZ.at(t).Fill(Ntp->Vertex_Signal_KF_pos(final_idx,true).Z(),1);
     PVPosX.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).X(),1);
     PVPosY.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Y(),1);
     PVPosZ.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true).Z(),1);
     SVPVEta.at(t).Fill(SVPV.Eta(),1);
     SVPVPhi.at(t).Fill(SVPV.Phi(),1);
     SVPVTransverse.at(t).Fill(SVPV.Perp(),1);


    

      for (unsigned int isigp=0; isigp<Ntp->NMCSignalParticles(); isigp++){
        for (int is=0; is<Ntp->NMCSignalParticleSources(isigp); is++){
          if (abs(Ntp->MCSignalParticle_Sourcepdgid(isigp,is))>500){
            isPrompt=false;
          }
        }
      }

      if (isPrompt) {DecayLength_prompt.at(t).Fill(DecayL,w*w_peak);}
      else if (!isPrompt) {DecayLength_non_prompt.at(t).Fill(DecayL,w*w_peak);}
      DecayLength.at(t).Fill(DecayL,w*w_peak);
    }


  }

}

void  DsPhiPeak_2018::Finish(){

  int id(Ntp->GetMCID());
  if (id==1) {

    std::vector<double> scaleRun;

    //if(RunB){scaleRun.push_back(2265.07);scaleRun.push_back(4586.64);} 
    //if(RunC){scaleRun.push_back(15965.4);scaleRun.push_back(24220.4);}
    //if(RunD){scaleRun.push_back(5952.73);scaleRun.push_back(11303.6);}
    //if(RunE){scaleRun.push_back(10661.2);scaleRun.push_back(19461.1);}
    //if(RunF){scaleRun.push_back(10093.0);scaleRun.push_back(19046.7);}
    //scaleRun.push_back(1122.09);scaleRun.push_back(2246.3); //For 2017
    scaleRun.push_back(543.064);scaleRun.push_back(4377.86); //For 2018

    for ( unsigned int j=0; j<validationCollection.size(); ++j){
      (validationCollection.at(j)->at(0)).Add(&(peakCollection.at(j)->at(0)),1.0);
      (validationCollection.at(j)->at(0)).Add(&(sidebandCollection.at(j)->at(0)),-(scaleRun[0]/scaleRun[1]));
    }

    DecayLength_prompt.at(0).Add(&DecayLength_peak.at(0));
    DecayLength_prompt.at(0).Add(&DecayLength_sideband.at(0),-(scaleRun[0]/scaleRun[1]));
    DecayLength_non_prompt.at(0).Add(&DecayLength_peak.at(0));
    DecayLength_non_prompt.at(0).Add(&DecayLength_sideband.at(0),-(scaleRun[0]/scaleRun[1]));
 
      }

  if(mode == RECONSTRUCT){
  //  for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
  //    double scale(1.);
  //    if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/1;
  //    ScaleAllHistOfType(i,scale);
  //  }

    for ( unsigned int j=0; j<validationCollection.size(); ++j){
      std::cout << "Plot " << j << " has data_Int = " << validationCollection.at(j)->at(0).Integral() << ", MC_Int = " << validationCollection.at(j)->at(1).Integral() << ", and ";
      std::cout << "Scale = " << validationCollection.at(j)->at(0).Integral()/validationCollection.at(j)->at(1).Integral() << std::endl;
      validationCollection.at(j)->at(1).Scale(validationCollection.at(j)->at(0).Integral()/validationCollection.at(j)->at(1).Integral());
    }

  }
  Selection::Finish();

}



