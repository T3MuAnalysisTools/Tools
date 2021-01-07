#include "CommonSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace std;

CommonSelector::CommonSelector(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0),
  phiVetoCut1(0.96),
  phiVetoCut2(1.07),
  rmgCutVeto1(0.77),  // rmg = rho&omega
  rmgCutVeto2(0.812), // rmg = rho&omega
  PEMassResolutionCut1_(0.007),
  PEMassResolutionCut2_(0.01),
  mvaA1_(0.0927636), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2_(0.160007),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1_(0.113621),
  mvaB2_(0.206321),
  mvaC1_(0.139706),
  mvaC2_(0.219904)
{
  // This is a class constructor;
}



CommonSelector::~CommonSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  CommonSelector::Configure(){


  //  This mini tree is for limit extraction

  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("bdt",&bdt);
  T3MMiniTree->Branch("category",&category);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  T3MMiniTree->Branch("LumiScale",&LumiScale);
  T3MMiniTree->Branch("A1",&mvaA1);
  T3MMiniTree->Branch("A2",&mvaA2);
  T3MMiniTree->Branch("B1",&mvaB1);
  T3MMiniTree->Branch("B2",&mvaB2);
  T3MMiniTree->Branch("C1",&mvaC1);
  T3MMiniTree->Branch("C2",&mvaC2);




  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  //*** defined the bdt reader for event selection; readerA- category A, readerB - category B ...
  readerA = new TMVA::Reader( "!Color:!Silent" );
  readerA->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2);
  readerA->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerA->AddVariable( "var_flightLenSig", &var_flightLenSig);
  readerA->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerA->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerA->AddVariable( "var_maxMuonsDca", &var_maxMuonsDca);
  readerA->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerA->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerA->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerA->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);

  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->BookMVA( "BDT", basedir+"weights/August_A_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerB = new TMVA::Reader( "!Color:!Silent" );
  readerB->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2);
  readerB->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerB->AddVariable( "var_flightLenSig", &var_flightLenSig);
  readerB->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerB->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerB->AddVariable( "var_maxMuonsDca", &var_maxMuonsDca);
  readerB->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerB->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerB->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerB->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerB->AddSpectator("var_tauMass",&var_tauMass);
  readerB->BookMVA( "BDT", basedir+"weights/August_B_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerC = new TMVA::Reader( "!Color:!Silent" );
  readerC->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2);
  readerC->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerC->AddVariable( "var_flightLenSig", &var_flightLenSig);
  readerC->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerC->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerC->AddVariable( "var_maxMuonsDca", &var_maxMuonsDca);
  readerC->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerC->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerC->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerC->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerC->AddSpectator("var_tauMass",&var_tauMass);
  readerC->BookMVA( "BDT", basedir+"weights/August_C_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles




  //*** Muon MVA ID  readers

  readerMuIDBarrel= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  readerMuIDBarrel->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDBarrel->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
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
  readerMuIDBarrel->BookMVA( "BDT", basedir+"MuonMVA_02may_barrel/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles


  readerMuIDEndcap= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
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
  readerMuIDEndcap->BookMVA( "BDT", basedir+"MuonMVA_02may_endcap/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles




  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1T)                cut.at(L1T)=1;
    if(i==HLT)                cut.at(HLT)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto1)           cut.at(PhiVeto1)=0; // defined below
    if(i==OmegaVeto1)         cut.at(OmegaVeto1)=0; // defined below
    if(i==PhiVeto2)           cut.at(PhiVeto2)=0; // defined below
    if(i==OmegaVeto2)         cut.at(OmegaVeto2)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;

    if(i==L1T){
      title.at(i)="L1T trigger ";
      hlabel="Level 1 Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1T_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1T_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

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
      title.at(i)="Muons GL and PF";
      hlabel="gl,gl,gl";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==PhiVeto1){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
    }

    else if(i==OmegaVeto1){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
    }

    else if(i==PhiVeto2){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
    }

    else if(i==OmegaVeto2){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
    }

    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
    }

    else if(i==TauMassCut){
      title.at(i)="$\\tau$ mass 1.6 - 2 GeV ";
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
  TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVASV =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVASV","Mu1TrackInvariantMassBeforeMVASV",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVASV =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVASV","Mu2TrackInvariantMassBeforeMVASV",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVASV =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVASV","Mu3TrackInvariantMassBeforeMVASV",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVASVAngle0 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVASVAngle0","Mu1TrackInvariantMassBeforeMVASVAngle0",100,0,2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVASVAngle0 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVASVAngle0","Mu2TrackInvariantMassBeforeMVASVAngle0",100,0,2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVASVAngle0 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVASVAngle0","Mu3TrackInvariantMassBeforeMVASVAngle0",100,0,2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVASVAngle1 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVASVAngle1","Mu1TrackInvariantMassBeforeMVASVAngle1",50,0,0.1,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVASVAngle1 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVASVAngle1","Mu2TrackInvariantMassBeforeMVASVAngle1",50,0,0.1,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVASVAngle1 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVASVAngle1","Mu3TrackInvariantMassBeforeMVASVAngle1",50,0,0.1,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVASVAngle2 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVASVAngle2","Mu1TrackInvariantMassBeforeMVASVAngle2",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVASVAngle2 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVASVAngle2","Mu2TrackInvariantMassBeforeMVASVAngle2",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVASVAngle2 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVASVAngle2","Mu3TrackInvariantMassBeforeMVASVAngle2",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVASVAngle3 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVASVAngle3","Mu1TrackInvariantMassBeforeMVASVAngle3",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVASVAngle3 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVASVAngle3","Mu2TrackInvariantMassBeforeMVASVAngle3",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVASVAngle3 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVASVAngle3","Mu3TrackInvariantMassBeforeMVASVAngle3",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarSV =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarSV","Mu1TrackInvariantMassBeforeMVAKStarSV",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarSV =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarSV","Mu2TrackInvariantMassBeforeMVAKStarSV",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarSV =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarSV","Mu3TrackInvariantMassBeforeMVAKStarSV",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarSVAngle0 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarSVAngle0","Mu1TrackInvariantMassBeforeMVAKStarSVAngle0",100,0,2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarSVAngle0 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarSVAngle0","Mu2TrackInvariantMassBeforeMVAKStarSVAngle0",100,0,2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarSVAngle0 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarSVAngle0","Mu3TrackInvariantMassBeforeMVAKStarSVAngle0",100,0,2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarSVAngle1 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarSVAngle1","Mu1TrackInvariantMassBeforeMVAKStarSVAngle1",50,0,0.1,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarSVAngle1 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarSVAngle1","Mu2TrackInvariantMassBeforeMVAKStarSVAngle1",50,0,0.1,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarSVAngle1 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarSVAngle1","Mu3TrackInvariantMassBeforeMVAKStarSVAngle1",50,0,0.1,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarSVAngle2 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarSVAngle2","Mu1TrackInvariantMassBeforeMVAKStarSVAngle2",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarSVAngle2 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarSVAngle2","Mu2TrackInvariantMassBeforeMVAKStarSVAngle2",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarSVAngle2 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarSVAngle2","Mu3TrackInvariantMassBeforeMVAKStarSVAngle2",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarSVAngle3 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarSVAngle3","Mu1TrackInvariantMassBeforeMVAKStarSVAngle3",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarSVAngle3 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarSVAngle3","Mu2TrackInvariantMassBeforeMVAKStarSVAngle3",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarSVAngle3 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarSVAngle3","Mu3TrackInvariantMassBeforeMVAKStarSVAngle3",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVA","Mu1TrackInvariantMassBeforeMVA",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVA","Mu2TrackInvariantMassBeforeMVA",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVA","Mu3TrackInvariantMassBeforeMVA",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStar =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStar","Mu1TrackInvariantMassBeforeMVAKStar",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStar =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStar","Mu2TrackInvariantMassBeforeMVAKStar",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStar =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStar","Mu3TrackInvariantMassBeforeMVAKStar",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAFiner =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAFiner","Mu1TrackInvariantMassBeforeMVAFiner",80,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAFiner =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAFiner","Mu2TrackInvariantMassBeforeMVAFiner",80,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAFiner =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAFiner","Mu3TrackInvariantMassBeforeMVAFiner",80,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVACoarser =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVACoarser","Mu1TrackInvariantMassBeforeMVACoarser",35,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVACoarser =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVACoarser","Mu2TrackInvariantMassBeforeMVACoarser",35,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVACoarser =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVACoarser","Mu3TrackInvariantMassBeforeMVACoarser",35,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  MuOSSS1InvariantMassBeforeMVA =HConfig.GetTH1D(Name+"_MuOSSS1InvariantMassBeforeMVA","MuOSSS1InvariantMassBeforeMVA",200,0.95,1.2,"OS - SS1 #mu Invariant Mass, GeV","Events");
  MuOSSS2InvariantMassBeforeMVA =HConfig.GetTH1D(Name+"_MuOSSS2InvariantMassBeforeMVA","MuOSSS2InvariantMassBeforeMVA",200,0.95,1.2,"OS - SS2 #mu Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdR =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdR","Mu1TrackInvariantMassBeforeMVABestdR",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdR =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdR","Mu2TrackInvariantMassBeforeMVABestdR",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdR =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdR","Mu3TrackInvariantMassBeforeMVABestdR",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdR =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdR","Mu1TrackInvariantMassBeforeMVAKStarBestdR",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdR =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdR","Mu2TrackInvariantMassBeforeMVAKStarBestdR",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdR =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdR","Mu3TrackInvariantMassBeforeMVAKStarBestdR",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRIncrease =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRIncrease","Mu1TrackInvariantMassBeforeMVABestdRIncrease",100,0,0.1,"SV dR SS2","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRDecrease =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRDecrease","Mu1TrackInvariantMassBeforeMVABestdRDecrease",100,0,0.2,"#Delta #theta - #pi SS2","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd","Mu1TrackInvariantMassBeforeMVABestdRLtd",100,0,30,"OS #mu - Isolation Track Chi2","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd","Mu2TrackInvariantMassBeforeMVABestdRLtd",100,0,30,"SS1 #mu - Isolation Track Chi2","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd","Mu3TrackInvariantMassBeforeMVABestdRLtd",100,0,30,"SS2 #mu - Isolation Track Chi2","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd1 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd1","Mu1TrackInvariantMassBeforeMVABestdRLtd1",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd1 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd1","Mu2TrackInvariantMassBeforeMVABestdRLtd1",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd1 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd1","Mu3TrackInvariantMassBeforeMVABestdRLtd1",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd2 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd2","Mu1TrackInvariantMassBeforeMVABestdRLtd2",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd2 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd2","Mu2TrackInvariantMassBeforeMVABestdRLtd2",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd2 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd2","Mu3TrackInvariantMassBeforeMVABestdRLtd2",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd3 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd3","Mu1TrackInvariantMassBeforeMVABestdRLtd3",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd3 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd3","Mu2TrackInvariantMassBeforeMVABestdRLtd3",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd3 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd3","Mu3TrackInvariantMassBeforeMVABestdRLtd3",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd4 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd4","Mu1TrackInvariantMassBeforeMVABestdRLtd4",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd4 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd4","Mu2TrackInvariantMassBeforeMVABestdRLtd4",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd4 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd4","Mu3TrackInvariantMassBeforeMVABestdRLtd4",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd5 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd5","Mu1TrackInvariantMassBeforeMVABestdRLtd5",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd5 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd5","Mu2TrackInvariantMassBeforeMVABestdRLtd5",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd5 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd5","Mu3TrackInvariantMassBeforeMVABestdRLtd5",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd6 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd6","Mu1TrackInvariantMassBeforeMVABestdRLtd6",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd6 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd6","Mu2TrackInvariantMassBeforeMVABestdRLtd6",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd6 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd6","Mu3TrackInvariantMassBeforeMVABestdRLtd6",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestdRLtd7 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRLtd7","Mu1TrackInvariantMassBeforeMVABestdRLtd7",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestdRLtd7 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestdRLtd7","Mu2TrackInvariantMassBeforeMVABestdRLtd7",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestdRLtd7 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestdRLtd7","Mu3TrackInvariantMassBeforeMVABestdRLtd7",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7 =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7","Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7 =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7","Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7 =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7","Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVABestMass =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestMass","Mu1TrackInvariantMassBeforeMVABestMass",40,0.95,1.07,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVABestMass =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVABestMass","Mu2TrackInvariantMassBeforeMVABestMass",40,0.95,1.07,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVABestMass =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVABestMass","Mu3TrackInvariantMassBeforeMVABestMass",40,0.95,1.07,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassBeforeMVAKStarBestMass =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVAKStarBestMass","Mu1TrackInvariantMassBeforeMVAKStarBestMass",100,0.6,1.6,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassBeforeMVAKStarBestMass =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassBeforeMVAKStarBestMass","Mu2TrackInvariantMassBeforeMVAKStarBestMass",100,0.6,1.6,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassBeforeMVAKStarBestMass =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassBeforeMVAKStarBestMass","Mu3TrackInvariantMassBeforeMVAKStarBestMass",100,0.6,1.6,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  IsolationTrackCount =HConfig.GetTH1D(Name+"_IsolationTrackCount","IsolationTrackCount",21,-0.5,20.5,"No of tracks","Events");
  
  Mu1TrackInvariantMassAfterA1MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterA1MVA","Mu1TrackInvariantMassAfterA1MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterA1MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterA1MVA","Mu2TrackInvariantMassAfterA1MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterA1MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterA1MVA","Mu3TrackInvariantMassAfterA1MVA",150,0.9,1.2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassAfterB1MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterB1MVA","Mu1TrackInvariantMassAfterB1MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterB1MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterB1MVA","Mu2TrackInvariantMassAfterB1MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterB1MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterB1MVA","Mu3TrackInvariantMassAfterB1MVA",150,0.9,1.2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassAfterC1MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterC1MVA","Mu1TrackInvariantMassAfterC1MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterC1MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterC1MVA","Mu2TrackInvariantMassAfterC1MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterC1MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterC1MVA","Mu3TrackInvariantMassAfterC1MVA",150,0.9,1.2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassAfterA2MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterA2MVA","Mu1TrackInvariantMassAfterA2MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterA2MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterA2MVA","Mu2TrackInvariantMassAfterA2MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterA2MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterA2MVA","Mu3TrackInvariantMassAfterA2MVA",150,0.9,1.2,"SS2 #mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassAfterB2MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterB2MVA","Mu1TrackInvariantMassAfterB2MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterB2MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterB2MVA","Mu2TrackInvariantMassAfterB2MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterB2MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterB2MVA","Mu3TrackInvariantMassAfterB2MVA",150,0.9,1.2,"SS2#mu - Isolation Track Invariant Mass, GeV","Events");
  
  Mu1TrackInvariantMassAfterC2MVA =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassAfterC2MVA","Mu1TrackInvariantMassAfterC2MVA",150,0.9,1.2,"OS #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu2TrackInvariantMassAfterC2MVA =HConfig.GetTH1D(Name+"_Mu2TrackInvariantMassAfterC2MVA","Mu2TrackInvariantMassAfterC2MVA",150,0.9,1.2,"SS1 #mu - Isolation Track Invariant Mass, GeV","Events");
  Mu3TrackInvariantMassAfterC2MVA =HConfig.GetTH1D(Name+"_Mu3TrackInvariantMassAfterC2MVA","Mu3TrackInvariantMassAfterC2MVA",150,0.9,1.2,"SS2#mu - Isolation Track Invariant Mass, GeV","Events");
  
  TauAngleTest  =HConfig.GetTH1D(Name+"_TauAngleTest","3#mu  mass",60,1.65,1.95,"  M_{#tau} , GeV","Events");
  
  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");
  TauMass_all =HConfig.GetTH1D(Name+"_TauMass_all","3#mu  mass",60,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMass_allVsBDTA=HConfig.GetTH2D(Name+"_TauMass_allVsBDTA","3#mu mass vs BDTa",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTB=HConfig.GetTH2D(Name+"_TauMass_allVsBDTB","3#mu mass vs BDTb",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTC=HConfig.GetTH2D(Name+"_TauMass_allVsBDTC","3#mu mass vs BDTc",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");

  EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

  TauMassA1 =HConfig.GetTH1D(Name+"_TauMassA1","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitA1 =HConfig.GetTH1D(Name+"_TauMassRefitA1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (A1)","Events");
  TauMassRefitA1MassCut =HConfig.GetTH1D(Name+"_TauMassRefitA1MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (A1)","Events");
  TauMassRefitA2MassCut =HConfig.GetTH1D(Name+"_TauMassRefitA2MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (A2)","Events");

  TauMassRefitA1HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitA1HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (A1)","Events");
  TauMassRefitA2HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitA2HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (A2)","Events");

  TauMassRefitA1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitA1FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (A1)","Events");
  TauMassRefitA2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitA2FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (A2)","Events");



  TauMassRefitABC1 =HConfig.GetTH1D(Name+"_TauMassRefitABC1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2 =HConfig.GetTH1D(Name+"_TauMassRefitABC2","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (inclusive ABC2)","Events");


  TauMassB1 =HConfig.GetTH1D(Name+"_TauMassB1","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitB1 =HConfig.GetTH1D(Name+"_TauMassRefitB1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (B1)","Events");
  TauMassRefitB1MassCut =HConfig.GetTH1D(Name+"_TauMassRefitB1MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B1)","Events");
  TauMassRefitB2MassCut =HConfig.GetTH1D(Name+"_TauMassRefitB2MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B2)","Events");


  TauMassRefitB1HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitB1HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (B1)","Events");
  TauMassRefitB2HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitB2HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (B2)","Events");


  TauMassRefitB1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitB1FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B1)","Events");
  TauMassRefitB2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitB2FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B2)","Events");

  TauMassC1 =HConfig.GetTH1D(Name+"_TauMassC1","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitC1 =HConfig.GetTH1D(Name+"_TauMassRefitC1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (C1)","Events");
  TauMassRefitC1MassCut =HConfig.GetTH1D(Name+"_TauMassRefitC1MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (C1)","Events");
  TauMassRefitC2MassCut =HConfig.GetTH1D(Name+"_TauMassRefitC2MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (C2)","Events");

  TauMassRefitC1HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitC1HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (C1)","Events");
  TauMassRefitC2HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitC2HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (C2)","Events");


  TauMassRefitC1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitC1FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (C1)","Events");
  TauMassRefitC2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitC2FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (C2)","Events");

  TauMassA2 =HConfig.GetTH1D(Name+"_TauMassA2","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV (A2)","Events");
  TauMassRefitA2 =HConfig.GetTH1D(Name+"_TauMassRefitA2","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (A2)","Events");


  TauMassB2 =HConfig.GetTH1D(Name+"_TauMassB2","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV (B2)","Events");
  TauMassRefitB2 =HConfig.GetTH1D(Name+"_TauMassRefitB2","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (B2)","Events");


  TauMassC2 =HConfig.GetTH1D(Name+"_TauMassC2","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV (C2)","Events");
  TauMassRefitC2 =HConfig.GetTH1D(Name+"_TauMassRefitC2","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (C2)","Events");
  //TauMassRefitC2FullEtaVetoCut

  EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");

  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");

  BDTOutputA = HConfig.GetTH1D(Name+"_BDTOutputA","BDTOutputA",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputB = HConfig.GetTH1D(Name+"_BDTOutputB","BDTOutputB",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputC = HConfig.GetTH1D(Name+"_BDTOutputC","BDTOutputC",50,-0.4,0.4,"BDT Output","Events");



  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");
  PairMass=HConfig.GetTH2D(Name+"_PairMass","PairMass",100,0.2,1.8,100,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassFinalSel=HConfig.GetTH2D(Name+"_PairMassFinalSel","PairMassFinalSel",60,0.2,1.8,60,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1=HConfig.GetTH1D(Name+"_PairMass1","PairMass1",80,0.2,1.777,"M_{1}, GeV","");
  PairMass2=HConfig.GetTH1D(Name+"_PairMass2","PairMass2",80,0.2,1.777,"M_{2}, GeV","");





  PairMassWithCut=HConfig.GetTH2D(Name+"_PairMassWithCut","PairMassWithCut",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEta=HConfig.GetTH2D(Name+"_PairMassEta","PairMassEta",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEtaPrime=HConfig.GetTH2D(Name+"_PairMassEtaPrime","PairMassEtaPrime",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");

  IDOriginOfOSMuon =HConfig.GetTH1D(Name+"_IDOriginOfOSMuon","IDOriginOfOSMuon",400,200,600,"PDGID of OS muon origin","Events");

  Muon1MVAID=HConfig.GetTH1D(Name+"_Muon1MVAID","Muon1MVAID",50,0.0,1.0,"#mu_{1} MVA","Events");
  Muon2MVAID=HConfig.GetTH1D(Name+"_Muon2MVAID","Muon2MVAID",50,0.0,1.0,"#mu_{2} MVA","Events");
  Muon3MVAID=HConfig.GetTH1D(Name+"_Muon3MVAID","Muon3MVAID",50,0.0,1.0,"#mu_{3} MVA","Events");







  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  CommonSelector::Store_ExtraDist(){ 


  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  
  Extradist1d.push_back(&TauAngleTest);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVASV);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVASV);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVASV);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVASVAngle0);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVASVAngle0);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVASVAngle0);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVASVAngle1);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVASVAngle1);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVASVAngle1);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVASVAngle2);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVASVAngle2);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVASVAngle2);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVASVAngle3);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVASVAngle3);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVASVAngle3);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarSV);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarSV);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarSV);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarSVAngle0);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarSVAngle0);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarSVAngle0);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarSVAngle1);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarSVAngle1);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarSVAngle1);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarSVAngle2);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarSVAngle2);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarSVAngle2);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarSVAngle3);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarSVAngle3);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarSVAngle3);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStar);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStar);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStar);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAFiner);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAFiner);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAFiner);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVACoarser);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVACoarser);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVACoarser);
  
  Extradist1d.push_back(&MuOSSS1InvariantMassBeforeMVA);
  Extradist1d.push_back(&MuOSSS2InvariantMassBeforeMVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdR);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdR);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdR);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdR);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdR);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdR);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRIncrease);
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRDecrease);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd1);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd1);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd1);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd2);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd2);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd2);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd3);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd3);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd3);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd4);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd4);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd4);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd5);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd5);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd5);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd6);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd6);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd6);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestdRLtd7);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestdRLtd7);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestdRLtd7);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7);
  
  Extradist1d.push_back(&IsolationTrackCount);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVABestMass);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVABestMass);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVABestMass);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassBeforeMVAKStarBestMass);
  Extradist1d.push_back(&Mu2TrackInvariantMassBeforeMVAKStarBestMass);
  Extradist1d.push_back(&Mu3TrackInvariantMassBeforeMVAKStarBestMass);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterA1MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterA1MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterA1MVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterB1MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterB1MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterB1MVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterC1MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterC1MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterC1MVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterA2MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterA2MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterA2MVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterB2MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterB2MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterB2MVA);
  
  Extradist1d.push_back(&Mu1TrackInvariantMassAfterC2MVA);
  Extradist1d.push_back(&Mu2TrackInvariantMassAfterC2MVA);
  Extradist1d.push_back(&Mu3TrackInvariantMassAfterC2MVA);

  Extradist1d.push_back(&TauMassRefitABC1);  
  Extradist1d.push_back(&TauMassRefitABC2);

  Extradist1d.push_back(&TauMassRefitA1);
  Extradist1d.push_back(&TauMassRefitB1);
  Extradist1d.push_back(&TauMassRefitC1);
  Extradist1d.push_back(&TauMassRefitA2);
  Extradist1d.push_back(&TauMassRefitB2);
  Extradist1d.push_back(&TauMassRefitC2);





  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist2d.push_back(&TauMass_all_nophiVeto);
  Extradist1d.push_back(&TauMass_all);
  Extradist2d.push_back(&TauMass_allVsBDTA);
  Extradist2d.push_back(&TauMass_allVsBDTB);
  Extradist2d.push_back(&TauMass_allVsBDTC);

  Extradist2d.push_back(&EMR_tau_eta);



  Extradist1d.push_back(&SVPVTauDirAngle);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);


  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);


  // Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&EventMassResolution_PtEtaPhi);
  Extradist1d.push_back(&NSignalCandidates);


  Extradist1d.push_back(&BDTOutputA);
  Extradist1d.push_back(&BDTOutputB);
  Extradist1d.push_back(&BDTOutputC);


  Extradist2d.push_back(&PairMass);
  Extradist2d.push_back(&PairMassFinalSel);
  Extradist1d.push_back(&PairMass1);
  Extradist1d.push_back(&PairMass2);
  Extradist2d.push_back(&PairMassWithCut);
  Extradist2d.push_back(&PairMassEta);
  Extradist2d.push_back(&PairMassEtaPrime);



  Extradist1d.push_back(&Muon1MVAID);
  Extradist1d.push_back(&Muon2MVAID);
  Extradist1d.push_back(&Muon3MVAID);

}


void  CommonSelector::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  //    std::cout<<" id   "<< id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}


  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu0Fired(false);
  bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);

  random_num = rndm.Rndm();


  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    if(HLTName.Contains("DoubleMu3_TkMu_DsTau3Mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}
  }

  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;

  if (HLTOk) value.at(HLT) = true;
  else value.at(HLT) = false;

  if (L1Ok) value.at(L1T) = true;
  else value.at(L1T) = false;


  if(DoubleMuFired) value.at(L1T)=1;

  pass.at(L1T)= (value.at(L1T)==cut.at(L1T));
  pass.at(HLT)= (value.at(HLT)==cut.at(HLT));



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


  NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);
  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();

    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);

    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);
    //*** Muon ID cut
    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) &&  Ntp->Muon_isPFMuon(mu1_pt_idx) &&
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&  Ntp->Muon_isPFMuon(mu2_pt_idx) &&
			Ntp->Muon_isGlobalMuon(mu3_pt_idx) &&  Ntp->Muon_isPFMuon(mu3_pt_idx));



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

    value.at(PhiVeto1) =  M_osss1;//fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(PhiVeto2) =  M_osss2;//fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto1) = M_osss1;//fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;
    value.at(OmegaVeto2) = M_osss2;//fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;


    vector<TLorentzVector> trigobjTriplet;
    for (int i=0; i<Ntp->NTriggerObjects(); i++){
      TString name = Ntp->TriggerObject_name(i);
      if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
      trigobjTriplet.push_back(tmp);
    }

    std::vector<TLorentzVector> muonTriplet;
    muonTriplet.push_back(Ntp->Muon_P4(mu1_pt_idx));
    muonTriplet.push_back(Ntp->Muon_P4(mu2_pt_idx));
    muonTriplet.push_back(Ntp->Muon_P4(mu3_pt_idx));

    bool triggerCheck = false;
    if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
    value.at(TriggerMatch) = triggerCheck;

    value.at(TauMassCut) = TauRefittedLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));
  pass.at(PhiVeto1) = true;//(value.at(PhiVeto1) < 0.98 || value.at(PhiVeto1) > 1.06 );
  pass.at(OmegaVeto1) = true;//(value.at(OmegaVeto1) < 0.742 || value.at(OmegaVeto1) > 0.822 );
  pass.at(PhiVeto2) = true;//(value.at(PhiVeto2) < 0.98 || value.at(PhiVeto2) > 1.06 );
  pass.at(OmegaVeto2) = true;//(value.at(OmegaVeto2) < 0.742 || value.at(OmegaVeto2) > 0.822 );
  pass.at(TauMassCut) =( (value.at(TauMassCut) > tauMinSideBand_)  &&   (value.at(TauMassCut) < tauMaxSideBand_ ));

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(PhiVeto1);
  exclude_cuts.push_back(OmegaVeto1);


  exclude_cuts.push_back(PhiVeto2);
  exclude_cuts.push_back(OmegaVeto2);

  if(passAllBut(exclude_cuts)){


    TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);

    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
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


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);

  if(status){

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

    EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);

    //*** Rapidity sorted muons
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

    //*** With such sorting pT(ss1) > pT(ss2)
    
    TLorentzVector MuonOS = Ntp->Muon_P4(os_mu_idx);  
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);
    
    TLorentzVector TauAngleHalf1=MuonOS+MuonSS1;
    TLorentzVector TauAngleHalf2=MuonSS2;
    TLorentzVector TauCombination=TauAngleHalf1+TauAngleHalf2;
    TauCombination.SetE(sqrt(TauCombination.Px()*TauCombination.Px()+TauCombination.Py()*TauCombination.Py()+TauCombination.Pz()*TauCombination.Pz()+1.77686*1.77686));
    TauAngleHalf1.Boost(-1*TauCombination.BoostVector());
    TauAngleHalf2.Boost(-1*TauCombination.BoostVector());
    
    if(fabs(TauAngleHalf1.Vect().Angle(TauAngleHalf2.Vect())-3.141592653589793)<0.02){
      TauAngleTest.at(t).Fill((MuonOS+MuonSS1+MuonSS2).M(),1);
    }
    
    
    
    //TLorentzVector MuonOS = Muon1LV;
    
    TLorentzVector MuonOSReassigned = MuonOS;//reassign the masses to be similar to a kaon
    MuonOSReassigned.SetE(sqrt(MuonOS.Px()*MuonOS.Px()+MuonOS.Py()*MuonOS.Py()+MuonOS.Pz()*MuonOS.Pz()+0.493677*0.493677));
    TLorentzVector MuonSS1Reassigned = MuonSS1;
    MuonSS1Reassigned.SetE(sqrt(MuonSS1.Px()*MuonSS1.Px()+MuonSS1.Py()*MuonSS1.Py()+MuonSS1.Pz()*MuonSS1.Pz()+0.493677*0.493677));
    TLorentzVector MuonSS2Reassigned = MuonSS2;
    MuonSS2Reassigned.SetE(sqrt(MuonSS2.Px()*MuonSS2.Px()+MuonSS2.Py()*MuonSS2.Py()+MuonSS2.Pz()*MuonSS2.Pz()+0.493677*0.493677));
    
    MuOSSS1InvariantMassBeforeMVA.at(t).Fill((MuonOSReassigned+MuonSS1Reassigned).M(),1);
    MuOSSS2InvariantMassBeforeMVA.at(t).Fill((MuonOSReassigned+MuonSS2Reassigned).M(),1);
    
    double dR1 = 99.0;
    double dR2 = 99.0;
    double dR3 = 99.0;
    
    double M1dR = 99.0;
    double M2dR = 99.0;
    double M3dR = 99.0;
    
    double M1best = 99.0;
    double M2best = 99.0;
    double M3best = 99.0;
    
    double KdR1 = 99.0;
    double KdR2 = 99.0;
    double KdR3 = 99.0;
    
    double KM1dR = 99.0;
    double KM2dR = 99.0;
    double KM3dR = 99.0;
    
    double KM1best = 99.0;
    double KM2best = 99.0;
    double KM3best = 99.0;
    
    double OSdRSV=99.0*sqrt(2);
    double SS1dRSV=99.0*sqrt(2);
    double SS2dRSV=99.0*sqrt(2);
    
    double M1bestSV = 99.0;
    double M2bestSV = 99.0;
    double M3bestSV = 99.0;
    
    double M1bestSVKStar = 99.0;
    double M2bestSVKStar = 99.0;
    double M3bestSVKStar = 99.0;
    
    double bestpTSV1 = -1.0;
    double bestpTSV2 = -1.0;
    double bestpTSV3 = -1.0;
    
    double bestpTSVKStar1 = -1.0;
    double bestpTSVKStar2 = -1.0;
    double bestpTSVKStar3 = -1.0;
    
    double bestpSV1 = -1.0;
    double bestpSV2 = -1.0;
    double bestpSV3 = -1.0;
    
    double bestpSVKStar1 = -1.0;
    double bestpSVKStar2 = -1.0;
    double bestpSVKStar3 = -1.0;
    
    double BestOSSVAngle = 99.0;
    double BestSS1SVAngle = 99.0;
    double BestSS2SVAngle = 99.0;
        
    double BestOSSVAngleMass = 99.0;
    double BestSS1SVAngleMass = 99.0;
    double BestSS2SVAngleMass = 99.0;
    
    double BestOSSVAngleKStar = 99.0;
    double BestSS1SVAngleKStar = 99.0;
    double BestSS2SVAngleKStar = 99.0;
        
    double BestOSSVAngleMassKStar = 99.0;
    double BestSS1SVAngleMassKStar = 99.0;
    double BestSS2SVAngleMassKStar = 99.0;
    
    double bestOSAng = 99.0;
    double bestOSAngM = 99.0;
    double bestOSAngMK = 99.0;
    
    double bestSS1Ang = 99.0;
    double bestSS1AngM = 99.0;
    double bestSS1AngMK = 99.0;
    
    double bestSS2Ang = 99.0;
    double bestSS2AngM = 99.0;
    double bestSS2AngMK = 99.0;
    
    double M1bestSVQual1 = 99.0;
    double M2bestSVQual1 = 99.0;
    double M3bestSVQual1 = 99.0;
    
    double Mu1bestChi2 = 199.0;
    double Mu2bestChi2 = 199.0;
    double Mu3bestChi2 = 199.0;
    
    double M1bestChi2 = 99.0;
    double M2bestChi2 = 99.0;
    double M3bestChi2 = 99.0;
    
    double M1KStarbestChi2 = 99.0;
    double M2KStarbestChi2 = 99.0;
    double M3KStarbestChi2 = 99.0;
    
    double Mu1bestChi2AngleDiff = 99.0;//angle between SVPV and muon-track combination
    double Mu2bestChi2AngleDiff = 99.0;
    double Mu3bestChi2AngleDiff = 99.0;
    
    double Mu1bestChi2DistRatio = 99.0;//ratio of SVPV dist and c t gamma
    double Mu2bestChi2DistRatio = 99.0;
    double Mu3bestChi2DistRatio = 99.0;
    
    double Mu1KSbestChi2DistRatio = 99.0;//ratio of SVPV dist and c t gamma
    double Mu2KSbestChi2DistRatio = 99.0;
    double Mu3KSbestChi2DistRatio = 99.0;
    
    for(int i=0;i<Ntp->NIsolationTrack(signal_idx);i++){
      for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
        if(Ntp->IsolationTrack_VertexWithSignalMuon1IsValid(signal_idx,i)&&Ntp->IsolationTrack_VertexWithSignalMuon1IsValid(signal_idx,j)){
          Mu1TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx,i)+Ntp->IsolationTrack_p4(signal_idx,j)).M(),1);
        }
        if(Ntp->IsolationTrack_VertexWithSignalMuon2IsValid(signal_idx,i)&&Ntp->IsolationTrack_VertexWithSignalMuon2IsValid(signal_idx,j)){
          Mu2TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx,i)+Ntp->IsolationTrack_p4(signal_idx,j)).M(),1);
        }
        if(Ntp->IsolationTrack_VertexWithSignalMuon3IsValid(signal_idx,i)&&Ntp->IsolationTrack_VertexWithSignalMuon3IsValid(signal_idx,j)){
          Mu3TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx,i)+Ntp->IsolationTrack_p4(signal_idx,j)).M(),1);
        }
      }
    }
    
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      TLorentzVector CombOSPhi = ParticleLVReassigned+MuonOSReassigned;
      TLorentzVector CombOSKS = ParticleLV+MuonOSReassigned;
      
      TLorentzVector CombSS1Phi = ParticleLVReassigned+MuonSS1Reassigned;
      TLorentzVector CombSS1KS = ParticleLV+MuonSS1Reassigned;
      
      TLorentzVector CombSS2Phi = ParticleLVReassigned+MuonSS2Reassigned;
      TLorentzVector CombSS2KS = ParticleLV+MuonSS2Reassigned;
      
      //MassDifference.at(t).Fill(ParticleLV.M()-0.493677,1);
      
      //std::cout<<"The charge of OS is "<<Ntp->Muon_charge(os_mu_idx)<<" SS1 is "<<Ntp->Muon_charge(ss1_mu_idx)<<" SS2 is "<<Ntp->Muon_charge(ss2_mu_idx)<< std::endl;
      //std::cout<<"The charge of Particle is "<<Ntp->IsolationTrack_charge(signal_idx,j)<<std::endl;
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))&&Ntp->IsolationTrack_VertexWithSignalMuon1IsValid(signal_idx,j)){
        Mu1TrackInvariantMassBeforeMVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
        Mu1TrackInvariantMassBeforeMVAFiner.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
        Mu1TrackInvariantMassBeforeMVACoarser.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
        
        TVector3 SVPVVect =  Ntp->IsolationTrack_VertexWithSignalMuon1Position(signal_idx,j)-Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
        TVector3 CombDirOS = CombOSPhi.Vect();
        
        Mu1bestChi2DistRatio =  (4.6345*pow(10,-12)*CombOSPhi.Gamma())/(SVPVVect.Mag());
        Mu1KSbestChi2DistRatio =  (4.1594*pow(10,-13)*CombOSKS.Gamma())/(SVPVVect.Mag());
        
        if(Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,j)<Mu1bestChi2){
          Mu1bestChi2=Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,j);
          M1bestChi2=(ParticleLVReassigned+MuonOSReassigned).M();
          M1KStarbestChi2=(ParticleLV+MuonOSReassigned).M();
          
          Mu1bestChi2AngleDiff=fabs(SVPVVect.Angle(CombDirOS));
          
          Mu1bestChi2DistRatio =  (4.6345*pow(10,-12)*CombOSPhi.Gamma())/(SVPVVect.Mag());
          Mu1KSbestChi2DistRatio =  (4.1594*pow(10,-13)*CombOSKS.Gamma())/(SVPVVect.Mag());
        }
        
        Mu1TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,j),1);
        
        
        
        if(fabs(ParticleLVReassigned.DeltaR(MuonOSReassigned))<dR1){//try to find the pair with the lowest dR
          dR1=fabs(ParticleLVReassigned.DeltaR(MuonOSReassigned));
          M1dR=(ParticleLVReassigned+MuonOSReassigned).M();
        }
        
        if(fabs((ParticleLVReassigned+MuonOSReassigned).M()-1.019461)<fabs(M1best-1.019461)){//find the pair that is closest in mass
          M1best=(ParticleLVReassigned+MuonOSReassigned).M();
        }
        
        Mu1TrackInvariantMassBeforeMVAKStar.at(t).Fill((ParticleLV+MuonOSReassigned).M(),1);//For the Kstars
        
        if(fabs(ParticleLV.DeltaR(MuonOSReassigned))<KdR1){//try to find the pair with the lowest dR
          KdR1=fabs(ParticleLV.DeltaR(MuonOSReassigned));
          KM1dR=(ParticleLV+MuonOSReassigned).M();
        }
        
        if(fabs((ParticleLV+MuonOSReassigned).M()-0.892)<fabs(KM1best-0.892)){//find the pair that is closest in mass
          KM1best=(ParticleLV+MuonOSReassigned).M();
        }
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))&&Ntp->IsolationTrack_VertexWithSignalMuon2IsValid(signal_idx,j)){
        Mu2TrackInvariantMassBeforeMVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
        Mu2TrackInvariantMassBeforeMVAFiner.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
        Mu2TrackInvariantMassBeforeMVACoarser.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
        
        TVector3 SVPVVect =  Ntp->IsolationTrack_VertexWithSignalMuon2Position(signal_idx,j)-Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
        TVector3 CombDirSS1 = CombSS1Phi.Vect();
        
        if(Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,j)<Mu2bestChi2){
          Mu2bestChi2=Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,j);
          M2bestChi2=(ParticleLVReassigned+MuonSS1Reassigned).M();
          M2KStarbestChi2=(ParticleLV+MuonSS1Reassigned).M();
          
          Mu2bestChi2AngleDiff=fabs(SVPVVect.Angle(CombDirSS1));
          
          Mu2bestChi2DistRatio =  (4.6345*pow(10,-12)*CombSS1Phi.Gamma())/(SVPVVect.Mag());
          Mu2KSbestChi2DistRatio =  (4.1594*pow(10,-13)*CombSS1KS.Gamma())/(SVPVVect.Mag());
        }
        
        Mu2TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,j),1);
        
        if(fabs(ParticleLVReassigned.DeltaR(MuonSS1Reassigned))<dR2){//try to find the pair with the lowest dR
          dR2=fabs(ParticleLVReassigned.DeltaR(MuonSS1Reassigned));
          M2dR=(ParticleLVReassigned+MuonSS1Reassigned).M();
        }
        
        if(fabs((ParticleLVReassigned+MuonSS1Reassigned).M()-1.019461)<fabs(M2best-1.019461)){//find the pair that is closest in mass
          M2best=(ParticleLVReassigned+MuonSS1Reassigned).M();
        }
        
        Mu2TrackInvariantMassBeforeMVAKStar.at(t).Fill((ParticleLV+MuonSS1Reassigned).M(),1);
        
        if(fabs(ParticleLV.DeltaR(MuonSS1Reassigned))<KdR2){//try to find the pair with the lowest dR
          KdR2=fabs(ParticleLV.DeltaR(MuonSS1Reassigned));
          KM2dR=(ParticleLV+MuonSS1Reassigned).M();
        }
        
        if(fabs((ParticleLV+MuonSS1Reassigned).M()-0.892)<fabs(KM2best-0.892)){//find the pair that is closest in mass
          KM2best=(ParticleLV+MuonSS1Reassigned).M();
        }
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))&&Ntp->IsolationTrack_VertexWithSignalMuon3IsValid(signal_idx,j)){
        Mu3TrackInvariantMassBeforeMVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
        Mu3TrackInvariantMassBeforeMVAFiner.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
        Mu3TrackInvariantMassBeforeMVACoarser.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
        
        TVector3 SVPVVect =  Ntp->IsolationTrack_VertexWithSignalMuon3Position(signal_idx,j)-Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
        TVector3 CombDirSS2 = CombSS2Phi.Vect();
        
        if(Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,j)<Mu3bestChi2){
          Mu3bestChi2=Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,j);
          M3bestChi2=(ParticleLVReassigned+MuonSS2Reassigned).M();
          M3KStarbestChi2=(ParticleLV+MuonSS2Reassigned).M();
          
          Mu3bestChi2AngleDiff=fabs(SVPVVect.Angle(CombDirSS2));
          
          Mu3bestChi2DistRatio =  (4.6345*pow(10,-12)*CombSS2Phi.Gamma())/(SVPVVect.Mag());
          Mu3KSbestChi2DistRatio =  (4.1594*pow(10,-13)*CombSS2KS.Gamma())/(SVPVVect.Mag());
        }
        
        Mu3TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,j),1);
        
        if(fabs(ParticleLVReassigned.DeltaR(MuonSS2Reassigned))<dR3){//try to find the pair with the lowest dR
          dR3=fabs(ParticleLVReassigned.DeltaR(MuonSS2Reassigned));
          M3dR=(ParticleLVReassigned+MuonSS2Reassigned).M();
        }
        
        if(fabs((ParticleLVReassigned+MuonSS2Reassigned).M()-1.019461)<fabs(M3best-1.019461)){//find the pair that is closest in mass
          M3best=(ParticleLVReassigned+MuonSS2Reassigned).M();
        }
        
        Mu3TrackInvariantMassBeforeMVAKStar.at(t).Fill((ParticleLV+MuonSS2Reassigned).M(),1);
        
        if(fabs(ParticleLV.DeltaR(MuonSS2Reassigned))<KdR3){//try to find the pair with the lowest dR
          KdR3=fabs(ParticleLV.DeltaR(MuonSS2Reassigned));
          KM3dR=(ParticleLV+MuonSS2Reassigned).M();
        }
        
        if(fabs((ParticleLV+MuonSS2Reassigned).M()-0.892)<fabs(KM3best-0.892)){//find the pair that is closest in mass
          KM3best=(ParticleLV+MuonSS2Reassigned).M();
        }
      }
      
      
      /*
      if(dRmu1<dRmu2){
        if(dRmu1<dRmu3){
          Mu1TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
          if(dRmu2<dRmu3){
            Mu2TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
          }
          else {
            Mu3TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
          }
        }
        else {
          Mu1TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
          if(dRmu1<dRmu2){
            Mu2TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
          }
          else {
            Mu3TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
          }
        }
      }
      else {
        if(dRmu2<dRmu3){
          Mu1TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
          if(dRmu1<dRmu3){
            Mu2TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
          }
          else {
            Mu3TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
          }
        }
        else {
          Mu1TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
          if(dRmu1<dRmu2){
            Mu2TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
          }
          else {
            Mu3TrackInvariantMassAfterMVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
          }
        }
      }
      */
      
      
      //MassDifference.at(t).Fill((ParticleLV+Muon2LV).M(),1);
      //MassDifference.at(t).Fill((ParticleLV+Muon2LV).M(),1);
      //std::cout<<"Mass of particle "<<j<<" is "<<ParticleLV.M()<<" GeV"<< std::endl;
      for(int m=0;m<Ntp->NSecondaryVertices();m++){
        
        double OSdRSV1=99.0*sqrt(2);
        double OSdRSV1mu=99.0;
        double OSdRSV1track=99.0;
        
        double SS1dRSV1=99.0*sqrt(2);
        double SS1dRSV1mu=99.0;
        double SS1dRSV1track=99.0;
        
        double SS2dRSV1=99.0*sqrt(2);
        double SS2dRSV1mu=99.0;
        double SS2dRSV1track=99.0;
        
        double M1bestSV1 = 99.0;
        double M2bestSV1 = 99.0;
        double M3bestSV1 = 99.0;
        
        
        TLorentzVector MuLVSV1(0.,0.,0.,0.);
        TLorentzVector PtLVSV1(0.,0.,0.,0.);
        TLorentzVector MuLVSV2(0.,0.,0.,0.);
        TLorentzVector PtLVSV2(0.,0.,0.,0.);
        TLorentzVector MuLVSV3(0.,0.,0.,0.);
        TLorentzVector PtLVSV3(0.,0.,0.,0.);
        
        
        for(int n=0;n<Ntp->NTracksAtSecondaryVertex(m);n++){//for each MuLV, ParticleLV pair, try to find the tracks with common SV that closely match
          
          //std::cout<<"The mass of the track is "<< Ntp->SecondaryVertexTrack_P4(m,n).M() <<"GeV"<<std::endl;
          
         if(Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->Muon_charge(os_mu_idx)&&sqrt(pow(MuonOSReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+OSdRSV1track*OSdRSV1track)<OSdRSV1){
           OSdRSV1=sqrt(pow(MuonOSReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+OSdRSV1track*OSdRSV1track);
           OSdRSV1mu=MuonOSReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           MuLVSV1=Ntp->SecondaryVertexTrack_P4(m,n);
           MuLVSV1.SetE(sqrt(MuLVSV1.Px()*MuLVSV1.Px()+MuLVSV1.Py()*MuLVSV1.Py()+MuLVSV1.Pz()*MuLVSV1.Pz()+0.493677*0.493677));
           M1bestSV1=(MuLVSV1+PtLVSV1).M();
         }
         if(Ntp->SecondaryVertexTrackCharge(m,n)==(-1*Ntp->Muon_charge(os_mu_idx))&&Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->IsolationTrack_charge(signal_idx,j)&&sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+OSdRSV1mu*OSdRSV1mu)<OSdRSV1){//only tracks with opposite charge of muon consiered
           OSdRSV1=sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+OSdRSV1mu*OSdRSV1mu);
           OSdRSV1track=ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           PtLVSV1=Ntp->SecondaryVertexTrack_P4(m,n);
           PtLVSV1.SetE(sqrt(PtLVSV1.Px()*PtLVSV1.Px()+PtLVSV1.Py()*PtLVSV1.Py()+PtLVSV1.Pz()*PtLVSV1.Pz()+0.493677*0.493677));
           M1bestSV1=(MuLVSV1+PtLVSV1).M();
         }
         
         
         
         if(Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->Muon_charge(ss1_mu_idx)&&sqrt(pow(MuonSS1Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS1dRSV1track*SS1dRSV1track)<SS1dRSV1){
           SS1dRSV1=sqrt(pow(MuonSS1Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS1dRSV1track*SS1dRSV1track);
           SS1dRSV1mu=MuonSS1Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           MuLVSV2=Ntp->SecondaryVertexTrack_P4(m,n);
           MuLVSV2.SetE(sqrt(MuLVSV2.Px()*MuLVSV2.Px()+MuLVSV2.Py()*MuLVSV2.Py()+MuLVSV2.Pz()*MuLVSV2.Pz()+0.493677*0.493677));
           M2bestSV1=(MuLVSV2+PtLVSV2).M();
         }
         if(Ntp->SecondaryVertexTrackCharge(m,n)==(-1*Ntp->Muon_charge(ss1_mu_idx))&&Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->IsolationTrack_charge(signal_idx,j)&&sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS1dRSV1mu*SS1dRSV1mu)<SS1dRSV1){
           SS1dRSV1=sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS1dRSV1mu*SS1dRSV1mu);
           SS1dRSV1track=ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           PtLVSV2=Ntp->SecondaryVertexTrack_P4(m,n);
           PtLVSV2.SetE(sqrt(PtLVSV2.Px()*PtLVSV2.Px()+PtLVSV2.Py()*PtLVSV2.Py()+PtLVSV2.Pz()*PtLVSV2.Pz()+0.493677*0.493677));
           M2bestSV1=(MuLVSV2+PtLVSV2).M();
         }
         
         
         
         if(Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->Muon_charge(ss2_mu_idx)&&sqrt(pow(MuonSS2Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS2dRSV1track*SS2dRSV1track)<SS2dRSV1){
           SS2dRSV1=sqrt(pow(MuonSS2Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS2dRSV1track*SS2dRSV1track);
           SS2dRSV1mu=MuonSS2Reassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           MuLVSV3=Ntp->SecondaryVertexTrack_P4(m,n);
           MuLVSV3.SetE(sqrt(MuLVSV3.Px()*MuLVSV3.Px()+MuLVSV3.Py()*MuLVSV3.Py()+MuLVSV3.Pz()*MuLVSV3.Pz()+0.493677*0.493677));
           M3bestSV1=(MuLVSV3+PtLVSV3).M();
         }
         if(Ntp->SecondaryVertexTrackCharge(m,n)==(-1*Ntp->Muon_charge(ss2_mu_idx))&&Ntp->SecondaryVertexTrackCharge(m,n)==Ntp->IsolationTrack_charge(signal_idx,j)&&sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS2dRSV1mu*SS2dRSV1mu)<SS2dRSV1){
           SS2dRSV1=sqrt(pow(ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n)),2)+SS2dRSV1mu*SS2dRSV1mu);
           SS2dRSV1track=ParticleLVReassigned.DeltaR(Ntp->SecondaryVertexTrack_P4(m,n));
           PtLVSV3=Ntp->SecondaryVertexTrack_P4(m,n);
           PtLVSV3.SetE(sqrt(PtLVSV3.Px()*PtLVSV3.Px()+PtLVSV3.Py()*PtLVSV3.Py()+PtLVSV3.Pz()*PtLVSV3.Pz()+0.493677*0.493677));
           M3bestSV1=(MuLVSV3+PtLVSV3).M();
         }
         
         
         
        }//end of n loop
        
        TVector3 SVPVVect =  Ntp->SecondaryVertexPosition(m)-Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
        TVector3 CombDirOS = (MuonOSReassigned+ParticleLV).Vect();
        TVector3 CombDirSS1 = (MuonSS1Reassigned+ParticleLV).Vect();
        TVector3 CombDirSS2 = (MuonSS2Reassigned+ParticleLV).Vect();
        
        if(OSdRSV1<OSdRSV){
          OSdRSV=OSdRSV1;
          M1bestSV=M1bestSV1;
          M1bestSVKStar=(MuonOSReassigned+ParticleLV).M();
          
          bestpTSV1=(MuonOSReassigned+ParticleLVReassigned).Pt();
          bestpTSVKStar1=(MuonOSReassigned+ParticleLV).Pt();
          
          bestpSV1=(MuonOSReassigned+ParticleLVReassigned).P();
          bestpSVKStar1=(MuonOSReassigned+ParticleLV).P();
          
        }
        
               
        if(SS1dRSV1<SS1dRSV){
          SS1dRSV=SS1dRSV1;
          M2bestSV=M2bestSV1;
          M2bestSVKStar=(MuonSS1Reassigned+ParticleLV).M();
          
          bestpTSV2=(MuonSS1Reassigned+ParticleLVReassigned).Pt();
          bestpTSVKStar2=(MuonSS1Reassigned+ParticleLV).Pt();
          
          bestpSV2=(MuonSS1Reassigned+ParticleLVReassigned).P();
          bestpSVKStar2=(MuonSS1Reassigned+ParticleLV).P();
        }
        
        if(SS2dRSV1<SS2dRSV){
          SS2dRSV=SS2dRSV1;
          M3bestSV=M3bestSV1;
          M3bestSVKStar=(MuonSS2Reassigned+ParticleLV).M();
          
          bestpTSV3=(MuonSS2Reassigned+ParticleLVReassigned).Pt();
          bestpTSVKStar3=(MuonSS2Reassigned+ParticleLV).Pt();
          
          bestpSV3=(MuonSS2Reassigned+ParticleLVReassigned).P();
          bestpSVKStar3=(MuonSS2Reassigned+ParticleLV).P();
        }
        
        if(OSdRSV1<0.00001){
        }
        
        
        if(OSdRSV1<0.00001){
          //Mu1TrackInvariantMassBeforeMVASV.at(t).Fill(M1bestSV,1);
          //Mu1TrackInvariantMassBeforeMVAKStarSV.at(t).Fill((MuonOSReassigned+ParticleLV).M(),1);
          
          if(fabs(SVPVVect.Angle(CombDirOS))<bestOSAng){
            bestOSAng=fabs(SVPVVect.Angle(CombDirOS));
            bestOSAngM=(MuonOSReassigned+ParticleLVReassigned).M();
            bestOSAngMK=(MuonOSReassigned+ParticleLV).M();
          }
          
          TLorentzVector MuTest=MuonOSReassigned;
          TLorentzVector PtTest=ParticleLVReassigned;
          TLorentzVector SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+1.019461*1.019461));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestOSSVAngle){//0.03
            BestOSSVAngle=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestOSSVAngleMass=M1bestSV1;
          }
          
          
          
          
          MuTest=MuonOSReassigned;//For K star 896, neutral
          PtTest=ParticleLV;
          SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+0.896*0.896));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestOSSVAngleKStar){//0.03
            BestOSSVAngleKStar=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestOSSVAngleMassKStar=(MuonOSReassigned+ParticleLV).M();
          }
        }
        
        
        if(SS1dRSV1<0.00001){
          //Mu2TrackInvariantMassBeforeMVASV.at(t).Fill(M2bestSV,1);
          //Mu2TrackInvariantMassBeforeMVAKStarSV.at(t).Fill((MuonSS1Reassigned+ParticleLV).M(),1);
          
          if(fabs(SVPVVect.Angle(CombDirSS1))<bestSS1Ang){
            bestSS1Ang=fabs(SVPVVect.Angle(CombDirSS1));
            bestSS1AngM=(MuonSS1Reassigned+ParticleLVReassigned).M();
            bestSS1AngMK=(MuonSS1Reassigned+ParticleLV).M();
          }
          
          TLorentzVector MuTest=MuonSS1Reassigned;
          TLorentzVector PtTest=ParticleLVReassigned;
          TLorentzVector SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+1.019461*1.019461));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestSS1SVAngle){//0.03
            BestSS1SVAngle=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestSS1SVAngleMass=M2bestSV1;
          }
          
          
          
          MuTest=MuonSS1Reassigned;//For K star 896, neutral
          PtTest=ParticleLV;
          SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+0.896*0.896));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestSS1SVAngleKStar){//0.03
            BestSS1SVAngleKStar=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestSS1SVAngleMassKStar=(MuonSS1Reassigned+ParticleLV).M();
          }
        }
        
        
        
        if(SS2dRSV1<0.00001){
        
          //std::cout<<"The dR is, for SS2dRSV1, : "<< SS2dRSV1 <<std::endl; it's around 10^-8
          //Mu3TrackInvariantMassBeforeMVASV.at(t).Fill(M3bestSV,1);
          //Mu3TrackInvariantMassBeforeMVAKStarSV.at(t).Fill((MuonSS2Reassigned+ParticleLV).M(),1);
          
          if(fabs(SVPVVect.Angle(CombDirSS2))<bestSS2Ang){
            bestSS2Ang=fabs(SVPVVect.Angle(CombDirSS2));
            bestSS2AngM=(MuonSS2Reassigned+ParticleLVReassigned).M();
            bestSS2AngMK=(MuonSS2Reassigned+ParticleLV).M();
          }
          
          TLorentzVector MuTest=MuonSS2Reassigned;
          TLorentzVector PtTest=ParticleLVReassigned;
          TLorentzVector SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+1.019461*1.019461));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          Mu1TrackInvariantMassBeforeMVABestdRDecrease.at(t).Fill(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793),1);
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestSS2SVAngle){//0.03
            BestSS2SVAngle=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestSS2SVAngleMass=M3bestSV1;
          }
          
          
          
          
          Mu1TrackInvariantMassBeforeMVABestdRIncrease.at(t).Fill(MuLVSV3.DeltaR(PtLVSV3),1);
          
          MuTest=MuonSS2Reassigned;//For K star 896, neutral
          PtTest=ParticleLV;
          SumLV=MuTest+PtTest;
          SumLV.SetE(sqrt(SumLV.Px()*SumLV.Px()+SumLV.Py()*SumLV.Py()+SumLV.Pz()*SumLV.Pz()+0.896*0.896));
          MuTest.Boost(-1*SumLV.BoostVector());
          PtTest.Boost(-1*SumLV.BoostVector());
          
          if(fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793)<BestSS2SVAngleKStar){//0.03
            BestSS2SVAngleKStar=fabs(MuTest.Vect().Angle(PtTest.Vect())-3.141592653589793);
            BestSS2SVAngleMassKStar=(MuonSS2Reassigned+ParticleLV).M();
          }
        }
        
        
      }//end of m loop
    }
    
    
    /*
    if(OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVASV.at(t).Fill(bestOSAngM,1);
      Mu1TrackInvariantMassBeforeMVAKStarSV.at(t).Fill(bestOSAngMK,1);
      
      if(M1bestSV<1.03&&M1bestSV>1.01){
        Mu1TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(bestpTSV1,1);
        Mu1TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(bestOSAng,1);
      }
      if(M1bestSVKStar<0.95&&M1bestSVKStar>0.8){
        Mu1TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(bestpTSVKStar1,1);
        Mu1TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(bestOSAng,1);
      }
    }
    
    if(SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVASV.at(t).Fill(bestSS1AngM,1);
      Mu2TrackInvariantMassBeforeMVAKStarSV.at(t).Fill(bestSS1AngMK,1);
      
      if(M2bestSV<1.03&&M2bestSV>1.01){
        Mu2TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(bestpTSV2,1);
        Mu2TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(bestSS1Ang,1);
      }
      if(M2bestSVKStar<0.95&&M2bestSVKStar>0.8){
        Mu2TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(bestpTSVKStar2,1);
        Mu2TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(bestSS1Ang,1);
      }
    }
    
    if(SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVASV.at(t).Fill(bestSS2AngM,1);
      Mu3TrackInvariantMassBeforeMVAKStarSV.at(t).Fill(bestSS2AngMK,1);
      
      if(M3bestSV<1.03&&M3bestSV>1.01){
        Mu3TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(bestpTSV3,1);
        Mu3TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(bestSS2Ang,1);
      }
      if(M3bestSVKStar<0.95&&M3bestSVKStar>0.8){
        Mu3TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(bestpTSVKStar3,1);
        Mu3TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(bestSS2Ang,1);
      }
    }
    */
    
    
    
    
    
    IsolationTrackCount.at(t).Fill(Ntp->NIsolationTrack(signal_idx),1);
    
    Mu1TrackInvariantMassBeforeMVABestdR.at(t).Fill(M1dR,1);
    Mu2TrackInvariantMassBeforeMVABestdR.at(t).Fill(M2dR,1);
    Mu3TrackInvariantMassBeforeMVABestdR.at(t).Fill(M3dR,1);
    
    if(bestOSAng<0.02&&bestpTSV1>6){
      Mu1TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(bestOSAngM,1);
    }
    if(bestSS1Ang<0.02&&bestpTSV2>6){
      Mu2TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(bestSS1AngM,1);
    }
    if(bestSS2Ang<0.02&&bestpTSV3>6){
      Mu3TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(bestSS2AngM,1);
    }
    
    if(bestOSAng<0.05&&bestpTSV1>6){
      Mu1TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(bestOSAngM,1);
    }
    if(bestSS1Ang<0.05&&bestpTSV2>6){
      Mu2TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(bestSS1AngM,1);
    }
    if(bestSS2Ang<0.05&&bestpTSV3>6){
      Mu3TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(bestSS2AngM,1);
    }
    
    //KStar angles
    
    if(bestOSAng<0.02&&bestpTSV1>6){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(bestOSAngMK,1);
    }
    if(bestSS1Ang<0.02&&bestpTSV2>6){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(bestSS1AngMK,1);
    }
    if(bestSS2Ang<0.02&&bestpTSV3>6){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(bestSS2AngMK,1);
    }
    
    if(bestOSAng<0.05&&bestpTSV1>6){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(bestOSAngMK,1);
    }
    if(bestSS1Ang<0.05&&bestpTSV2>6){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(bestSS1AngMK,1);
    }
    if(bestSS2Ang<0.05&&bestpTSV3>6){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(bestSS2AngMK,1);
    }
    
    
    if(1==1){
      Mu1TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(Mu1bestChi2DistRatio,1);
    }
    if(1==1){
      Mu2TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(Mu2bestChi2DistRatio,1);
    }
    if(1==1){
      Mu3TrackInvariantMassBeforeMVASVAngle0.at(t).Fill(Mu3bestChi2DistRatio,1);
    }
    
    /*
    
    if(bestpTSV1>3&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(M1bestSV,1);
    }
    if(bestpTSV2>3&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(M2bestSV,1);
    }
    if(bestpTSV3>3&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVASVAngle1.at(t).Fill(M3bestSV,1);
    }
    */
    /*
    if(bestpTSV1>3&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(M1bestSV,1);
    }
    if(bestpTSV2>3&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(M2bestSV,1);
    }
    if(bestpTSV3>3&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVASVAngle2.at(t).Fill(M3bestSV,1);
    }
    
    
    if(bestpTSV1>10&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(M1bestSV,1);
    }
    if(bestpTSV2>10&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(M2bestSV,1);
    }
    if(bestpTSV3>10&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVASVAngle3.at(t).Fill(M3bestSV,1);
    }
    */
    
    
    //KStar
    
    if(1==1){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(Mu1KSbestChi2DistRatio,1);
    }
    if(1==1){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(Mu2KSbestChi2DistRatio,1);
    }
    if(1==1){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle0.at(t).Fill(Mu3KSbestChi2DistRatio,1);
    }
    
    /*
    if(bestpTSV1>3&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(M1bestSVKStar,1);
    }
    if(bestpTSV2>3&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(M2bestSVKStar,1);
    }
    if(bestpTSV3>3&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle1.at(t).Fill(M3bestSVKStar,1);
    }
    */
    /*
    if(bestpTSV1>3&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(M1bestSVKStar,1);
    }
    if(bestpTSV2>3&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(M2bestSVKStar,1);
    }
    if(bestpTSV3>3&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle2.at(t).Fill(M3bestSVKStar,1);
    }
    
    
    if(bestpTSV1>10&&OSdRSV<0.00001){
      Mu1TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(M1bestSVKStar,1);
    }
    if(bestpTSV2>10&&SS1dRSV<0.00001){
      Mu2TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(M2bestSVKStar,1);
    }
    if(bestpTSV3>10&&SS2dRSV<0.00001){
      Mu3TrackInvariantMassBeforeMVAKStarSVAngle3.at(t).Fill(M3bestSVKStar,1);
    }
    */
    
    
    
    /*
    if(OSdRSV<0.002){
      Mu1TrackInvariantMassBeforeMVABestdR.at(t).Fill(M1bestSV,1);
    }
    if(SS1dRSV<0.002){
      Mu2TrackInvariantMassBeforeMVABestdR.at(t).Fill(M2bestSV,1);
    }
    if(SS2dRSV<0.002){
      Mu3TrackInvariantMassBeforeMVABestdR.at(t).Fill(M3bestSV,1);
    }
    */
    
    Mu1TrackInvariantMassBeforeMVAKStarBestdR.at(t).Fill(KM1dR,1);
    Mu2TrackInvariantMassBeforeMVAKStarBestdR.at(t).Fill(KM2dR,1);
    Mu3TrackInvariantMassBeforeMVAKStarBestdR.at(t).Fill(KM3dR,1);
    
    /*
    Mu1TrackInvariantMassBeforeMVABestdRIncrease.at(t).Fill(OSdRSV,1);
    
    
    Mu1TrackInvariantMassBeforeMVABestdRDecrease.at(t).Fill(SS2dRSV,1);
    
    */
    
    //Mu1TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(dR1,1);
    //Mu2TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(dR2,1);
    //Mu3TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(dR3,1);
    
    /*
    
    if(dR1<0.125){
      Mu1TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(M1dR,1);
    }
    if(dR2<0.125){
      Mu2TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(M2dR,1);
    }
    if(dR3<0.125){
      Mu3TrackInvariantMassBeforeMVABestdRLtd.at(t).Fill(M3dR,1);
    }
    
    */
    
    Mu1TrackInvariantMassBeforeMVABestdRLtd1.at(t).Fill(M1bestChi2,1);
    
    
    Mu2TrackInvariantMassBeforeMVABestdRLtd1.at(t).Fill(M2bestChi2,1);
    
    
    Mu3TrackInvariantMassBeforeMVABestdRLtd1.at(t).Fill(M3bestChi2,1);
    
    /*
    if(Mu1bestChi2<100){
      Mu1TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<100){
      Mu2TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<100){
      Mu3TrackInvariantMassBeforeMVABestdRLtd2.at(t).Fill(M3bestChi2,1);
    }
    */
    
    if(Mu1bestChi2<20){
      Mu1TrackInvariantMassBeforeMVABestdRLtd3.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<20){
      Mu2TrackInvariantMassBeforeMVABestdRLtd3.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<20){
      Mu3TrackInvariantMassBeforeMVABestdRLtd3.at(t).Fill(M3bestChi2,1);
    }
    
    if(Mu1bestChi2<10){
      Mu1TrackInvariantMassBeforeMVABestdRLtd4.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<10){
      Mu2TrackInvariantMassBeforeMVABestdRLtd4.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<10){
      Mu3TrackInvariantMassBeforeMVABestdRLtd4.at(t).Fill(M3bestChi2,1);
    }
    
    if(Mu1bestChi2<3){
      Mu1TrackInvariantMassBeforeMVABestdRLtd5.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<3){
      Mu2TrackInvariantMassBeforeMVABestdRLtd5.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<3){
      Mu3TrackInvariantMassBeforeMVABestdRLtd5.at(t).Fill(M3bestChi2,1);
    }
    
    if(Mu1bestChi2<1){
      Mu1TrackInvariantMassBeforeMVABestdRLtd6.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<1){
      Mu2TrackInvariantMassBeforeMVABestdRLtd6.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<1){
      Mu3TrackInvariantMassBeforeMVABestdRLtd6.at(t).Fill(M3bestChi2,1);
    }
    
    if(Mu1bestChi2<0.5){
      Mu1TrackInvariantMassBeforeMVABestdRLtd7.at(t).Fill(M1bestChi2,1);
    }
    if(Mu2bestChi2<0.5){
      Mu2TrackInvariantMassBeforeMVABestdRLtd7.at(t).Fill(M2bestChi2,1);
    }
    if(Mu3bestChi2<0.5){
      Mu3TrackInvariantMassBeforeMVABestdRLtd7.at(t).Fill(M3bestChi2,1);
    }
    
    Mu1TrackInvariantMassBeforeMVABestMass.at(t).Fill(M1best,1);
    Mu2TrackInvariantMassBeforeMVABestMass.at(t).Fill(M2best,1);
    Mu3TrackInvariantMassBeforeMVABestMass.at(t).Fill(M3best,1);
    
    
    if(KdR1<0.125){//for kstars
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd.at(t).Fill(KM1dR,1);
    }
    if(KdR2<0.125){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd.at(t).Fill(KM2dR,1);
    }
    if(KdR3<0.125){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd.at(t).Fill(KM3dR,1);
    }
    
    
    
    Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd1.at(t).Fill(M1KStarbestChi2,1);
    
    
    Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd1.at(t).Fill(M2KStarbestChi2,1);
    
    
    Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd1.at(t).Fill(M3KStarbestChi2,1);
    
    
    if(Mu1bestChi2<100){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd2.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<100){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd2.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<100){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd2.at(t).Fill(M3KStarbestChi2,1);
    }
    
    if(Mu1bestChi2<20){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd3.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<20){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd3.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<20){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd3.at(t).Fill(M3KStarbestChi2,1);
    }
    
    if(Mu1bestChi2<10){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd4.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<10){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd4.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<10){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd4.at(t).Fill(M3KStarbestChi2,1);
    }
    
    if(Mu1bestChi2<3){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd5.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<3){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd5.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<3){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd5.at(t).Fill(M3KStarbestChi2,1);
    }
    
    if(Mu1bestChi2<1){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd6.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<1){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd6.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<1){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd6.at(t).Fill(M3KStarbestChi2,1);
    }
    
    if(Mu1bestChi2<0.5){
      Mu1TrackInvariantMassBeforeMVAKStarBestdRLtd7.at(t).Fill(M1KStarbestChi2,1);
    }
    if(Mu2bestChi2<0.5){
      Mu2TrackInvariantMassBeforeMVAKStarBestdRLtd7.at(t).Fill(M2KStarbestChi2,1);
    }
    if(Mu3bestChi2<0.5){
      Mu3TrackInvariantMassBeforeMVAKStarBestdRLtd7.at(t).Fill(M3KStarbestChi2,1);
    }
    
    Mu1TrackInvariantMassBeforeMVAKStarBestMass.at(t).Fill(KM1best,1);
    Mu2TrackInvariantMassBeforeMVAKStarBestMass.at(t).Fill(KM2best,1);
    Mu3TrackInvariantMassBeforeMVAKStarBestMass.at(t).Fill(KM3best,1);
    
    

    PairMass.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1); 

  
    bool RemoveEta(false);
    bool RemoveHalfEta(false);


    bool phiVeto(false);
    bool rmgVeto(false);


    double    m12v = (MuonOS+MuonSS1).M();
    double    m13v = (MuonOS+MuonSS2).M();


    if(( m12v < phiVetoCut1  || m12v > phiVetoCut2 )  && (m13v < phiVetoCut1 || m13v > phiVetoCut2)  )  phiVeto=true;
    if(( m12v < rmgCutVeto1 || m12v > rmgCutVeto2 )  && (m13v < rmgCutVeto1 || m13v > rmgCutVeto2))  rmgVeto=true;


    if((MuonOS+MuonSS1).M() > 0.549 && (MuonOS+MuonSS2).M() > 0.549) RemoveEta = true;
    if((MuonOS+MuonSS2).M() > 0.549) RemoveHalfEta = true;
    if(RemoveEta && phiVeto && rmgVeto)    PairMassWithCut.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);

    

    //*** 3-mu mass after the KF vertex constrain

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);
      
    


    //*** uncomment if you want to have print outs
    /*
    if(id ==120 ){// or id == 40){
      std::cout<<"-------------- All categoris ----------------"<< std::endl;
      std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Ntp->printMCDecayChainOfEvent(true, true, true, true);
      std::cout<< "\n\n\n\n\n\n";
    }
    */
    
    

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

    EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());  // Event Mass resolution

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


    //***  define the mva varables used for evaluation of BDT weights for selection
    
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

    var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});

    var_MaxMuon_chi2LocalPosition = std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });
    var_MaxtrkKink = std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});

    var_MaxD0SigSV=    std::max({Ntp->Vertex_d0sigSV_reco(signal_idx,0),
	  Ntp->Vertex_d0sigSV_reco(signal_idx,1),
	  Ntp->Vertex_d0sigSV_reco(signal_idx,2)});

    var_MindcaTrackSV=    Ntp->Isolation_MinDist(signal_idx);

    var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});



    var_MaxVertexPairQuality =   std::max({Ntp->Vertex_pair_quality(signal_idx,0),Ntp->Vertex_pair_quality(signal_idx,1),Ntp->Vertex_pair_quality(signal_idx,2)});
    var_MaxMuon_chi2LocalMomentum = std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });


    var_MuonglbkinkSum    = (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));



    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(signal_idx,0),
	  Ntp->Vertex_d0sig_reco(signal_idx,1),
	  Ntp->Vertex_d0sig_reco(signal_idx,2)});

    var_MaxD0Significance = MaxD0Significance;
    var_IsolationMinDist = Ntp->Isolation_MinDist(signal_idx);


    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      var_mass12_dRsorting = (MuonOS+MuonSS2).M();
      var_mass13_drSorting = (MuonOS+MuonSS1).M();
    }else{
      var_mass12_dRsorting = (MuonOS+MuonSS1).M();
      var_mass13_drSorting = (MuonOS+MuonSS2).M();
    }


    var_tauMass=TauRefitLV.M();
    TauMass_all.at(t).Fill(TauRefitLV.M(),1);


    //*** define variables for Mu ID and evaluate the BDT

    for(unsigned int imu=0; imu<3;imu++){

      mu_combinedQuality_chi2LocalMomentum=Ntp->Muon_combinedQuality_chi2LocalMomentum(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_chi2LocalPosition=Ntp->Muon_combinedQuality_chi2LocalPosition(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_staRelChi2=Ntp->Muon_combinedQuality_staRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_trkRelChi2=Ntp->Muon_combinedQuality_trkRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_globalDeltaEtaPhi=Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_trkKink=log(Ntp->Muon_combinedQuality_glbKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)));
      mu_combinedQuality_glbKink=log(Ntp->Muon_combinedQuality_trkKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)));
      mu_combinedQuality_glbTrackProbability=Ntp->Muon_combinedQuality_glbTrackProbability(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_Numberofvalidtrackerhits=Ntp->Muon_numberofValidPixelHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_Numberofvalidpixelhits=Ntp->Muon_innerTrack_numberOfValidTrackerHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_validMuonHitComb=Ntp->Muon_hitPattern_numberOfValidMuonHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_numberOfMatchedStations=Ntp->Muon_numberOfMatchedStations(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_segmentCompatibility=Ntp->Muon_segmentCompatibility(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_timeAtIpInOutErr=Ntp->Muon_timeAtIpInOutErr(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_GLnormChi2=Ntp->Muon_normChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      
      mu_innerTrack_normalizedChi2=Ntp->Muon_innerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_outerTrack_normalizedChi2= Ntp->Muon_outerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_innerTrack_validFraction=Ntp->Muon_innerTrack_validFraction(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));

      if(fabs(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)).Eta()) < 1.2    )
	{
	  if(imu==0)
	    {
	      Muon1DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }
	  if(imu==1)
	    {
	      Muon2DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }

	  if(imu==2)
	    {
	      Muon3DetID = readerMuIDBarrel->EvaluateMVA("BDT");
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
	  if(imu==2)
	    {
	      Muon3DetID= readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	}
    }

    Muon1MVAID.at(t).Fill(Muon1DetID);
    Muon2MVAID.at(t).Fill(Muon2DetID);
    Muon3MVAID.at(t).Fill(Muon3DetID);

    var_Muon1DetID = Muon1DetID; 
    var_Muon2DetID = Muon2DetID;
    var_Muon3DetID = Muon3DetID;



    double dRSortedMass;

    //*** define per event resolution categroies 
    //*** Category A1
    //    if(phiVeto && rhoVeto)
      {
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {

	    TauMass_allVsBDTA.at(t).Fill(TauRefitLV.M(),readerA->EvaluateMVA("BDT"));
	    BDTOutputA.at(t).Fill(    readerA->EvaluateMVA("BDT"),1 );

	    if(readerA->EvaluateMVA("BDT") > mvaA2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
		

		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassA1.at(t).Fill(TauLV.M(),1);                 // three mu mass 
		  TauMassRefitA1.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
		  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		  if(RemoveEta)	TauMassRefitA1MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta) TauMassRefitA1HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitA1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);    //  this is to be checked
        
        //std::cout<<"is it working?"<<std::endl;
        
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterA1MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterA1MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterA1MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    }
        
        
        
        
        //stuff
		
		}
	    }
	  }

	//Category B1
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {

	    TauMass_allVsBDTB.at(t).Fill(TauRefitLV.M(),readerB->EvaluateMVA("BDT"));
	    BDTOutputB.at(t).Fill(readerB->EvaluateMVA("BDT"), 1);

	    if(readerB->EvaluateMVA("BDT") > mvaB2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	  
		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassB1.at(t).Fill(TauLV.M(),1);                  // three mu mass 
		  TauMassRefitB1.at(t).Fill(TauRefitLV.M(),1);        // three mu KF reffited mass
		  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);      // fill up all categories inclusive
		  if(RemoveEta)	TauMassRefitB1MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta)	TauMassRefitB1HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitB1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);    //  this is to be checked
        
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterB1MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterB1MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterB1MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    }
    
    
		
		}
	    }
	  }

	//Category C1
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {

	    TauMass_allVsBDTC.at(t).Fill(TauRefitLV.M(),readerC->EvaluateMVA("BDT"));
	    BDTOutputC.at(t).Fill(    readerC->EvaluateMVA("BDT") );

	    if(readerC->EvaluateMVA("BDT") > mvaC2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);

		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassC1.at(t).Fill(TauLV.M(),1);	          // three mu mass 
		  TauMassRefitC1.at(t).Fill(TauRefitLV.M(),1);      // three mu KF reffited mass
		  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);    // fill up all categories inclusive
		  if(RemoveEta)	TauMassRefitC1MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta)	TauMassRefitC1HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitC1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1); 
      
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterC1MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterC1MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterC1MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    }
    
       
		
		}
	    }
	  }

	//Category A2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {

	    if(readerA->EvaluateMVA("BDT") > mvaA1_ && readerA->EvaluateMVA("BDT") < mvaA2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	  
		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassA2.at(t).Fill(TauLV.M(),1);
		  TauMassRefitA2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveEta)	TauMassRefitA2MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta)	TauMassRefitA2HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitA2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);    
        
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterA2MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterA2MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterA2MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    }
    
    

		}
	    }
	  }

	//Category B2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(readerB->EvaluateMVA("BDT") > mvaB1_ && readerB->EvaluateMVA("BDT") < mvaB2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	  
		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassB2.at(t).Fill(TauLV.M(),1);
		  TauMassRefitB2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveEta)	TauMassRefitB2MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta)	TauMassRefitB2HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitB2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);  
      
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterB2MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterB2MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterB2MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    } 
    
     

		}
	    }
	  }
    
	//Category C2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {
	    if(readerC->EvaluateMVA("BDT") > mvaC1_ && readerC->EvaluateMVA("BDT")< mvaC2_){
	      if(phiVeto && rmgVeto)
		{
		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	       
		  //*** defined the pair with SS best alligned to OS
		  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		    dRSortedMass = (MuonOS+MuonSS2).M();
		  }else{
		    dRSortedMass = (MuonOS+MuonSS1).M();
		  }
		  //***

		  TauMassC2.at(t).Fill(TauLV.M(),1);	      
		  TauMassRefitC2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    

		  if(RemoveEta)	TauMassRefitC2MassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(RemoveHalfEta)	TauMassRefitC2HalfMassCut.at(t).Fill(TauRefitLV.M(),1);    
		  if(dRSortedMass < 0.549) TauMassRefitC2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);    
        
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){
      TLorentzVector ParticleLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TLorentzVector ParticleLVReassigned = ParticleLV;//reassign the masses to be similar to a kaon
      ParticleLVReassigned.SetE(sqrt(ParticleLV.Px()*ParticleLV.Px()+ParticleLV.Py()*ParticleLV.Py()+ParticleLV.Pz()*ParticleLV.Pz()+0.493677*0.493677));
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(os_mu_idx))){
        Mu1TrackInvariantMassAfterC2MVA.at(t).Fill((ParticleLVReassigned+MuonOSReassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss1_mu_idx))){
        Mu2TrackInvariantMassAfterC2MVA.at(t).Fill((ParticleLVReassigned+MuonSS1Reassigned).M(),1);
      }
      
      if(Ntp->IsolationTrack_charge(signal_idx,j)==(-1*Ntp->Muon_charge(ss2_mu_idx))){
        Mu3TrackInvariantMassAfterC2MVA.at(t).Fill((ParticleLVReassigned+MuonSS2Reassigned).M(),1);
      }
      
      
    }
    
    
        
		}

	    }
	  }


	//*** below are basic  purity and resolution plots
	if(id==40 || id == 60 || id ==90){
	  if(Ntp->MCEventIsReconstructed()){
	    TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	    TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	    TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	    TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
	
	    TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	    TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

	    Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	    Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	    Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
	  }
	}
    
    
	//*** fill up the T3MMiniTree.root for statistical analysis
    
	m3m = TauRefitLV.M();
    
	dataMCtype = id;
	event_weight =1; // 1 for data
	if(dataMCtype == 1){event_weight =1;}
	else if(dataMCtype == 40){event_weight =0.00128;} // event_weight is a value Lumi Scale 
	else if(dataMCtype == 60){event_weight =0.000497;}
	else if(dataMCtype == 90){event_weight =0.00149;}
    
    
	mvaA1 = mvaA1_;
	mvaA2 = mvaA2_;
	mvaB1 = mvaB1_;
	mvaB2 = mvaB2_;
	mvaC1 = mvaC1_;
	mvaC2 = mvaC2_;
    

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_ ){
	  category = 1;
	  bdt = readerA->EvaluateMVA("BDT");
	}
    
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_){
	  category = 2;
	  bdt = readerB->EvaluateMVA("BDT");
	}

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_){
	  category = 3;
	  bdt = readerC->EvaluateMVA("BDT");
	}

	m12 = (MuonOS+MuonSS1).M();
	m13 = (MuonOS+MuonSS2).M();
	LumiScale = 1.;
	T3MMiniTree->Fill();
      }
      
   
  }// end of if status
}


void  CommonSelector::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis
  T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();


  //*** extra actions
  if(mode == RECONSTRUCT){

    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    //    int id(Ntp->GetMCID());
    //    double scale(1.);
    //    double scaleDsTau(0.637);
    //    double scaleBpTau(0.262);
    //    double scaleB0Tau(0.099);
    //total xsection of producing taus is 12.848 ub 
    // if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    // ScaleAllHistOfType(2,scale*scaleDsTau);
    // if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    // ScaleAllHistOfType(3,scale*scaleB0Tau);
    // if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    // ScaleAllHistOfType(4,scale*scaleBpTau);
    //    }
  }
  Selection::Finish();
}





