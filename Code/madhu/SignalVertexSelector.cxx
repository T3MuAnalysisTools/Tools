#include "SignalVertexSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TH3.h"


using namespace std;

SignalVertexSelector::SignalVertexSelector(TString Name_, TString id_):
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
  PEMassResolutionCut2_(0.0105),
  mvaA1_(0.132), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2_(0.233),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1_(0.148),
  mvaB2_(0.279),
  mvaC1_(0.170),
  mvaC2_(0.257),
  mvaBTrainA1_(0.090),
  mvaBTrainA2_(0.182),
  mvaBTrainB1_(0.121),
  mvaBTrainB2_(0.211),
  mvaBTrainC1_(0.127),
  mvaBTrainC2_(0.211),
  mvaDTrainA1_(0.153),
  mvaDTrainA2_(0.212),
  mvaDTrainB1_(0.154),
  mvaDTrainB2_(0.224),
  mvaDTrainC1_(0.129),
  mvaDTrainC2_(0.210)
{
  // This is a class constructor;
}



SignalVertexSelector::~SignalVertexSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SignalVertexSelector::Configure(){

  gErrorIgnoreLevel = kFatal;
  
  //  This mini tree is for limit extraction

  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("xv",&xv);
  T3MMiniTree->Branch("phiv",&phiv);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("bdt",&bdt);
  T3MMiniTree->Branch("category",&category);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  T3MMiniTree->Branch("mDr1",&mDr1);
  T3MMiniTree->Branch("mDr2",&mDr2);
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
  readerA->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerA->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerA->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerA->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerA->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerA->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerA->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerA->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerA->AddVariable("var_NtracksClose", &var_NtracksClose);
  readerA->AddVariable("var_Iso08", &var_Iso08);
  

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_0_A/weights/TMVAClassification_BDT.weights.xml"); 




  readerB = new TMVA::Reader( "!Color:!Silent" );
  readerB->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerB->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerB->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerB->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerB->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerB->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerB->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerB->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerB->AddVariable( "var_NtracksClose", &var_NtracksClose);
  readerB->AddVariable("var_Iso08", &var_Iso08);
  

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerB->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_0_B/weights/TMVAClassification_BDT.weights.xml"); 


  readerC = new TMVA::Reader( "!Color:!Silent" );
  readerC->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerC->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerC->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerC->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerC->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerC->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerC->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerC->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerC->AddVariable( "var_NtracksClose", &var_NtracksClose);
  readerC->AddVariable("var_Iso08", &var_Iso08);
  

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerC->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_0_C/weights/TMVAClassification_BDT.weights.xml"); 





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


  readerBvsD= new TMVA::Reader( "!Color:!Silent" );
  readerBvsD->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBvsD->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBvsD->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBvsD->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBvsD->AddVariable("var_nsv",&var_nsv);
  readerBvsD->AddVariable("var_MaxD0SigBS",&var_MaxD0SigBS);
  readerBvsD->AddVariable("var_MinD0SigBS",&var_MinD0SigBS);
  readerBvsD->AddVariable("var_Iso08",&var_Iso08);
  readerBvsD->AddVariable("var_dcaTrackPV",&var_dcaTrackPV);
  readerBvsD->AddVariable("var_MinMuonImpactAngle",&var_MinMuonImpactAngle);
  readerBvsD->AddVariable("var_flightLenDist",&var_flightLenDist);
  readerBvsD->BookMVA( "BDTG", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_0_MCTrainA/weights/TMVAClassification_BDTG.weights.xml" );



  readerBTrainA= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainA->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainA->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainA->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainA->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainA->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainA->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainA->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainA->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainA->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainA->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainA->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_A_B/weights/TMVAClassification_BDT.weights.xml");


  readerBTrainB= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainB->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainB->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainB->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainB->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainB->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainB->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainB->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainB->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainB->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainB->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainB->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_B_B/weights/TMVAClassification_BDT.weights.xml");



  readerBTrainC= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainC->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainC->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainC->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainC->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainC->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainC->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainC->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainC->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainC->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainC->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainC->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_C_B/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainA= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainA->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainA->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainA->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainA->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainA->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainA->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainA->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainA->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainA->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainA->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainA->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_A_DS/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainB= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainB->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainB->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainB->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainB->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainB->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainB->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainB->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainB->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainB->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainB->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainB->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_B_DS/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainC= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainC->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainC->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainC->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainC->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainC->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainC->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainC->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainC->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainC->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainC->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainC->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVATree_Verify1Apr_27_2021/Code/CommonUtils/IterativeTrain/output_7_C_DS/weights/TMVAClassification_BDT.weights.xml");



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
  
  Mu1TrackInvariantMassBeforeMVABestdRDecrease =HConfig.GetTH1D(Name+"_Mu1TrackInvariantMassBeforeMVABestdRDecrease","Mu1TrackInvariantMassBeforeMVABestdRDecrease",100,0,0.8,"#Delta #theta - #pi SS2","Events");
  
  
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
  TauMassResolutionHelixRefit=HConfig.GetTH1D(Name+"_TauMassResolutionHelixRefit","TauMassResolutionHelixRefit",50,-0.2,0.2,"Helix refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

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



  TauMassRefitABC1 =HConfig.GetTH1D(Name+"_TauMassRefitABC1","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2 =HConfig.GetTH1D(Name+"_TauMassRefitABC2","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2)","Events");


  TauMassRefitABC1_BDSeparateTrain=HConfig.GetTH1D(Name+"_TauMassRefitABC1_BDSeparateTrain","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2_BDSeparateTrain=HConfig.GetTH1D(Name+"_TauMassRefitABC2_BDSeparateTrain","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2)","Events");



  TauMassRefitABC1_eta =HConfig.GetTH2D(Name+"_TauMassRefitABC1_eta","Refit #tau lepton mass vs eta",30,1.5,2.1,30,0,2.5,"M_{3#mu} , GeV (inclusive ABC1)","#eta_{#tau}");
  TauMassRefitABC2_eta =HConfig.GetTH2D(Name+"_TauMassRefitABC2_eta","Refit #tau lepton mass vs eta",30,1.5,2.1,30,0,2.5,"M_{3#mu} , GeV (inclusive ABC2)","#eta_{#tau}");


  TauMassB1 =HConfig.GetTH1D(Name+"_TauMassB1","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitB1 =HConfig.GetTH1D(Name+"_TauMassRefitB1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (B1)","Events");
  TauMassRefitB1MassCut =HConfig.GetTH1D(Name+"_TauMassRefitB1MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B1)","Events");
  TauMassRefitB2MassCut =HConfig.GetTH1D(Name+"_TauMassRefitB2MassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B2)","Events");


  TauMassRefitB1HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitB1HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (B1)","Events");
  TauMassRefitB2HalfMassCut =HConfig.GetTH1D(Name+"_TauMassRefitB2HalfMassCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (half #eta veto) (B2)","Events");


  TauMassRefitB1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitB1FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B1)","Events");
  TauMassRefitB2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitB2FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (B2)","Events");


  TauMassRefitABC1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitABC1FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (ABC1)","Events");
  TauMassRefitABC2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitABC2FullEtaVetoCut","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (ABC2)","Events");





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
  VertexChi2KF_vs_HelixFit=HConfig.GetTH2D(Name+"_VertexChi2KF_vs_HelixFit","VertexChi2KF_vs_HelixFit",100,0,100,100,0,100,"Kalman Vertex #chi^{2}","Helix Vertex  Fitter #chi^{2}");

  KF_Helix_deltaX=HConfig.GetTH1D(Name+"_KF_Helix_deltaX","KF_Helix_deltaX",50,-0.05,0.05,"#Delta X, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaY=HConfig.GetTH1D(Name+"_KF_Helix_deltaY","KF_Helix_deltaY",50,-0.05,0.05,"#Delta Y, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaZ=HConfig.GetTH1D(Name+"_KF_Helix_deltaZ","KF_Helix_deltaZ",50,-0.05,0.05,"#Delta Z, cm (Helix Fitter - Kalman Fitter)","Events");

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

  BvsDBDTG  = HConfig.GetTH1D(Name+"_BvsDBDTG","BvsDBDTG",50,-1.0,1.0," B vs D BDTG","Events");

  BvsDBDTG_ABC1  = HConfig.GetTH1D(Name+"_BvsDBDTG_ABC1","BvsDBDTG_ABC1",50,-1.0,1.0," B vs D BDTG (ABC1 inclusive)","Events");
  BvsDBDTG_ABC2  = HConfig.GetTH1D(Name+"_BvsDBDTG_ABC2","BvsDBDTG_ABC2",50,-1.0,1.0," B vs D BDTG (ABC2 inclusive)","Events");

  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");
  PairMass=HConfig.GetTH2D(Name+"_PairMass","PairMass",100,0.2,1.8,100,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");


  KKMass_dR_sort=HConfig.GetTH2D(Name+"_KKMass_dR_sort","KKMass_dR_sort",110,0.9,2.0,110,0.9,2.0,"M_{1}(K^{+}K^{-}), GeV (dR sort)","M_{2}(K^{+}K^{-}), GeV (dR sort)");
  KKMass_dR_sort1=HConfig.GetTH1D(Name+"_KKMass_dR_sort1","KKMass_dR_sort1",100,0.9,1.8,"M_{1}(K^{+}K^{-}), GeV (dR sort)","");
  KKMass_dR_sort2=HConfig.GetTH1D(Name+"_KKMass_dR_sort2","KKMass_dR_sort2",100,0.9,1.8,"M_{2}(K^{+}K^{-}), GeV (dR sort)","");


  KKMass_pt_sort=HConfig.GetTH2D(Name+"_KKMass_pt_sort","KKMass_pt_sort",110,0.9,2.0,110,0.9,2.0,"M_{1}(K^{+}K^{-}), GeV (pt sort)","M_{2}(K^{+}K^{-}), GeV (pt sort)");
  KKMass_pt_sort1=HConfig.GetTH1D(Name+"_KKMass_pt_sort1","KKMass_pt_sort1",100,0.9,1.8,"M_{1}(K^{+}K^{-}), GeV (pt sort)","");
  KKMass_pt_sort2=HConfig.GetTH1D(Name+"_KKMass_pt_sort2","KKMass_pt_sort2",100,0.9,1.8,"M_{2}(K^{+}K^{-}), GeV (pt sort)","");


  KKMass_dR_sort_XVeto=HConfig.GetTH2D(Name+"_KKMass_dR_sort_XVeto","KKMass_dR_sort_XVeto",110,0.9,2.0,110,0.9,2.0,"M_{1}(K^{+}K^{-}), GeV (dR sort) XV","M_{2}(K^{+}K^{-}), GeV (dR sort) XV");
  KKMass_dR_sort1_XVeto=HConfig.GetTH1D(Name+"_KKMass_dR_sort1_XVeto","KKMass_dR_sort1_XVeto",100,0.9,1.8,"M_{1}(K^{+}K^{-}), GeV (dR sort) XV","");
  KKMass_dR_sort2_XVeto=HConfig.GetTH1D(Name+"_KKMass_dR_sort2_XVeto","KKMass_dR_sort2_XVeto",100,0.9,1.8,"M_{2}(K^{+}K^{-}), GeV (dR sort) XV","");

  KKMass_pt_sort_XVeto=HConfig.GetTH2D(Name+"_KKMass_pt_sort_XVeto","KKMass_pt_sort_XVeto",110,0.9,2.0,110,0.9,2.0,"M_{1}(K^{+}K^{-}), GeV (dR sort)","M_{2}(K^{+}K^{-}), GeV (dR sort)");
  KKMass_pt_sort1_XVeto=HConfig.GetTH1D(Name+"_KKMass_pt_sort1_XVeto","KKMass_pt_sort1_XVeto",100,0.9,1.8,"M_{1}(K^{+}K^{-}), GeV (dR sort)","");
  KKMass_pt_sort2_XVeto=HConfig.GetTH1D(Name+"_KKMass_pt_sort2_XVeto","KKMass_pt_sort2_XVeto",100,0.9,1.8,"M_{2}(K^{+}K^{-}), GeV (dR sort)","");




 
  KpiIsolationMass_OS=HConfig.GetTH1D(Name+"_KpiIsolationMass_OS","KpiIsolationMass_OS",100,0.6,1.8,"M_{1}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS1=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS1","KpiIsolationMass_SS1",100,0.6,1.8,"M_{2}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS2=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS2","KpiIsolationMass_SS2",100,0.6,1.8,"M_{3}(K#pi), GeV (comb. iso #pi)","");
 

  PairMass1NoSorting=HConfig.GetTH1D(Name+"_PairMass1NoSorting","PairMass1NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","Events");
  PairMass2NoSorting=HConfig.GetTH1D(Name+"_PairMass2NoSorting","PairMass2NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2 no sorting), GeV","Events");
  MuMuMassNoSorting=HConfig.GetTH2D(Name+"_MuMuMassNoSorting","MuMuMassNoSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","M_{#mu#mu} (OS-SS, 2 no sorting sorting), GeV");


  PairMass1PTSorting=HConfig.GetTH1D(Name+"_PairMass1PTSorting","PairMass1PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","Events");
  PairMass2PTSorting=HConfig.GetTH1D(Name+"_PairMass2PTSorting","PairMass2PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV","Events");
  MuMuMassPTSorting=HConfig.GetTH2D(Name+"_MuMuMassPTSorting","MuMuMassPTSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV");


  PairMass1AllignedSorting=HConfig.GetTH1D(Name+"_PairMass1AllignedSorting","PairMass1AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV","Events");
  PairMass2AllignedSorting=HConfig.GetTH1D(Name+"_PairMass2AllignedSorting","PairMass2AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","Events");
  MuMuMassAllignedSorting=HConfig.GetTH2D(Name+"_MuMuMassAllignedSorting","MuMuMassAllignedSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","M_{#mu#mu} (OS-SS, 1st collimated pair) GeV");



  PairMassDRSorted1A=HConfig.GetTH1D(Name+"_PairMassDRSorted1A","PairMassDRSorted1A",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV A","Events");
  PairMassDRSorted2A=HConfig.GetTH1D(Name+"_PairMassDRSorted2A","PairMassDRSorted2A",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV A","Events");


  PairMassDRSorted1B=HConfig.GetTH1D(Name+"_PairMassDRSorted1B","PairMassDRSorted1B",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV B","Events");
  PairMassDRSorted2B=HConfig.GetTH1D(Name+"_PairMassDRSorted2B","PairMassDRSorted2B",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV B","Events");

  PairMassDRSorted1C=HConfig.GetTH1D(Name+"_PairMassDRSorted1C","PairMassDRSorted1C",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV C","Events");
  PairMassDRSorted2C=HConfig.GetTH1D(Name+"_PairMassDRSorted2C","PairMassDRSorted2C",80,0.5,1.2,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV C","Events");




  PairMassdRSorted=HConfig.GetTH2D(Name+"_PairMassdRSorted","PairMassdRSorted",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassVertexSorted=HConfig.GetTH2D(Name+"_PairMassVertexSorted","PairMassVertexSorted",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1VertexSorting=HConfig.GetTH1D(Name+"_PairMass1VertexSorting","PairMass1VertexSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st vertex sorting pair), GeV","Events");
  PairMass2VertexSorting=HConfig.GetTH1D(Name+"_PairMass2VertexSorting","PairMass2VertexSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd vertex sorting pair), GeV","Events");




  PairMassPhiMassSorting=HConfig.GetTH2D(Name+"_PairMassPhiMassSorting","PairMassPhiMassSorting",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1PhiMassSorting=HConfig.GetTH1D(Name+"_PairMass1PhiMassSorting","PairMass1PhiMassSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st #phi mass sorting pair), GeV","Events");
  PairMass2PhiMassSorting=HConfig.GetTH1D(Name+"_PairMass2PhiMassSorting","PairMass2PhiMassSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd #phi mass sorting pair), GeV","Events");



  PairMass1TauPhiMassSorting=HConfig.GetTH2D(Name+"_PairMass1TauPhiMassSorting","PairMass1TauPhiMassSorting",200,0.2,1.8,100,1.6,2.0,"M_{OS}, GeV","M_{#tau}, GeV");
  PairMass2TauPhiMassSorting=HConfig.GetTH2D(Name+"_PairMass2TauPhiMassSorting","PairMass2TauPhiMassSorting",200,0.2,1.8,100,1.6,2.0,"M_{OS}, GeV","M_{#tau}, GeV");


  PairMassdRSortedXVeto=HConfig.GetTH2D(Name+"_PairMassdRSortedXVeto","PairMassdRSortedXVeto",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");



  PairMassFinalSel=HConfig.GetTH2D(Name+"_PairMassFinalSel","PairMassFinalSel",60,0.2,1.8,60,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1=HConfig.GetTH1D(Name+"_PairMass1","PairMass1",80,0.2,1.777,"M_{1}, GeV","");
  PairMass2=HConfig.GetTH1D(Name+"_PairMass2","PairMass2",80,0.2,1.777,"M_{2}, GeV","");



  AllignSortMass1=HConfig.GetTH1D(Name+"_AllignSortMass1","AllignSortMass1",80,0.2,1.777,"M_{1} (#Delta R OS sorted), GeV","");
  AllignSortMass2=HConfig.GetTH1D(Name+"_AllignSortMass2","AllignSortMass2",80,0.2,1.777,"M_{2} (#Delta R OS sorted), GeV","");





  //  AllignSortMass1XVeto=HConfig.GetTH1D(Name+"_AllignSortMass1XVeto","AllignSortMass1XVeto",80,0.2,1.777,"M_{1} (#Delta R OS sorted), GeV","");
  //  AllignSortMass2XVeto=HConfig.GetTH1D(Name+"_AllignSortMass2XVeto","AllignSortMass2XVeto",80,0.2,1.777,"M_{2} (#Delta R OS sorted), GeV","");




  PairMassWithCut=HConfig.GetTH2D(Name+"_PairMassWithCut","PairMassWithCut",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEta=HConfig.GetTH2D(Name+"_PairMassEta","PairMassEta",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEtaPrime=HConfig.GetTH2D(Name+"_PairMassEtaPrime","PairMassEtaPrime",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");

  IDOriginOfOSMuon =HConfig.GetTH1D(Name+"_IDOriginOfOSMuon","IDOriginOfOSMuon",400,200,600,"PDGID of OS muon origin","Events");

  Muon1MVAID=HConfig.GetTH1D(Name+"_Muon1MVAID","Muon1MVAID",50,0.0,1.0,"#mu_{1} MVA","Events");
  Muon2MVAID=HConfig.GetTH1D(Name+"_Muon2MVAID","Muon2MVAID",50,0.0,1.0,"#mu_{2} MVA","Events");
  Muon3MVAID=HConfig.GetTH1D(Name+"_Muon3MVAID","Muon3MVAID",50,0.0,1.0,"#mu_{3} MVA","Events");



  BetterMuMuVertex=HConfig.GetTH1D(Name+"_BetterMuMuVertex","BetterMuMuVertex",30,0,5,"vertex pair quality (close)","");
  WorseMuMuVertex=HConfig.GetTH1D(Name+"_WorseMuMuVertex","WorseMuMuVertex",30,0,5,"vertex pair quality (far)","");
  
  WhetherTau3Mu=HConfig.GetTH1D(Name+"_WhetherTau3Mu","WhetherTau3Mu",2,-0.5,1.5,"whether tau with 3 mu was found","Entries");
  WhetherdRMatch=HConfig.GetTH1D(Name+"_WhetherdRMatch","WhetherdRMatch",2,-0.5,1.5,"whether dR match is found","Entries");
  
  IsoTrackToMCdR01 =HConfig.GetTH1D(Name+"_IsoTrackToMCdR01","IsoTrackToMCdR01",100,0,0.1,"Reconstructed Track to MC dR","Entries");
  IsoTrackToMCdR08 =HConfig.GetTH1D(Name+"_IsoTrackToMCdR08","IsoTrackToMCdR08",100,0,0.8,"Reconstructed Track to MC dR","Entries");
  IsoTrackToMCAngle01 =HConfig.GetTH1D(Name+"_IsoTrackToMCAngle01","IsoTrackToMCAngle01",50,0,0.1,"Reconstructed Track to MC Angle","Entries");
  
  TrackToTauDr2Prong =HConfig.GetTH1D(Name+"_TrackToTauDr2Prong","TrackToTauDr2Prong",100,0,1.0,"Reconstructed Track to Tau dR 2 Prong","Entries");
  TrackToTauDr3Prong =HConfig.GetTH1D(Name+"_TrackToTauDr3Prong","TrackToTauDr3Prong",100,0,1.0,"Reconstructed Track to Tau dR 3 Prong","Entries");
  TrackToTauDrAll =HConfig.GetTH1D(Name+"_TrackToTauDrAll","TrackToTauDrAll",100,0,1.5,"All Reconstructed Track to Tau dR","Entries");
  
  NumberOfFS_ChargedParticles =HConfig.GetTH1D(Name+"_NumberOfFS_ChargedParticles","NumberOfFS_ChargedParticles",8,-0.5,7.5,"No Of final state particles","Entries");
  NumberOfFS_ChargedParticles_RecoMatch =HConfig.GetTH1D(Name+"_NumberOfFS_ChargedParticles_RecoMatch","NumberOfFS_ChargedParticles_RecoMatch",8,-0.5,7.5,"No Of final state particles reconstructed","Entries");
  
  NumberOfRecoChargedParticlesIfMC1 =HConfig.GetTH1D(Name+"_NumberOfRecoChargedParticlesIfMC1","NumberOfRecoChargedParticlesIfMC1",8,-0.5,7.5,"No Of final state particles reconstructed","Entries");
  NumberOfRecoChargedParticlesIfMC2 =HConfig.GetTH1D(Name+"_NumberOfRecoChargedParticlesIfMC2","NumberOfRecoChargedParticlesIfMC2",8,-0.5,7.5,"No Of final state particles reconstructed","Entries");
  NumberOfRecoChargedParticlesIfMC3 =HConfig.GetTH1D(Name+"_NumberOfRecoChargedParticlesIfMC3","NumberOfRecoChargedParticlesIfMC3",8,-0.5,7.5,"No Of final state particles reconstructed","Entries");
  
  TwoProngInvariantMassReco =HConfig.GetTH1D(Name+"_TwoProngInvariantMassReco","TwoProngInvariantMassReco",50,0,5.0,"Invariant Mass of two tracks, reco","Events");
  TwoProngInvariantMassMC =HConfig.GetTH1D(Name+"_TwoProngInvariantMassMC","TwoProngInvariantMassMC",50,0,5.0,"Invariant Mass of two tracks, MC","Events");
  
  ThreeProngInvariantMassReco =HConfig.GetTH1D(Name+"_ThreeProngInvariantMassReco","ThreeProngInvariantMassReco",50,0,5.0,"Invariant Mass of three tracks, reco","Events");
  ThreeProngInvariantMassMC =HConfig.GetTH1D(Name+"_ThreeProngInvariantMassMC","ThreeProngInvariantMassMC",50,0,5.0,"Invariant Mass of three tracks, MC","Events");
  
  TwoProngInvariantMassReco005 =HConfig.GetTH1D(Name+"_TwoProngInvariantMassReco005","TwoProngInvariantMassReco005",50,0,5.0,"Invariant Mass of two tracks, reco","Events");
  TwoProngInvariantMassMC005 =HConfig.GetTH1D(Name+"_TwoProngInvariantMassMC005","TwoProngInvariantMassMC005",50,0,5.0,"Invariant Mass of two tracks, MC","Events");
  
  ThreeProngInvariantMassReco005 =HConfig.GetTH1D(Name+"_ThreeProngInvariantMassReco005","ThreeProngInvariantMassReco005",50,0,5.0,"Invariant Mass of three tracks, reco","Events");
  ThreeProngInvariantMassMC005 =HConfig.GetTH1D(Name+"_ThreeProngInvariantMassMC005","ThreeProngInvariantMassMC005",50,0,5.0,"Invariant Mass of three tracks, MC","Events");
  
  TwoProngTrackPt =HConfig.GetTH1D(Name+"_TwoProngTrackPt","TwoProngTrackPt",50,0,5.0,"Track Pt","Events");
  TwoProngTrack2Pt =HConfig.GetTH1D(Name+"_TwoProngTrack2Pt","TwoProngTrack2Pt",50,0,5.0,"Track Pt","Events");
  
  ThreeProngTrackPt =HConfig.GetTH1D(Name+"_ThreeProngTrackPt","ThreeProngTrackPt",50,0,5.0,"Track Pt","Events");
  ThreeProngTrack2Pt =HConfig.GetTH1D(Name+"_ThreeProngTrack2Pt","ThreeProngTrack2Pt",50,0,5.0,"Track Pt","Events");
  ThreeProngTrack3Pt =HConfig.GetTH1D(Name+"_ThreeProngTrack3Pt","ThreeProngTrack3Pt",50,0,5.0,"Track Pt","Events");
  
  TwoProngTrackEta =HConfig.GetTH1D(Name+"_TwoProngTrackEta","TwoProngTrackEta",50,0,5.0,"Track Eta","Events");
  TwoProngTrack2Eta =HConfig.GetTH1D(Name+"_TwoProngTrack2Eta","TwoProngTrack2Eta",50,0,5.0,"Track Eta","Events");
  
  dR_vs_dP=HConfig.GetTH2D(Name+"_dR_vs_dP","dR_vs_dP",40,0,0.1,40,0,2.0,"dR","dP");
  
  InvMass2_vs_pdgid=HConfig.GetTH2D(Name+"_InvMass2_vs_pdgid","InvMass2_vs_pdgid",50,0,5.0,1000,-0.5,999.5,"Invariant Mass MC 2 Prong","pdgid");
  InvMass3_vs_pdgid=HConfig.GetTH2D(Name+"_InvMass3_vs_pdgid","InvMass3_vs_pdgid",50,0,5.0,1000,-0.5,999.5,"Invariant Mass MC 3 Prong","pdgid");
  
  NoOfIsoTracks2Prong =HConfig.GetTH1D(Name+"_NoOfIsoTracks2Prong","NoOfIsoTracks2Prong",40,-0.5,39.5,"No of isolation tracks for 2 prong cat","Events");
  NoOfIsoTracks3Prong =HConfig.GetTH1D(Name+"_NoOfIsoTracks3Prong","NoOfIsoTracks3Prong",40,-0.5,39.5,"No of isolation tracks for 3 prong cat","Events");
  
  dRmin_sum_vs_InvariantMass_2prong=HConfig.GetTH2D(Name+"_dRmin_sum_vs_InvariantMass_2prong","dRmin_sum_vs_InvariantMass_2prong",40,0,0.03,40,0,4.0,"dR","Mass");
  dRmin_sum_vs_InvariantMass_3prong=HConfig.GetTH2D(Name+"_dRmin_sum_vs_InvariantMass_3prong","dRmin_sum_vs_InvariantMass_3prong",40,0,0.03,40,0,4.0,"dR","Mass");
  
  RankMatchedTrackpT =HConfig.GetTH1D(Name+"_RankMatchedTrackpT","RankMatchedTrackpT",40,-0.5,8.5,"Rank of matched iso track in pT sorted tracks after dR sorting","Events");
  RankMatchedTrackdR =HConfig.GetTH1D(Name+"_RankMatchedTrackdR","RankMatchedTrackdR",40,-0.5,39.5,"Rank of matched iso trk among dR sorted tracks","Events");
  RankMatchedTrackdR_trim =HConfig.GetTH1D(Name+"_RankMatchedTrackdR_trim","RankMatchedTrackdR_trim",7,-0.5,6.5,"Rank of matched iso trk among dR sorted tracks","Events");
  
  RankMatchedTrackAvgDiff =HConfig.GetTH1D(Name+"_RankMatchedTrackAvgDiff","RankMatchedTrackAvgDiff",40,-0.5,8.5,"Rank of matched iso trk in tracks sorted using vtx with muon","Events");
  RankMatchedTrackCombn =HConfig.GetTH1D(Name+"_RankMatchedTrackCombn","RankMatchedTrackCombn",40,-0.5,39.5,"Rank of matched iso trk among vtx position sorted tracks","Events");
  
  var_All_7_Iso_dR =HConfig.GetTH1D(Name+"_var_All_7_Iso_dR","var_All_7_Iso_dR",100,0,1.0,"Isolation Track to Tau dR","Entries");
  var_Correct_Iso_dR =HConfig.GetTH1D(Name+"_var_Correct_Iso_dR","var_Correct_Iso_dR",100,0,1.0,"Matched Track to Tau dR","Entries");
  
  var_All_7_Iso_AvgDiff =HConfig.GetTH1D(Name+"_var_All_7_Iso_AvgDiff","var_All_7_Iso_AvgDiff",100,0,2.5,"Isolation Track, Average Vertex Distance","Entries");
  var_Correct_Iso_AvgDiff =HConfig.GetTH1D(Name+"_var_Correct_Iso_AvgDiff","var_Correct_Iso_AvgDiff",100,0,2.5,"Matched Track, Average Vertex Distance","Entries");
  
  var_All_7_Iso_Avg =HConfig.GetTH1D(Name+"_var_All_7_Iso_Avg","var_All_7_Iso_Avg",100,0,5.0,"Isolation Track, Average Vertex Distance to SV","Entries");
  var_Correct_Iso_Avg =HConfig.GetTH1D(Name+"_var_Correct_Iso_Avg","var_Correct_Iso_Avg",100,0,5.0,"Matched Track, Average Vertex Distance to SV","Entries");
  
  var_All_7_Iso_Combn =HConfig.GetTH1D(Name+"_var_All_7_Iso_Combn","var_All_7_Iso_Combn",200,0,5.0,"Isolation Track, Avg Distance x dR x sqrt(chi2)","Entries");
  var_Correct_Iso_Combn =HConfig.GetTH1D(Name+"_var_Correct_Iso_Combn","var_Correct_Iso_Combn",200,0,5.0,"Matched Track, Avg Distance x dR x sqrt(chi2)","Entries");
  
  var_All_7_Iso_Trial =HConfig.GetTH1D(Name+"_var_All_7_Iso_Trial","var_All_7_Iso_Trial",100,0,10.0,"Isolation Track, Trial","Entries");
  var_Correct_Iso_Trial =HConfig.GetTH1D(Name+"_var_Correct_Iso_Trial","var_Correct_Iso_Trial",100,0,10.0,"Matched Track, Trial","Entries");
  
  RankMatchedTrackdR_cut =HConfig.GetTH1D(Name+"_RankMatchedTrackdR_cut","RankMatchedTrackdR_cut",7,-0.5,6.5,"Rank of matched iso trk in dR sorted tracks","Events");
  
  RankMatchedTrackPairdR =HConfig.GetTH1D(Name+"_RankMatchedTrackPairdR","RankMatchedTrackPairdR",21,-0.5,20.5,"Rank of matched pair among dR sorted pair","Events");
  
  TrackPairdR_Crt =HConfig.GetTH1D(Name+"_TrackPairdR_Crt","TrackPairdR_Crt",100,0,1.0,"Pair dR - Matched Tracks","Entries");
  TrackPairdR_Bkg =HConfig.GetTH1D(Name+"_TrackPairdR_Bkg","TrackPairdR_Bkg",100,0,1.0,"Pair dR - Isolation Tracks","Entries");
  
  IsoTrackMatchedToSV_1=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_1","IsoTrackMatchedToSV_1",2,-0.5,1.5,"Whether atleast 2 iso track matched to SV track","Entries");
  IsoTrackMatchedToSV_TwoMatched=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_TwoMatched","IsoTrackMatchedToSV_TwoMatched",2,-0.5,1.5,"Whether 2 correct iso track matched to SV","Entries");
  IsoTrackMatchedToSV_TwoCrtIso=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_TwoCrtIso","IsoTrackMatchedToSV_TwoCrtIso",2,-0.5,1.5,"Whether 2 correct iso track matched to SV","Entries");
  IsoTrackMatchedToSV_MassMatch=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_MassMatch","IsoTrackMatchedToSV_MassMatch",2,-0.5,1.5,"Whether 2 correct iso track matched to SV and invariant mass matches","Entries");
  IsoTrackMatchedToSV_MassMatch1=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_MassMatch1","IsoTrackMatchedToSV_MassMatch1",2,-0.5,1.5,"Whether 2 correct iso track matched to SV and invariant mass matches","Entries");
  IsoTrackMatchedToSV_CombMatch=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_CombMatch","IsoTrackMatchedToSV_CombMatch",2,-0.5,1.5,"Whether 2 sorted iso track matched to SV","Entries");
  
  IsoTrackMatchedToSV_ThreeMassMatch=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_ThreeMassMatch","IsoTrackMatchedToSV_ThreeMassMatch",2,-0.5,1.5,"Whether 3 correct iso track matched to SV and invariant mass matches","Entries");
  IsoTrackMatchedToSV_ThreeMassMatch1=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_ThreeMassMatch1","IsoTrackMatchedToSV_ThreeMassMatch1",2,-0.5,1.5,"Whether 3 correct iso track matched to SV and invariant mass matches","Entries");
  
  IsoTrackMatchedToSV_Count=HConfig.GetTH1D(Name+"_IsoTrackMatchedToSV_Count","IsoTrackMatchedToSV_Count",8,-0.5,7.5,"No of iso tracks matched to SV","Entries");
  
  CombMatch_Avg1 =HConfig.GetTH1D(Name+"_CombMatch_Avg1","CombMatch_Avg1",100,0,10.0,"Value of combination var after matching","Entries");
  CombMatch_Avg2 =HConfig.GetTH1D(Name+"_CombMatch_Avg2","CombMatch_Avg2",100,0,10.0,"Value of combination var after matching","Entries");
  
  Angle_SVPV_iSVSV =HConfig.GetTH1D(Name+"_Angle_SVPV_iSVSV","Angle_SVPV_iSVSV",100,0,3.14159265,"Angle between SVPV and iSVSV","Events");
  Angle_SVPV_isvSV =HConfig.GetTH1D(Name+"_Angle_SVPV_isvSV","Angle_SVPV_isvSV",100,0,3.14159265,"Angle between SVPV and isvSV","Events");
  
  iSVSV_Distance =HConfig.GetTH1D(Name+"_iSVSV_Distance","iSVSV_Distance",100,0,2,"Dist between SV and iSV","Events");
  iSVSV_Distance_Sig =HConfig.GetTH1D(Name+"_iSVSV_Distance_Sig","iSVSV_Distance_Sig",100,0,10,"Significance of dist between SV and iSV","Events");
  
  isvSV_Distance =HConfig.GetTH1D(Name+"_isvSV_Distance","isvSV_Distance",100,0,2,"Dist between SV and isv","Events");
  isvSV_Distance_Sig =HConfig.GetTH1D(Name+"_isvSV_Distance_Sig","isvSV_Distance_Sig",100,0,10,"Significance of dist between SV and isv","Events");
  
  InvMass2ProngMatched =HConfig.GetTH1D(Name+"_InvMass2ProngMatched","InvMass2ProngMatched",50,0,5.0,"Invariant Mass of two tracks","Events");
  InvMass2ProngMatchedSV =HConfig.GetTH1D(Name+"_InvMass2ProngMatchedSV","InvMass2ProngMatchedSV",50,0,5.0,"Invariant Mass of SV","Events");
  InvMass2ProngNotMatched =HConfig.GetTH1D(Name+"_InvMass2ProngNotMatched","InvMass2ProngNotMatched",50,0,5.0,"Invariant Mass of SV, unmatched","Events");
  
  InvMass3ProngMatched =HConfig.GetTH1D(Name+"_InvMass3ProngMatched","InvMass3ProngMatched",50,0,5.0,"Invariant Mass of three tracks","Events");
  InvMass3ProngMatchedSV =HConfig.GetTH1D(Name+"_InvMass3ProngMatchedSV","InvMass3ProngMatchedSV",50,0,5.0,"Invariant Mass of SV","Events");
  InvMass3ProngNotMatched =HConfig.GetTH1D(Name+"_InvMass3ProngNotMatched","InvMass3ProngNotMatched",50,0,5.0,"Invariant Mass of SV, unmatched","Events");
  
  InvMassTotal =HConfig.GetTH1D(Name+"_InvMassTotal","InvMassTotal",50,0,5.0,"Invariant Mass of tau+2 tracks","Events");
  InvMassTotal1 =HConfig.GetTH1D(Name+"_InvMassTotal1","InvMassTotal1",50,0,5.0,"Invariant Mass of tau+2 tracks","Events");
  InvMassTotal2 =HConfig.GetTH1D(Name+"_InvMassTotal2","InvMassTotal2",50,0,5.0,"Invariant Mass of 5 tracks minus 2","Events");
  
  SVSize=HConfig.GetTH1D(Name+"_SVSize","SVSize",11,-0.5,10.5,"No of SV","Entries");
  SVNoOfTracksMatched=HConfig.GetTH1D(Name+"_SVNoOfTracksMatched","SVNoOfTracksMatched",11,-0.5,10.5,"No of tracks in matched SV","Entries");
  SVNoOfTracksMatchedThree=HConfig.GetTH1D(Name+"_SVNoOfTracksMatchedThree","SVNoOfTracksMatchedThree",11,-0.5,10.5,"No of tracks in matched SV","Entries");
  SVNoOfTracksUnmatched=HConfig.GetTH1D(Name+"_SVNoOfTracksUnmatched","SVNoOfTracksUnmatched",11,-0.5,10.5,"No of tracks in unmatched SV","Entries");
  
  SVCollectionNoOfSignalMu=HConfig.GetTH1D(Name+"_SVCollectionNoOfSignalMu","SVCollectionNoOfSignalMu",11,-0.5,10.5,"No of signal muons in SV","Entries");
  SVCollectionNoOfNeither=HConfig.GetTH1D(Name+"_SVCollectionNoOfNeither","SVCollectionNoOfNeither",11,-0.5,10.5,"No of tracks in SV; neither signal nor crt track","Entries");
  
  SVCollectionNoOfSignalMu_if1=HConfig.GetTH1D(Name+"_SVCollectionNoOfSignalMu_if1","SVCollectionNoOfSignalMu_if1",11,-0.5,10.5,"No of signal muons in SV","Entries");
  SVCollectionNoOfNeither_if1=HConfig.GetTH1D(Name+"_SVCollectionNoOfNeither_if1","SVCollectionNoOfNeither_if1",11,-0.5,10.5,"No of tracks in SV; neither signal nor crt track","Entries");
  
  SVCollectionNoOfSignalMu_ifmore1=HConfig.GetTH1D(Name+"_SVCollectionNoOfSignalMu_ifmore1","SVCollectionNoOfSignalMu_ifmore1",11,-0.5,10.5,"No of signal muons in SV","Entries");
  SVCollectionNoOfNeither_ifmore1=HConfig.GetTH1D(Name+"_SVCollectionNoOfNeither_ifmore1","SVCollectionNoOfNeither_ifmore1",11,-0.5,10.5,"No of tracks in SV; neither signal nor crt track","Entries");
  
  SVCollectionNoOfCrt=HConfig.GetTH1D(Name+"_SVCollectionNoOfCrt","SVCollectionNoOfCrt",11,-0.5,10.5,"No of correct tracks in SV","Entries");
  SVCollectionNoOfCrt_if3=HConfig.GetTH1D(Name+"_SVCollectionNoOfCrt_if3","SVCollectionNoOfCrt_if3",11,-0.5,10.5,"No of correct tracks in SV","Entries");
  
  Whether_Lowest_Chi2_is_Correct_2iso=HConfig.GetTH1D(Name+"_Whether_Lowest_Chi2_is_Correct_2iso","Whether_Lowest_Chi2_is_Correct_2iso",2,-0.5,1.5,"Whether 2 vertex with lowest chi2 is correct","Entries");
  Rank_Correct_2iso_Chi2 =HConfig.GetTH1D(Name+"_Rank_Correct_2iso_Chi2","Rank_Correct_2iso_Chi2",22,-0.5,21.5,"Rank of matched iso pair by chi2","Events");
  Rank_Correct_2iso3mu_Chi2 =HConfig.GetTH1D(Name+"_Rank_Correct_2iso3mu_Chi2","Rank_Correct_2iso3mu_Chi2",22,-0.5,21.5,"Rank of matched iso pair and 3mu by chi2","Events");
  Rank_Correct_2iso_PairAngle =HConfig.GetTH1D(Name+"_Rank_Correct_2iso_PairAngle","Rank_Correct_2iso_PairAngle",22,-0.5,21.5,"Rank of matched iso pair by angle","Events");
  
  Rank_Correct_1iso3mu_Chi2 =HConfig.GetTH1D(Name+"_Rank_Correct_1iso3mu_Chi2","Rank_Correct_1iso3mu_Chi2",8,-0.5,7.5,"Rank of matched iso and 3mu by chi2","Events");
  
  Rank_Correct_3iso_Chi2 =HConfig.GetTH1D(Name+"_Rank_Correct_3iso_Chi2","Rank_Correct_3iso_Chi2",36,-0.5,35.5,"Rank of matched iso triplet by chi2","Events");
  Rank_Correct_3iso3mu_Chi2 =HConfig.GetTH1D(Name+"_Rank_Correct_3iso3mu_Chi2","Rank_Correct_3iso3mu_Chi2",36,-0.5,35.5,"Rank of matched iso triplet and 3mu by chi2","Events");
  
  TauEnergyVsPairdR=HConfig.GetTH2D(Name+"_TauEnergyVsPairdR","TauEnergyVsPairdR",110,10.0,65.0,40,0,0.2,"E","dR");
  EnergyVsPairdR=HConfig.GetTH2D(Name+"_EnergyVsPairdR","EnergyFractionVsPairdR",100,0.0,1.0,80,0,0.8,"E fraction, Tau","dR");
  EnergyVsPairdR_1=HConfig.GetTH2D(Name+"_EnergyVsPairdR_1","EnergyFractionVsPairdR",100,0.0,1.0,80,0,0.8,"E fraction, Iso Pair","dR");
  EnergyVsPairdR_2=HConfig.GetTH2D(Name+"_EnergyVsPairdR_2","dR vs Angle Sum",100,0.0,6.0,80,0,0.8,"Sum of Angles","dR");
  EnergyVsPairdR_Circular=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Circular","Angular Distribution in COM frame, wrt B Meson",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Tau=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Tau","Angular Distribution in COM frame, wrt Tau",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Circular_reco=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Circular_reco","Reco angular Distribution in COM frame, wrt B Meson",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Tau_reco=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Tau_reco","Reco angular Distribution in COM frame, wrt Tau",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Circular_dR=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Circular_dR","Angular Distribution in COM frame, Vtx",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Circular_Incorrect=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Circular_Incorrect","EnergyVsPairdR_Circular_Incorrect",40,-4.,4.,40,0.,4.,"#phi","#theta");
  EnergyVsPairdR_Circular_Incorrect_reco=HConfig.GetTH2D(Name+"_EnergyVsPairdR_Circular_Incorrect_reco","EnergyVsPairdR_Circular_Incorrect_reco",40,-4.,4.,40,0.,4.,"#phi","#theta");
  
  GammaVs_Tau_Energy=HConfig.GetTH2D(Name+"_GammaVs_Tau_Energy","B Meson Gamma vs Tau Energy",100,0.0,100.0,80,0,80.0,"E, Tau","#Gamma");
  GammaVs_Pair_Energy=HConfig.GetTH2D(Name+"_GammaVs_Pair_Energy","B Meson Gamma vs Pair Energy",100,0.0,100.0,80,0,80.0,"E, Iso Pair","#Gamma");
  
  PairAngle_2Cor =HConfig.GetTH1D(Name+"_PairAngle_2Cor","PairAngle_2Cor",80,0,1.6,"Pair dR, 2 Correct","Events");
  PairAngle_1Cor =HConfig.GetTH1D(Name+"_PairAngle_1Cor","PairAngle_1Cor",80,0,1.6,"Pair dR, 1 Correct","Events");
  PairAngle_0Cor =HConfig.GetTH1D(Name+"_PairAngle_0Cor","PairAngle_0Cor",80,0,1.6,"Pair dR, 0 Correct","Events");
  
  PairMass_2Cor =HConfig.GetTH1D(Name+"_PairMass_2Cor","PairMass_2Cor",50,0,5.0,"Pair Mass, 2 Correct","Events");
  PairMass_1Cor =HConfig.GetTH1D(Name+"_PairMass_1Cor","PairMass_1Cor",50,0,5.0,"Pair Mass, 1 Correct","Events");
  PairMass_0Cor =HConfig.GetTH1D(Name+"_PairMass_0Cor","PairMass_0Cor",50,0,5.0,"Pair Mass, 0 Correct","Events");
  
  Angle_SVPV_BMeson =HConfig.GetTH1D(Name+"_Angle_SVPV_BMeson","Angle_SVPV_BMeson",100,0.,0.25,"Angle between SVPV and BMeson","Events");
  Angle_Reco_Comparison =HConfig.GetTH1D(Name+"_Angle_Reco_Comparison","Angle_Reco_Comparison",100,0.,2.5,"Angle between boosted tracks","Events");
  Angle_Tau_Vtx =HConfig.GetTH1D(Name+"_Angle_Tau_Vtx","Angle_Tau_Vtx",100,0,3.14159265,"Angle between tau and vtx dirn.","Events");
  
  Angle_dVtx_Pair_Crt =HConfig.GetTH1D(Name+"_Angle_dVtx_Pair_Crt","Angle_dVtx_Pair_Crt",100,0,3.14159265,"Angle between dVertex and Pair dirn.","Events");
  Angle_dVtx_Pair_InCrt =HConfig.GetTH1D(Name+"_Angle_dVtx_Pair_InCrt","Angle_dVtx_Pair_InCrt",100,0,3.14159265,"Angle between dVertex and Pair dirn.","Events");
  
  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  SignalVertexSelector::Store_ExtraDist(){ 

  Extradist1d.push_back(&WhetherTau3Mu);
  Extradist1d.push_back(&WhetherdRMatch);
  Extradist1d.push_back(&IsoTrackToMCdR01);
  Extradist1d.push_back(&IsoTrackToMCdR08);
  Extradist1d.push_back(&IsoTrackToMCAngle01);
  
  Extradist1d.push_back(&NumberOfFS_ChargedParticles);
  Extradist1d.push_back(&NumberOfFS_ChargedParticles_RecoMatch);
  
  Extradist1d.push_back(&NumberOfRecoChargedParticlesIfMC1);
  Extradist1d.push_back(&NumberOfRecoChargedParticlesIfMC2);
  Extradist1d.push_back(&NumberOfRecoChargedParticlesIfMC3);
  
  Extradist1d.push_back(&TwoProngInvariantMassReco);
  Extradist1d.push_back(&TwoProngInvariantMassMC);
  
  Extradist1d.push_back(&ThreeProngInvariantMassReco);
  Extradist1d.push_back(&ThreeProngInvariantMassMC);
  
  Extradist1d.push_back(&TwoProngInvariantMassReco005);
  Extradist1d.push_back(&TwoProngInvariantMassMC005);
  
  Extradist1d.push_back(&ThreeProngInvariantMassReco005);
  Extradist1d.push_back(&ThreeProngInvariantMassMC005);
  
  Extradist1d.push_back(&TwoProngTrackPt);
  Extradist1d.push_back(&TwoProngTrack2Pt);
  
  Extradist1d.push_back(&ThreeProngTrackPt);
  Extradist1d.push_back(&ThreeProngTrack2Pt);
  Extradist1d.push_back(&ThreeProngTrack3Pt);
  
  Extradist1d.push_back(&TrackToTauDr2Prong);
  Extradist1d.push_back(&TrackToTauDr3Prong);
  Extradist1d.push_back(&TrackToTauDrAll);
  
  Extradist1d.push_back(&TwoProngTrackEta);
  Extradist1d.push_back(&TwoProngTrack2Eta);
  
  Extradist2d.push_back(&dR_vs_dP);
  Extradist2d.push_back(&dRmin_sum_vs_InvariantMass_2prong);
  Extradist2d.push_back(&dRmin_sum_vs_InvariantMass_3prong);
  
  Extradist1d.push_back(&NoOfIsoTracks2Prong);
  Extradist1d.push_back(&NoOfIsoTracks3Prong);
  
  Extradist2d.push_back(&InvMass2_vs_pdgid);
  Extradist2d.push_back(&InvMass3_vs_pdgid);
  
  Extradist1d.push_back(&RankMatchedTrackpT);
  Extradist1d.push_back(&RankMatchedTrackdR);
  Extradist1d.push_back(&RankMatchedTrackdR_trim);
  
  Extradist1d.push_back(&RankMatchedTrackAvgDiff);
  Extradist1d.push_back(&RankMatchedTrackCombn);
  
  
  Extradist1d.push_back(&var_All_7_Iso_dR);
  Extradist1d.push_back(&var_Correct_Iso_dR);
  
  Extradist1d.push_back(&var_All_7_Iso_AvgDiff);
  Extradist1d.push_back(&var_Correct_Iso_AvgDiff);
  
  Extradist1d.push_back(&var_All_7_Iso_Avg);
  Extradist1d.push_back(&var_Correct_Iso_Avg);
  
  Extradist1d.push_back(&var_All_7_Iso_Combn);
  Extradist1d.push_back(&var_Correct_Iso_Combn);
  
  Extradist1d.push_back(&var_All_7_Iso_Trial);
  Extradist1d.push_back(&var_Correct_Iso_Trial);
  
  Extradist1d.push_back(&RankMatchedTrackdR_cut);
  
  Extradist1d.push_back(&RankMatchedTrackPairdR);
  Extradist1d.push_back(&TrackPairdR_Crt);
  Extradist1d.push_back(&TrackPairdR_Bkg);
  
  Extradist1d.push_back(&IsoTrackMatchedToSV_1);
  Extradist1d.push_back(&IsoTrackMatchedToSV_TwoMatched);
  Extradist1d.push_back(&IsoTrackMatchedToSV_TwoCrtIso);
  Extradist1d.push_back(&IsoTrackMatchedToSV_MassMatch);
  Extradist1d.push_back(&IsoTrackMatchedToSV_MassMatch1);
  Extradist1d.push_back(&IsoTrackMatchedToSV_CombMatch);
  
  Extradist1d.push_back(&IsoTrackMatchedToSV_Count);
  
  Extradist1d.push_back(&CombMatch_Avg1);
  Extradist1d.push_back(&CombMatch_Avg2);
  
  Extradist1d.push_back(&Angle_SVPV_iSVSV);
  Extradist1d.push_back(&Angle_SVPV_isvSV);
  
  Extradist1d.push_back(&IsoTrackMatchedToSV_ThreeMassMatch);
  Extradist1d.push_back(&IsoTrackMatchedToSV_ThreeMassMatch1);
  
  Extradist1d.push_back(&iSVSV_Distance);
  Extradist1d.push_back(&iSVSV_Distance_Sig);
  
  Extradist1d.push_back(&isvSV_Distance);
  Extradist1d.push_back(&isvSV_Distance_Sig);
  
  Extradist1d.push_back(&InvMass2ProngMatched);
  Extradist1d.push_back(&InvMass2ProngMatchedSV);
  Extradist1d.push_back(&InvMass2ProngNotMatched);
  
  Extradist1d.push_back(&SVSize);
  Extradist1d.push_back(&SVNoOfTracksMatched);
  Extradist1d.push_back(&SVNoOfTracksMatchedThree);
  Extradist1d.push_back(&SVNoOfTracksUnmatched);
  
  Extradist1d.push_back(&InvMass3ProngMatched);
  Extradist1d.push_back(&InvMass3ProngMatchedSV);
  Extradist1d.push_back(&InvMass3ProngNotMatched);
  
  Extradist1d.push_back(&InvMassTotal);
  Extradist1d.push_back(&InvMassTotal1);
  Extradist1d.push_back(&InvMassTotal2);
  
  Extradist1d.push_back(&SVCollectionNoOfSignalMu);
  Extradist1d.push_back(&SVCollectionNoOfNeither);
  
  Extradist1d.push_back(&SVCollectionNoOfSignalMu_if1);
  Extradist1d.push_back(&SVCollectionNoOfNeither_if1);
  
  Extradist1d.push_back(&SVCollectionNoOfSignalMu_ifmore1);
  Extradist1d.push_back(&SVCollectionNoOfNeither_ifmore1);
  
  Extradist1d.push_back(&SVCollectionNoOfCrt);
  Extradist1d.push_back(&SVCollectionNoOfCrt_if3);
  
  Extradist1d.push_back(&Whether_Lowest_Chi2_is_Correct_2iso);
  
  Extradist1d.push_back(&Rank_Correct_1iso3mu_Chi2);
  Extradist1d.push_back(&Rank_Correct_2iso_Chi2);
  Extradist1d.push_back(&Rank_Correct_2iso3mu_Chi2);
  Extradist1d.push_back(&Rank_Correct_2iso_PairAngle);
  Extradist1d.push_back(&Rank_Correct_3iso_Chi2);
  Extradist1d.push_back(&Rank_Correct_3iso3mu_Chi2);
  
  Extradist2d.push_back(&TauEnergyVsPairdR);
  Extradist2d.push_back(&EnergyVsPairdR);
  Extradist2d.push_back(&EnergyVsPairdR_1);
  Extradist2d.push_back(&EnergyVsPairdR_2);
  Extradist2d.push_back(&EnergyVsPairdR_Circular);
  Extradist2d.push_back(&EnergyVsPairdR_Tau);
  Extradist2d.push_back(&EnergyVsPairdR_Circular_reco);
  Extradist2d.push_back(&EnergyVsPairdR_Tau_reco);
  Extradist2d.push_back(&EnergyVsPairdR_Circular_dR);
  Extradist2d.push_back(&EnergyVsPairdR_Circular_Incorrect);
  Extradist2d.push_back(&EnergyVsPairdR_Circular_Incorrect_reco);
  
  Extradist2d.push_back(&GammaVs_Tau_Energy);
  Extradist2d.push_back(&GammaVs_Pair_Energy);
  
  Extradist1d.push_back(&PairAngle_2Cor);
  Extradist1d.push_back(&PairAngle_1Cor);
  Extradist1d.push_back(&PairAngle_0Cor);
  
  Extradist1d.push_back(&PairMass_2Cor);
  Extradist1d.push_back(&PairMass_1Cor);
  Extradist1d.push_back(&PairMass_0Cor);
  
  Extradist1d.push_back(&Angle_SVPV_BMeson);
  Extradist1d.push_back(&Angle_Reco_Comparison);
  Extradist1d.push_back(&Angle_Tau_Vtx);
  
  Extradist1d.push_back(&Angle_dVtx_Pair_Crt);
  Extradist1d.push_back(&Angle_dVtx_Pair_InCrt);


}


void  SignalVertexSelector::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  //  std::cout<<" id   "<< id << std::endl;
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
    //    std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("DoubleMu3_TkMu_DsTau3Mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}
  }

  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //    std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
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


  //  if(DoubleMuFired) value.at(L1T)=1;

  //  std::cout<<"  "<< value.at(L1T) << "  "<<value.at(HLT)  << std::endl;


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
  
  unsigned int final_idx=signal_idx;


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


  //  if(id!=1)  std::cout<<" id:   "<< id << "  NMCSignalParticles  "<< Ntp->NMCSignalParticles() << "  NMCTaus   "<< Ntp->NMCTaus() << std::endl;


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
    
    

    //    std::cout <<"  " << Muon1LV.M() <<"   " <<Muon2LV.M() <<"   "  <<  Muon3LV.M() <<std::endl;

    vector<unsigned int> idx_vec;
    idx_vec.push_back(Muon_index_1);
    idx_vec.push_back(Muon_index_2);
    idx_vec.push_back(Muon_index_3);

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    //*** With such sorting pT(ss1) > pT(ss2)
    TLorentzVector MuonOS  = Ntp->Muon_P4(os_mu_idx);  
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
    
    
    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    
    TVector3 Tau_Vector = TauLV.Vect();
    TVector3 SVPV_Vector = Ntp->Vertex_Signal_KF_pos(signal_idx) - Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
    TVector3 SVPV_Vector_Mag = SVPV_Vector;
    SVPV_Vector_Mag.SetMag(Tau_Vector.Mag());
    TVector3 Tau_Vector_Reflection = Tau_Vector;
    Tau_Vector_Reflection.Rotate(TMath::Pi()/2, SVPV_Vector);
    TLorentzVector TauLV_Reflection(Tau_Vector_Reflection,TauLV.E());
    TLorentzVector SVPV_LV(SVPV_Vector_Mag,TauLV.E());
    
    
    //TLorentzVector MuonOS = Muon1LV;
    
    TLorentzVector MuonOSReassigned = MuonOS;//reassign the masses to be similar to a kaon
    MuonOSReassigned.SetE(sqrt(MuonOS.Px()*MuonOS.Px()+MuonOS.Py()*MuonOS.Py()+MuonOS.Pz()*MuonOS.Pz()+0.493677*0.493677));
    TLorentzVector MuonSS1Reassigned = MuonSS1;
    MuonSS1Reassigned.SetE(sqrt(MuonSS1.Px()*MuonSS1.Px()+MuonSS1.Py()*MuonSS1.Py()+MuonSS1.Pz()*MuonSS1.Pz()+0.493677*0.493677));
    TLorentzVector MuonSS2Reassigned = MuonSS2;
    MuonSS2Reassigned.SetE(sqrt(MuonSS2.Px()*MuonSS2.Px()+MuonSS2.Py()*MuonSS2.Py()+MuonSS2.Pz()*MuonSS2.Pz()+0.493677*0.493677));
    
    MuOSSS1InvariantMassBeforeMVA.at(t).Fill((MuonOSReassigned+MuonSS1Reassigned).M(),1);
    MuOSSS2InvariantMassBeforeMVA.at(t).Fill((MuonOSReassigned+MuonSS2Reassigned).M(),1);
    
    // This is to print out selected events content
    /*
    if(id ==60 ||  id ==90){// or id == 40){
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
    
    int NoOfTracksAfterdRCut = 7;
    
    std::vector<int> b_meson_full_childidx;//All descendents of B meson except for the tau and the three muons
    std::vector<int> b_meson_full_childidx_FS;//All descendents of B meson which are charged final state particles (except t3mu)
    std::vector<int> b_meson_full_childidx_FS_RecoMC;//All descendents of B meson which are charged final state particles and have been reconstructed within a certain dR (except t3mu)
    
    std::vector<float> pts_reco; // Transverse momenta
    std::vector<float> phis_reco;// signed values of phi
    std::vector<float> etas_reco;// signed values of eta
    std::vector<int> charges_reco;// indirectly calculated charges
    
    std::vector<float> es;// particle energy
    std::vector<float> p1;// particle px
    std::vector<float> p2;
    std::vector<float> p3;
    
    std::vector<float> mc_pdgid;
    
    std::vector<float> es_reco;// particle energy
    std::vector<float> p1_reco;// particle px
    std::vector<float> p2_reco;
    std::vector<float> p3_reco;
    
    std::vector<float> dR_min_reco;
    
    //Do everything with dR<0.005
    std::vector<int> b_meson_full_childidx_FS_RecoMC005;
    
    std::vector<float> es005;// particle energy
    std::vector<float> p1005;// particle px
    std::vector<float> p2005;
    std::vector<float> p3005;
    
    std::vector<float> es_reco005;// particle energy
    std::vector<float> p1_reco005;// particle px
    std::vector<float> p2_reco005;
    std::vector<float> p3_reco005;
    
    std::vector<float> dR_min_reco005;
    
    std::vector<float> dR_To_Tau_Prong;
    
    std::vector<float> MatchedIsoTrackNo; // Just the index, 'j' of the matched isolation track
    
    TLorentzVector LV_BMeson; //LV of B Meson from MC
    
    //Trying to figure out if there are final state particles coming from b signal that match the isolation tracks
    for(int gen_part_index=0; gen_part_index < Ntp->NMCParticles(); gen_part_index++){        
          //std::cout<<"Particle PDGID is:"<< Ntp->MCParticle_pdgid(gen_part_index) << std::endl;
          //std::cout<<"Index of the mother particle is:"<< Ntp->MCParticle_midx(gen_part_index) << std::endl;
          if(Ntp->MCParticle_midx(gen_part_index)>-0.5){
            
            //Trying to figure out if there are final state particles coming from b signal that match the isolation tracks
            
            if(abs(Ntp->MCParticle_pdgid(gen_part_index)) == 15){
              bool Three_Children(false);
              
              std::vector<int> tau_childpdgid = Ntp->MCParticle_childpdgid(gen_part_index);
              int count_muon(0);
              
              //for(int i : tau_childpdgid){std::cout<< i << std::endl;}
              
              for(int i : tau_childpdgid){if(abs(i)==13){count_muon=count_muon+1;}}//counts number of muons for which tau is a parent
              
              if(count_muon>2){// Tau with three muon children found
                Three_Children=true;
                
                std::vector<int> b_meson_childidx = Ntp->MCParticle_childidx(Ntp->MCParticle_midx(gen_part_index)); //all 'id's of direct children of parent of tau. Kind of useless
                
                //Create a vector of indices of all children of parent of tau
                int tau_parent_idx = Ntp->MCParticle_midx(gen_part_index);
                LV_BMeson=Ntp->MCParticle_p4(tau_parent_idx);
                
                for(int all_index=0; all_index < Ntp->NMCParticles(); all_index++){//used to get indices of all children of parent of tau, excluding tau and 3 muons
                  
                  if(Ntp->MCParticle_midx(all_index)>-0.5&&all_index!=gen_part_index&&Ntp->MCParticle_midx(all_index)!=gen_part_index){
                    if(Ntp->MCParticle_midx(all_index)==tau_parent_idx){
                      b_meson_full_childidx.push_back(all_index);
                    }
                    if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index))>-0.5){
                      if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index))==tau_parent_idx){
                        b_meson_full_childidx.push_back(all_index);
                      }
                      if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index)))>-0.5){
                        if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index)))==tau_parent_idx){
                          b_meson_full_childidx.push_back(all_index);
                        }
                        if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index))))>-0.5){
                          if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index))))==tau_parent_idx){
                            b_meson_full_childidx.push_back(all_index);
                          }
                          if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index)))))>-0.5){
                            if(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(Ntp->MCParticle_midx(all_index)))))==tau_parent_idx){
                              b_meson_full_childidx.push_back(all_index);
                            }
                          }
                        }
                      }
                    }
                  }
                  
                }// end of all_index for loop
                
                //std::cout<<"Children of parent of tau with index: "<< gen_part_index << std::endl;
                //for(int i : b_meson_full_childidx){std::cout<<"Particle index: "<< i <<" with pdgid "<< Ntp->MCParticle_pdgid(i)<< " charge "<<Ntp->MCParticle_charge(i)<<" status "<< Ntp->MCParticle_status(i)<< std::endl;}
                
                
                for(int i : b_meson_full_childidx){if(Ntp->MCParticle_status(i)==1&&abs(Ntp->MCParticle_charge(i))>0){//here, try to match charged MC final state particles to isolation tracks
                  b_meson_full_childidx_FS.push_back(i);// Fill all charged final state particles we're interested in
                  double dR_min(199.0);
                  double Ang_min(99.0);
                  bool dR_Match(false);
                  bool dR_Match005(false);
                  
                  
                  TLorentzVector TrackLV_min;
                  TLorentzVector FinalStateParticleLV_min;
                  
                  float TrackNoMin;// Just the track index
                  
                  for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){//loop over isolation tracks
                    TLorentzVector TrackLV = Ntp->IsolationTrack_p4(signal_idx,j);
                    TLorentzVector FinalStateParticleLV = Ntp->MCParticle_p4(i);
                    
                    TVector3 TrackLV1 = TrackLV.Vect();
                    TVector3 FinalStateParticleLV1 = FinalStateParticleLV.Vect();
                    
                    double dR1=fabs(FinalStateParticleLV.DeltaR(TrackLV));
                    if(dR1<dR_min&&(std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), j))==0){
                      dR_min=dR1;
                      TrackLV_min=TrackLV;
                      FinalStateParticleLV_min=FinalStateParticleLV;
                      TrackNoMin=j;
                    }
                    
                    double AngleDiff=fabs(TrackLV1.Angle(FinalStateParticleLV1));
                    (AngleDiff<Ang_min)?(Ang_min=AngleDiff):(Ang_min=Ang_min);
                    
                  }// end of j for loop
                  
                  IsoTrackToMCdR01.at(t).Fill(dR_min,1);
                  IsoTrackToMCdR08.at(t).Fill(dR_min,1);
                  IsoTrackToMCAngle01.at(t).Fill(Ang_min,1);
                  
                  dR_vs_dP.at(t).Fill(dR_min,abs((TrackLV_min-FinalStateParticleLV_min).Perp()),1);
                  
                  (dR_min<0.015)?(dR_Match=true):(dR_Match=false);
                  (dR_min<0.005)?(dR_Match005=true):(dR_Match005=false);
                  
                  if(dR_Match005){// Fill all charged final state particles we're interested in that have been reconstructed
                    b_meson_full_childidx_FS_RecoMC005.push_back(i);
                    
                    es005.push_back(FinalStateParticleLV_min.E());// particle energy
                    p1005.push_back(FinalStateParticleLV_min.Px());// particle px
                    p2005.push_back(FinalStateParticleLV_min.Py());
                    p3005.push_back(FinalStateParticleLV_min.Pz());
                    
                    es_reco005.push_back(TrackLV_min.E());// particle energy
                    p1_reco005.push_back(TrackLV_min.Px());// particle px
                    p2_reco005.push_back(TrackLV_min.Py());
                    p3_reco005.push_back(TrackLV_min.Pz());
                    
                    dR_min_reco005.push_back(dR_min);
                    //charges_reco.push_back(TrackLV_min);
                  }
                  
                  if(dR_Match){// Fill all charged final state particles we're interested in that have been reconstructed
                    b_meson_full_childidx_FS_RecoMC.push_back(i);
                    
                    es.push_back(FinalStateParticleLV_min.E());// particle energy
                    p1.push_back(FinalStateParticleLV_min.Px());// particle px
                    p2.push_back(FinalStateParticleLV_min.Py());
                    p3.push_back(FinalStateParticleLV_min.Pz());
                    
                    es_reco.push_back(TrackLV_min.E());// particle energy
                    p1_reco.push_back(TrackLV_min.Px());// particle px
                    p2_reco.push_back(TrackLV_min.Py());
                    p3_reco.push_back(TrackLV_min.Pz());
                    
                    pts_reco.push_back(TrackLV_min.Perp()); // Transverse momenta
                    phis_reco.push_back(TrackLV_min.Phi());// signed values of phi
                    etas_reco.push_back(TrackLV_min.Eta());// signed values of eta
                    
                    dR_min_reco.push_back(dR_min);
                    
                    mc_pdgid.push_back(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i))));
                    
                    dR_To_Tau_Prong.push_back(TrackLV_min.DeltaR(TauLV));
                    //charges_reco.push_back(TrackLV_min);
                    
                    MatchedIsoTrackNo.push_back(TrackNoMin);
                  }
                  
                  WhetherdRMatch.at(t).Fill(dR_Match,1);
                }}//end of i for loop?
                
              }
              //if(!Three_Children){std::cout<<"Tau not found with index: "<< gen_part_index << std::endl;}
              WhetherTau3Mu.at(t).Fill(Three_Children,1);
            }//if(abs(Ntp->MCParticle_pdgid(gen_part_index)) == 15)
            
            
            
          }//end of (gen_part_index)>-0.5
          
          
          
    }//end of gen_part_index for loop
    
    // Attempt at reconstructing B-Meson LV from estimate of gamma and direction
    double B_Gamma = 0.341265 * TauLV.E() + 0.883522;
    double B_Beta = sqrt(1-(1/(B_Gamma*B_Gamma)));
    TVector3 B_Vector_reco = ((5.2795*B_Gamma*B_Beta)/SVPV_Vector.Mag())*SVPV_Vector;
    TLorentzVector LV_BMeson_reco(B_Vector_reco,5.2795*B_Gamma);
    
    using Row = vector<double>;
    using Matrix = vector<Row>;
    
    Matrix dR_No;// n x 2 matrix where first column is dR and second column is the index of the isolation track
    Matrix pT_No;// n x 2 matrix where first column is pT and second column is the index of the isolation track
    
    Matrix pT_No_7;// n x 2 matrix where first column is pT and second column is the index of the isolation track. These are the 7 tracks closest to tau by dR
    
    for(int j=0;j<Ntp->NIsolationTrack(signal_idx);j++){//loop over isolation tracks
      TLorentzVector TrackLV = Ntp->IsolationTrack_p4(signal_idx,j);
      TrackToTauDrAll.at(t).Fill(TrackLV.DeltaR(TauLV));
      
      dR_No.push_back({TrackLV.DeltaR(TauLV),j});
      pT_No.push_back({TrackLV.Perp(),j});
    }
    
    sort( dR_No.begin(), dR_No.end() ); // sort based on first column, lowest first
    sort( pT_No.rbegin(), pT_No.rend() ); // sort based on first column, highest first
    
    //Adding conditions for cuts applied on iso tracks:
    Matrix Selection_Status_OneProng;
    //for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>0;i++){
    //}
    Matrix Selection_Status_TwoProng;
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>1;i++){
      for(int j=0;j<dR_No.size()&&j<NoOfTracksAfterdRCut&&j<i;j++){// combinations of tracks with indices from 0 - 6
        TLorentzVector TrackLV1 = Ntp->IsolationTrack_p4(signal_idx,dR_No[i][1]);
        TLorentzVector TrackLV2 = Ntp->IsolationTrack_p4(signal_idx,dR_No[j][1]);
        
        bool Mass_Cut = ((TrackLV1+TrackLV2).M()<2.0)?true:false;
        bool Angle_Cut = (fabs((TrackLV1.Vect()).Angle(TrackLV2.Vect()))<1.0)?true:false;
        if(Mass_Cut&&Angle_Cut){
          Selection_Status_TwoProng.push_back({true,i,j});
        }
        else{
          Selection_Status_TwoProng.push_back({false,i,j});
        }
      }
    }
    /*
    Matrix Selection_Status_ThreeProng;
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>2;i++){
      for(int j=0;j<dR_No.size()&&j<NoOfTracksAfterdRCut&&j<i;j++){
        for(int k=0;k<dR_No.size()&&k<NoOfTracksAfterdRCut&&k<j;k++){// combinations of tracks with indices from 0 - 6
          if(){
            Selection_Status_TwoProng.push_back({true,i,j,k});
          }
          else{
            Selection_Status_TwoProng.push_back({false,i,j,k});
          }
        }
      }
    }
    */
    
    /*
    for(int i=0;i<pT_No.size();i++){
      for(int j=0;j<MatchedIsoTrackNo.size();j++){// to get the rank
        if(pT_No[i][1]==MatchedIsoTrackNo[j]){
          RankMatchedTrackpT.at(t).Fill(i,1);
        }
      }
    }
    */
    
    //std::cout<<" The dR sorted indices of the isolation tracks are: " << std::endl;
    
    for(int i=0;i<dR_No.size();i++){
      //std::cout<<" dR : " << dR_No[i][0] << " index : " << dR_No[i][1] << std::endl;
      
      for(int j=0;j<MatchedIsoTrackNo.size();j++){// to get the rank
        if(dR_No[i][1]==MatchedIsoTrackNo[j]){
          RankMatchedTrackdR.at(t).Fill(i,1);
          RankMatchedTrackdR_trim.at(t).Fill(i,1);
        }
      }
      
    }
    
    Matrix Distance_SV_Avg;// Average of distance of common vertex (btw iso track and signal muons) from the signal vertex
    Matrix Distance_SV_Sum;// Sum of positions common vertex (btw iso track and signal muons) from the signal vertex
    Matrix Distance_Difference_Avg;
    Matrix Chi2_Avg;
    Matrix Distance_Difference_Avg_Over_Distance_SV_Avg;
    Matrix Angles3;
    Matrix PosSV_Avg_1D_Perp;// projecting the sum of Pos1SV,Pos2SV,Pos3SV onto the tau axis
    Matrix PosSV_Avg_1D_SVPV;
    Matrix Testing;
    Matrix Testing1;
    Matrix dRtoTauRefl;
    Matrix dRtoSVPV;
    Matrix all_vars;
    
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut;i++){// Selecting pTs of the 7 tracks (indices 0 - 6) that are closest to tau
      TLorentzVector TrackLV = Ntp->IsolationTrack_p4(signal_idx,dR_No[i][1]);
      
      bool Whether_Matched = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[i][1]))==1?true:false;
      
      // Transform everything so that Track1 is along Z-axis (beam axis). +X points from CMS to Saint Genis (phi is measured from that)
      TLorentzVector TrackLV_mod = TrackLV;
      TLorentzVector TauLV_mod = TauLV;
      TrackLV_mod.Boost(-1*LV_BMeson.BoostVector());
      TauLV_mod.Boost(-1*LV_BMeson.BoostVector());
      
      TVector3 Bmeson_mod = LV_BMeson.Vect();
      TVector3 Tau_mod = TauLV_mod.Vect();
      TVector3 Track_mod = TrackLV_mod.Vect();
      
      //Making Bmeson point towards Z axis (and making phi = 0 for the tau: commented)
      double Phi_init = Bmeson_mod.Phi();
      double Theta_init = Bmeson_mod.Theta();
      if(Phi_init >= TMath::Pi()) Phi_init = Phi_init-2*TMath::Pi();
      if(Phi_init <=-TMath::Pi()) Phi_init = Phi_init+2*TMath::Pi();
      Bmeson_mod.RotateZ(-Phi_init);
      Bmeson_mod.RotateY(-Theta_init);
      Tau_mod.RotateZ(-Phi_init);
      Tau_mod.RotateY(-Theta_init);
      Track_mod.RotateZ(-Phi_init);
      Track_mod.RotateY(-Theta_init);
      
      //double Phi_Second = Tau_mod.Phi();
      //if(Phi_Second >= TMath::Pi()) Phi_Second = Phi_Second-2*TMath::Pi();
      //if(Phi_Second <=-TMath::Pi()) Phi_Second = Phi_Second+2*TMath::Pi();
      //Bmeson_mod.RotateZ(-Phi_Second);
      //Tau_mod.RotateZ(-Phi_Second);
      //Track_mod.RotateZ(-Phi_Second);
      
      
      //Making tau point towards Z axis
      TVector3 Bmeson_mod1 = LV_BMeson.Vect();
      TVector3 Tau_mod1 = TauLV_mod.Vect();
      TVector3 Track_mod1 = TrackLV_mod.Vect();
      double Phi_init1 = Tau_mod1.Phi();
      double Theta_init1 = Tau_mod1.Theta();
      if(Phi_init1 >= TMath::Pi()) Phi_init1 = Phi_init1-2*TMath::Pi();
      if(Phi_init1 <=-TMath::Pi()) Phi_init1 = Phi_init1+2*TMath::Pi();
      Tau_mod1.RotateZ(-Phi_init1);
      Tau_mod1.RotateY(-Theta_init1);
      Track_mod1.RotateZ(-Phi_init1);
      Track_mod1.RotateY(-Theta_init1);
      
      // Trying everything Fully Reco
      // Transform everything so that Track1 is along Z-axis (beam axis). +X points from CMS to Saint Genis (phi is measured from that) : fully reco
      TLorentzVector TrackLV_mod_reco = TrackLV;
      TLorentzVector TauLV_mod_reco = TauLV;
      TrackLV_mod_reco.Boost(-1*LV_BMeson_reco.BoostVector());
      TauLV_mod_reco.Boost(-1*LV_BMeson_reco.BoostVector());
      TVector3 Bmeson_mod_reco = LV_BMeson_reco.Vect();
      TVector3 Tau_mod_reco = TauLV_mod_reco.Vect();
      TVector3 Track_mod_reco = TrackLV_mod_reco.Vect();
      
      //Making Bmeson point towards Z axis (and making phi = 0 for the tau: commented): fully reco
      double Phi_init_reco = Bmeson_mod_reco.Phi();
      double Theta_init_reco = Bmeson_mod_reco.Theta();
      if(Phi_init_reco >= TMath::Pi()) Phi_init_reco = Phi_init_reco-2*TMath::Pi();
      if(Phi_init_reco <=-TMath::Pi()) Phi_init_reco = Phi_init_reco+2*TMath::Pi();
      Bmeson_mod_reco.RotateZ(-Phi_init_reco);
      Bmeson_mod_reco.RotateY(-Theta_init_reco);
      Tau_mod_reco.RotateZ(-Phi_init_reco);
      Tau_mod_reco.RotateY(-Theta_init_reco);
      Track_mod_reco.RotateZ(-Phi_init_reco);
      Track_mod_reco.RotateY(-Theta_init_reco);
      
      //Making tau point towards Z axis: fully reco
      TVector3 Bmeson_mod_reco1 = LV_BMeson_reco.Vect();
      TVector3 Tau_mod_reco1 = TauLV_mod_reco.Vect();
      TVector3 Track_mod_reco1 = TrackLV_mod_reco.Vect();
      double Phi_init_reco1 = Tau_mod_reco1.Phi();
      double Theta_init_reco1 = Tau_mod_reco1.Theta();
      if(Phi_init_reco1 >= TMath::Pi()) Phi_init_reco1 = Phi_init_reco1-2*TMath::Pi();
      if(Phi_init_reco1 <=-TMath::Pi()) Phi_init_reco1 = Phi_init_reco1+2*TMath::Pi();
      Tau_mod_reco1.RotateZ(-Phi_init_reco1);
      Tau_mod_reco1.RotateY(-Theta_init_reco1);
      Track_mod_reco1.RotateZ(-Phi_init_reco1);
      Track_mod_reco1.RotateY(-Theta_init_reco1);
      
      
      if(Whether_Matched){
        EnergyVsPairdR_Circular.at(t).Fill(Track_mod.Phi(),Track_mod.Theta());
        EnergyVsPairdR_Tau.at(t).Fill(Track_mod1.Phi(),Track_mod1.Theta());
        
        EnergyVsPairdR_Circular_reco.at(t).Fill(Track_mod_reco.Phi(),Track_mod_reco.Theta());
        EnergyVsPairdR_Tau_reco.at(t).Fill(Track_mod_reco1.Phi(),Track_mod_reco1.Theta());
        
        Angle_Reco_Comparison.at(t).Fill(fabs((Track_mod).Angle(Track_mod_reco)),1);
        
        /*
        if(Track_mod1.Theta()<0.5){//print out events where track is in a similar dirn. as tau
          if(id ==60 ||  id ==90){// or id == 40){
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
        }
        */
      }
      if(!Whether_Matched){
        EnergyVsPairdR_Circular_Incorrect.at(t).Fill(Track_mod.Phi(),Track_mod.Theta());
        EnergyVsPairdR_Circular_Incorrect_reco.at(t).Fill(Track_mod_reco.Phi(),Track_mod_reco.Theta());
      }
      
      //pT_No_7.push_back({TrackLV.Perp(TauLV.Vect()),dR_No[i][1]});
      pT_No_7.push_back({(TrackLV.Vect()).Dot(TauLV.Vect()),dR_No[i][1]});
      
      TVector3 Position1=Ntp->IsolationTrack_VertexWithSignalMuon1Position(signal_idx,dR_No[i][1]);
      TVector3 Position2=Ntp->IsolationTrack_VertexWithSignalMuon2Position(signal_idx,dR_No[i][1]);
      TVector3 Position3=Ntp->IsolationTrack_VertexWithSignalMuon3Position(signal_idx,dR_No[i][1]);
      TVector3 SV_Position = Ntp->Vertex_Signal_KF_pos(signal_idx);// Signal vertex
      TVector3 Tau_Dir = TauLV.Vect();
      TVector3 SVPV_Dir = Ntp->Vertex_Signal_KF_pos(signal_idx) - Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
      
      TVector3 Pos1SV=Position1-SV_Position;// position of common vertex (btw iso track and signal muons), wrt signal vertex
      TVector3 Pos2SV=Position2-SV_Position;
      TVector3 Pos3SV=Position3-SV_Position;
      
      double Distance_Difference_Avg_var = ((Position1-Position2).Mag()+(Position2-Position3).Mag()+(Position3-Position1).Mag())/3.0;
      double Distance_SV_Avg_var = ((Pos1SV).Mag()+(Pos2SV).Mag()+(Pos3SV).Mag())/3.0;
      double Chi2_Avg_var = (Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,dR_No[i][1])+Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,dR_No[i][1])+Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,dR_No[i][1]))/3.0;
      double PosSV_Avg_1D_var = (Tau_Dir.Dot(Pos1SV+Pos2SV+Pos3SV))/Tau_Dir.Mag();
      double PosSV_Avg_1D_SVPV_var = (SVPV_Dir.Dot(Pos1SV+Pos2SV+Pos3SV))/SVPV_Dir.Mag();
      double PosSV_Avg_1D_Perp_var = (SVPV_Dir.Cross(Pos1SV+Pos2SV+Pos3SV)).Mag()/SVPV_Dir.Mag();
      
      Distance_Difference_Avg.push_back({Distance_Difference_Avg_var,dR_No[i][1]});
      
      Distance_SV_Avg.push_back({Distance_SV_Avg_var,dR_No[i][1]});
      
      Chi2_Avg.push_back({Chi2_Avg_var,dR_No[i][1]});
      
      Distance_Difference_Avg_Over_Distance_SV_Avg.push_back({Distance_Difference_Avg_var/Distance_SV_Avg_var,dR_No[i][1]});
      
      Angles3.push_back({Tau_Dir.Angle(Pos1SV),Tau_Dir.Angle(Pos2SV),Tau_Dir.Angle(Pos3SV),dR_No[i][1]});
      
      PosSV_Avg_1D_Perp.push_back({PosSV_Avg_1D_Perp_var,dR_No[i][1]});
      
      PosSV_Avg_1D_SVPV.push_back({PosSV_Avg_1D_SVPV_var,dR_No[i][1]});// dot product doesn't work for some reason, but cross product does
      
      Distance_SV_Sum.push_back({(Pos1SV+Pos2SV+Pos3SV).Mag(),dR_No[i][1]});
      
      Testing.push_back({dR_No[i][0]*Distance_SV_Avg_var*pow(Chi2_Avg_var,0.5),dR_No[i][1]});
      
      Testing1.push_back({dR_No[i][0]*Distance_SV_Avg_var*pow(Chi2_Avg_var,0.5),dR_No[i][1]});
      
      dRtoTauRefl.push_back({TrackLV.DeltaR(TauLV_Reflection),dR_No[i][1]});
      
      dRtoSVPV.push_back({TrackLV.DeltaR(SVPV_LV),dR_No[i][1]});
      
      all_vars.push_back({dR_No[i][0],dR_No[i][1],Distance_SV_Avg_var,dR_No[i][0]*Distance_SV_Avg_var*pow(Chi2_Avg_var,0.5)});
      
      //Distance_Difference_Avg and Distance_SV_Difference_Avg are the same because you're just looking at differences after coordinate transformation
      
      
    }
    
    sort( Testing1.begin(), Testing1.end() ); // sort based on first column, lowest first
    
    int OneProngCount(0);
    Matrix OneProngChi2;// n x 3 matrix where first column is chi2 of fit of four tracks, second column is the index of one isolation track, third column is the index of one isolation track
    TVector3 tau_vtx = Ntp->Vertex_Signal_KF_pos(signal_idx);
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>0;i++){
      OneProngCount+=1;
      
      std::vector<TrackParticle> TrackPair;
      TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[i][1]));
      TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
      TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
      TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_3));
      TVector3 SV_Position = Ntp->Vertex_Signal_KF_pos(signal_idx);
      TVector3 FirstGuess(0.1,0.1,0.1);
      
      //using signal SV_Position as the first guess
      Chi2VertexFitter  PairFittedVertex(TrackPair,FirstGuess);
      //PairFittedVertex.Fit();
      
      bool Whether_Matched = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[i][1]))==1?true:false;
      
      
      OneProngChi2.push_back({PairFittedVertex.ChiSquare(),dR_No[i][1],Whether_Matched,dR_No[i][0]});
      
    }
    sort( OneProngChi2.begin(), OneProngChi2.end() ); // sort based on first column, lowest first
    
    
    for(int i=0;i<OneProngChi2.size()&&MatchedIsoTrackNo.size()>0;i++){// Selecting Chi squares of 7 track combinations with lowest Chi2
      std::cout<<"Chi2 of 1-prong track with 3mu: "<< OneProngChi2[i][0] << " index1 : " << OneProngChi2[i][1] <<std::endl;
      
      if((std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), OneProngChi2[i][1]))==1){
        Rank_Correct_1iso3mu_Chi2.at(t).Fill(i,1);
      }
    }
    if(MatchedIsoTrackNo.size()>0){
      //std::cout<<" Mass of B meson is: " << LV_BMeson.M() << std::endl;
      
      std::cout<<" The required isolation track indices are: " << std::endl;
      for(int i=0;i<MatchedIsoTrackNo.size();i++){
        std::cout<<" index : " << MatchedIsoTrackNo[i] << std::endl;
      }
    }
    
    GammaVs_Tau_Energy.at(t).Fill(TauLV.E(),LV_BMeson.Gamma());
    Angle_SVPV_BMeson.at(t).Fill(fabs((SVPV_Vector).Angle(LV_BMeson.Vect())),1);
    //std::cout<<" Size of Selection_Status_TwoProng: " << Selection_Status_TwoProng.size() << std::endl;
    int TwoProngCount(-1);
    Matrix TwoProngChi2;// n x 3 matrix where first column is chi2 of fit of two tracks, second column is the index of one isolation track, third column is the index of one isolation track
    Matrix TwoProngChi2_With3mu;
    Matrix TwoProngVtxDist;
    Matrix TwoProngChi2_OtherVar;
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>1;i++){
      for(int j=0;j<dR_No.size()&&j<NoOfTracksAfterdRCut&&j<i;j++){// combinations of tracks with indices from 0 - 6
        //TwoProngCount+=1;
        
        //if (!Selection_Status_TwoProng[TwoProngCount][0]) {
        //  continue;
        //}
        
        std::vector<TrackParticle> TrackPair;
        TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[i][1]));
        TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[j][1]));
        TVector3 SV_Position = Ntp->Vertex_Signal_KF_pos(signal_idx);
        TVector3 FirstGuess(0.1,0.1,0.1);
        
        TLorentzVector TrackLV1 = Ntp->IsolationTrack_p4(signal_idx,dR_No[i][1]);
        TLorentzVector TrackLV2 = Ntp->IsolationTrack_p4(signal_idx,dR_No[j][1]);
        
        int IsoTrack_Charge1 = Ntp->IsolationTrack_charge(signal_idx,dR_No[i][1]);
        int IsoTrack_Charge2 = Ntp->IsolationTrack_charge(signal_idx,dR_No[j][1]);
        
        double TauCharge = Ntp->Muon_charge(Muon_index_1) + Ntp->Muon_charge(Muon_index_2) + Ntp->Muon_charge(Muon_index_3);
        
        //Removed Chi2 temporarily to get rid of errors
        //using signal SV_Position as the first guess
        Chi2VertexFitter  PairFittedVertex(TrackPair,FirstGuess);
        PairFittedVertex.Fit();
        
        bool Pair_Matched = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[j][1]))==2?true:false;
        bool One_Matched = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[j][1]))==1?true:false;
        
        double Pair_dR = (Ntp->IsolationTrack_p4(signal_idx,dR_No[i][1])).DeltaR(Ntp->IsolationTrack_p4(signal_idx,dR_No[j][1]));
        
        double Product1 = (Pair_dR*Testing[i][0]*Testing[j][0]*100000)/3;
        double Product2 = (Testing[i][0]*Testing[j][0]*1000)/2;
        
        
        double AngleDiff=fabs((TrackLV1.Vect()).Angle(TrackLV2.Vect()));
        
        
        TLorentzVector LVhalf1=TauLV;
        TLorentzVector LVhalf2=TrackLV1+TrackLV2;
        TLorentzVector LVmod1=TrackLV1;
        TLorentzVector LVmod2=TrackLV2;
        TLorentzVector LVCombination=LVhalf1+LVhalf2;
        LVCombination.SetE(sqrt(LVCombination.Px()*LVCombination.Px()+LVCombination.Py()*LVCombination.Py()+LVCombination.Pz()*LVCombination.Pz()+5.2795*5.2795));
        //LVhalf1.Boost(-1*LVCombination.BoostVector());
        //LVhalf2.Boost(-1*LVCombination.BoostVector());
        
        LVhalf1.Boost(-1*LV_BMeson.BoostVector());
        LVhalf2.Boost(-1*LV_BMeson.BoostVector());
        
        LVmod1.Boost(-1*LV_BMeson.BoostVector());
        LVmod2.Boost(-1*LV_BMeson.BoostVector());
        
        //double AngleCheck=fabs(LVhalf1.Vect().Angle(LVhalf2.Vect())-3.141592653589793);
        double AngleCheck=LVhalf1.Vect().Angle(LVhalf2.Vect());
        
        double AngleWithBMesonDir=LVmod1.Vect().Angle(LV_BMeson.Vect())+LVmod2.Vect().Angle(LV_BMeson.Vect());
        
        //if(fabs(LVhalf1.Vect().Angle(LVhalf2.Vect())-3.141592653589793)<0.02){
        //  TauAngleTest.at(t).Fill((MuonOS+MuonSS1+MuonSS2).M(),1);
        //}
        
        
        TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
        TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
        TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_3));
        
        Chi2VertexFitter  PairFitted_WithMus_Vertex(TrackPair,FirstGuess);
        PairFitted_WithMus_Vertex.Fit();
        
        TwoProngChi2.push_back({PairFittedVertex.ChiSquare(),dR_No[i][1],dR_No[j][1],Product1,Pair_Matched,dR_No[i][0],dR_No[j][0],Pair_dR,Product2,PairFitted_WithMus_Vertex.ChiSquare(),AngleDiff,AngleCheck,IsoTrack_Charge1,IsoTrack_Charge2,TauCharge,TrackLV1.E()+TrackLV2.E(),AngleWithBMesonDir,One_Matched,(TrackLV1+TrackLV2).M()});
        TwoProngChi2_OtherVar.push_back({AngleDiff,dR_No[i][1],dR_No[j][1],Product1,Pair_Matched,dR_No[i][0],dR_No[j][0],Pair_dR,Product2,PairFitted_WithMus_Vertex.ChiSquare(),PairFittedVertex.ChiSquare(),AngleCheck,IsoTrack_Charge1,IsoTrack_Charge2,TauCharge});
        
        TVector3 PrVtx = PairFittedVertex.GetVertex();
        TVector3 PrVtx_WithMu = PairFitted_WithMus_Vertex.GetVertex();
        
        TwoProngChi2_With3mu.push_back({PairFitted_WithMus_Vertex.ChiSquare(),dR_No[i][1],dR_No[j][1],Product1,Pair_Matched,dR_No[i][0],dR_No[j][0],Pair_dR,Product2,PairFittedVertex.ChiSquare()});
        TwoProngVtxDist.push_back({PairFittedVertex.ChiSquare(),(PrVtx-tau_vtx).Mag()});
        
        
        TVector3 Vector_A = (PrVtx-PrVtx_WithMu);
        TVector3 Vector_B = (tau_vtx-PrVtx_WithMu);
        bool Whether_Tau_Dirn_Matches_dVtx(true);
        if(fabs((Vector_A).Angle(TauLV.Vect()))>TMath::Pi()/2){
          Vector_A= -1 * Vector_A;
          Vector_B= -1 * Vector_B;
          Whether_Tau_Dirn_Matches_dVtx=false;
        }
        
        if(Pair_Matched){
          
          Angle_Tau_Vtx.at(t).Fill(fabs((tau_vtx-PrVtx_WithMu).Angle(TauLV.Vect())),1);
          
          TVector3 Vector_1 = Vector_A*(((LV_BMeson_reco-TauLV).Vect()).Mag()/(Vector_A).Mag());
          TVector3 Vector_2 = Vector_B*(((TauLV).Vect()).Mag()/(Vector_B).Mag());
          TLorentzVector PseudoLV_Iso(Vector_1,0.);
          TLorentzVector PseudoLV_Tau(Vector_2,0.);
          PseudoLV_Iso.SetE(sqrt(Vector_1.Px()*Vector_1.Px()+Vector_1.Py()*Vector_1.Py()+Vector_1.Pz()*Vector_1.Pz()+2.006*2.006));
          PseudoLV_Tau.SetE(sqrt(Vector_2.Px()*Vector_2.Px()+Vector_2.Py()*Vector_2.Py()+Vector_2.Pz()*Vector_2.Pz()+1.776*1.776));
          
          PseudoLV_Iso.Boost(-1*LV_BMeson.BoostVector());
          PseudoLV_Tau.Boost(-1*LV_BMeson.BoostVector());
          
          //Making tau point towards Z axis
          TVector3 Tau_mod1 = PseudoLV_Tau.Vect();
          TVector3 Track_mod1 = PseudoLV_Iso.Vect();
          double Phi_init1 = Tau_mod1.Phi();
          double Theta_init1 = Tau_mod1.Theta();
          if(Phi_init1 >= TMath::Pi()) Phi_init1 = Phi_init1-2*TMath::Pi();
          if(Phi_init1 <=-TMath::Pi()) Phi_init1 = Phi_init1+2*TMath::Pi();
          Track_mod1.RotateZ(-Phi_init1);
          Track_mod1.RotateY(-Theta_init1);
          
          EnergyVsPairdR_Circular_dR.at(t).Fill(Track_mod1.Phi(),Track_mod1.Theta());
          
          GammaVs_Pair_Energy.at(t).Fill((TrackLV1+TrackLV2).E(),LV_BMeson.Gamma());
          
          //check what (PrVtx-PrVtx_WithMu) matches to
          double Ang_min(99.0);
          int Min_idx(20);
          for(int m : b_meson_full_childidx){
            TLorentzVector ProductParticleLV = Ntp->MCParticle_p4(m);
            TVector3 ProductVect = ProductParticleLV.Vect();
            double Angle1 = fabs((Vector_A).Angle(ProductVect));
            //if(Angle1 >= TMath::Pi()/2) Angle1 = TMath::Pi()-Angle1;
            if(Angle1<Ang_min){
              Ang_min=Angle1;
              Min_idx=m;
            }
          }
          if(Ang_min<0.25){
            std::cout<< "For indices: "<< i <<" and "<< j <<" PDGID of MC particle is: " <<Ntp->MCParticle_pdgid(Min_idx)<<std::endl;
          }
          
          // check the angle from d(Vertex) to charged pair
          double AngleCh = fabs((Vector_A).Angle((TrackLV1+TrackLV2).Vect()));
          //if(AngleCh >= TMath::Pi()/2) AngleCh = TMath::Pi()-AngleCh;
          //if(!Whether_Tau_Dirn_Matches_dVtx){
            Angle_dVtx_Pair_Crt.at(t).Fill(AngleCh,1);
          //}
        }//if(Pair_Matched)
        if(!Pair_Matched){
          // check the angle from d(Vertex) to charged pair
          double AngleCh = fabs((Vector_A).Angle((TrackLV1+TrackLV2).Vect()));
          //if(AngleCh >= TMath::Pi()/2) AngleCh = TMath::Pi()-AngleCh;
          //if(!Whether_Tau_Dirn_Matches_dVtx){
            Angle_dVtx_Pair_InCrt.at(t).Fill(AngleCh,1);
          //}
        }
      }
    }// end i and j loops
    Matrix TwoProngChi2_AfterCuts;// n x 3 matrix where first column is chi2 of fit of two tracks, second column is the index of one isolation track, third column is the index of one isolation track
    Matrix TwoProngChi2_With3mu_AfterCuts;
    for(int i=0;i<TwoProngChi2.size();i++){
      if(TwoProngChi2[i][4]){
        TwoProngChi2_AfterCuts.push_back(TwoProngChi2[i]);
        TwoProngChi2_With3mu_AfterCuts.push_back(TwoProngChi2_With3mu[i]);
      }
    }
    
    sort( TwoProngChi2.begin(), TwoProngChi2.end() ); // sort based on first column, lowest first
    sort( TwoProngVtxDist.begin(), TwoProngVtxDist.end() ); // sort based on first column, lowest first
    sort( TwoProngChi2_With3mu.begin(), TwoProngChi2_With3mu.end() ); // sort based on first column, lowest first
    sort( TwoProngChi2_OtherVar.begin(), TwoProngChi2_OtherVar.end() ); // sort based on first column, lowest first
    
    
    if(MatchedIsoTrackNo.size()>1){
      
      std::cout<<" The dR sorted indices of the isolation tracks are: " << std::endl;
      
      for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>0;i++){
        std::cout<< "Index : " << dR_No[i][1] << " dR : " << dR_No[i][0] << " dR to refl: " << dRtoTauRefl[i][0] << " dR to SVPV: " << dRtoSVPV[i][0] << " Distance_Difference_Avg: " << Distance_Difference_Avg[i][0] << " Distance_SV_Avg: " << Distance_SV_Avg[i][0] << " Chi2_Avg: " << Chi2_Avg[i][0] << "  Angle1:" << Angles3[i][0] << "Angle2:" << Angles3[i][1] << "Angle1:" << Angles3[i][2] << " Cross:" << PosSV_Avg_1D_Perp[i][0] << " Dot:" << PosSV_Avg_1D_SVPV[i][0] << " Testing1:" << Testing1[i][0] << std::endl;
      }
      std::cout<<" The Testing1 sorted indices of the isolation tracks are: " << std::endl;
      
      for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut;i++){
        std::cout<<" Testing1 : " << Testing1[i][0] << " index : " << Testing1[i][1] << std::endl;
      }
      
      std::cout<<"Tau energy: " << TauLV.E() << std::endl;
      
    }
    
    
    for(int i=0;i<TwoProngChi2.size()&&MatchedIsoTrackNo.size()>1;i++){// Selecting Chi squares of 7 track combinations with lowest Chi2
      std::cout<<"Chi2 of 2-prong track: "<< TwoProngChi2[i][0] <<" Chi2 with 3mu: "<< TwoProngChi2[i][9] <<" avg dR: "<<TwoProngChi2[i][3]<<" avg dR2: "<<TwoProngChi2[i][8]<< " TwoProngVtxDist: " << TwoProngVtxDist[i][1] <<" pair dr: "<< TwoProngChi2[i][7]<<" pair angle: "<< TwoProngChi2[i][10]<< " Angle after boosting: "<< TwoProngChi2[i][11] << " Sum of charges: " << TwoProngChi2[i][12] + TwoProngChi2[i][13] << " Tau charge: " << TwoProngChi2[i][14] << " index1 and 2: " << TwoProngChi2[i][1] << " , " << TwoProngChi2[i][2] <<std::endl;
      
      if(TwoProngChi2[i][4]){//if 2 correct iso tracks
        Rank_Correct_2iso_Chi2.at(t).Fill(i,1);
        
        TauEnergyVsPairdR.at(t).Fill(TauLV.E(),TwoProngChi2[i][7]);
        EnergyVsPairdR.at(t).Fill(1-(TauLV.E()/LV_BMeson.E()),TwoProngChi2[i][7]);
        EnergyVsPairdR_1.at(t).Fill(TwoProngChi2[i][15]/LV_BMeson.E(),TwoProngChi2[i][7]);
        EnergyVsPairdR_2.at(t).Fill(TwoProngChi2[i][16],TwoProngChi2[i][7]);
        
        PairAngle_2Cor.at(t).Fill(TwoProngChi2[i][7],1);
        PairMass_2Cor.at(t).Fill(TwoProngChi2[i][18],1);
      }
      if(TwoProngChi2[i][17]){//if 1 correct iso track
        PairAngle_1Cor.at(t).Fill(TwoProngChi2[i][7],1);
        PairMass_1Cor.at(t).Fill(TwoProngChi2[i][18],1);
      }
      if(!TwoProngChi2[i][17]&&!TwoProngChi2[i][4]){//if no correct iso tracks
        PairAngle_0Cor.at(t).Fill(TwoProngChi2[i][7],1);
        PairMass_0Cor.at(t).Fill(TwoProngChi2[i][18],1);
      }
      if((std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), TwoProngChi2_With3mu[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), TwoProngChi2_With3mu[i][2]))==2){
        Rank_Correct_2iso3mu_Chi2.at(t).Fill(i,1);
      }
      if(TwoProngChi2_OtherVar[i][4]){
        Rank_Correct_2iso_PairAngle.at(t).Fill(i,1);
      }
    }
    
    
    if(dR_No.size()>1&&TwoProngChi2.size()>0){
      bool Lowest_Chi2_is_Correct_2iso(false);//whether the lowest Chi2 is the correct one
      if(TwoProngChi2[0][4]){
        Lowest_Chi2_is_Correct_2iso=true;
      }
      Whether_Lowest_Chi2_is_Correct_2iso.at(t).Fill(Lowest_Chi2_is_Correct_2iso,1);
    }
    
    
    // Three iso track vertex
    
    int ThreeProngCount(0);
    Matrix ThreeProngChi2;// n x 3 matrix where first column is chi2 of fit of two tracks, second column is the index of one isolation track, third column is the index of one isolation track
    Matrix ThreeProngChi2_With3mu;
    for(int i=0;i<dR_No.size()&&i<NoOfTracksAfterdRCut&&dR_No.size()>2;i++){
      for(int j=0;j<dR_No.size()&&j<NoOfTracksAfterdRCut&&j<i;j++){
        for(int k=0;k<dR_No.size()&&k<NoOfTracksAfterdRCut&&k<j;k++){// combinations of tracks with indices from 0 - 6
          ThreeProngCount+=1;
          
          std::vector<TrackParticle> TrackPair;
          TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[i][1]));
          TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[j][1]));
          TrackPair.push_back(Ntp->IsolationTrack_TrackParticle(dR_No[k][1]));
          
          TVector3 SV_Position = Ntp->Vertex_Signal_KF_pos(signal_idx);
          TVector3 FirstGuess(0.1,0.1,0.1);
          
          //using signal SV_Position as the first guess
          Chi2VertexFitter  PairFittedVertex(TrackPair,FirstGuess);
          //PairFittedVertex.Fit();
          
          TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
          TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
          TrackPair.push_back(Ntp->Muon_TrackParticle(Muon_index_3));
          
          Chi2VertexFitter  PairFitted_WithMus_Vertex(TrackPair,FirstGuess);
          //PairFitted_WithMus_Vertex.Fit();
          
          bool Whether_Matched = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[j][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), dR_No[k][1]))==3?true:false;
          
          
          ThreeProngChi2.push_back({PairFittedVertex.ChiSquare(),dR_No[i][1],dR_No[j][1],dR_No[k][1],Whether_Matched,dR_No[i][0],dR_No[j][0],dR_No[k][0]});
          ThreeProngChi2_With3mu.push_back({PairFitted_WithMus_Vertex.ChiSquare(),dR_No[i][1],dR_No[j][1],dR_No[k][1],Whether_Matched,dR_No[i][0],dR_No[j][0],dR_No[k][0]});
          
        }
      }
    }
    sort( ThreeProngChi2.begin(), ThreeProngChi2.end() ); // sort based on first column, lowest first
    sort( ThreeProngChi2_With3mu.begin(), ThreeProngChi2_With3mu.end() ); // sort based on first column, lowest first
    
    for(int i=0;i<ThreeProngChi2.size()&&MatchedIsoTrackNo.size()>2;i++){// Selecting Chi squares of 7 track combinations with lowest Chi2
      std::cout<<"Chi2 of 3-prong track: "<< ThreeProngChi2[i][0] << " index1 : " << ThreeProngChi2[i][1] << " index2 : " << ThreeProngChi2[i][2] << " index3 : " << ThreeProngChi2[i][3] <<std::endl;
      
      if((std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2[i][2])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2[i][3]))==3){
        Rank_Correct_3iso_Chi2.at(t).Fill(i,1);
      }
      if((std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2_With3mu[i][1])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2_With3mu[i][2])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), ThreeProngChi2_With3mu[i][3]))==3){
        Rank_Correct_3iso3mu_Chi2.at(t).Fill(i,1);
      }
    }
    
    
    
    Matrix IsoTrackMatchedToSV_MassMatch_Index;// contains SV index, SV track indices. Input is the correct iso tracks
    Matrix IsoTrackMatchedToSV_ThreeMassMatch_Index;
    Matrix IsoTrackMatchedToSV_Index;// contains SV index, iso track indices. Input is the 7 dR sorted tracks
    Matrix IsoTrackMatchedToSV_Combination_Index;// contains SV index, iso track indices. Input is the 2 TwoProngChi2 sorted tracks
    Matrix SV_Detailed;// for each SV in the SV collection, count the no of 3mu candidates in it and the no of correct iso tracks in it and the no. of tracks which are neither
    bool var_IsoTrackMatchedToSV_1(false);
    bool var_IsoTrackMatchedToSV_TwoMatched(false);
    bool var_IsoTrackMatchedToSV_TwoCrtIso(false);
    bool var_IsoTrackMatchedToSV_MassMatch(false);//TwoMassMatch: matching two correct iso tracks
    bool var_IsoTrackMatchedToSV_ThreeMassMatch(false);
    for(int i=0;i<Ntp->NSecondaryVertices();i++){// see if we can get matched tracks in SecondaryVertexTrack_P4 std::count
      std::vector<int> matched_track_from_SV_to_iso;// iso track index
      std::vector<int> matched_track_from_SV_to_iso_SVIndex;// SVIndex = SV track index
      std::vector<int> matched_track_from_SV_to_correct_SVIndex;
      std::vector<int> matched_track_from_SV_to_correct_IsoIndex;//to get rid of duplicates
      std::vector<int> matched_track_from_SV_to_Sorted_Combination;// iso track index
      std::vector<int> matched_track_from_SV_to_Sorted_Combination_SVIndex;
      std::vector<int> matched_track_from_SV_to_SignalMu_SVIndex;
      std::vector<int> matched_track_from_SV_to_SignalMu_SortedChargeMuIndex;//to get rid of duplicates
      for(int j=0;j<Ntp->NTracksAtSecondaryVertex(i);j++){
        TLorentzVector SV_Track_LV = Ntp->SecondaryVertexTrack_P4(i,j);
        for(int k=0;k<dR_No.size()&&k<NoOfTracksAfterdRCut;k++){// Matching dR sorted iso tracks to SV tracks
          TLorentzVector Track_LV = Ntp->IsolationTrack_p4(signal_idx,dR_No[k][1]);
          if(SV_Track_LV.DeltaR(Track_LV)<0.001&&(std::count(matched_track_from_SV_to_iso.begin(), matched_track_from_SV_to_iso.end(), dR_No[k][1]))==0){// making sure same iso track doesn't get matched to different SV tracks
            matched_track_from_SV_to_iso_SVIndex.push_back(j);
            matched_track_from_SV_to_iso.push_back(dR_No[k][1]);
          }
        }
        for(int k=0;k<MatchedIsoTrackNo.size();k++){//store correct iso tracks in this SV
          TLorentzVector Track_LV = Ntp->IsolationTrack_p4(signal_idx,MatchedIsoTrackNo[k]);
          if(SV_Track_LV.DeltaR(Track_LV)<0.001&&(std::count(matched_track_from_SV_to_correct_IsoIndex.begin(), matched_track_from_SV_to_correct_IsoIndex.end(), MatchedIsoTrackNo[k]))==0){
            matched_track_from_SV_to_correct_SVIndex.push_back(j);
            matched_track_from_SV_to_correct_IsoIndex.push_back(MatchedIsoTrackNo[k]);
          }
        }
        for(int k=0;k<3;k++){//store signal mu in this SV
          TLorentzVector Track_LV = Ntp->Muon_P4(Ntp->SortedChargeMuons(idx_vec).at(k));
          if(SV_Track_LV.DeltaR(Track_LV)<0.0005&&(std::count(matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.begin(), matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.end(), k))==0){
            matched_track_from_SV_to_SignalMu_SVIndex.push_back(j);
            matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.push_back(k);
          }
        }
        for(int k=0;k<2&&dR_No.size()>1&&TwoProngChi2.size()>0;k++){// for sorted TwoProngChi2
          TLorentzVector Track_LV = Ntp->IsolationTrack_p4(signal_idx,TwoProngChi2[0][k+1]);
          if(SV_Track_LV.DeltaR(Track_LV)<0.001&&(std::count(matched_track_from_SV_to_Sorted_Combination.begin(), matched_track_from_SV_to_Sorted_Combination.end(), TwoProngChi2[0][k+1]))==0){
            matched_track_from_SV_to_Sorted_Combination_SVIndex.push_back(j);
            matched_track_from_SV_to_Sorted_Combination.push_back(TwoProngChi2[0][k+1]);
          }
        }
      }// end of track loop
      
      
      if(matched_track_from_SV_to_correct_SVIndex.size()==2){//If two of the correct iso tracks have been reconstructed in the SV collection. Then do an invariant mass check to verify we've got the right tracks
        TLorentzVector LV3 = Ntp->SecondaryVertexTrack_P4(i,matched_track_from_SV_to_correct_SVIndex[0]);
        TLorentzVector LV4 = Ntp->SecondaryVertexTrack_P4(i,matched_track_from_SV_to_correct_SVIndex[1]);
        
        
        if(MatchedIsoTrackNo.size()>=2){//should be MatchedIsoTrackNo.size()>=2
          var_IsoTrackMatchedToSV_TwoCrtIso=true;
          
          
          TLorentzVector LV1 = Ntp->IsolationTrack_p4(signal_idx,matched_track_from_SV_to_correct_IsoIndex[0]);
          TLorentzVector LV2 = Ntp->IsolationTrack_p4(signal_idx,matched_track_from_SV_to_correct_IsoIndex[1]);
          
          //std::cout<<"Before Invariant mass (crt): "<<(LV1+LV2).M()<<" Invariant mass (from SV): "<<(LV3+LV4).M()<<std::endl;
          if(abs((LV1+LV2).M()-(LV3+LV4).M())<0.01&&IsoTrackMatchedToSV_MassMatch_Index.size()==0){
            //::cout<<"Invariant mass (crt): "<<(LV1+LV2).M()<<" Invariant mass (from SV): "<<(LV3+LV4).M()<<std::endl;
            var_IsoTrackMatchedToSV_MassMatch=true;
            IsoTrackMatchedToSV_MassMatch_Index.push_back({i,matched_track_from_SV_to_correct_SVIndex[0],matched_track_from_SV_to_correct_IsoIndex[0]});
            IsoTrackMatchedToSV_MassMatch_Index.push_back({i,matched_track_from_SV_to_correct_SVIndex[1],matched_track_from_SV_to_correct_IsoIndex[1]});
          }
          IsoTrackMatchedToSV_MassMatch.at(t).Fill(var_IsoTrackMatchedToSV_MassMatch,1);
        }
        IsoTrackMatchedToSV_TwoCrtIso.at(t).Fill(var_IsoTrackMatchedToSV_TwoCrtIso,1);
        
        
      }// if(matched_track_from_SV_to_correct_SVIndex.size()==2)
      
      if(matched_track_from_SV_to_correct_SVIndex.size()==3){//if three of the correct iso tracks have been reconstructed in the SV collection
        TLorentzVector LV4 = Ntp->SecondaryVertexTrack_P4(i,matched_track_from_SV_to_correct_SVIndex[0]);
        TLorentzVector LV5 = Ntp->SecondaryVertexTrack_P4(i,matched_track_from_SV_to_correct_SVIndex[1]);
        TLorentzVector LV6 = Ntp->SecondaryVertexTrack_P4(i,matched_track_from_SV_to_correct_SVIndex[2]);
        
        if(MatchedIsoTrackNo.size()>=3){//should be MatchedIsoTrackNo.size()>=3
          
          TLorentzVector LV1 = Ntp->IsolationTrack_p4(signal_idx,matched_track_from_SV_to_correct_IsoIndex[0]);
          TLorentzVector LV2 = Ntp->IsolationTrack_p4(signal_idx,matched_track_from_SV_to_correct_IsoIndex[1]);
          TLorentzVector LV3 = Ntp->IsolationTrack_p4(signal_idx,matched_track_from_SV_to_correct_IsoIndex[2]);
          
          //std::cout<<"Before Invariant mass-3 (crt): "<<(LV1+LV2+LV3).M()<<" Invariant mass-3 (from SV): "<<(LV4+LV5+LV6).M()<<std::endl;
          
          if(abs((LV1+LV2+LV3).M()-(LV4+LV5+LV6).M())<0.01&&IsoTrackMatchedToSV_ThreeMassMatch_Index.size()==0){
            //std::cout<<"Invariant mass-3 (crt): "<<(LV1+LV2+LV3).M()<<" Invariant mass-3 (from SV): "<<(LV4+LV5+LV6).M()<<std::endl;
            var_IsoTrackMatchedToSV_ThreeMassMatch=true;
            IsoTrackMatchedToSV_ThreeMassMatch_Index.push_back({i,matched_track_from_SV_to_correct_SVIndex[0],matched_track_from_SV_to_correct_IsoIndex[0]});
            IsoTrackMatchedToSV_ThreeMassMatch_Index.push_back({i,matched_track_from_SV_to_correct_SVIndex[1],matched_track_from_SV_to_correct_IsoIndex[1]});
            IsoTrackMatchedToSV_ThreeMassMatch_Index.push_back({i,matched_track_from_SV_to_correct_SVIndex[2],matched_track_from_SV_to_correct_IsoIndex[2]});
          }
          
          IsoTrackMatchedToSV_ThreeMassMatch.at(t).Fill(var_IsoTrackMatchedToSV_ThreeMassMatch,1);
        }
      }
      
      //Here, for each SV in the SV collection, count the no of 3mu candidates in it and the no of correct iso tracks in it and the no. of tracks which are neither. In order to account for mass matches, do (no of correct iso tracks in it) in two columns where in one column we have 2, 3 or -1 (if mass doesn't match or if the number is different from 2 or 3). In the other column, put matched_track_from_SV_to_correct_SVIndex.size(). Also store the sum of all p4 so that we can do dR from SV to tau.
      if(var_IsoTrackMatchedToSV_MassMatch||var_IsoTrackMatchedToSV_ThreeMassMatch){
        SV_Detailed.push_back({i,matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.size(),matched_track_from_SV_to_correct_SVIndex.size(),matched_track_from_SV_to_correct_SVIndex.size(),(Ntp->NTracksAtSecondaryVertex(i)-matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.size()-matched_track_from_SV_to_correct_SVIndex.size()),matched_track_from_SV_to_iso_SVIndex.size()});
      }
      else{
        SV_Detailed.push_back({i,matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.size(),-1,matched_track_from_SV_to_correct_SVIndex.size(),(Ntp->NTracksAtSecondaryVertex(i)-matched_track_from_SV_to_SignalMu_SortedChargeMuIndex.size()-matched_track_from_SV_to_correct_SVIndex.size()),matched_track_from_SV_to_iso_SVIndex.size()});
      }
      
      
      if(matched_track_from_SV_to_iso.size()>1){
        var_IsoTrackMatchedToSV_1=true;
        //std::cout<<"SV No. "<<i+1<<" has matched to "<<matched_track_from_SV_to_iso.size()<<" iso tracks of which "<<matched_track_from_SV_to_correct_SVIndex.size()<<" are correct. Iso tracks: "<<dR_No.size()<<" and correct tracks: "<<MatchedIsoTrackNo.size()<<std::endl;
        
        if(matched_track_from_SV_to_correct_SVIndex.size()==2){
          var_IsoTrackMatchedToSV_TwoMatched=true;
        }
        IsoTrackMatchedToSV_TwoMatched.at(t).Fill(var_IsoTrackMatchedToSV_TwoMatched,1);
        
        if(IsoTrackMatchedToSV_Index.size()==0){
          for(int j=0;j<matched_track_from_SV_to_iso.size();j++){
            IsoTrackMatchedToSV_Index.push_back({i,matched_track_from_SV_to_iso_SVIndex[j],matched_track_from_SV_to_iso[j]});//note: this has both SV and iso track indices
          }
        }
      }
      
      if(IsoTrackMatchedToSV_Combination_Index.size()==0&&matched_track_from_SV_to_Sorted_Combination_SVIndex.size()>1){// for matching with sorted combinations
        for(int j=0;j<matched_track_from_SV_to_Sorted_Combination_SVIndex.size()&&TwoProngChi2.size()>0;j++){
          IsoTrackMatchedToSV_Combination_Index.push_back({i,matched_track_from_SV_to_Sorted_Combination_SVIndex[j],matched_track_from_SV_to_Sorted_Combination[j],TwoProngChi2[0][3]});//note: this has both SV and iso track indices
        }
      }
      IsoTrackMatchedToSV_1.at(t).Fill(var_IsoTrackMatchedToSV_1,1);
    }// i
    
    // SV: signal SV, iSV: secondary vertex from 2 prong, isv: all the other SVs from the collection
    
    if(var_IsoTrackMatchedToSV_MassMatch){
      TVector3 iSVSV_Vector = Ntp->SecondaryVertexPosition(IsoTrackMatchedToSV_MassMatch_Index[0][0]) - Ntp->Vertex_Signal_KF_pos(signal_idx);// IsoTrackMatchedToSV_MassMatch_Index[0][0] = i, index of iSV
      Angle_SVPV_iSVSV.at(t).Fill(SVPV_Vector.Angle(iSVSV_Vector),1);
      iSVSV_Distance_Sig.at(t).Fill(Ntp->FlightLength_significance(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx),Ntp->SecondaryVertexPosition(IsoTrackMatchedToSV_MassMatch_Index[0][0]),Ntp->SecondaryVertexCovariance(IsoTrackMatchedToSV_MassMatch_Index[0][0])),1);
      iSVSV_Distance.at(t).Fill(iSVSV_Vector.Mag(),1);
    }
    
    if(!var_IsoTrackMatchedToSV_MassMatch){
      for(int i=0;i<Ntp->NSecondaryVertices();i++){
        TVector3 isvSV_Vector = Ntp->SecondaryVertexPosition(i) - Ntp->Vertex_Signal_KF_pos(signal_idx);// IsoTrackMatchedToSV_MassMatch_Index[0][0] = i, index of iSV
        Angle_SVPV_isvSV.at(t).Fill(SVPV_Vector.Angle(isvSV_Vector),1);
        if(!(isvSV_Vector.Mag()<0.00001)){
          isvSV_Distance_Sig.at(t).Fill(Ntp->FlightLength_significance(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx),Ntp->SecondaryVertexPosition(i),Ntp->SecondaryVertexCovariance(i)),1);
          isvSV_Distance.at(t).Fill(isvSV_Vector.Mag(),1);
        }
      }
    }
    
    /*
    for(int i=0;i<SV_Detailed.size();i++){
      if(SV_Detailed[i][3]>0){
        SVCollectionNoOfSignalMu.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][3]>1){
        SVCollectionNoOfSignalMu_if1.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither_if1.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][3]==1){
        SVCollectionNoOfSignalMu_ifmore1.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither_ifmore1.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][1]>0&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]){//&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]: to make sure if we have signal muons, then there should be other tracks as well
        SVCollectionNoOfCrt.at(t).Fill(SV_Detailed[i][3],1);
      }
      if(SV_Detailed[i][1]==3&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]){
        SVCollectionNoOfCrt_if3.at(t).Fill(SV_Detailed[i][3],1);
      }
    }
    */
    
    for(int i=0;i<SV_Detailed.size();i++){
      if(SV_Detailed[i][3]>0&&SV_Detailed[i][5]==(Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])-SV_Detailed[i][1])){
        SVCollectionNoOfSignalMu.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][3]>1&&SV_Detailed[i][5]==(Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])-SV_Detailed[i][1])){
        SVCollectionNoOfSignalMu_if1.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither_if1.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][3]==1&&SV_Detailed[i][5]==(Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])-SV_Detailed[i][1])){
        SVCollectionNoOfSignalMu_ifmore1.at(t).Fill(SV_Detailed[i][1],1);
        SVCollectionNoOfNeither_ifmore1.at(t).Fill(SV_Detailed[i][4],1);
      }
      if(SV_Detailed[i][1]>0&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]&&SV_Detailed[i][5]==(Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])-SV_Detailed[i][1])){//&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]: to make sure if we have signal muons, then there should be other tracks as well
        SVCollectionNoOfCrt.at(t).Fill(SV_Detailed[i][3],1);
      }
      if(SV_Detailed[i][1]==3&&Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])>SV_Detailed[i][1]&&SV_Detailed[i][5]==(Ntp->NTracksAtSecondaryVertex(SV_Detailed[i][0])-SV_Detailed[i][1])){
        SVCollectionNoOfCrt_if3.at(t).Fill(SV_Detailed[i][3],1);
      }
    }
    
    
    if(IsoTrackMatchedToSV_Index.size()>1){
      IsoTrackMatchedToSV_Count.at(t).Fill(IsoTrackMatchedToSV_Index.size(),1);
      /*
      for(int i=0;i<IsoTrackMatchedToSV_Index.size();i++){
        for(int j=0;j<IsoTrackMatchedToSV_Index.size();j++){
          //try to get all combinations of tracks and sort them by dR
        }
      }
      */
    }
    
    /*
    if(MatchedIsoTrackNo.size()==2){
      
      TLorentzVector LV1 = Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_MassMatch_Index[0][0],IsoTrackMatchedToSV_MassMatch_Index[0][1]);
      TLorentzVector LV2 = Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_MassMatch_Index[1][0],IsoTrackMatchedToSV_MassMatch_Index[1][1]);
      
      if(IsoTrackMatchedToSV_MassMatch_Index.size()==2&&IsoTrackMatchedToSV_Combination_Index.size()==2){
        std::cout<<"MatchedIsoTrackNo: "<<MatchedIsoTrackNo[0]<<","<<MatchedIsoTrackNo[1]<<" IsoTrackMatchedToSV_MassMatch_Index: "<<IsoTrackMatchedToSV_MassMatch_Index[0][2]<<","<<IsoTrackMatchedToSV_MassMatch_Index[1][2]<<" IsoTrackMatchedToSV_Combination_Index: "<<IsoTrackMatchedToSV_Combination_Index[0][2]<<","<<IsoTrackMatchedToSV_Combination_Index[1][2]<<std::endl;
      }
      
    }
    */
    
    if(MatchedIsoTrackNo.size()>=2){
      IsoTrackMatchedToSV_MassMatch1.at(t).Fill(var_IsoTrackMatchedToSV_MassMatch,1);// Efficiency check
      
      
      if(var_IsoTrackMatchedToSV_MassMatch){
        TLorentzVector LV_2Prong = Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_MassMatch_Index[0][0],IsoTrackMatchedToSV_MassMatch_Index[0][1])+Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_MassMatch_Index[0][0],IsoTrackMatchedToSV_MassMatch_Index[1][1]);
        float LV_2ProngFullSV = Ntp->SecondaryVertexMass(IsoTrackMatchedToSV_MassMatch_Index[0][0]);
        InvMass2ProngMatched.at(t).Fill(LV_2Prong.M(),1);
        InvMass2ProngMatchedSV.at(t).Fill(LV_2ProngFullSV,1);
        SVNoOfTracksMatched.at(t).Fill(Ntp->NTracksAtSecondaryVertex(IsoTrackMatchedToSV_MassMatch_Index[0][0]),1);
        
        if(Ntp->NTracksAtSecondaryVertex(IsoTrackMatchedToSV_MassMatch_Index[0][0])==5){
          InvMassTotal1.at(t).Fill(Ntp->SecondaryVertexMass(IsoTrackMatchedToSV_MassMatch_Index[0][0]),1);
          
          //Here, I'm checking to see if the 5 bin is 3 mu + 2 correct tracks
          TLorentzVector Possible_Tau_LV;
          for(int i=0;i<Ntp->NTracksAtSecondaryVertex(IsoTrackMatchedToSV_MassMatch_Index[0][0]);i++){
            Possible_Tau_LV = Possible_Tau_LV + Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_MassMatch_Index[0][0],i);
          }
          InvMassTotal2.at(t).Fill((Possible_Tau_LV-LV_2Prong).M(),1);
        }
        
        InvMassTotal.at(t).Fill((LV_2Prong+TauLV).M(),1);
      }
      if(!var_IsoTrackMatchedToSV_MassMatch&&!var_IsoTrackMatchedToSV_ThreeMassMatch){
        for(int i=0;i<Ntp->NSecondaryVertices();i++){
          InvMass2ProngNotMatched.at(t).Fill(Ntp->SecondaryVertexMass(i),1);
          SVNoOfTracksUnmatched.at(t).Fill(Ntp->NTracksAtSecondaryVertex(i),1);
        }
        
      }
    }
    
    SVSize.at(t).Fill(Ntp->NSecondaryVertices(),1);
    
    if(MatchedIsoTrackNo.size()>=3){
      IsoTrackMatchedToSV_ThreeMassMatch1.at(t).Fill(var_IsoTrackMatchedToSV_ThreeMassMatch,1);// Efficiency check
      
      
      if(var_IsoTrackMatchedToSV_ThreeMassMatch){
        TLorentzVector LV_3Prong = Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_ThreeMassMatch_Index[0][0],IsoTrackMatchedToSV_ThreeMassMatch_Index[0][1])+Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_ThreeMassMatch_Index[0][0],IsoTrackMatchedToSV_ThreeMassMatch_Index[1][1])+Ntp->SecondaryVertexTrack_P4(IsoTrackMatchedToSV_ThreeMassMatch_Index[0][0],IsoTrackMatchedToSV_ThreeMassMatch_Index[2][1]);
        float LV_3ProngFullSV = Ntp->SecondaryVertexMass(IsoTrackMatchedToSV_ThreeMassMatch_Index[0][0]);
        InvMass3ProngMatched.at(t).Fill(LV_3Prong.M(),1);
        InvMass3ProngMatchedSV.at(t).Fill(LV_3ProngFullSV,1);
        SVNoOfTracksMatchedThree.at(t).Fill(Ntp->NTracksAtSecondaryVertex(IsoTrackMatchedToSV_ThreeMassMatch_Index[0][0]),1);
      }
      if(!var_IsoTrackMatchedToSV_MassMatch&&!var_IsoTrackMatchedToSV_ThreeMassMatch){
        for(int i=0;i<Ntp->NSecondaryVertices();i++){
          InvMass3ProngNotMatched.at(t).Fill(Ntp->SecondaryVertexMass(i),1);
        }
      }
      
    }
    
    
    if(IsoTrackMatchedToSV_Combination_Index.size()==2){// results of matching iso tracks sorted by combination var to SV tracks
      //std::cout<<"IsoTrackMatchedToSV_Combination_Index: "<<IsoTrackMatchedToSV_Combination_Index[0][2]<<","<<IsoTrackMatchedToSV_Combination_Index[1][2]<<std::endl;
      if(IsoTrackMatchedToSV_MassMatch_Index.size()==2&&MatchedIsoTrackNo.size()==2){
        //std::cout<<"MatchedIsoTrackNo: "<<MatchedIsoTrackNo[0]<<","<<MatchedIsoTrackNo[1]<<" IsoTrackMatchedToSV_MassMatch_Index: "<<IsoTrackMatchedToSV_MassMatch_Index[0][2]<<","<<IsoTrackMatchedToSV_MassMatch_Index[1][2]<<std::endl;
      }
      bool Whether_Comb_Correct = (std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), IsoTrackMatchedToSV_Combination_Index[0][2])+std::count(MatchedIsoTrackNo.begin(), MatchedIsoTrackNo.end(), IsoTrackMatchedToSV_Combination_Index[1][2]))==2;
      IsoTrackMatchedToSV_CombMatch.at(t).Fill(Whether_Comb_Correct,1);
      if(Whether_Comb_Correct){
        CombMatch_Avg1.at(t).Fill(IsoTrackMatchedToSV_Combination_Index[0][3],1);
      }
      if(!Whether_Comb_Correct){
        CombMatch_Avg2.at(t).Fill(IsoTrackMatchedToSV_Combination_Index[0][3],1);
      }
    }
    
    
    
    
    
    NumberOfFS_ChargedParticles.at(t).Fill(b_meson_full_childidx_FS.size(),1);// Number of charged particles that can leave tracks
    NumberOfFS_ChargedParticles_RecoMatch.at(t).Fill(b_meson_full_childidx_FS_RecoMC.size(),1);// Number of charged particles reconstructed
    
    if(b_meson_full_childidx_FS.size()==1){
      NumberOfRecoChargedParticlesIfMC1.at(t).Fill(b_meson_full_childidx_FS_RecoMC.size(),1);
    }
    
    if(b_meson_full_childidx_FS.size()==2){
      NumberOfRecoChargedParticlesIfMC2.at(t).Fill(b_meson_full_childidx_FS_RecoMC.size(),1);
    }
    
    if(b_meson_full_childidx_FS.size()==3){
      NumberOfRecoChargedParticlesIfMC3.at(t).Fill(b_meson_full_childidx_FS_RecoMC.size(),1);
    }
    
    // Two-Prong Events
    
    
    if(b_meson_full_childidx_FS_RecoMC.size()==2){
      
      for(int i=0;i<Ntp->NIsolationTrack(signal_idx);i++){
      
      }
    
     TrackToTauDr2Prong.at(t).Fill(dR_To_Tau_Prong.at(0),1);
     TrackToTauDr2Prong.at(t).Fill(dR_To_Tau_Prong.at(1),1);
     
     NoOfIsoTracks2Prong.at(t).Fill(Ntp->NIsolationTrack(signal_idx),1);
     
     TLorentzVector L1_Track(p1_reco.at(0),p2_reco.at(0),p3_reco.at(0),es_reco.at(0));//Track
     TLorentzVector L2_Track(p1_reco.at(1),p2_reco.at(1),p3_reco.at(1),es_reco.at(1));//Track
     
     TwoProngInvariantMassReco.at(t).Fill((L1_Track+L2_Track).M(),1);
     
     TLorentzVector L1_MC(p1.at(0),p2.at(0),p3.at(0),es.at(0));//MC
     TLorentzVector L2_MC(p1.at(1),p2.at(1),p3.at(1),es.at(1));//MC
     
     TwoProngInvariantMassMC.at(t).Fill((L1_MC+L2_MC).M(),1);
     
     InvMass2_vs_pdgid.at(t).Fill((L1_MC+L2_MC).M(),mc_pdgid.at(0),1);
     InvMass2_vs_pdgid.at(t).Fill((L1_MC+L2_MC).M(),mc_pdgid.at(1),1);
     
     if(pts_reco.at(0)>pts_reco.at(1)){
       TwoProngTrackPt.at(t).Fill(pts_reco.at(0),1);
       TwoProngTrack2Pt.at(t).Fill(pts_reco.at(1),1);
       
       TwoProngTrackEta.at(t).Fill(etas_reco.at(0),1);
       TwoProngTrack2Eta.at(t).Fill(etas_reco.at(1),1);
       
     }
     else{
       TwoProngTrackPt.at(t).Fill(pts_reco.at(1),1);
       TwoProngTrack2Pt.at(t).Fill(pts_reco.at(0),1);
       
       TwoProngTrackEta.at(t).Fill(etas_reco.at(1),1);
       TwoProngTrack2Eta.at(t).Fill(etas_reco.at(0),1);
     }
     
     dRmin_sum_vs_InvariantMass_2prong.at(t).Fill((dR_min_reco.at(0)+dR_min_reco.at(1)),(L1_Track+L2_Track).M(),1);
     /*
     if((L1_MC+L2_MC).M()<0.6){
       if(id ==60 ||  id ==90){// or id == 40){
        std::cout<<"--------------  Low mass ----------------"<< (L1_MC+L2_MC).M() <<std::endl;
        std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Ntp->printMCDecayChainOfEvent(true, true, true, true);
        std::cout<< "\n\n\n\n\n\n";
       }
     }
     */
     
    }// if(b_meson_full_childidx_FS_RecoMC.size()==2)
    //---------------------- Fit signal vertex
    
    /*
    TVector3 vguess(0.,0.,0.);
    std::vector<TrackParticle> MuonsTrackParticles;
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_3));

    Chi2VertexFitter  Fitter(MuonsTrackParticles,vguess);

    if(Fitter.Fit())std::cout<<" 3mu tau vertex fit  "<< Fitter.ChiSquare() << std::endl;
    
    */


    //---------------------- Fit signal vertex
    
    
    
    sort( pT_No_7.rbegin(), pT_No_7.rend() ); // sort based on first column, highest first. Moving this here due to errors when running with vertex finder
    
    if(MatchedIsoTrackNo.size()>0){
      //std::cout<<" The pT sorted indices of the 7 isolation tracks are: " << std::endl;
      
      for(int i=0;i<pT_No_7.size();i++){
        //std::cout<<" pT : " << pT_No_7[i][0] << " index : " << pT_No_7[i][1] << std::endl;
        
        for(int j=0;j<MatchedIsoTrackNo.size();j++){// to get the rank
          if(pT_No_7[i][1]==MatchedIsoTrackNo[j]){
            RankMatchedTrackpT.at(t).Fill(i,1);
          }
        }
      }// i<pT_No_7.size()
      
      std::cout<<" The required isolation track indices are: " << std::endl;
      
      for(int i=0;i<MatchedIsoTrackNo.size();i++){
        std::cout<<" index : " << MatchedIsoTrackNo[i] << std::endl;
      }
      
      std::cout<<"\n\n\n" << std::endl;
      std::cout<<"No of secondary vertices is: " <<Ntp->NSecondaryVertices()<< std::endl;
      std::cout<<"Size of TwoProngChi2 is: "<< TwoProngChi2.size() << std::endl;
      std::cout<<"\n\n\n" << std::endl;
      std::cout<<"End of event with MatchedIsoTrackNo size > 0." << std::endl;
      std::cout<<"\n\n\n" << std::endl;
      
    }// end if MatchedIsoTrackNo.size()>0
    
    sort( Testing1.begin(), Testing1.end() ); // sort based on first column, lowest first
    sort( Distance_Difference_Avg.begin(), Distance_Difference_Avg.end() ); // sort based on first column, lowest first
    
    for(int i=0;i<Testing1.size();i++){ // Filling some track ranks, ranked according to certain variables. Also some comparisons
      
      bool Track_Combn_NotReco(true);
      bool Track_AvgDiff_NotReco(true);
      bool Track_dR_NotReco(true);
      bool Track_Avg_NotReco(true);
      
      for(int j=0;j<MatchedIsoTrackNo.size();j++){// to get the rank
        if(Testing1[i][1]==MatchedIsoTrackNo[j]){
          Track_Combn_NotReco=false;
          RankMatchedTrackCombn.at(t).Fill(i,1);
          var_Correct_Iso_Combn.at(t).Fill(Testing1[i][0]);
        }
        if(Distance_Difference_Avg[i][1]==MatchedIsoTrackNo[j]){
          Track_AvgDiff_NotReco=false;
          RankMatchedTrackAvgDiff.at(t).Fill(i,1);
          var_Correct_Iso_AvgDiff.at(t).Fill(Distance_Difference_Avg[i][0]);
        }
        if(dR_No[i][1]==MatchedIsoTrackNo[j]){
          Track_dR_NotReco=false;
          var_Correct_Iso_dR.at(t).Fill(dR_No[i][0]);
        }
        if(Distance_SV_Avg[i][1]==MatchedIsoTrackNo[j]){
          Track_Avg_NotReco=false;
          var_Correct_Iso_Avg.at(t).Fill(Distance_SV_Avg[i][0]);
        }
      }
      
      if(Track_Combn_NotReco){
        var_All_7_Iso_Combn.at(t).Fill(Testing1[i][0]);
      }
      if(Track_AvgDiff_NotReco){
        var_All_7_Iso_AvgDiff.at(t).Fill(Distance_Difference_Avg[i][0]);
      }
      if(Track_dR_NotReco){
        var_All_7_Iso_dR.at(t).Fill(dR_No[i][0]);
      }
      if(Track_Avg_NotReco){
        var_All_7_Iso_Avg.at(t).Fill(Distance_SV_Avg[i][0]);
      }
      
    }
    
    for(int i=0;i<TwoProngChi2.size();i++){ // Filling some track pair ranks, ranked according to certain variables. Also some comparisons
      
      if(TwoProngChi2[i][4]){
        RankMatchedTrackPairdR.at(t).Fill(i,1);
        TrackPairdR_Crt.at(t).Fill(TwoProngChi2[i][3]);
      }
      
      if(!TwoProngChi2[i][4]){
        TrackPairdR_Bkg.at(t).Fill(TwoProngChi2[i][8]);
      }
      
    }
    
    
    
    // Three prong
    
    if(b_meson_full_childidx_FS_RecoMC.size()==3){
     
     TrackToTauDr3Prong.at(t).Fill(dR_To_Tau_Prong.at(0),1);
     TrackToTauDr3Prong.at(t).Fill(dR_To_Tau_Prong.at(1),1);
     TrackToTauDr3Prong.at(t).Fill(dR_To_Tau_Prong.at(2),1);
     
     std::vector<float> vect = pts_reco;
     std::sort (vect.begin(), vect.end());
     
     ThreeProngTrackPt.at(t).Fill(vect.at(2),1);
     ThreeProngTrack2Pt.at(t).Fill(vect.at(1),1);
     ThreeProngTrack3Pt.at(t).Fill(vect.at(0),1);
     
     
     NoOfIsoTracks3Prong.at(t).Fill(Ntp->NIsolationTrack(signal_idx),1);
     
     TLorentzVector L1_Track(p1_reco.at(0),p2_reco.at(0),p3_reco.at(0),es_reco.at(0));//Track
     TLorentzVector L2_Track(p1_reco.at(1),p2_reco.at(1),p3_reco.at(1),es_reco.at(1));//Track
     TLorentzVector L3_Track(p1_reco.at(2),p2_reco.at(2),p3_reco.at(2),es_reco.at(2));//Track
     
     ThreeProngInvariantMassReco.at(t).Fill((L1_Track+L2_Track+L3_Track).M(),1);
     
     TLorentzVector L1_MC(p1.at(0),p2.at(0),p3.at(0),es.at(0));//MC
     TLorentzVector L2_MC(p1.at(1),p2.at(1),p3.at(1),es.at(1));//MC
     TLorentzVector L3_MC(p1.at(2),p2.at(2),p3.at(2),es.at(2));//MC
     
     ThreeProngInvariantMassMC.at(t).Fill((L1_MC+L2_MC+L3_MC).M(),1);
     
     InvMass3_vs_pdgid.at(t).Fill((L1_MC+L2_MC+L3_MC).M(),mc_pdgid.at(0),1);
     InvMass3_vs_pdgid.at(t).Fill((L1_MC+L2_MC+L3_MC).M(),mc_pdgid.at(1),1);
     InvMass3_vs_pdgid.at(t).Fill((L1_MC+L2_MC+L3_MC).M(),mc_pdgid.at(2),1);
     
     dRmin_sum_vs_InvariantMass_2prong.at(t).Fill((dR_min_reco.at(0)+dR_min_reco.at(1)+dR_min_reco.at(2)),(L1_Track+L2_Track+L3_Track).M(),1);
     
     
    }//b_meson_full_childidx_FS_RecoMC.size()==3
    
    
    /*
    
    
    if(MatchedIsoTrackNo.size()>1){
      std::cout<<" The required isolation track indices are: " << std::endl;
      for(int i=0;i<MatchedIsoTrackNo.size();i++){
        std::cout<<" index : " << MatchedIsoTrackNo[i] << std::endl;
      }
    }// end if MatchedIsoTrackNo.size()>0
    
    */
    
    
    
    // Event Printouts
    /*
    if(b_meson_full_childidx_FS.size()==2){
      if(id ==60 ||  id ==90){// or id == 40){
        std::cout<<"--------------  Two Prong MC ----------------" << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Ntp->printMCDecayChainOfEvent(true, true, true, true);
        std::cout<< "\n\n\n\n\n\n";
      }
    }
    if(b_meson_full_childidx_FS.size()==3){
      if(id ==60 ||  id ==90){// or id == 40){
        std::cout<<"--------------  Three Prong MC ----------------" << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        std::cout<<" idx = "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
        Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
        Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
        Ntp->printMCDecayChainOfEvent(true, true, true, true);
        std::cout<< "\n\n\n\n\n\n";
      }
    }
    */
    
    
    
    //Doing the same with dR<0.005
    
    if(b_meson_full_childidx_FS_RecoMC005.size()==2){
    
     TLorentzVector L1_Track(p1_reco005.at(0),p2_reco005.at(0),p3_reco005.at(0),es_reco005.at(0));//Track
     TLorentzVector L2_Track(p1_reco005.at(1),p2_reco005.at(1),p3_reco005.at(1),es_reco005.at(1));//Track
     
     TwoProngInvariantMassReco005.at(t).Fill((L1_Track+L2_Track).M(),1);
     
     TLorentzVector L1_MC(p1005.at(0),p2005.at(0),p3005.at(0),es005.at(0));//MC
     TLorentzVector L2_MC(p1005.at(1),p2005.at(1),p3005.at(1),es005.at(1));//MC
     
     TwoProngInvariantMassMC005.at(t).Fill((L1_MC+L2_MC).M(),1);
    }
    
    if(b_meson_full_childidx_FS_RecoMC005.size()==3){
     TLorentzVector L1_Track(p1_reco005.at(0),p2_reco005.at(0),p3_reco005.at(0),es_reco005.at(0));//Track
     TLorentzVector L2_Track(p1_reco005.at(1),p2_reco005.at(1),p3_reco005.at(1),es_reco005.at(1));//Track
     TLorentzVector L3_Track(p1_reco005.at(2),p2_reco005.at(2),p3_reco005.at(2),es_reco005.at(2));//Track
     
     ThreeProngInvariantMassReco005.at(t).Fill((L1_Track+L2_Track+L3_Track).M(),1);
     
     TLorentzVector L1_MC(p1005.at(0),p2005.at(0),p3005.at(0),es005.at(0));//MC
     TLorentzVector L2_MC(p1005.at(1),p2005.at(1),p3005.at(1),es005.at(1));//MC
     TLorentzVector L3_MC(p1005.at(2),p2005.at(2),p3005.at(2),es005.at(2));//MC
     
     ThreeProngInvariantMassMC005.at(t).Fill((L1_MC+L2_MC+L3_MC).M(),1);
    }
    
    
    
    // Apply cuts to all 7 iso tracks close to tau
    // all_vars variables:  dR_No[i][0],dR_No[i][1],Distance_SV_Avg_var,dR_No[i][0]*Distance_SV_Avg_var*pow(Chi2_Avg_var,0.5)
    
    Matrix dR_No_cut;
    std::vector<float> MatchedIsoTrackNo_cut; // Just the index, 'j' of the matched isolation tracks after cut
    
    for(int i=0;i<all_vars.size()&&i<NoOfTracksAfterdRCut;i++){
      if(all_vars[i][0]<0.7&&all_vars[i][2]<1.0&&all_vars[i][3]<0.5){
        dR_No_cut.push_back({all_vars[i][0],all_vars[i][1]});
        
        for(int j=0;j<MatchedIsoTrackNo.size();j++){// to get the rank
          if(all_vars[i][1]==MatchedIsoTrackNo[j]){
            MatchedIsoTrackNo_cut.push_back(MatchedIsoTrackNo[j]);
          }
        }
        
      }
    }
    
    for(int i=0;i<dR_No_cut.size();i++){
      //std::cout<<" dR : " << dR_No[i][0] << " index : " << dR_No[i][1] << std::endl;
      
      for(int j=0;j<MatchedIsoTrackNo_cut.size();j++){// to get the rank
        if(dR_No_cut[i][1]==MatchedIsoTrackNo_cut[j]){
          RankMatchedTrackdR_cut.at(t).Fill(i,1);
        }
      }
      
    }
    
    // Two-Prong after cuts
    
    if(MatchedIsoTrackNo_cut.size()==2){
      
      
    }// if(MatchedIsoTrackNo_cut.size()==2)
    
    
    
  }// end of if status
}


void  SignalVertexSelector::Finish(){

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





