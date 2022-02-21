#include "VertexFits.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


//using namespace std;

VertexFits::VertexFits(TString Name_, TString id_):
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
  mvaA1_(0.117), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2_(0.221),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1_(0.134),
  mvaB2_(0.223),
  mvaC1_(0.143),
  mvaC2_(0.227)
{
  // This is a class constructor;
}



VertexFits::~VertexFits(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  VertexFits::Configure(){


  //  This mini tree is for limit extraction
  gErrorIgnoreLevel = kFatal;





  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  //*** defined the bdt reader for event selection; readerA- category A, readerB - category B ...
  readerA = new TMVA::Reader( "!Color:!Silent" );
  readerA->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2);
  readerA->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerA->AddVariable( "var_flightLenSig", &var_flightLenSig);
  readerA->AddVariable( "var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerA->AddVariable( "var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerA->AddVariable( "var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);
  readerA->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerA->AddVariable( "var_maxMuonsDca", &var_maxMuonsDca);
  readerA->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerA->AddVariable( "var_IsoPhiKKMass_Mu1", &var_IsoPhiKKMass_Mu1);
  readerA->AddVariable( "var_IsoPhiKKMass_Mu2", &var_IsoPhiKKMass_Mu2);
  readerA->AddVariable( "var_IsoPhiKKMass_Mu3", &var_IsoPhiKKMass_Mu3);
  readerA->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerA->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerA->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerA->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);


  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_2_A/weights/TMVAClassification_BDT.weights.xml"); 




  readerB = new TMVA::Reader( "!Color:!Silent" );
  readerB->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2);
  readerB->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerB->AddVariable( "var_flightLenSig", &var_flightLenSig);
  readerB->AddVariable( "var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerB->AddVariable( "var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerB->AddVariable( "var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);
  readerB->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerB->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerB->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerB->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerB->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerB->AddVariable( "var_IsoPhiKKMass_Mu1", &var_IsoPhiKKMass_Mu1);
  readerB->AddVariable( "var_IsoPhiKKMass_Mu2", &var_IsoPhiKKMass_Mu2);
  readerB->AddVariable( "var_IsoPhiKKMass_Mu3", &var_IsoPhiKKMass_Mu3);
  readerB->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerB->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_2_B/weights/TMVAClassification_BDT.weights.xml");



  readerC = new TMVA::Reader( "!Color:!Silent" );

  readerC->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle);
  readerC->AddVariable( "var_vertexKFChi2", &var_flightLenSig);
  readerC->AddVariable( "var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerC->AddVariable( "var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerC->AddVariable( "var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);
  readerC->AddVariable( "var_MaxD0SigSV", &var_MaxD0SigSV);
  readerC->AddVariable( "var_maxMuonsDca", &var_maxMuonsDca);
  readerC->AddVariable( "var_MindcaTrackSV", &var_MindcaTrackSV);
  readerC->AddVariable( "var_IsoPhiKKMass_Mu1", &var_IsoPhiKKMass_Mu1);
  readerC->AddVariable( "var_IsoPhiKKMass_Mu2", &var_IsoPhiKKMass_Mu2);
  readerC->AddVariable( "var_IsoPhiKKMass_Mu3", &var_IsoPhiKKMass_Mu3);
  readerC->AddVariable( "var_Muon1DetID", &var_Muon1DetID);
  readerC->AddVariable( "var_Muon2DetID", &var_Muon2DetID);
  readerC->AddVariable( "var_Muon3DetID", &var_Muon3DetID);
  readerC->AddVariable( "var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerC->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_2_C/weights/TMVAClassification_BDT.weights.xml"); 





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


  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionHelixRefit=HConfig.GetTH1D(Name+"_TauMassResolutionHelixRefit","TauMassResolutionHelixRefit",50,-0.2,0.2,"Helix refit #Delta M_{#tau}  (reco - mc)/mc ","Events");



  TauMassA1 =HConfig.GetTH1D(Name+"_TauMassA1","#tau lepton mass",30,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitA1 =HConfig.GetTH1D(Name+"_TauMassRefitA1","Refit #tau lepton mass",30,1.5,2.1,"KF refit  M_{#tau} , GeV (A1)","Events");


  TauMassRefitABC1 =HConfig.GetTH1D(Name+"_TauMassRefitABC1","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2 =HConfig.GetTH1D(Name+"_TauMassRefitABC2","Refit #tau lepton mass",30,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2)","Events");



  MassPairFit =HConfig.GetTH1D(Name+"_MassPairFit","MassPairFit",80,0.5,4,"Vertex Mass, GeV (2 tracks)","Events");
  Chi2SquarePairFit=HConfig.GetTH1D(Name+"_Chi2SquarePairFit","Chi2SquarePairFit",100,0,2,"Vertex Chi2, (2 tracks)","Events");
  DistanceToSignalVertexPairFit=HConfig.GetTH1D(Name+"_DistanceToSignalVertexPairFit","DistanceToSignalVertexPairFit",50,0,5,"Distance between D vertex and Signal vertex pos (2 tracks)","Events");

  MassPairFit_Chi2  =HConfig.GetTH2D(Name+"_MassPairFit_Chi2","MassPairFit_Chi2",50,0.2,4,50,0,2,"Vertex Mass, GeV (2 tracks)","Vertex #chi^{2}");
  MassTriplFit_Chi2 =HConfig.GetTH2D(Name+"_MassTriplFit_Chi2","MassTriplFit_Chi2",50,0.2,4,50,0,2,"Vertex Mass, GeV (3 tracks)","Vertex #chi^{2}");


  IsoTrack_mcDr_vs_drToTau=HConfig.GetTH2D(Name+"_IsoTrack_mcDr_vs_drToTau","IsoTrack_mcDr_vs_drToTau",100,0.,1.,100,0.,1,"#Delta R wrt to mc","#Delta R wrt to #tau");

  DeltaR_D_SV_3tracks=HConfig.GetTH1D(Name+"_DeltaR_D_SV_3tracks","DeltaR_D_SV_3tracks",30,0,2.5,"#DeltaR(#vec{D} - #vec{SV}) (D and SV are w.r.t #vec{PV} (3 tracks)","Events");
  DeltaR_D_SV_2tracks=HConfig.GetTH1D(Name+"_DeltaR_D_SV_2tracks","DeltaR_D_SV_2tracks",30,0,2.5,"#DeltaR(#vec{D} - #vec{SV}) (D and SV are w.r.t #vec{PV} (2 tracks)","Events");


  Angle_D_SV_3tracks=HConfig.GetTH1D(Name+"_Angle_D_SV_3tracks","Angle_D_SV_3tracks",40,0,1.5,"Angle (#vec{D} - #vec{SV}) (D and SV are w.r.t #vec{PV} (3 tracks), rad","Events");
  Angle_D_SV_2tracks=HConfig.GetTH1D(Name+"_Angle_D_SV_2tracks","Angle_D_SV_2tracks",40,0,1.5,"Angle (#vec{D} - #vec{SV}) (D and SV are w.r.t #vec{PV} (2 tracks), rad","Events");


  MassTripleFit =HConfig.GetTH1D(Name+"_MassTripleFit","MassTripleFit",80,0.5,5,"Vertex Mass, GeV (3 tracks)","Events");
  Chi2SquareTripleFit=HConfig.GetTH1D(Name+"_Chi2SquareTripleFit","Chi2SquareTripleFit",100,0,5,"Vertex Chi2, (3 tracks)","Events");
  DistanceToSignalVertexTripleFit=HConfig.GetTH1D(Name+"_DistanceToSignalVertexTripleFit","DistanceToSignalVertexTripleFit",50,0,5,"Distance between D vertex and Signal vertex pos (3 tracks)","Events");





  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  VertexChi2KF_vs_HelixFit=HConfig.GetTH2D(Name+"_VertexChi2KF_vs_HelixFit","VertexChi2KF_vs_HelixFit",50,0,10,50,0,10,"Kalman Vertex #chi^{2}","Helix Vertex  Fitter #chi^{2}");


  DisplacementFromPV_vs_mass3tracks=HConfig.GetTH2D(Name+"_DisplacementFromPV_vs_mass3tracks","DisplacementFromPV_vs_mass3tracks",80,0,3,80,0.5,2.5,"Displacement from PV 3 tracks, cm","SV Mass, GeV");
  DisplacementFromPV_vs_mass2tracks=HConfig.GetTH2D(Name+"_DisplacementFromPV_vs_mass2tracks","DisplacementFromPV_vs_mass2tracks",80,0,3,80,0.5,2.5,"Displacement from PV 2 tracks, cm","SV Mass, GeV");


  KF_Helix_deltaX=HConfig.GetTH1D(Name+"_KF_Helix_deltaX","KF_Helix_deltaX",50,-0.05,0.05,"#Delta X, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaY=HConfig.GetTH1D(Name+"_KF_Helix_deltaY","KF_Helix_deltaY",50,-0.05,0.05,"#Delta Y, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaZ=HConfig.GetTH1D(Name+"_KF_Helix_deltaZ","KF_Helix_deltaZ",50,-0.05,0.05,"#Delta Z, cm (Helix Fitter - Kalman Fitter)","Events");

  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");


 
  KpiIsolationMass_OS=HConfig.GetTH1D(Name+"_KpiIsolationMass_OS","KpiIsolationMass_OS",100,0.6,1.8,"M_{1}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS1=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS1","KpiIsolationMass_SS1",100,0.6,1.8,"M_{2}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS2=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS2","KpiIsolationMass_SS2",100,0.6,1.8,"M_{3}(K#pi), GeV (comb. iso #pi)","");
 

  BbkgVertexHypothesis=HConfig.GetTH1D(Name+"_BbkgVertexHypothesis","BbkgVertexHypothesis",2,-0.5,1.5,"1st bin - 2#mu+track, 2nd bin 3#mu","Events");

  BetterMuMuVertex=HConfig.GetTH1D(Name+"_BetterMuMuVertex","BetterMuMuVertex",30,0,5,"vertex pair quality (close)","");
  WorseMuMuVertex=HConfig.GetTH1D(Name+"_WorseMuMuVertex","WorseMuMuVertex",30,0,5,"vertex pair quality (far)","");


  MassFit12 =HConfig.GetTH1D(Name+"_MassFit12","MassFit12",80,0.2,4,"Vertex Mass Fit12, GeV ","Events");
  MassFit13 =HConfig.GetTH1D(Name+"_MassFit13","MassFit13",80,0.2,4,"Vertex Mass Fit13, GeV ","Events");
  MassFit23 =HConfig.GetTH1D(Name+"_MassFit23","MassFit23",80,0.2,4,"Vertex Mass Fit23, GeV ","Events");
  MassFit123 =HConfig.GetTH1D(Name+"_MassFit123","MassFit123",80,0.2,4,"Vertex Mass Fit123, GeV ","Events");
  MassFitLeastChi2 =HConfig.GetTH1D(Name+"_MassFitLeastChi2","MassFitLeastChi2",80,0.2,4,"Vertex Mass 2(tracks) least #chi^{2}, GeV ","Events");


  PairMassHighestPlusBestVertexMass=HConfig.GetTH1D(Name+"_PairMassHighestPlusBestVertexMass","PairMassHighestPlusBestVertexMass",80,0.2,4,"PairMass Highest+BestVertex, GeV ","Events");
  PairMassHighestPlusSubHighest=HConfig.GetTH1D(Name+"_PairMassHighestPlusSubHighest","PairMassHighestPlusSubHighest",80,0.2,4,"PairMass Highest+SubHighest, GeV ","Events");
  TripleMassHighestPlusBestVertexMass=HConfig.GetTH1D(Name+"_TripleMassHighestPlusBestVertexMass","TripleMassHighestPlusBestVertexMass",80,0.2,4,"TripleMass Highest+BestVertex, GeV ","Events");

  PairMassHighestPlusBestVertexMassVsChi2=HConfig.GetTH2D(Name+"_PairMassHighestPlusBestVertexMassVsChi2","PairMassHighestPlusBestVertexMassVsChi2",80,0.2,4,100,0,5,"PairMass Highest+BestVertex, GeV ","#chi^{2}");
  TripleMassHighestPlusBestVertexMassVsChi2=HConfig.GetTH2D(Name+"_TripleMassHighestPlusBestVertexMassVsChi2","TripleMassHighestPlusBestVertexMassVsChi2",80,0.2,4,100,0,5,"TripleMass Highest+BestVertex, GeV ","#chi^{2}");


  PairMassHighestPlusBestVertexMassDRToTruth=HConfig.GetTH1D(Name+"_PairMassHighestPlusBestVertexMassDRToTruth","PairMassHighestPlusBestVertexMassDRToTruth",80,0.2,4,"PairMass Highest+BestVertex, GeV ","Events");
  PairMassHighestPlusSubHighestDRToTruth=HConfig.GetTH1D(Name+"_PairMassHighestPlusSubHighestDRToTruth","PairMassHighestPlusSubHighestDRToTruth",80,0.2,4,"PairMass Highest+SubHighest, GeV ","Events");




  HighestPtTrackDRToTruth=HConfig.GetTH1D(Name+"_HighestPtTrackDRToTruth","HighestPtTrackDRToTruth",40,0,0.05,"Highest iso pT matched to MC?  #DeltaR ","Events");
  SecondTrackByBestFitIndexDRToTruth=HConfig.GetTH1D(Name+"_SecondTrackByBestFitIndexDRToTruth","SecondTrackByBestFitIndexDRToTruth",40,0,0.5,"Best vertex iso track with highest pT matched to MC?  #DeltaR ","Events");
  SubHighestPtTrackDRToTruth=HConfig.GetTH1D(Name+"_SubHighestPtTrackDRToTruth","SubHighestPtTrackDRToTruth",40,0,0.05,"NextHighest iso pT matched to MC?  #DeltaR ","Events");
  TripleMassHighestPlusBestVertexMassDRToTruth=HConfig.GetTH1D(Name+"_TripleMassHighestPlusBestVertexMassDRToTruth","TripleMassHighestPlusBestVertexMassDRToTruth",80,0.2,4,"TripleMass Highest+BestVertex, GeV ","Events");
    




  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  VertexFits::Store_ExtraDist(){ 


  Extradist1d.push_back(&HighestPtTrackDRToTruth);
  Extradist1d.push_back(&SecondTrackByBestFitIndexDRToTruth);
  Extradist1d.push_back(&SubHighestPtTrackDRToTruth);


  Extradist1d.push_back(&PairMassHighestPlusBestVertexMass);
  Extradist1d.push_back(&PairMassHighestPlusSubHighest);
  Extradist1d.push_back(&TripleMassHighestPlusBestVertexMass);
  Extradist1d.push_back(&TripleMassHighestPlusBestVertexMassDRToTruth);

  Extradist1d.push_back(&PairMassHighestPlusBestVertexMassDRToTruth);
  Extradist1d.push_back(&PairMassHighestPlusSubHighestDRToTruth);


  Extradist2d.push_back(&PairMassHighestPlusBestVertexMassVsChi2);
  Extradist2d.push_back(&TripleMassHighestPlusBestVertexMassVsChi2);


  Extradist1d.push_back(&BbkgVertexHypothesis);

  Extradist1d.push_back(&MassPairFit);
  Extradist1d.push_back(&Chi2SquarePairFit);
  Extradist1d.push_back(&DistanceToSignalVertexPairFit);

  Extradist1d.push_back(&MassTripleFit);
  Extradist1d.push_back(&Chi2SquareTripleFit);
  Extradist1d.push_back(&DistanceToSignalVertexTripleFit);


  Extradist1d.push_back(&MassFit12);
  Extradist1d.push_back(&MassFit13);
  Extradist1d.push_back(&MassFit23);
  Extradist1d.push_back(&MassFit123);

  Extradist1d.push_back(&MassFitLeastChi2);

  Extradist1d.push_back(&DeltaR_D_SV_3tracks);
  Extradist1d.push_back(&DeltaR_D_SV_2tracks);

  Extradist1d.push_back(&Angle_D_SV_3tracks);
  Extradist1d.push_back(&Angle_D_SV_2tracks);


  Extradist2d.push_back(&MassPairFit_Chi2);
  Extradist2d.push_back(&MassTriplFit_Chi2);
  Extradist2d.push_back(&IsoTrack_mcDr_vs_drToTau);

  Extradist1d.push_back(&TauMassRefitABC1);  
  Extradist1d.push_back(&TauMassRefitABC2);

  Extradist1d.push_back(&BetterMuMuVertex);
  Extradist1d.push_back(&WorseMuMuVertex);


  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);
  Extradist1d.push_back(&TauMassResolutionHelixRefit);

  Extradist1d.push_back(&VertexChi2KF);
  Extradist2d.push_back(&VertexChi2KF_vs_HelixFit);
  Extradist1d.push_back(&KF_Helix_deltaX);
  Extradist1d.push_back(&KF_Helix_deltaY);
  Extradist1d.push_back(&KF_Helix_deltaZ);


  Extradist2d.push_back(&DisplacementFromPV_vs_mass3tracks);
  Extradist2d.push_back(&DisplacementFromPV_vs_mass2tracks);

  Extradist1d.push_back(&KpiIsolationMass_OS);
  Extradist1d.push_back(&KpiIsolationMass_SS1);
  Extradist1d.push_back(&KpiIsolationMass_SS2);



}


void  VertexFits::doEvent(){ 

  
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
  pass.at(PhiVeto1) = (value.at(PhiVeto1) < 0.98 || value.at(PhiVeto1) > 1.06 );
  pass.at(OmegaVeto1) =(value.at(OmegaVeto1) < 0.742 || value.at(OmegaVeto1) > 0.822 );
  pass.at(PhiVeto2) = (value.at(PhiVeto2) < 0.98 || value.at(PhiVeto2) > 1.06 );
  pass.at(OmegaVeto2) =(value.at(OmegaVeto2) < 0.742 || value.at(OmegaVeto2) > 0.822 );
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


  }


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);

  if(status){
    //    std::cout<<"id ---- "<< id << std::endl;
    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);


    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);  
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);

    std::vector<unsigned int> EtaSortedIndices;
    
    EtaSortedIndices.push_back(Muon_Eta_index_1);
    EtaSortedIndices.push_back(Muon_Eta_index_2);
    EtaSortedIndices.push_back(Muon_Eta_index_3);


    //    std::cout<<" id --------- "<< id << std::endl;
    float DChi2(9999999.);
    float DChi2_Triple(9999999.);
    int TrackIndex1(-1);
    int TrackIndex2(-1);
    int TrackIndex3(-1);

    std::cout<<"all iso tracks:   "<< std::endl;
    int DToKPiPi_index(-1);
    if(id!=1){


      //------------------ study the MC content ------------
      /*
      for(unsigned int imc=0; imc< Ntp->NMCParticles(); imc++){
	if( abs(Ntp->MCParticle_pdgid(imc)) == 511 ){
	  std::cout<<" B0/antiB0  found  "<<   Ntp->MCSignalParticle_Nchilds(imc) <<std::endl;

	  for(unsigned int dau_mc = 0; dau_mc < Ntp->MCSignalParticle_Nchilds(imc); dau_mc++){
	    std::cout<<" Childs of B  "<< std::endl;
	    std::cout<<"  "<< Ntp->MCSignalParticle_childpdgid(imc, dau_mc) << std::endl;

	  }
	}
      }

      */






    for(unsigned int imc=0; imc< Ntp->NMCParticles(); imc++){

      if(abs(Ntp->MCParticle_pdgid(imc)) == 411 or abs(Ntp->MCParticle_pdgid(imc)) == 421){
	std::cout<<" found particle  "<< Ntp->MCParticle_pdgid(imc) << " size  " <<Ntp->MCParticle_childpdgid(imc).size() <<std::endl;
	if(Ntp->MCParticle_childpdgid(imc).size()>=2){
	  
	  for (int i = 0; i < Ntp->MCParticle_childpdgid(imc).size(); i++) {
	    std::cout << Ntp->MCParticle_childpdgid(imc).at(i) << ' ';
	  }

	  //	  if((std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , 211)!=Ntp->MCParticle_childpdgid(imc).end()   or 
	  //	      std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , -211)!=Ntp->MCParticle_childpdgid(imc).end()   ) and 
	  //	     (std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , 321)!=Ntp->MCParticle_childpdgid(imc).end()   or 
	  //	      std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , -321)!=Ntp->MCParticle_childpdgid(imc).end()))

	  if((std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , 321)!=Ntp->MCParticle_childpdgid(imc).end()   or 
	      std::find(Ntp->MCParticle_childpdgid(imc).begin() , Ntp->MCParticle_childpdgid(imc).end() , -321)!=Ntp->MCParticle_childpdgid(imc).end()))
	    {
	      DToKPiPi_index = imc;	       
	      //	      std::cout<<"====================================================== "<< DToKPiPi_index  <<std::endl;
	    }
	  
	}
      }
    }

    /*    if(DToKPiPi_index>0){
      
      //      std::cout<<"------------------      DToKPiPi found " << std::endl;
      Ntp->printMCDecayChainOfParticle(DToKPiPi_index,true,true,true,true);

      }*/
    }

    for(unsigned int i1=0; i1 < Ntp->NIsolationTrack(signal_idx); i1++){// first loop

      if(id!=1)     std::cout<<"i:   "<<  i1 << "   pt:   "<<Ntp->IsolationTrack_p4(signal_idx, i1).Pt() << "  eta:    " << Ntp->IsolationTrack_p4(signal_idx, i1).Eta() 
         	       <<  "   truth pt:   "<<Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1))).Pt() 
          	       <<  "   truth eta:  "<<Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1))).Eta() <<"   id   " <<  Ntp->MCParticle_pdgid(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1)))  <<"  dR   " <<   Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1)) <<" dR to tau   " <<  TauLV.DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1))<<" delta Phi  " << Ntp->DeltaPhi(TauLV.Phi(), Ntp->IsolationTrack_p4(signal_idx, i1).Phi()) <<std::endl;

      

      if(id!=1)      IsoTrack_mcDr_vs_drToTau.at(t).Fill(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1)), TauLV.DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1)));




      //      if(id!=1)      if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i1))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1)) > 0.003) continue;// lets take only matched

      if(Ntp->IsolationTrack_p4(signal_idx, i1).Pt() < 0.4 ) continue;

      for(unsigned int i2=i1+1; i2 < Ntp->NIsolationTrack(signal_idx); i2++){// second loop
	if(Ntp->IsolationTrack_p4(signal_idx, i2).Pt() < 0.4 ) continue;
	//	if(Ntp->IsolationTrack_p4(signal_idx, i2).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i1)) > 1.5) continue;
	//	if(id!=1)	if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i2))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i2)) > 0.003) continue;// lets take only matched
	if(Ntp->IsolationTrack_Helcharge(i1)!=Ntp->IsolationTrack_Helcharge(i2)){

	  std::vector<TrackParticle> FirstPairCollection;
	  FirstPairCollection.push_back(Ntp->IsolationTrack_TrackParticle(i1));
	  FirstPairCollection.push_back(Ntp->IsolationTrack_TrackParticle(i2));
	  TVector3 vguess(0,0,0);
	  Chi2VertexFitter  PairFit(FirstPairCollection,vguess);
	  PairFit.Fit();
	  if(PairFit.ChiSquare() > 0 && PairFit.ChiSquare() < 999999. ){
	    if(PairFit.ChiSquare() < DChi2){
	      DChi2  = PairFit.ChiSquare();
	      TrackIndex1 = i1;
	      TrackIndex2 = i2;
	    }
	  }
	}
      }
    }
    if(TrackIndex1!=-1 && TrackIndex2!=-1){
      if(Ntp->IsolationTrack_NTracks() > 2){
	for(unsigned int i3=0; i3 < Ntp->NIsolationTrack(signal_idx); i3++){
	  if(i3==TrackIndex1) continue;
	  if(i3==TrackIndex2) continue;
	  if(Ntp->IsolationTrack_p4(signal_idx, i3).Pt() < 0.4 ) continue;
	  //	  if(Ntp->IsolationTrack_p4(signal_idx, i3).DeltaR(Ntp->IsolationTrack_p4(signal_idx, TrackIndex1)) > 1.5
	  //	     or Ntp->IsolationTrack_p4(signal_idx, i3).DeltaR(Ntp->IsolationTrack_p4(signal_idx, TrackIndex2)) > 1.5) continue;
	  //	  if(id!=1)	  if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, i3))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, i3)) > 0.003) continue;// lets take only matched
	  std::vector<TrackParticle> TripleCollection;
	  TripleCollection.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex1));
	  TripleCollection.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex2));
	  Tripl.push_back(Ntp->IsolationTrack_TrackParticle(i3));


	  TVector3 vguess(0,0,0);
	  Chi2VertexFitter  TripleFit(TripleCollection,vguess);
	  TripleFit.Fit();

	  if(TripleFit.ChiSquare() > 0 && TripleFit.ChiSquare() < 999999. ){
	    if(TripleFit.ChiSquare() < DChi2_Triple){
	      DChi2_Triple  = TripleFit.ChiSquare();
	      TrackIndex3 = i3;
	    }
	  }
	}
      }
    }

    //  ----------------------------- secondary vertices ----------------
    int NumberOfPrimaryVertices(0);
    for(unsigned int iVertex=0; iVertex < Ntp->NSecondaryVertices(); iVertex++){

      bool nonSignalMuons(true);
      for(unsigned int iTrack = 0; iTrack < Ntp->NTracksAtSecondaryVertex(iVertex); iTrack++){
        if(Muon1LV.DeltaR(Ntp->SecondaryVertexTrack_P4(iVertex, iTrack)) < 0.01 or
           Muon2LV.DeltaR(Ntp->SecondaryVertexTrack_P4(iVertex, iTrack)) < 0.01 or
           Muon3LV.DeltaR(Ntp->SecondaryVertexTrack_P4(iVertex, iTrack)) < 0.01) {nonSignalMuons = false; continue; }
	std::cout<<"Vertex "<< iVertex << " TrackPt  "<< Ntp->SecondaryVertexTrack_P4(iVertex, iTrack).Pt() << " eta  "<< Ntp->SecondaryVertexTrack_P4(iVertex, iTrack).Eta() <<std::endl;
      }


      if(nonSignalMuons){
	if(Ntp->NTracksAtSecondaryVertex(iVertex)==3){
	  if(Ntp->SecondaryVertexTrackCharge(iVertex,0) * Ntp->SecondaryVertexTrackCharge(iVertex,1)*Ntp->SecondaryVertexTrackCharge(iVertex,2)==-1 ){
	    DisplacementFromPV_vs_mass3tracks.at(t).Fill(  (Ntp->Vertex_MatchedPrimaryVertex(signal_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag(), Ntp->SecondaryVertexMass(iVertex)  );
	  }
	}
	if( Ntp->NTracksAtSecondaryVertex(iVertex)==2){
	  if(Ntp->SecondaryVertexTrackCharge(iVertex,0) * Ntp->SecondaryVertexTrackCharge(iVertex,1)==-1){
	    DisplacementFromPV_vs_mass2tracks.at(t).Fill(  (Ntp->Vertex_MatchedPrimaryVertex(signal_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag(), Ntp->SecondaryVertexMass(iVertex)  );
	  }
	}
      }
    }


    if(TrackIndex1!=-1 && TrackIndex2 !=-1){
      std::vector<TrackParticle> TwoTracks;
      TwoTracks.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex1));
      TwoTracks.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex2));
      TVector3 vguess(0,0,0);
      Chi2VertexFitter  TwoTracksFit(TwoTracks,vguess);
      TwoTracksFit.Fit();


      TVector3 ReffitedVertex2Particles = TwoTracksFit.GetVertex();
      LorentzVectorParticle  MotherParticle2Particles = TwoTracksFit.GetMother(888);


      MassPairFit.at(t).Fill(MotherParticle2Particles.LV().M(),1);
      Chi2SquarePairFit.at(t).Fill(TwoTracksFit.ChiSquare(),1);
      //      DistanceToSignalVertexPairFit.at(t).Fill( (Ntp->Vertex_Signal_KF_pos(signal_idx) - ReffitedVertex2Particles).Mag(),1);

      TVector3 pvTosv_2tracks = Ntp->Vertex_Signal_KF_pos(signal_idx) -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
      TVector3 pvToD_2tracks  = ReffitedVertex2Particles -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
      //      DeltaR_D_SV_2tracks.at(t).Fill(pvTosv_2tracks.DeltaR(pvToD_2tracks),1);
      //      Angle_D_SV_2tracks.at(t).Fill(pvTosv_2tracks.Angle(pvToD_2tracks),1);

      if(pvTosv_2tracks.Angle(pvToD_2tracks) > 0.5*TMath::Pi()){

	std::cout<<"iso 2 - 1: "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Pt()<< "  eta:    "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Eta() << " phi:  " << Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Phi() <<std::endl;
	std::cout<<"iso 2 - 2: "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Pt()<< "  eta:    "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Eta() << " phi:  " << Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Phi()  <<std::endl;
	
	/*	if(id == 40 or id == 60 or id == 90 ){// or id == 40){
	    std::cout<<"-------------- All categoris ----------------"<< std::endl;
	    std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	    std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	    std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	    Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	    Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	    Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	    Ntp->printMCDecayChainOfEvent(true, true, true, true);
	    std::cout<< "\n\n\n\n\n\n";
	    }*/

	if(DToKPiPi_index!=-1 && id !=1){
	  for (unsigned int i_dau = 0; i_dau < Ntp->MCParticle_childidx(DToKPiPi_index).size(); i_dau++){
	    std::cout<<"  delta Eta to Tau   "<<Ntp->MCParticle_p4(Ntp->MCParticle_childidx(DToKPiPi_index).at(i_dau)).Eta()  - TauLV.Eta() <<"  pdgid  "  
		     << PDGInfo::pdgIdToName(Ntp->MCParticle_pdgid(Ntp->MCParticle_childidx(DToKPiPi_index).at(i_dau)))<<" pT  " << Ntp->MCParticle_p4(Ntp->MCParticle_childidx(DToKPiPi_index).at(i_dau)).Pt() << "  eta  " << Ntp->MCParticle_p4(Ntp->MCParticle_childidx(DToKPiPi_index).at(i_dau)).Eta() << std::endl;
	  }
	}
      
      }
    
      MassPairFit_Chi2.at(t).Fill(MotherParticle2Particles.LV().M(),TwoTracksFit.ChiSquare(),1);
    }

    if(TrackIndex1!=-1 && TrackIndex2 !=-1 && TrackIndex3!=-1){
      std::vector<TrackParticle> ThreeTracks;
      ThreeTracks.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex1));
      ThreeTracks.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex2));
      ThreeTracks.push_back(Ntp->IsolationTrack_TrackParticle(TrackIndex3));
      TVector3 vguess(0,0,0);
      Chi2VertexFitter  ThreeTracksFit(ThreeTracks,vguess);
      ThreeTracksFit.Fit();


      TVector3 ReffitedVertex3Particles = ThreeTracksFit.GetVertex();
      LorentzVectorParticle  MotherParticle3Particles = ThreeTracksFit.GetMother(888);


      TVector3 pvTosv_3tracks = Ntp->Vertex_Signal_KF_pos(signal_idx) -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
      TVector3 pvToD_3tracks  = ReffitedVertex3Particles -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);

      //      DeltaR_D_SV_3tracks.at(t).Fill(pvTosv_3tracks.DeltaR(pvToD_3tracks),1);
      //      Angle_D_SV_3tracks.at(t).Fill(pvTosv_3tracks.Angle(pvToD_3tracks),1);


      if(pvTosv_3tracks.Angle(pvToD_3tracks) > 0.5*TMath::Pi()){


	std::cout<<"iso 3 - 1: "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Pt()<< "  eta:    "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Eta() << " phi:   " <<Ntp->IsolationTrack_p4(signal_idx, TrackIndex1).Phi() <<std::endl;
	std::cout<<"iso 3 - 2: "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Pt()<< "  eta:    "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Eta() << " phi:   "<<Ntp->IsolationTrack_p4(signal_idx, TrackIndex2).Phi() <<std::endl;
       	std::cout<<"iso 3x - 3: "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex3).Pt()<< "  eta:    "<< Ntp->IsolationTrack_p4(signal_idx, TrackIndex3).Eta() << " phi:   "<<Ntp->IsolationTrack_p4(signal_idx, TrackIndex3).Phi() <<std::endl;

	/*	  if(id == 40 or id == 60 or id == 90 ){// or id == 40){

	    std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	    std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	    std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	    Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	    Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	    Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	    Ntp->printMCDecayChainOfEvent(true, true, true, true);
	    std::cout<< "\n\n\n\n\n\n";
	    }*/
      }



      MassTripleFit.at(t).Fill(MotherParticle3Particles.LV().M(),1);
      Chi2SquareTripleFit.at(t).Fill(ThreeTracksFit.ChiSquare(),1);
      DistanceToSignalVertexTripleFit.at(t).Fill( (Ntp->Vertex_Signal_KF_pos(signal_idx) - ReffitedVertex3Particles).Mag(),1);
      MassTriplFit_Chi2.at(t).Fill(MotherParticle3Particles.LV().M(),ThreeTracksFit.ChiSquare());
    }




    if(id == 40 or id == 60 or id == 90 ){// or id == 40){
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
    //*********  Tracking Other Approach

    double deltaR1(999.);
    double deltaR2(999.);
    double deltaR3(999.);
    int FirstTrack_index;
    int SecondTrack_index;
    int ThirdTrack_index;
    bool FirstTrack_found(false);
    bool SecondTrack_found(false);
    bool ThirdTrack_found(false);
    int HighestPtTrack(-1);
    int SubHighestPtTrack(-1);
    double MaxPT(.01);
    for(unsigned int i1=0; i1 < Ntp->NIsolationTrack(signal_idx); i1++){// first loop
      //first find the track closest to tau by dR
      if(Ntp->IsolationTrack_p4(signal_idx, i1).DeltaR(TauLV) < deltaR1){
	deltaR1 = Ntp->IsolationTrack_p4(signal_idx, i1).DeltaR(TauLV);
	FirstTrack_index = i1;
	FirstTrack_found = true;
      }  
      if(Ntp->IsolationTrack_p4(signal_idx, i1).DeltaR(TauLV) < 0.5){
	if(Ntp->IsolationTrack_p4(signal_idx, i1).Pt() > MaxPT){ 
	  MaxPT = Ntp->IsolationTrack_p4(signal_idx, i1).Pt();
	  HighestPtTrack = i1;
	}
      }
    }
    //    std::cout<<"deb4"<<std::endl;
    MaxPT = .01;
    double deltaChi2(9e+9);
    int SecondTrackByBestFitIndex(-1);
    if(HighestPtTrack!=-1 ){
      for(unsigned int i2=0; i2 < Ntp->NIsolationTrack(signal_idx); i2++){// 
	if(i2==HighestPtTrack) continue;
	if( Ntp->IsolationTrack_p4(signal_idx, i2).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))  < 0.35 ){
	  TVector3 v(0,0,0);
	  std::vector<TrackParticle> TracksToFit;
	  TracksToFit.push_back(Ntp->IsolationTrack_TrackParticle(FirstTrack_index));
	  TracksToFit.push_back(Ntp->IsolationTrack_TrackParticle(i2));
	  Chi2VertexFitter  vfit(TracksToFit,v);
	  bool fitOk = vfit.Fit();

	  //	  std::cout<<" pT  "<< Ntp->IsolationTrack_p4(signal_idx, i2).Pt() << "  dR to first  found   "
	  //		   <<  Ntp->IsolationTrack_p4(signal_idx, i2).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack)) <<"   chi2  " << vfit.ChiSquare() <<"  ok?   " << fitOk <<std::endl;
	  
	  if(fitOk){
	    if( vfit.ChiSquare()< deltaChi2){

	      deltaChi2 = vfit.ChiSquare();
	      SecondTrackByBestFitIndex = i2;

	    }
	  }
  
	  if(Ntp->IsolationTrack_p4(signal_idx, i2).Pt() > MaxPT){
	    MaxPT = Ntp->IsolationTrack_p4(signal_idx, i2).Pt();
	    SubHighestPtTrack = i2;
	  }
	}
      }
    }
      
    // ----------------------------------   debugging fit -----------------

    std::cout<<" ========================  "<<std::endl;
    std::cout<<"SubHighestPtTrack   "<< SubHighestPtTrack << "   HighestPtTrack   " << HighestPtTrack << std::endl;

    if(SubHighestPtTrack!=-1 && HighestPtTrack!=-1){


      std::vector<TrackParticle> DebugFitTracks;
      DebugFitTracks.push_back(Ntp->IsolationTrack_TrackParticle(HighestPtTrack));
      DebugFitTracks.push_back(Ntp->IsolationTrack_TrackParticle(SubHighestPtTrack));
      TVector3 vguess(0.1,0.1,0.1);
      Chi2VertexFitter  DebugFitter(DebugFitTracks,vguess);
      DebugFitter.Fit();

      std::cout<<"  Debugger chi2    "<< DebugFitter.ChiSquare() << std::endl;
      std::cout<<"  Vertex Position    "<< std::endl;
      DebugFitter.GetVertex().Print();

    }


    //--------------------------------------------------------------------



    if(FirstTrack_found){
      for(unsigned int i2=0; i2 < Ntp->NIsolationTrack(signal_idx); i2++){// second loop
	if(i2 == FirstTrack_index) continue;
	if(Ntp->IsolationTrack_p4(signal_idx, i2).DeltaR(Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index)) < deltaR2){
	  
	  deltaR2 = Ntp->IsolationTrack_p4(signal_idx, i2).DeltaR(Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index));
	  SecondTrack_index=i2;
	  SecondTrack_found = true;
	}
      }
      
      if(SecondTrack_found){
	for(unsigned int i3=0; i3 < Ntp->NIsolationTrack(signal_idx); i3++){// second loop
	  if(i3 == FirstTrack_index) continue;
	  if(i3 == SecondTrack_index) continue;
	  if(Ntp->IsolationTrack_p4(signal_idx, i3).DeltaR(Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index) + Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index) ) < deltaR3){
	    deltaR3 = Ntp->IsolationTrack_p4(signal_idx, i3).DeltaR(Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index) + Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index) );
	    ThirdTrack_index = i3;
	    ThirdTrack_found = true;
	  }
	}

	if(ThirdTrack_found && SubHighestPtTrack!=-1 && HighestPtTrack!=-1){
	  /*	  std::cout<<"index1  "<< FirstTrack_index << " dR to tau "<< Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index).DeltaR(TauLV) << " pt  " << Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index).Pt()  <<"  eta  "<<Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index).Eta() <<std::endl;
	  std::cout<<"index2  "<< SecondTrack_index << " dR to tau  "<< Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index).DeltaR(TauLV) << " pt  " << Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index).Pt() << "  eta  "<<Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index).Eta() <<std::endl;
	  std::cout<<"index3  "<< ThirdTrack_index << " dR to tau  "<< Ntp->IsolationTrack_p4(signal_idx, ThirdTrack_index).DeltaR(TauLV) << " pt  " << Ntp->IsolationTrack_p4(signal_idx, ThirdTrack_index).Pt()<< "  eta  "<<Ntp->IsolationTrack_p4(signal_idx, ThirdTrack_index).Eta() <<std::endl;

	  std::cout<<"Highest pT  "<< HighestPtTrack << " dR to tau  "<< Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack).DeltaR(TauLV) << " pt  " <<Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack).Pt() <<std::endl;
	  std::cout<<"Second Highest pT  "<< SubHighestPtTrack << " dR to tau  "<< Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack).DeltaR(TauLV) << " pt  " <<Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack).Pt() <<"  dR to first track "  << Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))<<std::endl;*/

	}
      }
    }

    if(SecondTrackByBestFitIndex!=-1 && HighestPtTrack!=-1){


      if(id!=1){
	HighestPtTrackDRToTruth.at(t).Fill(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack)),1);
	SecondTrackByBestFitIndexDRToTruth.at(t).Fill(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex)),1);

	if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack)) < 0.005 && 
	   Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex)) < 0.005)
	  {
	    PairMassHighestPlusBestVertexMassDRToTruth.at(t).Fill( (Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack) + Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex)).M(),1);
	  }
      }
      PairMassHighestPlusBestVertexMass.at(t).Fill( (Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack) + Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex)).M(),1);
      PairMassHighestPlusBestVertexMassVsChi2.at(t).Fill( (Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack) + Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex)).M(),deltaChi2);
      
      //-----------------



      std::vector<TrackParticle> TwoTrackToFit;
      TwoTrackToFit.push_back(Ntp->IsolationTrack_TrackParticle(HighestPtTrack));
      TwoTrackToFit.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrackByBestFitIndex));
      TVector3 vguess(0,0,0);
      Chi2VertexFitter  TwoTracksFit(TwoTrackToFit,vguess);
      TwoTracksFit.Fit();
      TVector3 ReffitedVertex2Particles = TwoTracksFit.GetVertex();


      TVector3 pvTosv_2tracks = Ntp->Vertex_Signal_KF_pos(signal_idx) -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
      TVector3 pvToD_2tracks  = ReffitedVertex2Particles -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);

      DeltaR_D_SV_2tracks.at(t).Fill(pvTosv_2tracks.DeltaR(pvToD_2tracks),1);
      Angle_D_SV_2tracks.at(t).Fill(pvTosv_2tracks.Angle(pvToD_2tracks),1);

      DistanceToSignalVertexPairFit.at(t).Fill( (Ntp->Vertex_Signal_KF_pos(signal_idx) - ReffitedVertex2Particles).Mag(),1);



      //-----------------
      double TempChi2(9e+9);
	int ThirdTrackIndexByBestVertex(-1);
	double TripletMassFit;
        for(unsigned int i3=0; i3 < Ntp->NIsolationTrack(signal_idx); i3++){// third loop
          if(i3 == HighestPtTrack) continue;
          if(i3 == SecondTrackByBestFitIndex) continue;


          TVector3 v(0,0,0);
	  std::vector<TrackParticle> TracksToFit;
          TracksToFit.push_back(Ntp->IsolationTrack_TrackParticle(HighestPtTrack));
          TracksToFit.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrackByBestFitIndex));
          TracksToFit.push_back(Ntp->IsolationTrack_TrackParticle(i3));
          Chi2VertexFitter  vfit(TracksToFit,v);
          bool fitOk = vfit.Fit();

          if(fitOk){
            if( vfit.ChiSquare()< TempChi2){
              TempChi2 = vfit.ChiSquare();
              ThirdTrackIndexByBestVertex = i3;
	      TripletMassFit = vfit.GetMother(888).LV().M();
            }
          }
	}

	if(ThirdTrackIndexByBestVertex!=-1) {
	  TripleMassHighestPlusBestVertexMass.at(t).Fill( TripletMassFit ,1);
	  TripleMassHighestPlusBestVertexMassVsChi2.at(t).Fill( TripletMassFit ,TempChi2);



	  //-----------------------------------
	  std::vector<TrackParticle> ThreeTrackToFit;
	  ThreeTrackToFit.push_back(Ntp->IsolationTrack_TrackParticle(HighestPtTrack));
	  ThreeTrackToFit.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrackByBestFitIndex));
	  ThreeTrackToFit.push_back(Ntp->IsolationTrack_TrackParticle(ThirdTrackIndexByBestVertex));

	  TVector3 vguess(0,0,0);
	  Chi2VertexFitter  ThreeTracksFit(ThreeTrackToFit,vguess);
	  ThreeTracksFit.Fit();
	  TVector3 ReffitedVertex3Particles = ThreeTracksFit.GetVertex();


	  TVector3 pvTosv_3tracks = Ntp->Vertex_Signal_KF_pos(signal_idx) -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);
	  TVector3 pvToD_3tracks  = ReffitedVertex3Particles -Ntp->Vertex_MatchedPrimaryVertex(signal_idx);

	  DeltaR_D_SV_3tracks.at(t).Fill(pvTosv_3tracks.DeltaR(pvToD_3tracks),1);
	  Angle_D_SV_3tracks.at(t).Fill(pvTosv_3tracks.Angle(pvToD_3tracks),1);


	  //-------------------------------------



	  if(ThirdTrackIndexByBestVertex!=-1)
	    {
	    if(id!=1)
	      {
		if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, ThirdTrackIndexByBestVertex))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, ThirdTrackIndexByBestVertex))<0.005&& 
		   Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, SecondTrackByBestFitIndex))<0.005&&
		   Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))<0.005)
		  {
		    TripleMassHighestPlusBestVertexMassDRToTruth.at(t).Fill(TripletMassFit,1);
		  }
	      }
	    }

	}
    }

    


    if(SubHighestPtTrack!=-1 && HighestPtTrack!=-1){


      if(id!=1){      
	SubHighestPtTrackDRToTruth.at(t).Fill(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack)),1);
	if(Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack))<0.005 && 
	   Ntp->MCParticle_p4(Ntp->getMatchTruthIndex(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))).DeltaR(Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack))<0.005)
	  {
	    PairMassHighestPlusSubHighestDRToTruth.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack) + Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack)).M(),1);
	  }
      }

      PairMassHighestPlusSubHighest.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx, HighestPtTrack) + Ntp->IsolationTrack_p4(signal_idx, SubHighestPtTrack)).M(),1);
    }


    if(FirstTrack_found && SecondTrack_found && ThirdTrack_found){
      TVector3 v(0,0,0);
      std::vector<TrackParticle> Tracks12;
      Tracks12.push_back(Ntp->IsolationTrack_TrackParticle(FirstTrack_index));
      Tracks12.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrack_index));
      
      std::vector<TrackParticle> Tracks13;
      Tracks13.push_back(Ntp->IsolationTrack_TrackParticle(FirstTrack_index));
      Tracks13.push_back(Ntp->IsolationTrack_TrackParticle(ThirdTrack_index));

      std::vector<TrackParticle> Tracks23;
      Tracks23.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrack_index));
      Tracks23.push_back(Ntp->IsolationTrack_TrackParticle(ThirdTrack_index));
            

     std::vector<TrackParticle> Tracks123;
     Tracks123.push_back(Ntp->IsolationTrack_TrackParticle(FirstTrack_index));
     Tracks123.push_back(Ntp->IsolationTrack_TrackParticle(SecondTrack_index));
     Tracks123.push_back(Ntp->IsolationTrack_TrackParticle(ThirdTrack_index));
 

      Chi2VertexFitter  Fit12(Tracks12,v);
      Chi2VertexFitter  Fit13(Tracks13,v);
      Chi2VertexFitter  Fit23(Tracks23,v);
      Chi2VertexFitter  Fit123(Tracks123,v);
   
      
      if(Fit12.Fit() &&  Fit13.Fit() &&  Fit23.Fit()){
	//	std::cout<<"12 chi2  "<< Fit12.ChiSquare() << "  13 chi2  "<< Fit13.ChiSquare() << "  23 chi2   " << Fit23.ChiSquare()<< std::endl;
	MassFit12.at(t).Fill(Fit12.GetMother(888).LV().M(),1);
	MassFit13.at(t).Fill(Fit13.GetMother(888).LV().M(),1);
	MassFit23.at(t).Fill(Fit23.GetMother(888).LV().M(),1);
	double LeastChi2 = std::min({Fit12.ChiSquare(), Fit13.ChiSquare(),Fit23.ChiSquare()});
	if(LeastChi2 == Fit12.ChiSquare()) MassFitLeastChi2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index) + Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index) ).M(),1);
	if(LeastChi2 == Fit13.ChiSquare()) MassFitLeastChi2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx, FirstTrack_index) + Ntp->IsolationTrack_p4(signal_idx, ThirdTrack_index) ).M(),1);
	if(LeastChi2 == Fit23.ChiSquare()) MassFitLeastChi2.at(t).Fill((Ntp->IsolationTrack_p4(signal_idx, SecondTrack_index) + Ntp->IsolationTrack_p4(signal_idx, ThirdTrack_index) ).M(),1);

      }
      if(Fit123.Fit()) {	MassFit123.at(t).Fill(Fit123.GetMother(888).LV().M(),1); }
      
    }
  



    TVector3 vguess(0,0,0);
    std::vector<TrackParticle> MuonsTrackParticles;
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_3));

    Chi2VertexFitter  Fitter(MuonsTrackParticles,vguess);
    LorentzVectorParticle  MotherParticle;// = Fitter.GetMother(888);

    if(Fitter.Fit()){

      VertexChi2KF_vs_HelixFit.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(signal_idx),Fitter.ChiSquare());
      TVector3 ReffitedVertex = Fitter.GetVertex();


      KF_Helix_deltaX.at(t).Fill(ReffitedVertex.X()  - Ntp->Vertex_Signal_KF_pos(signal_idx).X(),1);
      KF_Helix_deltaY.at(t).Fill(ReffitedVertex.Y()  - Ntp->Vertex_Signal_KF_pos(signal_idx).Y(),1);
      KF_Helix_deltaZ.at(t).Fill(ReffitedVertex.Z()  - Ntp->Vertex_Signal_KF_pos(signal_idx).Z(),1);
      
      std::vector<LorentzVectorParticle> ReffitedLVParticles = Fitter.GetReFitLorentzVectorParticles();
      //    std::cout<<" pt1  "<< Ntp->Muon_P4(Muon_index_1).Pt() << "  "<< ReffitedLVParticles.at(0).LV().Pt() << std::endl;
      MotherParticle = Fitter.GetMother(888);
    }

    std::vector<TrackParticle> MuonsTrackParticles12;
    MuonsTrackParticles12.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
    MuonsTrackParticles12.push_back(Ntp->Muon_TrackParticle(Muon_index_2));

    std::vector<TrackParticle> MuonsTrackParticles13;
    MuonsTrackParticles13.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
    MuonsTrackParticles13.push_back(Ntp->Muon_TrackParticle(Muon_index_3));


    std::vector<TrackParticle> MuonsTrackParticles23;
    MuonsTrackParticles23.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
    MuonsTrackParticles23.push_back(Ntp->Muon_TrackParticle(Muon_index_3));


    Chi2VertexFitter  Fitter12(MuonsTrackParticles12,vguess);
    Chi2VertexFitter  Fitter13(MuonsTrackParticles13,vguess);
    Chi2VertexFitter  Fitter23(MuonsTrackParticles23,vguess);


    if(Fitter12.Fit() && Fitter13.Fit()  && Fitter23.Fit()){
      std::vector<TrackParticle> MuonsTrackParticles_BestPair;

      double MinPairMuVertexChi2 = std::min({Fitter12.ChiSquare(), Fitter13.ChiSquare(),Fitter23.ChiSquare()});

      //      std::cout<<"    "<<Fitter12.ChiSquare()<<"  " << Fitter13.ChiSquare()<< "  "<< Fitter23.ChiSquare() << std::endl;
      //      std::cout<<"min:  "<< MinPairMuVertexChi2 << std::endl;

      if(MinPairMuVertexChi2 == Fitter12.ChiSquare()) {
	MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
	MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
      }
      if(MinPairMuVertexChi2 == Fitter13.ChiSquare()){
	MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
	MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_3));
      }
      if(MinPairMuVertexChi2 == Fitter23.ChiSquare()){
        MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
        MuonsTrackParticles_BestPair.push_back(Ntp->Muon_TrackParticle(Muon_index_3));
       }


      
      Chi2VertexFitter FitterBestMuonPair(MuonsTrackParticles_BestPair,vguess);
      if(FitterBestMuonPair.Fit()){

	int BestIsoTrackAddedToMuonPair = -1;
	double tempChi2(10e+7);
	for(unsigned int iso_Track=0; iso_Track < Ntp->NIsolationTrack(signal_idx); iso_Track++){
	  TVector3 vguess(0,0,0);
	  MuonsTrackParticles_BestPair.push_back(Ntp->IsolationTrack_TrackParticle(iso_Track)); // adding a track to best muon pair for vertex fit
	  Chi2VertexFitter FitterBestMuonPair_plusTrack(MuonsTrackParticles_BestPair,vguess);
	  if(FitterBestMuonPair_plusTrack.Fit() && FitterBestMuonPair_plusTrack.ChiSquare() < 10e+6){
	    //	    std::cout<<"all chi2   "<<FitterBestMuonPair_plusTrack.ChiSquare() << std::endl;
	    if(FitterBestMuonPair_plusTrack.ChiSquare() < tempChi2){
	      tempChi2 = FitterBestMuonPair_plusTrack.ChiSquare();
	      BestIsoTrackAddedToMuonPair = iso_Track;
	    }
	  }
	  MuonsTrackParticles_BestPair.pop_back();
	}
	if(Fitter.Fit()){ // if also three mu vertex is valid, compare chi2's
	  //	  std::cout<< "tempChi2  "<< tempChi2 << "  3mu chi2   "  << Fitter.ChiSquare() <<std::endl;

	  if(tempChi2< Fitter.ChiSquare()) BbkgVertexHypothesis.at(t).Fill(0.,1);
	  if(tempChi2> Fitter.ChiSquare()) BbkgVertexHypothesis.at(t).Fill(1.,1);

	}
      }
    }



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


    double VertexQuality_OS_SS1;    
    double VertexQuality_OS_SS2;    


    //    std::cout<<"lets check "<< std::endl;
    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,0);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,2);
	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,2);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,0);

	  }
      }


    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,0);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,1);
	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,1);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,0);

	  }
      }


    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,2);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,1);

	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,1);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,2);
	  }
      }



    float VertexQualitySortedMass1,VertexQualitySortedMass2;
    float BetterPhiVertex, WorsePhiVertex;
    if(VertexQuality_OS_SS1 < VertexQuality_OS_SS2){
      VertexQualitySortedMass1 = (MuonOS+MuonSS1).M();
      VertexQualitySortedMass2 = (MuonOS+MuonSS2).M();

      BetterPhiVertex=VertexQuality_OS_SS1;
      WorsePhiVertex=VertexQuality_OS_SS2;
    }else{
      VertexQualitySortedMass1 = (MuonOS+MuonSS2).M();
      VertexQualitySortedMass2 = (MuonOS+MuonSS1).M();
      BetterPhiVertex=VertexQuality_OS_SS2;
      WorsePhiVertex=VertexQuality_OS_SS1;

    }
    BetterMuMuVertex.at(t).Fill(BetterPhiVertex,1);
    WorseMuMuVertex.at(t).Fill(WorsePhiVertex,1);


    std::vector<unsigned int> Indices;
    Indices.push_back(ss1_mu_idx);
    Indices.push_back(ss2_mu_idx);



    unsigned int SS1RandomIndex(0);
    unsigned int SS2RandomIndex(0);


    float random_muon_index = rndm.Uniform();
    if(random_muon_index >= 0.5 ){SS1RandomIndex =  Indices.at(0); SS2RandomIndex = Indices.at(1) ; }
    if(random_muon_index <  0.5 ){SS1RandomIndex =  Indices.at(1); SS2RandomIndex = Indices.at(0) ; }

    TLorentzVector MuonLV_RandomSS1 = Ntp->Muon_P4(SS1RandomIndex);
    TLorentzVector MuonLV_RandomSS2 = Ntp->Muon_P4(SS2RandomIndex);
    
      
    float dRSortedMassPair1,dRSortedMassPair2;
    unsigned int ss1_mu_idx_dr, ss2_mu_idx_dr;

    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      ss1_mu_idx_dr=ss2_mu_idx;
      ss2_mu_idx_dr=ss1_mu_idx;
      dRSortedMassPair1 = (MuonOS+MuonSS2).M();
      dRSortedMassPair2 = (MuonOS+MuonSS1).M();
    }else{
      ss1_mu_idx_dr=ss1_mu_idx;
      ss2_mu_idx_dr=ss2_mu_idx;
      dRSortedMassPair1 = (MuonOS+MuonSS1).M();
      dRSortedMassPair2 = (MuonOS+MuonSS2).M();
    }

    TLorentzVector KaonOS(0,0,0,0);  KaonOS.SetXYZM(MuonOS.Px(),MuonOS.Py(),MuonOS.Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS1(0,0,0,0); KaonSS1.SetXYZM(Ntp->Muon_P4(ss1_mu_idx_dr).Px(),Ntp->Muon_P4(ss1_mu_idx_dr).Py(),Ntp->Muon_P4(ss1_mu_idx_dr).Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS2(0,0,0,0); KaonSS2.SetXYZM(Ntp->Muon_P4(ss2_mu_idx_dr).Px(),Ntp->Muon_P4(ss2_mu_idx_dr).Py(),Ntp->Muon_P4(ss2_mu_idx_dr).Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS1pTsort(0,0,0,0); KaonSS1pTsort.SetXYZM(MuonSS1.Px(),MuonSS1.Py(),MuonSS1.Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS2pTsort(0,0,0,0); KaonSS2pTsort.SetXYZM(MuonSS2.Px(),MuonSS2.Py(),MuonSS2.Pz(),PDG_Var::Kp_mass());
    for(unsigned int iIsoTrack=0; iIsoTrack < Ntp->NIsolationTrack(signal_idx); iIsoTrack++){
      if(Ntp->Muon_charge(os_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_OS.at(t).Fill((KaonOS + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
      if(Ntp->Muon_charge(ss1_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_SS1.at(t).Fill((KaonSS1pTsort + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
      if(Ntp->Muon_charge(ss2_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_SS2.at(t).Fill((KaonSS2pTsort + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
    }



    
    
    float Mass_osss1 = (Ntp->Muon_P4(os_mu_idx)+Ntp->Muon_P4(ss1_mu_idx)).M();
    float Mass_osss2 = (Ntp->Muon_P4(os_mu_idx)+Ntp->Muon_P4(ss2_mu_idx)).M();
    
    float CloserToPhiMassPair  = fabs(Mass_osss1-PDG_Var::Phi_mass())  < fabs(Mass_osss2-PDG_Var::Phi_mass()) ? Mass_osss1 : Mass_osss2;
    
        
    
    bool phiVeto(false);
    
    bool CrossVeto1(true);
    bool CrossVeto2(true);
    bool CrossVeto3(true);
    bool CrossVeto(true);
    double    m12v = (MuonOS+MuonSS1).M();
    double    m13v = (MuonOS+MuonSS2).M();
    
    
    
    //    if((dRSortedMassPair1 < phiVetoCut1 || dRSortedMassPair1 > phiVetoCut2 ) && (dRSortedMassPair2 < 0.65 || dRSortedMassPair2 > 1.6) )CrossVeto=true;
    if((dRSortedMassPair1 > phiVetoCut1 &&  dRSortedMassPair2 > 0.65)  &&  (dRSortedMassPair1 < phiVetoCut2 &&  dRSortedMassPair2 < 1.6))CrossVeto1=false;
    if((dRSortedMassPair2 > phiVetoCut1 &&  dRSortedMassPair1 > 0.2)  &&  (dRSortedMassPair2 < phiVetoCut2 &&  dRSortedMassPair1 < 1.4))CrossVeto2=false;
    if((dRSortedMassPair1 > rmgCutVeto1 &&  dRSortedMassPair2 > 0.95)  &&  (dRSortedMassPair1 < rmgCutVeto2 &&  dRSortedMassPair2 < 1.5))CrossVeto3=false;





    //*** 3-mu mass after the KF vertex constrain

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);






    VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(signal_idx),w);
    FLSignificance.at(t).Fill(sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(signal_idx),Ntp->Vertex_PrimaryVertex_Covariance(signal_idx),
								   Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx))),w);
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    SVPVTauDirAngle.at(t).Fill(SVPV.Angle(TauLV.Vect()),w);

    
    //***  define the mva varables used for evaluation of BDT weights for selection
    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      var_mass12_dRsorting = (MuonOS+MuonSS2).M();
      var_mass13_drSorting = (MuonOS+MuonSS1).M();
    }else{
      var_mass12_dRsorting = (MuonOS+MuonSS1).M();
      var_mass13_drSorting = (MuonOS+MuonSS2).M();
    }
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


    var_tauMass=TauRefitLV.M();

    // ** Isolation mass spectra
    double Chi2IsoTrackVertexToMuon3(999.);
    TLorentzVector IsoTrack_P4_closestToMu3(0,0,0,0);

    double Chi2IsoTrackVertexToMuon2(999.);
    TLorentzVector IsoTrack_P4_closestToMu2(0,0,0,0);

    double Chi2IsoTrackVertexToMuon1(999.);
    TLorentzVector IsoTrack_P4_closestToMu1(0,0,0,0);

    for(int i =0; i< Ntp->NIsolationTrack(signal_idx); i++){

      if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(2)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon3)
	       {
		 Chi2IsoTrackVertexToMuon3= Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu3 = Ntp->IsolationTrack_p4(signal_idx,i);
	       }
	   }


	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(1)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon2)
	       {
		 Chi2IsoTrackVertexToMuon2= Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu2 = Ntp->IsolationTrack_p4(signal_idx,i);

	       }

	   }


	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(0)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon1)
	       {
		 Chi2IsoTrackVertexToMuon1= Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu1 = Ntp->IsolationTrack_p4(signal_idx,i);

	       }

	   }
    }




    IsoTrack_P4_closestToMu3.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
				       IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
				       IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.493677*0.493677));

    TLorentzVector IsoTrack_P4_closestToMu3_PiMass = IsoTrack_P4_closestToMu3;
    IsoTrack_P4_closestToMu3_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
					      IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
					      IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.135*0.135));

    TLorentzVector IsoTrack_P4_closestToMu3_MuMass = IsoTrack_P4_closestToMu3;
    IsoTrack_P4_closestToMu3_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
					      IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
					      IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.106*0.106));


    TLorentzVector Mu3_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
    
    Mu3_WithKMass.SetE(sqrt(Mu3_WithKMass.Px()*Mu3_WithKMass.Px() + 
			    Mu3_WithKMass.Py()*Mu3_WithKMass.Py() + 
			    Mu3_WithKMass.Pz()*Mu3_WithKMass.Pz() + 0.493677*0.493677));
    
  
    
    var_IsoPhiKKMass_Mu3 = (IsoTrack_P4_closestToMu3+Mu3_WithKMass).M();
    var_IsoKStarMass_Mu3 = (IsoTrack_P4_closestToMu3_PiMass+Mu3_WithKMass).M();
    var_IsoMuMuMass_Mu3 = (IsoTrack_P4_closestToMu3_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2))).M();
    
    
    
    
    IsoTrack_P4_closestToMu2.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
				       IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
				       IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.493677*0.493677));
    TLorentzVector Mu2_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1));
    
    Mu2_WithKMass.SetE(sqrt(Mu2_WithKMass.Px()*Mu2_WithKMass.Px() + 
			    Mu2_WithKMass.Py()*Mu2_WithKMass.Py() + 
			    Mu2_WithKMass.Pz()*Mu2_WithKMass.Pz() + 0.493677*0.493677));
    
    TLorentzVector IsoTrack_P4_closestToMu2_PiMass = IsoTrack_P4_closestToMu2;
    IsoTrack_P4_closestToMu2_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
					      IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
					      IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.135*0.135));
    
    TLorentzVector IsoTrack_P4_closestToMu2_MuMass = IsoTrack_P4_closestToMu2;
    IsoTrack_P4_closestToMu2_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
					      IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
					      IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.106*0.106));
    
   
    
    var_IsoPhiKKMass_Mu2 = (IsoTrack_P4_closestToMu2+Mu2_WithKMass).M();
    var_IsoKStarMass_Mu2 = (IsoTrack_P4_closestToMu2_PiMass+Mu2_WithKMass).M();
    var_IsoMuMuMass_Mu2 = (IsoTrack_P4_closestToMu2_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))).M();
    


    IsoTrack_P4_closestToMu1.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
				       IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
				       IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.493677*0.493677));
    
    TLorentzVector IsoTrack_P4_closestToMu1_PiMass = IsoTrack_P4_closestToMu1;
    IsoTrack_P4_closestToMu1_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
					      IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
					      IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.135*0.135));
    
    TLorentzVector IsoTrack_P4_closestToMu1_MuMass = IsoTrack_P4_closestToMu1;
    IsoTrack_P4_closestToMu1_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
					      IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
					      IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.106*0.106));
    
    
    
    TLorentzVector Mu1_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0));
    
    Mu1_WithKMass.SetE(sqrt(Mu1_WithKMass.Px()*Mu1_WithKMass.Px() + 
			    Mu1_WithKMass.Py()*Mu1_WithKMass.Py() + 
			    Mu1_WithKMass.Pz()*Mu1_WithKMass.Pz() + 0.493677*0.493677));
    
   
    
    var_IsoPhiKKMass_Mu1 = (IsoTrack_P4_closestToMu1+Mu1_WithKMass).M();
    var_IsoKStarMass_Mu1 = (IsoTrack_P4_closestToMu1_PiMass+Mu1_WithKMass).M();
    var_IsoMuMuMass_Mu1 = (IsoTrack_P4_closestToMu1_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))).M();
    



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


    var_Muon1DetID = Muon1DetID; 
    var_Muon2DetID = Muon2DetID;
    var_Muon3DetID = Muon3DetID;

    bool KeepSignalRegionForMC(true);
    if(id!=1) KeepSignalRegionForMC = true;
    if(id==1 && (TauRefitLV.M() > tauMinSideBand_ && TauRefitLV.M() < tauMinMass_) or (TauRefitLV.M() > tauMaxMass_ && TauRefitLV.M() < tauMaxSideBand_) ) KeepSignalRegionForMC=true;


    //  define also phi veto here

    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_){


      if((dRSortedMassPair1<0.994 ||  dRSortedMassPair1 > 1.044) && (dRSortedMassPair2<0.994 || dRSortedMassPair2> 1.044) ) phiVeto = true;


    } else if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_){

      if((dRSortedMassPair1<0.985 ||  dRSortedMassPair1 > 1.053) && (dRSortedMassPair2<0.985 || dRSortedMassPair2> 1.053) ) phiVeto = true;



    }else if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_){

      if((dRSortedMassPair1<0.974 ||  dRSortedMassPair1 > 1.064) && (dRSortedMassPair2<0.974 || dRSortedMassPair2> 1.064) ) phiVeto = true;


    }



    double dRSortedMass;

    if(phiVeto){
      //**** Categories ABC1
      if(  (Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_ && readerA->EvaluateMVA("BDT") > mvaA2_)  or 
	   (Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_ && readerB->EvaluateMVA("BDT") > mvaB2_)  or 
	   (Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_ && readerC->EvaluateMVA("BDT") > mvaC2_) )
	{
	  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1); 
	  
	}


   //**** Categories ABC2
    if(  (Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_ && readerA->EvaluateMVA("BDT") > mvaA1_ && readerA->EvaluateMVA("BDT") < mvaA2_)  or 
	 (Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_ && 
	  readerB->EvaluateMVA("BDT") > mvaB1_ && readerB->EvaluateMVA("BDT") < mvaB2_)  or 
	 (Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_ && readerC->EvaluateMVA("BDT") > mvaC1_ && readerC->EvaluateMVA("BDT")< mvaC2_) )
	{
	  
	  TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    

	  
	}
    }
    /*
    if(id == 40 or id == 60 or id == 90 ){// or id == 40){
      std::cout<<"-------------- All categoris ----------------"<< std::endl;
      std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Ntp->printMCDecayChainOfEvent(true, true, true, true);
      std::cout<< "\n\n\n\n\n\n";
      }*/


    //*** below are basic  purity and resolution plots
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){
	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
	
	TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	if(Fitter.Fit())	    TauMassResolutionHelixRefit.at(t).Fill((MotherParticle.LV().M() - MCTauLV.M())/MCTauLV.M(),1);


	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }
    
    
    
  }
}



void  VertexFits::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis


  //*** extra actions
    if(mode == RECONSTRUCT){


      MassPairFit.at(1).Scale( MassPairFit.at(0).Integral()/ MassPairFit.at(1).Integral());
      Chi2SquarePairFit.at(1).Scale( Chi2SquarePairFit.at(0).Integral()/ Chi2SquarePairFit.at(1).Integral());
    DistanceToSignalVertexPairFit.at(1).Scale( DistanceToSignalVertexPairFit.at(0).Integral()/ DistanceToSignalVertexPairFit.at(1).Integral());


    MassTripleFit.at(1).Scale( MassTripleFit.at(0).Integral()/ MassTripleFit.at(1).Integral());
    Chi2SquareTripleFit.at(1).Scale( Chi2SquareTripleFit.at(0).Integral()/ Chi2SquareTripleFit.at(1).Integral());
    DistanceToSignalVertexTripleFit.at(1).Scale( DistanceToSignalVertexTripleFit.at(0).Integral()/ DistanceToSignalVertexTripleFit.at(1).Integral());

    BbkgVertexHypothesis.at(1).Scale( BbkgVertexHypothesis.at(0).Integral()/BbkgVertexHypothesis.at(1).Integral());



    DeltaR_D_SV_3tracks.at(1).Scale( DeltaR_D_SV_3tracks.at(0).Integral()/DeltaR_D_SV_3tracks.at(1).Integral());
    DeltaR_D_SV_2tracks.at(1).Scale( DeltaR_D_SV_2tracks.at(0).Integral()/DeltaR_D_SV_2tracks.at(1).Integral());


    Angle_D_SV_3tracks.at(1).Scale( Angle_D_SV_3tracks.at(0).Integral()/Angle_D_SV_3tracks.at(1).Integral());
    Angle_D_SV_2tracks.at(1).Scale( Angle_D_SV_2tracks.at(0).Integral()/Angle_D_SV_2tracks.at(1).Integral());

    PairMassHighestPlusBestVertexMass.at(1).Scale( PairMassHighestPlusBestVertexMass.at(0).Integral()/ PairMassHighestPlusBestVertexMass.at(1).Integral());

    PairMassHighestPlusSubHighest.at(1).Scale( PairMassHighestPlusSubHighest.at(0).Integral()/ PairMassHighestPlusSubHighest.at(1).Integral());
    TripleMassHighestPlusBestVertexMass.at(1).Scale( TripleMassHighestPlusBestVertexMass.at(0).Integral()/ TripleMassHighestPlusBestVertexMass.at(1).Integral());



    

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


    //    if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    //    ScaleAllHistOfType(3,scale*scaleB0Tau);
    }
  Selection::Finish();
}





