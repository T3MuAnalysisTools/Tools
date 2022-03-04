#include "DebugFit.h"
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

DebugFit::DebugFit(TString Name_, TString id_):
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



DebugFit::~DebugFit(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  DebugFit::Configure(){


  //  This mini tree is for limit extraction
  gErrorIgnoreLevel = kFatal;





  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  //*** defined the bdt reader for event selection; readerA- category A, readerB - category B ...
   readerA = new TMVA::Reader( "!Color:!Silent" );
  readerA->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerA->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerA->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerA->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerA->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerA->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerA->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerA->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);



  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_24_A/weights/TMVAClassification_BDT.weights.xml"); 




  readerB = new TMVA::Reader( "!Color:!Silent" );
  readerB->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerB->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerB->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerB->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerB->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerB->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerB->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerB->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);


  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerB->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_24_B/weights/TMVAClassification_BDT.weights.xml"); 


  readerC = new TMVA::Reader( "!Color:!Silent" );
  readerC->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerC->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerC->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerC->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerC->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerC->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerC->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerC->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);


  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerC->BookMVA( "BDT", "/afs/cern.ch/work/m/mmadhu/Analysis/workdirMakeMVAWithFilterMay_09_2021/Code/CommonUtils/IterativeTrain/output_24_C/weights/TMVAClassification_BDT.weights.xml");





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
    



  TwoTracksFitOk=HConfig.GetTH1D(Name+"_TwoTracksFitOk","TwoTracksFitOk",2,-0.5,1.5,"Fit is Ok ","Events");
  SignalVertexFitOk=HConfig.GetTH1D(Name+"_SignalVertexFitOk","SignalVertexFitOk",2,-0.5,1.5,"Signal Vertex Fit is Ok","Events");
  SignalFitVertexChi2=HConfig.GetTH1D(Name+"_SignalFitVertexChi2","SignalFitVertexChi2",50,-0.5,49.5,"Helix Fitter Vertex #chi^{2} ","Events");
  SignalFitVertexChi2VsKalman=HConfig.GetTH2D(Name+"_SignalFitVertexChi2VsKalman","SignalFitVertexChi2VsKalman",50,-0.5,49.5,50,-0.5,49.5,"Helix Fitter Vertex #chi^{2}","Kalman Vertex Fit #chi^{2}");
  TwoLargestPtIsolationTracks=HConfig.GetTH1D(Name+"_TwoLargestPtIsolationTracks","TwoLargestPtIsolationTracks",50,-0.5,49.5,"Two largest pT isolation tracks  #chi^{2} ","Events");

  HelixKalmanDeltaX=HConfig.GetTH1D(Name+"_HelixKalmanDeltaX","HelixKalmanDeltaX",50,-0.1,0.1,"3#mu vertex #Delta X, cm (Helix - Kalman) ","Events");
  HelixKalmanDeltaY=HConfig.GetTH1D(Name+"_HelixKalmanDeltaY","HelixKalmanDeltaY",50,-0.1,0.1,"3#mu vertex #Delta Y, cm (Helix - Kalman) ","Events");
  HelixKalmanDeltaZ=HConfig.GetTH1D(Name+"_HelixKalmanDeltaZ","HelixKalmanDeltaZ",50,-0.1,0.1,"3#mu vertex #Delta Z, cm (Helix - Kalman) ","Events");

  ReffitedMotherMass = HConfig.GetTH1D(Name+"_ReffitedMotherMass","ReffitedMotherMass",50,-0.1,0.1,"#Delta M_{3#mu}, GeV (Helix - Kalman)","Events");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  DebugFit::Store_ExtraDist(){ 

  Extradist1d.push_back(&TwoTracksFitOk);
  Extradist1d.push_back(&SignalVertexFitOk);
  Extradist1d.push_back(&SignalFitVertexChi2);
  Extradist2d.push_back(&SignalFitVertexChi2VsKalman);
  Extradist1d.push_back(&TwoLargestPtIsolationTracks);

  Extradist1d.push_back(&ReffitedMotherMass);

  Extradist1d.push_back(&HelixKalmanDeltaX);
  Extradist1d.push_back(&HelixKalmanDeltaY);
  Extradist1d.push_back(&HelixKalmanDeltaZ);





}


void  DebugFit::doEvent(){ 

  
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

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    std::vector<TrackParticle> MuonsTrackParticles;

    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_1));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_2));
    MuonsTrackParticles.push_back(Ntp->Muon_TrackParticle(Muon_index_3));



    TVector3 vg(0.1,0.1,0.1);
    Chi2VertexFitter  SignalVertex(MuonsTrackParticles,vg);
    bool SignalFitOk  =    SignalVertex.Fit();
    SignalVertexFitOk.at(t).Fill(SignalFitOk,1);
    SignalFitVertexChi2.at(t).Fill(SignalVertex.ChiSquare(),1);
    SignalFitVertexChi2VsKalman.at(t).Fill(SignalVertex.ChiSquare(),Ntp->Vertex_signal_KF_Chi2(signal_idx) );

    HelixKalmanDeltaX.at(t).Fill(SignalVertex.GetVertex().X() - Ntp->Vertex_Signal_KF_pos(signal_idx).X(),1);
    HelixKalmanDeltaY.at(t).Fill(SignalVertex.GetVertex().Y() - Ntp->Vertex_Signal_KF_pos(signal_idx).Y(),1);
    HelixKalmanDeltaZ.at(t).Fill(SignalVertex.GetVertex().Z() - Ntp->Vertex_Signal_KF_pos(signal_idx).Z(),1);

    TLorentzVector ReffitedTau = SignalVertex.GetMother(444).LV();


    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);  
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    ReffitedMotherMass.at(t).Fill(ReffitedTau.M() - TauLV.M(),1);
    std::vector<unsigned int> EtaSortedIndices;
    
    EtaSortedIndices.push_back(Muon_Eta_index_1);
    EtaSortedIndices.push_back(Muon_Eta_index_2);
    EtaSortedIndices.push_back(Muon_Eta_index_3);



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
	  
	if(Ntp->IsolationTrack_p4(signal_idx, i2).Pt() > MaxPT){
	  MaxPT = Ntp->IsolationTrack_p4(signal_idx, i2).Pt();
	  SubHighestPtTrack = i2;
	}
      }
    }
  
      
    // ----------------------------------   debugging fit -----------------

    //    std::cout<<" ========================  "<<std::endl;
    //    std::cout<<"SubHighestPtTrack   "<< SubHighestPtTrack << "   HighestPtTrack   " << HighestPtTrack << std::endl;

    if(SubHighestPtTrack!=-1 && HighestPtTrack!=-1){


	std::vector<TrackParticle> DebugFitTracks;
	TMatrixT<double> vala;
	TMatrixTSym<double> cova;
	TMatrixTSym<double> cova_inv;

	TMatrixTSym<double> cova_inflated;
	vala.ResizeTo(TrackParticle::NHelixPar,1);
	cova.ResizeTo(TrackParticle::NHelixPar,TrackParticle::NHelixPar);
	cova_inflated.ResizeTo(TrackParticle::NHelixPar,TrackParticle::NHelixPar);
	cova_inv.ResizeTo(TrackParticle::NHelixPar,TrackParticle::NHelixPar);

	double cova_det = cova.Determinant();


	for(unsigned int j=0; j<TrackParticle::NHelixPar;j++){
	  vala(j,0)=Ntp->IsolationTrack_TrackParticle(HighestPtTrack).Parameter(j);
	  for(unsigned int k=0; k<TrackParticle::NHelixPar;k++){
	    cova(j,k)=Ntp->IsolationTrack_TrackParticle(HighestPtTrack).Covariance(j,k);

	  }
	}





	//	std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  DebugFit  :  "<<std::endl; 
	//	std::cout<<"   cova  inv  and eigen values  "<< std::endl;
	//	cova.Print(); Ntp->EigenValues(cova).Print();
	//	TMatrixTSym<double>   Inverted1 = Ntp->InvertMatrix(cova);
	//std::cout<<" ///////////////////////////////////////  Inverted eigen values 1 "<< std::endl;Ntp->EigenValues(Inverted1).Print();
	//std::cout<<"cova_infl_SIMILARITY  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<=====================  "<< std::endl;
	//	Inverted1.SimilarityT(vala).Print();
	/*
	std::cout<<" cova inv  inflamed   and eigen values"<< std::endl;

	TMatrixTSym<double> cova_infl= 	Ntp->RegulariseCovariance(cova, 1.0);
	cova_infl.Print();   Ntp->EigenValues(cova_infl).Print();
	TMatrixTSym<double>   Inverted2 = Ntp->InvertMatrix(cova_infl);

	std::cout<<" ///////////////////////////////////////  Inverted eigen values 1 "<< std::endl;Ntp->EigenValues(Inverted1).Print();
	std::cout<<" \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  Inverted eigen values 2 "<< std::endl;Ntp->EigenValues(Inverted2).Print();





	std::cout<<"cova_SIMILARITY  =========================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<< std::endl;
	Inverted1.SimilarityT(vala).Print();


	std::cout<<"cova_infl_SIMILARITY  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<=====================  "<< std::endl;
	Inverted2.SimilarityT(vala).Print();
	*/


	DebugFitTracks.push_back(Ntp->IsolationTrack_TrackParticle(HighestPtTrack));
	DebugFitTracks.push_back(Ntp->IsolationTrack_TrackParticle(SubHighestPtTrack));
	//	DebugFitTracks.push_back(Ntp->IsolationTrack_TrackParticle(4));
	TVector3 vguess(0.1,0.1,0.1);
	Chi2VertexFitter  DebugFitter(DebugFitTracks,vguess);
	bool fitOk =   DebugFitter.Fit();
	TwoTracksFitOk.at(t).Fill(fitOk,1);
	TwoLargestPtIsolationTracks.at(t).Fill(DebugFitter.ChiSquare(),1);
	//	std::cout<<"  Debugger chi2    "<< DebugFitter.ChiSquare() << " isFit ok    "  << fitOk <<std::endl;
	//	std::cout<<"  Vertex Position    "<< std::endl;
	//      DebugFitter.GetVertex().Print();
      }
    }


    //--------------------------------------------------------------------
}
  



void  DebugFit::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis


  //*** extra actions
   Selection::Finish();
}





