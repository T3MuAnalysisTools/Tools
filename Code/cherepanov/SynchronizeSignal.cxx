#include "SynchronizeSignal.h"
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

SynchronizeSignal::SynchronizeSignal(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.62),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

SynchronizeSignal::~SynchronizeSignal(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SynchronizeSignal::Configure(){

  Sync_tree= new TTree("tree","tree");

  Sync_tree->Branch("evt",&evt);
  Sync_tree->Branch("run",&run);
  Sync_tree->Branch("lumi",&lumi);

  Sync_tree->Branch("sync_pt_1",&sync_pt_1);
  Sync_tree->Branch("sync_pt_2",&sync_pt_2);
  Sync_tree->Branch("sync_pt_3",&sync_pt_3);
  Sync_tree->Branch("sync_chi2",&sync_chi2);



  /*  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("bdt",&bdt);
  T3MMiniTree->Branch("category",&category);
  T3MMiniTree->Branch("rapidity",&rapidity);
  T3MMiniTree->Branch("LumiScale",&LumiScale);
  */
  readerA = new TMVA::Reader( "!Color:!Silent" );

  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonUtils/tmva/dataset/weights/";

  readerA->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerA->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerA->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerA->AddVariable( "var_sumMuTrkKinkChi2", &var_sumMuTrkKinkChi2 );
  readerA->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerA->AddVariable( "var_MinMIPLikelihood", &var_MinMIPLikelihood );
  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->AddVariable( "var_maxdca", &var_maxdca );
  readerA->AddVariable( "var_RelPt_Mu1Tau", &var_RelPt_Mu1Tau );
  readerA->AddVariable( "var_MaxD0Significance", &var_MaxD0Significance );
  readerA->AddVariable( "var_IsolationMinDist", &var_IsolationMinDist);


  readerA->BookMVA( "BDT", basedir+"Run2017Classification_unbiased_A_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerB = new TMVA::Reader( "!Color:!Silent" );

  readerB->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerB->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerB->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerB->AddVariable( "var_sumMuTrkKinkChi2", &var_sumMuTrkKinkChi2 );
  readerB->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerB->AddVariable( "var_MinMIPLikelihood", &var_MinMIPLikelihood );
  readerB->AddSpectator("var_tauMass",&var_tauMass);
  readerB->AddVariable( "var_maxdca", &var_maxdca );
  readerB->AddVariable( "var_RelPt_Mu1Tau", &var_RelPt_Mu1Tau );
  readerB->AddVariable( "var_MaxD0Significance", &var_MaxD0Significance );
  readerB->AddVariable( "var_IsolationMinDist", &var_IsolationMinDist);
  readerB->BookMVA( "BDT", basedir+"Run2017Classification_unbiased_B_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerC = new TMVA::Reader( "!Color:!Silent" );
  readerC->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerC->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerC->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerC->AddVariable( "var_sumMuTrkKinkChi2", &var_sumMuTrkKinkChi2 );
  readerC->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerC->AddVariable( "var_MinMIPLikelihood", &var_MinMIPLikelihood );
  readerC->AddSpectator("var_tauMass",&var_tauMass);
  readerC->AddVariable( "var_maxdca", &var_maxdca );
  readerC->AddVariable( "var_RelPt_Mu1Tau", &var_RelPt_Mu1Tau );
  readerC->AddVariable( "var_MaxD0Significance", &var_MaxD0Significance );
  readerC->AddVariable( "var_IsolationMinDist", &var_IsolationMinDist);


  readerC->BookMVA( "BDT", basedir+"Run2017Classification_unbiased_C_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles


  readerBarrel = new TMVA::Reader( "!Color:!Silent" );
  readerBarrel->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerBarrel->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerBarrel->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerBarrel->AddVariable( "var_sumMuTrkKinkChi2", &var_sumMuTrkKinkChi2 );
  readerBarrel->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerBarrel->AddVariable( "var_MinMIPLikelihood", &var_MinMIPLikelihood );
  readerBarrel->AddSpectator("var_tauMass",&var_tauMass);
  readerBarrel->AddVariable( "var_maxdca", &var_maxdca );
  readerBarrel->AddVariable( "var_RelPt_Mu1Tau", &var_RelPt_Mu1Tau );
  readerBarrel->AddVariable( "var_MaxD0Significance", &var_MaxD0Significance );
  readerBarrel->AddVariable( "var_IsolationMinDist", &var_IsolationMinDist);


  readerBarrel->BookMVA( "BDT", basedir+"Run2017Classification_unbiased_Barrel_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles

  readerEndcap = new TMVA::Reader( "!Color:!Silent" );
  readerEndcap->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerEndcap->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerEndcap->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerEndcap->AddVariable( "var_sumMuTrkKinkChi2", &var_sumMuTrkKinkChi2 );
  readerEndcap->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerEndcap->AddVariable( "var_MinMIPLikelihood", &var_MinMIPLikelihood );
  readerEndcap->AddSpectator("var_tauMass",&var_tauMass);
  readerEndcap->AddVariable( "var_maxdca", &var_maxdca );
  readerEndcap->AddVariable( "var_RelPt_Mu1Tau", &var_RelPt_Mu1Tau );
  readerEndcap->AddVariable( "var_MaxD0Significance", &var_MaxD0Significance );
  readerEndcap->AddVariable( "var_IsolationMinDist", &var_IsolationMinDist);


  readerEndcap->BookMVA( "BDT", basedir+"Run2017Classification_unbiased_Endcap_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



   // dataloader->AddVariable("var_vertexKFChi2","Variable vertexKFChi2","units", 'F' );
   // dataloader->AddVariable("var_svpvTauAngle","Variable svpvTauAngle","units", 'F' );
   // dataloader->AddVariable("var_flightLenSig","Variable flightLenSig","units", 'F' );
   // dataloader->AddVariable("var_sumMuTrkKinkChi2","Variable sumMuTrkKinkChi2","units", 'F' );
   // dataloader->AddVariable("var_segCompMuMin","Variable segCompMuMin","units", 'F' );
   // dataloader->AddVariable("var_MinMIPLikelihood","Variable MinMIPLikelihood","units", 'F' );
   // dataloader->AddSpectator("var_tauMass","Variable tauMass","units", 'F' );
   // dataloader->AddVariable("var_maxdca","Variable maxdca","units", 'F' );
   // dataloader->AddVariable("var_RelPt_Mu1Tau","Variable RelPt_Mu1Tau","units", 'F' );
   // dataloader->AddVariable("var_MaxD0Significance","Variable MaxD0Significance","units", 'F' );
   // dataloader->AddVariable("var_IsolationMinDist","Variable IsolationMinDist","units", 'F' );




  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==VertChi2)           cut.at(VertChi2)=15.0;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=2.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch1)      cut.at(TriggerMatch1)=0.03;
    if(i==TriggerMatch2)      cut.at(TriggerMatch2)=0.03;
    if(i==TriggerMatch3)      cut.at(TriggerMatch3)=0.03;
    if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==TriggerOk){
      title.at(i)="Pass HLT and L1";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="Signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==VertChi2){
      title.at(i)="Secondary Vertex $\\chi^{2} <$ 15";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Chi squared of triple vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VertChi2_",htitle,50,0,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VertChi2_",htitle,50,0,25,hlabel,"Events"));
    }
    
    else if(i==Mu1PtCut){
      title.at(i)="$p_{T}(\\mu_{1}) >$ 2.0 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Muon1 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
    }
    else if(i==Mu2PtCut){
      title.at(i)="$p_{T}(\\mu_{2}) >$ 2.0 GeV";
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
      title.at(i)="All muons are global";
      hlabel="gl,gl,gl";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PhiVeto){
      title.at(i)="0.997 GeV $ < $  $\\phi$ mass  $ < $ 1.041 GeV ";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Phi mass Veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto){
      title.at(i)="0.765 GeV $ < $  $\\omega$ mass $ < $ 0.799 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Rho mass veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
    }
    else if(i==TriggerMatch1){
      title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{1}} <$ 0.03";
      hlabel="Sum of dR_{reco-trigger}";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch1_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch1_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==TriggerMatch2){
      title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{2}} <$ 0.03";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch2_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch2_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==TriggerMatch3){
      title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{3}} <$ 0.03";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch3_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch3_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==TauMassCut){
      title.at(i)="$1.62 < \\tau $ mass $ < 2.0$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,60,1.4,2.2,hlabel,"Events"));
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

  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");
  TauMass_all =HConfig.GetTH1D(Name+"_TauMass_all","3#mu  mass",60,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMass_allVsBDTA=HConfig.GetTH2D(Name+"_TauMass_allVsBDTA","3#mu mass vs BDTa",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTB=HConfig.GetTH2D(Name+"_TauMass_allVsBDTB","3#mu mass vs BDTb",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTC=HConfig.GetTH2D(Name+"_TauMass_allVsBDTC","3#mu mass vs BDTc",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTBarrel=HConfig.GetTH2D(Name+"_TauMass_allVsBDTBarrel","3#mu mass vs BDTBarre;",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTEndcap=HConfig.GetTH2D(Name+"_TauMass_allVsBDTEndcap","3#mu mass vs BDTEndcap",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");

  EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

  TauMassA1 =HConfig.GetTH1D(Name+"_TauMassA1","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitA1 =HConfig.GetTH1D(Name+"_TauMassRefitA1","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassB1 =HConfig.GetTH1D(Name+"_TauMassB1","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitB1 =HConfig.GetTH1D(Name+"_TauMassRefitB1","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassC1 =HConfig.GetTH1D(Name+"_TauMassC1","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitC1 =HConfig.GetTH1D(Name+"_TauMassRefitC1","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassA2 =HConfig.GetTH1D(Name+"_TauMassA2","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitA2 =HConfig.GetTH1D(Name+"_TauMassRefitA2","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassB2 =HConfig.GetTH1D(Name+"_TauMassB2","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitB2 =HConfig.GetTH1D(Name+"_TauMassRefitB2","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassC2 =HConfig.GetTH1D(Name+"_TauMassC2","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitC2 =HConfig.GetTH1D(Name+"_TauMassRefitC2","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");

  TauMassBarrel1 =HConfig.GetTH1D(Name+"_TauMassBarrel1","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitBarrel1 =HConfig.GetTH1D(Name+"_TauMassRefitBarrel1","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");

  TauMassEndcap1 =HConfig.GetTH1D(Name+"_TauMassEndcap1","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitEndcap1 =HConfig.GetTH1D(Name+"_TauMassRefitEndcap1","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  TauMassBarrel2 =HConfig.GetTH1D(Name+"_TauMassBarrel2","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitBarrel2 =HConfig.GetTH1D(Name+"_TauMassRefitBarrel2","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");

  TauMassEndcap2 =HConfig.GetTH1D(Name+"_TauMassEndcap2","#tau lepton mass",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitEndcap2 =HConfig.GetTH1D(Name+"_TauMassRefitEndcap2","Refit #tau lepton mass",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");


  EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");


  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,3.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");

  BDTOutputA = HConfig.GetTH1D(Name+"_BDTOutputA","BDTOutputA",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputB = HConfig.GetTH1D(Name+"_BDTOutputB","BDTOutputB",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputC = HConfig.GetTH1D(Name+"_BDTOutputC","BDTOutputC",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputBarrel = HConfig.GetTH1D(Name+"_BDTOutputBarrel","BDTOutputBarrel",50,-0.4,0.4,"BDT Output","Events");
  BDTOutputEndcap = HConfig.GetTH1D(Name+"_BDTOutputEndcap","BDTOutputENdcap",50,-0.4,0.4,"BDT Output","Events");

  L1Triggers=HConfig.GetTH2D(Name+"_L1Triggers","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  L1TriggersB=HConfig.GetTH2D(Name+"_L1TriggersB","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersC=HConfig.GetTH2D(Name+"_L1TriggersC","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersD=HConfig.GetTH2D(Name+"_L1TriggersD","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersE=HConfig.GetTH2D(Name+"_L1TriggersE","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1TriggersF=HConfig.GetTH2D(Name+"_L1TriggersF","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  SynchronizeSignal::Store_ExtraDist(){ 


  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  //Extradist1d.push_back(&Muon1isGlob);
  //Extradist1d.push_back(&Muon2isGlob);
  //Extradist1d.push_back(&Muon3isGlob);


  //Extradist1d.push_back(&Muon1isStand);
  //Extradist1d.push_back(&Muon2isStand);
  //Extradist1d.push_back(&Muon3isStand);
    

  //Extradist1d.push_back(&Muon1isTrack);
  //Extradist1d.push_back(&Muon2isTrack);
  //Extradist1d.push_back(&Muon3isTrack);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  // Extradist1d.push_back(&TauP);


  Extradist1d.push_back(&TauMassA1);
  Extradist1d.push_back(&TauMassRefitA1);

  Extradist1d.push_back(&TauMassB1);
  Extradist1d.push_back(&TauMassRefitB1);

  Extradist1d.push_back(&TauMassC1);
  Extradist1d.push_back(&TauMassRefitC1);


  Extradist1d.push_back(&TauMassA2);
  Extradist1d.push_back(&TauMassRefitA2);

  Extradist1d.push_back(&TauMassB2);
  Extradist1d.push_back(&TauMassRefitB2);

  Extradist1d.push_back(&TauMassC2);
  Extradist1d.push_back(&TauMassRefitC2);

  Extradist1d.push_back(&TauMassBarrel1);
  Extradist1d.push_back(&TauMassRefitBarrel1);
  Extradist1d.push_back(&TauMassEndcap1);
  Extradist1d.push_back(&TauMassRefitEndcap1);

  Extradist1d.push_back(&TauMassBarrel2);
  Extradist1d.push_back(&TauMassRefitBarrel2);
  Extradist1d.push_back(&TauMassEndcap2);
  Extradist1d.push_back(&TauMassRefitEndcap2);

  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist2d.push_back(&TauMass_all_nophiVeto);
  Extradist1d.push_back(&TauMass_all);
  Extradist2d.push_back(&TauMass_allVsBDTA);
  Extradist2d.push_back(&TauMass_allVsBDTB);
  Extradist2d.push_back(&TauMass_allVsBDTC);
  Extradist2d.push_back(&TauMass_allVsBDTBarrel);
  Extradist2d.push_back(&TauMass_allVsBDTEndcap);
  Extradist2d.push_back(&EMR_tau_eta);



  Extradist1d.push_back(&SVPVTauDirAngle);

  //  Extradist1d.push_back(&Muon1DRToTruth);
  //  Extradist1d.push_back(&Muon2DRToTruth);
  //  Extradist1d.push_back(&Muon3DRToTruth);


  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);


  // Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&EventMassResolution_PtEtaPhi);

  //Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&NSignalCandidates);


  Extradist1d.push_back(&BDTOutputA);
  Extradist1d.push_back(&BDTOutputB);
  Extradist1d.push_back(&BDTOutputC);
  Extradist1d.push_back(&BDTOutputBarrel);
  Extradist1d.push_back(&BDTOutputEndcap);

  Extradist2d.push_back(&L1Triggers);

  Extradist2d.push_back(&L1TriggersB);
  Extradist2d.push_back(&L1TriggersC);
  Extradist2d.push_back(&L1TriggersD);
  Extradist2d.push_back(&L1TriggersE);
  Extradist2d.push_back(&L1TriggersF);

}


void  SynchronizeSignal::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  bool HLTOk(false);
  bool L1Ok(false);
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);

    if(id==1){
      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    }

    //    if(id==1 && Ntp->WhichEra(2017).Contains("RunF")){
    //      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    //    }	 

    if(id!=1){
      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    }
    
  }


  bool DoubleMuFired(0);
  bool TripleMuFired(0);
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17"))      TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }
  }



  if(DoubleMuFired  or TripleMuFired) L1Ok = true;
  value.at(TriggerOk)=(HLTOk and L1Ok);
  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));



  value.at(SignalCandidate)=0;
  unsigned int  signal_idx=0;
  value.at(TriggerMatch1)=0;
  value.at(TriggerMatch2)=0;
  value.at(TriggerMatch3)=0;

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
    value.at(VertChi2) = Ntp->Vertex_Signal_KF_Chi2(signal_idx);
    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
    //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
    //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);
    //
    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
    			Ntp->Muon_isGlobalMuon(mu3_pt_idx));
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
    TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);


    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

    value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

    //    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx)){
    //      value.at(TriggerMatch)+=i; 
    //    }
    value.at(TriggerMatch1) = Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx).at(0);
    value.at(TriggerMatch2) = Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx).at(1);
    value.at(TriggerMatch3) = Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx).at(2);
    //    std::cout<<"Trigger Match  "<< value.at(TriggerMatch1) << std::endl;
    //    value.at(TauMassCut) = TauLV.M();
    value.at(TauMassCut) = TauRefittedLV.M();
  }


  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  pass.at(VertChi2) = (value.at(VertChi2) <= cut.at(VertChi2));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   = (value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch1) = (value.at(TriggerMatch1)  <  cut.at(TriggerMatch1));
  pass.at(TriggerMatch2) = (value.at(TriggerMatch2)  <  cut.at(TriggerMatch2));
  pass.at(TriggerMatch3) = (value.at(TriggerMatch3)  <  cut.at(TriggerMatch3));
  pass.at(PhiVeto)       = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*0.011);
  pass.at(OmegaVeto)     = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*0.0085);
  pass.at(TauMassCut)    = ( (value.at(TauMassCut) > tauMinSideBand_)  &&   (value.at(TauMassCut)< tauMaxSideBand_ ));


    //  if(id!=1) pass.at(TauMassCut) = true;
    //  else  pass.at(TauMassCut) =true;//( (value.at(TauMassCut) < tauMinMass_)  ||   (value.at(TauMassCut)> tauMaxMass_ ));

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(PhiVeto);
  exclude_cuts.push_back(OmegaVeto);


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);

  if(status){


    //    std::cout<<" mindist   "<< Ntp->Isolation_MinDist(signal_idx) << std::endl;
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

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    sync_pt_1 =  Muon1LV.Pt();
    sync_pt_2 =  Muon2LV.Pt();
    sync_pt_3 =  Muon3LV.Pt();
    evt = Ntp->EventNumber();
    run = Ntp->RunNumber();
    lumi = Ntp->LuminosityBlock();
    sync_chi2 = Ntp->Vertex_signal_KF_Chi2(signal_idx);
    if(id==1) Sync_tree->Fill();


    vector<unsigned int> idx_vec;

    idx_vec.push_back(Muon_index_1);
    idx_vec.push_back(Muon_index_2);
    idx_vec.push_back(Muon_index_3);

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    TLorentzVector MuonOS = Ntp->Muon_P4(os_mu_idx);
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);


    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);



    if(id==1){
      if((TauRefitLV.M() > 1.65 && TauRefitLV.M() < 1.73 )  ||   (TauRefitLV.M() > 1.82 && TauRefitLV.M() < 1.9)){
	L1Triggers.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunB"))L1TriggersB.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunC"))L1TriggersC.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunD"))L1TriggersD.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunE"))L1TriggersE.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunF"))L1TriggersF.at(t).Fill(DoubleMuFired,TripleMuFired);
      }
    }
    if(id!=1)L1Triggers.at(t).Fill(DoubleMuFired,TripleMuFired);




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

    EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());
    //    TauMass.at(t).Fill(TauLV.M(),1);
    //    TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    

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

    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(signal_idx,0),
	  Ntp->Vertex_d0sig_reco(signal_idx,1),
	  Ntp->Vertex_d0sig_reco(signal_idx,2)});

    var_MaxD0Significance = MaxD0Significance;
    var_IsolationMinDist = Ntp->Isolation_MinDist(signal_idx);



    var_tauMass=TauRefitLV.M();
    TauMass_all.at(t).Fill(TauRefitLV.M(),1);


  }
}


void  SynchronizeSignal::Finish(){

  file= new TFile("Sync_signal_tree_UF.root","recreate");
  Sync_tree->SetDirectory(file);

  file->Write();
  file->Close();

  Selection::Finish();
}





