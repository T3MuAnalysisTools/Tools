#include "BParkingSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TRegexp.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace std;

BParkingSelector::BParkingSelector(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

BParkingSelector::~BParkingSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  BParkingSelector::Configure(){


  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("bdt",&bdt);
  T3MMiniTree->Branch("category",&category);
  T3MMiniTree->Branch("rapidity",&rapidity);
  T3MMiniTree->Branch("LumiScale",&LumiScale);


  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  readerA_3glb = new TMVA::Reader( "!Color:!Silent" );
  readerA_3glb->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerA_3glb->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerA_3glb->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerA_3glb->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerA_3glb->AddVariable( "var_pmin", &var_pmin);
  readerA_3glb->AddVariable( "var_max_cLP", &var_max_cLP);
  readerA_3glb->AddVariable( "var_max_tKink", &var_max_tKink);
  readerA_3glb->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerA_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerA_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
  readerA_3glb->AddSpectator("var_tauMass",&var_tauMass);
  readerA_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerA_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
  readerA_3glb->BookMVA( "BDT", basedir+"TMVAClassification_ThreeGlbMu_Veto_2016vars_A/weights/TMVAClassification_ThreeGlbMu_Veto_2016vars_A_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles
  
  readerB_3glb = new TMVA::Reader( "!Color:!Silent" );
  readerB_3glb->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerB_3glb->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerB_3glb->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerB_3glb->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerB_3glb->AddVariable( "var_pmin", &var_pmin);
  readerB_3glb->AddVariable( "var_max_cLP", &var_max_cLP);
  readerB_3glb->AddVariable( "var_max_tKink", &var_max_tKink);
  readerB_3glb->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerB_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerB_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
  readerB_3glb->AddSpectator("var_tauMass",&var_tauMass);
  readerB_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerB_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);

  readerB_3glb->BookMVA( "BDT", basedir+"TMVAClassification_ThreeGlbMu_Veto_2016vars_B/weights/TMVAClassification_ThreeGlbMu_Veto_2016vars_B_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles


  readerC_3glb = new TMVA::Reader( "!Color:!Silent" );
  readerC_3glb->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerC_3glb->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerC_3glb->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerC_3glb->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerC_3glb->AddVariable( "var_pmin", &var_pmin);
  readerC_3glb->AddVariable( "var_max_cLP", &var_max_cLP);
  readerC_3glb->AddVariable( "var_max_tKink", &var_max_tKink);
  readerC_3glb->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerC_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerC_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
  readerC_3glb->AddSpectator("var_tauMass",&var_tauMass);
  readerC_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerC_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);




  readerC_3glb->BookMVA( "BDT", basedir+"TMVAClassification_ThreeGlbMu_Veto_2016vars_C/weights/TMVAClassification_ThreeGlbMu_Veto_2016vars_C_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

  readerA_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
  readerA_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerA_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerA_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerA_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerA_2glbTrk->AddVariable( "var_pmin", &var_pmin);
  readerA_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
  readerA_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
  readerA_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerA_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerA_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
  readerA_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
  readerA_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerA_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);



  readerA_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_A/weights/TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_A_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

  readerB_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
  readerB_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerB_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerB_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerB_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerB_2glbTrk->AddVariable( "var_pmin", &var_pmin);
  readerB_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
  readerB_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
  readerB_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerB_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerB_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
  readerB_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
  readerB_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerB_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);





  readerB_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_B/weights/TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_B_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles


  readerC_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
  readerC_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
  readerC_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
  readerC_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
  readerC_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
  readerC_2glbTrk->AddVariable( "var_pmin", &var_pmin);
  readerC_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
  readerC_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
  readerC_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
  //readerC_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
  readerC_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
  readerC_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
  readerC_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
  readerC_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);

  readerC_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_C/weights/TMVAClassification_TwoGlbMuOneTrk_Veto_2016vars_C_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

   // dataloader->AddVariable("var_vertexKFChi2","Variable vertexKFChi2","units", 'F' );
   // dataloader->AddVariable("var_svpvTauAngle","Variable svpvTauAngle","units", 'F' );
   // dataloader->AddVariable("var_flightLenSig","Variable flightLenSig","units", 'F' );
   // dataloader->AddVariable("var_sumMuTrkKinkChi2","Variable sumMuTrkKinkChi2","units", 'F' );
   // dataloader->AddVariable("var_segCompMuMin","Variable segCompMuMin","units", 'F' );
   // dataloader->AddVariable("var_MinMIPLikelihood","Variable MinMIPLikelihood","units", 'F' );
   // dataloader->AddSpectator("var_tauMass","Variable tauMass","units", 'F' );
   // dataloader->AddVariable("var_maxdca","Variable maxdca","units", 'F' );
   // dataloader->AddVariable("var_RelPt_Mu1Tau","Variable RelPt_Mu1Tau","units", 'F' );
   // dataloader->AddVariable("var_MinD0Significance","Variable MinD0Significance","units", 'F' );
   // dataloader->AddVariable("var_IsolationMinDist","Variable IsolationMinDist","units", 'F' );




  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    //    if(i==VertChi2)           cut.at(VertChi2)=20.0;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch1)       cut.at(TriggerMatch1)=0.03;
    if(i==TriggerMatch2)       cut.at(TriggerMatch2)=0.03;
    if(i==TriggerMatch3)       cut.at(TriggerMatch3)=0.03;
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
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    /*  else if(i==VertChi2){
      title.at(i)="Triple Vertex Chi Squared $<$ 15";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Chi squared of triple vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VertChi2_",htitle,50,0,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VertChi2_",htitle,50,0,25,hlabel,"Events"));
      }*/





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
      title.at(i)="$p_{T}(\\mu_{3}) >$ 2 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Muon3 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
    }
    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="gl,gl,gl/tr";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PhiVeto){
      title.at(i)="$\\phi$ mass veto";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Phi mass Veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,60,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto){
      title.at(i)="$\\omega$ mass veto";
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
      title.at(i)="$\\tau$ mass (sideband in data)";
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

  TriggersDecision = HConfig.GetTH2D(Name+"_TriggersDecision","TriggersDecision",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu trigger","SingleMu");

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

  EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");
// category-wise Tau mass plots

  TauMass_allVsBDT_3glb_A=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_3glb_A","3#mu mass vs BDTa (3 glb)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDT_3glb_B=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_3glb_B","3#mu mass vs BDTb (3 glb)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDT_3glb_C=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_3glb_C","3#mu mass vs BDTc (3 glb)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_3glb_A1 =HConfig.GetTH1D(Name+"_TauMass_3glb_A1","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_A1 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_A1","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_3glb_B1 =HConfig.GetTH1D(Name+"_TauMass_3glb_B1","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_B1 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_B1","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_3glb_C1 =HConfig.GetTH1D(Name+"_TauMass_3glb_C1","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_C1 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_C1","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_3glb_A2 =HConfig.GetTH1D(Name+"_TauMass_3glb_A2","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_A2 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_A2","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_3glb_B2 =HConfig.GetTH1D(Name+"_TauMass_3glb_B2","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_B2 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_B2","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_3glb_C2 =HConfig.GetTH1D(Name+"_TauMass_3glb_C2","#tau lepton mass (3 glb)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_3glb_C2 =HConfig.GetTH1D(Name+"_TauMassRefit_3glb_C2","Refit #tau lepton mass (3 glb)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");

  TauMass_allVsBDT_2glbTrk_A=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_2glbTrk_A","3#mu mass vs BDTa (2 glb + trk)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDT_2glbTrk_B=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_2glbTrk_B","3#mu mass vs BDTb (2 glb + trk)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDT_2glbTrk_C=HConfig.GetTH2D(Name+"_TauMass_allVsBDT_2glbTrk_C","3#mu mass vs BDTc (2 glb + trk)",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_2glbTrk_A1 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_A1","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_A1 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_A1","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_2glbTrk_B1 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_B1","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_B1 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_B1","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_2glbTrk_C1 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_C1","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_C1 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_C1","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_2glbTrk_A2 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_A2","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_A2 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_A2","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_2glbTrk_B2 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_B2","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_B2 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_B2","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");
  TauMass_2glbTrk_C2 =HConfig.GetTH1D(Name+"_TauMass_2glbTrk_C2","#tau lepton mass (2 glb + trk)",50,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefit_2glbTrk_C2 =HConfig.GetTH1D(Name+"_TauMassRefit_2glbTrk_C2","Refit #tau lepton mass (2 glb + trk)",50,1.5,2.1,"KF refit  M_{#tau} , GeV","Events");

// End of category-wise Tau mass plots
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

  BDTOutput_3glb_A = HConfig.GetTH1D(Name+"_BDTOutput_3glb_A","BDTOutput_3glb_A",50,-0.4,0.4,"BDT Output","Events");
  BDTOutput_3glb_B = HConfig.GetTH1D(Name+"_BDTOutput_3glb_B","BDTOutpu_3glb_B",50,-0.4,0.4,"BDT Output","Events");
  BDTOutput_3glb_C = HConfig.GetTH1D(Name+"_BDTOutput_3glb_C","BDTOutpu_3glb_C",50,-0.4,0.4,"BDT Output","Events");

  BDTOutput_2glbTrk_A = HConfig.GetTH1D(Name+"_BDTOutputA_2glbTrk_","BDTOutput_2glbTrk_A",50,-0.4,0.4,"BDT Output","Events");
  BDTOutput_2glbTrk_B = HConfig.GetTH1D(Name+"_BDTOutputB_2glbTrk_","BDTOutput_2glbTrk_B",50,-0.4,0.4,"BDT Output","Events");
  BDTOutput_2glbTrk_C = HConfig.GetTH1D(Name+"_BDTOutputC_2glbTrk_","BDTOutput_2glbTrk_C",50,-0.4,0.4,"BDT Output","Events");

  L1Triggers=HConfig.GetTH2D(Name+"_L1Triggers","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  L1Triggers2017B=HConfig.GetTH2D(Name+"_L1Triggers2017B","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2017C=HConfig.GetTH2D(Name+"_L1Triggers2017C","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2017D=HConfig.GetTH2D(Name+"_L1Triggers2017D","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2017E=HConfig.GetTH2D(Name+"_L1Triggers2017E","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2017F=HConfig.GetTH2D(Name+"_L1Triggers2017F","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2018A=HConfig.GetTH2D(Name+"_L1Triggers2018A","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2018B=HConfig.GetTH2D(Name+"_L1Triggers2018B","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2018C=HConfig.GetTH2D(Name+"_L1Triggers2018C","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");
  L1Triggers2018D=HConfig.GetTH2D(Name+"_L1Triggers2018D","DoubleMuVsTripleMu",2,-0.5,1.5,2,-0.5,1.5,"DoubleMu","TripleMu");

  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  BParkingSelector::Store_ExtraDist(){ 


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


// category-wise Tau mass plots
// 3 global muons
  Extradist1d.push_back(&TauMass_3glb_A1);
  Extradist1d.push_back(&TauMassRefit_3glb_A1);

  Extradist1d.push_back(&TauMass_3glb_B1);
  Extradist1d.push_back(&TauMassRefit_3glb_B1);

  Extradist1d.push_back(&TauMass_3glb_C1);
  Extradist1d.push_back(&TauMassRefit_3glb_C1);


  Extradist1d.push_back(&TauMass_3glb_A2);
  Extradist1d.push_back(&TauMassRefit_3glb_A2);

  Extradist1d.push_back(&TauMass_3glb_B2);
  Extradist1d.push_back(&TauMassRefit_3glb_B2);

  Extradist1d.push_back(&TauMass_3glb_C2);
  Extradist1d.push_back(&TauMassRefit_3glb_C2);
  Extradist2d.push_back(&TauMass_allVsBDT_3glb_A);
  Extradist2d.push_back(&TauMass_allVsBDT_3glb_B);
  Extradist2d.push_back(&TauMass_allVsBDT_3glb_C);
  
// 2 global muons + one tracker
  Extradist1d.push_back(&TauMass_2glbTrk_A1);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_A1);

  Extradist1d.push_back(&TauMass_2glbTrk_B1);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_B1);

  Extradist1d.push_back(&TauMass_2glbTrk_C1);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_C1);


  Extradist1d.push_back(&TauMass_2glbTrk_A2);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_A2);

  Extradist1d.push_back(&TauMass_2glbTrk_B2);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_B2);

  Extradist1d.push_back(&TauMass_2glbTrk_C2);
  Extradist1d.push_back(&TauMassRefit_2glbTrk_C2);

  Extradist2d.push_back(&TauMass_allVsBDT_2glbTrk_A);
  Extradist2d.push_back(&TauMass_allVsBDT_2glbTrk_B);
  Extradist2d.push_back(&TauMass_allVsBDT_2glbTrk_C);
// End of category-wise Tau mass plots

  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist2d.push_back(&TauMass_all_nophiVeto);
  Extradist1d.push_back(&TauMass_all);

  Extradist2d.push_back(&EMR_tau_eta);



  Extradist1d.push_back(&SVPVTauDirAngle);
  Extradist2d.push_back(&TriggersDecision);
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


  Extradist1d.push_back(&BDTOutput_2glbTrk_A);
  Extradist1d.push_back(&BDTOutput_2glbTrk_B);
  Extradist1d.push_back(&BDTOutput_2glbTrk_C);

  Extradist1d.push_back(&BDTOutput_3glb_A);
  Extradist1d.push_back(&BDTOutput_3glb_B);
  Extradist1d.push_back(&BDTOutput_3glb_C);
 

  Extradist2d.push_back(&L1Triggers);

  Extradist2d.push_back(&L1Triggers2017B);
  Extradist2d.push_back(&L1Triggers2017C);
  Extradist2d.push_back(&L1Triggers2017D);
  Extradist2d.push_back(&L1Triggers2017E);
  Extradist2d.push_back(&L1Triggers2017F);

  Extradist2d.push_back(&L1Triggers2018A);
  Extradist2d.push_back(&L1Triggers2018B);
  Extradist2d.push_back(&L1Triggers2018C);
  Extradist2d.push_back(&L1Triggers2018D);
}


void  BParkingSelector::doEvent(){ 

  
  unsigned int t;
  bool threeGlobal = false;

  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu_fired(false);
  bool Parked_fired(false);
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
	 TRegexp parkingHLT("HLT_Mu*_IP*_v");
       if(id==1){
	 if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ){ doubleMu_counter++; HLTOk=true; DoubleMu_fired = true;}
         if(HLT.Contains("HLT_Mu") && HLT.Contains("IP") && Ntp->HLTDecision(iTrigger)  ) {singleMu_counter++; HLTOk = true; Parked_fired = true;}
		 }

        if(id!=1){

	  if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) { doubleMu_counter++; HLTOk=true; DoubleMu_fired =true;}
	  if(HLT.Contains("HLT_Mu") && HLT.Contains("IP") && Ntp->HLTDecision(iTrigger)  ) { singleMu_counter++; HLTOk = true; Parked_fired = true;}
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



  //if(DoubleMuFired  or TripleMuFired) L1Ok = true;
  L1Ok = true;
  value.at(TriggerOk)=(HLTOk && L1Ok);
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
    //    value.at(VertChi2) = Ntp->Vertex_Signal_KF_Chi2(signal_idx);
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
    			(Ntp->Muon_isGlobalMuon(mu3_pt_idx) || Ntp->Muon_isTrackerMuon(mu3_pt_idx)));

	 if (Ntp->Muon_isGlobalMuon(mu3_pt_idx)) threeGlobal = true;
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
    //    value.at(TauMassCut) = TauLV.M();
    value.at(TauMassCut) = TauRefittedLV.M();
  
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  //  pass.at(VertChi2) = (value.at(VertChi2) <= cut.at(VertChi2));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch1) = true;//(value.at(TriggerMatch1)  <  cut.at(TriggerMatch1));
  pass.at(TriggerMatch2) = true;//(value.at(TriggerMatch2)  <  cut.at(TriggerMatch2));
  pass.at(TriggerMatch3) = true;//(value.at(TriggerMatch3)  <  cut.at(TriggerMatch3));
  pass.at(PhiVeto)      = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 8*PDG_Var::Phi_width());
  pass.at(OmegaVeto)    = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 3*PDG_Var::Omega_width());

  if(id!=1) pass.at(TauMassCut) = true;
  else  pass.at(TauMassCut) =( (value.at(TauMassCut) < tauMinMass_)  ||   (value.at(TauMassCut)> tauMaxMass_ ));

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(PhiVeto);
  exclude_cuts.push_back(OmegaVeto);

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
    
    TriggersDecision.at(t).Fill(DoubleMu_fired, Parked_fired);

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
	if(Ntp->WhichEra(2017).Contains("RunB"))L1Triggers2017B.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunC"))L1Triggers2017C.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunD"))L1Triggers2017D.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunE"))L1Triggers2017E.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2017).Contains("RunF"))L1Triggers2017F.at(t).Fill(DoubleMuFired,TripleMuFired);
   if(Ntp->WhichEra(2018).Contains("RunA"))L1Triggers2018A.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2018).Contains("RunB"))L1Triggers2018B.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2018).Contains("RunC"))L1Triggers2018C.at(t).Fill(DoubleMuFired,TripleMuFired);
	if(Ntp->WhichEra(2018).Contains("RunD"))L1Triggers2018D.at(t).Fill(DoubleMuFired,TripleMuFired);
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
  if (threeGlobal){ 
     var_max_tKink = std::max(Ntp->Muon_combinedQuality_trkKink(mu1_pt_idx), std::max(Ntp->Muon_combinedQuality_trkKink(mu2_pt_idx), Ntp->Muon_combinedQuality_trkKink(mu3_pt_idx)));
     }
  else {
    var_max_tKink = std::max(Ntp->Muon_combinedQuality_trkKink(mu1_pt_idx), Ntp->Muon_combinedQuality_trkKink(mu2_pt_idx));
    }

    var_max_cLP = std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1), std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2), Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)));
    var_pmin = std::min(Ntp->Muon_P4(Muon_index_1).P(), std::min(Ntp->Muon_P4(Muon_index_2).P(), Ntp->Muon_P4(Muon_index_3).P()));
//    var_mindca_iso = Ntp->Isolation_MinDist(signal_idx);


    float MinD0Significance = std::max({Ntp->Vertex_d0sig_reco(signal_idx,0),
	  Ntp->Vertex_d0sig_reco(signal_idx,1),
	  Ntp->Vertex_d0sig_reco(signal_idx,2)});

    var_MinD0Significance = MinD0Significance;


    var_tauMass=TauRefitLV.M();
    TauMass_all.at(t).Fill(TauRefitLV.M(),1);

    m3m = TauRefitLV.M();
    dataMCtype = id;
    event_weight = w;
    rapidity = TauLV.Eta();
    LumiScale=1;

    //---------------- define per event resolution categroies 
    if (threeGlobal){
      //Category A1(3 global)
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007 ){
        TauMass_allVsBDT_3glb_A.at(t).Fill(TauRefitLV.M(),readerA_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_A.at(t).Fill(    readerA_3glb->EvaluateMVA("BDT") );
        if(readerA_3glb->EvaluateMVA("BDT") > 0.0851216952639){
  		    TauMassRefit_3glb_A1.at(t).Fill(TauRefitLV.M(),1);    
	       TauMass_3glb_A1.at(t).Fill(TauLV.M(),1);
	       category=1;
	       bdt = readerA_3glb->EvaluateMVA("BDT");
        }
      }
  
      //Category B1(3 global)
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.01 && threeGlobal){
        TauMass_allVsBDT_3glb_B.at(t).Fill(TauRefitLV.M(),readerB_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_B.at(t).Fill(    readerB_3glb->EvaluateMVA("BDT") );
        if(readerB_3glb->EvaluateMVA("BDT") > 0.0127777643976){
 	     	 TauMassRefit_3glb_B1.at(t).Fill(TauRefitLV.M(),1);    
 	     	 TauMass_3glb_B1.at(t).Fill(TauLV.M(),1);
  		    category =2 ;
			 bdt = readerB_3glb->EvaluateMVA("BDT");
        }
      }
  
      //Category C1(3 global)
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01){
        TauMass_allVsBDT_3glb_C.at(t).Fill(TauRefitLV.M(),readerC_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_C.at(t).Fill(    readerC_3glb->EvaluateMVA("BDT") );
        if(readerC_3glb->EvaluateMVA("BDT") > 0.0996578295929){
 	     	 TauMassRefit_3glb_C1.at(t).Fill(TauRefitLV.M(),1);    
  	    	 TauMass_3glb_C1.at(t).Fill(TauLV.M(),1);
  	    	 category = 3;
  		    bdt = readerC_3glb->EvaluateMVA("BDT");
        }
      }
  
      //Category A2(3 global)  
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007){
        TauMass_allVsBDT_3glb_A.at(t).Fill(TauRefitLV.M(),readerA_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_A.at(t).Fill(    readerA_3glb->EvaluateMVA("BDT") );
        if(readerA_3glb->EvaluateMVA("BDT") > 0.0245207610239 && readerA_3glb->EvaluateMVA("BDT") < 0.0851216952639){
  		    TauMassRefit_3glb_A2.at(t).Fill(TauRefitLV.M(),1);    
  		    TauMass_3glb_A2.at(t).Fill(TauLV.M(),1);
  		    category=4;
  		    bdt = readerA_3glb->EvaluateMVA("BDT");
        }
      }
  
   
     //Category B2(3 global)
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.01){
        TauMass_allVsBDT_3glb_B.at(t).Fill(TauRefitLV.M(),readerB_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_B.at(t).Fill(    readerB_3glb->EvaluateMVA("BDT") );
        if(readerB_3glb->EvaluateMVA("BDT") > 0.0668544925093 && readerB_3glb->EvaluateMVA("BDT") < 0.0127777643976){
  		    TauMassRefit_3glb_B2.at(t).Fill(TauRefitLV.M(),1);    
  		    TauMass_3glb_B2.at(t).Fill(TauLV.M(),1);
  		    category =5 ;
  		    bdt = readerB_3glb->EvaluateMVA("BDT");
        }
      }
  
      //Category C2(3 global)
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01){
        TauMass_allVsBDT_3glb_C.at(t).Fill(TauRefitLV.M(),readerC_3glb->EvaluateMVA("BDT"));
        BDTOutput_3glb_C.at(t).Fill(    readerC_3glb->EvaluateMVA("BDT") );
        if(readerC_3glb->EvaluateMVA("BDT") >  0.0199134813988 && readerC_3glb->EvaluateMVA("BDT")< 0.0996578295929){
  		    TauMassRefit_3glb_C2.at(t).Fill(TauRefitLV.M(),1);    
  		    TauMass_3glb_C2.at(t).Fill(TauLV.M(),1);
  		    category = 6;
  		    bdt = readerC_3glb->EvaluateMVA("BDT");
        }
      }
    }

   else{
     //Category A1(2 global)
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007){
      TauMass_allVsBDT_2glbTrk_A.at(t).Fill(TauRefitLV.M(),readerA_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_A.at(t).Fill(    readerA_2glbTrk->EvaluateMVA("BDT") );
      if(readerA_2glbTrk->EvaluateMVA("BDT") > 0.147629725692){
		  TauMassRefit_2glbTrk_A1.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_A1.at(t).Fill(TauLV.M(),1);
		  category=7;
		  bdt = readerA_2glbTrk->EvaluateMVA("BDT");
      }
    }

    //Category B1(2 global)
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.01){
      TauMass_allVsBDT_2glbTrk_B.at(t).Fill(TauRefitLV.M(),readerB_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_B.at(t).Fill(    readerB_2glbTrk->EvaluateMVA("BDT") );
      if(readerB_2glbTrk->EvaluateMVA("BDT") > 0.0745167711418){
		  TauMassRefit_2glbTrk_B1.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_B1.at(t).Fill(TauLV.M(),1);
		  category =8 ;
		  bdt = readerB_2glbTrk->EvaluateMVA("BDT");
      }
    }

    //Category C1(2 global)
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01){
      TauMass_allVsBDT_2glbTrk_C.at(t).Fill(TauRefitLV.M(),readerC_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_C.at(t).Fill(    readerC_2glbTrk->EvaluateMVA("BDT") );
      if(readerC_2glbTrk->EvaluateMVA("BDT") > 0.121108058178){
		  TauMassRefit_2glbTrk_C1.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_C1.at(t).Fill(TauLV.M(),1);
		  category = 9;
		  bdt = readerC_2glbTrk->EvaluateMVA("BDT");
      }
    }

    //Category A2(2 global)  
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007){
      TauMass_allVsBDT_2glbTrk_A.at(t).Fill(TauRefitLV.M(),readerA_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_A.at(t).Fill(    readerA_2glbTrk->EvaluateMVA("BDT") );
      if(readerA_2glbTrk->EvaluateMVA("BDT") > -0.00190483084183 && readerA_2glbTrk->EvaluateMVA("BDT") < 0.147629725692){
		  TauMassRefit_2glbTrk_A2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_A2.at(t).Fill(TauLV.M(),1);
		  category=10;
		  bdt = readerA_2glbTrk->EvaluateMVA("BDT");
      }
    }

 
   //Category B2(2 global)
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.01){
      TauMass_allVsBDT_2glbTrk_B.at(t).Fill(TauRefitLV.M(),readerB_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_B.at(t).Fill(    readerB_2glbTrk->EvaluateMVA("BDT") );
      if(readerB_2glbTrk->EvaluateMVA("BDT") >  -0.040602412133 && readerB_2glbTrk->EvaluateMVA("BDT") < 0.0745167711418){
		  TauMassRefit_2glbTrk_B2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_B2.at(t).Fill(TauLV.M(),1);
		  category =11 ;
		  bdt = readerB_2glbTrk->EvaluateMVA("BDT");
      }
    }

    //Category C2(2 global)
    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01){
      TauMass_allVsBDT_2glbTrk_C.at(t).Fill(TauRefitLV.M(),readerC_2glbTrk->EvaluateMVA("BDT"));
      BDTOutput_2glbTrk_C.at(t).Fill(    readerC_2glbTrk->EvaluateMVA("BDT") );
      if(readerC_2glbTrk->EvaluateMVA("BDT") > -0.0241806939129 && readerC_2glbTrk->EvaluateMVA("BDT")< 0.121108058178){
		  TauMassRefit_2glbTrk_C2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMass_2glbTrk_C2.at(t).Fill(TauLV.M(),1);
		  category = 12;
		  bdt = readerC_2glbTrk->EvaluateMVA("BDT");
      }
    }
    }
	 T3MMiniTree->Fill();

    //---------------  Fill MC plots 
	 /*    if(id==40 || id == 60 || id ==90){
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
      }*/


  }
 }
}


void  BParkingSelector::Finish(){
  cout<<"Number of singleMu HLTs fired: "<<singleMu_counter<<endl;
  cout<<"Number of doubleMu HLTs fired: "<<doubleMu_counter<<endl;

  T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();

  if(mode == RECONSTRUCT){
    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    //    int id(Ntp->GetMCID());
    //    double scale(1.);
    //    double scaleDsTau(0.637);
    //    double scaleBpTau(0.262);
    //    double scaleB0Tau(0.099);

    //total xsection of producing taus is 12.848 ub 


    //    if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    // ScaleAllHistOfType(2,scale*scaleDsTau);
    
    //if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    // ScaleAllHistOfType(3,scale*scaleB0Tau);

    //if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    //ScaleAllHistOfType(4,scale*scaleBpTau);

    //    }
  }
  Selection::Finish();
}





