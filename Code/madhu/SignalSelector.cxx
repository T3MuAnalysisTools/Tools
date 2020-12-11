#include "SignalSelector.h"
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

SignalSelector::SignalSelector(TString Name_, TString id_):
   Selection(Name_,id_),
   tauMinMass_(1.73),
   tauMaxMass_(1.82),
   tauMinSideBand_(1.65),
   tauMaxSideBand_(2.02),
   tauMassResCutLow(0.007),
   tauMassResCutHigh(0.01),
   phiVetoSigma(0.03),
   omegaVetoSigma(0.03),
   bdt_cutA2_3glb(0.0178036087),
   bdt_cutA1_3glb(0.13855843498),
   bdt_cutA2_2glbTrk(-0.00931573728943),
   bdt_cutA1_2glbTrk(0.0900170166134),
   bdt_cutB2_3glb(0.0364416637001),
   bdt_cutB1_3glb(0.1669982497),
   bdt_cutB2_2glbTrk(0.00258006234741),
   bdt_cutB1_2glbTrk(0.117157105615),
   bdt_cutC2_3glb(0.032225350338),
   bdt_cutC1_3glb(0.121164828298),
   bdt_cutC2_2glbTrk(-0.0483397164297),
   bdt_cutC1_2glbTrk(0.0750508898802)
{
   // This is a class constructor;
   T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
}

SignalSelector::~SignalSelector(){
   for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
         << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
         << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
   }
   Logger(Logger::Info) << "complete." << std::endl;
}

void  SignalSelector::Configure(){

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
   readerA_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerA_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
   readerA_3glb->AddSpectator("var_tauMass",&var_tauMass);
   readerA_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerA_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerA_3glb->AddSpectator("threeGlobal",&threeGlobal);



   readerA_3glb->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_A_threeGlobal/weights/TMVAClassification_2016vars_A_threeGlobal_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

   readerB_3glb = new TMVA::Reader( "!Color:!Silent" );
   readerB_3glb->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
   readerB_3glb->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
   readerB_3glb->AddVariable( "var_flightLenSig", &var_flightLenSig );
   readerB_3glb->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
   readerB_3glb->AddVariable( "var_pmin", &var_pmin);
   readerB_3glb->AddVariable( "var_max_cLP", &var_max_cLP);
   readerB_3glb->AddVariable( "var_max_tKink", &var_max_tKink);
   readerB_3glb->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
   readerB_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerB_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
   readerB_3glb->AddSpectator("var_tauMass",&var_tauMass);
   readerB_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerB_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerB_3glb->AddSpectator("threeGlobal",&threeGlobal);




   readerB_3glb->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_B_threeGlobal/weights/TMVAClassification_2016vars_B_threeGlobal_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles


   readerC_3glb = new TMVA::Reader( "!Color:!Silent" );
   readerC_3glb->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
   readerC_3glb->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
   readerC_3glb->AddVariable( "var_flightLenSig", &var_flightLenSig );
   readerC_3glb->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
   readerC_3glb->AddVariable( "var_pmin", &var_pmin);
   readerC_3glb->AddVariable( "var_max_cLP", &var_max_cLP);
   readerC_3glb->AddVariable( "var_max_tKink", &var_max_tKink);
   readerC_3glb->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
   readerC_3glb->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerC_3glb->AddVariable( "var_trk_relPt", &var_trk_relPt );
   readerC_3glb->AddSpectator("var_tauMass",&var_tauMass);
   readerC_3glb->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerC_3glb->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerC_3glb->AddSpectator("threeGlobal",&threeGlobal);




   readerC_3glb->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_C_threeGlobal/weights/TMVAClassification_2016vars_C_threeGlobal_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

   readerA_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
   readerA_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
   readerA_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
   readerA_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
   readerA_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
   readerA_2glbTrk->AddVariable( "var_pmin", &var_pmin);
   readerA_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
   readerA_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
   readerA_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
   readerA_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerA_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
   readerA_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
   readerA_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerA_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerA_2glbTrk->AddSpectator("threeGlobal",&threeGlobal);



   readerA_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_A_twoGlobaltrk/weights/TMVAClassification_2016vars_A_twoGlobaltrk_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

   readerB_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
   readerB_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
   readerB_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
   readerB_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
   readerB_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
   readerB_2glbTrk->AddVariable( "var_pmin", &var_pmin);
   readerB_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
   readerB_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
   readerB_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
   readerB_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerB_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
   readerB_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
   readerB_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerB_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerB_2glbTrk->AddSpectator("threeGlobal",&threeGlobal);



   readerB_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_B_twoGlobaltrk/weights/TMVAClassification_2016vars_B_twoGlobaltrk_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles


   readerC_2glbTrk = new TMVA::Reader( "!Color:!Silent" );
   readerC_2glbTrk->AddVariable( "var_vertexKFChi2", &var_vertexKFChi2 );
   readerC_2glbTrk->AddVariable( "var_svpvTauAngle", &var_svpvTauAngle );
   readerC_2glbTrk->AddVariable( "var_flightLenSig", &var_flightLenSig );
   readerC_2glbTrk->AddVariable( "var_segCompMuMin", &var_segCompMuMin );
   readerC_2glbTrk->AddVariable( "var_pmin", &var_pmin);
   readerC_2glbTrk->AddVariable("var_max_cLP", &var_max_cLP);
   readerC_2glbTrk->AddVariable("var_max_tKink", &var_max_tKink);
   readerC_2glbTrk->AddVariable( "var_MinD0Significance", &var_MinD0Significance );
   readerC_2glbTrk->AddVariable( "var_mindca_iso", &var_mindca_iso );
   readerC_2glbTrk->AddVariable("var_trk_relPt", &var_trk_relPt );
   readerC_2glbTrk->AddSpectator("var_tauMass",&var_tauMass);
   readerC_2glbTrk->AddSpectator("var_tauMassRes",&var_tauMassRes);
   readerC_2glbTrk->AddSpectator("var_Eta_Tau",&var_Eta_Tau);
   readerC_2glbTrk->AddSpectator("threeGlobal",&threeGlobal);

   readerC_2glbTrk->BookMVA( "BDT", basedir+"TMVAClassification_2016vars_C_twoGlobaltrk/weights/TMVAClassification_2016vars_C_twoGlobaltrk_BDT.weights.xml" ); // weights weights.xml file after training, place it to CommonFiles

   for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==TriggerOk)          cut.at(TriggerOk)=1;
      if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
      if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
      if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
      if(i==MuonID)             cut.at(MuonID)=1;
      if(i==PVRefit) 			  cut.at(PVRefit)=1;
      if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
      if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
      if(i==TriggerMatchMu1)    cut.at(TriggerMatchMu1)=0.03;
      if(i==TriggerMatchMu2)    cut.at(TriggerMatchMu2)=0.03;
      if(i==TriggerMatchMu3)    cut.at(TriggerMatchMu3)=0.03;
      if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
      if(i==GenMatch)	  		  cut.at(GenMatch)=0.03;
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
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
      }
      else if(i==Mu1PtCut){
         title.at(i)="$p_{T}(\\mu_{1}) >$";
         title.at(i)+=cut.at(Mu1PtCut);
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="$p_{T}(\\mu_{2}) >$";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      }
      else if(i==Mu3PtCut){
         title.at(i)="$p_{T}(\\mu_{3}) >$";
         title.at(i)+=cut.at(Mu3PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon3 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="gl,gl,(gl/trk)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==PVRefit){
         title.at(i)="PV refit valid";
         hlabel="Primary Vertex refit valid";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
         hlabel="Omega mass veto, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu1){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{1}} <$ 0.03";
         hlabel="dR muon 1";
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");

         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu1_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu1_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu2){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{2}} <$ 0.03";
         hlabel="dR muon 2";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu2_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu2_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu3){
         title.at(i)="$\\Delta R(reco-trigger)_{\\mu_{3}} <$ 0.03";
         hlabel="dR muon 3";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu3_",htitle,40,0,0.05,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu3_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==TauMassCut){
         title.at(i)="$\\tau$ mass (sideband in data)";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");

         hlabel="three mu mass, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,80,1.4,2.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,80,1.4,2.2,hlabel,"Events"));
      }
      else if(i==GenMatch){
         title.at(i)="GEN matching";
         hlabel="GEN match";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
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

   Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"#mu_{1} p_{T}, GeV","Events");
   Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"#mu_{2} p_{T}, GeV","Events");
   Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"#mu_{3} p_{T}, GeV","Events");

   Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#eta(#mu_{1})","Events");
   Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#eta(#mu_{2})","Events");
   Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#eta(#mu_{3})","Events");

   TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
   TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
   TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");

   TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
   TauMassResolutionRefit_3glb=HConfig.GetTH1D(Name+"_TauMassResolutionRefit_3glb","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
   TauMassResolutionRefit_2glbtrk=HConfig.GetTH1D(Name+"_TauMassResolutionRefit_2glbtrk","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
   TauMassResolution_3glb=HConfig.GetTH1D(Name+"_TauMassResolution_3glb","TauMassResolution",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
   TauMassResolution_2glbtrk=HConfig.GetTH1D(Name+"_TauMassResolution_2glbtrk","TauMassResolution",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

   TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");
   TauMass_all =HConfig.GetTH1D(Name+"_TauMass_all","3#mu  mass",60,1.5,2.1,"  M_{#tau} , GeV","Events");

   TauMass_allVsBDTBarrel=HConfig.GetTH2D(Name+"_TauMass_allVsBDTBarrel","3#mu mass vs BDTBarre;",60,1.6,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
   TauMass_allVsBDTEndcap=HConfig.GetTH2D(Name+"_TauMass_allVsBDTEndcap","3#mu mass vs BDTEndcap",60,1.6,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");

   EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

   // Invariant mass (3 global category)

   TauMass_allVsBDTA_3glb=HConfig.GetTH2D(Name+"_TauMass_allVsBDTA_3glb","3#mu mass vs BDTa",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
   TauMass_allVsBDTB_3glb=HConfig.GetTH2D(Name+"_TauMass_allVsBDTB_3glb","3#mu mass vs BDTb",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
   TauMass_allVsBDTC_3glb=HConfig.GetTH2D(Name+"_TauMass_allVsBDTC_3glb","3#mu mass vs BDTc",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");

   TauMassA1_3glb =HConfig.GetTH1D(Name+"_TauMassA1_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitA1_3glb =HConfig.GetTH1D(Name+"_TauMassRefitA1_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassB1_3glb =HConfig.GetTH1D(Name+"_TauMassB1_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitB1_3glb =HConfig.GetTH1D(Name+"_TauMassRefitB1_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassC1_3glb =HConfig.GetTH1D(Name+"_TauMassC1_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitC1_3glb =HConfig.GetTH1D(Name+"_TauMassRefitC1_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassA2_3glb =HConfig.GetTH1D(Name+"_TauMassA2_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitA2_3glb =HConfig.GetTH1D(Name+"_TauMassRefitA2_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassB2_3glb =HConfig.GetTH1D(Name+"_TauMassB2_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitB2_3glb =HConfig.GetTH1D(Name+"_TauMassRefitB2_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassC2_3glb =HConfig.GetTH1D(Name+"_TauMassC2_3glb","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitC2_3glb =HConfig.GetTH1D(Name+"_TauMassRefitC2_3glb","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau}, GeV","Events");

   // Invariant mass (2 global + trk category)

   TauMass_allVsBDTA_2glbtrk=HConfig.GetTH2D(Name+"_TauMass_allVsBDTA_2glbtrk","3#mu mass vs BDTa",25,1.65,1.9,50,-0.3,0.3,"3#mu mass, GeV","BDT");
   TauMass_allVsBDTB_2glbtrk=HConfig.GetTH2D(Name+"_TauMass_allVsBDTB_2glbtrk","3#mu mass vs BDTb",25,1.65,1.9,50,-0.3,0.3,"3#mu mass, GeV","BDT");
   TauMass_allVsBDTC_2glbtrk=HConfig.GetTH2D(Name+"_TauMass_allVsBDTC_2glbtrk","3#mu mass vs BDTc",25,1.65,1.9,50,-0.3,0.3,"3#mu mass, GeV","BDT");

   TauMassA1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassA1_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau}, GeV","Events");
   TauMassRefitA1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitA1_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassB1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassB1_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitB1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitB1_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassC1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassC1_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitC1_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitC1_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassA2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassA2_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitA2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitA2_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassB2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassB2_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitB2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitB2_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassC2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassC2_2glbtrk","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitC2_2glbtrk =HConfig.GetTH1D(Name+"_TauMassRefitC2_2glbtrk","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassBarrel1 =HConfig.GetTH1D(Name+"_TauMassBarrel1","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitBarrel1 =HConfig.GetTH1D(Name+"_TauMassRefitBarrel1","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassEndcap1 =HConfig.GetTH1D(Name+"_TauMassEndcap1","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitEndcap1 =HConfig.GetTH1D(Name+"_TauMassRefitEndcap1","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassBarrel2 =HConfig.GetTH1D(Name+"_TauMassBarrel2","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitBarrel2 =HConfig.GetTH1D(Name+"_TauMassRefitBarrel2","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");

   TauMassEndcap2 =HConfig.GetTH1D(Name+"_TauMassEndcap2","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassRefitEndcap2 =HConfig.GetTH1D(Name+"_TauMassRefitEndcap2","Refit #tau lepton mass",25,1.65,1.9,"KF refit  M_{#tau} , GeV","Events");


   EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
   VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
   FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

   SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

   Muon1DRToTruth_3glb=HConfig.GetTH1D(Name+"_Muon1DRToTruth_3glb","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
   Muon2DRToTruth_3glb=HConfig.GetTH1D(Name+"_Muon2DRToTruth_3glb","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
   Muon3DRToTruth_3glb=HConfig.GetTH1D(Name+"_Muon3DRToTruth_3glb","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

   Muon1DRToTruth_2glbtrk=HConfig.GetTH1D(Name+"_Muon1DRToTruth_2glbtrk","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
   Muon2DRToTruth_2glbtrk=HConfig.GetTH1D(Name+"_Muon2DRToTruth_2glbtrk","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
   Muon3DRToTruth_2glbtrk=HConfig.GetTH1D(Name+"_Muon3DRToTruth_2glbtrk","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

   TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,0.1,"trigger match #Delta R 1","Events");
   TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,0.1,"trigger match #Delta R 2","Events");
   TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,0.1,"trigger match #Delta R 3","Events");

   BDTOutputA_3glb = HConfig.GetTH1D(Name+"_BDTOutputA_3glb","BDTOutputA",50,-0.4,0.4,"BDT Output","Events");
   BDTOutputB_3glb = HConfig.GetTH1D(Name+"_BDTOutputB_3glb","BDTOutputB",50,-0.4,0.4,"BDT Output","Events");
   BDTOutputC_3glb = HConfig.GetTH1D(Name+"_BDTOutputC_3glb","BDTOutputC",50,-0.4,0.4,"BDT Output","Events");

   BDTOutputA_2glbtrk = HConfig.GetTH1D(Name+"_BDTOutputA_2glbtrk","BDTOutputA",50,-0.4,0.4,"BDT Output","Events");
   BDTOutputB_2glbtrk = HConfig.GetTH1D(Name+"_BDTOutputB_2glbtrk","BDTOutputB",50,-0.4,0.4,"BDT Output","Events");
   BDTOutputC_2glbtrk = HConfig.GetTH1D(Name+"_BDTOutputC_2glbtrk","BDTOutputC",50,-0.4,0.4,"BDT Output","Events");

   BDTOutputBarrel = HConfig.GetTH1D(Name+"_BDTOutputBarrel","BDTOutputBarrel",50,-0.4,0.4,"BDT Output","Events");
   BDTOutputEndcap = HConfig.GetTH1D(Name+"_BDTOutputEndcap","BDTOutputENdcap",50,-0.4,0.4,"BDT Output","Events");

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


void  SignalSelector::Store_ExtraDist(){ 

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


   // 3 global category
   Extradist2d.push_back(&TauMass_allVsBDTA_3glb);
   Extradist2d.push_back(&TauMass_allVsBDTB_3glb);
   Extradist2d.push_back(&TauMass_allVsBDTC_3glb);

   Extradist1d.push_back(&TauMassA1_3glb);
   Extradist1d.push_back(&TauMassRefitA1_3glb);

   Extradist1d.push_back(&TauMassB1_3glb);
   Extradist1d.push_back(&TauMassRefitB1_3glb);

   Extradist1d.push_back(&TauMassC1_3glb);
   Extradist1d.push_back(&TauMassRefitC1_3glb);


   Extradist1d.push_back(&TauMassA2_3glb);
   Extradist1d.push_back(&TauMassRefitA2_3glb);

   Extradist1d.push_back(&TauMassB2_3glb);
   Extradist1d.push_back(&TauMassRefitB2_3glb);

   Extradist1d.push_back(&TauMassC2_3glb);
   Extradist1d.push_back(&TauMassRefitC2_3glb);

   Extradist1d.push_back(&BDTOutputA_3glb);
   Extradist1d.push_back(&BDTOutputB_3glb);
   Extradist1d.push_back(&BDTOutputC_3glb);

   // 2 global + trk category
   Extradist2d.push_back(&TauMass_allVsBDTA_2glbtrk);
   Extradist2d.push_back(&TauMass_allVsBDTB_2glbtrk);
   Extradist2d.push_back(&TauMass_allVsBDTC_2glbtrk);

   Extradist1d.push_back(&TauMassA1_2glbtrk);
   Extradist1d.push_back(&TauMassRefitA1_2glbtrk);

   Extradist1d.push_back(&TauMassB1_2glbtrk);
   Extradist1d.push_back(&TauMassRefitB1_2glbtrk);

   Extradist1d.push_back(&TauMassC1_2glbtrk);
   Extradist1d.push_back(&TauMassRefitC1_2glbtrk);


   Extradist1d.push_back(&TauMassA2_2glbtrk);
   Extradist1d.push_back(&TauMassRefitA2_2glbtrk);

   Extradist1d.push_back(&TauMassB2_2glbtrk);
   Extradist1d.push_back(&TauMassRefitB2_2glbtrk);

   Extradist1d.push_back(&TauMassC2_2glbtrk);
   Extradist1d.push_back(&TauMassRefitC2_2glbtrk);

   Extradist1d.push_back(&BDTOutputA_2glbtrk);
   Extradist1d.push_back(&BDTOutputB_2glbtrk);
   Extradist1d.push_back(&BDTOutputC_2glbtrk);

   Extradist1d.push_back(&TauMassBarrel1);
   Extradist1d.push_back(&TauMassRefitBarrel1);
   Extradist1d.push_back(&TauMassEndcap1);
   Extradist1d.push_back(&TauMassRefitEndcap1);

   Extradist1d.push_back(&TauMassBarrel2);
   Extradist1d.push_back(&TauMassRefitBarrel2);
   Extradist1d.push_back(&TauMassEndcap2);
   Extradist1d.push_back(&TauMassRefitEndcap2);


   Extradist1d.push_back(&TauMassResolution_3glb);
   Extradist1d.push_back(&TauMassResolutionRefit_3glb);

   Extradist1d.push_back(&TauMassResolution_2glbtrk);
   Extradist1d.push_back(&TauMassResolutionRefit_2glbtrk);

   Extradist2d.push_back(&TauMass_all_nophiVeto);
   Extradist1d.push_back(&TauMass_all);

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


   Extradist1d.push_back(&BDTOutputBarrel);
   Extradist1d.push_back(&BDTOutputEndcap);

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


void  SignalSelector::doEvent(){ 

   value.at(TriggerOk)=0;
   value.at(SignalCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(Mu3PtCut)=0;
   value.at(MuonID)=0;
   value.at(TriggerMatchMu1)=1.0;
   value.at(TriggerMatchMu2)=1.0;
   value.at(TriggerMatchMu3)=1.0;
   value.at(PVRefit)=0;
   value.at(PhiVeto)=99.0;
   value.at(OmegaVeto)=99.0;
   value.at(TauMassCut)=0;
   value.at(GenMatch)=1;

   bool threeGlobal = false;
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
      if( (id==1 && Ntp->WhichEra(2017).Contains("RunF")) || (id==1 && Ntp->WhichEra(2018).Contains("Run")) ){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }
      if(id!=1){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }

   }

   // ------------ HLT selection -----------------
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
      if (id==1 && Ntp->WhichEra(2018).Contains("Run")){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))							DoubleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))	TripleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2"))								DoubleMuFired = Ntp->L1Decision(il1);
      }
      if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
      }
   }



   if(DoubleMuFired  or TripleMuFired) L1Ok = true;
   value.at(TriggerOk)=(HLTOk && L1Ok);
   pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));


   // ------------ Define new variables -----------------
   unsigned int  final_idx=0;

   double min_chi2(100.);
   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_tau = 99.0;


   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_tau = 0.;
   double sumPtTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_tau = 0;
   // ---------------------------------------------------

   // ------------ Chose candidate with minimum chi square -----------------
   for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
      if(Ntp->Vertex_Signal_KF_Chi2(i_idx,false) < min_chi2){
         min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx,false);
         final_idx = i_idx;
      }
   }
   // ---------------------------------------------------


   NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);
   value.at(SignalCandidate) = Ntp->NThreeMuons();

   if(Ntp->NThreeMuons()>0){
      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
      //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
      //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
      //
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) && 
            Ntp->Muon_isGlobalMuon(Muon_index_2) &&
            ( Ntp->Muon_isGlobalMuon(Muon_index_3) || Ntp->Muon_isTrackerMuon(Muon_index_3)));	 

      if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal = true;
      //------------------------------------------------------------------------------------------------------

      value.at(PVRefit) = (Ntp->Vertex_RefitPVisValid(final_idx)==1);
      value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();

      vector<unsigned int> idx_vec;

      idx_vec.push_back(mu1_idx);
      idx_vec.push_back(mu2_idx);
      idx_vec.push_back(mu3_idx);

      unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
      unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
      unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

      TLorentzVector TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);
      TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);


      double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
      double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

      value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
      value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

      value.at(TriggerMatchMu1) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(0);
      value.at(TriggerMatchMu2) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(1);
      value.at(TriggerMatchMu3) = Ntp-> ThreeMuons_TriggerMatch_dR(final_idx).at(2);
      //value.at(TauMassCut) = TauRefittedLV.M();
      value.at(TauMassCut) = TauLV.M();
   }
   if (id!=1) value.at(GenMatch) = Ntp->TauGenMatch(final_idx);
   else value.at(GenMatch) = 0;

   pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
   pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
   pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
   pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
   pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
   pass.at(PVRefit) = (value.at(PVRefit) == cut.at(PVRefit));
   pass.at(TriggerMatchMu1) = (value.at(TriggerMatchMu1)  <  cut.at(TriggerMatchMu1));
   pass.at(TriggerMatchMu2) = (value.at(TriggerMatchMu2)  <  cut.at(TriggerMatchMu2));
   pass.at(TriggerMatchMu3) = (value.at(TriggerMatchMu3)  <  cut.at(TriggerMatchMu3));
   //pass.at(TriggerMatchMu1) = true;
   //pass.at(TriggerMatchMu2) = true;
   //pass.at(TriggerMatchMu3) = true;
   pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 8*PDG_Var::Phi_width());
   pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 3*PDG_Var::Omega_width()); 

   if(id!=1) pass.at(TauMassCut) = true;
   else  pass.at(TauMassCut) = (value.at(TauMassCut)>tauMinSideBand_  && value.at(TauMassCut) < tauMaxSideBand_);


   //pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
   pass.at(GenMatch) = true;
   std::vector<unsigned int> exclude_cuts;
   exclude_cuts.push_back(PhiVeto);
   exclude_cuts.push_back(OmegaVeto);

   if(passAllBut(exclude_cuts)){


      TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);

      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
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
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

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
      
      TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
      TauEta.at(t).Fill(TauLV.Eta(),1);
      TauPt.at(t).Fill(TauLV.Pt(),1);
      TauP.at(t).Fill(TauLV.P(),1);
      TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false);


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

      TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0),1);
      TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1),1);
      TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2),1);

      VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,false),w);
      FLSignificance.at(t).Fill(sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,false),
                  Ntp->Vertex_PrimaryVertex_Covariance(final_idx,false),
                  Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx,false))),w);
      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,false),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));
      SVPVTauDirAngle.at(t).Fill(SVPV.Angle(TauLV.Vect()),w);

      double tauMassRes = Ntp->TauMassResolution(EtaSortedIndices,1,false);
      float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,false),
            Ntp->Vertex_d0sig_reco(final_idx,1,false),
            Ntp->Vertex_d0sig_reco(final_idx,2,false)});
      float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0),
            Ntp->Vertex_d0sig_reco(final_idx,1,false),
            Ntp->Vertex_d0sig_reco(final_idx,2,false)});

      double maxMuondR = std::max({Muon1LV.DeltaR(TauLV), Muon2LV.DeltaR(TauLV), Muon3LV.DeltaR(TauLV)});
      double minMuonPt = std::min({Muon1LV.Pt(), Muon2LV.Pt(), Muon3LV.Pt()});

      // Isolation algorithm
      for (int it=0; it<Ntp->NIsolationTrack(final_idx, false); it++){
         double dxy_track = Ntp->IsolationTrack_dxySV(final_idx, it, false);
         double dz_track = Ntp->IsolationTrack_dzSV(final_idx, it, false);
         TLorentzVector TrackLV = Ntp->IsolationTrack_p4(final_idx, it, false);
         double dca_fv = TMath::Sqrt(pow(dxy_track, 2)+ pow(dz_track, 2));

         double dr_tau = TauLV.DeltaR(TrackLV); 
         double dr_mu1 = Muon1LV.DeltaR(TrackLV);
         double dr_mu2 = Muon2LV.DeltaR(TrackLV);
         double dr_mu3 = Muon3LV.DeltaR(TrackLV);

         // Isolation 1
         if ( dca_fv<0.5 && TrackLV.Pt()<0.33*minMuonPt && dr_tau<3*maxMuondR ){
            sumPtTracks_tau += TrackLV.Pt();
            nTracks_tau++;
            if (dca_fv < mindca_tau) mindca_tau = dca_fv;
         }

         // Isolation 2
         if (TrackLV.Pt()<1.0) {
            continue;
         }

         if (dca_fv < mindca_iso) mindca_iso = dca_fv;

         // Isolation 3 (within dR = 0.5 of tau)
         if (dr_tau<0.5 && dca_fv<0.5){
            sumPtTracks_iso05 += TrackLV.Pt();
            nTracks_iso05++;
            if(dca_fv<mindca_iso05) mindca_iso05 = dca_fv;
         }

         // Isolation 4 (Muon isolation)
         if (dr_mu1 < 0.3 && Ntp->IsolationTrack_DocaMu1(final_idx, it, false) < 0.1 ) sumPtTracks_mu1 += TrackLV.Pt();
         if (dr_mu2 < 0.3 && Ntp->IsolationTrack_DocaMu2(final_idx, it, false) < 0.1 ) sumPtTracks_mu2 += TrackLV.Pt();
         if (dr_mu3 < 0.3 && Ntp->IsolationTrack_DocaMu3(final_idx, it, false) < 0.1 ) sumPtTracks_mu3 += TrackLV.Pt();
      }
      // Relative Pt calculation
      double mu1_relPt = sumPtTracks_mu1/Muon1LV.Pt();
      double mu2_relPt = sumPtTracks_mu2/Muon2LV.Pt();
      double mu3_relPt = sumPtTracks_mu3/Muon3LV.Pt();
      double relPt_iso05 = sumPtTracks_tau/TauLV.Pt();


      // Categorization variables
      var_Eta_Tau = TauLV.Eta();
      var_tauMassRes = tauMassRes;

      // Muon/Tau kinematic variables
      var_pmin = std::min(Ntp->Muon_P4(Muon_index_1).P(), std::min(Ntp->Muon_P4(Muon_index_2).P(), Ntp->Muon_P4(Muon_index_3).P()));
      var_RelPt_Mu1Tau = Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt();
      var_MuMu_mindR = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
      var_MuTau_maxdR = std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});

      // Muon ID variables
      var_max_cLP = std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1), std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2), Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)));
      var_max_tKink = std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_1), std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_2), Ntp->Muon_combinedQuality_trkKink(Muon_index_3)));
      var_segCompMuMin  = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      var_MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      var_sumMuTrkKinkChi2= (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));

      // vertex variables
      var_maxdca = std::max({Ntp->Vertex_DCA12(final_idx, false),Ntp->Vertex_DCA23(final_idx, false),Ntp->Vertex_DCA31(final_idx, false)});
      var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(final_idx, 0, false), Ntp->Vertex_pair_quality(final_idx, 1, false), Ntp->Vertex_pair_quality(final_idx, 2, false)});
      var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
      var_flightLenSig = sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx, false),Ntp->Vertex_PrimaryVertex_Covariance(final_idx, false),
               Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx, false)));
      var_vertexKFChi2 =Ntp->Vertex_signal_KF_Chi2(final_idx, false);
      // Isolation variables
      var_MinD0Significance = MinD0Significance;
      var_MaxD0Significance = MaxD0Significance;
      var_mindca_iso = mindca_iso;
      var_trk_relPt = std::max({mu1_relPt, mu2_relPt, mu3_relPt});

      var_tauMass=TauRefitLV.M();
      TauMass_all.at(t).Fill(TauRefitLV.M(),1);
      m3m = TauRefitLV.M();
      dataMCtype = id;
      event_weight = w;
      rapidity = TauLV.Eta();
      LumiScale=1;
      bdt = -99.0;


      //---------------- define per event resolution categroies 
      if (threeGlobal){
         //Category A1(3 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutLow ){
            TauMass_allVsBDTA_3glb.at(t).Fill(TauRefitLV.M(),readerA_3glb->EvaluateMVA("BDT"));
            BDTOutputA_3glb.at(t).Fill(    readerA_3glb->EvaluateMVA("BDT") );
            if(readerA_3glb->EvaluateMVA("BDT") > bdt_cutA1_3glb){
               TauMassRefitA1_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassA1_3glb.at(t).Fill(TauLV.M(),1);
               category=1;
               bdt = readerA_3glb->EvaluateMVA("BDT");
            }
         }

         //Category B1(3 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutLow && Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutHigh){
            TauMass_allVsBDTB_3glb.at(t).Fill(TauRefitLV.M(),readerB_3glb->EvaluateMVA("BDT"));
            BDTOutputB_3glb.at(t).Fill(    readerB_3glb->EvaluateMVA("BDT") );
            if(readerB_3glb->EvaluateMVA("BDT") > bdt_cutB1_3glb){
               TauMassRefitB1_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassB1_3glb.at(t).Fill(TauLV.M(),1);
               category =2 ;
               bdt = readerB_3glb->EvaluateMVA("BDT");
            }
         }

         //Category C1(3 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutHigh){
            TauMass_allVsBDTC_3glb.at(t).Fill(TauRefitLV.M(),readerC_3glb->EvaluateMVA("BDT"));
            BDTOutputC_3glb.at(t).Fill(    readerC_3glb->EvaluateMVA("BDT") );
            if(readerC_3glb->EvaluateMVA("BDT") > bdt_cutC1_3glb){
               TauMassRefitC1_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassC1_3glb.at(t).Fill(TauLV.M(),1);
               category = 3;
               bdt = readerC_3glb->EvaluateMVA("BDT");
            }
         }

         //Category A2(3 global)  
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutLow){
            if(readerA_3glb->EvaluateMVA("BDT") > bdt_cutA2_3glb && readerA_3glb->EvaluateMVA("BDT") <= bdt_cutA1_3glb){
               TauMassRefitA2_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassA2_3glb.at(t).Fill(TauLV.M(),1);
               category=4;
               bdt = readerA_3glb->EvaluateMVA("BDT");
            }
         }


         //Category B2(3 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutLow && Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutHigh){
            if(readerB_3glb->EvaluateMVA("BDT") > bdt_cutB2_3glb && readerB_3glb->EvaluateMVA("BDT") <= bdt_cutB1_3glb){
               TauMassRefitB2_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassB2_3glb.at(t).Fill(TauLV.M(),1);
               category =5 ;
               bdt = readerB_3glb->EvaluateMVA("BDT");
            }
         }

         //Category C2(3 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutHigh){
            if(readerC_3glb->EvaluateMVA("BDT") >  bdt_cutC2_3glb && readerC_3glb->EvaluateMVA("BDT")< bdt_cutC1_3glb){
               TauMassRefitC2_3glb.at(t).Fill(TauRefitLV.M(),1);    
               TauMassC2_3glb.at(t).Fill(TauLV.M(),1);
               category = 6;
               bdt = readerC_3glb->EvaluateMVA("BDT");
            }
         }
      }

      else{
         //Category A1(2 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutLow){
            TauMass_allVsBDTA_2glbtrk.at(t).Fill(TauRefitLV.M(),readerA_2glbTrk->EvaluateMVA("BDT"));
            BDTOutputA_2glbtrk.at(t).Fill(    readerA_2glbTrk->EvaluateMVA("BDT") );
            if(readerA_2glbTrk->EvaluateMVA("BDT") > bdt_cutA1_2glbTrk){
               TauMassRefitA1_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassA1_2glbtrk.at(t).Fill(TauLV.M(),1);
               category=7;
               bdt = readerA_2glbTrk->EvaluateMVA("BDT");
            }
         }

         //Category B1(2 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutLow && Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutHigh){
            TauMass_allVsBDTB_2glbtrk.at(t).Fill(TauRefitLV.M(),readerB_2glbTrk->EvaluateMVA("BDT"));
            BDTOutputB_2glbtrk.at(t).Fill(    readerB_2glbTrk->EvaluateMVA("BDT") );
            if(readerB_2glbTrk->EvaluateMVA("BDT") > bdt_cutB1_2glbTrk){
               TauMassRefitB1_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassB1_2glbtrk.at(t).Fill(TauLV.M(),1);
               category =8 ;
               bdt = readerB_2glbTrk->EvaluateMVA("BDT");
            }
         }

         //Category C1(2 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutHigh){
            TauMass_allVsBDTC_2glbtrk.at(t).Fill(TauRefitLV.M(),readerC_2glbTrk->EvaluateMVA("BDT"));
            BDTOutputC_2glbtrk.at(t).Fill(    readerC_2glbTrk->EvaluateMVA("BDT") );
            if(readerC_2glbTrk->EvaluateMVA("BDT") > bdt_cutC1_2glbTrk){
               TauMassRefitC1_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassC1_2glbtrk.at(t).Fill(TauLV.M(),1);
               category = 9;
               bdt = readerC_2glbTrk->EvaluateMVA("BDT");
            }
         }

         //Category A2(2 global)  
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutLow){
            if(readerA_2glbTrk->EvaluateMVA("BDT") > bdt_cutA2_2glbTrk && readerA_2glbTrk->EvaluateMVA("BDT") <= bdt_cutA1_2glbTrk){
               TauMassRefitA2_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassA2_2glbtrk.at(t).Fill(TauLV.M(),1);
               category=10;
               bdt = readerA_2glbTrk->EvaluateMVA("BDT");
            }
         }


         //Category B2(2 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutLow && Ntp->TauMassResolution(EtaSortedIndices,1,false) < tauMassResCutHigh){
            if(readerB_2glbTrk->EvaluateMVA("BDT") > bdt_cutB2_2glbTrk && readerB_2glbTrk->EvaluateMVA("BDT") <= bdt_cutB1_2glbTrk){
               TauMassRefitB2_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassB2_2glbtrk.at(t).Fill(TauLV.M(),1);
               category =11 ;
               bdt = readerB_2glbTrk->EvaluateMVA("BDT");
            }
         }

         //Category C2(2 global)
         if(Ntp->TauMassResolution(EtaSortedIndices,1,false) >= tauMassResCutHigh){
            if(readerC_2glbTrk->EvaluateMVA("BDT") >bdt_cutC2_2glbTrk && readerC_2glbTrk->EvaluateMVA("BDT")<= bdt_cutC1_2glbTrk){
               TauMassRefitC2_2glbtrk.at(t).Fill(TauRefitLV.M(),1);    
               TauMassC2_2glbtrk.at(t).Fill(TauLV.M(),1);
               category = 12;
               bdt = readerC_2glbTrk->EvaluateMVA("BDT");
            }
         }
      }
      T3MMiniTree->Fill();


      //---------------  Fill MC plots 
      if(id==40 || id == 60 || id ==90){
         if(Ntp->MCEventIsReconstructed()){
            TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0)));
            TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1)));
            TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2)));
            TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;

            if (threeGlobal){
               TauMassResolution_3glb.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
               TauMassResolutionRefit_3glb.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);
               Muon1DRToTruth_3glb.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
               Muon2DRToTruth_3glb.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
               Muon3DRToTruth_3glb.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);

            }
            else{
               TauMassResolution_2glbtrk.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
               TauMassResolutionRefit_2glbtrk.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);
               Muon1DRToTruth_2glbtrk.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
               Muon2DRToTruth_2glbtrk.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
               Muon3DRToTruth_2glbtrk.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
            }
         }
      }
   }
}


void  SignalSelector::Finish(){


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





