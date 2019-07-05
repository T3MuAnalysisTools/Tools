#include "MCStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

MCStudy::MCStudy(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.70),
  tauMaxMass_(1.82),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

MCStudy::~MCStudy(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  MCStudy::Configure(){
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=2.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=1.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=0.03;
    if(i==ThreeMuMass)        cut.at(ThreeMuMass)=1;// true for MC and mass side band for data
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==TriggerOk){
      title.at(i)="Pass HLT";
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
      title.at(i)="$p_{T}(\\mu_{3}) >$ 1 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Muon3 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
    }
    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="pass MuonID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PhiVeto){
      title.at(i)="$\\phi$ mass veto";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="Phi mass Veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,50,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,50,0.8,1.2,hlabel,"Events"));
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
    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==ThreeMuMass){
      title.at(i)="$\\tau$ mass (sideband in data)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuMass_",htitle,50,1.4,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuMass_",htitle,50,1.4,2.2,hlabel,"Events"));
    }


  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
  Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
  Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");


  Muon1isGlob =HConfig.GetTH1D(Name+"_Muon1isGlob","Muon1isGlob",2,-0.5,1.5,"  #mu_{1} is global muon","Events");
  Muon2isGlob =HConfig.GetTH1D(Name+"_Muon2isGlob","Muon2isGlob",2,-0.5,1.5,"  #mu_{2} is global muon","Events");
  Muon3isGlob =HConfig.GetTH1D(Name+"_Muon3isGlob","Muon3isGlob",2,-0.5,1.5,"  #mu_{3} is global muon","Events");

  Muon1isStand =HConfig.GetTH1D(Name+"_Muon1isStand","Muon1isStand",2,-0.5,1.5,"  #mu_{1} is a standalone muon","Events");
  Muon2isStand =HConfig.GetTH1D(Name+"_Muon2isStand","Muon2isStand",2,-0.5,1.5,"  #mu_{2} is a standalone muon","Events");
  Muon3isStand =HConfig.GetTH1D(Name+"_Muon3isStand","Muon3isStand",2,-0.5,1.5,"  #mu_{3} is a standalone muon","Events");

  Muon1isTrack =HConfig.GetTH1D(Name+"_Muon1isTrack","Muon1isTrack",2,-0.5,1.5,"  #mu_{1} is tracker muon","Events");
  Muon2isTrack =HConfig.GetTH1D(Name+"_Muon2isTrack","Muon2isTrack",2,-0.5,1.5,"  #mu_{2} is tracker muon","Events");
  Muon3isTrack =HConfig.GetTH1D(Name+"_Muon3isTrack","Muon3isTrack",2,-0.5,1.5,"  #mu_{3} is tracker muon","Events");




  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#eta(#mu_{1})","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#eta(#mu_{2})","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#eta(#mu_{3})","Events");

  Muon1Eta_EtaSort=HConfig.GetTH1D(Name+"_Muon1Eta_EtaSort","Muon1Eta_EtaSort",26,-2.6,2.6,"#eta(#mu_{1})( #eta sorting)","Events");
  Muon2Eta_EtaSort=HConfig.GetTH1D(Name+"_Muon2Eta_EtaSort","Muon2Eta_EtaSort",26,-2.6,2.6,"#eta(#mu_{2})( #eta sorting)","Events");
  Muon3Eta_EtaSort=HConfig.GetTH1D(Name+"_Muon3Eta_EtaSort","Muon3Eta_EtaSort",26,-2.6,2.6,"#eta(#mu_{3})( #eta sorting)","Events");


  Muon1Pt_EtaSort =HConfig.GetTH1D(Name+"_Muon1Pt_EtaSort","Muon1Pt_EtaSort",50,0,25,"  #mu_{1} p_{T}( #eta sorting), GeV","Events");
  Muon2Pt_EtaSort =HConfig.GetTH1D(Name+"_Muon2Pt_EtaSort","Muon2Pt_EtaSort",50,0,20,"  #mu_{2} p_{T}( #eta sorting), GeV","Events");
  Muon3Pt_EtaSort =HConfig.GetTH1D(Name+"_Muon3Pt_EtaSort","Muon3Pt_EtaSort",50,0,15,"  #mu_{3} p_{T}( #eta sorting), GeV","Events");

  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",100,-0.1,0.1,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",100,-0.1,0.1,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
  
  TauMassResolutionVsEta=HConfig.GetTH2D(Name+"_TauMassResolutionVsEta","TauMassResolutionVsEta",50,-2.5,2.5,50, -0.1,0.1,"#eta","#Delta M_{#tau}  (reco - mc)/mc");
  TauMassResolutionVsPt=HConfig.GetTH2D(Name+"_TauMassResolutionVsPt","TauMassResolutionVsPt",50,5,40,50,-0.1,0.1,"p_{T}, GeV","#Delta M_{#tau}  (reco - mc)/mc ");

  TauMassResolutionVsMu1Eta=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu1Eta","TauMassResolutionVsMu1Eta",50,-2.5,2.5,50, -0.1,0.1,"#mu_{1} #eta","#Delta M_{#tau}  (reco - mc)/mc");
  TauMassResolutionVsMu1Pt=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu1Pt","TauMassResolutionVsMu1Pt",50,5,25,50,-0.1,0.1,"#mu_{1} p_{T}, GeV","#Delta M_{#tau}  (reco - mc)/mc ");

  TauMassResolutionVsMu2Eta=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu2Eta","TauMassResolutionVsMu2Eta",50,-2.5,2.5,50, -0.1,0.1,"#mu_{2} #eta","#Delta M_{#tau}  (reco - mc)/mc");
  TauMassResolutionVsMu2Pt=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu2Pt","TauMassResolutionVsMu2Pt",50,5,20,50,-0.1,0.1,"#mu_{2} p_{T}, GeV","#Delta M_{#tau}  (reco - mc)/mc ");

  TauMassResolutionVsMu3Eta=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu3Eta","TauMassResolutionVsMu3Eta",50,-2.5,2.5,50, -0.1,0.1,"#mu_{3} #eta","#Delta M_{#tau}  (reco - mc)/mc");
  TauMassResolutionVsMu3Pt=HConfig.GetTH2D(Name+"_TauMassResolutionVsMu3Pt","TauMassResolutionVsMu3Pt",50,5,15,50,-0.1,0.1,"#mu_{3} p_{T}, GeV","#Delta M_{#tau}  (reco - mc)/mc ");


  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"reco - mc #mu_{3} #Delta R","Events");

  PETauMassResolution_Pt = HConfig.GetTH1D(Name+"_PETauMassResolution_Pt","PETauMassResolution_Pt",50,0,0.02,"#frac{#Delta m}{m} (pt)","Events");
  PETauMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_PETauMassResolution_PtEtaPhi","PETauMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
  
  PETauMassResolution_PtEtaPhi_RefitTracks = HConfig.GetTH1D(Name+"_PETauMassResolution_PtEtaPhi_RefitTracks",
							     "PETauMassResolution_PtEtaPhi_RefitTracks",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi) (refit tracks)","Events");

  PETauMassResVsMu1Pt=HConfig.GetTH2D(Name+"_PETauMassResVsMu1Pt","PETauMassResVsMu1Pt",50,0,25,50,0,0.02,"#mu_{1} p_{T}, GeV","#frac{#Delta m}{m} (ptEtaPhi)");
  PETauMassResVsMu2Pt=HConfig.GetTH2D(Name+"_PETauMassResVsMu2Pt","PETauMassResVsMu2Pt",50,0,20,50,0,0.02,"#mu_{2} p_{T}, GeV","#frac{#Delta m}{m} (ptEtaPhi)");
  PETauMassResVsMu3Pt=HConfig.GetTH2D(Name+"_PETauMassResVsMu3Pt","PETauMassResVsMu3Pt",50,0,15,50,0,0.02,"#mu_{3} p_{T}, GeV","#frac{#Delta m}{m} (ptEtaPhi)");



  Mu1PtvsPtError=HConfig.GetTH2D(Name+"_Mu1PtvsPtError","Mu1PtvsPtError",50,0,25,50,0.01,0.2,"#mu_{1} p_{T} (p_{T} sorted), GeV","#delta p_{T}, GeV");
  Mu2PtvsPtError=HConfig.GetTH2D(Name+"_Mu2PtvsPtError","Mu2PtvsPtError",50,0,20,50,0.01,0.2,"#mu_{2} p_{T} (p_{T} sorted), GeV","#delta p_{T}, GeV");
  Mu3PtvsPtError=HConfig.GetTH2D(Name+"_Mu3PtvsPtError","Mu3PtvsPtError",50,0,15,50,0.01,0.2,"#mu_{3} p_{T} (p_{T} sorted), GeV","#delta p_{T}, GeV");
  
  Mu1PtvsEta=HConfig.GetTH2D(Name+"_Mu1PtvsEta","Mu1PtvsEta",50,0,25,50,-2.5,2.5,"#mu_{1} p_{T} (p_{T} sorted), GeV","#eta");
  Mu2PtvsEta=HConfig.GetTH2D(Name+"_Mu2PtvsEta","Mu2PtvsEta",50,0,20,50,-2.5,2.5,"#mu_{2} p_{T} (p_{T} sorted), GeV","#eta");
  Mu3PtvsEta=HConfig.GetTH2D(Name+"_Mu3PtvsEta","Mu3PtvsEta",50,0,15,50,-2.5,2.5,"#mu_{3} p_{T} (p_{T} sorted), GeV","#eta");
  TauPtVsEta=HConfig.GetTH2D(Name+"_TauPtVsEta","TauPtVsEta",50,5,45,50,-2.5,2.5,"#tau p_{T}, GeV","#eta");




  Mu1EtavsEtaError=HConfig.GetTH2D(Name+"_Mu1EtavsEtaError","Mu1EtavsEtaError",50,-2.5,2.5,50,0,0.004,"#mu_{1} #eta (#eta sorted)","#delta#eta");
  Mu2EtavsEtaError=HConfig.GetTH2D(Name+"_Mu2EtavsEtaError","Mu2EtavsEtaError",50,-2.5,2.5,50,0,0.004,"#mu_{2} #eta (#eta sorted)","#delta#eta");
  Mu3EtavsPetError=HConfig.GetTH2D(Name+"_Mu3EtavsEtaError","Mu3EtavsEtaError",50,-2.5,2.5,50,0,0.004,"#mu_{3} #eta (#eta sorted)","#delta#eta");

  Mu1PhivsPhiError=HConfig.GetTH2D(Name+"_Mu1PhivsPhiError","Mu1PhivsPhiError",50,-3.14,3.14,50,0,0.003,"#mu_{1} #phi (p_{T} sorted)","#delta#phi, rad");
  Mu2PhivsPhiError=HConfig.GetTH2D(Name+"_Mu2PhivsPhiError","Mu2PhivsPhiError",50,-3.14,3.14,50,0,0.003,"#mu_{2} #phi (p_{T} sorted)","#delta#phi, rad");
  Mu3PhivsPhiError=HConfig.GetTH2D(Name+"_Mu3PhivsPhiError","Mu3PhivsPhiError",50,-3.14,3.14,50,0,0.003,"#mu_{3} #phi (p_{T} sorted)","#delta#phi, rad");


 

  PETauMassResVsMu1Eta=HConfig.GetTH2D(Name+"_PETauMassResVsMu1Eta","PETauMassResVsMu1Eta",50,-2.5,2.5,50,0.,0.02,"#mu_{1} #eta","#frac{#Delta m}{m} (ptEtaPhi)");
  PETauMassResVsMu2Eta=HConfig.GetTH2D(Name+"_PETauMassResVsMu2Eta","PETauMassResVsMu2Eta",50,-2.5,2.5,50,0.,0.02,"#mu_{2} #eta","#frac{#Delta m}{m} (ptEtaPhi)");
  PETauMassResVsMu3Eta=HConfig.GetTH2D(Name+"_PETauMassResVsMu3Eta","PETauMassResVsMu3Eta",50,-2.5,2.5,50,0.,0.02,"#mu_{3} #eta","#frac{#Delta m}{m} (ptEtaPhi)");

  TauMassResolutionVsPtEta = HConfig.GetTH3F(Name+"_TauMassResolutionVsPtEta","TauMassResolutionVsPtEta",50,0,25,50,-2.5,2.5,50,0.,0.02,"p_{T}(#tau)","#eta(#tau)","#frac{#Delta m}{m}");


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  MCStudy::Store_ExtraDist(){ 


  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  Extradist1d.push_back(&Muon1isGlob);
  Extradist1d.push_back(&Muon2isGlob);
  Extradist1d.push_back(&Muon3isGlob);


  Extradist1d.push_back(&Muon1isStand);
  Extradist1d.push_back(&Muon2isStand);
  Extradist1d.push_back(&Muon3isStand);
    

  Extradist1d.push_back(&Muon1isTrack);
  Extradist1d.push_back(&Muon2isTrack);
  Extradist1d.push_back(&Muon3isTrack);

  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);

  Extradist2d.push_back(&TauMassResolutionVsEta);
  Extradist2d.push_back(&TauMassResolutionVsPt);

  Extradist2d.push_back(&TauMassResolutionVsMu1Eta);
  Extradist2d.push_back(&TauMassResolutionVsMu1Pt);

  Extradist2d.push_back(&TauMassResolutionVsMu2Eta);
  Extradist2d.push_back(&TauMassResolutionVsMu2Pt);

  Extradist2d.push_back(&TauMassResolutionVsMu3Eta);
  Extradist2d.push_back(&TauMassResolutionVsMu3Pt);

  Extradist1d.push_back(&Muon1Eta_EtaSort);
  Extradist1d.push_back(&Muon2Eta_EtaSort);
  Extradist1d.push_back(&Muon3Eta_EtaSort);
  Extradist1d.push_back(&Muon1Pt_EtaSort);
  Extradist1d.push_back(&Muon2Pt_EtaSort);
  Extradist1d.push_back(&Muon3Pt_EtaSort);

  Extradist1d.push_back(&PETauMassResolution_Pt);
  Extradist1d.push_back(&PETauMassResolution_PtEtaPhi);
  Extradist1d.push_back(&PETauMassResolution_PtEtaPhi_RefitTracks);

  Extradist2d.push_back(&PETauMassResVsMu1Pt);
  Extradist2d.push_back(&PETauMassResVsMu2Pt);
  Extradist2d.push_back(&PETauMassResVsMu3Pt);

  Extradist2d.push_back(&PETauMassResVsMu1Eta);
  Extradist2d.push_back(&PETauMassResVsMu2Eta);
  Extradist2d.push_back(&PETauMassResVsMu3Eta);

  Extradist2d.push_back(&Mu1PtvsPtError);
  Extradist2d.push_back(&Mu2PtvsPtError);
  Extradist2d.push_back(&Mu3PtvsPtError);
  
  Extradist2d.push_back(&Mu1EtavsEtaError);
  Extradist2d.push_back(&Mu2EtavsEtaError);
  Extradist2d.push_back(&Mu3EtavsPetError);

  Extradist2d.push_back(&Mu1PhivsPhiError);
  Extradist2d.push_back(&Mu2PhivsPhiError);
  Extradist2d.push_back(&Mu3PhivsPhiError);


  Extradist2d.push_back(&Mu1PtvsEta);
  Extradist2d.push_back(&Mu2PtvsEta);
  Extradist2d.push_back(&Mu3PtvsEta);
  Extradist2d.push_back(&TauPtVsEta);

  Extradist3d.push_back(&TauMassResolutionVsPtEta);



}


void  MCStudy::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu"))) value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);
  }

  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));
  value.at(SignalCandidate)=0;
  unsigned int  final_idx=0;
  value.at(TriggerMatch)=0;


  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();
    unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
    //    value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
    //    			 Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
    //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
    //


    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
    			Ntp->Muon_isTrackerMuon(mu3_pt_idx));
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

    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();


    value.at(PhiVeto)   = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(final_idx)){
      value.at(TriggerMatch)+=i; 
    }
    
    value.at(ThreeMuMass) = TauLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) == cut.at(SignalCandidate));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
  pass.at(MuonID) =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = (value.at(TriggerMatch)  <  cut.at(TriggerMatch));
  pass.at(PhiVeto) = true;//(fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*PDG_Var::Phi_width());
  pass.at(OmegaVeto) = true;//(fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*PDG_Var::Omega_width());

  if(id!=1) pass.at(ThreeMuMass) = true;
  else  pass.at(ThreeMuMass) = ( (value.at(ThreeMuMass) > tauMinSideBand_ && value.at(ThreeMuMass) < tauMinMass_)  ||   (value.at(ThreeMuMass)> tauMaxMass_ && value.at(ThreeMuMass) < tauMaxSideBand_));



  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);

  if(status){

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);


    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

    
    TLorentzVector Muon1EtaLV = Ntp->Muon_P4(Muon_Eta_index_1);
    TLorentzVector Muon2EtaLV = Ntp->Muon_P4(Muon_Eta_index_2);
    TLorentzVector Muon3EtaLV = Ntp->Muon_P4(Muon_Eta_index_3);




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


    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);

    std::vector<TLorentzVector > list;
    list.push_back(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0));
    list.push_back(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1));
    list.push_back(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2));

    TLorentzVector Mu1Refit_EtaSorted = Ntp->MatchedLV(list,Muon_Eta_index_1);
    TLorentzVector Mu2Refit_EtaSorted = Ntp->MatchedLV(list,Muon_Eta_index_2);
    TLorentzVector Mu3Refit_EtaSorted = Ntp->MatchedLV(list,Muon_Eta_index_3);


    Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
    Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
    Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);

    Muon1Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Eta(),1);
    Muon2Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Eta(),1);
    Muon3Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Eta(),1);

    Muon1Eta_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_1).Eta(),1);
    Muon2Eta_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_2).Eta(),1);
    Muon3Eta_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_3).Eta(),1);

    Muon1Pt_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_1).Pt(),1);
    Muon2Pt_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_2).Pt(),1);
    Muon3Pt_EtaSort.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_3).Pt(),1);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);


    Mu1PtvsEta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(), Ntp->Muon_P4(Muon_index_1).Eta());
    Mu2PtvsEta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(), Ntp->Muon_P4(Muon_index_2).Eta());
    Mu3PtvsEta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(), Ntp->Muon_P4(Muon_index_3).Eta());
    TauPtVsEta.at(t).Fill(TauLV.Pt(),TauLV.Eta());


    Muon1isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_1),1);
    Muon2isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_2),1);
    Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);


    Muon1isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_1),w);
    Muon2isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_2),w);
    Muon3isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_3),w);


    Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
    Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
    Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);

    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));


    Mu1PtvsPtError.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(), Ntp->Muon_ptError(Muon_index_1),w);
    Mu2PtvsPtError.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(), Ntp->Muon_ptError(Muon_index_2),w);
    Mu3PtvsPtError.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(), Ntp->Muon_ptError(Muon_index_3),w);

    Mu1EtavsEtaError.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_1).Eta(), Ntp->Muon_etaError(Muon_Eta_index_1),w);
    Mu2EtavsEtaError.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_2).Eta(), Ntp->Muon_etaError(Muon_Eta_index_2),w);
    Mu3EtavsPetError.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_3).Eta(), Ntp->Muon_etaError(Muon_Eta_index_3),w);
    
    Mu1PhivsPhiError.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Phi(), Ntp->Muon_phiError(Muon_index_1),w);
    Mu2PhivsPhiError.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Phi(), Ntp->Muon_phiError(Muon_index_2),w);
    Mu3PhivsPhiError.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Phi(), Ntp->Muon_phiError(Muon_index_3),w);



    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){
	
	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;



	std::vector<unsigned int> EtaSortedIndices;
	
	EtaSortedIndices.push_back(Muon_Eta_index_1);
	EtaSortedIndices.push_back(Muon_Eta_index_2);
	EtaSortedIndices.push_back(Muon_Eta_index_3);
	

	PETauMassResolution_Pt.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,0,false),w);
	PETauMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);
	PETauMassResolution_PtEtaPhi_RefitTracks.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,true),w);


	TauMassResolutionVsPtEta.at(t).Fill(TauLV.Pt(), TauLV.Eta(), Ntp->TauMassResolution(EtaSortedIndices,1,false));

	PETauMassResVsMu1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),Ntp->TauMassResolution(EtaSortedIndices,1,false));
	PETauMassResVsMu2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),Ntp->TauMassResolution(EtaSortedIndices,1,false));
	PETauMassResVsMu3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),Ntp->TauMassResolution(EtaSortedIndices,1,false));

	PETauMassResVsMu1Eta.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_1).Eta(),Ntp->TauMassResolution(EtaSortedIndices,1,false));
	PETauMassResVsMu2Eta.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_2).Eta(),Ntp->TauMassResolution(EtaSortedIndices,1,false));
	PETauMassResVsMu3Eta.at(t).Fill(Ntp->Muon_P4(Muon_Eta_index_3).Eta(),Ntp->TauMassResolution(EtaSortedIndices,1,false));

	TauMassResolutionVsMu1Eta.at(t).Fill(MCMuon1LV.Eta(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());
	TauMassResolutionVsMu1Pt.at(t).Fill(MCMuon1LV.Pt(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());

	TauMassResolutionVsMu2Eta.at(t).Fill(MCMuon2LV.Eta(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());
	TauMassResolutionVsMu2Pt.at(t).Fill(MCMuon2LV.Pt(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());

	TauMassResolutionVsMu3Eta.at(t).Fill(MCMuon3LV.Eta(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());
	TauMassResolutionVsMu3Pt.at(t).Fill(MCMuon3LV.Pt(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());

	TauMassResolutionVsEta.at(t).Fill(MCTauLV.Eta(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());
	TauMassResolutionVsPt.at(t).Fill(MCTauLV.Pt(),(TauLV.M() - MCTauLV.M())/MCTauLV.M());

	TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }
    
  }
}


void  MCStudy::Finish(){

  if(mode == RECONSTRUCT){
    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/3;
      //      ScaleAllHistOfType(i,scale);
    }
  }
  Selection::Finish();
}





