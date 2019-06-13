#include "ThreeMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

ThreeMu::ThreeMu(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.70),
  tauMaxMass_(1.82),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

ThreeMu::~ThreeMu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ThreeMu::Configure(){
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
      hlabel="Omega mass veto, GeV";
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

  Muon1kink =HConfig.GetTH1D(Name+"_Muon1kink","Muon1kink",50,0.,50,"  #mu_{1} trk kink #chi^{2}","Events");
  Muon2kink =HConfig.GetTH1D(Name+"_Muon2kink","Muon2kink",50,0.,50,"  #mu_{2} trk kink #chi^{2}","Events");
  Muon3kink =HConfig.GetTH1D(Name+"_Muon3kink","Muon3kink",50,0.,50,"  #mu_{3} trk kink #chi^{2}","Events");
  MuonkinkMax =HConfig.GetTH1D(Name+"_MuonkinkMax","MuonkinkMax",50,0.,50,"  #mu max trk kink #chi^{2}","Events");


  Muon1glbkink =HConfig.GetTH1D(Name+"_Muon1glbkink","Muon1glbkink",50,0.,50,"  #mu_{1} glb kink #chi^{2}","Events");
  Muon2glbkink =HConfig.GetTH1D(Name+"_Muon2glbkink","Muon2glbkink",50,0.,50,"  #mu_{2} glb kink #chi^{2}","Events");
  Muon3glbkink =HConfig.GetTH1D(Name+"_Muon3glbkink","Muon3glbkink",50,0.,50,"  #mu_{3} glb kink #chi^{2}","Events");
  MuonglbkinkMax =HConfig.GetTH1D(Name+"_MuonglbkinkMax","MuonglbkinkMax",50,0.,50,"  #mu max glb kink #chi^{2}","Events");

  MuonkinkMin  =HConfig.GetTH1D(Name+"_MuonkinkMin","MuonkinkMin",50,0.,50,"  #mu min trk kink #chi^{2}","Events"); 
  MuonglbkinkMin  =HConfig.GetTH1D(Name+"_MuonglbkinkMin","MuonglbkinkMin",50,0.,50,"  #mu min glb kink #chi^{2}","Events"); 
  MuonglbkinkSum  =HConfig.GetTH1D(Name+"_MuonglbkinkSum","MuonglbkinkSum",50,0.,50," #sum  #mu glb kink #chi^{2}","Events");
  MuonkinkSum  =HConfig.GetTH1D(Name+"_MuonkinkSum","MuonkinkSum",50,0.,50," #sum  #mu  trk kink #chi^{2}","Events");

  Muon_segmentCompatibility_mu1  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu1","Muon_segmentCompatibility_mu1",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
  Muon_segmentCompatibility_mu2  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu2","Muon_segmentCompatibility_mu2",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
  Muon_segmentCompatibility_mu3  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu3","Muon_segmentCompatibility_mu3",50,0.,1,"Inner Track and muon segment match  #mu_{3} ","Events");

  Muon1InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon1InOutTrackMatch","Muon1InOutTrackMatch",50,0.,10,"  #mu_{1} inner and outer tracker match","Events");
  Muon2InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon2InOutTrackMatch","Muon2InOutTrackMatch",50,0.,10,"  #mu_{2} inner and outer tracker match","Events");
  Muon3InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon3InOutTrackMatch","Muon3InOutTrackMatch",50,0.,10,"  #mu_{3} inner and outer tracker match","Events");
  MuonInOutTrackMatchMax =HConfig.GetTH1D(Name+"_MuonInOutTrackMatchMax","MuonInOutTrackMatchMax",50,0.,10,"  #mu max inner and outer tracker match","Events");


  Muon1_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_Muon1_globalDeltaEtaPhi","Muon1_globalDeltaEtaPhi",40,0,0.2,"#mu_{1} #Delta #eta-#phi","Events");
  Muon2_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_Muon2_globalDeltaEtaPhi","Muon2_globalDeltaEtaPhi",40,0,0.2,"#mu_{2} #Delta #eta-#phi","Events");
  Muon3_globalDeltaEtaPhi=HConfig.GetTH1D(Name+"_Muon2_globalDeltaEtaPhi","Muon3_globalDeltaEtaPhi",40,0,0.2,"#mu_{3} #Delta #eta-#phi","Events");

  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#eta(#mu_{1})","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#eta(#mu_{2})","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#eta(#mu_{3})","Events");

  Muon1StandardSelector=HConfig.GetTH1D(Name+"_Muon1StandardSelector","Muon1StandardSelector",23,-0.5,22.5,"#mu_{1} standard selector; bin 0 - no ID","Events");
  Muon2StandardSelector=HConfig.GetTH1D(Name+"_Muon2StandardSelector","Muon2StandardSelector",23,-0.5,22.5,"#mu_{2} standard selector; bin 0 - no ID","Events");
  Muon3StandardSelector=HConfig.GetTH1D(Name+"_Muon3StandardSelector","Muon3StandardSelector",23,-0.5,22.5,"#mu_{3} standard selector; bin 0 - no ID","Events");

  Muon1PtResolution=HConfig.GetTH1D(Name+"_Muon1PtResolution","Muon1PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{1})  (reco - mc)/mc ","Events");
  Muon2PtResolution=HConfig.GetTH1D(Name+"_Muon2PtResolution","Muon2PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{2})  (reco - mc)/mc ","Events");
  Muon3PtResolution=HConfig.GetTH1D(Name+"_Muon3PtResolution","Muon3PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{3})  (reco - mc)/mc  ","Events");

  Muon1EtaResolution=HConfig.GetTH1D(Name+"_Muon1EtaResolution","Muon1EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc  ","Events");
  Muon2EtaResolution=HConfig.GetTH1D(Name+"_Muon2EtaResolution","Muon2EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
  Muon3EtaResolution=HConfig.GetTH1D(Name+"_Muon3EtaResolution","Muon3EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
  TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");
  TauMass =HConfig.GetTH1D(Name+"_TauMass","#tau lepton mass",50,1.5,1.9,"  M_{#tau} , GeV","Events");
  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassRefit =HConfig.GetTH1D(Name+"_TauMassRefit","Refit #tau lepton mass",50,1.5,1.9,"KF refit  M_{#tau} , GeV","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  VertexDCA12=HConfig.GetTH1D(Name+"_VertexDCA12","VertexDCA12",40,0,0.10,"Closest distance (#mu_{1}#mu_{2})","Events");
  VertexDCA23=HConfig.GetTH1D(Name+"_VertexDCA23","VertexDCA23",40,0,0.10,"Closest distance (#mu_{2}#mu_{3})","Events");
  VertexDCA31=HConfig.GetTH1D(Name+"_VertexDCA31","VertexDCA31",40,0,0.10,"Closest distance (#mu_{3}#mu_{1})","Events");
  VertexDCAMax=HConfig.GetTH1D(Name+"_VertexDCAMax","VertexDCAMax",40,0,0.15,"Max closest distance between muons","Events");
  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  VertexChi2AF=HConfig.GetTH1D(Name+"_VertexChi2AF","VertexChi2AF",50,0,10,"AF vertex #chi^{2}","Events");

  FLSignificance2=HConfig.GetTH1D(Name+"_FLSignificance2","FLSignificance2",50,0,50,"PV - SV distance significance (squared)","Events");
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  VertexSignalKFRefittedMu1P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1P","VertexSignalKFRefittedMu1P",50,0,20,"KF refitted #mu_{1} track p (GeV)","Events");
  VertexSignalKFRefittedMu1Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Pt","VertexSignalKFRefittedMu1P",50,0,20,"KF refitted #mu_{1} track p_{T} (GeV)","Events");
  VertexSignalKFRefittedMu1Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Eta","VertexSignalKFRefittedMu1Eta",25,-2.5,2.5,"KF refitted #mu_{1} track #eta","Events");
  VertexSignalKFRefittedMu1Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Phi","VertexSignalKFRefittedMu1Phi",32,-3.2,3.2,"KF refitted #mu_{1} track #phi","Events");

  VertexSignalKFRefittedMu2P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2P","VertexSignalKFRefittedMu2P",50,0,20,"KF refitted #mu_{2} track p (GeV)","Events");
  VertexSignalKFRefittedMu2Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Pt","VertexSignalKFRefittedMu2Pt",50,0,20,"KF refitted #mu_{2} track p_{T} (GeV)","Events");
  VertexSignalKFRefittedMu2Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Eta","VertexSignalKFRefittedMu2Eta",25,-2.5,2.5,"KF refitted #mu_{2} track #eta","Events");
  VertexSignalKFRefittedMu2Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Phi","VertexSignalKFRefittedMu2Phi",32,-3.2,3.2,"KF refitted #mu_{2} track #phi","Events");

  VertexSignalKFRefittedMu3P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3P","VertexSignalKFRefittedMu3P",50,0,20,"KF refitted Mu3 p (GeV)","Events");
  VertexSignalKFRefittedMu3Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Pt","VertexSignalKFRefittedMu3Pt",50,0,20,"KF refitted Mu3 p_{T} (GeV)","Events");
  VertexSignalKFRefittedMu3Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Eta","VertexSignalKFRefittedMu3Eta",25,-2.5,2.5,"KF refitted Mu3 #eta","Events");
  VertexSignalKFRefittedMu3Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Phi","VertexSignalKFRefittedMu3Phi",32,-3.2,3.2,"KF refitted Mu3 #phi","Events");

  VertexMu1D0Reco=HConfig.GetTH1D(Name+"_VertexMu1D0Reco","VertexMu1D0Reco",40,0,0.01,"#mu_{1} - PV transverse distance","Events");
  VertexMu1D0SigReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigReco","VertexMu1D0SigReco",40,0,3,"#mu_{1} - PV transverse distance significance","Events");
  VertexMu2D0Reco=HConfig.GetTH1D(Name+"_VertexMu2D0Reco","VertexMu2D0Reco",40,0,0.01,"#mu_{2} - PV transverse distance","Events");
  VertexMu2D0SigReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigReco","VertexMu2D0SigReco",40,0,3,"#mu_{2} - PV transverse distance significance","Events");
  VertexMu3D0Reco=HConfig.GetTH1D(Name+"_VertexMu3D0Reco","VertexMu3D0Reco",40,0,0.01,"#mu_{3} - PV transverse distance","Events");
  VertexMu3D0SigReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigReco","VertexMu3D0SigReco",50,0,3,"#mu_{3} - PV transverse distance significance","Events");
  VertexMuMinD0SigReco=HConfig.GetTH1D(Name+"_VertexMuMinD0SigReco","VertexMuMinD0SigReco",50,0,3,"#mu - PV min transverse distance significance","Events");
  VertexMuMaxD0SigReco=HConfig.GetTH1D(Name+"_VertexMuMaxD0SigReco","VertexMuMaxD0SigReco",50,0,4,"#mu - PV max transverse distance significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}","Events");

  Vertex2DDisplacementSignificance=HConfig.GetTH1D(Name+"_Vertex2DDisplacementSignificance","Vertex2DDisplacementSignificance",20,0,4,"PV-SV  transverse distance  significance","Events");
  Vertex3DDisplacementSignificance=HConfig.GetTH1D(Name+"_Vertex3DDisplacementSignificance","Vertex3DDisplacementSignificance",20,0,4,"PV-SV  distance significance","Events");

  Vertex2DDisplacement=HConfig.GetTH1D(Name+"_Vertex2DDisplacement","Vertex2DDisplacement",10,0,0.1,"PV-SV transverse displacement","Events");
  Vertex3DDisplacement=HConfig.GetTH1D(Name+"_Vertex3DDisplacement","Vertex3DDisplacement",10,0,0.1,"PV-SV displacement","Events");

  VertexPairQuality=HConfig.GetTH1D(Name+"_VertexPairQuality","VertexPairQuality",10,0,10,"vertex pair quality","Events");
  VertexPairfitStatus=HConfig.GetTH1D(Name+"_VertexPairfitStatus","VertexPairfitStatus",2,0,2,"vertex pair fit status","Events");

  VertexSignalKFChi2=HConfig.GetTH1D(Name+"_VertexSignalKFChi2","VertexSignalKFChi2",40,0,20,"vertex KF #chi^{2}","Events");
  VertexSignalAFX=HConfig.GetTH1D(Name+"_VertexSignalAFX","VertexSignalAFX",100,0,10,"AF vertex x (cm)","Events");
  VertexSignalAFY=HConfig.GetTH1D(Name+"_VertexSignalAFY","VertexSignalAFY",100,0,10,"AF vertex y (cm)","Events");
  VertexSignalAFZ=HConfig.GetTH1D(Name+"_VertexSignalAFZ","VertexSignalAFZ",100,0,10,"AF vertex z (cm)","Events");
  VertexSignalAFChi2=HConfig.GetTH1D(Name+"_VertexSignalChi2","VertexSignalChi2",100,0,10,"AF vertex #chi^{2}","Events");
  VertexSignalAFNdf=HConfig.GetTH1D(Name+"_VertexSignalAFNdf","VertexSignalAFNdf",10,0,10,"AF vertex ndf","Events");

  VertexMatchedPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexX","VertexMatchedPrimaryVertexX",50,-0.1,0.1,"vertex matched pv x (cm)","Events");
  VertexMatchedPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexY","VertexMatchedPrimaryVertexY",50,-0.1,0.1,"vertex matched pv y (cm)","Events");
  VertexMatchedPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexZ","VertexMatchedPrimaryVertexZ",50,-0.5,0.5,"vertex matched pv z (cm)","Events");
  VertexRefitPVisValid=HConfig.GetTH1D(Name+"_VertexRefitPVisValid","VertexRefitPVisValid",2,0,2,"vertex refit pv is valid","Events");
  VertexMatchedRefitPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexX","VertexMatchedRefitPrimaryVertexX",50,-0.1,0.1,"vertex matched refit pv x (cm)","Events");
  VertexMatchedRefitPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexY","VertexMatchedRefitPrimaryVertexY",50,-0.1,0.1,"vertex matched refit pv y (cm)","Events");
  VertexMatchedRefitPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexZ","VertexMatchedRefitPrimaryVertexZ",50,-0.1,0.1,"vertex matched refit pv z (cm)","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");

  MuPair1_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair1_vertex_chi2","MuPair1_vertex_chi2",50,0,3,"KF  #chi^{2} of first #mu pair","Events");
  MuPair2_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair2_vertex_chi2","MuPair2_vertex_chi2",50,0,3,"KF  #chi^{2} of second #mu pair","Events");
  MuPair3_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair3_vertex_chi2","MuPair3_vertex_chi2",50,0,3,"KF  #chi^{2} of third #mu pair","Events");
  MuPair_vertex_chi2_min=HConfig.GetTH1D(Name+"_MuPair_vertex_chi2_min","MuPair_vertex_chi2_min",50,0,1.5,"KF min #chi^{2} of #mu pair","Events");

  Pair1Mass =HConfig.GetTH1D(Name+"_Pair1Mass","Pair1Mass",50,0,2," mass of #mu pair (p_{T} sorted  1-2), GeV","Events");
  Pair2Mass =HConfig.GetTH1D(Name+"_Pair2Mass","Pair2Mass",50,0,2," mass of #mu pair (p_{T} sorted  2-3), GeV","Events");
  Pair3Mass =HConfig.GetTH1D(Name+"_Pair3Mass","Pair3Mass",50,0,2," mass of #mu pair (p_{T} sorted  1-3), GeV","Events");

  Pair1Mass_OS1  =HConfig.GetTH1D(Name+"_Pair1Mass_OS1","Pair1Mass_OS1",50,0,2," mass of OS #mu#mu pair, GeV","Events");
  Pair2Mass_OS2  =HConfig.GetTH1D(Name+"_Pair2Mass_OS2","Pair2Mass_OS2",50,0,2," mass of #mu#mu pair SS1, GeV","Events");
  Pair3Mass_SS   =HConfig.GetTH1D(Name+"_Pair3Mass_SS","Pair3Mass_SS",50,0,2," mass of #mu#mu pair SS2,  GeV","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");

  dR12 = HConfig.GetTH1D(Name+"_dR12","dR12",50,0,0.6,"#DeltaR(#mu_{1}#mu_{2})","Events");
  dR23 = HConfig.GetTH1D(Name+"_dR23","dR23",50,0,0.6,"#DeltaR(#mu_{2}#mu_{3})","Events");
  dR31 = HConfig.GetTH1D(Name+"_dR31","dR31",50,0,0.6,"#DeltaR(#mu_{3}#mu_{1})","Events");

  maxMudR = HConfig.GetTH1D(Name+"_maxMudR","maxMudR",50,0,0.6,"max #DeltaR(#mu#mu)","Events");
  minMudR = HConfig.GetTH1D(Name+"_minMudR","minMudR",50,0,0.6,"min #DeltaR(#mu#mu)","Events");

  Mu1TauPTRatio = HConfig.GetTH1D(Name+"_Mu1TauPTRatio","Mu1TauPTRatio",50,0,1,"p_{T}(#mu_{1})/p_{T}(#tau)","Events");
  Mu2TauPTRatio = HConfig.GetTH1D(Name+"_Mu2TauPTRatio","Mu2TauPTRatio",50,0,1,"p_{T}(#mu_{2})/p_{T}(#tau)","Events");
  Mu3TauPTRatio = HConfig.GetTH1D(Name+"_Mu3TauPTRatio","Mu3TauPTRatio",50,0,1,"p_{T}(#mu_{3})/p_{T}(#tau)","Events");

  Mu1TauPRatio = HConfig.GetTH1D(Name+"_Mu1TauPRatio","Mu1TauPRatio",50,0,1,"p(#mu_{1})/p(#tau)","Events");
  Mu2TauPRatio = HConfig.GetTH1D(Name+"_Mu2TauPRatio","Mu2TauPRatio",50,0,1,"p(#mu_{2})/p(#tau)","Events");
  Mu3TauPRatio = HConfig.GetTH1D(Name+"_Mu3TauPRatio","Mu3TauPRatio",50,0,1,"p(#mu_{3})/p(#tau)","Events");

  dR1Tau = HConfig.GetTH1D(Name+"_dR1Tau","dR1Tau",50,0,0.5,"#DeltaR(#mu_{1}#tau)","Events");
  dR2Tau = HConfig.GetTH1D(Name+"_dR2Tau","dR2Tau",50,0,0.5,"#DeltaR(#mu_{2}#tau)","Events");
  dR3Tau = HConfig.GetTH1D(Name+"_dR3Tau","dR3Tau",50,0,0.5,"#DeltaR(#mu_{3}#tau)","Events");
  dRMaxMuTau = HConfig.GetTH1D(Name+"_dRMaxMuTau","dRMaxMuTau",50,0,0.5,"max #DeltaR(#mu#tau)","Events");
  dRMinMuTau = HConfig.GetTH1D(Name+"_dRMinMuTau","dRMinMuTau",50,0,0.5,"min #DeltaR(#mu#tau)","Events");

  Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","Isolation_NTracks",10,-0.5,9.5,"N tracks","Events");
  Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","Isolation_RelPt",50,0,1,"relative p_{T}","Events");
  Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","Isolation_MinDist",50,0,0.1,"Iso MinDist","Events");
  Isolation05_RelPt=HConfig.GetTH1D(Name+"_Isolation05_RelPt","Isolation05_RelPt",50,0,1,"relative  rel p_{T} in 0.5 cone","Events");
  Isolation05_NTracks=HConfig.GetTH1D(Name+"_Isolation05_NTracks","Isolation05_NTracks",10,-0.5,9.5,"N tracks in 0.5 cone","Events");
  Isolation05_MinDist=HConfig.GetTH1D(Name+"_Isolation05_MinDist","Isolation05_MinDist",50,0,0.1,"Iso05 MinDist","Events");

  Isolation_Ntrk1=HConfig.GetTH1D(Name+"_Isolation_Ntrk1","Isolation_Ntrk1",20,-0.5,19.5,"N tracks wihtin a 0.3 cone","Events");
  Isolation_Ntrk2=HConfig.GetTH1D(Name+"_Isolation_Ntrk2","Isolation_Ntrk2",20,-0.5,19.5,"N tracks wihtin a 0.3 cone","Events");
  Isolation_Ntrk3=HConfig.GetTH1D(Name+"_Isolation_Ntrk3","Isolation_Ntrk3",20,-0.5,19.5,"N tracks wihtin a 0.3 cone","Events");
  Isolation_Ntrk0p1=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p1","Isolation_Ntrk0p1",20,-0.5,19.5,"N tracks transverse displaced from PV > 0.1","Events");
  Isolation_Ntrk0p2=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p2","Isolation_Ntrk0p2",20,-0.5,19.5,"N tracks transverse displaced from PV > 0.2 ","Events");
  Isolation_Ntrk0p5=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p5","Isolation_Ntrk0p5",20,-0.5,19.5,"N tracks transverse displaced from PV > 0.5","Events");

  Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","Isolation_maxdxy",40,0,15,"Iso maximum transversely displaced track","Events");

  Isolation_RelativePtMu1=HConfig.GetTH1D(Name+"_Isolation_RelativePtMu1","Isolation_RelativePtMu1",50,0,5,"p_{T} ratio  track/#mu_{1} (0.3 cone)","Events");
  Isolation_RelativePtMu2=HConfig.GetTH1D(Name+"_Isolation_RelativePtMu2","Isolation_RelativePtMu2",50,0,5,"p_{T} ratio  track/#mu_{2} (0.3 cone)","Events");
  Isolation_RelativePtMu3=HConfig.GetTH1D(Name+"_Isolation_RelativePtMu3","Isolation_RelativePtMu3",50,0,5,"p_{T} ratio  track/#mu_{3} (0.3 cone)","Events");
  Isolation_RelativePtMaxMu=HConfig.GetTH1D(Name+"_Isolation_RelativePtMaxMu","Isolation_RelativePtMaxMu",50,0,5,"p_{T} ratio  track/#mu Max (0.3 cone)","Events");


  Isolation_Muon_hadVetoEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_hadVetoEt03","Isolation_Muon_hadVetoEt03",50,0,10,"Leading #mu iso HCAL E_{T} Veto (0.3 cone)","Events");
  Isolation_Muon_hadEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_hadEt03","Isolation_Muon_hadEt03",50,0,10,"Leading #mu iso HCAL E_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_emVetoEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_emVetoEt03","Isolation_Muon_emVetoEt03",50,0,4,"Leading #mu iso ECAL E_{T} Veto (0.3 cone)","Events");
  Isolation_Muon_emEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_emEt03","Isolation_Muon_emEt03",50,0,10,"Leading #mu iso ECAL E_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_nJets03   =HConfig.GetTH1D(Name+"_Isolation_Muon_nJets03","Isolation_Muon_nJets03",4,-0.5,3.5,"Leading #mu  nJets (0.3 cone)","Events");
  Isolation_Muon_nTracks03   =HConfig.GetTH1D(Name+"_Isolation_Muon_nTracks03","Isolation_Muon_nTracks03",20,-0.5,19.5,"Leading #mu  nTracks (0.3 cone)","Events");
  Isolation_Muon_sumPt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPt03","Isolation_Muon_sumPt03",50,0,25,"Leading #mu #sum tracks  p_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_trackerVetoPt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_trackerVetoPt03","Isolation_Muon_trackerVetoPt03",50,3,20,"Leading #mu tracker veto p_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_hadVetoEt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_hadVetoEt05","Isolation_Muon_hadVetoEt05",50,0,10,"Leading #mu iso HCAL E_{T} Veto (0.5 cone)","Events");
  Isolation_Muon_hadEt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_hadEt05","Isolation_Muon_hadEt05",50,0,20,"Leading #mu iso HCAL E_{T} (0.5 cone), GeV","Events");
  Isolation_Muon_emVetoEt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_emVetoEt05","Isolation_Muon_emVetoEt05",50,0,4,"Leading #mu iso ECAL E_{T} Veto (0.5 cone)","Events");
  Isolation_Muon_emEt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_emEt05","Isolation_Muon_emEt05",50,0,15,"Leading #mu iso ECAL E_{T} (0.5 cone), GeV","Events");
  Isolation_Muon_nJets05   =HConfig.GetTH1D(Name+"_Isolation_Muon_nJets05","Isolation_Muon_nJets05",4,-0.5,3.5,"Leading #mu  nJets (0.5 cone)","Events");
  Isolation_Muon_nTracks05   =HConfig.GetTH1D(Name+"_Isolation_Muon_nTracks05","Isolation_Muon_nTracks05",20,-0.5,19.5,"Leading #mu  nTracks (0.5 cone)","Events");
  Isolation_Muon_sumPt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPt05","Isolation_Muon_sumPt05",50,4,25,"Leading #mu #sum tracks  p_{T} (0.5 cone), GeV","Events");
  Isolation_Muon_trackerVetoPt05   =HConfig.GetTH1D(Name+"_Isolation_Muon_trackerVetoPt05","Isolation_Muon_trackerVetoPt05",50,3,20,"Leading #mu tracker veto p_{T} (0.5 cone), GeV","Events");
  Isolation_Muon_sumChargedHadronPt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumChargedHadronPt03","Isolation_Muon_sumChargedHadronPt03",50,0,10,"Leading #mu #sum ChHad p_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_sumChargedParticlePt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumChargedParticlePt03","Isolation_Muon_sumChargedParticlePt03",50,0,25,"Leading #mu #sum ChParticle p_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_sumNeutralHadronEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumNeutralHadronEt03","Isolation_Muon_sumNeutralHadronEt03",50,0,10,"Leading #mu #sum Neutral Hadron E_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_sumPhotonEt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPhotonEt03","Isolation_Muon_sumPhotonEt03",50,0,5,"Leading #mu #sum Photon E_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_sumPUPt03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPUPt03","Isolation_Muon_sumPUPt03",50,0,10,"Leading #mu #sum PU  p_{T} (0.3 cone), GeV","Events");
  Isolation_Muon_sumNeutralHadronEtHighThreshold03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumNeutralHadronEtHighThreshold03","Isolation_Muon_sumNeutralHadronEtHighThreshold03",50,0,6,"Lead #mu #sum Neutral Had E_{T}(High Threshold) (0.3 cone), GeV","Events");
  Isolation_Muon_sumPhotonEtHighThreshold03   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPhotonEtHighThreshold03","Isolation_Muon_sumPhotonEtHighThreshold03",50,0,5,"Leading #mu #sum Photon E_{T} (High Threshold) (0.3 cone), GeV","Events");
  
  Isolation_Muon_sumChargedHadronPt04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumChargedHadronPt04","Isolation_Muon_sumChargedHadronPt04",50,0,10,"Leading #mu #sum ChHad p_{T} (0.4 cone), GeV","Events");
  Isolation_Muon_sumChargedParticlePt04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumChargedParticlePt04","Isolation_Muon_sumChargedParticlePt04",50,5,35,"Leading #mu #sum ChParticle p_{T} (0.4 cone), GeV","Events");
  Isolation_Muon_sumNeutralHadronEt04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumNeutralHadronEt04","Isolation_Muon_sumNeutralHadronEt04",50,0,6,"Leading #mu #sum Neutral Had E_{T} (0.4 cone), GeV","Events");
  Isolation_Muon_sumPhotonEt04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPhotonEt04","Isolation_Muon_sumPhotonEt04",50,0,5,"Leading #mu #sum Photon E_{T} (0.4 cone), GeV","Events");
  Isolation_Muon_sumPUPt04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPUPt04","Isolation_Muon_sumPUPt04",50,0,20,"Leading #mu #sum PU  p_{T} (0.4 cone), GeV","Events");
  Isolation_Muon_sumNeutralHadronEtHighThreshold04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumNeutralHadronEtHighThreshold04","Isolation_Muon_sumNeutralHadronEtHighThreshold04",50,0,6,"Lead #mu #sum Neutral Had E_{T}(High Threshold) (0.4 cone), GeV","Events");
  Isolation_Muon_sumPhotonEtHighThreshold04   =HConfig.GetTH1D(Name+"_Isolation_Muon_sumPhotonEtHighThreshold04","Isolation_Muon_sumPhotonEtHighThreshold04",50,0,5,"Leading #mu #sum Photon E_{T} (High Threshold) (0.4 cone), GeV","Events");
  




  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  ThreeMu::Store_ExtraDist(){ 


  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  Extradist1d.push_back(&Muon1StandardSelector);
  Extradist1d.push_back(&Muon2StandardSelector);
  Extradist1d.push_back(&Muon3StandardSelector);

  Extradist1d.push_back(&Muon1isGlob);
  Extradist1d.push_back(&Muon2isGlob);
  Extradist1d.push_back(&Muon3isGlob);


  Extradist1d.push_back(&Muon1isStand);
  Extradist1d.push_back(&Muon2isStand);
  Extradist1d.push_back(&Muon3isStand);
    

  Extradist1d.push_back(&Muon1isTrack);
  Extradist1d.push_back(&Muon2isTrack);
  Extradist1d.push_back(&Muon3isTrack);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauP);
  Extradist1d.push_back(&TauMass);
  Extradist1d.push_back(&TauMassRefit);
  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist1d.push_back(&SVPVTauDirAngle);

  Extradist1d.push_back(&Muon1kink);
  Extradist1d.push_back(&Muon2kink);
  Extradist1d.push_back(&Muon3kink);

  Extradist1d.push_back(&MuonkinkMax);


  Extradist1d.push_back(&Muon1glbkink);
  Extradist1d.push_back(&Muon2glbkink);
  Extradist1d.push_back(&Muon3glbkink);
  Extradist1d.push_back(&MuonglbkinkMax);

  Extradist1d.push_back(&Muon1InOutTrackMatch);
  Extradist1d.push_back(&Muon2InOutTrackMatch);
  Extradist1d.push_back(&Muon3InOutTrackMatch);

  Extradist1d.push_back(&MuonInOutTrackMatchMax);

  Extradist1d.push_back(&maxMudR);
  Extradist1d.push_back(&minMudR);

  Extradist1d.push_back(&Mu1TauPTRatio);
  Extradist1d.push_back(&Mu2TauPTRatio);
  Extradist1d.push_back(&Mu3TauPTRatio);

  Extradist1d.push_back(&Mu1TauPRatio);
  Extradist1d.push_back(&Mu2TauPRatio);
  Extradist1d.push_back(&Mu3TauPRatio);

  Extradist1d.push_back(&dRMaxMuTau);
  Extradist1d.push_back(&dRMinMuTau);

  Extradist1d.push_back(&Muon1PtResolution);
  Extradist1d.push_back(&Muon2PtResolution);
  Extradist1d.push_back(&Muon3PtResolution);

  Extradist1d.push_back(&Muon1EtaResolution);
  Extradist1d.push_back(&Muon2EtaResolution);
  Extradist1d.push_back(&Muon3EtaResolution);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);

  Extradist1d.push_back(&MuPair1_vertex_chi2);
  Extradist1d.push_back(&MuPair2_vertex_chi2);
  Extradist1d.push_back(&MuPair3_vertex_chi2);
  Extradist1d.push_back(&MuPair_vertex_chi2_min);

  Extradist1d.push_back(&Pair1Mass);
  Extradist1d.push_back(&Pair2Mass);
  Extradist1d.push_back(&Pair3Mass);

  Extradist1d.push_back(&Pair1Mass_OS1);
  Extradist1d.push_back(&Pair2Mass_OS2);
  Extradist1d.push_back(&Pair3Mass_SS);

  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);

  Extradist1d.push_back(&dR12);
  Extradist1d.push_back(&dR23);
  Extradist1d.push_back(&dR31);
  Extradist1d.push_back(&dR1Tau);
  Extradist1d.push_back(&dR2Tau);
  Extradist1d.push_back(&dR3Tau);

  Extradist1d.push_back(&FLSignificance2);
  Extradist1d.push_back(&FLSignificance);

  Extradist1d.push_back(&Isolation_NTracks);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_MinDist);
  Extradist1d.push_back(&Isolation05_RelPt);
  Extradist1d.push_back(&Isolation05_NTracks);
  Extradist1d.push_back(&Isolation05_MinDist);
  Extradist1d.push_back(&Isolation_Ntrk1);
  //  Extradist1d.push_back(&Isolation_Ntrk2);
  //  Extradist1d.push_back(&Isolation_Ntrk3);
  Extradist1d.push_back(&Isolation_Ntrk0p1);
  Extradist1d.push_back(&Isolation_Ntrk0p2);
  Extradist1d.push_back(&Isolation_Ntrk0p5);
  Extradist1d.push_back(&Isolation_maxdxy);

  Extradist1d.push_back(&Isolation_RelativePtMu1);
  Extradist1d.push_back(&Isolation_RelativePtMu2);
  Extradist1d.push_back(&Isolation_RelativePtMu3);
  Extradist1d.push_back(&Isolation_RelativePtMaxMu);

  Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&VertexChi2AF);
  Extradist1d.push_back(&VertexDCA12);
  Extradist1d.push_back(&VertexDCA23);
  Extradist1d.push_back(&VertexDCA31);
  Extradist1d.push_back(&VertexDCAMax);
  Extradist1d.push_back(&VertexSignalKFRefittedMu1P);
  Extradist1d.push_back(&VertexSignalKFRefittedMu1Pt);
  Extradist1d.push_back(&VertexSignalKFRefittedMu1Eta);
  Extradist1d.push_back(&VertexSignalKFRefittedMu1Phi);
  Extradist1d.push_back(&VertexSignalKFRefittedMu2P);
  Extradist1d.push_back(&VertexSignalKFRefittedMu2Pt);
  Extradist1d.push_back(&VertexSignalKFRefittedMu2Eta);
  Extradist1d.push_back(&VertexSignalKFRefittedMu2Phi);
  Extradist1d.push_back(&VertexSignalKFRefittedMu3P);
  Extradist1d.push_back(&VertexSignalKFRefittedMu3Pt);
  Extradist1d.push_back(&VertexSignalKFRefittedMu3Eta);
  Extradist1d.push_back(&VertexSignalKFRefittedMu3Phi);

  Extradist1d.push_back(&VertexMu1D0Reco);
  Extradist1d.push_back(&VertexMu1D0SigReco);
  Extradist1d.push_back(&VertexMu2D0Reco);
  Extradist1d.push_back(&VertexMu2D0SigReco);
  Extradist1d.push_back(&VertexMu3D0Reco);
  Extradist1d.push_back(&VertexMu3D0SigReco);

  Extradist1d.push_back(&VertexMuMinD0SigReco);
  Extradist1d.push_back(&VertexMuMaxD0SigReco);

  Extradist1d.push_back(&Vertex2DDisplacement);
  Extradist1d.push_back(&Vertex3DDisplacement);

  Extradist1d.push_back(&Vertex2DDisplacementSignificance);
  Extradist1d.push_back(&Vertex3DDisplacementSignificance);
  Extradist1d.push_back(&VertexPairQuality);
  Extradist1d.push_back(&VertexPairfitStatus);
  Extradist1d.push_back(&VertexSignalKFChi2);
  //  Extradist1d.push_back(&VertexSignalAFX);
  //  Extradist1d.push_back(&VertexSignalAFY);
  //  Extradist1d.push_back(&VertexSignalAFZ);
  //  Extradist1d.push_back(&VertexSignalAFChi2);
  //  Extradist1d.push_back(&VertexSignalAFNdf);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexX);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexY);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexZ);
  Extradist1d.push_back(&VertexRefitPVisValid);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexX);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexY);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexZ);

  Extradist1d.push_back(&Muon_segmentCompatibility_mu1);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu2);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu3);

  Extradist1d.push_back(&Muon1_globalDeltaEtaPhi);
  Extradist1d.push_back(&Muon2_globalDeltaEtaPhi);
  Extradist1d.push_back(&Muon3_globalDeltaEtaPhi);

  Extradist1d.push_back(&MuonkinkMin);
  Extradist1d.push_back(&MuonglbkinkMin);
  Extradist1d.push_back(&MuonglbkinkSum);
  Extradist1d.push_back(&MuonkinkSum);


  Extradist1d.push_back(&Isolation_Muon_hadVetoEt03);
  Extradist1d.push_back(&Isolation_Muon_hadEt03);
  Extradist1d.push_back(&Isolation_Muon_emVetoEt03);
  Extradist1d.push_back(&Isolation_Muon_emEt03);
  Extradist1d.push_back(&Isolation_Muon_nJets03);
  Extradist1d.push_back(&Isolation_Muon_nTracks03);
  Extradist1d.push_back(&Isolation_Muon_sumPt03);
  Extradist1d.push_back(&Isolation_Muon_trackerVetoPt03);
  Extradist1d.push_back(&Isolation_Muon_hadVetoEt05);
  Extradist1d.push_back(&Isolation_Muon_hadEt05);
  Extradist1d.push_back(&Isolation_Muon_emVetoEt05);
  Extradist1d.push_back(&Isolation_Muon_emEt05);
  Extradist1d.push_back(&Isolation_Muon_nJets05);
  Extradist1d.push_back(&Isolation_Muon_nTracks05);
  Extradist1d.push_back(&Isolation_Muon_sumPt05);
  Extradist1d.push_back(&Isolation_Muon_trackerVetoPt05);
  Extradist1d.push_back(&Isolation_Muon_sumChargedHadronPt03);
  Extradist1d.push_back(&Isolation_Muon_sumChargedParticlePt03);
  Extradist1d.push_back(&Isolation_Muon_sumNeutralHadronEt03);
  Extradist1d.push_back(&Isolation_Muon_sumPhotonEt03);
  Extradist1d.push_back(&Isolation_Muon_sumPUPt03);
  Extradist1d.push_back(&Isolation_Muon_sumNeutralHadronEtHighThreshold03);
  Extradist1d.push_back(&Isolation_Muon_sumPhotonEtHighThreshold03);
  Extradist1d.push_back(&Isolation_Muon_sumChargedHadronPt04);
  Extradist1d.push_back(&Isolation_Muon_sumChargedParticlePt04);
  Extradist1d.push_back(&Isolation_Muon_sumNeutralHadronEt04);
  Extradist1d.push_back(&Isolation_Muon_sumPhotonEt04);
  Extradist1d.push_back(&Isolation_Muon_sumPUPt04);
  Extradist1d.push_back(&Isolation_Muon_sumNeutralHadronEtHighThreshold04);
  Extradist1d.push_back(&Isolation_Muon_sumPhotonEtHighThreshold04);

}


void  ThreeMu::doEvent(){ 

  
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

    Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
    Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
    Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);

    Muon1Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Eta(),1);
    Muon2Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Eta(),1);
    Muon3Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Eta(),1);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    Mu1TauPTRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt(),w);
    Mu2TauPTRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt(),w);
    Mu3TauPTRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt(),w);

    Mu1TauPRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).P()/TauLV.P(),w);
    Mu2TauPRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).P()/TauLV.P(),w);
    Mu3TauPRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).P()/TauLV.P(),w);
    
    dR12.at(t).Fill(Muon1LV.DeltaR(Muon2LV),1);
    dR23.at(t).Fill(Muon2LV.DeltaR(Muon3LV),1);
    dR31.at(t).Fill(Muon1LV.DeltaR(Muon3LV),1);

    float maxDrMuons = std::max({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
    float minDrMuons = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
 
    maxMudR.at(t).Fill(maxDrMuons,w);
    minMudR. at(t).Fill(minDrMuons,w);

    dR1Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
    dR2Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
    dR3Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);

    float maxDrTauMu =  std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});
    float minDrTauMu =  std::min({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});

    dRMaxMuTau.at(t).Fill(maxDrTauMu);
    dRMinMuTau.at(t).Fill(minDrTauMu);


    TauEta.at(t).Fill(TauLV.Eta(),1);
    TauPt.at(t).Fill(TauLV.Pt(),1);
    TauP.at(t).Fill(TauLV.P(),1);
    TauMass.at(t).Fill(TauLV.M(),1);
    TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    


    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_1).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)==1)  Muon1StandardSelector.at(t).Fill(iMuSelector,1);
    }

    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_2).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)==1)  Muon2StandardSelector.at(t).Fill(iMuSelector,1);
    }

    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_3).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)==1)  Muon3StandardSelector.at(t).Fill(iMuSelector,1);
    }


    Muon_segmentCompatibility_mu1.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_1),w);
    Muon_segmentCompatibility_mu2.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_2),w);
    Muon_segmentCompatibility_mu3.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_3),w);


    Muon1isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_1),1);
    Muon2isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_2),1);
    Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);


    Muon1isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_1),w);
    Muon2isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_2),w);
    Muon3isStand.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_index_3),w);


    Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
    Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
    Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);

    Muon1kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1),1);
    Muon2kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_2),1);
    Muon3kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_3),1);

    Muon1glbkink.at(t).Fill(Ntp->Muon_combinedQuality_glbKink(Muon_index_1),1);
    Muon2glbkink.at(t).Fill(Ntp->Muon_combinedQuality_glbKink(Muon_index_2),1);
    Muon3glbkink.at(t).Fill(Ntp->Muon_combinedQuality_glbKink(Muon_index_3),1);

    float MaxMuonKink =  std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),
	  Ntp->Muon_combinedQuality_trkKink(Muon_index_2),
	  Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});

    float MaxMuonglbKink =  std::max({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),
	  Ntp->Muon_combinedQuality_glbKink(Muon_index_2),
	  Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});

    float MinMuonKink =  std::min({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),
	  Ntp->Muon_combinedQuality_trkKink(Muon_index_2),
	  Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});

    float MinMuonglbKink =  std::min({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),
	  Ntp->Muon_combinedQuality_glbKink(Muon_index_2),
	  Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});


    MuonkinkMin.at(t).Fill(MinMuonKink,w);
    MuonglbkinkMin.at(t).Fill(MinMuonglbKink,w);


    MuonglbkinkSum.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+
			      Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+
			      Ntp->Muon_combinedQuality_trkKink(Muon_index_3),w);

    MuonkinkSum.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+
			   Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+
			   Ntp->Muon_combinedQuality_trkKink(Muon_index_3),w);




    MuonkinkMax.at(t).Fill(MaxMuonKink,1);
    MuonglbkinkMax.at(t).Fill(MaxMuonglbKink,1);

    Muon1InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),1);
    Muon2InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),1);
    Muon3InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),1);

    float maxInOuterTrackerMatch = std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),
	  Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),
	  Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)});

    MuonInOutTrackMatchMax.at(t).Fill(maxInOuterTrackerMatch,1);

    MuPair1_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,0),1);
    MuPair2_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,1),1);
    MuPair3_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,2),1);
    float MinMuPairVertex_chi2 = std::min({Ntp->Vertex_pair_quality(final_idx,0),
	  Ntp->Vertex_pair_quality(final_idx,1),
	  Ntp->Vertex_pair_quality(final_idx,2)});

    MuPair_vertex_chi2_min.at(t).Fill(MinMuPairVertex_chi2,1);


    Muon1_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_1),w);
    Muon2_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_2),w);
    Muon3_globalDeltaEtaPhi.at(t).Fill(Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_3),w);



    Pair1Mass.at(t).Fill((Muon1LV + Muon2LV).M(),1);
    Pair2Mass.at(t).Fill((Muon2LV + Muon3LV).M(),1);
    Pair3Mass.at(t).Fill((Muon1LV + Muon3LV).M(),1);

    Pair1Mass_OS1.at(t).Fill((MuonOS+MuonSS1).M(),1);
    Pair2Mass_OS2.at(t).Fill((MuonOS+MuonSS1).M(),1);
    Pair3Mass_SS.at(t).Fill((MuonSS1+MuonSS2).M(),1);


    TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0),1);
    TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1),1);
    TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2),1);

    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(final_idx),w);
    Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(final_idx),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(final_idx),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(final_idx),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(final_idx),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(final_idx),w);


    //    std::cout<<" min dist  "<< Ntp->Isolation_MinDist(final_idx) << std::endl;
    // -------------------------------- isolation 
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(final_idx),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(final_idx),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(final_idx),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(final_idx),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(final_idx),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(final_idx),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(final_idx),w);

    Isolation_RelativePtMu1.at(t).Fill(Ntp->Isolation_Mu1RelIso(final_idx),w);
    Isolation_RelativePtMu2.at(t).Fill(Ntp->Isolation_Mu2RelIso(final_idx),w);
    Isolation_RelativePtMu3.at(t).Fill(Ntp->Isolation_Mu3RelIso(final_idx),w);
    Isolation_RelativePtMaxMu.at(t).Fill(Ntp->Isolation_MuMaxRelIso(final_idx),w);


    Isolation_Muon_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(Muon_index_1)  ,w);
    Isolation_Muon_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(Muon_index_1)  ,w);
    Isolation_Muon_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(Muon_index_1)  ,w);
    Isolation_Muon_emEt03.at(t).Fill(Ntp->Muon_emEt03(Muon_index_1)  ,w);
    Isolation_Muon_nJets03.at(t).Fill(Ntp->Muon_nJets03(Muon_index_1)  ,w);
    Isolation_Muon_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(Muon_index_1)  ,w);
    Isolation_Muon_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(Muon_index_1)  ,w);
    Isolation_Muon_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(Muon_index_1)  ,w);
    Isolation_Muon_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(Muon_index_1)  ,w);
    Isolation_Muon_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(Muon_index_1)  ,w);
    Isolation_Muon_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(Muon_index_1)  ,w);
    Isolation_Muon_emEt05.at(t).Fill(Ntp->Muon_emEt05(Muon_index_1)  ,w);
    Isolation_Muon_nJets05.at(t).Fill(Ntp->Muon_nJets05(Muon_index_1)  ,w);
    Isolation_Muon_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(Muon_index_1)  ,w);
    Isolation_Muon_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(Muon_index_1)  ,w);
    Isolation_Muon_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(Muon_index_1)  ,w);
    Isolation_Muon_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(Muon_index_1)  ,w);
    Isolation_Muon_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(Muon_index_1)  ,w);
    Isolation_Muon_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(Muon_index_1)  ,w);
    Isolation_Muon_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(Muon_index_1)  ,w);
    Isolation_Muon_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(Muon_index_1)  ,w);
    Isolation_Muon_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_1)  ,w);
    Isolation_Muon_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_1)  ,w);
    Isolation_Muon_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(Muon_index_1)  ,w);
    Isolation_Muon_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(Muon_index_1)  ,w);
    Isolation_Muon_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(Muon_index_1)  ,w);
    Isolation_Muon_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(Muon_index_1)  ,w);
    Isolation_Muon_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(Muon_index_1)  ,w);
    Isolation_Muon_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(Muon_index_1)  ,w);
    Isolation_Muon_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(Muon_index_1)  ,w);


    VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx),w);
    VertexChi2AF.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(final_idx),w);
    VertexDCA12.at(t).Fill(Ntp->Vertex_DCA12(final_idx),w);
    VertexDCA23.at(t).Fill(Ntp->Vertex_DCA23(final_idx),w);
    VertexDCA31.at(t).Fill(Ntp->Vertex_DCA31(final_idx),w);

    float MaxDCA = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
    VertexDCAMax.at(t).Fill(MaxDCA,w);


    VertexSignalKFRefittedMu1P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0).P(),w);
    VertexSignalKFRefittedMu1Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0).Pt(),w);
    VertexSignalKFRefittedMu1Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0).Eta(),w);
    VertexSignalKFRefittedMu1Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0).Phi(),w);
    VertexSignalKFRefittedMu2P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1).P(),w);
    VertexSignalKFRefittedMu2Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1).Pt(),w);
    VertexSignalKFRefittedMu2Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1).Eta(),w);
    VertexSignalKFRefittedMu2Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1).Phi(),w);
    VertexSignalKFRefittedMu3P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2).P(),w);
    VertexSignalKFRefittedMu3Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2).Pt(),w);
    VertexSignalKFRefittedMu3Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2).Eta(),w);
    VertexSignalKFRefittedMu3Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2).Phi(),w);
    VertexMu1D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,0),w);
    VertexMu1D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,0),w);
    VertexMu2D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,1),w);
    VertexMu2D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,1),w);
    VertexMu3D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,2),w);
    VertexMu3D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,2),w);
    Vertex2DDisplacement.at(t).Fill(Ntp->Vertex_2Ddisplacement(final_idx,0),w);
    Vertex3DDisplacement.at(t).Fill(Ntp->Vertex_3Ddisplacement(final_idx,0),w);

    float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0),
	  Ntp->Vertex_d0sig_reco(final_idx,1),
	  Ntp->Vertex_d0sig_reco(final_idx,2)});


    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0),
	  Ntp->Vertex_d0sig_reco(final_idx,1),
	  Ntp->Vertex_d0sig_reco(final_idx,2)});


    VertexMuMinD0SigReco.at(t).Fill(MinD0Significance,1);
    VertexMuMaxD0SigReco.at(t).Fill(MaxD0Significance,1);

    Vertex2DDisplacementSignificance.at(t).Fill(Ntp->Vertex_2Ddisplacement(final_idx,1),w);
    Vertex3DDisplacementSignificance.at(t).Fill(Ntp->Vertex_3Ddisplacement(final_idx,1),w);



    VertexPairQuality.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,0),w);
    VertexPairfitStatus.at(t).Fill(Ntp->Vertex_pairfit_status(final_idx,0),w);
    VertexSignalKFChi2.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx),w);
    /*
      VertexSignalAFX.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).X(),w);
      VertexSignalAFY.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).Y(),w);
      VertexSignalAFZ.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).Z(),w);
      VertexSignalAFChi2.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(0),w);
      VertexSignalAFNdf.at(t).Fill(Ntp->Vertex_signal_AF_Ndf(0),w);
    */

    VertexMatchedPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx).x(),w);
    VertexMatchedPrimaryVertexY.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx).y(),w);
    VertexMatchedPrimaryVertexZ.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(final_idx).z(),w);


    FLSignificance2.at(t).Fill( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx),Ntp->Vertex_PrimaryVertex_Covariance(final_idx),
							       Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx)),w);

    FLSignificance.at(t).Fill(sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx),Ntp->Vertex_PrimaryVertex_Covariance(final_idx),
								   Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx))),w);

    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
    SVPVTauDirAngle.at(t).Fill(SVPV.Angle(TauLV.Vect()),w);

    VertexRefitPVisValid.at(t).Fill(Ntp->Vertex_RefitPVisValid(final_idx),w);
    if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
      VertexMatchedRefitPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx).x(),w);
      VertexMatchedRefitPrimaryVertexY.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx).y(),w);
      VertexMatchedRefitPrimaryVertexZ.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx).z(),w);
    }

    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){

	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
	Muon1PtResolution.at(t).Fill((Muon1LV.Pt() - MCMuon1LV.Pt())/MCMuon1LV.Pt(), 1);
	Muon2PtResolution.at(t).Fill((Muon2LV.Pt() - MCMuon2LV.Pt())/MCMuon2LV.Pt(), 1);
	Muon3PtResolution.at(t).Fill((Muon3LV.Pt() - MCMuon3LV.Pt())/MCMuon3LV.Pt(), 1);

	Muon1EtaResolution.at(t).Fill((Muon1LV.Eta() - MCMuon1LV.Eta())/MCMuon1LV.Eta(), 1);
	Muon2EtaResolution.at(t).Fill((Muon2LV.Eta() - MCMuon2LV.Eta())/MCMuon2LV.Eta(), 1);
	Muon3EtaResolution.at(t).Fill((Muon3LV.Eta() - MCMuon3LV.Eta())/MCMuon3LV.Eta(), 1);

	TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }
    
  }
}


void  ThreeMu::Finish(){

  if(mode == RECONSTRUCT){
    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/3;
      ScaleAllHistOfType(i,scale);
    }
  }
  Selection::Finish();
}





