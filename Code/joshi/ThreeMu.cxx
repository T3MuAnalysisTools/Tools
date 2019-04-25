#include "ThreeMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

ThreeMu::ThreeMu(TString Name_, TString id_):
  Selection(Name_,id_)
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
    if(i==HLTOk)               cut.at(HLTOk)=1;
    if(i==ThreeMuCandidate)    cut.at(ThreeMuCandidate)=1;
  }
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    
	 if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==ThreeMuCandidate){
      title.at(i)="ThreeMuCandidate is found ";
      hlabel="is Three Mu Candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention

  Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",75,0,25,"  #mu_{1} p_{T}, GeV","Events");
  Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",75,0,20,"  #mu_{2} p_{T}, GeV","Events");
  Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",75,0,15,"  #mu_{3} p_{T}, GeV","Events");

  Muon1isGlob =HConfig.GetTH1D(Name+"_Muon1isGlob","Muon1isGlob",2,-0.5,1.5,"  #mu_{1} is global muon","Events");
  Muon2isGlob =HConfig.GetTH1D(Name+"_Muon2isGlob","Muon2isGlob",2,-0.5,1.5,"  #mu_{2} is global muon","Events");
  Muon3isGlob =HConfig.GetTH1D(Name+"_Muon3isGlob","Muon3isGlob",2,-0.5,1.5,"  #mu_{3} is global muon","Events");

  Muon1isTrack =HConfig.GetTH1D(Name+"_Muon1isTrack","Muon1isTrack",2,-0.5,1.5,"  #mu_{1} is tracker muon","Events");
  Muon2isTrack =HConfig.GetTH1D(Name+"_Muon2isTrack","Muon2isTrack",2,-0.5,1.5,"  #mu_{2} is tracker muon","Events");
  Muon3isTrack =HConfig.GetTH1D(Name+"_Muon3isTrack","Muon3isTrack",2,-0.5,1.5,"  #mu_{3} is tracker muon","Events");

  Muon1kink =HConfig.GetTH1D(Name+"_Muon1kink","Muon1kink",50,0.,50,"  #mu_{1} kink","Events");
  Muon2kink =HConfig.GetTH1D(Name+"_Muon2kink","Muon2kink",50,0.,50,"  #mu_{2} kink","Events");
  Muon3kink =HConfig.GetTH1D(Name+"_Muon3kink","Muon3kink",50,0.,50,"  #mu_{3} kink","Events");

  Muon1InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon1InOutTrackMatch","Muon1InOutTrackMatch",50,0.,10,"  #mu_{1} inner and outer tracker match","Events");
  Muon2InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon2InOutTrackMatch","Muon2InOutTrackMatch",50,0.,10,"  #mu_{2} inner and outer tracker match","Events");
  Muon3InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon3InOutTrackMatch","Muon3InOutTrackMatch",50,0.,10,"  #mu_{3} inner and outer tracker match","Events");

  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",25,-2.6,2.6,"#mu_{1}  rapidity","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",25,-2.6,2.6,"#mu_{2}  rapidity","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",25,-2.6,2.6,"#mu_{3}  rapidity","Events");

  Muon1StandardSelector=HConfig.GetTH1D(Name+"_Muon1StandardSelector","Muon1StandardSelector",23,-0.5,22.5,"#mu_{1} standard selector; bin 0 - no ID","Events");
  Muon2StandardSelector=HConfig.GetTH1D(Name+"_Muon2StandardSelector","Muon2StandardSelector",23,-0.5,22.5,"#mu_{2} standard selector; bin 0 - no ID","Events");
  Muon3StandardSelector=HConfig.GetTH1D(Name+"_Muon3StandardSelector","Muon3StandardSelector",23,-0.5,22.5,"#mu_{3} standard selector; bin 0 - no ID","Events");

  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#tau rapidity","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"  #tau p_{T}, GeV","Events");
  TauP =HConfig.GetTH1D(Name+"_TauP","TauP",30,0,50,"  #tau |p|, GeV","Events");
  TauMass =HConfig.GetTH1D(Name+"_TauMass","#tau lepton mass",50,1.5,1.9,"  M_{#tau} , GeV","Events");
  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassRefit =HConfig.GetTH1D(Name+"_TauMassRefit","Refit #tau lepton mass",50,1.5,1.9,"KF refit  M_{#tau} , GeV","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  VertexDCA12=HConfig.GetTH1D(Name+"_VertexDCA12","VertexDCA12",100,0,0.5,"dca (#mu_{1}#mu_{2})","Events");
  VertexDCA23=HConfig.GetTH1D(Name+"_VertexDCA23","VertexDCA23",100,0,0.5,"dca (#mu_{1}trk)","Events");
  VertexDCA31=HConfig.GetTH1D(Name+"_VertexDCA31","VertexDCA31",100,0,0.5,"dca (#mu_{2}trk)","Events");
  VertexDCAMax=HConfig.GetTH1D(Name+"_VertexDCAMax","VertexDCAMax",100,0,0.5,"max dca","Events");
  
  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",100,0,10,"KF vertex #chi^{2}","Events");
  VertexChi2AF=HConfig.GetTH1D(Name+"_VertexChi2AF","VertexChi2AF",100,0,10,"AF vertex #chi^{2}","Events");

  VertexSignalKFRefittedMu1P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1P","VertexSignalKFRefittedMu1P",200,0,20,"KF refitted #mu_{1} track p (GeV)","Events");
  VertexSignalKFRefittedMu1Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Pt","VertexSignalKFRefittedMu1P",200,0,20,"KF refitted #mu_{1} track p_{T} (GeV)","Events");
  VertexSignalKFRefittedMu1Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Eta","VertexSignalKFRefittedMu1Eta",250,-2.5,2.5,"KF refitted #mu_{1} track #eta","Events");
  VertexSignalKFRefittedMu1Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Phi","VertexSignalKFRefittedMu1Phi",320,-3.2,3.2,"KF refitted #mu_{1} track #phi","Events");

  VertexSignalKFRefittedMu2P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2P","VertexSignalKFRefittedMu2P",200,0,20,"KF refitted #mu_{2} track p (GeV)","Events");
  VertexSignalKFRefittedMu2Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Pt","VertexSignalKFRefittedMu2Pt",200,0,20,"KF refitted #mu_{2} track p_{T} (GeV)","Events");
  VertexSignalKFRefittedMu2Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Eta","VertexSignalKFRefittedMu2Eta",250,-2.5,2.5,"KF refitted #mu_{2} track #eta","Events");
  VertexSignalKFRefittedMu2Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Phi","VertexSignalKFRefittedMu2Phi",320,-3.2,3.2,"KF refitted #mu_{2} track #phi","Events");

  VertexSignalKFRefittedTrackP=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedTrackP","VertexSignalKFRefittedTrackP",200,0,20,"KF refitted Track p (GeV)","Events");
  VertexSignalKFRefittedTrackPt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedTrackPt","VertexSignalKFRefittedTrackPt",200,0,20,"KF refitted Track p_{T} (GeV)","Events");
  VertexSignalKFRefittedTrackEta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedTrackEta","VertexSignalKFRefittedTrackEta",250,-2.5,2.5,"KF refitted Track #eta","Events");
  VertexSignalKFRefittedTrackPhi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedTrackPhi","VertexSignalKFRefittedTrackPhi",320,-3.2,3.2,"KF refitted Track #phi","Events");
  
  VertexMu1D0Reco=HConfig.GetTH1D(Name+"_VertexMu1D0Reco","VertexMu1D0Reco",100,0,5,"#mu_{1} d0 vertex","Events");
  VertexMu1D0SigReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigReco","VertexMu1D0SigReco",100,0,5,"#mu_{1} d0 vertex significance","Events");
  VertexMu2D0Reco=HConfig.GetTH1D(Name+"_VertexMu2D0Reco","VertexMu2D0Reco",100,0,5,"#mu_{2} d0 vertex","Events");
  VertexMu2D0SigReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigReco","VertexMu2D0SigReco",100,0,5,"#mu_{2} d0 vertex significance","Events");
  VertexTrackD0Reco=HConfig.GetTH1D(Name+"_VertexTrackD0Reco","VertexTrackD0Reco",100,0,5,"track d0 vertex","Events");
  VertexTrackD0SigReco=HConfig.GetTH1D(Name+"_VertexTrackD0SigReco","VertexTrackD0SigReco",100,0,5,"track d0 vertex significance","Events");
  Vertex2DDisplacement=HConfig.GetTH1D(Name+"_Vertex2DDisplacement","Vertex2DDisplacement",10,0,1,"vertex 2d displacement","Events");
  Vertex3DDisplacement=HConfig.GetTH1D(Name+"_Vertex3DDisplacement","Vertex3DDisplacement",10,0,1,"vertex 3d displacement","Events");
  VertexPairQuality=HConfig.GetTH1D(Name+"_VertexPairQuality","VertexPairQuality",10,0,10,"vertex pair quality","Events");
  VertexPairfitStatus=HConfig.GetTH1D(Name+"_VertexPairfitStatus","VertexPairfitStatus",100,0,10,"vertex pair fit status","Events");
  VertexSignalKFChi2=HConfig.GetTH1D(Name+"_VertexSignalKFChi2","VertexSignalKFChi2",10,0,10,"vertex KF #chi^{2}","Events");
  VertexSignalAFX=HConfig.GetTH1D(Name+"_VertexSignalAFX","VertexSignalAFX",100,0,10,"AF vertex x (cm)","Events");
  VertexSignalAFY=HConfig.GetTH1D(Name+"_VertexSignalAFY","VertexSignalAFY",100,0,10,"AF vertex y (cm)","Events");
  VertexSignalAFZ=HConfig.GetTH1D(Name+"_VertexSignalAFZ","VertexSignalAFZ",100,0,10,"AF vertex z (cm)","Events");
  VertexSignalAFChi2=HConfig.GetTH1D(Name+"_VertexSignalChi2","VertexSignalChi2",100,0,10,"AF vertex #chi^{2}","Events");
  VertexSignalAFNdf=HConfig.GetTH1D(Name+"_VertexSignalAFNdf","VertexSignalAFNdf",10,0,10,"AF vertex ndf","Events");
  VertexMatchedPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexX","VertexMatchedPrimaryVertexX",100,0,0.5,"vertex matched pv x (cm)","Events");
  VertexMatchedPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexY","VertexMatchedPrimaryVertexY",100,0,0.5,"vertex matched pv y (cm)","Events");
  VertexMatchedPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexZ","VertexMatchedPrimaryVertexZ",100,0,0.5,"vertex matched pv z (cm)","Events");
  VertexRefitPVisValid=HConfig.GetTH1D(Name+"_VertexRefitPVisValid","VertexRefitPVisValid",2,0,2,"vertex refit pv is valid","Events");
  VertexMatchedRefitPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexX","VertexMatchedRefitPrimaryVertexX",100,0,0.5,"vertex matched refit pv x (cm)","Events");
  VertexMatchedRefitPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexY","VertexMatchedRefitPrimaryVertexY",100,0,0.5,"vertex matched refit pv y (cm)","Events");
  VertexMatchedRefitPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexZ","VertexMatchedRefitPrimaryVertexZ",100,0,0.5,"vertex matched refit pv z (cm)","Events");
  
  Muon1PtResolution=HConfig.GetTH1D(Name+"_Muon1PtResolution","Muon1PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{1})  (reco - mc)/mc ","Events");
  Muon2PtResolution=HConfig.GetTH1D(Name+"_Muon2PtResolution","Muon2PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{2})  (reco - mc)/mc ","Events");
  Muon3PtResolution=HConfig.GetTH1D(Name+"_Muon3PtResolution","Muon3PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{3})  (reco - mc)/mc  ","Events");

  Muon1EtaResolution=HConfig.GetTH1D(Name+"_Muon1EtaResolution","Muon1EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc  ","Events");
  Muon2EtaResolution=HConfig.GetTH1D(Name+"_Muon2EtaResolution","Muon2EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
  Muon3EtaResolution=HConfig.GetTH1D(Name+"_Muon3EtaResolution","Muon3EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");


  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");

  MuPair1_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair1_vertex_chi2","MuPair1_vertex_chi2",50,0,5,"KF  #chi^{2} of first #mu pair","Events");
  MuPair2_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair2_vertex_chi2","MuPair2_vertex_chi2",50,0,5,"KF  #chi^{2} of second #mu pair","Events");
  MuPair3_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair3_vertex_chi2","MuPair3_vertex_chi2",50,0,5,"KF  #chi^{2} of third #mu pair","Events");

  Pair1Mass =HConfig.GetTH1D(Name+"_Pair1Mass","Pair1Mass",50,0,2," mass of #mu pair (12), GeV","Events");
  Pair2Mass =HConfig.GetTH1D(Name+"_Pair2Mass","Pair2Mass",50,0,2," mass of #mu pair (23), GeV","Events");
  Pair3Mass =HConfig.GetTH1D(Name+"_Pair3Mass","Pair3Mass",50,0,2," mass of #mu pair (31), GeV","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match dR 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match dR 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match dR 3","Events");

  dR12 =HConfig.GetTH1D(Name+"_dR12","dR12",50,0,1,"dR(#mu_{1}#mu_{2})","Events");
  dR23 =HConfig.GetTH1D(Name+"_dR23","dR23",50,0,1,"dR(#mu_{2}#mu_{3})","Events");
  dR31 = HConfig.GetTH1D(Name+"_dR31","dR31",50,0,1,"dR(#mu_{3}#mu_{1})","Events");
  dR1Tau = HConfig.GetTH1D(Name+"_dR1Tau","dR1Tau",50,0,1,"dR(#mu_{1}#tau)","Events");
  dR2Tau = HConfig.GetTH1D(Name+"_dR2Tau","dR2Tau",50,0,1,"dR(#mu_{2}#tau)","Events");
  dR3Tau = HConfig.GetTH1D(Name+"_dR3Tau","dR3Tau",50,0,1,"dR(#mu_{3}#tau)","Events");

  Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","Isolation_NTracks",20,0,10,"N tracks","Events");
  Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","Isolation_RelPt",20,0,10,"relative p_{T}","Events");
  Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","Isolation_MinDist",10,0,1,"Iso MinDist","Events");
  Isolation05_RelPt=HConfig.GetTH1D(Name+"_Isolation05_RelPt","Isolation05_RelPt",10,0,3,"relative  rel p_{T} in 0.5 cone","Events");
  Isolation05_NTracks=HConfig.GetTH1D(Name+"_Isolation05_NTracks","Isolation05_NTracks",20,0,10,"N tracks in 0.5 cone","Events");
  Isolation05_MinDist=HConfig.GetTH1D(Name+"_Isolation05_MinDist","Isolation05_MinDist",10,0,1,"Iso05 MinDist","Events");
  Isolation_Ntrk1=HConfig.GetTH1D(Name+"_Isolation_Ntrk1","Isolation_Ntrk1",10,0,10,"Iso ntrk 1","Events");
  Isolation_Ntrk2=HConfig.GetTH1D(Name+"_Isolation_Ntrk2","Isolation_Ntrk2",10,0,10,"Iso ntrk 2","Events");
  Isolation_Ntrk3=HConfig.GetTH1D(Name+"_Isolation_Ntrk3","Isolation_Ntrk3",10,0,10,"Iso ntrk 3","Events");
  Isolation_Ntrk0p1=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p1","Isolation_Ntrk0p1",10,0,10,"Iso ntrk0p1","Events");
  Isolation_Ntrk0p2=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p2","Isolation_Ntrk0p2",10,0,10,"Iso ntrk0p2","Events");
  Isolation_Ntrk0p5=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p5","Isolation_Ntrk0p5",10,0,10,"Iso ntrk0p5","Events");
  Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","Isolation_maxdxy",40,0,20,"Iso max(dxy)","Events");

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

  Extradist1d.push_back(&Muon1kink);
  Extradist1d.push_back(&Muon2kink);
  Extradist1d.push_back(&Muon3kink);

  Extradist1d.push_back(&Muon1InOutTrackMatch);
  Extradist1d.push_back(&Muon2InOutTrackMatch);
  Extradist1d.push_back(&Muon3InOutTrackMatch);

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

  Extradist1d.push_back(&Pair1Mass);
  Extradist1d.push_back(&Pair2Mass);
  Extradist1d.push_back(&Pair3Mass);

  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);

  Extradist1d.push_back(&dR12);
  Extradist1d.push_back(&dR23);
  Extradist1d.push_back(&dR31);
  Extradist1d.push_back(&dR1Tau);
  Extradist1d.push_back(&dR2Tau);
  Extradist1d.push_back(&dR3Tau);

  Extradist1d.push_back(&Isolation_NTracks);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_MinDist);
  Extradist1d.push_back(&Isolation05_RelPt);
  Extradist1d.push_back(&Isolation05_NTracks);
  Extradist1d.push_back(&Isolation05_MinDist);
  Extradist1d.push_back(&Isolation_Ntrk1);
  Extradist1d.push_back(&Isolation_Ntrk2);
  Extradist1d.push_back(&Isolation_Ntrk3);
  Extradist1d.push_back(&Isolation_Ntrk0p1);
  Extradist1d.push_back(&Isolation_Ntrk0p2);
  Extradist1d.push_back(&Isolation_Ntrk0p5);
  Extradist1d.push_back(&Isolation_maxdxy);

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
  Extradist1d.push_back(&VertexSignalKFRefittedTrackP);
  Extradist1d.push_back(&VertexSignalKFRefittedTrackPt);
  Extradist1d.push_back(&VertexSignalKFRefittedTrackEta);
  Extradist1d.push_back(&VertexSignalKFRefittedTrackPhi);

  Extradist1d.push_back(&VertexMu1D0Reco);
  Extradist1d.push_back(&VertexMu1D0SigReco);
  Extradist1d.push_back(&VertexMu2D0Reco);
  Extradist1d.push_back(&VertexMu2D0SigReco);
  Extradist1d.push_back(&VertexTrackD0Reco);
  Extradist1d.push_back(&VertexTrackD0SigReco);
  Extradist1d.push_back(&Vertex2DDisplacement);
  Extradist1d.push_back(&Vertex3DDisplacement);
  Extradist1d.push_back(&VertexPairQuality);
  Extradist1d.push_back(&VertexPairfitStatus);
  Extradist1d.push_back(&VertexSignalKFChi2);
  Extradist1d.push_back(&VertexSignalAFX);
  Extradist1d.push_back(&VertexSignalAFY);
  Extradist1d.push_back(&VertexSignalAFZ);
  Extradist1d.push_back(&VertexSignalAFChi2);
  Extradist1d.push_back(&VertexSignalAFNdf);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexX);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexY);
  Extradist1d.push_back(&VertexMatchedPrimaryVertexZ);
  Extradist1d.push_back(&VertexRefitPVisValid);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexX);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexY);
  Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexZ);

}


void  ThreeMu::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  //  std::cout<<"id  "<<id<< std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  value.at(HLTOk) = 0;
  value.at(ThreeMuCandidate) = 0;
  if(Ntp->NThreeMuons()>0)  value.at(ThreeMuCandidate)=1;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  pass.at(ThreeMuCandidate)= (value.at(ThreeMuCandidate)==cut.at(ThreeMuCandidate)); 


  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
    
    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(0,0)+Ntp->Vertex_signal_KF_refittedTracksP4(0,1)+Ntp->Vertex_signal_KF_refittedTracksP4(0,2);

    Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
    Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
    Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);

    Muon1Eta.at(t).Fill( Ntp->Muon_P4(Muon_index_1).Eta(),1);
    Muon2Eta.at(t).Fill(  Ntp->Muon_P4(Muon_index_2).Eta(),1);
    Muon3Eta.at(t).Fill(  Ntp->Muon_P4(Muon_index_3).Eta(),1);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    
    dR12.at(t).Fill(Muon1LV.DeltaR(Muon2LV),1);
    dR23.at(t).Fill(Muon2LV.DeltaR(Muon3LV),1);
    dR31.at(t).Fill(Muon1LV.DeltaR(Muon3LV),1);
    dR1Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
    dR2Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
    dR3Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);

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

    Muon1isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_1),1);
    Muon2isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_2),1);
    Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);

    Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
    Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
    Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);

    Muon1kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1),1);
    Muon2kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_2),1);
    Muon3kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_3),1);

    Muon1InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),1);
    Muon2InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),1);
    Muon3InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),1);

    MuPair1_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,0),1);
    MuPair2_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,1),1);
    MuPair3_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,2),1);

    Pair1Mass.at(t).Fill((Muon1LV + Muon2LV).M(),1);
    Pair2Mass.at(t).Fill((Muon2LV + Muon3LV).M(),1);
    Pair3Mass.at(t).Fill((Muon1LV + Muon3LV).M(),1);

    TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(0).at(0),1);
    TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(0).at(1),1);
    TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(0).at(2),1);
    
	 Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(0),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(0),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(0),w);
    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(0),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(0),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(0),w);
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(0),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(0),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(0),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(0),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(0),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(0),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(0),w); 
    Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(0),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(0),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(0),w);
    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(0),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(0),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(0),w);
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(0),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(0),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(0),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(0),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(0),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(0),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(0),w);
    
	 VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(0),w);
    VertexChi2AF.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(0),w);
    VertexDCA12.at(t).Fill(Ntp->Vertex_DCA12(0),w);
    VertexDCA23.at(t).Fill(Ntp->Vertex_DCA23(0),w);
    VertexDCA31.at(t).Fill(Ntp->Vertex_DCA31(0),w);
    VertexDCAMax.at(t).Fill(Ntp->Vertex_DCAMax(0),w);
    VertexSignalKFRefittedMu1P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,0).P(),w);
    VertexSignalKFRefittedMu1Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,0).Pt(),w);
    VertexSignalKFRefittedMu1Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,0).Eta(),w);
    VertexSignalKFRefittedMu1Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,0).Phi(),w);
	 VertexSignalKFRefittedMu2P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,1).P(),w);
    VertexSignalKFRefittedMu2Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,1).Pt(),w);
    VertexSignalKFRefittedMu2Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,1).Eta(),w);
    VertexSignalKFRefittedMu2Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,1).Phi(),w);
	 VertexSignalKFRefittedTrackP.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,2).P(),w);
    VertexSignalKFRefittedTrackPt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,2).Pt(),w);
    VertexSignalKFRefittedTrackEta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,2).Eta(),w);
    VertexSignalKFRefittedTrackPhi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(0,2).Phi(),w);
    VertexMu1D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(0,0),w);
    VertexMu1D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(0,0),w);
    VertexMu2D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(0,1),w);
    VertexMu2D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(0,1),w);
    VertexTrackD0Reco.at(t).Fill(Ntp->Vertex_d0_reco(0,2),w);
    VertexTrackD0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(0,2),w);
    Vertex2DDisplacement.at(t).Fill(Ntp->Vertex_2Ddisplacement(0,0),w);
    Vertex3DDisplacement.at(t).Fill(Ntp->Vertex_3Ddisplacement(0,0),w);
    VertexPairQuality.at(t).Fill(Ntp->Vertex_pair_quality(0,0),w);
    VertexPairfitStatus.at(t).Fill(Ntp->Vertex_pairfit_status(0,0),w);
    VertexSignalKFChi2.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(0),w);
    /*
	 VertexSignalAFX.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).X(),w);
    VertexSignalAFY.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).Y(),w);
    VertexSignalAFZ.at(t).Fill(Ntp->Vertex_signal_AF_pos(0).Z(),w);
    VertexSignalAFChi2.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(0),w);
    VertexSignalAFNdf.at(t).Fill(Ntp->Vertex_signal_AF_Ndf(0),w);
	 */
    VertexMatchedPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(0).x(),w);
    VertexMatchedPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(0).y(),w);
    VertexMatchedPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedPrimaryVertex(0).z(),w);
	 VertexRefitPVisValid.at(t).Fill(Ntp->Vertex_RefitPVisValid(0),w);
	 if (Ntp->Vertex_RefitPVisValid(0)==1){
	 VertexMatchedRefitPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(0).x(),w);
	 VertexMatchedRefitPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(0).y(),w);
	 VertexMatchedRefitPrimaryVertexX.at(t).Fill(Ntp->Vertex_MatchedRefitPrimaryVertex(0).z(),w);
	 }

    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){

      TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0)));
      TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1)));
      TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2)));
cout<<"error 16"<<endl;
      TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
      Muon1PtResolution.at(t).Fill((Muon1LV.Pt() - MCMuon1LV.Pt())/MCMuon1LV.Pt(), 1);
      Muon2PtResolution.at(t).Fill((Muon2LV.Pt() - MCMuon2LV.Pt())/MCMuon2LV.Pt(), 1);
      Muon3PtResolution.at(t).Fill((Muon3LV.Pt() - MCMuon3LV.Pt())/MCMuon3LV.Pt(), 1);

      Muon1EtaResolution.at(t).Fill((Muon1LV.Eta() - MCMuon1LV.Eta())/MCMuon1LV.Eta(), 1);
      Muon2EtaResolution.at(t).Fill((Muon2LV.Eta() - MCMuon2LV.Eta())/MCMuon2LV.Eta(), 1);
      Muon3EtaResolution.at(t).Fill((Muon3LV.Eta() - MCMuon3LV.Eta())/MCMuon3LV.Eta(), 1);

      TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
cout<<"error 17"<<endl;
      TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

      Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
      Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
      Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }

  }
}



void  ThreeMu::Finish(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 

  if(mode == RECONSTRUCT){
    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/3;
      ScaleAllHistOfType(i,scale);
    }
  }


  Selection::Finish();

}





