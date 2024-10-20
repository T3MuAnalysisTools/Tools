#include "MCBackgroundStudy.h"
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

MCBackgroundStudy::MCBackgroundStudy(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

MCBackgroundStudy::~MCBackgroundStudy(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  MCBackgroundStudy::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto)            cut.at(PhiVeto)=0; // defined below
    if(i==OmegaVeto)          cut.at(OmegaVeto)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=0.03;
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
      title.at(i)="All mu pass ID";
      hlabel="gl,gl,gl";
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
    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="Sum of dR_{reco-trigger}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
    }
    else if(i==TauMassCut){
      title.at(i)="$\\tau$ mass (sideband in data)";
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



  Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
  Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
  Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");


  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#eta(#mu_{1})","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#eta(#mu_{2})","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#eta(#mu_{3})","Events");

  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
  TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");



  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");


  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");


  IDOriginOfOSMuon =HConfig.GetTH1D(Name+"_IDOriginOfOSMuon","IDOriginOfOSMuon",400,200,600,"PDGID of OS muon origin","Events");
  
  MassEta =HConfig.GetTH1D(Name+"_MassEta","MassEta",50,0,2,"Mass of eta meson","Events");
  MassEtaPrime =HConfig.GetTH1D(Name+"_MassEtaPrime","MassEtaPrime",50,0,2,"Mass of #eta' meson","Events");
  MassPhi =HConfig.GetTH1D(Name+"_MassPhi","MassPhi",50,0,2,"Mass of #phi meson","Events");
  MassOmega =HConfig.GetTH1D(Name+"_MassOmega","MassOmega",50,0,2,"Mass of #omega meson","Events");
  
  MassEtaReco =HConfig.GetTH1D(Name+"_MassEtaReco","MassEtaReco",50,0,2,"Mass of reco #eta meson","Events");
  MassEtaPrimeReco =HConfig.GetTH1D(Name+"_MassEtaPrimeReco","MassEtaPrimeReco",50,0,2,"Mass of reco #eta' meson","Events");
  MassPhiReco =HConfig.GetTH1D(Name+"_MassPhiReco","MassPhiReco",50,0,2,"Mass of reco #phi meson","Events");
  
  MassEtaMixed =HConfig.GetTH1D(Name+"_MassEtaMixed","MassEtaMixed",50,0,2,"Mass of mixed-reco #eta meson","Events");
  MassEtaPrimeMixed =HConfig.GetTH1D(Name+"_MassEtaPrimeMixed","MassEtaPrimeMixed",50,0,2,"Mass of mixed-reco #eta' meson","Events");
  MassPhiMixed =HConfig.GetTH1D(Name+"_MassPhiMixed","MassPhiMixed",50,0,2,"Mass of mixed-reco #phi meson","Events");
  
  PhotonDRToTruth_Eta=HConfig.GetTH1D(Name+"_PhotonDRToTruth_Eta","PhotonDRToTruth_Eta",40,0,1,"reco - mc #gamma #Delta R","Events");
  PhotonDRToTruth_EtaPrime=HConfig.GetTH1D(Name+"_PhotonDRToTruth_EtaPrime","PhotonDRToTruth_EtaPrime",40,0,1,"reco - mc #gamma #Delta R","Events");
  PhotonDRToTruth_Phi=HConfig.GetTH1D(Name+"_PhotonDRToTruth_Phi","PhotonDRToTruth_Phi",40,0,1,"reco - mc #gamma #Delta R","Events");
  
  Photon_hasPixelSeed_Eta=HConfig.GetTH1D(Name+"_Photon_hasPixelSeed_Eta","Photon_hasPixelSeed_Eta",2,-0.5,1.5,"Whether photon has pixel seed","Events");
  Photon_hasPixelSeed_EtaPrime=HConfig.GetTH1D(Name+"_Photon_hasPixelSeed_EtaPrime","Photon_hasPixelSeed_EtaPrime",2,-0.5,1.5,"Whether photon has pixel seed","Events");
  Photon_hasPixelSeed_Phi=HConfig.GetTH1D(Name+"_Photon_hasPixelSeed_Phi","Photon_hasPixelSeed_Phi",2,-0.5,1.5,"Whether photon has pixel seed","Events");
  
  Photon_hasPixelSeed=HConfig.GetTH1D(Name+"_Photon_hasPixelSeed","Photon_hasPixelSeed",2,-0.5,1.5,"Whether photon has pixel seed","Events");
  Photon_hasConversionTracks=HConfig.GetTH1D(Name+"_Photon_hasConversionTracks","Photon_hasConversionTracks",2,-0.5,1.5,"Whether photon has conversion tracks","Events");
  Photon_isPF=HConfig.GetTH1D(Name+"_Photon_isPF","Photon_isPF",2,-0.5,1.5,"Whether photon is particle flow","Events");
  
  DeltaRPhotontoTau=HConfig.GetTH1D(Name+"_DeltaRPhotontoTau","DeltaRPhotontoTau",4,0.5,4.5,"DeltaR to Tau of all reco #gamma 0.1,0.5,1,2","Events");
  PhotonRecoSuccess=HConfig.GetTH1D(Name+"_PhotonRecoSuccess","PhotonRecoSuccess",2,-2,2,"Whether atleast 1 photon is found","Events");
  DeltaEnergyPhoton=HConfig.GetTH1D(Name+"_DeltaEnergyPhoton","DeltaEnergyPhoton",20,0,20,"reco - mc #gamma #Delta E","Events");
  
  
  PassedCount =HConfig.GetTH1D(Name+"_PassedCount","PassedCount",4,-0.5,3.5,"How many indices changed value from -1","Events");
  ChildCount =HConfig.GetTH1D(Name+"_ChildCount","ChildCount",10,-0.5,9.5,"Number of children of particle with pdgid 221","Events");
  ChildCountPrime =HConfig.GetTH1D(Name+"_ChildCountPrime","ChildCountPrime",10,-0.5,9.5,"Number of children of particle with pdgid 331","Events"); 
   
  PhotonSpectrum =HConfig.GetTH2D(Name+"_PhotonSpectrum","Energy/Rapidity of Photon from #eta",40,0,20,40,-3.0,3.0,"Energy, Gev","Rapidity");
  PhotonSpectrumPrime =HConfig.GetTH2D(Name+"_PhotonSpectrumPrime","Energy/Rapidity of Photon from #eta'",40,0,20,40,-3.0,3.0,"Energy, Gev","Rapidity");
  PhotonSpectrumPhi =HConfig.GetTH2D(Name+"_PhotonSpectrumPhi","Energy/Rapidity of Photon from #phi'",40,0,20,40,-3.0,3.0,"Energy, Gev","Rapidity");

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  MCBackgroundStudy::Store_ExtraDist(){ 


  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauP);


  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);


  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);
  
  Extradist1d.push_back(&NSignalCandidates);
  Extradist1d.push_back(&IDOriginOfOSMuon);
  
  Extradist1d.push_back(&MassEta);
  Extradist1d.push_back(&MassEtaPrime);
  Extradist1d.push_back(&MassPhi);
  Extradist1d.push_back(&MassOmega);
  
  Extradist1d.push_back(&MassEtaMixed);
  Extradist1d.push_back(&MassEtaPrimeMixed);
  Extradist1d.push_back(&MassPhiMixed);
  
  Extradist1d.push_back(&MassEtaReco);
  Extradist1d.push_back(&MassEtaPrimeReco);
  Extradist1d.push_back(&MassPhiReco);
  
  Extradist1d.push_back(&PhotonDRToTruth_Eta);
  Extradist1d.push_back(&PhotonDRToTruth_EtaPrime);
  Extradist1d.push_back(&PhotonDRToTruth_Phi);
  
  Extradist1d.push_back(&PassedCount);
  Extradist1d.push_back(&ChildCount);
  Extradist1d.push_back(&ChildCountPrime);
  
  Extradist2d.push_back(&PhotonSpectrum);
  Extradist2d.push_back(&PhotonSpectrumPrime);
  Extradist2d.push_back(&PhotonSpectrumPhi);
  
  Extradist1d.push_back(&DeltaRPhotontoTau);
  Extradist1d.push_back(&PhotonRecoSuccess);
  Extradist1d.push_back(&DeltaEnergyPhoton);
  
  Extradist1d.push_back(&Photon_hasPixelSeed_Eta);
  Extradist1d.push_back(&Photon_hasPixelSeed_EtaPrime);
  Extradist1d.push_back(&Photon_hasPixelSeed_Phi);
  
  Extradist1d.push_back(&Photon_hasPixelSeed);
  Extradist1d.push_back(&Photon_hasConversionTracks);
  Extradist1d.push_back(&Photon_isPF);


}


void  MCBackgroundStudy::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  //  std::cout<<" id   "<< id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  bool HLTOk(false);
  bool L1Ok(false);
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);

    if(id==1){
    if( HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") && Ntp->HLTDecision(iTrigger)==1) HLTOk = true;
    }

    //    if(id==1 && Ntp->WhichEra(2017).Contains("RunF")){
    //      if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
    //    }	 

    if(id!=1){
    if( HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") && Ntp->HLTDecision(iTrigger)==1) HLTOk = true;
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
  value.at(TriggerOk)=(HLTOk);// and L1Ok);
  pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));



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
			//    			Ntp->Muon_isGlobalMuon(mu3_pt_idx) or Ntp->Muon_isTrackerMuon(mu3_pt_idx));
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

    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx)){
      value.at(TriggerMatch)+=i; 
    }
    
    //    value.at(TauMassCut) = TauLV.M();
    value.at(TauMassCut) = TauRefittedLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  //  pass.at(VertChi2) = (value.at(VertChi2) <= cut.at(VertChi2));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = true;//(value.at(TriggerMatch)  <  cut.at(TriggerMatch));
  pass.at(PhiVeto) = true;//(fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 8*PDG_Var::Phi_width());
  pass.at(OmegaVeto) = true;//(fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 3*PDG_Var::Omega_width());

  if(id!=1) pass.at(TauMassCut) = true;
  else  pass.at(TauMassCut) =( (value.at(TauMassCut) > tauMinSideBand_)  ||   (value.at(TauMassCut) < tauMaxSideBand_ ));


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



    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    /*    // This is to print out selected events content
    std::cout<<"------------------------------- "<< std::endl;
    std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
    std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
    std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
    Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
    Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
    Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
    Ntp->printMCDecayChainOfEvent(true, true, true, true);
    std::cout<< "\n\n\n\n\n\n";
    */


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

    // ----  MC 
    if(id == 122 || id == 128 || id == 132 || id == 142){ // id == 122 is Ds-> eta(mu mu gamma) mu nu   + eta' (mu mu gamma) mu nu
    
            double Dr1(99.), Dr2(99.), Dr3(99.);
      int mc_index1,mc_index2, mc_index3;
      for(int igen_particle=0; igen_particle < Ntp->NMCParticles(); igen_particle++){
	
	      if(Ntp->MCParticle_p4(igen_particle).DeltaR(MuonOS) < Dr1){
	        Dr1 = Ntp->MCParticle_p4(igen_particle).DeltaR(MuonOS);
      	  mc_index1 = igen_particle;
	      }
	
      	if(Ntp->MCParticle_p4(igen_particle).DeltaR(MuonSS1) < Dr2){
      	  Dr2 = Ntp->MCParticle_p4(igen_particle).DeltaR(MuonSS1);
	        mc_index2 = igen_particle;
      	}
	
    	  if(Ntp->MCParticle_p4(igen_particle).DeltaR(MuonSS2) < Dr3){
	        Dr3 = Ntp->MCParticle_p4(igen_particle).DeltaR(MuonSS2);
      	  mc_index3 = igen_particle;
	      }
      }

      Muon1DRToTruth.at(t).Fill(Dr1,1);
      Muon2DRToTruth.at(t).Fill(Dr2,1);
      Muon3DRToTruth.at(t).Fill(Dr3,1);


      // Truth muons momenta
      TLorentzVector GenMu_OS   = Ntp->MCParticle_p4(mc_index1); 
      TLorentzVector GenMu_SS1  = Ntp->MCParticle_p4(mc_index2);
      TLorentzVector GenMu_SS2  = Ntp->MCParticle_p4(mc_index3);


      if(Ntp->MCParticle_midx(mc_index1)>-0.5){
        IDOriginOfOSMuon.at(t).Fill(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(mc_index1)),1);

        if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(mc_index1))) == 221 );
        //	std::cout<<"This isEta  "<<std::endl;	Ntp->printMCDecayChainOfEvent(true, true, true, true);
        //	std::cout<< "\n\n\n\n\n\n";
      
        if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(mc_index1))) == 331 );
        //	std::cout<<"This isEtaPrime  "<<std::endl;	Ntp->printMCDecayChainOfEvent(true, true, true, true);
        //	std::cout<< "\n\n\n\n\n\n";
      }
      
      //To get lorentz vector for eta and eta prime
      
      TLorentzVector EtaLV(0,0,0,0);
      TLorentzVector EtaPrimeLV(0,0,0,0);
      TLorentzVector PhiLV(0,0,0,0);
      TLorentzVector OmegaLV(0,0,0,0);
      
      TLorentzVector EtaLV_reco(0,0,0,0);
      TLorentzVector EtaPrimeLV_reco(0,0,0,0);
      TLorentzVector PhiLV_reco(0,0,0,0);
      TLorentzVector OmegaLV_reco(0,0,0,0);
      
      TLorentzVector EtaLV_mixed(0,0,0,0);
      TLorentzVector EtaPrimeLV_mixed(0,0,0,0);
      TLorentzVector PhiLV_mixed(0,0,0,0);
      TLorentzVector OmegaLV_mixed(0,0,0,0);
      
      int mc_muon1idx(-1),mc_muon2idx(-1), mc_photon_idx(-1);
      int mc_prime_muon1idx(-1),mc_prime_muon2idx(-1), mc_prime_photon_idx(-1);
      int mc_phi_muon1idx(-1),mc_phi_muon2idx(-1), mc_phi_photon_idx(-1);
      int mc_omega_muon1idx(-1),mc_omega_muon2idx(-1), mc_omega_pion_idx(-1);
      
      int idx_fail_count(0);
      int idx_success_count(0);
      std::vector<int> child_vector;
      std::vector<int> child_prime_vector;
      
      for(int gen_part_index=0; gen_part_index < Ntp->NMCParticles(); gen_part_index++){        
          //std::cout<<"Particle PDGID is:"<< Ntp->MCParticle_pdgid(gen_part_index) << std::endl;
          //std::cout<<"Index of the mother particle is:"<< Ntp->MCParticle_midx(gen_part_index) << std::endl;
          if(Ntp->MCParticle_midx(gen_part_index)>-0.5){
          
            if(Ntp->MCParticle_pdgid(gen_part_index) == 22 || abs(Ntp->MCParticle_pdgid(gen_part_index)) == 13){
            
              if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(gen_part_index))) == 221 ){      //for eta
                //std::cout<<"Index of mother of"<< Ntp->MCParticle_pdgid(gen_part_index) <<" is " << Ntp->MCParticle_midx(gen_part_index) << std::endl;
                mc_photon_idx = (Ntp->MCParticle_pdgid(gen_part_index)==22)?gen_part_index:mc_photon_idx;
                mc_muon1idx = (Ntp->MCParticle_pdgid(gen_part_index)==13)?gen_part_index:mc_muon1idx;
                mc_muon2idx = (Ntp->MCParticle_pdgid(gen_part_index)<0)?gen_part_index:mc_muon2idx;
                idx_success_count=idx_success_count+1;
              }
              
              if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(gen_part_index))) == 331 ){       //for eta'
                //std::cout<<"Index of mother of "<< Ntp->MCParticle_pdgid(gen_part_index) <<" is " << Ntp->MCParticle_midx(gen_part_index) << std::endl;
                mc_prime_photon_idx = (Ntp->MCParticle_pdgid(gen_part_index)==22)?gen_part_index:mc_prime_photon_idx;
                mc_prime_muon1idx = (Ntp->MCParticle_pdgid(gen_part_index)==13)?gen_part_index:mc_prime_muon1idx;
                mc_prime_muon2idx = (Ntp->MCParticle_pdgid(gen_part_index)<0)?gen_part_index:mc_prime_muon2idx;
              }
              
              if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(gen_part_index))) == 333 ){        //for phi
                //std::cout<<"Index of mother of "<< Ntp->MCParticle_pdgid(gen_part_index) <<" is " << Ntp->MCParticle_midx(gen_part_index) << std::endl;
                mc_phi_muon1idx = (Ntp->MCParticle_pdgid(gen_part_index)==22)?gen_part_index:mc_phi_muon1idx;
                mc_phi_muon2idx = (Ntp->MCParticle_pdgid(gen_part_index)==13)?gen_part_index:mc_phi_muon2idx;
                mc_phi_photon_idx = (Ntp->MCParticle_pdgid(gen_part_index)<0)?gen_part_index:mc_phi_photon_idx;
              }
            }//end of MuMuGamma
            
            if(Ntp->MCParticle_pdgid(gen_part_index) == -111){
              std::cout << "There is pi0 antiparticle" << std::endl;
            }
            if(Ntp->MCParticle_pdgid(gen_part_index) == 111 || abs(Ntp->MCParticle_pdgid(gen_part_index)) == 13){
            
              if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(gen_part_index))) == 223 ){      //for Omega meson (not baryon)
                //std::cout<<"Index of mother of"<< Ntp->MCParticle_pdgid(gen_part_index) <<" is " << Ntp->MCParticle_midx(gen_part_index) << std::endl;
                mc_omega_muon1idx = (Ntp->MCParticle_pdgid(gen_part_index)==111)?gen_part_index:mc_omega_muon1idx;
                mc_omega_muon2idx = (Ntp->MCParticle_pdgid(gen_part_index)==13)?gen_part_index:mc_omega_muon2idx;
                mc_omega_pion_idx = (Ntp->MCParticle_pdgid(gen_part_index)<0)?gen_part_index:mc_omega_pion_idx;
              }
              
            }//end of MuMuPi0
            
          }//end of (gen_part_index)>-0.5
          
          
          //if(Ntp->MCParticle_midx(gen_part_index)<0){
            //std::cout << "The energy of the failed photon is " << Ntp->MCParticle_p4(gen_part_index).E() << std::endl;
            //idx_fail_count=idx_fail_count+1;
          //}
          
          
      }//end of for loop
      
      
      TLorentzVector eta_ph(0,0,0,0);
      TLorentzVector etaprime_ph(0,0,0,0);
      TLorentzVector phi_ph(0,0,0,0);
      TLorentzVector omega_pi(0,0,0,0);
      
      TLorentzVector eta_ph_reco(0,0,0,0);
      TLorentzVector etaprime_ph_reco(0,0,0,0);
      TLorentzVector phi_ph_reco(0,0,0,0);
      TLorentzVector omega_pi_reco(0,0,0,0);
      
      
      std::vector<unsigned int> pass_list{0,0,0,0};//stores eta_pass(0), etaprime_pass(0), phi_pass(0), omega_pass(0)
      
      if(mc_photon_idx>-0.5 && mc_muon1idx>-0.5 && mc_muon2idx>-0.5){//fill eta mass
        if((Ntp->MCParticle_midx(mc_photon_idx)==Ntp->MCParticle_midx(mc_muon1idx))&&(Ntp->MCParticle_midx(mc_muon1idx)==Ntp->MCParticle_midx(mc_muon2idx))){
          //std::cout<<"Successful 221"<< std::endl;
          EtaLV = Ntp->MCParticle_p4(mc_photon_idx) + Ntp->MCParticle_p4(mc_muon1idx) + Ntp->MCParticle_p4(mc_muon2idx);
          if(EtaLV.M()>0.52&&EtaLV.M()<0.57){
            PhotonSpectrum.at(t).Fill(Ntp->MCParticle_p4(mc_photon_idx).E(),Ntp->MCParticle_p4(mc_photon_idx).Eta());
            ChildCount.at(t).Fill(Ntp->MCParticle_childpdgid(Ntp->MCParticle_midx(mc_muon1idx)).size(),1);
            MassEta.at(t).Fill(EtaLV.M(),1);
            pass_list.at(0)=1;
            eta_ph=Ntp->MCParticle_p4(mc_photon_idx);
            std::cout << "The ID of the event is " << id << std::endl;
          }
        }
      }
      
      if(mc_prime_photon_idx>-0.5 && mc_prime_muon1idx>-0.5 && mc_prime_muon2idx>-0.5){//fill eta' mass
        if((Ntp->MCParticle_midx(mc_prime_photon_idx)==Ntp->MCParticle_midx(mc_prime_muon1idx))&&(Ntp->MCParticle_midx(mc_prime_muon1idx)==Ntp->MCParticle_midx(mc_prime_muon2idx))){
          //std::cout<<"Successful 331"<< std::endl;
          EtaPrimeLV = Ntp->MCParticle_p4(mc_prime_photon_idx) + Ntp->MCParticle_p4(mc_prime_muon1idx) + Ntp->MCParticle_p4(mc_prime_muon2idx);
          //std::cout<<"Mass of Eta Prime is: "<< EtaPrimeLV.M() << std::endl;
          if(EtaPrimeLV.M()>0.93&&EtaPrimeLV.M()<0.97){
            PhotonSpectrumPrime.at(t).Fill(Ntp->MCParticle_p4(mc_prime_photon_idx).E(),Ntp->MCParticle_p4(mc_prime_photon_idx).Eta());
            ChildCountPrime.at(t).Fill(Ntp->MCParticle_childpdgid(Ntp->MCParticle_midx(mc_prime_muon1idx)).size(),1);
            MassEtaPrime.at(t).Fill(EtaPrimeLV.M(),1);
            pass_list.at(1)=1;
            etaprime_ph=Ntp->MCParticle_p4(mc_prime_photon_idx);
            std::cout << "The ID of the event is " << id << std::endl;
          }
        }
      }
      
      if(mc_phi_photon_idx>-0.5 && mc_phi_muon1idx>-0.5 && mc_phi_muon2idx>-0.5){//fill phi mass
        if((Ntp->MCParticle_midx(mc_phi_photon_idx)==Ntp->MCParticle_midx(mc_phi_muon1idx))&&(Ntp->MCParticle_midx(mc_phi_muon1idx)==Ntp->MCParticle_midx(mc_phi_muon2idx))){
          PhiLV = Ntp->MCParticle_p4(mc_phi_photon_idx) + Ntp->MCParticle_p4(mc_phi_muon1idx) + Ntp->MCParticle_p4(mc_phi_muon2idx);
          if(PhiLV.M()>1.005&&PhiLV.M()<1.030){
            PhotonSpectrumPhi.at(t).Fill(Ntp->MCParticle_p4(mc_phi_photon_idx).E(),Ntp->MCParticle_p4(mc_phi_photon_idx).Eta());
            MassPhi.at(t).Fill(PhiLV.M(),1);
            pass_list.at(2)=1;
            phi_ph=Ntp->MCParticle_p4(mc_phi_photon_idx);
            std::cout << "The ID of the event is " << id << std::endl;
          }
        }
      }
      
      if(mc_omega_pion_idx>-0.5 && mc_omega_muon1idx>-0.5 && mc_omega_muon2idx>-0.5){//fill omega mass
        if((Ntp->MCParticle_midx(mc_omega_pion_idx)==Ntp->MCParticle_midx(mc_omega_muon1idx))&&(Ntp->MCParticle_midx(mc_omega_muon1idx)==Ntp->MCParticle_midx(mc_omega_muon2idx))){
          OmegaLV = Ntp->MCParticle_p4(mc_omega_pion_idx) + Ntp->MCParticle_p4(mc_omega_muon1idx) + Ntp->MCParticle_p4(mc_omega_muon2idx);
          if(OmegaLV.M()>0.770&&OmegaLV.M()<0.795){
            MassOmega.at(t).Fill(OmegaLV.M(),1);
            pass_list.at(3)=1;
            omega_pi=Ntp->MCParticle_p4(mc_omega_pion_idx);
            std::cout << "The ID of the event is " << id << std::endl;
          }
        }
      }
      
      
      
      int mc_idx_sum = (mc_photon_idx>-0.5) + (mc_muon1idx>-0.5) + (mc_muon2idx>-0.5);
      PassedCount.at(t).Fill(mc_idx_sum,1);      
      
      
      //To find the reco photons:
      
      
      TLorentzVector TauLV1 = Ntp->Muon_P4(os_mu_idx)  + Ntp->Muon_P4(ss1_mu_idx) + Ntp->Muon_P4(ss2_mu_idx);
      
      std::vector<int> recoph_index{-1,-1,-1,-1};//reco indices: photon from eta, eta', phi, and then pion_index
      
      if(pass_list.at(0)==1||pass_list.at(1)==1||pass_list.at(2)==1){
        double Dr1(99.),Dr2(99.),Dr3(99.);
        for(int reco_photon=0; reco_photon < Ntp->NPhotons(); reco_photon++){
          //std::cout<<"Energy of photon no. "<<reco_photon+1<<" is "<<Ntp->Photon_p4(reco_photon).E()<<std::endl;
          if(Ntp->Photon_p4(reco_photon).DeltaR(TauLV1) < 2){
            DeltaRPhotontoTau.at(t).Fill(4,1);
          }
          if(Ntp->Photon_p4(reco_photon).DeltaR(TauLV1) < 1){
            DeltaRPhotontoTau.at(t).Fill(3,1);
          }
          if(Ntp->Photon_p4(reco_photon).DeltaR(TauLV1) < 0.5){
            DeltaRPhotontoTau.at(t).Fill(2,1);
          }
          if(Ntp->Photon_p4(reco_photon).DeltaR(TauLV1) < 0.1){
            DeltaRPhotontoTau.at(t).Fill(1,1);
          }
          
          
          if(Ntp->Photon_p4(reco_photon).DeltaR(eta_ph) < Dr1){
            Dr1 = Ntp->Photon_p4(reco_photon).DeltaR(eta_ph);
            recoph_index.at(0) = reco_photon;
          }
          if(Ntp->Photon_p4(reco_photon).DeltaR(etaprime_ph) < Dr2){
            Dr2 = Ntp->Photon_p4(reco_photon).DeltaR(etaprime_ph);
            recoph_index.at(1) = reco_photon;
          }
          if(Ntp->Photon_p4(reco_photon).DeltaR(phi_ph) < Dr3){
            Dr3 = Ntp->Photon_p4(reco_photon).DeltaR(phi_ph);
            recoph_index.at(2) = reco_photon;
          }
        }
      }
      
      //std::cout<<"pass_list is:"<< pass_list.at(0)<< " " << pass_list.at(1)<< " " <<pass_list.at(2)<< " " <<pass_list.at(3)<< " " <<std::endl;
      //std::cout<<"recoph_index is:"<< recoph_index.at(0)<< " " << recoph_index.at(1)<< " " <<recoph_index.at(2)<< " " <<recoph_index.at(3)<< " " <<std::endl;
      
      std::vector<unsigned int> genidx{mc_index1,mc_index2,mc_index3}, recomu_idx{os_mu_idx,ss1_mu_idx,ss2_mu_idx};
      
      if(pass_list.at(0)==1&&recoph_index.at(0)>-0.5){
        std::cout<<"Whether eta photon has pixel seed: "<< Ntp->Photon_hasPixelSeed(recoph_index.at(0)) << std::endl;
        std::cout<<"Whether eta photon has Conversion Tracks: "<< Ntp->Photon_hasConversionTracks(recoph_index.at(0)) << std::endl;
        std::cout<<"Whether eta photon is PF: "<< Ntp->Photon_isPF(recoph_index.at(0)) << std::endl;
        Photon_hasPixelSeed_Eta.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(0)),1);
        Photon_hasPixelSeed.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(0)),1);
        Photon_hasConversionTracks.at(t).Fill(Ntp->Photon_hasConversionTracks(recoph_index.at(0)),1);
        Photon_isPF.at(t).Fill(Ntp->Photon_isPF(recoph_index.at(0)),1);
        
        eta_ph_reco=Ntp->Photon_p4(recoph_index.at(0));
        DeltaEnergyPhoton.at(t).Fill(abs(eta_ph_reco.E()-eta_ph.E()),1);
        int match1(-1), match2(-1);//finding the corresponding muons
        for(int idx_match =0; idx_match < 3; idx_match++){
          match1=genidx.at(idx_match)==mc_muon1idx?idx_match:match1;
          match2=genidx.at(idx_match)==mc_muon2idx?idx_match:match2;
        }
        if(match1>-0.5&&match2>-0.5){
          PhotonDRToTruth_Eta.at(t).Fill(eta_ph_reco.DeltaR(eta_ph),1);
          EtaLV_reco = Ntp->Muon_P4(recomu_idx.at(match1)) + Ntp->Muon_P4(recomu_idx.at(match2)) + eta_ph_reco;
          MassEtaReco.at(t).Fill(EtaLV_reco.M(),1);
          std::cout<<"Mass of reco eta is "<< EtaLV_reco.M()<< std::endl;
        }
        EtaLV_mixed = Ntp->MCParticle_p4(mc_muon1idx) + Ntp->MCParticle_p4(mc_muon2idx) + eta_ph_reco;
        MassEtaMixed.at(t).Fill(EtaLV_mixed.M(),1);
        std::cout<<"Mass of mixed-reco eta is "<< EtaLV_mixed.M()<< std::endl;
      }
      
      
      if(pass_list.at(1)==1&&recoph_index.at(1)>-0.5){
        std::cout<<"Whether eta' photon has pixel seed: "<< Ntp->Photon_hasPixelSeed(recoph_index.at(1)) << std::endl;
        std::cout<<"Whether eta' photon has Conversion Tracks: "<< Ntp->Photon_hasConversionTracks(recoph_index.at(1)) << std::endl;
        std::cout<<"Whether eta' photon is PF: "<< Ntp->Photon_isPF(recoph_index.at(1)) << std::endl;
        Photon_hasPixelSeed_EtaPrime.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(1)),1);
        Photon_hasPixelSeed.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(1)),1);
        Photon_hasConversionTracks.at(t).Fill(Ntp->Photon_hasConversionTracks(recoph_index.at(1)),1);
        Photon_isPF.at(t).Fill(Ntp->Photon_isPF(recoph_index.at(1)),1);
        
        etaprime_ph_reco=Ntp->Photon_p4(recoph_index.at(1));
        DeltaEnergyPhoton.at(t).Fill(abs(etaprime_ph_reco.E()-etaprime_ph.E()),1);
        int match1(-1), match2(-1);//finding the corresponding muons
        for(int idx_match =0; idx_match < 3; idx_match++){
          match1=genidx.at(idx_match)==mc_prime_muon1idx?idx_match:match1;
          match2=genidx.at(idx_match)==mc_prime_muon2idx?idx_match:match2;
        }
        if(match1>-0.5&&match2>-0.5){
          PhotonDRToTruth_EtaPrime.at(t).Fill(etaprime_ph_reco.DeltaR(etaprime_ph),1);
          EtaPrimeLV_reco = Ntp->Muon_P4(recomu_idx.at(match1)) + Ntp->Muon_P4(recomu_idx.at(match2)) + etaprime_ph_reco;
          MassEtaPrimeReco.at(t).Fill(EtaPrimeLV_reco.M(),1);
          std::cout<<"Mass of reco eta prime is "<< EtaPrimeLV_reco.M()<< std::endl;
        }
        EtaPrimeLV_mixed = Ntp->MCParticle_p4(mc_prime_muon1idx) + Ntp->MCParticle_p4(mc_prime_muon2idx) + etaprime_ph_reco;
        MassEtaPrimeMixed.at(t).Fill(EtaPrimeLV_mixed.M(),1);
        std::cout<<"Mass of mixed-reco eta' is "<< EtaPrimeLV_mixed.M()<< std::endl;
      }
      
      
      if(pass_list.at(2)==1&&recoph_index.at(2)>-0.5){
        std::cout<<"Whether phi photon has pixel seed: "<< Ntp->Photon_hasPixelSeed(recoph_index.at(2)) << std::endl;
        std::cout<<"Whether phi photon has Conversion Tracks: "<< Ntp->Photon_hasConversionTracks(recoph_index.at(2)) << std::endl;
        std::cout<<"Whether phi photon is PF: "<< Ntp->Photon_isPF(recoph_index.at(2)) << std::endl;
        Photon_hasPixelSeed_Phi.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(2)),1);
        Photon_hasPixelSeed.at(t).Fill(Ntp->Photon_hasPixelSeed(recoph_index.at(2)),1);
        Photon_hasConversionTracks.at(t).Fill(Ntp->Photon_hasConversionTracks(recoph_index.at(2)),1);
        Photon_isPF.at(t).Fill(Ntp->Photon_isPF(recoph_index.at(2)),1);
        
        phi_ph_reco=Ntp->Photon_p4(recoph_index.at(2));
        DeltaEnergyPhoton.at(t).Fill(abs(phi_ph_reco.E()-phi_ph.E()),1);
        int match1(-1), match2(-1);//finding the corresponding muons
        for(int idx_match =0; idx_match < 3; idx_match++){
          match1=genidx.at(idx_match)==mc_phi_muon1idx?idx_match:match1;
          match2=genidx.at(idx_match)==mc_phi_muon2idx?idx_match:match2;
        }
        if(match1>-0.5&&match2>-0.5){
          PhotonDRToTruth_Phi.at(t).Fill(phi_ph_reco.DeltaR(phi_ph),1);
          PhiLV_reco = Ntp->Muon_P4(recomu_idx.at(match1)) + Ntp->Muon_P4(recomu_idx.at(match2)) + phi_ph_reco;
          MassPhiReco.at(t).Fill(PhiLV_reco.M(),1);
          std::cout<<"Mass of reco phi is "<< PhiLV_reco.M()<< std::endl;
        }
        PhiLV_mixed = Ntp->MCParticle_p4(mc_phi_muon1idx) + Ntp->MCParticle_p4(mc_phi_muon2idx) + phi_ph_reco;
        MassPhiMixed.at(t).Fill(PhiLV_mixed.M(),1);
        std::cout<<"Mass of mixed-reco phi is "<< PhiLV_mixed.M()<< std::endl;
      }
      
      if((pass_list.at(0)==1||pass_list.at(1)==1||pass_list.at(2)==1)&&recoph_index.at(0)==-1&&recoph_index.at(1)==-1&&recoph_index.at(2)==-1){
        PhotonRecoSuccess.at(t).Fill(-1,1);
      }
      if((pass_list.at(0)==1||pass_list.at(1)==1||pass_list.at(2)==1)&&(recoph_index.at(0)>-0.5||recoph_index.at(1)>-0.5||recoph_index.at(2)>-0.5)){
        PhotonRecoSuccess.at(t).Fill(1,1);
      }
      
      
    }

    //----



    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);







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



    TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(0),1);
    TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(1),1);
    TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(2),1);

    //    std::cout<<"  Match 1 "<< Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(0) << std::endl;
    //    std::cout<<"  Match 2 "<< Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(1) << std::endl;
    //    std::cout<<"  Match 3 "<< Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(2) << std::endl;



    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){
	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;

	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
      }
    }
  }
}


void  MCBackgroundStudy::Finish(){
    std::cout<<" 3 "<< std::endl;
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