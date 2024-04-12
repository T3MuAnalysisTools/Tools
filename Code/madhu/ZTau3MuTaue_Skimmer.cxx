#include "ZTau3MuTaue_Skimmer.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTaue_Skimmer::ZTau3MuTaue_Skimmer(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTaue_Skimmer::~ZTau3MuTaue_Skimmer(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTaue_Skimmer::Configure(){

  //  Mini tree for limit extraction

  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  T3MMiniTree->Branch("LumiScale",&LumiScale);
  
  T3MMiniTree->Branch("var_TripletPT",&var_TripletPT);
  T3MMiniTree->Branch("var_Tau3MuIsolation",&var_Tau3MuIsolation);
  T3MMiniTree->Branch("var_Electron_pT",&var_Electron_pT);
  T3MMiniTree->Branch("var_VisMass",&var_VisMass);
  T3MMiniTree->Branch("var_mu1_pT",&var_mu1_pT);
  T3MMiniTree->Branch("var_mu2_pT",&var_mu2_pT);
  T3MMiniTree->Branch("var_mu3_pT",&var_mu3_pT);
  
  
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==WhetherZTTDecayFound) cut.at(WhetherZTTDecayFound)=1;
    if(i==L1_TriggerOk)       cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)      cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==L1_TriggerOk){
      title.at(i)="L1 Trigger ";
      hlabel="L1 Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT_TriggerOk){
      title.at(i)="HLT Trigger ";
      hlabel="HLT Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==WhetherZTTDecayFound){
      title.at(i)="Passed fiducial cuts";
      hlabel="Passed fiducial cuts ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_WhetherZTTDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_WhetherZTTDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate (3,3,2 GeV,  $|\\eta| < 2.4$, dz($\\mu_{i} , \\mu_{j}$)$<$0.5, dR($\\mu_{i} , \\mu_{j}$)$<$0.8, $\\Sigma \\mu_{charge}$ = +-1)";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }



  }
  // Setup NPassed Histogams




  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  ElectronSumIsolation=HConfig.GetTH1D(Name+"_ElectronSumIsolation","ElectronSumIsolation",50,0.,10,"I= neutralH + chargedH + photon Iso, GeV","Events ");
  
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",70,0.,140,"M_{#tau(#mu) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",40,1.1,2.2,"M_{3#mu}, GeV","Events");
  PairMass_OppositeSign_dR12=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR12","PairMass_OppositeSign_dR12",40,0.2,2.,"M_{1}, GeV (OS - SS dR sorted)","Events");
  PairMass_OppositeSign_dR13=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR13","PairMass_OppositeSign_dR13",40,0.2,2.,"M_{2}, GeV (OS - SS dR sorted)","Events");





  matched_pdgId=HConfig.GetTH1D(Name+"_matched_pdgId","matched_pdgId",25,-0.5,24.5,"pdgID MC matched","Events");
  matched_dR=HConfig.GetTH1D(Name+"_matched_dR","matched_dR",50,-0.1,0.5,"#Delta R(MC-RECO) Object opposite to #tau_{3#mu}","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{1} ","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{2} ","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{3} ","Events");
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,0.5,"#Delta R Truth #tau's prods ","Events");
  
  dR_betweenTruth_NeutrinoGuess=HConfig.GetTH1D(Name+"_dR_betweenTruth_NeutrinoGuess","dR_betweenTruth_NeutrinoGuess",20,0,0.5,"#Delta R: Truth to Guessed #nu","Events");
  dR_betweenTruth_Tau=HConfig.GetTH1D(Name+"_dR_betweenTruth_Tau","dR_betweenTruth_Tau",20,0,0.5,"#Delta R: Truth to Guessed #tau","Events");
  Z_Pt=HConfig.GetTH1D(Name+"_Z_Pt","Z_Pt",50,0,70,"Z_{pT}","Events");
  OS_vs_3mu_trigger=HConfig.GetTH2D(Name+"_OS_vs_3mu_trigger","OS_vs_3mu_trigger",3,-0.5,2.5,2,-0.5,1.5,"Whether 3mu Triggered","Whether OS #tau Triggered");
  
  MET_Et=HConfig.GetTH1D(Name+"_MET_Et","MET_Et",100,0.0,100.0,"MET Et, GeV","Events");
  MET_Phi=HConfig.GetTH1D(Name+"_MET_Phi","MET_Phi",20,-3.2,3.2,"MET #Phi ","Events");
  
  Selection_Cut_3mu_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Pt","Selection_Cut_3mu_Pt",100,0,50.0,"3#mu p_{T}, GeV","Events");
  Selection_Cut_3mu_Rel_Iso=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Rel_Iso","Selection_Cut_3mu_Rel_Iso",50,0,1.1,"3 #mu Relative Isolation, p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","Events");
  Selection_Cut_elect_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_elect_Pt","Selection_Cut_elect_Pt",100,0,50.0,"e p_{T}, GeV","Events");
  Selection_Cut_elect_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_elect_Eta","Selection_Cut_elect_Eta",30,0,3.0,"e |#eta|","Events");
  Selection_Cut_elect_DeltaR_3mu=HConfig.GetTH1D(Name+"_Selection_Cut_elect_DeltaR_3mu","Selection_Cut_elect_DeltaR_3mu",60,0,1.2,"#Delta R (e-3#mu)","Events");
  Selection_Cut_Vis_InvM=HConfig.GetTH1D(Name+"_Selection_Cut_Vis_InvM","Selection_Cut_Vis_InvM",75,0,150.0,"M_{e + #tau(3#mu)}, GeV (visible mass)","Events");
  
  Selection_Cut_Mu1_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_P","Selection_Cut_Mu1_P",160,0.0,80.0,"#mu_{1} p, GeV","Events");
  Selection_Cut_Mu1_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_Eta","Selection_Cut_Mu1_Eta",30,0,3.14,"#mu_{1} |#eta|","Events");
  Selection_Cut_Mu1_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_before","Selection_Cut_Mu1_p_eta_before",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after","Selection_Cut_Mu1_p_eta_after",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after_reco","Selection_Cut_Mu1_p_eta_after_reco",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  
  Selection_Cut_Mu2_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_P","Selection_Cut_Mu2_P",160,0.0,80.0,"#mu_{2} p, GeV","Events");
  Selection_Cut_Mu2_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_Eta","Selection_Cut_Mu2_Eta",30,0,3.14,"#mu_{2} |#eta|","Events");
  Selection_Cut_Mu2_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_before","Selection_Cut_Mu2_p_eta_before",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after","Selection_Cut_Mu2_p_eta_after",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after_reco","Selection_Cut_Mu2_p_eta_after_reco",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  
  Selection_Cut_Mu3_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_P","Selection_Cut_Mu3_P",160,0.0,80.0,"#mu_{3} p, GeV","Events");
  Selection_Cut_Mu3_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_Eta","Selection_Cut_Mu3_Eta",30,0,3.14,"#mu_{3} |#eta|","Events");
  Selection_Cut_Mu3_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_before","Selection_Cut_Mu3_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after","Selection_Cut_Mu3_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after_reco","Selection_Cut_Mu3_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  
  Selection_Cut_El_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_El_Pt","Selection_Cut_El_Pt",40,0.0,80.0,"e p_{T}, GeV","Events");
  Selection_Cut_El_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_El_Eta","Selection_Cut_El_Eta",30,0,3.14,"e |#eta|","Events");
  Selection_Cut_El_pt_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_before","Selection_Cut_El_pt_eta_before",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  Selection_Cut_El_pt_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_after","Selection_Cut_El_pt_eta_after",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  Selection_Cut_El_pt_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_after_reco","Selection_Cut_El_pt_eta_after_reco",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  
  Selection_Cut_Mu1_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_dR","Selection_Cut_Mu1_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu2_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_dR","Selection_Cut_Mu2_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu3_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_dR","Selection_Cut_Mu3_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_El_dR=HConfig.GetTH1D(Name+"_Selection_Cut_El_dR","Selection_Cut_El_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu1_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu1_dR_large_scale","Selection_Cut_Mu1_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu2_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu2_dR_large_scale","Selection_Cut_Mu2_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu3_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu3_dR_large_scale","Selection_Cut_Mu3_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_El_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_El_dR_large_scale","Selection_Cut_El_dR_large_scale",100,0,2.0,"#Delta R","Events");
  
  Selection_Cut_RecoMu_P=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_P","Selection_Cut_RecoMu_P",100,0.0,5.0,"#mu p, GeV","Events");
  Selection_Cut_RecoMu_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_Eta","Selection_Cut_RecoMu_Eta",50,2,3.0,"#mu |#eta|","Events");
  Selection_Cut_RecoEl_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_RecoEl_Pt","Selection_Cut_RecoEl_Pt",100,0.0,10.0,"e p_{T}, GeV","Events");
  Selection_Cut_RecoEl_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoEl_Eta","Selection_Cut_RecoEl_Eta",50,2,3.0,"e |#eta|","Events");
  
  Electron_Isolation_relative=HConfig.GetTH1D(Name+"_Electron_Isolation_relative","Electron_Isolation_relative",50,0.0,1.0,"e relative isolation","Events");
  Electron_Isolation_trackIso=HConfig.GetTH1D(Name+"_Electron_Isolation_trackIso","Electron_Isolation_trackIso",50,0.0,5.0,"e track isolation","Events");
  Electron_Isolation_puppiPhotonIso=HConfig.GetTH1D(Name+"_Electron_Isolation_puppiPhotonIso","Electron_Isolation_puppiPhotonIso",50,0.0,5.0,"e puppiPhoton isolation","Events");
  Electron_Isolation_puppiNeutralHadronIso=HConfig.GetTH1D(Name+"_Electron_Isolation_puppiNeutralHadronIso","Electron_Isolation_puppiNeutralHadronIso",50,0.0,5.0,"e puppiNeutralHadron isolation","Events");
  Electron_Isolation_puppiChargedHadronIso=HConfig.GetTH1D(Name+"_Electron_Isolation_puppiChargedHadronIso","Electron_Isolation_puppiChargedHadronIso",100,0.0,3.0,"e puppiChargedHadron isolation","Events");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTaue_Skimmer::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output


  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&ElectronSumIsolation);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);
  Extradist1d.push_back(&matched_pdgId);
  Extradist1d.push_back(&matched_dR);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);
  Extradist1d.push_back(&dR_betweenTruth_VisibleTaus);
  Extradist1d.push_back(&PairMass_OppositeSign_dR12);
  Extradist1d.push_back(&PairMass_OppositeSign_dR13);
  
  Extradist1d.push_back(&dR_betweenTruth_NeutrinoGuess);
  Extradist1d.push_back(&dR_betweenTruth_Tau);
  Extradist1d.push_back(&Z_Pt);
  Extradist2d.push_back(&OS_vs_3mu_trigger);
  
  Extradist1d.push_back(&MET_Et);
  Extradist1d.push_back(&MET_Phi);
  
  Extradist1d.push_back(&Selection_Cut_3mu_Pt);
  Extradist1d.push_back(&Selection_Cut_3mu_Rel_Iso);
  Extradist1d.push_back(&Selection_Cut_elect_Pt);
  Extradist1d.push_back(&Selection_Cut_elect_Eta);
  Extradist1d.push_back(&Selection_Cut_elect_DeltaR_3mu);
  Extradist1d.push_back(&Selection_Cut_Vis_InvM);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_P);
  Extradist1d.push_back(&Selection_Cut_Mu1_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu2_P);
  Extradist1d.push_back(&Selection_Cut_Mu2_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu3_P);
  Extradist1d.push_back(&Selection_Cut_Mu3_Eta);
  Extradist1d.push_back(&Selection_Cut_El_Pt);
  Extradist1d.push_back(&Selection_Cut_El_Eta);
  
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_before);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_after);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_after_reco);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR);
  Extradist1d.push_back(&Selection_Cut_El_dR);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_El_dR_large_scale);
  
  Extradist1d.push_back(&Selection_Cut_RecoMu_P);
  Extradist1d.push_back(&Selection_Cut_RecoMu_Eta);
  Extradist1d.push_back(&Selection_Cut_RecoEl_Pt);
  Extradist1d.push_back(&Selection_Cut_RecoEl_Eta);
  
  Extradist1d.push_back(&Electron_Isolation_relative);
  Extradist1d.push_back(&Electron_Isolation_trackIso);
  Extradist1d.push_back(&Electron_Isolation_puppiPhotonIso);
  Extradist1d.push_back(&Electron_Isolation_puppiNeutralHadronIso);
  Extradist1d.push_back(&Electron_Isolation_puppiChargedHadronIso);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTaue_Skimmer::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection



  bool HLTOk(false);
  bool L1Ok(false);
  bool SingleMuFired(false);
  //bool DoubleMu0Fired(false);
  //bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);
  bool HLT_OppositeSide(false);
  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
    if(HLTName.Contains("HLT_Ele32_WPTight_Gsf") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
    if(HLTName.Contains("HLT_Ele15_WPLoose_Gsf") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
  }
  
  random_num = rndm.Rndm();
  
  //std::cout<<"  Event:   "  << std::endl;
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    
    if( ( L1TriggerName.Contains("L1_SingleMu22") || L1TriggerName.Contains("L1_SingleMu25") )  && Ntp->L1Decision(il1)) { SingleMuFired = true; }
    if( ( L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") || L1TriggerName.Contains("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4") ) && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    if( ( L1TriggerName.Contains("L1_DoubleMu4p5_SQ_OS_dR_Max1p2") )  && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true;}
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
    if( ( L1TriggerName.Contains("L1_DoubleMu_15_7") )  && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    if( ( L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") || L1TriggerName.Contains("") ) && Ntp->L1Decision(il1)) { TripleMuFired = true; }
  }
  
  value.at(L1_TriggerOk)=0;  
  value.at(HLT_TriggerOk)= 0 ;
  //if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (SingleMuFired || DoubleMuFired || TripleMuFired) L1Ok = true;

  value.at(L1_TriggerOk)=(L1Ok);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  
  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));



  value.at(SignalCandidate) = Ntp->NThreeMuons();
  pass.at(SignalCandidate)=(value.at(SignalCandidate)>=cut.at(SignalCandidate));

  int  signal_idx=-1;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  TLorentzVector MC_NeutrinoSum_LV(0.,0.,0.,0.);
  
  bool WhetherSignalMC = id==210||id==210231||id==210232||id==210233;
  if(WhetherSignalMC){
  
  int Whether_decay_found(0);
  int TausFromZ_Count(0);
  int tau_3mu_idx(-1);
  int tau_mu_idx(-1);
  int tau_e_idx(-1);
  int tau_h_idx(-1);
  for(unsigned int imc =0; imc< Ntp->NMCParticles(); imc++){
    if(abs(Ntp->MCParticle_pdgid(imc)) == 15 && Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(imc) ) == 23){
      TausFromZ_Count++;
          
          //correcting for cases where tau decays multiple times
          int ChildIdx = imc;
          bool Whether_Another_Tau_Decay_In_Chain(false);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15) Whether_Another_Tau_Decay_In_Chain=true;
          }
          while (Whether_Another_Tau_Decay_In_Chain){
            Whether_Another_Tau_Decay_In_Chain=false;
            for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
              if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15){
                ChildIdx = Ntp->MCParticle_childidx(ChildIdx).at(i);
                Whether_Another_Tau_Decay_In_Chain=true;
              }
            }
          }
          
          int nmuons_temp(0);
          int nelectrons_temp(0);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            //std::cout<<"Child "<< i <<" of tau is :  "<< Ntp->MCParticle_childpdgid(ChildIdx).at(i) << std::endl;
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==13){
              nmuons_temp++;
            }
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==11){
              nelectrons_temp++;
            }
          }
          
          if(nmuons_temp==3){
            tau_3mu_idx=ChildIdx;
          }
          else if(nmuons_temp==1){
            tau_mu_idx=ChildIdx;
          }
          else if(nelectrons_temp==1){
            tau_e_idx=ChildIdx;
          }
          else{
            tau_h_idx=ChildIdx;
          }
          
          
      }
  }
    
  if(tau_3mu_idx>-0.5 && tau_e_idx>-0.5) Whether_decay_found=1;
  int Whether_decay_found_temp(0);
  if(tau_3mu_idx>-0.5) Whether_decay_found_temp=1;
  
  
  
  
  
  //Checking for alternate triggers
  bool HasSomeMuTrigger(false);
  bool HasIsoEleTrigger(false);
  bool HasSingleEle_L1Trigger(false);
  
  //L1T
  /*
  if(!HLTOk&&Whether_decay_found){
  
  std::cout<<" Events:  "<<std::endl;
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    
    if(Ntp->L1Decision(il1)&&L1TriggerName.Contains("L1_SingleEG")){
      //HasSingleEle_L1Trigger=true;
    }
    
  }
  
  //if(HasSingleEle_L1Trigger) std::cout<<" Single Ele L1 fired   "  << std::endl;
  //if(!HasSingleEle_L1Trigger) std::cout<<" Single Ele L1 Not fired   "  << std::endl;
  
  //If there's no single electron L1 trigger
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    
    if(Ntp->L1Decision(il1)&&!HasSingleEle_L1Trigger){
      std::cout<<"L1 fired:   "  << L1TriggerName  << std::endl;
    }
    
  }
  
  }
  */
  
  if(!HLTOk&&Whether_decay_found){
  
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(Ntp->L1Decision(il1)){
      //std::cout<<"L1 fired:   "  << L1TriggerName  << std::endl;
    }
    
  }
  
  }
  
  
  
  /*
  //HLT
  if(!HLTOk&&Whether_decay_found){
  //std::cout<<"No triple mu trigger for taue category:  " << std::endl;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    
    if(Ntp->HLTDecision(iTrigger)){
      HasSomeMuTrigger=true;
      if(HLTName.Contains("Ele")){
        HasIsoEleTrigger=true;
        //std::cout<<" HLTName:   "<< HLTName << std::endl;
      }
    }
  }
  
  
  OS_vs_3mu_trigger.at(t).Fill(2,HasSomeMuTrigger,1 );
  if(HasSomeMuTrigger){
    OS_vs_3mu_trigger.at(t).Fill(3,HasIsoEleTrigger,1 );
    //std::cout<<"Has some triggers: " << std::endl;
  }
  
  if(HasSomeMuTrigger&&!HasIsoEleTrigger){
    //std::cout<<" Has some mu trigger but no ele " << std::endl;
  }
  
  if(HasIsoEleTrigger){
    //std::cout<<" Has ele triggers " << std::endl;
  }
  
  
  }
  
  OS_vs_3mu_trigger.at(t).Fill(HLTOk,HLT_OppositeSide,1 );
  */
  
  
  
  
  
  
  //std::cout<<"  No of taus from Z:  "<< TausFromZ_Count << std::endl;
  //if(Whether_decay_found&&id==210233)      Ntp->printMCDecayChainOfEvent(true,true,true,true);
  
  TLorentzVector Mu1_LV;
  TLorentzVector Mu2_LV;
  TLorentzVector Mu3_LV;
  TLorentzVector Electron_LV;
  if(Whether_decay_found==1){
    std::vector<int> Sorted_MC_Indices = Ntp->SortedPtMuons_MC(Ntp->MCParticle_childidx(tau_3mu_idx));
    
    Mu1_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(0));
    Mu2_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(1));
    Mu3_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(2));
    
    for(int i = 0; i < Ntp->MCParticle_childpdgid(tau_e_idx).size(); i++){
      if(abs(Ntp->MCParticle_childpdgid(tau_e_idx).at(i))==11){
        Electron_LV=Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_e_idx).at(i));
      }
      if(abs(Ntp->MCParticle_childpdgid(tau_e_idx).at(i))==12||abs(Ntp->MCParticle_childpdgid(tau_e_idx).at(i))==14||abs(Ntp->MCParticle_childpdgid(tau_e_idx).at(i))==16){
        MC_NeutrinoSum_LV=MC_NeutrinoSum_LV+Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_e_idx).at(i));
      }
    }
    
    
    
    //std::cout << "New e Event" << std::endl;
    bool triggerCheck_os(false);
    //Trigger matching for the fourth leg, to be used in Skimmer
    if(signal_idx!=-1)
    {
              //Trigger Matching opposite side
              //if(fabs(Electron_LV.Eta())<2.41&&!HLTOk&&Electron_LV.Pt()>4.9)
              if(true)
                {
                  //std::cout << "tau e: No 3mu trig, yes OS trig" << std::endl;
                  vector<TLorentzVector> trigobjTriplet;
                  for (int i=0; i<Ntp->NTriggerObjects(); i++)
                    {
                      TString name = Ntp->TriggerObject_name(i);
                      //        if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
                      
                      TLorentzVector tmp;
                      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
                      trigobjTriplet.push_back(tmp);
                      
                      double dpT = fabs(Electron_LV.Pt()-tmp.Pt())/Electron_LV.Pt();
                      double dR = Electron_LV.DeltaR(tmp);
                      
                      if(dpT<0.1 && dR<0.05 && (name.Contains("hltEle15WPLoose1EcalIsoFilter")||name.Contains("hltEle5WPTightEcalIsoFilter")||name.Contains("hltL1sSingleEGor")||name.Contains("hltEle8EtFilter")) ){
                              //cout << " The trigger object is "<< name << " with dR: " << dR << " and dpT: "<< dpT << endl;
                              triggerCheck_os=true;
                              continue;
                      }
                    }
        	  
                }
                
    }
    //x-axis: 0: HLTOk not pass, 1: HLTOk pass, 2: what part of reconstructable object is trigger matched to an object
    
    //OS_vs_3mu_trigger.at(t).Fill(HLTOk,(fabs(Electron_LV.Eta())<2.41&&Electron_LV.Pt()>4.9),1 );
    //if(fabs(Electron_LV.Eta())<2.41&&!HLTOk&&Electron_LV.Pt()>4.9){
    //        OS_vs_3mu_trigger.at(t).Fill(2,triggerCheck_os,1 );
    //}
    
    OS_vs_3mu_trigger.at(t).Fill(HLTOk&&L1Ok,triggerCheck_os,1 );
    
  }
  
  double cut_Mu_Candidate_p=2.49;
  double cut_Mu_Candidate_eta=2.41;
  double cut_Tau_e_Candidate_p=4.9;
  double cut_Tau_e_Candidate_eta=2.41;
  
  bool var_Whether_decay_found = Whether_decay_found;
  bool var_Mu1_Candidate_p = (Mu1_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu1_Candidate_eta = (abs(Mu1_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Mu2_Candidate_p = (Mu2_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu2_Candidate_eta = (abs(Mu2_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Mu3_Candidate_p = (Mu3_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu3_Candidate_eta = (abs(Mu3_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Tau_e_Candidate_p = (Electron_LV.Pt()) > cut_Tau_e_Candidate_p;
  bool var_Tau_e_Candidate_eta = (abs(Electron_LV.Eta())) < cut_Tau_e_Candidate_eta;
  
  double dR1_max(99.0);
  double dR2_max(99.0);
  double dR3_max(99.0);
  double dR4_max(99.0);
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      if(Ntp->Muon_P4(imu).DeltaR(Mu1_LV)<dR1_max){
        dR1_max=Ntp->Muon_P4(imu).DeltaR(Mu1_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu2_LV)<dR2_max){
        dR2_max=Ntp->Muon_P4(imu).DeltaR(Mu2_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu3_LV)<dR3_max){
        dR3_max=Ntp->Muon_P4(imu).DeltaR(Mu3_LV);
      }
      
      Selection_Cut_RecoMu_P.at(t).Fill(Ntp->Muon_P4(imu).Vect().Mag(),1 );
      Selection_Cut_RecoMu_Eta.at(t).Fill(abs(Ntp->Muon_P4(imu).Eta()),1 );

    }
  
  for(unsigned int ie=0; ie < Ntp->NElectrons(); ie++)
    {
      if(Ntp->Electron_P4(ie).DeltaR(Electron_LV)<dR4_max){
        dR4_max=Ntp->Electron_P4(ie).DeltaR(Electron_LV);
      }
      
      Selection_Cut_RecoEl_Pt.at(t).Fill(Ntp->Electron_P4(ie).Vect().Mag(),1 );
      Selection_Cut_RecoEl_Eta.at(t).Fill(abs(Ntp->Electron_P4(ie).Eta()),1 );
      
    }
  
  bool var_Mu1_Candidate_recod = (dR1_max<0.01);
  bool var_Mu2_Candidate_recod = (dR2_max<0.01);
  bool var_Mu3_Candidate_recod = (dR3_max<0.01);
  bool var_Tau_e_Candidate_recod = (dR4_max<0.01);
  
  //value.at(WhetherZTTDecayFound)=var_Whether_decay_found&&var_Mu1_Candidate_p&&var_Mu1_Candidate_eta&&var_Mu2_Candidate_p&&var_Mu2_Candidate_eta&&var_Mu3_Candidate_p&&var_Mu3_Candidate_eta&&var_Tau_e_Candidate_p&&var_Tau_e_Candidate_eta&&var_Mu1_Candidate_recod&&var_Mu2_Candidate_recod&&var_Mu3_Candidate_recod&&var_Tau_e_Candidate_recod;
  value.at(WhetherZTTDecayFound)=var_Whether_decay_found;
  pass.at(WhetherZTTDecayFound)=(value.at(WhetherZTTDecayFound)==cut.at(WhetherZTTDecayFound));
  
  // This is to print out selected event content
  if(id==210231){
          
          //std::cout<<"Whether tau e decay found: "<< Whether_decay_found << ". Event id: "<<id<< std::endl;
          
          if(Whether_decay_found){
                  /*
                  std::cout<<"------------------------------- "<< std::endl;
                  std::cout<<"Event Content "<< std::endl;
                  std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Mu1_LV) << std::endl;
                  std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Mu2_LV) << std::endl;
                  std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Mu3_LV) << std::endl;
                  std::cout<<" idx4:  "<<Ntp->getMatchTruthIndex(Electron_LV) << std::endl;
                  Ntp->printMCDecayChainOfEvent(true, true, true, true);
                  std::cout<< "\n\n\n\n\n\n";
                  
                  for(unsigned int i=0; i < Ntp->NMCTaus(); i++){
                    for(unsigned int j=0; j < Ntp->NMCTauDecayProducts(i); j++){
                      std::cout<<"Product coming from index: "<<i<<" with product index: "<<j<<" : "<<Ntp->MCTauandProd_pdgid(i, j)<< std::endl;
                    }
                  }
                  */
                  /*
                  Selection_Cut_Mu1_P.at(t).Fill(Mu1_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu1_Eta.at(t).Fill(abs(Mu1_LV.Eta()),1 );
                  Selection_Cut_Mu2_P.at(t).Fill(Mu2_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu2_Eta.at(t).Fill(abs(Mu2_LV.Eta()),1 );
                  Selection_Cut_Mu3_P.at(t).Fill(Mu3_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu3_Eta.at(t).Fill(abs(Mu3_LV.Eta()),1 );
                  Selection_Cut_El_Pt.at(t).Fill(Electron_LV.Pt(),1 );
                  Selection_Cut_El_Eta.at(t).Fill(abs(Electron_LV.Eta()),1 );
                  
                  Selection_Cut_Mu1_p_eta_before.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                  if(pass.at(Mu1_Candidate_p)&&pass.at(Mu1_Candidate_eta)){
                    Selection_Cut_Mu1_p_eta_after.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    if(!pass.at(Mu1_Candidate_recod)){
                      Selection_Cut_Mu1_p_eta_after_reco.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_Mu2_p_eta_before.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                  if(pass.at(Mu2_Candidate_p)&&pass.at(Mu2_Candidate_eta)){
                    Selection_Cut_Mu2_p_eta_after.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    if(!pass.at(Mu2_Candidate_recod)){
                      Selection_Cut_Mu2_p_eta_after_reco.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_Mu3_p_eta_before.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                  if(pass.at(Mu3_Candidate_p)&&pass.at(Mu3_Candidate_eta)){
                    Selection_Cut_Mu3_p_eta_after.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    if(!pass.at(Mu3_Candidate_recod)){
                      Selection_Cut_Mu3_p_eta_after_reco.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_El_pt_eta_before.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                  if(pass.at(Tau_e_Candidate_p)&&pass.at(Tau_e_Candidate_eta)){
                    Selection_Cut_El_pt_eta_after.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                    if(!pass.at(Tau_e_Candidate_recod)){
                      Selection_Cut_El_pt_eta_after_reco.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                    }
                  }
                  }
                  */
                  Selection_Cut_Mu1_dR.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR.at(t).Fill(dR3_max,1 );
                  Selection_Cut_El_dR.at(t).Fill(dR4_max,1 );
                  
                  Selection_Cut_Mu1_dR_large_scale.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR_large_scale.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR_large_scale.at(t).Fill(dR3_max,1 );
                  Selection_Cut_El_dR_large_scale.at(t).Fill(dR4_max,1 );
                  
          }
          
  }
  
  }//if(id!=1)
  
  if(!WhetherSignalMC){
    value.at(WhetherZTTDecayFound)=1;
    pass.at(WhetherZTTDecayFound)=1;
  }
  

  

    double wobs=1;
    double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 

    

  }
}


void  ZTau3MuTaue_Skimmer::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis
  T3MFMiniTree = new TFile("T3MMiniTree_taue.root","recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();
  
  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





