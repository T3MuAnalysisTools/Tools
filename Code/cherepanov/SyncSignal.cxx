#include "SyncSignal.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

SyncSignal::SyncSignal(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.70),
  tauMaxMass_(1.82),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0)
{
  // This is a class constructor;
}

SyncSignal::~SyncSignal(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SyncSignal::Configure(){
  Sync_tree= new TTree("tree","tree");

  Sync_tree->Branch("evt",&evt);
  Sync_tree->Branch("run",&run);
  Sync_tree->Branch("lumi",&lumi);

  Sync_tree->Branch("Ptmu1",&Ptmu1);
  Sync_tree->Branch("Ptmu2",&Ptmu2);
  Sync_tree->Branch("Ptmu3",&Ptmu3);

  Sync_tree->Branch("Pmu1",&Pmu1);
  Sync_tree->Branch("Pmu2",&Pmu2);
  Sync_tree->Branch("Pmu3",&Pmu3);


  Sync_tree->Branch("Etamu1",&Etamu1);
  Sync_tree->Branch("Etamu2",&Etamu2);
  Sync_tree->Branch("Etamu3",&Etamu3);


  /*  Sync_tree->Branch("sync_eta_1",&sync_eta_1);
  Sync_tree->Branch("sync_eta_2",&sync_eta_2);
  Sync_tree->Branch("sync_eta_3",&sync_eta_3);

  Sync_tree->Branch("muon_1_isGlob",&muon_1_isGlob);
  Sync_tree->Branch("muon_2_isGlob",&muon_2_isGlob);
  Sync_tree->Branch("muon_3_isGlob",&muon_3_isGlob);

  Sync_tree->Branch("muon_1_isTrack",&muon_1_isTrack);
  Sync_tree->Branch("muon_2_isTrack",&muon_2_isTrack);
  Sync_tree->Branch("muon_3_isTrack",&muon_3_isTrack);
  */

  Sync_tree->Branch("SVx",&SVx);
  Sync_tree->Branch("SVy",&SVy);
  Sync_tree->Branch("SVz",&SVz);
  Sync_tree->Branch("fv_nC",&fv_nC);
  Sync_tree->Branch("fv_dphi3D",&fv_dphi3D);



  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2;
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
    else if(i==TauMassCut){
      title.at(i)="$\\tau$ mass (sideband in data)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");

      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,50,1.4,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,50,1.4,2.2,hlabel,"Events"));
    }


  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  NSignalCandidates=HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"N candidates","Events");
  KFChi2 =HConfig.GetTH2D(Name+"_KFChi2","KFChi2",20,0,0.1,20,0,0.1,"#Delta R (mc -first) ","#Delta R (mc - second)");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

  // Setup Extra Histograms


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  SyncSignal::Store_ExtraDist(){ 

  Extradist1d.push_back(&NSignalCandidates);
  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);

  Extradist2d.push_back(&KFChi2);

}


void  SyncSignal::doEvent(){ 



  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  value.at(TriggerOk)=0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    //    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu"))) value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);
    //    if((HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu")   or  HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu")  or HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass"))  and Ntp->HLTDecision(iTrigger)  )

    if(HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu")  && Ntp->HLTDecision(iTrigger)) value.at(TriggerOk)=true;//
    
    //DoubleMu3_Trk_Tau3mu
    
    //value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);

  }
  
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


  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();
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

    for (auto &i:Ntp-> ThreeMuons_TriggerMatch_dR(signal_idx)){
      value.at(TriggerMatch)+=i; 
    }
    
    value.at(TauMassCut) = TauLV.M();
  }
  pass.at(SignalCandidate) = (value.at(SignalCandidate) > 0);
  pass.at(Mu1PtCut) = true;//(value.at(Mu1PtCut) > cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = true;//(value.at(Mu2PtCut) > cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = true;//(value.at(Mu3PtCut) > cut.at(Mu3PtCut));
  pass.at(MuonID) =true;//(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = true;//(value.at(TriggerMatch)  <  cut.at(TriggerMatch));
  pass.at(PhiVeto) = true;//(fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*PDG_Var::Phi_width());
  pass.at(OmegaVeto) = true;//(fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*PDG_Var::Omega_width());


  pass.at(TauMassCut) = true;

  //  if(id!=1) pass.at(TauMassCut) = true;
  //  else  pass.at(TauMassCut) = ( (value.at(TauMassCut) > tauMinSideBand_ && value.at(TauMassCut) < tauMinMass_)  ||   (value.at(TauMassCut)> tauMaxMass_ && value.at(TauMassCut) < tauMaxSideBand_));



  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);


  //  if(Ntp->EventNumber()==  2132441755 or  Ntp->EventNumber()==2132442503 or  Ntp->EventNumber()==2216999905 or  Ntp->EventNumber()==2218018017 or  Ntp->EventNumber()==2218153836 or  Ntp->EventNumber()==2224367468 or  Ntp->EventNumber()==41571595 or  Ntp->EventNumber()==58447138 or  Ntp->EventNumber()==60544639 or  Ntp->EventNumber()==61679792 or  Ntp->EventNumber()==61764355 or  Ntp->EventNumber()==61953079 or  Ntp->EventNumber()==73515741 or  Ntp->EventNumber()==77540867 or  Ntp->EventNumber()==77688484 or  Ntp->EventNumber()==79478760 or  Ntp->EventNumber()==82606776 or  Ntp->EventNumber()==82622592 or  Ntp->EventNumber()==244383799 or  Ntp->EventNumber()==247314349 or  Ntp->EventNumber()==249332949 or  Ntp->EventNumber()==251563711 or  Ntp->EventNumber()==254191727 or  Ntp->EventNumber()==254932255 or  Ntp->EventNumber()==330051129 or  Ntp->EventNumber()==330446387 or  Ntp->EventNumber()==336362564 or  Ntp->EventNumber()==351321174 or  Ntp->EventNumber()==1885381355 or  Ntp->EventNumber()==1885637011 or  Ntp->EventNumber()==1908089549 or  Ntp->EventNumber()==1908162449 or  Ntp->EventNumber()==2175517428 or  Ntp->EventNumber()==2211403491 or  Ntp->EventNumber()==2213054018 or  Ntp->EventNumber()==2213935135 or  Ntp->EventNumber()==2214639779 or  Ntp->EventNumber()==2214844207 or  Ntp->EventNumber()==2215333651 or  Ntp->EventNumber()==2215692202 or  Ntp->EventNumber()==337004381 or  Ntp->EventNumber()==2836236016 or  Ntp->EventNumber()==2891478798 or  Ntp->EventNumber()==2892333356 or  Ntp->EventNumber()==2892485321 or  Ntp->EventNumber()==2905102182 or  Ntp->EventNumber()==2984100754 or  Ntp->EventNumber()==2984128642 or  Ntp->EventNumber()==2932977849 or  Ntp->EventNumber()==3666838055 or  Ntp->EventNumber()==3668909240 or  Ntp->EventNumber()==1044803226 or  Ntp->EventNumber()==1052984755 or  Ntp->EventNumber()==1065922170 or  Ntp->EventNumber()==1066066748 or  Ntp->EventNumber()==1066156777 or  Ntp->EventNumber()==1070543292 or  Ntp->EventNumber()==1091453148 or  Ntp->EventNumber()==1093002932 or  Ntp->EventNumber()==3762072241 or  Ntp->EventNumber()==3769424927 or  Ntp->EventNumber()==3774427758 or  Ntp->EventNumber()==3775170225 or  Ntp->EventNumber()==3775315149 or  Ntp->EventNumber()==3777258633 or  Ntp->EventNumber()==3795863837 or  Ntp->EventNumber()==3796093483 or  Ntp->EventNumber()==3796380512 or  Ntp->EventNumber()==3796611059 or  Ntp->EventNumber()==3805532781 or  Ntp->EventNumber()==3805600582 or  Ntp->EventNumber()==3805727172 or  Ntp->EventNumber()==3808944919 or  Ntp->EventNumber()==3812717658 or  Ntp->EventNumber()==3814520986 or  Ntp->EventNumber()==3815192506 or  Ntp->EventNumber()==3818120062 or  Ntp->EventNumber()==3818597616 or  Ntp->EventNumber()==3819472164 or  Ntp->EventNumber()==3821621232 or  Ntp->EventNumber()==3821929754 or  Ntp->EventNumber()==3821966898 or  Ntp->EventNumber()==3822260627 or  Ntp->EventNumber()==3830371384 or  Ntp->EventNumber()==3834572121 or  Ntp->EventNumber()==3837597884 or  Ntp->EventNumber()==3839158403 or  Ntp->EventNumber()==3847710495 or  Ntp->EventNumber()==1178627079 or  Ntp->EventNumber()==1247838103 or  Ntp->EventNumber()==2911636352 or  Ntp->EventNumber()==2912023250 or  Ntp->EventNumber()==2962227502 or  Ntp->EventNumber()==2965959367 or  Ntp->EventNumber()==2968060926 or  Ntp->EventNumber()==2975740897 or  Ntp->EventNumber()==2977185040 or  Ntp->EventNumber()==2978035514 or  Ntp->EventNumber()==2982321151 or  Ntp->EventNumber()==2985544485 or  Ntp->EventNumber()==1045688150 or  Ntp->EventNumber()==1092716992 or  Ntp->EventNumber()==1096339395 or  Ntp->EventNumber()==1097999561 or  Ntp->EventNumber()==1099035935 or  Ntp->EventNumber()==1100059313 or  Ntp->EventNumber()==1166692433 or  Ntp->EventNumber()==1166789026 or  Ntp->EventNumber()==1167061228 or  Ntp->EventNumber()==1168608430 or  Ntp->EventNumber()==1170668839 or  Ntp->EventNumber()==1170810737 or  Ntp->EventNumber()==1172887630){



    if(Ntp->EventNumber()==  41529489 or  Ntp->EventNumber()==82515382 or  Ntp->EventNumber()==250715270 or  Ntp->EventNumber()==255043546 or  Ntp->EventNumber()==2135865441 or  Ntp->EventNumber()==2892411351 or  Ntp->EventNumber()==2984144772 or  Ntp->EventNumber()==1065641651 or  Ntp->EventNumber()==1091471246 or  Ntp->EventNumber()==1094939126 or  Ntp->EventNumber()==3759433684 or  Ntp->EventNumber()==3774594172 or  Ntp->EventNumber()==3778240220 or  Ntp->EventNumber()==3797652855 or  Ntp->EventNumber()==3814552384 or  Ntp->EventNumber()==3819375874 or  Ntp->EventNumber()==3836387841 or  Ntp->EventNumber()==2977321169 or  Ntp->EventNumber()==1092217350 or  Ntp->EventNumber()==1092681120){



    //  if(Ntp->EventNumber()==3778090386 or  Ntp->EventNumber()==2988198475){

    std::cout<<"status (passed selection? ) "<< status << std::endl;
    std::cout<<"  EventNumber   "  << Ntp->EventNumber() <<"  N candidates  "<< Ntp-> NThreeMuons()<< " trigger:  "<< pass.at(TriggerOk) << " signal   "<< pass.at(SignalCandidate) << std::endl;
    
    std::cout<<" ----   Print Out HLT's for this event     "	<< std::endl;
    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      //      std::cout<<"   trigger  "<<   iTrigger << std::endl;
      if(Ntp->HLTDecision(iTrigger))      std::cout<<" Trigger Name:   "<<  Ntp->HLTName(iTrigger)<< "     "<< Ntp->HLTDecision(iTrigger) << std::endl;
    }

    std::cout<<" loop over signal candidates in this event: "<< std::endl;
    for(unsigned int icand=0; icand < Ntp-> NThreeMuons(); icand++){

      std::cout<<"candidate:  "<< icand<< "  Muon1_pt  "<< Ntp->Muon_P4(Ntp->ThreeMuonIndices(icand).at(0)).Pt()<< "  Muon2_pt  "<< Ntp->Muon_P4(Ntp->ThreeMuonIndices(icand).at(1)).Pt()<< " Muon3_pt   "<< Ntp->Muon_P4(Ntp->ThreeMuonIndices(icand).at(2)).Pt()<< std::endl;


    }

  }
  if(status){


    NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    TLorentzVector TauLV = Muon1LV+Muon2LV+Muon3LV;

    Ptmu1 =  Muon1LV.Pt();
    Ptmu2 =  Muon2LV.Pt();
    Ptmu3 =  Muon3LV.Pt();


    Pmu1 =  Muon1LV.P();
    Pmu2 =  Muon2LV.P();
    Pmu3 =  Muon3LV.P();

    Etamu1 =  Muon1LV.Eta();
    Etamu2 =  Muon2LV.Eta();
    Etamu3 =  Muon3LV.Eta();

    muon_1_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_1);
    muon_2_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_2);
    muon_3_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_3);

    muon_1_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_1);
    muon_2_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_2);
    muon_3_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_3);
    
    evt = Ntp->EventNumber();
    run = Ntp->RunNumber();
    lumi = Ntp->LuminosityBlock();
    //    std::cout<<"evt:  "<< evt<< std::endl;
    SVx = Ntp->Vertex_Signal_KF_pos(signal_idx).X();
    SVy = Ntp->Vertex_Signal_KF_pos(signal_idx).Y();
    SVz = Ntp->Vertex_Signal_KF_pos(signal_idx).Z();
    fv_nC = Ntp->Vertex_Signal_KF_Chi2(signal_idx);
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    fv_dphi3D=    SVPV.Angle(TauLV.Vect());

      //    if(Ntp->NThreeMuons()==1) Sync_tree->Fill();
    Sync_tree->Fill();

    /*
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){
	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;


	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);

		if(Ntp->NThreeMuons()==2){
	  int index1,index2;
	  if(signal_idx==0){index1 = 0; index2=1;}
	  if(signal_idx==1){index1 = 1; index2=0;}

	  TLorentzVector Tau1 = Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index1)).at(0)))+
	    Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index1)).at(1)))+
	    Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index1)).at(2)));
	  TLorentzVector Tau2 = Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index2)).at(0)))+
	    Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index2)).at(1)))+
	    Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(index2)).at(2)));
	  

	  KFChi2.at(t).Fill(Tau1.DeltaR(MCTauLV),Tau2.DeltaR(MCTauLV));
	  
	}
	

      }
    }
    */
    
  }
}


void  SyncSignal::Finish(){

  file= new TFile("Sync_signal_tree_UF.root","recreate");
  Sync_tree->SetDirectory(file);

  file->Write();
  file->Close();

  Selection::Finish();
}





