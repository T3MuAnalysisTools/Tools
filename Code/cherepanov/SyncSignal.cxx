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

  Sync_tree->Branch("sync_pt_1",&sync_pt_1);
  Sync_tree->Branch("sync_pt_2",&sync_pt_2);
  Sync_tree->Branch("sync_pt_3",&sync_pt_3);

  Sync_tree->Branch("sync_eta_1",&sync_eta_1);
  Sync_tree->Branch("sync_eta_2",&sync_eta_2);
  Sync_tree->Branch("sync_eta_3",&sync_eta_3);

  Sync_tree->Branch("muon_1_isGlob",&muon_1_isGlob);
  Sync_tree->Branch("muon_2_isGlob",&muon_2_isGlob);
  Sync_tree->Branch("muon_3_isGlob",&muon_3_isGlob);

  Sync_tree->Branch("muon_1_isTrack",&muon_1_isTrack);
  Sync_tree->Branch("muon_2_isTrack",&muon_2_isTrack);
  Sync_tree->Branch("muon_3_isTrack",&muon_3_isTrack);


  Sync_tree->Branch("sync_ThreeMuVtx_x",&sync_ThreeMuVtx_x);
  Sync_tree->Branch("sync_ThreeMuVtx_y",&sync_ThreeMuVtx_y);
  Sync_tree->Branch("sync_ThreeMuVtx_z",&sync_ThreeMuVtx_z);
  Sync_tree->Branch("sync_ThreeMuVtx_Chi2",&sync_ThreeMuVtx_Chi2);




  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=2.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=1.5;
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
  // Setup Extra Histograms


  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  SyncSignal::Store_ExtraDist(){ 




}


void  SyncSignal::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    //    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu"))) value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);
    if(HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu_v"))value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);

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
  if(Ntp->EventNumber()==6510923){
    std::cout<<"status "<< status << std::endl;
    std::cout<<"  N candidates  "<< Ntp-> NThreeMuons()<< " trigger:  "<< pass.at(TriggerOk) << " signal   "<< pass.at(SignalCandidate) << std::endl;

  }
  if(status){

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

    sync_pt_1 =  Muon1LV.Pt();
    sync_pt_2 =  Muon2LV.Pt();
    sync_pt_3 =  Muon3LV.Pt();

    sync_eta_1 =  Muon1LV.Eta();
    sync_eta_2 =  Muon2LV.Eta();
    sync_eta_3 =  Muon3LV.Eta();

    muon_1_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_1);
    muon_2_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_2);
    muon_3_isGlob = Ntp->Muon_isGlobalMuon(Muon_index_3);

    muon_1_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_1);
    muon_2_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_2);
    muon_3_isTrack = Ntp->Muon_isTrackerMuon(Muon_index_3);
    
    evt = Ntp->EventNumber();
    run = Ntp->RunNumber();
    lumi = Ntp->LuminosityBlock();
    std::cout<<"evt:  "<< evt<< std::endl;
    sync_ThreeMuVtx_x = Ntp->Vertex_Signal_KF_pos(signal_idx).X();
    sync_ThreeMuVtx_y = Ntp->Vertex_Signal_KF_pos(signal_idx).Y();
    sync_ThreeMuVtx_z = Ntp->Vertex_Signal_KF_pos(signal_idx).Z();
    sync_ThreeMuVtx_Chi2 = Ntp->Vertex_Signal_KF_Chi2(signal_idx);



    Sync_tree->Fill();
    
  }
}


void  SyncSignal::Finish(){

  file= new TFile("Sync_signal_tree_UF.root","recreate");
  Sync_tree->SetDirectory(file);

  file->Write();
  file->Close();

  Selection::Finish();
}





