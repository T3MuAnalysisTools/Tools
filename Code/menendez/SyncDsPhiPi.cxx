#include "SyncDsPhiPi.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"



SyncDsPhiPi::SyncDsPhiPi(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

SyncDsPhiPi::~SyncDsPhiPi(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  SyncDsPhiPi::Configure(){

  Sync_tree= new TTree("tree","tree");
  Sync_tree->Branch("sync_pt_1",&sync_pt_1);
  Sync_tree->Branch("sync_pt_2",&sync_pt_2);
  Sync_tree->Branch("sync_pt_3",&sync_pt_3);

  Sync_tree->Branch("sync_eta_1",&sync_eta_1);
  Sync_tree->Branch("sync_eta_2",&sync_eta_2);
  Sync_tree->Branch("sync_eta_3",&sync_eta_3);

  Sync_tree->Branch("phi_mass",&phi_mass);
  Sync_tree->Branch("ds_mass",&ds_mass);

  Sync_tree->Branch("evt",&evt);
  Sync_tree->Branch("run",&run);
  Sync_tree->Branch("lumi",&lumi);

  Sync_tree->Branch("sync_DsPhiPiVtx_x",&sync_DsPhiPiVtx_x);
  Sync_tree->Branch("sync_DsPhiPiVtx_y",&sync_DsPhiPiVtx_y);
  Sync_tree->Branch("sync_DsPhiPiVtx_z",&sync_DsPhiPiVtx_z);
  Sync_tree->Branch("sync_DsPhiPiVtx_Chi2",&sync_DsPhiPiVtx_Chi2);




  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==HLTOk)           cut.at(HLTOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==PhiMassCut)      cut.at(PhiMassCut)=1;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }	
    else if(i==PhiMassCut){
      title.at(i)= "0.95 $<$ $M_{\\mu\\mu}$ $<$ 1.1 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="OS muon pair mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiMassCut_",htitle,40,0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiMassCut_",htitle,40,0.5,1.5,hlabel,"Events"));
    }
    
  }
  // Track Candidate Information

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  SyncDsPhiPi::Store_ExtraDist(){ 
  
	 
}


void  SyncDsPhiPi::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
 
  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  int mu1=-1, mu2=-1, track=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999.0;
  
  value.at(is2MuTrk) = 0;
  value.at(PhiMassCut) = 0;

  if(Ntp->NTwoMuonsTrack()!=0 /*&& Ntp->NThreeMuons() == 0*/) value.at(is2MuTrk) = 1;
 
  if (value.at(is2MuTrk)==1){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
 
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);
      if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M)){
	tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
	mu1 = tmp_mu1;
	mu2 = tmp_mu2;
	track = tmp_track;
	tmp_idx = i2M;
      }
    }
  value.at(PhiMassCut) =  (Ntp->Muon_P4(mu1) + Ntp->Muon_P4(mu2)).M();
  }
  

  pass.at(PhiMassCut) = ( value.at(PhiMassCut) > 0.95 && value.at(PhiMassCut) < 1.1);
  pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  double wobs=1;
  double w;  
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);

  std::vector<unsigned int> misses{80586988,98764233,48390993,25573281,15034108,13262859,12};

  for (unsigned int i=0;i<misses.size();i++) {
    if (Ntp->EventNumber() == misses[i]) {
      std::cout << "///////////////////////////////////////////////////////////" << std::endl;
      std::cout << "For Event Number " << Ntp->EventNumber() << " , passes are: PhiMassCut=" << pass.at(PhiMassCut) << ", is2MuTrk=" << pass.at(is2MuTrk) << ", HLTOk=" << pass.at(HLTOk) << ". PhiMassCut value = " << value.at(PhiMassCut) << std::endl;
    }
  }

  if(status){


    TLorentzVector Mu1LV;
    TLorentzVector Mu2LV;
    TLorentzVector TrackLV = Ntp->Track_P4(track);

    if(Ntp->Muon_P4(mu1).Pt() > Ntp->Muon_P4(mu2).Pt()){
      Mu1LV = Ntp->Muon_P4(mu1);
      Mu2LV = Ntp->Muon_P4(mu2);

    }
    else {

      Mu1LV = Ntp->Muon_P4(mu2);
      Mu2LV = Ntp->Muon_P4(mu1);
    }

    sync_pt_1 = Mu1LV.Pt();
    sync_pt_2 = Mu2LV.Pt();
    sync_pt_3 = TrackLV.Pt();

    sync_eta_1 = Mu1LV.Eta();
    sync_eta_2 = Mu2LV.Eta();
    sync_eta_3 = TrackLV.Eta();

    phi_mass  = (Mu1LV+Mu2LV).M();
    ds_mass  = (Mu1LV+Mu2LV  + TrackLV).M();

    evt = Ntp->EventNumber();
    run = Ntp->RunNumber();
    lumi = Ntp->LuminosityBlock();

//    sync_DsPhiPiVtx_x = Ntp->Vertex_Signal_KF_pos(tmp_idx).X();
//    sync_DsPhiPiVtx_y = Ntp->Vertex_Signal_KF_pos(tmp_idx).Y();
//    sync_DsPhiPiVtx_z = Ntp->Vertex_Signal_KF_pos(tmp_idx).Z();
    sync_DsPhiPiVtx_Chi2 = Ntp->TwoMuonsTrack_SV_Chi2(tmp_idx);



    Sync_tree->Fill();

  }
}

void  SyncDsPhiPi::Finish(){


  file= new TFile("Sync_dsphipi_tree_UF.root","recreate");
  Sync_tree->SetDirectory(file);

  file->Write();
  file->Close();

  Selection::Finish();

}



