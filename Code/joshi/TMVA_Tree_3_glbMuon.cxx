#include "TMVA_Tree_3_glbMuon.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

TMVA_Tree_3_glbMuon::TMVA_Tree_3_glbMuon(TString Name_, TString id_):
    Selection(Name_,id_),
    tauMinMass_(1.731),
    tauMaxMass_(1.823),
    tauMinSideBand_(1.6),
    tauMaxSideBand_(2.0),
	 massRes_(0.02)
{
    // This is a class constructor;
}


TMVA_Tree_3_glbMuon::~TMVA_Tree_3_glbMuon(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }
    Logger(Logger::Info) << "complete." << std::endl;
}

void  TMVA_Tree_3_glbMuon::Configure(){
    for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==TriggerOk) cut.at(TriggerOk)=1;
      if(i==SignalCandidate) cut.at(SignalCandidate)=1;
      if(i==Mu1PtCut)  cut.at(Mu1PtCut)=2.5;
      if(i==Mu2PtCut)  cut.at(Mu2PtCut)=2.5;
      if(i==Mu3PtCut)  cut.at(Mu3PtCut)=2.5;
      if(i==MuonID) cut.at(MuonID)=1;
      if(i==OmegaVeto) cut.at(OmegaVeto)=0; // defined below
      if(i==TriggerMatch) cut.at(TriggerMatch)=0.03;
      if(i==ThreeMuMass)  cut.at(ThreeMuMass)=1;// true for MC and mass side band for data
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
        title.at(i)="is signal candidate";
        hlabel="is 3mu candidate";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==Mu1PtCut){
        title.at(i)="Mu1 Pt";
        hlabel="Muon1 PT, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
        title.at(i)="Mu2 Pt";
        hlabel="Muon2 PT, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      }
      else if(i==Mu3PtCut){
        title.at(i)="Mu3 Pt";
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
      else if(i==TriggerMatch){
        title.at(i)="Trigger Matching";
        hlabel="Sum of dR_{reco-trigger}";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      }
      else if(i==ThreeMuMass){
        title.at(i)="Tau Mass";
        hlabel="three mu mass, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
      }
    } 

    Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
    TMVA_Tree= new TTree("tree","tree");

	/*
    //----------------------------------------------
    //   2016 TMVA variables
    //----------------------------------------------
    TMVA_Tree->Branch("var_muMinPt",&var_muMinPt);
    TMVA_Tree->Branch("var_inOutTrackMatch_Chi2", &var_inOutTrackMatch_Chi2);
    TMVA_Tree->Branch("var_mu3kink", &var_mu3kink);
    TMVA_Tree->Branch("var_vertexKFChi2", &var_vertexKFChi2);
    TMVA_Tree->Branch("var_VertexMu3D0Sig", &var_VertexMu3D0Sig);
    TMVA_Tree->Branch("var_VertexMu3D0", &var_VertexMu3D0);
    TMVA_Tree->Branch("var_mindca_iso", &var_mindca_iso);
    TMVA_Tree->Branch("var_iso_relpt", &var_iso_relpt);
    TMVA_Tree->Branch("var_fv_cosdphi3d", &var_fv_cosdphi3d);
	 */
  
  TMVA_Tree->Branch("var_KFV_chiSq", &var_KFV_chiSq);
  TMVA_Tree->Branch("var_flightLenSig", &var_flightLenSig);
  TMVA_Tree->Branch("var_muMinPt",&var_muMinPt);
  TMVA_Tree->Branch("var_inOutTrackMatch_Chi2", &var_inOutTrackMatch_Chi2);
  TMVA_Tree->Branch("var_relPt", &var_relPt);
  TMVA_Tree->Branch("var_mindca_iso", &var_mindca_iso);
  TMVA_Tree->Branch("var_VertexMu3D0", &var_VertexMu3D0);
  TMVA_Tree->Branch("var_osMmumu1", &var_osMmumu1);
  TMVA_Tree->Branch("var_osMmumu2", &var_osMmumu2);
  
    Selection::ConfigureHistograms(); //do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove

}




void  TMVA_Tree_3_glbMuon::Store_ExtraDist(){ 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
    ////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  TMVA_Tree_3_glbMuon::doEvent(){ 
    unsigned int t;
    int id(Ntp->GetMCID());
    if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu"))) value.at(TriggerOk)=Ntp->HLTDecision(iTrigger);
    }

    pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));
    unsigned int final_idx = 0;

    value.at(TriggerMatch)=0;
    value.at(SignalCandidate)=0;
    value.at(Mu1PtCut)=0;
    value.at(Mu2PtCut)=0;
    value.at(Mu3PtCut)=0;
    value.at(TriggerMatch)=0;
    value.at(MuonID)=0;
    value.at(ThreeMuMass)=0;

    if(Ntp->NThreeMuons()>0){
      value.at(SignalCandidate) = Ntp->NThreeMuons();
      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
      // value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
      //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
      unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
      //
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
            Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
            Ntp->Muon_isGlobalMuon(mu3_pt_idx));
      //------------------------------------------------------------------------------------------------------

      value.at(Mu1PtCut) = Ntp->Muon_P4(mu1_pt_idx).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(mu2_pt_idx).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(mu3_pt_idx).Pt();


      TLorentzVector TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);



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

    if(id!=1) pass.at(ThreeMuMass) = true;
    else  pass.at(ThreeMuMass) = ( (value.at(ThreeMuMass) > tauMinSideBand_ && value.at(ThreeMuMass) < tauMinMass_)  || (value.at(ThreeMuMass)> tauMaxMass_ && value.at(ThreeMuMass) < tauMaxSideBand_));


    double wobs=1;
    double w; 

    if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
    else{w=1;}

    bool status=AnalysisCuts(t,w,wobs);
    
	 if(status){

      unsigned int Muon_index_1 =  Ntp->ThreeMuonIndices(final_idx).at(0);
      unsigned int Muon_index_2 =  Ntp->ThreeMuonIndices(final_idx).at(1);
      unsigned int Muon_index_3 =  Ntp->ThreeMuonIndices(final_idx).at(2);

      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
      TLorentzVector TauLV = Muon1LV+Muon2LV+Muon3LV;

      /*//------------------- 2016 variables ----------------------------
      var_muMinPt = Ntp->Muon_P4(Muon_index_3).Pt();
      var_inOutTrackMatch_Chi2 = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3);
      var_mu3kink = Ntp->Muon_combinedQuality_trkKink(Muon_index_3);
      var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx);
      var_VertexMu3D0Sig = Ntp->Vertex_d0sig_reco(final_idx,2);
      var_VertexMu3D0 = Ntp->Vertex_d0_reco(final_idx,2);
      var_mindca_iso = Ntp->Isolation_MinDist(final_idx);
      var_iso_relpt = Ntp->Isolation_RelPt(final_idx);
      //------------------ calculate dv_cosdphi3d ---------------------
      TVector3 vec_sv = Ntp->Vertex_Signal_KF_pos(final_idx);
	   TVector3 vec_pv(0,0,0);
      if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
         vec_pv = Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx);
	      }
      TVector3 vec_tau = TauLV.Vect();
      TVector3 d_pv_sv = vec_sv - vec_pv;
      float fv_cosdphi3d = (d_pv_sv.Dot(vec_tau)/(d_pv_sv.Mag()*vec_tau.Mag()));
      //---------------------------------------------------------------
      var_fv_cosdphi3d = fv_cosdphi3d;
      //---------------------------------------------------------------*/

	 TVector3 fls_pv = Ntp->Vertex_MatchedPrimaryVertex(final_idx);
	 TVector3 fls_sv = Ntp->Vertex_Signal_KF_pos(final_idx);
	 TMatrixTSym<double> fls_PVcov = Ntp->Vertex_PrimaryVertex_Covariance(final_idx);
	 TMatrixTSym<double> fls_SVcov = Ntp->Vertex_Signal_KF_Covariance(final_idx);
	
   vector<unsigned int> idx_vec;

    idx_vec.push_back(Muon_index_1);
    idx_vec.push_back(Muon_index_2);
    idx_vec.push_back(Muon_index_3);

    unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

  
var_KFV_chiSq = Ntp->Vertex_signal_KF_Chi2(final_idx);
var_flightLenSig = Ntp->FlightLength_significance(fls_pv,fls_PVcov,fls_sv,fls_SVcov); // Add flight length significance
var_mindca_iso = Ntp->Isolation_MinDist(final_idx);
var_muMinPt = Ntp->Muon_P4(Muon_index_3).Pt();
var_inOutTrackMatch_Chi2 = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3);
var_VertexMu3D0 = Ntp->Vertex_d0_reco(final_idx,2);
var_relPt = Ntp->Isolation_RelPt(final_idx);
var_osMmumu1 = M_osss1;
var_osMmumu2 = M_osss2;
  
  if (Ntp->Vertex_RefitPVisValid(final_idx)==1)  TMVA_Tree->Fill();
    }

}



void  TMVA_Tree_3_glbMuon::Finish(){
file= new TFile("TMVA_Tree_3_glbMuonInput.root","recreate");
    TMVA_Tree->SetDirectory(file);

    file->Write();
    file->Close();

    Selection::Finish();

}
