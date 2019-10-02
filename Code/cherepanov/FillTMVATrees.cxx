#include "FillTMVATrees.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

FillTMVATrees::FillTMVATrees(TString Name_, TString id_):
    Selection(Name_,id_),
    tauMinMass_(1.731),
    tauMaxMass_(1.823),
    tauMinSideBand_(1.6),
    tauMaxSideBand_(2.0),
    tauMassResCutLow(0.007),
    tauMassResCutHigh(0.01)
{
    // This is a class constructor;
}


FillTMVATrees::~FillTMVATrees(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }
    Logger(Logger::Info) << "complete." << std::endl;
}

void  FillTMVATrees::Configure(){
  // Set tree branches
  TMVA_Tree= new TTree("tree","tree");
  TMVA_Tree->Branch("MC",&MC);
  TMVA_Tree->Branch("category",&category);
  TMVA_Tree->Branch("var_vertexKFChi2",&var_vertexKFChi2);
  TMVA_Tree->Branch("var_svpvTauAngle",&var_svpvTauAngle);
  TMVA_Tree->Branch("var_flightLenSig",&var_flightLenSig);
  TMVA_Tree->Branch("var_sumMuTrkKinkChi2",&var_sumMuTrkKinkChi2);
  TMVA_Tree->Branch("var_segCompMuMin",&var_segCompMuMin);
  TMVA_Tree->Branch("var_segCompMuMax",&var_segCompMuMax);
  TMVA_Tree->Branch("var_segCompMu1",&var_segCompMu1);
  TMVA_Tree->Branch("var_segCompMu2",&var_segCompMu2);
  TMVA_Tree->Branch("var_segCompMu3",&var_segCompMu3);
  TMVA_Tree->Branch("var_caloCompMin",&var_caloCompMin);
  TMVA_Tree->Branch("var_caloCompMax",&var_caloCompMax);
  TMVA_Tree->Branch("var_caloCompMu1",&var_caloCompMu1);
  TMVA_Tree->Branch("var_caloCompMu2",&var_caloCompMu2);
  TMVA_Tree->Branch("var_caloCompMu3",&var_caloCompMu3);
  TMVA_Tree->Branch("var_MinMIPLikelihood",&var_MinMIPLikelihood);
  TMVA_Tree->Branch("var_tauMass",&var_tauMass);
  TMVA_Tree->Branch("var_ntracks",&var_ntracks);
  TMVA_Tree->Branch("var_relPt",&var_relPt);
  TMVA_Tree->Branch("var_isoMax",&var_isoMax);
  TMVA_Tree->Branch("var_mu3d0VertexSig",&var_mu3d0VertexSig);
  TMVA_Tree->Branch("var_maxdca",&var_maxdca);
  
  // -----------------
  
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk) cut.at(TriggerOk)=1;
    if(i==SignalCandidate) cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)  cut.at(Mu1PtCut)=2.0;
    if(i==Mu2PtCut)  cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)  cut.at(Mu3PtCut)=1.0;
    if(i==MuonID) cut.at(MuonID)=1;
    if(i==PhiVeto) cut.at(PhiVeto)=0; // defined below
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
      title.at(i)="Mu1 Pt $>$ ";
      title.at(i)+=cut.at(Mu1PtCut);
      title.at(i)+=" GeV";
      hlabel="Muon1 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
    }
    else if(i==Mu2PtCut){
      title.at(i)="Mu2 Pt $>$ ";
      title.at(i)+=cut.at(Mu2PtCut);
      title.at(i)+=" GeV";
      hlabel="Muon2 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
    }
    else if(i==Mu3PtCut){
      title.at(i)="Mu3 Pt $>$ ";
      title.at(i)+=cut.at(Mu3PtCut);
      title.at(i)+=" GeV";
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
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto){
      title.at(i)="omega mass veto";
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
      title.at(i)="Tau Mass";
      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
    }
  } 
  
      Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
      VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
      MuonglbkinkSum  =HConfig.GetTH1D(Name+"_MuonglbkinkSum","MuonglbkinkSum",50,0.,50," #sum  #mu glb kink #chi^{2}","Events");
      FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");
      SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");
      Muon_segmentCompatibility_mu1  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu1","Muon_segmentCompatibility_mu1",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
      Muon_segmentCompatibility_mu2  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu2","Muon_segmentCompatibility_mu2",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
      Muon_segmentCompatibility_mu3  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu3","Muon_segmentCompatibility_mu3",50,0.,1,"Inner Track and muon segment match  #mu_{3} ","Events");
      
      Muon_segmentCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_min","Muon_segmentCompatibility_min",50,0.,1,"Inner Track and muon segment match min ","Events");
      Muon_segmentCompatibility_max  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_max","Muon_segmentCompatibility_max",50,0.,1,"Inner Track and muon segment match max ","Events");
      
      Muon_ECALCompatibility_mu1  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu1","Muon_ECALCompatibility_mu1",50,0.,1,"MIP Likelihood  #mu_{1} ","Events");
      Muon_ECALCompatibility_mu2  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu2","Muon_ECALCompatibility_mu2",50,0.,1,"MIP Likelihood  #mu_{2} ","Events");
      Muon_ECALCompatibility_mu3  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu3","Muon_ECALCompatibility_mu3",50,0.,1,"MIP Likelihood  #mu_{3} ","Events");
      
      Muon_ECALCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_min","Muon_ECALCompatibility_min",50,0.,1,"MIP Likelihood min ","Events");
      Muon_ECALCompatibility_max  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_max","Muon_ECALCompatibility_max",50,0.,1,"MIP Likelihood max ","Events");
      
      Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","Isolation_NTracks",10,-0.5,9.5,"N tracks","Events");
      Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","Isolation_RelPt",50,0,1,"relative p_{T}","Events");
      Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","Isolation_maxdxy",40,0,15,"Iso maximum transversely displaced track","Events");
      VertexMu3D0SigReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigReco","VertexMu3D0SigReco",50,0,3,"#mu_{3} - PV transverse distance significance","Events");
      
      VertexDCAMax=HConfig.GetTH1D(Name+"_VertexDCAMax","VertexDCAMax",40,0,0.15,"Max closest distance between muons","Events");
      
      Selection::ConfigureHistograms(); //do not remove
      HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
      
}

void  FillTMVATrees::Store_ExtraDist(){ 
  
  Extradist1d.push_back(&SVPVTauDirAngle);
  Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&MuonglbkinkSum);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu1);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu2);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu3);
  Extradist1d.push_back(&Muon_segmentCompatibility_min);
  Extradist1d.push_back(&Muon_segmentCompatibility_max);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu1);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu2);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu3);
  Extradist1d.push_back(&Muon_ECALCompatibility_min);
  Extradist1d.push_back(&Muon_ECALCompatibility_max);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_NTracks);
  Extradist1d.push_back(&VertexMu3D0SigReco);
  Extradist1d.push_back(&Isolation_maxdxy);
  Extradist1d.push_back(&VertexDCAMax);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  ////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  FillTMVATrees::doEvent(){ 
    unsigned int t;
    int id(Ntp->GetMCID());
    bool hlt_pass = false;

    if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if( HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger)==1) hlt_pass = true;
      if( HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") && Ntp->HLTDecision(iTrigger)==1) hlt_pass = true;
    }
    
    value.at(TriggerOk) = hlt_pass;
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
			  (Ntp->Muon_isGlobalMuon(mu3_pt_idx) || Ntp->Muon_isTrackerMuon(mu3_pt_idx)));
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

      value.at(PhiVeto) = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
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
    pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*PDG_Var::Phi_width());
    pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*PDG_Var::Omega_width());

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

      std::vector<unsigned int> indices;
      indices.push_back(Muon_index_1);
      indices.push_back(Muon_index_2);
      indices.push_back(Muon_index_3);
      
      float tauMassRes = Ntp->TauMassResolution(indices,1,1);
      
      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
      TLorentzVector TauLV = Muon1LV+Muon2LV+Muon3LV;
      
      float MinSegmentCompatibility = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      float MaxSegmentCompatibility = std::max({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      float MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      float MaxMIPLikelihood = std::max({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      //------------------ calculate var_svpvTauAngle ---------------------
      TVector3 vec_sv = Ntp->Vertex_Signal_KF_pos(final_idx);
      TVector3 vec_pv(0,0,0);
      if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
	vec_pv = Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx);
      }
      TVector3 vec_tau = TauLV.Vect();
      TVector3 d_pvsv = vec_sv - vec_pv;
      //---------------------------------------------------------------
      //------------------ calculate var_flightLenSig ---------------------
      TMatrixTSym<double> fls_PVcov = Ntp->Vertex_PrimaryVertex_Covariance(final_idx);
      TMatrixTSym<double> fls_SVcov = Ntp->Vertex_Signal_KF_Covariance(final_idx);
      
      var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx);
      var_svpvTauAngle = TMath::ACos(d_pvsv.Dot(vec_tau)/(d_pvsv.Mag()*vec_tau.Mag()));
      var_flightLenSig = sqrt(Ntp->FlightLength_significance(vec_pv,fls_PVcov,vec_sv,fls_SVcov)); // Add flight length significance
      var_sumMuTrkKinkChi2 = (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));
      var_segCompMuMin = MinSegmentCompatibility;
      var_segCompMuMax = MaxSegmentCompatibility;
      var_segCompMu1 = Ntp->Muon_segmentCompatibility(Muon_index_1);
      var_segCompMu2 = Ntp->Muon_segmentCompatibility(Muon_index_2);
      var_segCompMu3 = Ntp->Muon_segmentCompatibility(Muon_index_3);
      var_caloCompMin = MinMIPLikelihood;
      var_caloCompMax = MaxMIPLikelihood;
      var_caloCompMu1 = Ntp->Muon_caloCompatibility(Muon_index_1);
      var_caloCompMu2 = Ntp->Muon_caloCompatibility(Muon_index_2);
      var_caloCompMu3 = Ntp->Muon_caloCompatibility(Muon_index_3);
      var_MinMIPLikelihood = MinMIPLikelihood;
      var_tauMass = TauLV.M();
      var_ntracks = Ntp->Isolation_NTracks(final_idx);
      var_relPt = Ntp->Isolation_RelPt(final_idx);
      var_isoMax = Ntp->Isolation_maxdy(final_idx);
      var_mu3d0VertexSig = Ntp->Vertex_d0sig_reco(final_idx,2);
      var_maxdca = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
      
      if (id==1) MC=0;
      else  MC=1;
	 
      if (id==1){
	if (tauMassRes<tauMassResCutLow && value.at(ThreeMuMass)<PDG_Var::Tau_mass()) category = 1;
	if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh && value.at(ThreeMuMass)<PDG_Var::Tau_mass()) category = 2;
	if (tauMassRes>tauMassResCutHigh && value.at(ThreeMuMass)<PDG_Var::Tau_mass()) category = 3;
	if (tauMassRes<tauMassResCutLow && value.at(ThreeMuMass)>PDG_Var::Tau_mass()) category = 4;
	if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh && value.at(ThreeMuMass)>PDG_Var::Tau_mass()) category = 5;
	if (tauMassRes>tauMassResCutHigh && value.at(ThreeMuMass)>PDG_Var::Tau_mass()) category = 6;
      }
      else category = 0;
      if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
	
	TMVA_Tree->Fill();
	
	// ----- Fill the histograms -----
	VertexDCAMax.at(t).Fill(var_maxdca,w);
	SVPVTauDirAngle.at(t).Fill(var_svpvTauAngle);
	FLSignificance.at(t).Fill(var_flightLenSig);
	VertexChi2KF.at(t).Fill(var_vertexKFChi2);
	MuonglbkinkSum.at(t).Fill(var_sumMuTrkKinkChi2);
	Muon_segmentCompatibility_mu1.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_1),w);
	Muon_segmentCompatibility_mu2.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_2),w);
	Muon_segmentCompatibility_mu3.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_3),w);
	Muon_segmentCompatibility_min.at(t).Fill(MinSegmentCompatibility,w);
	Muon_segmentCompatibility_max.at(t).Fill(MaxSegmentCompatibility,w);
	Muon_ECALCompatibility_mu1.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_1),w);
	Muon_ECALCompatibility_mu2.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_2),w);
	Muon_ECALCompatibility_mu3.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_3),w);
	Muon_ECALCompatibility_min.at(t).Fill(MinMIPLikelihood,w);
	Muon_ECALCompatibility_max.at(t).Fill(MaxMIPLikelihood,w);
	Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(final_idx),w);
	Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(final_idx),w);
	Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(final_idx),w);
	VertexMu3D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,2),w);
      }
    }
}



void  FillTMVATrees::Finish(){
  
  if(mode == RECONSTRUCT){
    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    int id(Ntp->GetMCID());
    double scale(1.);
    double scaleDsTau(0.637);
    double scaleBpTau(0.262);
    double scaleB0Tau(0.099);

    if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    ScaleAllHistOfType(2,scale*scaleDsTau);
    
    if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    ScaleAllHistOfType(3,scale*scaleB0Tau);

    if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    ScaleAllHistOfType(4,scale*scaleBpTau);

    //    }
  }
    file= new TFile("TMVATreesInput.root","recreate");
    TMVA_Tree->SetDirectory(file);
    
    file->Write();
    file->Close();
    
    Selection::Finish();
}
