#include "FillMVATree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

FillMVATree::FillMVATree(TString Name_, TString id_):
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


FillMVATree::~FillMVATree(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }
    Logger(Logger::Info) << "complete." << std::endl;
}

void  FillMVATree::Configure(){
    // Set tree branches
	 TMVA_Tree= new TTree("tree","tree");
	 TMVA_Tree->Branch("MC",&MC);
	 TMVA_Tree->Branch("category",&category);
	 TMVA_Tree->Branch("var_vertexKFChi2",&var_vertexKFChi2);
    TMVA_Tree->Branch("var_svpvTauAngle",&var_svpvTauAngle);
  	 TMVA_Tree->Branch("var_flightLenSig",&var_flightLenSig);
	 TMVA_Tree->Branch("var_sumMuTrkKinkChi2",&var_sumMuTrkKinkChi2);
	 TMVA_Tree->Branch("var_segCompMuMin",&var_segCompMuMin);
	 TMVA_Tree->Branch("var_MinMIPLikelihood",&var_MinMIPLikelihood);
	 TMVA_Tree->Branch("var_tauMass",&var_tauMass);
	 TMVA_Tree->Branch("var_maxdca",&var_maxdca);
    TMVA_Tree->Branch("var_MuMu_mindR",&var_MuMu_mindR);
    TMVA_Tree->Branch("var_RelPt_Mu1Tau",&var_RelPt_Mu1Tau);
    TMVA_Tree->Branch("var_Eta_au",&var_Eta_Tau);
    TMVA_Tree->Branch("var_MuMu_minKFChi2",&var_MuMu_minKFChi2);
    TMVA_Tree->Branch("var_MuTau_maxdR",&var_MuTau_maxdR);
	 TMVA_Tree->Branch("var_MaxD0Significance", &var_MaxD0Significance);
	 TMVA_Tree->Branch("var_IsolationMinDist", &var_IsolationMinDist);
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

  Muon_segmentCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_min","Muon_segmentCompatibility_min",50,0.,1,"Inner Track and muon segment match min ","Events");
  Muon_HCALCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_HCALCompatibility_min","Muon_ECALCompatibility_min",50,0.,1,"MIP Likelihood min ","Events");
  //New variables
  minMudR = HConfig.GetTH1D(Name+"_minMudR","minMudR",50,0,0.6,"min #DeltaR(#mu#mu)","Events");
  Mu1TauPTRatio = HConfig.GetTH1D(Name+"_Mu1TauPTRatio","Mu1TauPTRatio",50,0,1,"p_{T}(#mu_{1})/p_{T}(#tau)","Events");
  dRMaxMuTau = HConfig.GetTH1D(Name+"_dRMaxMuTau","dRMaxMuTau",50,0,0.5,"max #DeltaR(#mu#tau)","Events");
  MuPair_vertex_chi2_min=HConfig.GetTH1D(Name+"_MuPair_vertex_chi2_min","MuPair_vertex_chi2_min",50,0,1.5,"KF min #chi^{2} of #mu pair","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
  VertexDCAMax=HConfig.GetTH1D(Name+"_VertexDCAMax","VertexDCAMax",40,0,0.15,"Max closest distance between muons","Events");
  Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","Isolation_MinDist",50,0,0.1,"Iso MinDist","Events"); 
  VertexMuMaxD0SigReco=HConfig.GetTH1D(Name+"_VertexMuMaxD0SigReco","VertexMuMaxD0SigReco",50,0,4,"#mu - PV max transverse distance significance","Events");
  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove

}

void  FillMVATree::Store_ExtraDist(){ 
  
  Extradist1d.push_back(&SVPVTauDirAngle);
  Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&MuonglbkinkSum);
  Extradist1d.push_back(&Muon_segmentCompatibility_min);
  Extradist1d.push_back(&Muon_HCALCompatibility_min);
  
  Extradist1d.push_back(&minMudR);
  Extradist1d.push_back(&Mu1TauPTRatio);
  Extradist1d.push_back(&dRMaxMuTau);
  Extradist1d.push_back(&MuPair_vertex_chi2_min);
  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&VertexDCAMax);
  Extradist1d.push_back(&Isolation_MinDist);
  Extradist1d.push_back(&VertexMuMaxD0SigReco);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  ////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  FillMVATree::doEvent(){ 
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
    //pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*PDG_Var::Phi_width());
    //pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*PDG_Var::Omega_width());
    pass.at(PhiVeto) = (fabs(value.at(PhiVeto)-PDG_Var::Phi_mass()) > 2*0.1);
    pass.at(OmegaVeto) = (fabs(value.at(OmegaVeto)-PDG_Var::Omega_mass())> 2*0.1);

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
      float MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
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

		float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0),
			 	  										Ntp->Vertex_d0sig_reco(final_idx,1),
				  	  									Ntp->Vertex_d0sig_reco(final_idx,2)});

    var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx);
	 var_MaxD0Significance = MaxD0Significance;
	 var_IsolationMinDist = Ntp->Isolation_MinDist(final_idx);
    var_svpvTauAngle = TMath::ACos(d_pvsv.Dot(vec_tau)/(d_pvsv.Mag()*vec_tau.Mag()));
  	 var_flightLenSig = sqrt(Ntp->FlightLength_significance(vec_pv,fls_PVcov,vec_sv,fls_SVcov)); // Add flight length significance
	 var_sumMuTrkKinkChi2 = (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));
	 var_segCompMuMin = MinSegmentCompatibility;
	 var_MinMIPLikelihood = MinMIPLikelihood;
	 var_tauMass = TauLV.M();
	 var_MuMu_mindR = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
	 var_RelPt_Mu1Tau = Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt();
	 var_Eta_Tau = TauLV.Eta();
	 var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(final_idx,0), Ntp->Vertex_pair_quality(final_idx,1), Ntp->Vertex_pair_quality(final_idx,2)});
 	 var_maxdca = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
 	 var_MuTau_maxdR = std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});

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
	 else {
		 	if (tauMassRes<tauMassResCutLow) category = 1;
		 	if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh) category = 2;
		 	if (tauMassRes>tauMassResCutHigh) category = 3;
			}
	 if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
	 	
		TMVA_Tree->Fill();

		 // ----- Fill the histograms -----
    	 VertexDCAMax.at(t).Fill(var_maxdca,w);
		 SVPVTauDirAngle.at(t).Fill(var_svpvTauAngle);
       FLSignificance.at(t).Fill(var_flightLenSig);
       VertexChi2KF.at(t).Fill(var_vertexKFChi2);
       MuonglbkinkSum.at(t).Fill(var_sumMuTrkKinkChi2);
       Muon_segmentCompatibility_min.at(t).Fill(MinSegmentCompatibility,w);
       Muon_HCALCompatibility_min.at(t).Fill(MinMIPLikelihood,w);
		 minMudR.at(t).Fill(std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)}));
		 Mu1TauPTRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt());
		 dRMaxMuTau.at(t).Fill(std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)}));
		 MuPair_vertex_chi2_min.at(t).Fill(std::min({Ntp->Vertex_pair_quality(final_idx,0), Ntp->Vertex_pair_quality(final_idx,1), Ntp->Vertex_pair_quality(final_idx,2)}));
		 TauEta.at(t).Fill(TauLV.Eta());
		 VertexDCAMax.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)}));
		 Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(final_idx),w);
		 VertexMuMaxD0SigReco.at(t).Fill(MaxD0Significance,1);
	    }
    }
}


void  FillMVATree::Finish(){
 
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
    file= new TFile("FillMVATreeInput.root","recreate");
    TMVA_Tree->SetDirectory(file);

    file->Write();
    file->Close();

    Selection::Finish();

}
