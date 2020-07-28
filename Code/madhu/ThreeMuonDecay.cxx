#include "ThreeMuonDecay.h"
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

ThreeMuonDecay::ThreeMuonDecay(TString Name_, TString id_):
    Selection(Name_,id_),
    tauMinMass_(1.731),
    tauMaxMass_(1.823),
    tauMinSideBand_(1.61),
    tauMaxSideBand_(2.05),
    tauMassResCutLow(0.007),
    tauMassResCutHigh(0.01)
{
  // This is a class constructor;
}

ThreeMuonDecay::~ThreeMuonDecay(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ThreeMuonDecay::Configure(){
  ////////////////////////////////////////////////////////////////////////
  // Here you can defined your cuts. There are three vector for cuts:
  // std::vector<double> cut, std::vector<double> value and std::vecto<bool> pass.  
  // For exmaple if you want to aplly a  selection to variables Var1,Var2,Var3
  // with selection values Val1, Val2, Val3, then you have: vector cut = (Val1,Val2,Val3),
  // vector value  = (actual_value1, actual_value2, actual_value3), where the actuala_value
  // is an actual value of Var1,2,3 in a given event (this vector will be filled later)
  // vector pass contains boolean of the cut status.

  for(int i=0; i<NCuts;i++){             //"value" store the values from MC/data which must be above the "cut" value. "pass" stores whether "value" passes "cut"
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==SignalCandidate)  cut.at(SignalCandidate)=1;              //Checks whether there's at least one 3-muon combination
    if(i==HLT)              cut.at(HLT)=1;
    if(i==L1)               cut.at(L1)=1;
    if(i==Mu1PtCut)         cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)         cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)         cut.at(Mu3PtCut)=2.0;
    if(i==TRKLWithM)        cut.at(TRKLWithM)=1;                    //Checks whether tracker layer with measurement is a layer > 7
    if(i==MuonID)           cut.at(MuonID)=1;                       //Checks whether they're all global muons (and particle flow muons)
    if(i==PhiVeto1)           cut.at(PhiVeto1)=0; // defined below
    if(i==OmegaVeto1)         cut.at(OmegaVeto1)=0; // defined below
    if(i==PhiVeto2)           cut.at(PhiVeto2)=0; // defined below
    if(i==OmegaVeto2)         cut.at(OmegaVeto2)=0; // defined below
    if(i==TriggerMatch)     cut.at(TriggerMatch)=0.03;
    if(i==ThreeMuMass)      cut.at(ThreeMuMass)=1;// true for MC and mass side band for data
    if(i==CutCategory)      cut.at(CutCategory)=2;// Depends on tau mass resolution

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==SignalCandidate){
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==L1){
      title.at(i)="Pass L1";
      hlabel="L1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT){
      title.at(i)="Pass HLT";
      hlabel="HLT";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
    else if(i==TRKLWithM){
      title.at(i)="TrkLayerWithM";
      hlabel="pass TRKLWithM";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TRKLWithM_",htitle,20,-0.5,19.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TRKLWithM_",htitle,20,-0.5,19.5,hlabel,"Events"));
    }
    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="gl,gl,gl";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PhiVeto1){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto1){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
    }
    else if(i==PhiVeto2){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto2){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
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

   else if(i==CutCategory){
      title.at(i)="CutCategory";
      hlabel="CutCategory";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_CutCategory_",htitle,3,0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_CutCategory_",htitle,3,0.5,3.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams

  LeadMuonPt=HConfig.GetTH1D(Name+"_LeadMuonPt","LeadMuonPt",40,0,20,"p_{T}(#mu_{1}), GeV","Events");
  LeadMuonEta=HConfig.GetTH1D(Name+"_LeadMuonEta","LeadMuonEta",40,-2.6,2.6,"#eta(#mu_{1})","Events");
  LeadMuonPhi=HConfig.GetTH1D(Name+"_LeadMuonPhi","LeadMuonPhi",40,-3.15,3.15,"#phi(#mu_{1})","Events");


  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ThreeMuonDecay::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&LeadMuonPt);
  Extradist1d.push_back(&LeadMuonEta);
  Extradist1d.push_back(&LeadMuonPhi);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ThreeMuonDecay::doEvent(){ 
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Here the index t belongs to sample type, this value originally filled at Ntuple filling
  // level: https://github.com/T3MuAnalysisTools/DsTau23Mu/blob/master/T3MNtuple/interface/DataMCType.h
  // but you can flexibly redefined this and make combinations like, Data, MC1, MC2, MC3+MC4, etc ...
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  
  
  bool HLTOk(false);
  bool L1Ok(false);
  
    bool DoubleMu0Fired(false);
    bool DoubleMu4Fired(false);
    bool DoubleMuFired(false);
    bool TripleMuFired(false);
    bool randomFailed(false);
  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);

    if(HLT.Contains("DoubleMu3_TkMu_DsTau3Mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}

  }

  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
      TString L1TriggerName = Ntp->L1Name(il1);
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
      if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
      if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
      if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
	randomFailed = true;
      }
  }
  
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (!DoubleMu0Fired && !TripleMuFired && randomFailed) l1FailedRandom++;
  
  if(DoubleMuFired  or TripleMuFired) L1Ok = true;
  
    if (HLTOk) value.at(HLT) = true;
    else value.at(HLT) = false;

    if (L1Ok) value.at(L1) = true;
    else value.at(L1) = false;

    //    if (HLTFired && !L1Fired && !randomFailed) cout<<"wrong hlt"<<endl;

    pass.at(HLT) = (value.at(HLT) == cut.at(HLT));
    pass.at(L1)  = (value.at(L1) == cut.at(L1));


  
  
  value.at(SignalCandidate)=0;
  value.at(TriggerMatch)=0;
  
  unsigned int  signal_idx=0;

  double min_chi2(99.);
  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }


  unsigned int  final_idx=0;
  
    value.at(TRKLWithM) = 0;
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


      //      std::cout<<"  "<<Ntp->Muon_isGlobalMuon(mu1_pt_idx)<<"  "<< Ntp->Muon_isTrackerMuon(mu1_pt_idx) << std::endl;

      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
			  Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
			  Ntp->Muon_isGlobalMuon(mu3_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu1_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu2_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu3_pt_idx));
      value.at(TRKLWithM) = (Ntp->Muon_trackerLayersWithMeasurement(mu1_pt_idx) >= 7 ? 1:0 &&
			     Ntp->Muon_trackerLayersWithMeasurement(mu2_pt_idx) >= 7 ? 1:0 &&
			     Ntp->Muon_trackerLayersWithMeasurement(mu3_pt_idx) >= 7 ? 1:0 );
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



      value.at(PhiVeto1) =  M_osss1;//M_osss1;
      value.at(PhiVeto2) =  M_osss2;
      value.at(OmegaVeto1) = M_osss1;
      value.at(OmegaVeto2) = M_osss2;


      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      std::vector<unsigned int> EtaSortedIndices;


      EtaSortedIndices.push_back(Muon_Eta_index_1);
      EtaSortedIndices.push_back(Muon_Eta_index_2);
      EtaSortedIndices.push_back(Muon_Eta_index_3);
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007) value.at(CutCategory)=1;
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false)< 0.01) value.at(CutCategory)=2;
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01) value.at(CutCategory)=3;


      vector<TLorentzVector> trigobjTriplet;
      for (int i=0; i<Ntp->NTriggerObjects(); i++){
	      TString name = Ntp->TriggerObject_name(i);
	      if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
    	  TLorentzVector tmp;
	      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
    	  trigobjTriplet.push_back(tmp);
      }

      std::vector<TLorentzVector> muonTriplet;
      muonTriplet.push_back(Ntp->Muon_P4(mu1_pt_idx));
      muonTriplet.push_back(Ntp->Muon_P4(mu2_pt_idx));
      muonTriplet.push_back(Ntp->Muon_P4(mu3_pt_idx));

      bool triggerCheck = false;
      if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet);
      value.at(TriggerMatch) = triggerCheck;


 
      //      value.at(TriggerMatchMu) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0);
      //      value.at(TriggerMatchMu2) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1);
      //      value.at(TriggerMatchMu3) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2);
     
      value.at(ThreeMuMass) = TauLV.M();
    }
  
  


    pass.at(SignalCandidate) = (value.at(SignalCandidate) == cut.at(SignalCandidate));
    pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
    pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
    pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
    pass.at(MuonID)  =  (value.at(MuonID)     == cut.at(MuonID));
    pass.at(TRKLWithM) = (value.at(TRKLWithM) == cut.at(TRKLWithM));
    pass.at(TriggerMatch) = (value.at(TriggerMatch) == cut.at(TriggerMatch));
    pass.at(PhiVeto1) = (value.at(PhiVeto1) < 0.98 || value.at(PhiVeto1) > 1.06 );
    pass.at(OmegaVeto1) = (value.at(OmegaVeto1) < 0.742 || value.at(OmegaVeto1) > 0.822 );
    pass.at(PhiVeto2) = (value.at(PhiVeto2) < 0.98 || value.at(PhiVeto2) > 1.06 );
    pass.at(OmegaVeto2) = (value.at(OmegaVeto2) < 0.742 || value.at(OmegaVeto2) > 0.822 );
    pass.at(CutCategory) = true;//( value.at(CutCategory) == cut.at(CutCategory) );
    pass.at(ThreeMuMass) =( (value.at(ThreeMuMass) > tauMinSideBand_) &&  (value.at(ThreeMuMass) < tauMaxSideBand_));


  double wobs=1;
  double w;  //  This is an event weights, one may intorduce any weights to the event, for exmaple PU. 
             //  there can be several weights, e.g. w = w1*w2*w3 ...
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  
  
  bool status=AnalysisCuts(t,w,wobs);
  
  std::cout<<"The status is: "<< status <<" with signal candidate, "<<value.at(SignalCandidate)<<" ThreeMuMass, "<<value.at(ThreeMuMass)<<" MuonID, "<<value.at(MuonID)<< std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  // The status boolean is true if all elements in pass are true
  // and false if at least one is false: status = true if 
  // pass = (true, true, true ..., true)  and status = false
  // if pass = (true, true, false, ..., true)

  if(status){ // Only selected events pass this if statement
    // Lets fill below some plots ...
    // All available get functions can be found in https://github.com/T3MuAnalysisTools/Tools/blob/master/Code/Ntuple_Controller.h
    // Lets plot the pT, phi and eta of the leading muon
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon
    LeadMuonPt.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Pt(),1);
    LeadMuonEta.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Eta(),1);
    LeadMuonPhi.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Phi(),1);
    
    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
    
    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
    
    std::cout<<"------------------------------- "<< std::endl;
    std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
    std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
    std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
    Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
    Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
    Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
    Ntp->printMCDecayChainOfEvent(true, true, true, true);
    std::cout<< "\n\n\n\n\n\n";


  }
}


void  ThreeMuonDecay::Finish(){
  Selection::Finish();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}






