#include "SignalSelectorArunTwo.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

SignalSelectorArunTwo::SignalSelectorArunTwo(TString Name_, TString id_):
   Selection(Name_,id_),
   tauMinMass_(1.73),
   tauMaxMass_(1.82),
   tauMinSideBand_(1.65),
   tauMaxSideBand_(2.02),
   tauMassResCutLow(0.007),
   tauMassResCutHigh(0.01),
   phiVetoSigma(0.03),
   omegaVetoSigma(0.04)
{
   // This is a class constructor;
}

SignalSelectorArunTwo::~SignalSelectorArunTwo(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	   << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
     << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SignalSelectorArunTwo::Configure(){


   // Set tree branches
   TMVA_Tree= new TTree("tree","tree");
   TMVA_Tree->Branch("MC",&MC);
   TMVA_Tree->Branch("category",&category);
   TMVA_Tree->Branch("threeGlobal",&threeGlobal);

   //commmon variables (2016 + 2017)
   TMVA_Tree->Branch("var_vertexKFChi2",&var_vertexKFChi2); // <= should be changed to normalized KF chi2
   TMVA_Tree->Branch("var_svpvTauAngle",&var_svpvTauAngle); 
   TMVA_Tree->Branch("var_flightLenSig",&var_flightLenSig);
   TMVA_Tree->Branch("var_segCompMuMin",&var_segCompMuMin);

   // 2016 variables
   TMVA_Tree->Branch("var_pmin", &var_pmin); // Minimum p of the three muons
   TMVA_Tree->Branch("var_max_cLP", &var_max_cLP); // Maximum chi square of the STA-TRK matching
   TMVA_Tree->Branch("var_max_tKink", &var_max_tKink); // Maximum of the track kink of the 3 muons
   TMVA_Tree->Branch("var_MinD0Significance", &var_MinD0Significance); // Minimum of the transverse IP significance of the 3 muons
   TMVA_Tree->Branch("var_mindca_iso", &var_mindca_iso); // Minimum DCA of tracks to muons with pT > 1 GeV (Maximum of three muons)
   TMVA_Tree->Branch("var_trk_relPt", &var_trk_relPt); // Ratio of sum of Pt of the tracks in muon isolation to muon (max value) [trk_pt>1 GeV, dR<0.03, dca<1 mm]
   TMVA_Tree->Branch("var_minMatchedStations", &var_minMatchedStations); // number of minimum matched stations

   // 2017 variables
   TMVA_Tree->Branch("var_MuMu_minKFChi2",&var_MuMu_minKFChi2);
   TMVA_Tree->Branch("var_MuTau_maxdR",&var_MuTau_maxdR);
   TMVA_Tree->Branch("var_sumMuTrkKinkChi2",&var_sumMuTrkKinkChi2); // sum of chi square of STA-TRK matching of 3 muons
   TMVA_Tree->Branch("var_MaxD0Significance", &var_MaxD0Significance); // Maximum of the transverse IP significance of the 3 muons
   TMVA_Tree->Branch("var_MinMIPLikelihood", &var_MinMIPLikelihood); //Calo compatibility 
   TMVA_Tree->Branch("var_maxdca", &var_maxdca); // max dca between the initial tracks of two muons
   TMVA_Tree->Branch("var_MuMu_mindR", &var_MuMu_mindR);
   TMVA_Tree->Branch("var_RelPt_Mu1Tau",&var_RelPt_Mu1Tau);
   TMVA_Tree->Branch("var_MuTau_maxdR",&var_MuTau_maxdR);

   // Include calo energy in the HCAL tower

   // Spectator variables
   TMVA_Tree->Branch("var_Eta_Tau",&var_Eta_Tau); 
   TMVA_Tree->Branch("var_tauMass",&var_tauMass);
   TMVA_Tree->Branch("var_tauMassRes", &var_tauMassRes);
   TMVA_Tree->Branch("var_tauMassRefit", &var_tauMassRefit);
   // -----------------
   
  ////////////////////////////////////////////////////////////////////////
  // Here you can defined your cuts. There are three vector for cuts:
  // std::vector<double> cut, std::vector<double> value and std::vecto<bool> pass.  
  // For exmaple if you want to aplly a  selection to variables Var1,Var2,Var3
  // with selection values Val1, Val2, Val3, then you have: vector cut = (Val1,Val2,Val3),
  // vector value  = (actual_value1, actual_value2, actual_value3), where the actuala_value
  // is an actual value of Var1,2,3 in a given event (this vector will be filled later)
  // vector pass contains boolean of the cut status.

  for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==TriggerOk)          cut.at(TriggerOk)=1;
      if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
      if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
      if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
      if(i==TriggerMatchMu1)    cut.at(TriggerMatchMu1)=0.03;
      if(i==TriggerMatchMu2)    cut.at(TriggerMatchMu2)=0.03;
      if(i==TriggerMatchMu3)    cut.at(TriggerMatchMu3)=0.03;
      if(i==MuonID)             cut.at(MuonID)=1;
      if(i==PVRefit)            cut.at(PVRefit)=1;
      if(i==PhiVetoOS1)         cut.at(PhiVetoOS1)=0; // defined below
      if(i==PhiVetoOS2)         cut.at(PhiVetoOS2)=0; // defined below
      if(i==OmegaVetoOS1)       cut.at(OmegaVetoOS1)=0; // defined below
      if(i==OmegaVetoOS2)       cut.at(OmegaVetoOS2)=0; // defined below
      if(i==TauMassCut)         cut.at(TauMassCut)=1; // true for MC and mass side band for data
      if(i==DsGenMatch)         cut.at(DsGenMatch)=0.03; // cut is applied to only Ds->Phi(MuMu)Pi channel
      if(i==GenMatch)	  		  cut.at(GenMatch)=0.03;

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
         title.at(i)="Pass HLT";
         hlabel="DoubleMu3_Trk_Tau3mu";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==SignalCandidate){
         title.at(i)="signal candidates";
         hlabel="3mu candidates";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidates_",htitle,9,1.0,10.0,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidates_",htitle,9,1.0,10.0,hlabel,"Events"));
      }
    else if(i==LeadingMuonPt){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==LeadingMuonPt1){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt1_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt1_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==LeadingMuonPt2){
      title.at(i)="$\\mu$ Pt $>$ .5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of the leading  muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeadingMuonPt2_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeadingMuonPt2_",htitle,80,0,20,hlabel,"Events"));
    }
      else if(i==Mu1PtCut){
         title.at(i)="$p_{T}(\\mu_{1}) >$";
         title.at(i)+=cut.at(Mu1PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="$p_{T}(\\mu_{2}) >$";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
      }
      else if(i==Mu3PtCut){
         title.at(i)="$p_{T}(\\mu_{3}) >$";
         title.at(i)+=cut.at(Mu3PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon3 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,30,0,15,hlabel,"Events"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="gl,gl,(gl/trk)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==PVRefit){
         title.at(i)="PV refit valid";
         hlabel="Primary Vertex refit valid";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PVRefitValid_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
      else if(i==PhiVetoOS1){
         title.at(i)="phi mass veto (OS1)";
         hlabel="Phi mass Veto (OS1), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVetoOS1_",htitle,60,0.8,1.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVetoOS1_",htitle,60,0.8,1.2,hlabel,"Events"));
      }
      else if(i==PhiVetoOS2){
         title.at(i)="phi mass veto (OS2)";
         hlabel="Phi mass Veto (OS2), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVetoOS2_",htitle,60,0.8,1.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVetoOS2_",htitle,60,0.8,1.2,hlabel,"Events"));
      }
      else if(i==OmegaVetoOS1){
         title.at(i)="omega mass veto (OS1)";
         hlabel="Omega mass veto (OS1), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVetoOS1_",htitle,50,0.4,0.9,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVetoOS1_",htitle,50,0.4,0.9,hlabel,"Events"));
      }
      else if(i==OmegaVetoOS2){
         title.at(i)="omega mass veto (OS2)";
         hlabel="Omega mass veto (OS2), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVetoOS2_",htitle,50,0.4,0.9,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVetoOS2_",htitle,50,0.4,0.9,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu1){
         title.at(i)="Trigger Matching (mu1)";
         hlabel="trigger matching dR (#mu1)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu1_",htitle,20,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu1_",htitle,20,0,0.1,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu2){
         title.at(i)="Trigger Matching (mu2)";
         hlabel="trigger matching dR (#mu2)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu2_",htitle,20,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu2_",htitle,20,0,0.1,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu3){
         title.at(i)="Trigger Matching (mu3)";
         hlabel="trigger matching dR (#mu3)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu3_",htitle,20,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu3_",htitle,20,0,0.1,hlabel,"Events"));
      }
      else if(i==TauMassCut){
         title.at(i)="Tau Mass";
         hlabel="three mu mass, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,50,1.60,2.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,50,1.60,2.1,hlabel,"Events"));
      }
      else if(i==DsGenMatch){
         title.at(i)="Ds GEN matching (only dsphipi)";
         hlabel="Ds GEN match (only dsphipi)";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DsGENMatch_",htitle,20,0.0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DsGENMatch_",htitle,20,0.0,0.1,hlabel,"Events"));  
      }
      else if(i==GenMatch){
         title.at(i)="GEN matching";
         hlabel="GEN match";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GENMatch_",htitle,20,0.0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GENMatch_",htitle,20,0.0,0.1,hlabel,"Events"));
      }
      else if(i==OmegaVeto){
         title.at(i)="$\\omega$ mass veto";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Omega mass veto, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto_",htitle,50,0.4,0.9,hlabel,"Events"));
      }
  } 
  // Setup NPassed Histogams

  LeadMuonPt=HConfig.GetTH1D(Name+"_LeadMuonPt","LeadMuonPt",40,0,20,"p_{T}(#mu_{1}), GeV","Events");
  LeadMuonEta=HConfig.GetTH1D(Name+"_LeadMuonEta","LeadMuonEta",40,-2.6,2.6,"#eta(#mu_{1})","Events");
  LeadMuonPhi=HConfig.GetTH1D(Name+"_LeadMuonPhi","LeadMuonPhi",40,-3.15,3.15,"#phi(#mu_{1})","Events");
  
  LeadMuonPt1=HConfig.GetTH1D(Name+"_LeadMuonPt1","LeadMuonPt1",40,0,20,"p_{T}(#mu_{2}), GeV","Events");
  LeadMuonEta1=HConfig.GetTH1D(Name+"_LeadMuonEta1","LeadMuonEta1",40,-2.6,2.6,"#eta(#mu_{2})","Events");
  LeadMuonPhi1=HConfig.GetTH1D(Name+"_LeadMuonPhi1","LeadMuonPhi1",40,-3.15,3.15,"#phi(#mu_{2})","Events");
  
  LeadMuonPt2=HConfig.GetTH1D(Name+"_LeadMuonPt2","LeadMuonPt2",40,0,20,"p_{T}(#mu_{3}), GeV","Events");
  LeadMuonEta2=HConfig.GetTH1D(Name+"_LeadMuonEta2","LeadMuonEta2",40,-2.6,2.6,"#eta(#mu_{3})","Events");
  LeadMuonPhi2=HConfig.GetTH1D(Name+"_LeadMuonPhi2","LeadMuonPhi2",40,-3.15,3.15,"#phi(#mu_{3})","Events");
  
  TauPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",40,0,20,"p_{T}(#tau_{1}), GeV","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",40,-2.6,2.6,"#eta(#tau_{1})","Events");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","TauPhi",40,-3.15,3.15,"#phi(#tau_{1})","Events");
  
  InvMu1G=HConfig.GetTH1D(Name+"_InvMu1G","First muon is global",40,0,0.5,"LEading Muon is Global","Events");
  InvMu1T=HConfig.GetTH1D(Name+"_InvMu1T","First muon is tracker",40,0,0.5,"LEading Muon is Tracker","Events");
  InvMu1S=HConfig.GetTH1D(Name+"_InvMu1S","First muon is standalone",40,0,0.5,"LEading Muon is Standalone","Events");
  
  InvMu2G=HConfig.GetTH1D(Name+"_InvMu2G","Second muon is global",40,0,0.5,"LEading Muon 1 is Global","Events");
  InvMu2T=HConfig.GetTH1D(Name+"_InvMu2T","Second muon is tracker",40,0,0.5,"LEading Muon 1 is Tracker","Events");
  InvMu2S=HConfig.GetTH1D(Name+"_InvMu2S","Second muon is standalone",40,0,0.5,"LEading Muon 1 is Standalone","Events");
  
  InvMu3G=HConfig.GetTH1D(Name+"_InvMu3G","Third muon is global",40,0,0.5,"LEading Muon 2 is Global","Events");
  InvMu3T=HConfig.GetTH1D(Name+"_InvMu3T","Third muon is tracker",40,0,0.5,"LEading Muon 3 is Tracker","Events");
  InvMu3S=HConfig.GetTH1D(Name+"_InvMu3S","Third muon is standalone",40,0,0.5,"LEading Muon 3 is Standalone","Events");
  
  InvTG=HConfig.GetTH1D(Name+"_InvTG","Tau is global",40,0,3,"Tau is Global","Events");
  InvTT=HConfig.GetTH1D(Name+"_InvTT","Tau is tracker",40,0,3,"Tau is Tracker","Events");
  InvTS=HConfig.GetTH1D(Name+"_InvTS","Tau is standalone",40,0,3,"Tau is Standalone","Events");  


  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");
  
  TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");
  
  
   //for Trigger
  
   Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
   L1Seed=HConfig.GetTH1D(Name+"_L1Seeds","L1Seed",3,-0.5,2.5,"DoubleMu/TripleMu","Events");
   VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
   MuonglbkinkSum  =HConfig.GetTH1D(Name+"_MuonglbkinkSum","MuonglbkinkSum",50,0.,50," #sum  #mu glb kink #chi^{2}","Events");
   FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");
   SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");
   Muon_segmentCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_min","Muon_segmentCompatibility_min",50,0.,1,"Inner Track and muon segment match min ","Events");
   Muon_HCALCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_HCALCompatibility_min","Muon_ECALCompatibility_min",50,0.,1,"MIP Likelihood min ","Events");
   //New variables
   minMudR = HConfig.GetTH1D(Name+"_minMudR","minMudR",30,0,0.3,"min #DeltaR(#mu#mu)","Events");
   Mu1TauPTRatio = HConfig.GetTH1D(Name+"_Mu1TauPTRatio","Mu1TauPTRatio",40,0.1,0.9,"p_{T}(#mu_{1})/p_{T}(#tau)","Events");
   dRMaxMuTau = HConfig.GetTH1D(Name+"dRMaxMuTau","dRMaxMuTau",50,0,0.5,"max #DeltaR(#mu#tau)","Events");
   MuPair_vertex_chi2_min=HConfig.GetTH1D(Name+"_MuPair_vertex_chi2_min","MuPair_vertex_chi2_min",50,0,1.5,"KF min #chi^{2} of #mu pair","Events");
   TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
   VertexDCAMax=HConfig.GetTH1D(Name+"_VertexDCAMax","VertexDCAMax",40,0,0.15,"Max closest distance between muons","Events");
   Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","Isolation_MinDist",50,0,0.1,"Iso MinDist","Events"); 
   VertexMuMaxD0SigReco=HConfig.GetTH1D(Name+"_VertexMuMaxD0SigReco","VertexMuMaxD0SigReco",50,0,4,"#mu - PV max transverse distance significance","Events");
   EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
   L1TriggerMatch=HConfig.GetTH2D(Name+"_L1TriggerMatch","L1 vs TriggerMatch",4,-0.5,3.5,4,-0.5,3.5,"Muon matched to trObj","DoubleMu/TripleMu");
   Test = HConfig.GetTH2D(Name+"_Test","L1 vs trigger match",4,-0.5,3.5,4,-0.5,3.5,"Muon matched to trObj","DoubleMu/TripleMu");
   TrackCorrelation= HConfig.GetTH2D(Name+"_TrackCorrelation","N tracks vs trigger match",4,-0.5,3.5,25,0,25,"Muon matched to trObj","N tracks");

   // Candidates Info
   NTwoMuonsTrack=HConfig.GetTH1D(Name+"_NTwoMuonsTrack","Two Muons and Track Candidates",10,0,10,"2 Muons + Track candidates","Events");
   NMuons=HConfig.GetTH1D(Name+"_NMuons","Number of Muons",10,0,10,"N Muons","Events");
   NTracks=HConfig.GetTH1D(Name+"_NTracks","Number of Tracks",30,0,30,"N Tracks","Events");



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SignalSelectorArunTwo::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&LeadMuonPt);
  Extradist1d.push_back(&LeadMuonEta);
  Extradist1d.push_back(&LeadMuonPhi);
  
  Extradist1d.push_back(&LeadMuonPt1);
  Extradist1d.push_back(&LeadMuonEta1);
  Extradist1d.push_back(&LeadMuonPhi1);
  
  Extradist1d.push_back(&LeadMuonPt2);
  Extradist1d.push_back(&LeadMuonEta2);
  Extradist1d.push_back(&LeadMuonPhi2);
  
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPhi);
  
  Extradist1d.push_back(&InvMu1G);
  Extradist1d.push_back(&InvMu1T);
  Extradist1d.push_back(&InvMu1S);
  
  Extradist1d.push_back(&InvMu2G);
  Extradist1d.push_back(&InvMu2T);
  Extradist1d.push_back(&InvMu2S);
  
  Extradist1d.push_back(&InvMu3G);
  Extradist1d.push_back(&InvMu3T);
  Extradist1d.push_back(&InvMu3S);
  
  Extradist1d.push_back(&InvTG);
  Extradist1d.push_back(&InvTT);
  Extradist1d.push_back(&InvTS);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  SignalSelectorArunTwo::doEvent(){ 

   TLorentzVector Tau_P4;
   TLorentzVector TauG_P4;
   TLorentzVector TauT_P4;
   TLorentzVector TauS_P4;

   value.at(TriggerOk)=0;
   value.at(SignalCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(Mu3PtCut)=0;
   value.at(MuonID)=0;
   value.at(TriggerMatchMu1)=1.0;
   value.at(TriggerMatchMu2)=1.0;
   value.at(TriggerMatchMu3)=1.0;
   value.at(PVRefit)=0;
   value.at(PhiVetoOS1)=99.0;
   value.at(PhiVetoOS2)=99.0;
   value.at(OmegaVetoOS1)=99.0;
   value.at(OmegaVetoOS2)=99.0;
   value.at(TauMassCut)=0;
   value.at(DsGenMatch)=1;
   value.at(GenMatch)=1;
      
   unsigned int t;
   int id(Ntp->GetMCID());
   
   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
   bool HLTOk(false);
   bool L1Ok(false);
   
   for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);

      if(id==1){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }
      if( (id==1 && Ntp->WhichEra(2017).Contains("RunF")) || (id==1 && Ntp->WhichEra(2018).Contains("Run")) ){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v" or HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v") ) && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }
      if(id!=1){
         if(HLT.Contains("DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) HLTOk=true;
      }

   }

   // ------------ HLT selection -----------------
   bool DoubleMuFired(false);
   bool TripleMuFired(false);
   
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
      if (id==1 && Ntp->WhichEra(2018).Contains("Run")){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))							DoubleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))	TripleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2"))								DoubleMuFired = Ntp->L1Decision(il1);
      }
      if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
      }
   }



   if (DoubleMuFired || TripleMuFired) L1Ok = true;
   if (L1Ok && HLTOk) value.at(TriggerOk) = true;
   else value.at(TriggerOk) = false;
   pass.at(TriggerOk) = (value.at(TriggerOk) == cut.at(TriggerOk));
   
   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_tau = 99.0;

   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_tau = 0.;
   double sumPtTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_tau = 0;
   
   // Selection of the best candidate

   vector<unsigned int> selectedIndices;
   vector<unsigned int> candidateRank;
   unsigned int final_idx = 0;
   double minChiSq = 999.0;
   value.at(SignalCandidate) = Ntp->NThreeMuons();

   if (Ntp->NThreeMuons()>0){
      for (size_t j=0; j<Ntp->NThreeMuons(); ++j){ 
         Muon_index_1 = Ntp->ThreeMuonIndices(j).at(0); 
         Muon_index_2 = Ntp->ThreeMuonIndices(j).at(1); 
         Muon_index_3 = Ntp->ThreeMuonIndices(j).at(2);

         // value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
         //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
         //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::CutBasedIdMedium));
         //----------------  alternatively require two leading muons to be global and trailing muon to be tracker

         // Older version of OS pairs
         //vector<unsigned int> idx_vec;
         //idx_vec.push_back(Muon_index_1);
         //idx_vec.push_back(Muon_index_2);
         //idx_vec.push_back(Muon_index_3);
         //unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
         //unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
         //unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);
         unsigned int os_idx=Muon_index_1, ss1_idx=Muon_index_1, ss2_idx=Muon_index_3;
         
         if (Ntp->Muon_charge(Muon_index_1)!=Ntp->Muon_charge(Muon_index_2) && Ntp->Muon_charge(Muon_index_1)!=Ntp->Muon_charge(Muon_index_3)){
            os_idx = Muon_index_1;
            if (Ntp->Vertex_pairfit_status(j,0,false) && Ntp->Vertex_pairfit_status(j,2,false)){
               if (Ntp->Vertex_pair_quality(j,0,false)<Ntp->Vertex_pair_quality(j,2,false)) { ss1_idx = Muon_index_2; ss2_idx = Muon_index_3; }
               else { ss1_idx = Muon_index_3; ss2_idx = Muon_index_2; }
            }
            else {
               ss1_idx = Muon_index_2;
               ss2_idx = Muon_index_3;
            }
         }
         else if (Ntp->Muon_charge(Muon_index_2)!=Ntp->Muon_charge(Muon_index_3) && Ntp->Muon_charge(Muon_index_2)!=Ntp->Muon_charge(Muon_index_1)){
            os_idx = Muon_index_2;
            if (Ntp->Vertex_pairfit_status(j,0,false) && Ntp->Vertex_pairfit_status(j,1,false)){
               if (Ntp->Vertex_pair_quality(j,0,false)<Ntp->Vertex_pair_quality(j,1,false)) { ss1_idx = Muon_index_1; ss2_idx = Muon_index_3; }
               else { ss1_idx = Muon_index_3; ss2_idx = Muon_index_1; }
            }
            else {
               ss1_idx = Muon_index_1;
               ss2_idx = Muon_index_3;
            }
         }
         else if (Ntp->Muon_charge(Muon_index_3)!=Ntp->Muon_charge(Muon_index_1) && Ntp->Muon_charge(Muon_index_3)!=Ntp->Muon_charge(Muon_index_2)){
            os_idx = Muon_index_3;            
            if (Ntp->Vertex_pairfit_status(j,1,false) && Ntp->Vertex_pairfit_status(j,2,false)){
               if (Ntp->Vertex_pair_quality(j,1,false)<Ntp->Vertex_pair_quality(j,2,false)) { ss1_idx = Muon_index_2; ss2_idx = Muon_index_1; }
               else { ss1_idx = Muon_index_1; ss2_idx = Muon_index_2; }
            }
            else {
               ss1_idx = Muon_index_2;
               ss2_idx = Muon_index_1;
            }
         }
         
         M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
         M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

         Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(j)).at(0);
         Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(j)).at(1);
         Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(j)).at(2);

         //
         value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) && 
               Ntp->Muon_isGlobalMuon(Muon_index_2) &&
               (!Ntp->Muon_isGlobalMuon(Muon_index_3) && Ntp->Muon_isTrackerMuon(Muon_index_3)));
         //------------------------------------------------------------------------------------------------------

         if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal = true;
         else threeGlobal = false;
         value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
         value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
         value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();
         value.at(PVRefit) = Ntp->Vertex_RefitPVisValid(j);

         TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)+Ntp->Muon_P4(Muon_index_2)+Ntp->Muon_P4(Muon_index_3);    

         //value.at(PhiVeto) = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
         //value.at(OmegaVeto) = fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;

         value.at(PhiVetoOS1) = M_osss1;
         value.at(PhiVetoOS2) = M_osss2;
         value.at(OmegaVetoOS1) = M_osss1;
         value.at(OmegaVetoOS2) = M_osss2;

         value.at(TriggerMatchMu1) = Ntp->ThreeMuons_TriggerMatch_dR(j).at(0);
         value.at(TriggerMatchMu2) = Ntp->ThreeMuons_TriggerMatch_dR(j).at(1);
         value.at(TriggerMatchMu3) = Ntp->ThreeMuons_TriggerMatch_dR(j).at(2);
         value.at(TauMassCut) = TauLV.M();

         if (id==30) value.at(DsGenMatch) = Ntp->DsGenMatch(j);
         else value.at(DsGenMatch) = 0;

         if (id!=1) value.at(GenMatch) = Ntp->TauGenMatch(j);
         else value.at(GenMatch) = 0;

         pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
         pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
         pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
         pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
         pass.at(MuonID) = (value.at(MuonID)  == cut.at(MuonID));
         pass.at(PVRefit) = (value.at(PVRefit) == cut.at(PVRefit));
         pass.at(TriggerMatchMu1) = (value.at(TriggerMatchMu1) <  cut.at(TriggerMatchMu1));
         pass.at(TriggerMatchMu2) = (value.at(TriggerMatchMu2) <  cut.at(TriggerMatchMu2));
         pass.at(TriggerMatchMu3) = (value.at(TriggerMatchMu3) <  cut.at(TriggerMatchMu3));
         //pass.at(TriggerMatchMu1) = true;
         //pass.at(TriggerMatchMu2) = true;
         //pass.at(TriggerMatchMu3) = true;
         if (threeGlobal) {
            pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma);
            pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma);
         }
         else {
            pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS1)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.02); // extend the veto region for two global + tracker category
            pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.02); // extend the veto region for two global + tracker category
         }
         pass.at(OmegaVetoOS1) = (fabs(value.at(OmegaVetoOS1)-PDG_Var::Omega_mass())> omegaVetoSigma);
         pass.at(OmegaVetoOS2) = (fabs(value.at(OmegaVetoOS2)-PDG_Var::Omega_mass())> omegaVetoSigma);
         pass.at(DsGenMatch) = ( value.at(DsGenMatch) < cut.at(DsGenMatch) );
         pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
         if(id!=1) pass.at(TauMassCut) = ( value.at(TauMassCut) > (tauMinMass_+0.02) && value.at(TauMassCut) < (tauMaxMass_-0.02));
         else  pass.at(TauMassCut) = ( (value.at(TauMassCut) > tauMinSideBand_ && value.at(TauMassCut) < tauMinMass_)  || (value.at(TauMassCut)> tauMaxMass_ && value.at(TauMassCut) < tauMaxSideBand_));
         unsigned int score = 0;
         for (unsigned int k=0; k<NCuts; ++k) if (pass.at(k)) score++;

         if (score==NCuts) selectedIndices.push_back(j);
         candidateRank.push_back(score);
      }

      if (selectedIndices.size()==0) final_idx=std::distance(candidateRank.begin(), std::max_element(candidateRank.begin(), candidateRank.end()));
      else{
         for (size_t j=0; j<selectedIndices.size(); ++j){
            int _ = selectedIndices.at(j);
            double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(_);
            if (tmpchi<minChiSq) {
               minChiSq = tmpchi;
               final_idx = _;
            }
         }
      }

      Muon_index_1 = Ntp->ThreeMuonIndices(final_idx).at(0); 
      Muon_index_2 = Ntp->ThreeMuonIndices(final_idx).at(1); 
      Muon_index_3 = Ntp->ThreeMuonIndices(final_idx).at(2);

      // value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::CutBasedIdMedium));
      //----------------  alternatively require two leading muons to be global and trailing muon to be tracker
      unsigned int os_idx=Muon_index_1, ss1_idx=Muon_index_2, ss2_idx=Muon_index_3;
      
      if (Ntp->Muon_charge(Muon_index_1)!=Ntp->Muon_charge(Muon_index_2) && Ntp->Muon_charge(Muon_index_1)!=Ntp->Muon_charge(Muon_index_3)){
         os_idx = Muon_index_1;
         if (Ntp->Vertex_pairfit_status(final_idx,0,false) && Ntp->Vertex_pairfit_status(final_idx,2,false)){
            if (Ntp->Vertex_pair_quality(final_idx,0,false)<Ntp->Vertex_pair_quality(final_idx,2,false)) { ss1_idx = Muon_index_2; ss2_idx = Muon_index_3; }
            else { ss1_idx = Muon_index_3; ss2_idx = Muon_index_2; }
         }
         else {
            ss1_idx = Muon_index_2;
            ss2_idx = Muon_index_3;
         }
      }
      else if (Ntp->Muon_charge(Muon_index_2)!=Ntp->Muon_charge(Muon_index_3) && Ntp->Muon_charge(Muon_index_2)!=Ntp->Muon_charge(Muon_index_1)){
         os_idx = Muon_index_2;
         if (Ntp->Vertex_pairfit_status(final_idx,0,false) && Ntp->Vertex_pairfit_status(final_idx,1,false)){
            if (Ntp->Vertex_pair_quality(final_idx,0,false)<Ntp->Vertex_pair_quality(final_idx,1,false)) { ss1_idx = Muon_index_1; ss2_idx = Muon_index_3; }
            else { ss1_idx = Muon_index_3; ss2_idx = Muon_index_1; }
         }
         else {
            ss1_idx = Muon_index_1;
            ss2_idx = Muon_index_3;
         }
      }
      else if (Ntp->Muon_charge(Muon_index_3)!=Ntp->Muon_charge(Muon_index_1) && Ntp->Muon_charge(Muon_index_3)!=Ntp->Muon_charge(Muon_index_2)){
         os_idx = Muon_index_3;            
         if (Ntp->Vertex_pairfit_status(final_idx,1,false) && Ntp->Vertex_pairfit_status(final_idx,2,false)){
            if (Ntp->Vertex_pair_quality(final_idx,1,false)<Ntp->Vertex_pair_quality(final_idx,2,false)) { ss1_idx = Muon_index_2; ss2_idx = Muon_index_1; }
            else { ss1_idx = Muon_index_1; ss2_idx = Muon_index_2; }
         }
         else {
            ss1_idx = Muon_index_2;
            ss2_idx = Muon_index_1;
         }
      }

      M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
      M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

      value.at(PhiVetoOS1) = M_osss1;
      value.at(PhiVetoOS2) = M_osss2;
      value.at(OmegaVetoOS1) = M_osss1;
      value.at(OmegaVetoOS2) = M_osss2;

      Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)+Ntp->Muon_P4(Muon_index_2)+Ntp->Muon_P4(Muon_index_3);

      //
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) && 
            Ntp->Muon_isGlobalMuon(Muon_index_2) &&
            (!Ntp->Muon_isGlobalMuon(Muon_index_3) && Ntp->Muon_isTrackerMuon(Muon_index_3)));
      //------------------------------------------------------------------------------------------------------

      if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal = true;
      else threeGlobal = false;
      value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();
      value.at(PVRefit) = Ntp->Vertex_RefitPVisValid(final_idx);

      /*
         vector<unsigned int> idx_vec;
         idx_vec.push_back(Muon_index_1);
         idx_vec.push_back(Muon_index_2);
         idx_vec.push_back(Muon_index_3);
         unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
         unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
         unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);
         M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
         M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();
         */

      value.at(TriggerMatchMu1) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0);
      value.at(TriggerMatchMu2) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1);
      value.at(TriggerMatchMu3) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2);
      value.at(TauMassCut) = TauLV.M();

      if (id==30) value.at(DsGenMatch) = Ntp->DsGenMatch(final_idx);
      else value.at(DsGenMatch) = 0;

      if (id!=1) value.at(GenMatch) = Ntp->TauGenMatch(final_idx);
      else value.at(GenMatch) = 0;

      pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
      pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
      pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
      pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
      pass.at(MuonID) = (value.at(MuonID)  == cut.at(MuonID));
      pass.at(PVRefit) = (value.at(PVRefit) == cut.at(PVRefit));
      pass.at(TriggerMatchMu1) = (value.at(TriggerMatchMu1) <  cut.at(TriggerMatchMu1));
      pass.at(TriggerMatchMu2) = (value.at(TriggerMatchMu2) <  cut.at(TriggerMatchMu2));
      pass.at(TriggerMatchMu3) = (value.at(TriggerMatchMu3) <  cut.at(TriggerMatchMu3));
      //pass.at(TriggerMatchMu1) = true;
      //pass.at(TriggerMatchMu2) = true;
      //pass.at(TriggerMatchMu3) = true;
      if (threeGlobal) {
         pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma);
         pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma);
      }
      else {
         pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS1)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.02); // extend the veto region for two global + tracker category
         pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.02); // extend the veto region for two global + tracker category
      }
      pass.at(OmegaVetoOS1) = (fabs(value.at(OmegaVetoOS1)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(OmegaVetoOS2) = (fabs(value.at(OmegaVetoOS2)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(DsGenMatch) = ( value.at(DsGenMatch) < cut.at(DsGenMatch) );
      pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );

      if(id!=1) pass.at(TauMassCut) = ( value.at(TauMassCut) > (tauMinMass_+0.02) && value.at(TauMassCut) < (tauMaxMass_-0.02));
      else  pass.at(TauMassCut) = ( (value.at(TauMassCut) > tauMinSideBand_ && value.at(TauMassCut) < tauMinMass_)  || (value.at(TauMassCut)> tauMaxMass_ && value.at(TauMassCut) < tauMaxSideBand_));

      // Find how many muons are matched to the trigger leg
      std::vector<unsigned int> triggerObjectMatches;
      triggerObjectMatches.push_back(TriggerMatchMu1);
      triggerObjectMatches.push_back(TriggerMatchMu2);
      triggerObjectMatches.push_back(TriggerMatchMu3);

      int nTriggerMatches = 0;
      int l1 = 0;

      l1 = 2*(TripleMuFired)+DoubleMuFired;

      for (size_t j=0; j<triggerObjectMatches.size(); ++j){
         if (pass.at(triggerObjectMatches.at(j))) nTriggerMatches++;
      }

      if (passAllBut(triggerObjectMatches)){
         TrackCorrelation.at(t).Fill(nTriggerMatches, Ntp->NTracks());
         L1TriggerMatch.at(t).Fill(nTriggerMatches, l1);
         Test.at(t).Fill(nTriggerMatches,l1);
      }
   }
   
   
   value.at(PrimeVtx)=Ntp->NVtx(); // Here the actual_value of a cut is set
   pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); // Here we check that the actuall value of PrimeVrtices is above 5.
  
   if(Ntp->NThreeMuons()>0){  // Check if this is a signal category (take the first triplet only in this example)

    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon 0
    value.at(LeadingMuonPt) = Ntp->Muon_P4(mu1_pt_idx).Pt();
    
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);  // leading pT muon 1
    value.at(LeadingMuonPt1) = Ntp->Muon_P4(mu2_pt_idx).Pt();
    
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);  // leading pT muon 2
    value.at(LeadingMuonPt2) = Ntp->Muon_P4(mu3_pt_idx).Pt();
    
    Tau_P4=Ntp->Muon_P4(mu1_pt_idx)+Ntp->Muon_P4(mu2_pt_idx)+Ntp->Muon_P4(mu3_pt_idx);

   }

   pass.at(LeadingMuonPt)   = (value.at(LeadingMuonPt)    > cut.at(LeadingMuonPt));
   pass.at(LeadingMuonPt1)   = (value.at(LeadingMuonPt1)    > cut.at(LeadingMuonPt1));
   pass.at(LeadingMuonPt2)   = (value.at(LeadingMuonPt2)    > cut.at(LeadingMuonPt2));
   
   double wobs=1;
   double w; 

   if(!Ntp->isData()) { w=1; } //  No weights to data
   else{w=1;}

   bool status=AnalysisCuts(t,w,wobs);
   
  
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
    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);  // leading pT muon 0
    LeadMuonPt.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Pt(),1);
    LeadMuonEta.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Eta(),1);
    LeadMuonPhi.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).Phi(),1);
    
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);  // leading pT muon 1
    LeadMuonPt1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Pt(),1);
    LeadMuonEta1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Eta(),1);
    LeadMuonPhi1.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).Phi(),1);
    
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);  // leading pT muon 2
    LeadMuonPt2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Pt(),1);
    LeadMuonEta2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Eta(),1);
    LeadMuonPhi2.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).Phi(),1); 
    
    TauPt.at(t).Fill(Tau_P4.Pt(),1);
    TauEta.at(t).Fill(Tau_P4.Eta(),1);
    TauPhi.at(t).Fill(Tau_P4.Phi(),1);
    
    if(Ntp->Muon_isGlobalMuon(mu1_pt_idx)){
    InvMu1G.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu1_pt_idx)){
    InvMu1T.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu1_pt_idx)){
    InvMu1S.at(t).Fill(Ntp->Muon_P4(mu1_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu2_pt_idx)){
    InvMu2G.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu2_pt_idx)){
    InvMu2T.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu2_pt_idx)){
    InvMu2S.at(t).Fill(Ntp->Muon_P4(mu2_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu3_pt_idx)){
    InvMu3G.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu3_pt_idx)){
    InvMu3T.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu3_pt_idx)){
    InvMu3S.at(t).Fill(Ntp->Muon_P4(mu3_pt_idx).M(),1);
    }
    
    
    
    if(Ntp->Muon_isGlobalMuon(mu1_pt_idx)&&Ntp->Muon_isGlobalMuon(mu2_pt_idx)&&Ntp->Muon_isGlobalMuon(mu3_pt_idx)){
    InvTG.at(t).Fill(Tau_P4.M(),1);
    }
    
    if(Ntp->Muon_isTrackerMuon(mu1_pt_idx)&&Ntp->Muon_isTrackerMuon(mu2_pt_idx)&&Ntp->Muon_isTrackerMuon(mu3_pt_idx)){
    InvTT.at(t).Fill(Tau_P4.M(),1);
    }
    
    if(Ntp->Muon_isStandAloneMuon(mu1_pt_idx)&&Ntp->Muon_isStandAloneMuon(mu2_pt_idx)&&Ntp->Muon_isStandAloneMuon(mu3_pt_idx)){
    InvTS.at(t).Fill(Tau_P4.M(),1);
    }


  }
}


void  SignalSelectorArunTwo::Finish(){
  Selection::Finish();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}






