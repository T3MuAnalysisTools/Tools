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
   tauMinMass_(1.75),
   tauMaxMass_(1.80),
   tauMinSideBand_(1.65),
   tauMaxSideBand_(1.90),
   tauMassResCutLow(0.007),
   tauMassResCutHigh(0.01),
   phiVetoSigma(0.03),
   omegaVetoSigma(0.04)
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
   TMVA_Tree->Branch("threeGlobal",&threeGlobal);
   TMVA_Tree->Branch("l1seed",&l1seed);

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
      if(i==TriggerOk){
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
   } 

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

   // Muon Kinematic plots
   Muon1PtEta=HConfig.GetTH2D(Name+"_Mu1_PtEta","Muon 1 Pt vs Eta",50,0,25,50,-2.5,2.5,"Pt (GeV)","#eta");
   Muon2PtEta=HConfig.GetTH2D(Name+"_Mu2_PtEta","Muon 2 Pt vs Eta",50,0,25,50,-2.5,2.5,"Pt (GeV)","#eta");
   Muon3PtEta=HConfig.GetTH2D(Name+"_Mu3_PtEta","Muon 3 Pt vs Eta",50,0,25,50,-2.5,2.5,"Pt (GeV)","#eta");

   Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
   Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
   Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");

   Muon3isGlob =HConfig.GetTH1D(Name+"_Muon3isGlob","Muon3isGlob",2,-0.5,1.5,"  #mu_{3} is global muon","Events");

   Muon1isTrack =HConfig.GetTH1D(Name+"_Muon1isTrack","Muon1isTrack",2,-0.5,1.5,"  #mu_{1} is tracker muon","Events");
   Muon2isTrack =HConfig.GetTH1D(Name+"_Muon2isTrack","Muon2isTrack",2,-0.5,1.5,"  #mu_{2} is tracker muon","Events");
   Muon3isTrack =HConfig.GetTH1D(Name+"_Muon3isTrack","Muon3isTrack",2,-0.5,1.5,"  #mu_{3} is tracker muon","Events");

   Muon1kink =HConfig.GetTH1D(Name+"_Muon1kink","Muon1kink",50,0.,50,"  #mu_{1} kink","Events");
   Muon2kink =HConfig.GetTH1D(Name+"_Muon2kink","Muon2kink",50,0.,50,"  #mu_{2} kink","Events");
   Muon3kink =HConfig.GetTH1D(Name+"_Muon3kink","Muon3kink",50,0.,50,"  #mu_{3} kink","Events");

   Muon1InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon1InOutTrackMatch","Muon1InOutTrackMatch",50,0.,10,"  #mu_{1} inner and outer tracker match","Events");
   Muon2InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon2InOutTrackMatch","Muon2InOutTrackMatch",50,0.,10,"  #mu_{2} inner and outer tracker match","Events");
   Muon3InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon3InOutTrackMatch","Muon3InOutTrackMatch",50,0.,10,"  #mu_{3} inner and outer tracker match","Events");

   Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",26,-2.6,2.6,"#mu_{1}  rapidity","Events");
   Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",26,-2.6,2.6,"#mu_{2}  rapidity","Events");
   Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",26,-2.6,2.6,"#mu_{3}  rapidity","Events");

   Muon1StandardSelector=HConfig.GetTH1D(Name+"_Muon1StandardSelector","Muon1StandardSelector",23,-0.5,22.5,"#mu_{1} standard selector; bin 0 - no ID","Events");
   Muon2StandardSelector=HConfig.GetTH1D(Name+"_Muon2StandardSelector","Muon2StandardSelector",23,-0.5,22.5,"#mu_{2} standard selector; bin 0 - no ID","Events");
   Muon3StandardSelector=HConfig.GetTH1D(Name+"_Muon3StandardSelector","Muon3StandardSelector",23,-0.5,22.5,"#mu_{3} standard selector; bin 0 - no ID","Events");
   Muon1PtResolution=HConfig.GetTH1D(Name+"_Muon1PtResolution","Muon1PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{1})  (reco - mc)/mc ","Events");
   Muon2PtResolution=HConfig.GetTH1D(Name+"_Muon2PtResolution","Muon2PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{2})  (reco - mc)/mc ","Events");
   Muon3PtResolution=HConfig.GetTH1D(Name+"_Muon3PtResolution","Muon3PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{3})  (reco - mc)/mc  ","Events");

   Muon1EtaResolution=HConfig.GetTH1D(Name+"_Muon1EtaResolution","Muon1EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc  ","Events");
   Muon2EtaResolution=HConfig.GetTH1D(Name+"_Muon2EtaResolution","Muon2EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
   Muon3EtaResolution=HConfig.GetTH1D(Name+"_Muon3EtaResolution","Muon3EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
   TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"  #tau p_{T}, GeV","Events");
   TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"  #tau |p|, GeV","Events");
   TauMass =HConfig.GetTH1D(Name+"_TauMass","#tau lepton mass",25,1.65,1.9,"  M_{#tau} , GeV","Events");
   TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
   TauMassRefit =HConfig.GetTH1D(Name+"_TauMassRefit","Refit #tau lepton mass",40,1.5,1.9,"KF refit  M_{#tau} , GeV","Events");
   TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

   VertexDCA12=HConfig.GetTH1D(Name+"_VertexDCA12","VertexDCA12",40,0,0.15,"dca (#mu_{1}#mu_{2})","Events");
   VertexDCA23=HConfig.GetTH1D(Name+"_VertexDCA23","VertexDCA23",40,0,0.15,"dca (#mu_{1}trk)","Events");
   VertexDCA31=HConfig.GetTH1D(Name+"_VertexDCA31","VertexDCA31",40,0,0.15,"dca (#mu_{2}trk)","Events");
   VertexChi2AF=HConfig.GetTH1D(Name+"_VertexChi2AF","VertexChi2AF",50,0,10,"AF vertex #chi^{2}","Events");

   VertexSignalKFRefittedMu1P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1P","VertexSignalKFRefittedMu1P",50,0,20,"KF refitted #mu_{1} track p (GeV)","Events");
   VertexSignalKFRefittedMu1Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Pt","VertexSignalKFRefittedMu1P",50,0,20,"KF refitted #mu_{1} track p_{T} (GeV)","Events");
   VertexSignalKFRefittedMu1Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Eta","VertexSignalKFRefittedMu1Eta",25,-2.5,2.5,"KF refitted #mu_{1} track #eta","Events");
   VertexSignalKFRefittedMu1Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu1Phi","VertexSignalKFRefittedMu1Phi",32,-3.2,3.2,"KF refitted #mu_{1} track #phi","Events");

   VertexSignalKFRefittedMu2P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2P","VertexSignalKFRefittedMu2P",50,0,20,"KF refitted #mu_{2} track p (GeV)","Events");
   VertexSignalKFRefittedMu2Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Pt","VertexSignalKFRefittedMu2Pt",50,0,20,"KF refitted #mu_{2} track p_{T} (GeV)","Events");
   VertexSignalKFRefittedMu2Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Eta","VertexSignalKFRefittedMu2Eta",25,-2.5,2.5,"KF refitted #mu_{2} track #eta","Events");
   VertexSignalKFRefittedMu2Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu2Phi","VertexSignalKFRefittedMu2Phi",32,-3.2,3.2,"KF refitted #mu_{2} track #phi","Events");

   VertexSignalKFRefittedMu3P=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3P","VertexSignalKFRefittedMu3P",50,0,20,"KF refitted Mu3 p (GeV)","Events");
   VertexSignalKFRefittedMu3Pt=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Pt","VertexSignalKFRefittedMu3Pt",50,0,20,"KF refitted Mu3 p_{T} (GeV)","Events");
   VertexSignalKFRefittedMu3Eta=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Eta","VertexSignalKFRefittedMu3Eta",25,-2.5,2.5,"KF refitted Mu3 #eta","Events");
   VertexSignalKFRefittedMu3Phi=HConfig.GetTH1D(Name+"_VertexSignalKFRefittedMu3Phi","VertexSignalKFRefittedMu3Phi",32,-3.2,3.2,"KF refitted Mu3 #phi","Events");

   VertexMu1D0Reco=HConfig.GetTH1D(Name+"_VertexMu1D0Reco","VertexMu1D0Reco",40,0,0.05,"#mu_{1} d0 vertex","Events");
   VertexMu1D0SigReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigReco","VertexMu1D0SigReco",40,0,5,"#mu_{1} d0 vertex significance","Events");
   VertexMu2D0Reco=HConfig.GetTH1D(Name+"_VertexMu2D0Reco","VertexMu2D0Reco",40,0,0.05,"#mu_{2} d0 vertex","Events");
   VertexMu2D0SigReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigReco","VertexMu2D0SigReco",40,0,5,"#mu_{2} d0 vertex significance","Events");
   VertexMu3D0Reco=HConfig.GetTH1D(Name+"_VertexMu3D0Reco","VertexMu3D0Reco",40,0,0.05,"#mu_{3} d0 vertex","Events");
   VertexMu3D0SigReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigReco","VertexMu3D0SigReco",50,0,5,"#mu_{3} d0 vertex significance","Events");

   Vertex2DDisplacement=HConfig.GetTH1D(Name+"_Vertex2DDisplacement","Vertex2DDisplacement",100,0,1.0,"vertex 2d displacement","Events");
   Vertex3DDisplacement=HConfig.GetTH1D(Name+"_Vertex3DDisplacement","Vertex3DDisplacement",100,0,1.0,"vertex 3d displacement","Events");

   VertexPairQuality=HConfig.GetTH1D(Name+"_VertexPairQuality","VertexPairQuality",10,0,10,"vertex pair quality","Events");
   VertexPairfitStatus=HConfig.GetTH1D(Name+"_VertexPairfitStatus","VertexPairfitStatus",2,0,2,"vertex pair fit status","Events");

   VertexSignalKFChi2=HConfig.GetTH1D(Name+"_VertexSignalKFChi2","VertexSignalKFChi2",40,0,20,"vertex KF #chi^{2}","Events");
   /*
      VertexSignalAFX=HConfig.GetTH1D(Name+"_VertexSignalAFX","VertexSignalAFX",100,0,10,"AF vertex x (cm)","Events");
      VertexSignalAFY=HConfig.GetTH1D(Name+"_VertexSignalAFY","VertexSignalAFY",100,0,10,"AF vertex y (cm)","Events");
      VertexSignalAFZ=HConfig.GetTH1D(Name+"_VertexSignalAFZ","VertexSignalAFZ",100,0,10,"AF vertex z (cm)","Events");
      VertexSignalAFChi2=HConfig.GetTH1D(Name+"_VertexSignalChi2","VertexSignalChi2",100,0,10,"AF vertex #chi^{2}","Events");
      VertexSignalAFNdf=HConfig.GetTH1D(Name+"_VertexSignalAFNdf","VertexSignalAFNdf",10,0,10,"AF vertex ndf","Events");
      */

   VertexMatchedPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexX","VertexMatchedPrimaryVertexX",50,0,0.15,"vertex matched pv x (cm)","Events");
   VertexMatchedPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexY","VertexMatchedPrimaryVertexY",50,0,0.15,"vertex matched pv y (cm)","Events");
   VertexMatchedPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedPrimaryVertexZ","VertexMatchedPrimaryVertexZ",50,0,0.15,"vertex matched pv z (cm)","Events");
   VertexRefitPVisValid=HConfig.GetTH1D(Name+"_VertexRefitPVisValid","VertexRefitPVisValid",2,0,2,"vertex refit pv is valid","Events");
   VertexMatchedRefitPrimaryVertexX=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexX","VertexMatchedRefitPrimaryVertexX",50,0,0.15,"vertex matched refit pv x (cm)","Events");
   VertexMatchedRefitPrimaryVertexY=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexY","VertexMatchedRefitPrimaryVertexY",50,0,0.15,"vertex matched refit pv y (cm)","Events");
   VertexMatchedRefitPrimaryVertexZ=HConfig.GetTH1D(Name+"_VertexMatchedRefitPrimaryVertexZ","VertexMatchedRefitPrimaryVertexZ",50,0,0.15,"vertex matched refit pv z (cm)","Events");

   Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
   Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
   Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");

   MuPair1_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair1_vertex_chi2","MuPair1_vertex_chi2",50,0,5,"KF  #chi^{2} of first #mu pair","Events");
   MuPair2_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair2_vertex_chi2","MuPair2_vertex_chi2",50,0,5,"KF  #chi^{2} of second #mu pair","Events");
   MuPair3_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair3_vertex_chi2","MuPair3_vertex_chi2",50,0,5,"KF  #chi^{2} of third #mu pair","Events");

   TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",15,0,0.03,"trigger match dR 1","Events");
   TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",15,0,0.03,"trigger match dR 2","Events");
   TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",15,0,0.03,"trigger match dR 3","Events");

   dR12 =HConfig.GetTH1D(Name+"_dR12","dR12",25,0,0.5,"dR(#mu_{1}#mu_{2})","Events");
   dR23 =HConfig.GetTH1D(Name+"_dR23","dR23",25,0,0.5,"dR(#mu_{2}#mu_{3})","Events");
   dR31 = HConfig.GetTH1D(Name+"_dR31","dR31",25,0,1.5,"dR(#mu_{3}#mu_{1})","Events");
   dR1Tau = HConfig.GetTH1D(Name+"_dR1Tau","dR1Tau",15,0,0.3,"dR(#mu_{1}#tau)","Events");
   dR2Tau = HConfig.GetTH1D(Name+"_dR2Tau","dR2Tau",15,0,0.3,"dR(#mu_{2}#tau)","Events");
   dR3Tau = HConfig.GetTH1D(Name+"_dR3Tau","dR3Tau",15,0,0.3,"dR(#mu_{3}#tau)","Events");

   Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","Isolation_NTracks",10,-0.5,9.5,"N tracks","Events");
   Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","Isolation_RelPt",20,0,1,"relative p_{T}","Events");
   Isolation05_RelPt=HConfig.GetTH1D(Name+"_Isolation05_RelPt","Isolation05_RelPt",8,0,0.8,"relative  rel p_{T} in 0.5 cone","Events");
   Isolation05_NTracks=HConfig.GetTH1D(Name+"_Isolation05_NTracks","Isolation05_NTracks",10,-0.5,9.5,"N tracks in 0.5 cone","Events");
   Isolation05_MinDist=HConfig.GetTH1D(Name+"_Isolation05_MinDist","Isolation05_MinDist",10,0,0.5,"Iso05 MinDist","Events");

   Isolation_Mu1RelPt=HConfig.GetTH1D(Name+"Isolation_Mu1RelPt","Isolation_Mu1RelPt",10,0,1,"relPt (>1 GeV) in #mu_{1} iso (dR=0.3)","Events");
   Isolation_Mu2RelPt=HConfig.GetTH1D(Name+"Isolation_Mu2RelPt","Isolation_Mu2RelPt",10,0,1,"relPt (>1 GeV) in #mu_{2} iso (dR=0.3)","Events");
   Isolation_Mu3RelPt=HConfig.GetTH1D(Name+"Isolation_Mu3RelPt","Isolation_Mu3RelPt",10,0,1,"relPt (>1 GeV) in #mu_{3} iso (dR=0.3)","Events");

   // Muon selection in BDT input variables
   MuonVsMinSC = HConfig.GetTH1D(Name+"_MuonVsMinSC","min SC vs Muon",3,0.0,3.0,"Muon (Decreasing Pt)","Events");
   MuonVsMaxCLP = HConfig.GetTH1D(Name+"_MuonvsMaxCLP","max cLP vs Kink",3,0.0,3.0,"Muon (Decreasing Pt)","Events");
   MuonVsMinP = HConfig.GetTH1D(Name+"_MuonVsMinP","min P vs Kink",3,0.0,3.0,"Muon (Decreasing Pt)","Events");
   MuonVstKink = HConfig.GetTH1D(Name+"_MuonVstKink","tKink vs Muon",3,0.0,3.0,"Muon (Decreasing Pt)","Events");
   MuonVsD0Sig = HConfig.GetTH1D(Name+"_MuonVsD0Sig","MinD0Sig vs Muon",3,0.0,3.0,"Muon (Decreasing Pt)","Events");
   MuonVsRelPt = HConfig.GetTH1D(Name+"_MuonVsRelPt","RelPt vs Muon",3,0.0,3.0,"Muon (Decreasing Pt)","Events");

   MuMuMass_OS1 = HConfig.GetTH1D(Name+"_MuMuMass_OS1","Invariant Mass OS Mu pair 1",20,0,2.0);
   MuMuMass_OS2 = HConfig.GetTH1D(Name+"_MuMuMass_OS2","Invariant Mass OS Mu pair 2",20,0,2.0);

   MinMatchedStations=HConfig.GetTH1D(Name+"_MinMatchedStations","MinMatchedStations",10,-0.5,9.5,"min matched stations","");
   MaxMatchedStations=HConfig.GetTH1D(Name+"_MaxMatchedStations","MaxMatchedStations",10,-0.5,9.5,"max matched stations ","");
   Mu1MatchedStations=HConfig.GetTH1D(Name+"_Mu1MatchedStations","Mu1MatchedStations",10,-0.5,9.5,"#mu_{1} matched stations","");
   Mu2MatchedStations=HConfig.GetTH1D(Name+"_Mu2MatchedStations","Mu2MatchedStations",10,-0.5,9.5,"#mu_{2} matched stations","");
   Mu3MatchedStations=HConfig.GetTH1D(Name+"_Mu3MatchedStations","Mu3MatchedStations",10,-0.5,9.5,"#mu_{3} matched stations","");

   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  FillMVATree::Store_ExtraDist(){ 

   Extradist1d.push_back(&L1Seed);
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
   Extradist2d.push_back(&L1TriggerMatch);
   Extradist2d.push_back(&Test);
   Extradist2d.push_back(&TrackCorrelation);
   Extradist1d.push_back(&EventMassResolution_PtEtaPhi);

   Extradist1d.push_back(&MuonVsMinSC);
   Extradist1d.push_back(&MuonVsMaxCLP);
   Extradist1d.push_back(&MuonVsMinP);

   Extradist1d.push_back(&MuMuMass_OS1);
   Extradist1d.push_back(&MuMuMass_OS2);

   // Muon Kinematic plots
   Extradist1d.push_back(&Muon1Pt);
   Extradist1d.push_back(&Muon2Pt);
   Extradist1d.push_back(&Muon3Pt);

   Extradist1d.push_back(&Muon1Eta);
   Extradist1d.push_back(&Muon2Eta);
   Extradist1d.push_back(&Muon3Eta);

   Extradist1d.push_back(&Muon1StandardSelector);
   Extradist1d.push_back(&Muon2StandardSelector);
   Extradist1d.push_back(&Muon3StandardSelector);

   Extradist1d.push_back(&Muon3isGlob);

   Extradist1d.push_back(&Muon1isTrack);
   Extradist1d.push_back(&Muon2isTrack);
   Extradist1d.push_back(&Muon3isTrack);

   Extradist1d.push_back(&TauPt);
   Extradist1d.push_back(&TauP);
   Extradist1d.push_back(&TauMass);
   Extradist1d.push_back(&TauMassRefit);
   Extradist1d.push_back(&TauMassResolution);
   Extradist1d.push_back(&TauMassResolutionRefit);

   Extradist1d.push_back(&Muon1kink);
   Extradist1d.push_back(&Muon2kink);
   Extradist1d.push_back(&Muon3kink);

   Extradist1d.push_back(&Muon1InOutTrackMatch);
   Extradist1d.push_back(&Muon2InOutTrackMatch);
   Extradist1d.push_back(&Muon3InOutTrackMatch);

   Extradist1d.push_back(&Muon1PtResolution);
   Extradist1d.push_back(&Muon2PtResolution);
   Extradist1d.push_back(&Muon3PtResolution);

   Extradist1d.push_back(&Muon1EtaResolution);
   Extradist1d.push_back(&Muon2EtaResolution);
   Extradist1d.push_back(&Muon3EtaResolution);

   Extradist1d.push_back(&Muon1DRToTruth);
   Extradist1d.push_back(&Muon2DRToTruth);
   Extradist1d.push_back(&Muon3DRToTruth);

   Extradist1d.push_back(&MuPair1_vertex_chi2);
   Extradist1d.push_back(&MuPair2_vertex_chi2);
   Extradist1d.push_back(&MuPair3_vertex_chi2);

   Extradist1d.push_back(&TriggerMatchdR1);
   Extradist1d.push_back(&TriggerMatchdR2);
   Extradist1d.push_back(&TriggerMatchdR3);

   Extradist1d.push_back(&dR12);
   Extradist1d.push_back(&dR23);
   Extradist1d.push_back(&dR31);
   Extradist1d.push_back(&dR1Tau);
   Extradist1d.push_back(&dR2Tau);
   Extradist1d.push_back(&dR3Tau);

   Extradist1d.push_back(&Isolation_NTracks);
   Extradist1d.push_back(&Isolation_RelPt);
   Extradist1d.push_back(&Isolation05_RelPt);
   Extradist1d.push_back(&Isolation05_NTracks);
   Extradist1d.push_back(&Isolation05_MinDist);

   Extradist1d.push_back(&VertexChi2AF);
   Extradist1d.push_back(&VertexDCA12);
   Extradist1d.push_back(&VertexDCA23);
   Extradist1d.push_back(&VertexDCA31);
   Extradist1d.push_back(&VertexSignalKFRefittedMu1P);
   Extradist1d.push_back(&VertexSignalKFRefittedMu1Pt);
   Extradist1d.push_back(&VertexSignalKFRefittedMu1Eta);
   Extradist1d.push_back(&VertexSignalKFRefittedMu1Phi);
   Extradist1d.push_back(&VertexSignalKFRefittedMu2P);
   Extradist1d.push_back(&VertexSignalKFRefittedMu2Pt);
   Extradist1d.push_back(&VertexSignalKFRefittedMu2Eta);
   Extradist1d.push_back(&VertexSignalKFRefittedMu2Phi);
   Extradist1d.push_back(&VertexSignalKFRefittedMu3P);
   Extradist1d.push_back(&VertexSignalKFRefittedMu3Pt);
   Extradist1d.push_back(&VertexSignalKFRefittedMu3Eta);
   Extradist1d.push_back(&VertexSignalKFRefittedMu3Phi);

   Extradist1d.push_back(&VertexMu1D0Reco);
   Extradist1d.push_back(&VertexMu1D0SigReco);
   Extradist1d.push_back(&VertexMu2D0Reco);
   Extradist1d.push_back(&VertexMu2D0SigReco);
   Extradist1d.push_back(&VertexMu3D0Reco);
   Extradist1d.push_back(&VertexMu3D0SigReco);
   Extradist1d.push_back(&Vertex2DDisplacement);
   Extradist1d.push_back(&Vertex3DDisplacement);
   Extradist1d.push_back(&VertexPairQuality);
   Extradist1d.push_back(&VertexPairfitStatus);
   Extradist1d.push_back(&VertexSignalKFChi2);

   /*
      Extradist1d.push_back(&VertexSignalAFX);
      Extradist1d.push_back(&VertexSignalAFY);
      Extradist1d.push_back(&VertexSignalAFZ);
      Extradist1d.push_back(&VertexSignalAFChi2);
      Extradist1d.push_back(&VertexSignalAFNdf);
      */

   Extradist1d.push_back(&VertexMatchedPrimaryVertexX);
   Extradist1d.push_back(&VertexMatchedPrimaryVertexY);
   Extradist1d.push_back(&VertexMatchedPrimaryVertexZ);
   Extradist1d.push_back(&VertexRefitPVisValid);
   Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexX);
   Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexY);
   Extradist1d.push_back(&VertexMatchedRefitPrimaryVertexZ);

   Extradist1d.push_back(&Isolation_Mu1RelPt);
   Extradist1d.push_back(&Isolation_Mu2RelPt);
   Extradist1d.push_back(&Isolation_Mu3RelPt); 

   Extradist2d.push_back(&Muon1PtEta);
   Extradist2d.push_back(&Muon2PtEta);
   Extradist2d.push_back(&Muon3PtEta);

   Extradist1d.push_back(&NMuons);
   Extradist1d.push_back(&NTwoMuonsTrack);
   Extradist1d.push_back(&NTracks);

   Extradist1d.push_back(&MinMatchedStations);
   Extradist1d.push_back(&MaxMatchedStations);
   Extradist1d.push_back(&Mu1MatchedStations);
   Extradist1d.push_back(&Mu2MatchedStations);
   Extradist1d.push_back(&Mu3MatchedStations);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
   ///t/////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  FillMVATree::doEvent(){ 

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
   bool DoubleMuFired(false);
   bool TripleMuFired(false);

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

   if (DoubleMuFired && !TripleMuFired) l1seed = 1;
   if (DoubleMuFired && TripleMuFired) l1seed = 2;
   if (!DoubleMuFired && TripleMuFired) l1seed = 3;

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
               Ntp->Muon_isGlobalMuon(Muon_index_3) );
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
            Ntp->Muon_isGlobalMuon(Muon_index_3));
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
   double wobs=1;
   double w; 

   if(!Ntp->isData()) { w=1; } //  No weights to data
   else{w=1;}

   bool status=AnalysisCuts(t,w,wobs);

   if(status){

      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      std::vector<unsigned int> EtaSortedIndices;

      EtaSortedIndices.push_back(Muon_Eta_index_1);
      EtaSortedIndices.push_back(Muon_Eta_index_2);
      EtaSortedIndices.push_back(Muon_Eta_index_3);

      EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);

      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);

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

      MuMuMass_OS1.at(t).Fill((MuonOS+MuonSS1).M());
      MuMuMass_OS2.at(t).Fill((MuonOS+MuonSS2).M());

      TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);
      TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);

      //    TauMass.at(t).Fill(TauLV.M(),1);
      //    TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    
      double tauMassRes = Ntp->TauMassResolution(EtaSortedIndices,1,false);
      float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,false),
            Ntp->Vertex_d0sig_reco(final_idx,1,false),
            Ntp->Vertex_d0sig_reco(final_idx,2,false)});
      float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0,false),
            Ntp->Vertex_d0sig_reco(final_idx,1,false),
            Ntp->Vertex_d0sig_reco(final_idx,2,false)});
      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,false),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));

      double maxMuondR = std::max({Muon1LV.DeltaR(TauLV), Muon2LV.DeltaR(TauLV), Muon3LV.DeltaR(TauLV)});
      double minMuonPt = std::min({Muon1LV.Pt(), Muon2LV.Pt(), Muon3LV.Pt()});

      // Isolation algorithm
      for (int it=0; it<Ntp->NIsolationTrack(final_idx, false); it++){
         double dxy_track = Ntp->IsolationTrack_dxySV(final_idx, it, false);
         double dz_track = Ntp->IsolationTrack_dzSV(final_idx, it, false);
         TLorentzVector TrackLV = Ntp->IsolationTrack_p4(final_idx, it, false);
         double dca_fv = TMath::Sqrt(pow(dxy_track, 2)+ pow(dz_track, 2));

         double dr_tau = TauLV.DeltaR(TrackLV); 
         double dr_mu1 = Muon1LV.DeltaR(TrackLV);
         double dr_mu2 = Muon2LV.DeltaR(TrackLV);
         double dr_mu3 = Muon3LV.DeltaR(TrackLV);

         // Isolation 1
         if ( dca_fv<0.5 && TrackLV.Pt()<0.33*minMuonPt && dr_tau<3*maxMuondR ){
            sumPtTracks_tau += TrackLV.Pt();
            nTracks_tau++;
            if (dca_fv < mindca_tau) mindca_tau = dca_fv;
         }

         // Isolation 2
         if (TrackLV.Pt()<1.0) {
            continue;
         }

         if (dca_fv < mindca_iso) mindca_iso = dca_fv;

         // Isolation 3 (within dR = 0.5 of tau)
         if (dr_tau<0.5 && dca_fv<0.5){
            sumPtTracks_iso05 += TrackLV.Pt();
            nTracks_iso05++;
            if(dca_fv<mindca_iso05) mindca_iso05 = dca_fv;
         }

         // Isolation 4 (Muon isolation)
         if (dr_mu1 < 0.3 && Ntp->IsolationTrack_DocaMu1(final_idx, it, false) < 0.1 ) sumPtTracks_mu1 += TrackLV.Pt();
         if (dr_mu2 < 0.3 && Ntp->IsolationTrack_DocaMu2(final_idx, it, false) < 0.1 ) sumPtTracks_mu2 += TrackLV.Pt();
         if (dr_mu3 < 0.3 && Ntp->IsolationTrack_DocaMu3(final_idx, it, false) < 0.1 ) sumPtTracks_mu3 += TrackLV.Pt();
      }
      // Relative Pt calculation
      double mu1_relPt = sumPtTracks_mu1/Muon1LV.Pt();
      double mu2_relPt = sumPtTracks_mu2/Muon2LV.Pt();
      double mu3_relPt = sumPtTracks_mu3/Muon3LV.Pt();
      double relPt_iso05 = sumPtTracks_iso05/TauLV.Pt();

      // Categorization variables
      var_Eta_Tau = TauLV.Eta();
      var_tauMassRes = tauMassRes;

      // Muon/Tau kinematic variables
      var_pmin = std::min(Ntp->Muon_P4(Muon_index_1).P(), std::min(Ntp->Muon_P4(Muon_index_2).P(), Ntp->Muon_P4(Muon_index_3).P()));
      var_RelPt_Mu1Tau = Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt();
      var_MuMu_mindR = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
      var_MuTau_maxdR = std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});
      var_minMatchedStations = std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});

      // Muon ID variables
      var_max_cLP = std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1), std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2), Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)));
      var_max_tKink = std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_1), std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_2), Ntp->Muon_combinedQuality_trkKink(Muon_index_3)));
      var_segCompMuMin  = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      var_MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      var_sumMuTrkKinkChi2= (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));

      // vertex variables
      var_maxdca = std::max({Ntp->Vertex_DCA12(final_idx, false),Ntp->Vertex_DCA23(final_idx, false),Ntp->Vertex_DCA31(final_idx, false)});
      var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(final_idx, 0, false), Ntp->Vertex_pair_quality(final_idx, 1, false), Ntp->Vertex_pair_quality(final_idx, 2, false)});
      var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
      var_flightLenSig = sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx, false),Ntp->Vertex_PrimaryVertex_Covariance(final_idx, false),
               Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx, false)));
      var_vertexKFChi2 =Ntp->Vertex_signal_KF_Chi2(final_idx, false);
      // Isolation variables
      var_MinD0Significance = MinD0Significance;
      var_MaxD0Significance = MaxD0Significance;
      var_mindca_iso = mindca_iso;
      var_trk_relPt = std::max({mu1_relPt, mu2_relPt, mu3_relPt});

      // Spectator variables
      var_tauMass = TauLV.M();
      var_tauMassRefit = TauRefitLV.M();
      if (id==1) MC=0;
      else  MC=1;

      if (tauMassRes<tauMassResCutLow) category = 1;
      if (tauMassRes>=tauMassResCutLow && tauMassRes<tauMassResCutHigh) category = 2;
      if (tauMassRes>=tauMassResCutHigh) category = 3;

      //    if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
      TMVA_Tree->Fill();
      // ----- Fill the histograms ----- 
      if (DoubleMuFired && !TripleMuFired) L1Seed.at(t).Fill(0.0,w);
      else if (!DoubleMuFired && TripleMuFired) L1Seed.at(t).Fill(1.0,w);
      else if (DoubleMuFired && TripleMuFired) L1Seed.at(t).Fill(2.0,w);
      VertexDCAMax.at(t).Fill(var_maxdca,w);
      SVPVTauDirAngle.at(t).Fill(var_svpvTauAngle,w);
      FLSignificance.at(t).Fill(var_flightLenSig,w);
      VertexChi2KF.at(t).Fill(var_vertexKFChi2,w);
      MuonglbkinkSum.at(t).Fill(var_sumMuTrkKinkChi2,w);
      Muon_segmentCompatibility_min.at(t).Fill(var_segCompMuMin,w);
      Muon_HCALCompatibility_min.at(t).Fill(var_MinMIPLikelihood,w);
      minMudR.at(t).Fill(std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)}));
      Mu1TauPTRatio.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt());
      dRMaxMuTau.at(t).Fill(std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)}));
      MuPair_vertex_chi2_min.at(t).Fill(std::min({Ntp->Vertex_pair_quality(final_idx, 0, false), Ntp->Vertex_pair_quality(final_idx, 1, false), Ntp->Vertex_pair_quality(final_idx, 2, false)}));
      TauEta.at(t).Fill(TauLV.Eta(),w);
      VertexMuMaxD0SigReco.at(t).Fill(MaxD0Significance,w);

      std::vector<double> vecIndices, vecSC, vecCLP, vecP, vecD0Sig, vecRelPt, vecKink;
      vecIndices.push_back(Muon_index_1);
      vecIndices.push_back(Muon_index_2);
      vecIndices.push_back(Muon_index_3);

      for (size_t k=0; k<3; ++k){
         int tmp_index = vecIndices.at(k);
         vecSC.push_back(Ntp->Muon_segmentCompatibility(tmp_index));
         vecCLP.push_back(Ntp->Muon_combinedQuality_chi2LocalPosition(tmp_index));
         vecP.push_back(Ntp->Muon_P4(tmp_index).P());
         vecKink.push_back(Ntp->Muon_combinedQuality_trkKink(tmp_index));
         for (size_t imu=0; imu<3; imu++){
            if (tmp_index==Ntp->ThreeMuonIndices(final_idx).at(imu)) vecD0Sig.push_back(Ntp->Vertex_d0sig_reco(final_idx,imu,false));
         }
      }
      vecRelPt.push_back(mu1_relPt);
      vecRelPt.push_back(mu2_relPt);
      vecRelPt.push_back(mu3_relPt);

      MuonVsMinSC.at(t).Fill(minQuantityIndex(vecSC));
      MuonVsMaxCLP.at(t).Fill(maxQuantityIndex(vecCLP));
      MuonVsMinP.at(t).Fill(minQuantityIndex(vecP));
      MuonVstKink.at(t).Fill(maxQuantityIndex(vecKink));
      MuonVsD0Sig.at(t).Fill(minQuantityIndex(vecD0Sig));
      MuonVsRelPt.at(t).Fill(maxQuantityIndex(vecRelPt));

      // Count number of muons, tracks and dsphipi candidates
      NMuons.at(t).Fill(Ntp->NMuons());
      NTracks.at(t).Fill(Ntp->NTracks());
      NTwoMuonsTrack.at(t).Fill(Ntp->NTwoMuonsTrack());

      // Muon kinematic plots
      Muon1PtEta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),Ntp->Muon_P4(Muon_index_1).Eta());
      Muon2PtEta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),Ntp->Muon_P4(Muon_index_2).Eta());
      Muon3PtEta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),Ntp->Muon_P4(Muon_index_3).Eta());

      Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
      Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
      Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);

      Muon1Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Eta(),1);
      Muon2Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Eta(),1);
      Muon3Eta.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Eta(),1);

      dR12.at(t).Fill(Muon1LV.DeltaR(Muon2LV),1);
      dR23.at(t).Fill(Muon2LV.DeltaR(Muon3LV),1);
      dR31.at(t).Fill(Muon1LV.DeltaR(Muon3LV),1);
      dR1Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
      dR2Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);
      dR3Tau.at(t).Fill(Muon1LV.DeltaR(TauLV),1);

      TauPt.at(t).Fill(TauLV.Pt(),1);
      TauP.at(t).Fill(TauLV.P(),1);
      TauMass.at(t).Fill(TauLV.M(),1);
      TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    

      Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);

      Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
      Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
      Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);

      Muon1kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1),1);
      Muon2kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_2),1);
      Muon3kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_3),1);

      Muon1InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),1);
      Muon2InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),1);
      Muon3InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),1);

      MuPair1_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,0),1);
      MuPair2_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,1),1);
      MuPair3_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,2),1);

      TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0),1);
      TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1),1);
      TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2),1);

      MinMatchedStations.at(t).Fill(var_minMatchedStations,w);
      MaxMatchedStations.at(t).Fill(std::max({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)}),w);
      Mu1MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_1),w);
      Mu2MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_2),w);
      Mu3MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_3),w);

      Isolation_NTracks.at(t).Fill(nTracks_tau);
      Isolation_RelPt.at(t).Fill(var_trk_relPt);
      Isolation_MinDist.at(t).Fill(mindca_iso);
      Isolation05_RelPt.at(t).Fill(relPt_iso05);
      Isolation05_NTracks.at(t).Fill(nTracks_iso05);
      Isolation05_MinDist.at(t).Fill(mindca_iso05);

      Isolation_Mu1RelPt.at(t).Fill(mu1_relPt);
      Isolation_Mu2RelPt.at(t).Fill(mu2_relPt);
      Isolation_Mu3RelPt.at(t).Fill(mu3_relPt);

      VertexChi2AF.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(final_idx),w);
      VertexDCA12.at(t).Fill(Ntp->Vertex_DCA12(final_idx),w);
      VertexDCA23.at(t).Fill(Ntp->Vertex_DCA23(final_idx),w);
      VertexDCA31.at(t).Fill(Ntp->Vertex_DCA31(final_idx),w);
      VertexSignalKFRefittedMu1P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false).P(),w);
      VertexSignalKFRefittedMu1Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false).Pt(),w);
      VertexSignalKFRefittedMu1Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false).Eta(),w);
      VertexSignalKFRefittedMu1Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false).Phi(),w);
      VertexSignalKFRefittedMu2P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false).P(),w);
      VertexSignalKFRefittedMu2Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false).Pt(),w);
      VertexSignalKFRefittedMu2Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false).Eta(),w);
      VertexSignalKFRefittedMu2Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false).Phi(),w);
      VertexSignalKFRefittedMu3P.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false).P(),w);
      VertexSignalKFRefittedMu3Pt.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false).Pt(),w);
      VertexSignalKFRefittedMu3Eta.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false).Eta(),w);
      VertexSignalKFRefittedMu3Phi.at(t).Fill(Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false).Phi(),w);
      VertexMu1D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,0,false),w);
      VertexMu1D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,0,false),w);
      VertexMu2D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,1,false),w);
      VertexMu2D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,1,false),w);
      VertexMu3D0Reco.at(t).Fill(Ntp->Vertex_d0_reco(final_idx,2,false),w);
      VertexMu3D0SigReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,2,false),w);
      Vertex2DDisplacement.at(t).Fill(Ntp->Vertex_2Ddisplacement(final_idx,0,false),w);
      Vertex3DDisplacement.at(t).Fill(Ntp->Vertex_3Ddisplacement(final_idx,0,false),w);
      VertexPairQuality.at(t).Fill(Ntp->Vertex_pair_quality(final_idx,0,false),w);
      VertexPairfitStatus.at(t).Fill(Ntp->Vertex_pairfit_status(final_idx,0,false),w);
      VertexSignalKFChi2.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,false),w);

      //---------------  Fill MC plots 
      if(id==40 || id == 60 || id ==90){
         if(Ntp->MCEventIsReconstructed()){

            TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0)));
            TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1)));
            TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2)));
            TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
            Muon1PtResolution.at(t).Fill((Muon1LV.Pt() - MCMuon1LV.Pt())/MCMuon1LV.Pt(), 1);
            Muon2PtResolution.at(t).Fill((Muon2LV.Pt() - MCMuon2LV.Pt())/MCMuon2LV.Pt(), 1);
            Muon3PtResolution.at(t).Fill((Muon3LV.Pt() - MCMuon3LV.Pt())/MCMuon3LV.Pt(), 1);

            Muon1EtaResolution.at(t).Fill((Muon1LV.Eta() - MCMuon1LV.Eta())/MCMuon1LV.Eta(), 1);
            Muon2EtaResolution.at(t).Fill((Muon2LV.Eta() - MCMuon2LV.Eta())/MCMuon2LV.Eta(), 1);
            Muon3EtaResolution.at(t).Fill((Muon3LV.Eta() - MCMuon3LV.Eta())/MCMuon3LV.Eta(), 1);

            TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
            TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

            Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
            Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
            Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
         }
      }

      //	}
   }
}

template <typename T>
int FillMVATree::minQuantityIndex(std::vector<T>& vec){
   if (vec.at(0)<=vec.at(1) && vec.at(0)<=vec.at(2)) return 0;
   if (vec.at(1)<=vec.at(2) && vec.at(1)<=vec.at(0)) return 1;
   if (vec.at(2)<=vec.at(0) && vec.at(2)<=vec.at(1)) return 2;
   return -1;
}

template <typename T>
int FillMVATree::maxQuantityIndex(std::vector<T>& vec){
   if (vec.at(0)>=vec.at(1) && vec.at(0)>=vec.at(2)) return 0;
   if (vec.at(1)>=vec.at(2) && vec.at(1)>=vec.at(0)) return 1;
   if (vec.at(2)>=vec.at(0) && vec.at(2)>=vec.at(1)) return 2;
   return -1;
}

void  FillMVATree::Finish(){
   /* 
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
   */
   file= new TFile("FillMVATreeInput.root","recreate");
   TMVA_Tree->SetDirectory(file);

   file->Write();
   file->Close();

   Selection::Finish();
}
