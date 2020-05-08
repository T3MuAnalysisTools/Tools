// T3MSelectionTree
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
   TString basedir = "";
   basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2017";
   PUWeightFile = new TFile(basedir+"/PUWeights_Run2017.root");
   puWeights = (TH1D*)PUWeightFile->Get("h1_weights");
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

   TMVA_Tree->Branch("var_MaxdeltaMuZ",&var_MaxdeltaMuZ);
   TMVA_Tree->Branch("var_MindeltaMuZ",&var_MindeltaMuZ);
   TMVA_Tree->Branch("var_maxMuonsDca",&var_maxMuonsDca);
   TMVA_Tree->Branch("var_minMuonsDca",&var_minMuonsDca);
   TMVA_Tree->Branch("var_nsv",&var_nsv);


   TMVA_Tree->Branch("var_VertexMu1D0SigPVReco",&var_VertexMu1D0SigPVReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigPVReco",&var_VertexMu2D0SigPVReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigPVReco",&var_VertexMu3D0SigPVReco);

   TMVA_Tree->Branch("var_MaxD0SigPV",&var_MaxD0SigPV);
   TMVA_Tree->Branch("var_MinD0SigPV",&var_MinD0SigPV);


   TMVA_Tree->Branch("var_VertexMu1D0SigBSReco",&var_VertexMu1D0SigBSReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigBSReco",&var_VertexMu2D0SigBSReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigBSReco",&var_VertexMu3D0SigBSReco);

   TMVA_Tree->Branch("var_MaxD0SigBS",&var_MaxD0SigBS);
   TMVA_Tree->Branch("var_MinD0SigBS",&var_MinD0SigBS);



   TMVA_Tree->Branch("var_VertexMu1D0SigSVReco",&var_VertexMu1D0SigSVReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigSVReco",&var_VertexMu2D0SigSVReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigSVReco",&var_VertexMu3D0SigSVReco);


   TMVA_Tree->Branch("var_MaxD0SigSV",&var_MaxD0SigSV);
   TMVA_Tree->Branch("var_MinD0SigSV",&var_MinD0SigSV);


   TMVA_Tree->Branch("var_MinMuon_chi2LocalPosition",&var_MinMuon_chi2LocalPosition);
   TMVA_Tree->Branch("var_MaxMuon_chi2LocalPosition",&var_MaxMuon_chi2LocalPosition);

   TMVA_Tree->Branch("var_Muon1_chi2LocalPosition",&var_Muon1_chi2LocalPosition);
   TMVA_Tree->Branch("var_Muon2_chi2LocalPosition",&var_Muon2_chi2LocalPosition);
   TMVA_Tree->Branch("var_Muon3_chi2LocalPosition",&var_Muon3_chi2LocalPosition);


   TMVA_Tree->Branch("var_MinMuon_chi2LocalMomentum",&var_MinMuon_chi2LocalMomentum);
   TMVA_Tree->Branch("var_MaxMuon_chi2LocalMomentum",&var_MaxMuon_chi2LocalMomentum);

   TMVA_Tree->Branch("var_Muon1_chi2LocalMomentum",&var_Muon1_chi2LocalMomentum);
   TMVA_Tree->Branch("var_Muon2_chi2LocalMomentum",&var_Muon2_chi2LocalMomentum);
   TMVA_Tree->Branch("var_Muon3_chi2LocalMomentum",&var_Muon3_chi2LocalMomentum);


   TMVA_Tree->Branch("var_MintrkKink",&var_MintrkKink);
   TMVA_Tree->Branch("var_MaxtrkKink",&var_MaxtrkKink);
   TMVA_Tree->Branch("var_MinglbKink",&var_MinglbKink);
   TMVA_Tree->Branch("var_MaxglbKink",&var_MaxglbKink);

   TMVA_Tree->Branch("var_MuonglbkinkSum",&var_MuonglbkinkSum);

   TMVA_Tree->Branch("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
   TMVA_Tree->Branch("var_MinVertexPairQuality",&var_MinVertexPairQuality);

   TMVA_Tree->Branch("var_Iso02",&var_Iso02);
   TMVA_Tree->Branch("var_Iso04",&var_Iso04);
   TMVA_Tree->Branch("var_Iso06",&var_Iso06);
   TMVA_Tree->Branch("var_Iso08",&var_Iso08);
   TMVA_Tree->Branch("var_Iso1",&var_Iso1);
   TMVA_Tree->Branch("var_Iso12",&var_Iso12);

   TMVA_Tree->Branch("var_Iso02Mu1",&var_Iso02Mu1);
   TMVA_Tree->Branch("var_Iso04Mu1",&var_Iso04Mu1);
   TMVA_Tree->Branch("var_Iso06Mu1",&var_Iso06Mu1);
   TMVA_Tree->Branch("var_Iso08Mu1",&var_Iso08Mu1);
   TMVA_Tree->Branch("var_Iso1Mu1",&var_Iso1Mu1);
   TMVA_Tree->Branch("var_Iso12Mu1",&var_Iso12Mu1);



   TMVA_Tree->Branch("var_Iso02Mu2",&var_Iso02Mu2);
   TMVA_Tree->Branch("var_Iso04Mu2",&var_Iso04Mu2);
   TMVA_Tree->Branch("var_Iso06Mu2",&var_Iso06Mu2);
   TMVA_Tree->Branch("var_Iso08Mu2",&var_Iso08Mu2);
   TMVA_Tree->Branch("var_Iso1Mu2",&var_Iso1Mu2);
   TMVA_Tree->Branch("var_Iso12Mu2",&var_Iso12Mu2);

   TMVA_Tree->Branch("var_Iso02Mu3",&var_Iso02Mu3);
   TMVA_Tree->Branch("var_Iso04Mu3",&var_Iso04Mu3);
   TMVA_Tree->Branch("var_Iso06Mu3",&var_Iso06Mu3);
   TMVA_Tree->Branch("var_Iso08Mu3",&var_Iso08Mu3);
   TMVA_Tree->Branch("var_Iso1Mu3",&var_Iso1Mu3);
   TMVA_Tree->Branch("var_Iso12Mu3",&var_Iso12Mu3);

   TMVA_Tree->Branch("var_Iso08MuMax",&var_Iso08MuMax);
   TMVA_Tree->Branch("var_Iso08MuMin",&var_Iso08MuMin);

   TMVA_Tree->Branch("var_NtracksClose",&var_NtracksClose);
   TMVA_Tree->Branch("var_Muon1Pt",&var_Muon1Pt);
   TMVA_Tree->Branch("var_Muon2Pt",&var_Muon2Pt);
   TMVA_Tree->Branch("var_Muon3Pt",&var_Muon3Pt);


   TMVA_Tree->Branch("var_MindcaTrackSV",&var_MindcaTrackSV);
   TMVA_Tree->Branch("var_Mu1TrackMass",&var_Mu1TrackMass);
   TMVA_Tree->Branch("var_Mu2TrackMass",&var_Mu2TrackMass);
   TMVA_Tree->Branch("var_Mu3TrackMass",&var_Mu3TrackMass);
   TMVA_Tree->Branch("var_dcaTrackPV",&var_dcaTrackPV);


   TMVA_Tree->Branch("var_MinMatchedStations",&var_MinMatchedStations);
   TMVA_Tree->Branch("var_MaxMatchedStations",&var_MaxMatchedStations);
   TMVA_Tree->Branch("var_Mu1MatchedStations",&var_Mu1MatchedStations);
   TMVA_Tree->Branch("var_Mu2MatchedStations",&var_Mu2MatchedStations);
   TMVA_Tree->Branch("var_Mu3MatchedStations",&var_Mu3MatchedStations);


   TMVA_Tree->Branch("var_MinMuon_numberOfChambers",&var_MinMuon_numberOfChambers);
   TMVA_Tree->Branch("var_MaxMuon_numberOfChambers",&var_MaxMuon_numberOfChambers);
   TMVA_Tree->Branch("var_Mu1Muon_numberOfChambers",&var_Mu1Muon_numberOfChambers);
   TMVA_Tree->Branch("var_Mu2Muon_numberOfChambers",&var_Mu2Muon_numberOfChambers);
   TMVA_Tree->Branch("var_Mu3Muon_numberOfChambers",&var_Mu3Muon_numberOfChambers);



   TMVA_Tree->Branch("var_MinMuon_numberOfMatches",&var_MinMuon_numberOfMatches);
   TMVA_Tree->Branch("var_MaxMuon_numberOfMatches",&var_MaxMuon_numberOfMatches);
   TMVA_Tree->Branch("var_Mu1Muon_numberOfMatches",&var_Mu1Muon_numberOfMatches);
   TMVA_Tree->Branch("var_Mu2Muon_numberOfMatches",&var_Mu2Muon_numberOfMatches);
   TMVA_Tree->Branch("var_Mu3Muon_numberOfMatches",&var_Mu3Muon_numberOfMatches);

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
      if(i==GenMatch)           cut.at(GenMatch)=0.03;
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


   Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
   Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
   Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");

   TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
   TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
   TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");

   VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
   MuonglbkinkSum  =HConfig.GetTH1D(Name+"_MuonglbkinkSum","MuonglbkinkSum",50,0.,50," #sum  #mu glb kink #chi^{2}","Events");
   FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",60,0,60,"PV - SV distance  significance","Events");
   FL=HConfig.GetTH1D(Name+"_FL","FL",60,0,1.5,"Flight length ,cm","Events");

   SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",60,0,0.25,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");
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



   EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
   EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

   VertexMu1D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigPVReco","VertexMu1D0SigPVReco",50,0,15,"#mu_{1} - PV transverse distance significance","Events");
   VertexMu2D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigPVReco","VertexMu2D0SigPVReco",50,0,15,"#mu_{2} - PV transverse distance significance","Events");
   VertexMu3D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigPVReco","VertexMu3D0SigPVReco",50,0,15,"#mu_{3} - PV transverse distance significance","Events");


   MaxD0SigPV=HConfig.GetTH1D(Name+"_MaxD0SigPV","MaxD0SigPV",50,0,15,"Max Transverse Impact significance w.r.t PV","");
   MinD0SigPV=HConfig.GetTH1D(Name+"_MinD0SigPV","MinD0SigPV",50,0,15,"Min Transverse Impact significance w.r.t PV","");


   VertexMu1D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigBSReco","VertexMu1D0SigBSReco",50,0,15,"#mu_{1} - BS transverse distance significance","Events");
   VertexMu2D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigBSReco","VertexMu2D0SigBSReco",50,0,15,"#mu_{2} - BS transverse distance significance","Events");
   VertexMu3D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigBSReco","VertexMu3D0SigBSReco",50,0,15,"#mu_{3} - BS transverse distance significance","Events");


   MaxD0SigBS=HConfig.GetTH1D(Name+"_MaxD0SigBS","MaxD0SigBS",50,0,15,"Max Transverse Impact significance w.r.t BS","");
   MinD0SigBS=HConfig.GetTH1D(Name+"_MinD0SigBS","MinD0SigBS",50,0,15,"Min Transverse Impact significance w.r.t BS","");


   VertexMu1D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigSVReco","VertexMu1D0SigSVReco",50,0,5,"#mu_{1} - SV transverse distance significance","Events");
   VertexMu2D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigSVReco","VertexMu2D0SigSVReco",50,0,5,"#mu_{2} - SV transverse distance significance","Events");
   VertexMu3D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigSVReco","VertexMu3D0SigSVReco",50,0,5,"#mu_{3} - SV transverse distance significance","Events");


   MaxD0SigSV=HConfig.GetTH1D(Name+"_MaxD0SigSV","MaxD0SigSV",50,0,5,"Max Transverse Impact significance w.r.t SV","");
   MinD0SigSV=HConfig.GetTH1D(Name+"_MinD0SigSV","MinD0SigSV",50,0,5,"Min Transverse Impact significance w.r.t SV","");

   MintrkKink= HConfig.GetTH1D(Name+"_MintrkKink","MintrkKink",30,0,30,"Min Tracker Kink","");
   MaxtrkKink= HConfig.GetTH1D(Name+"_MaxtrkKink","MaxtrkKink",30,0,30,"Max Tracker Kink","");
   MinglbKink= HConfig.GetTH1D(Name+"_MinglbKink","MinglbKink",30,0,30,"Min Global  Kink","");
   MaxglbKink= HConfig.GetTH1D(Name+"_MaxglbKink","MaxglbKink",30,0,30,"Max Global  Kink","");

   MinMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition","MinMuon_chi2LocalPosition",50,0,2,"Min Inner/Outer track #chi^{2}","");
   MaxMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MaxMuon_chi2LocalPosition","MaxMuon_chi2LocalPosition",50,0,30,"Max Inner/Outer track #chi^{2}","");

   Muon1_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon1_chi2LocalPosition","Muon1_chi2LocalPosition",50,0,30,"#mu_{1} Inner/Outer track #chi^{2}","");
   Muon2_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon2_chi2LocalPosition","Muon2_chi2LocalPosition",50,0,30,"#mu_{2} Inner/Outer track #chi^{2}","");
   Muon3_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon3_chi2LocalPosition","Muon3_chi2LocalPosition",50,0,30,"#mu_{3} Inner/Outer track #chi^{2}","");


   MinMuon_chi2LocalMomentum=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalMomentum","MinMuon_chi2LocalMomentum",50,0,5,"Min Inner/Outer track momentum #chi^{2}","");
   MaxMuon_chi2LocalMomentum=HConfig.GetTH1D(Name+"_MaxMuon_chi2LocalMomentum","MaxMuon_chi2LocalMomentum",50,0,50,"Max Inner/Outer track momentum #chi^{2}","");

   Muon1_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon1_chi2LocalMomentum","Muon1_chi2LocalMomentum",50,0,50,"#mu_{1} Inner/Outer track momentum #chi^{2}","");
   Muon2_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon2_chi2LocalMomentum","Muon2_chi2LocalMomentum",50,0,50,"#mu_{2} Inner/Outer track momentum #chi^{2}","");
   Muon3_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon3_chi2LocalMomentum","Muon3_chi2LocalMomentum",50,0,50,"#mu_{3} Inner/Outer track momentum #chi^{2}","");


   MinMatchedStations=HConfig.GetTH1D(Name+"_MinMatchedStations","MinMatchedStations",10,-0.5,9.5,"min matched stations","");
   MaxMatchedStations=HConfig.GetTH1D(Name+"_MaxMatchedStations","MaxMatchedStations",10,-0.5,9.5,"max matched stations ","");
   Mu1MatchedStations=HConfig.GetTH1D(Name+"_Mu1MatchedStations","Mu1MatchedStations",10,-0.5,9.5,"#mu_{1} matched stations","");
   Mu2MatchedStations=HConfig.GetTH1D(Name+"_Mu2MatchedStations","Mu2MatchedStations",10,-0.5,9.5,"#mu_{2} matched stations","");
   Mu3MatchedStations=HConfig.GetTH1D(Name+"_Mu3MatchedStations","Mu3MatchedStations",10,-0.5,9.5,"#mu_{3} matched stations","");


   MinMuon_numberOfChambers=HConfig.GetTH1D(Name+"_MinMuon_numberOfChambers","MinMuon_numberOfChambers",10,-0.5,9.5,"min numberOfChambers","");
   MaxMuon_numberOfChambers=HConfig.GetTH1D(Name+"_MaxMuon_numberOfChambers","MaxMuon_numberOfChambers",10,-0.5,9.5,"max numberOfChambers","");
   Mu1Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu1Muon_numberOfChambers","Mu1Muon_numberOfChambers",10,-0.5,9.5,"#mu_{1} numberOfChambers","");
   Mu2Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu2Muon_numberOfChambers","Mu2Muon_numberOfChambers",10,-0.5,9.5,"#mu_{2} numberOfChambers","");
   Mu3Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu3Muon_numberOfChambers","Mu3Muon_numberOfChambers",10,-0.5,9.5,"#mu_{3} numberOfChambers","");



   MinMuon_numberOfMatches=HConfig.GetTH1D(Name+"_MinMuon_numberOfMatches","MinMuon_numberOfMatches",10,-0.5,9.5,"min numberOfMatches","");
   MaxMuon_numberOfMatches=HConfig.GetTH1D(Name+"_MaxMuon_numberOfMatches","MaxMuon_numberOfMatches",10,-0.5,9.5,"max numberOfMatches","");
   Mu1Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu1Muon_numberOfMatches","Mu1Muon_numberOfMatches",10,-0.5,9.5,"#mu_{1} numberOfMatches","");
   Mu2Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu2Muon_numberOfMatches","Mu2Muon_numberOfMatches",10,-0.5,9.5,"#mu_{2} numberOfMatches","");
   Mu3Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu3Muon_numberOfMatches","Mu3Muon_numberOfMatches",10,-0.5,9.5,"#mu_{3} numberOfMatches","");

   MinMuon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_MinMuon_hitPattern_numberOfValidMuonHits","MinMuon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"min numberOfValidMuonHits","");
   MaxMuon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_MaxMuon_hitPattern_numberOfValidMuonHits","MaxMuon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"max numberOfValidMuonHits","");
   Mu1Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu1Muon_hitPattern_numberOfValidMuonHits","Mu1Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{1} numberOfValidMuonHits","");
   Mu2Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu2Muon_hitPattern_numberOfValidMuonHits","Mu2Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{2} numberOfValidMuonHits","");
   Mu3Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu3Muon_hitPattern_numberOfValidMuonHits","Mu3Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{3} numberOfValidMuonHits","");



   MaxVertexPairQuality=HConfig.GetTH1D(Name+"_MaxVertexPairQuality","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events");
   MinVertexPairQuality=HConfig.GetTH1D(Name+"_MinVertexPairQuality","MinVertexPairQuality",20,0,2,"minvertex pair quality","Events");


   deltaMuZ12 = HConfig.GetTH1D(Name+"_deltaMuZ12","deltaMuZ12",30,0,0.6,"#Delta z (#mu_{1}-#mu_{2}), cm","");
   deltaMuZ13 = HConfig.GetTH1D(Name+"_deltaMuZ13","deltaMuZ13",30,0,0.6,"#Delta z (#mu_{1}-#mu_{3}), cm","");
   deltaMuZ23 = HConfig.GetTH1D(Name+"_deltaMuZ23","deltaMuZ23",30,0,0.6,"#Delta z (#mu_{2}-#mu_{3}), cm","");
   MaxdeltaMuZ = HConfig.GetTH1D(Name+"_MaxdeltaMuZ","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","");
   MindeltaMuZ = HConfig.GetTH1D(Name+"_MindeltaMuZ","MindeltaMuZ",30,0,0.2,"Min #Delta z (#mu-#mu), cm","");
   MaxMuonsDca=HConfig.GetTH1D(Name+"_MaxMuonsDca","MaxMuonsDca",50,0,0.10,"Max distance between muons","");
   MinMuonsDca=HConfig.GetTH1D(Name+"_MinMuonsDca","MinMuonsDca",50,0,0.04,"Min distance between muons","");



   NSV=HConfig.GetTH1D(Name+"_NSV","NSV",8,-0.5,7.5,"N vertices in the tau cone","");

   SVDeltaR=HConfig.GetTH1D(Name+"_SVDeltaR","SVDeltaR",50,0,0.3,"#Delta R (#vec{#tau}  - #vec{Vertex-PV})","");
   SVDistance=HConfig.GetTH1D(Name+"_SVDistance","SVDistance",100,0,0.5,"Distance(SV  - Vertex),cm",""); 
   SV_Mass=HConfig.GetTH1D(Name+"_SV_Mass","SV_Mass",50,0.2,2.5,"VertexMass, GeV","");
   NtracksClose=HConfig.GetTH1D(Name+"_NtracksClose","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","");


   Mu1TrackMass=HConfig.GetTH1D(Name+"_Mu1TrackMass","Mu1TrackMass",60,0.2,4.5,"M_{#mu1-track}, GeV","");
   Mu2TrackMass=HConfig.GetTH1D(Name+"_Mu2TrackMass","Mu2TrackMass",60,0.2,4.5,"M_{#mu2-track}, GeV","");
   Mu3TrackMass=HConfig.GetTH1D(Name+"_Mu3TrackMass","Mu3TrackMass",60,0.2,4.5,"M_{#mu3-track}, GeV","");



   Iso02=HConfig.GetTH1D(Name+"_Iso02","Iso02",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
   Iso04=HConfig.GetTH1D(Name+"_Iso04","Iso04",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
   Iso06=HConfig.GetTH1D(Name+"_Iso06","Iso06",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
   Iso08=HConfig.GetTH1D(Name+"_Iso08","Iso08",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
   Iso1=HConfig.GetTH1D(Name+"_Iso1","Iso1",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.0");
   Iso12=HConfig.GetTH1D(Name+"_Iso12","Iso12",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.2");
   Iso14=HConfig.GetTH1D(Name+"_Iso14","Iso14",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.4");
   Iso16=HConfig.GetTH1D(Name+"_Iso16","Iso16",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
   Iso18=HConfig.GetTH1D(Name+"_Iso18","Iso18",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
   Iso2=HConfig.GetTH1D(Name+"_Iso2","Iso2",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 2");




   Iso02Mu1=HConfig.GetTH1D(Name+"_Iso02Mu1","Iso02Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
   Iso04Mu1=HConfig.GetTH1D(Name+"_Iso04Mu1","Iso04Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
   Iso06Mu1=HConfig.GetTH1D(Name+"_Iso06Mu1","Iso06Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
   Iso08Mu1=HConfig.GetTH1D(Name+"_Iso08Mu1","Iso08Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
   Iso1Mu1=HConfig.GetTH1D(Name+"_Iso1Mu1","Iso1Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0");
   Iso12Mu1=HConfig.GetTH1D(Name+"_Iso12Mu1","Iso12Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.2");
   Iso14Mu1=HConfig.GetTH1D(Name+"_Iso14Mu1","Iso14Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.4");
   Iso16Mu1=HConfig.GetTH1D(Name+"_Iso16Mu1","Iso16Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
   Iso18Mu1=HConfig.GetTH1D(Name+"_Iso18Mu1","Iso18Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
   Iso2Mu1=HConfig.GetTH1D(Name+"_Iso2Mu1","Iso2Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 2");


   Iso02Mu2=HConfig.GetTH1D(Name+"_Iso02Mu2","Iso02Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
   Iso04Mu2=HConfig.GetTH1D(Name+"_Iso04Mu2","Iso04Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
   Iso06Mu2=HConfig.GetTH1D(Name+"_Iso06Mu2","Iso06Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
   Iso08Mu2=HConfig.GetTH1D(Name+"_Iso08Mu2","Iso08Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
   Iso1Mu2=HConfig.GetTH1D(Name+"_Iso1Mu2","Iso1Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.0");
   Iso12Mu2=HConfig.GetTH1D(Name+"_Iso12Mu2","Iso12Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.2");
   Iso14Mu2=HConfig.GetTH1D(Name+"_Iso14Mu2","Iso14Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.4");
   Iso16Mu2=HConfig.GetTH1D(Name+"_Iso16Mu2","Iso16Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
   Iso18Mu2=HConfig.GetTH1D(Name+"_Iso18Mu2","Iso18Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
   Iso2Mu2=HConfig.GetTH1D(Name+"_Iso2Mu2","Iso2Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 2");

   Iso02Mu3=HConfig.GetTH1D(Name+"_Iso02Mu3","Iso02Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
   Iso04Mu3=HConfig.GetTH1D(Name+"_Iso04Mu3","Iso04Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
   Iso06Mu3=HConfig.GetTH1D(Name+"_Iso06Mu3","Iso06Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
   Iso08Mu3=HConfig.GetTH1D(Name+"_Iso08Mu3","Iso08Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
   Iso1Mu3=HConfig.GetTH1D(Name+"_Iso1Mu3","Iso1Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.0");
   Iso12Mu3=HConfig.GetTH1D(Name+"_Iso12Mu3","Iso12Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.2");
   Iso14Mu3=HConfig.GetTH1D(Name+"_Iso14Mu3","Iso14Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.4");
   Iso16Mu3=HConfig.GetTH1D(Name+"_Iso16Mu3","Iso16Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
   Iso18Mu3=HConfig.GetTH1D(Name+"_Iso18Mu3","Iso18Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
   Iso2Mu3=HConfig.GetTH1D(Name+"_Iso2Mu3","Iso2Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 2");


   Iso08MuMax  =HConfig.GetTH1D(Name+"_Iso08MuMax","Iso08MuMax",50,0,1.1,"Max I_{#mu}= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
   Iso08MuMin  =HConfig.GetTH1D(Name+"_Iso08MuMin","Iso08MuMin",50,0,1.1,"Min I_{#mu}= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T}(tracks))","#Delta R < 0.8 ");

   MindcaTrackSV=HConfig.GetTH1D(Name+"_MindcaTrackSV","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","");
   dcaTrackPV=HConfig.GetTH1D(Name+"_dcaTrackPV","dcaTrackPV",50,0,0.1,"distance of closest approach to PV","");

   NBJet4pi=HConfig.GetTH1D(Name+"_NBJet4pi","NBJet4pi",10,-0.5,9.5,"N Jets","");


   BTagCSV =HConfig.GetTH1D(Name+"_BTagCSV","BTagCSV",50,0.,1.,"BTagCSV ","");
   BTagMVA =HConfig.GetTH1D(Name+"_BTagMVA","BTagMVA",50,-1.,1.,"BTagMVA ","");
   BTagCVSB=HConfig.GetTH1D(Name+"_BTagCVSB","BTagCVSB",50,-1.,1.,"BTagCVSB ","");

   NBJet4piSH=HConfig.GetTH1D(Name+"_NBJet4piSH","NBJet4piSH",10,-0.5,9.5,"N Jets SH","");
   NBJet4piOH=HConfig.GetTH1D(Name+"_NBJet4piOH","NBJet4piOH",10,-0.5,9.5,"N Jets OH","");

   BTagCSVSH =HConfig.GetTH1D(Name+"_BTagCSVSH","BTagCSVSH",50,0.,1.,"BTagCSV SH","");
   BTagMVASH =HConfig.GetTH1D(Name+"_BTagMVASH","BTagMVASH",50,-1.,1.,"BTagMVA SH","");
   BTagCVSBSH=HConfig.GetTH1D(Name+"_BTagCVSBSH","BTagCVSBSH",50,-1.,1.,"BTagCVSB SH","");

   BTagCSVOH =HConfig.GetTH1D(Name+"_BTagCSVOH","BTagCSVOH",50,0.,1.,"BTagCSV OH ","");
   BTagMVAOH =HConfig.GetTH1D(Name+"_BTagMVAOH","BTagMVAOH",50,-1.,1.,"BTagMVA OH","");
   BTagCVSBOH=HConfig.GetTH1D(Name+"_BTagCVSBOH","BTagCVSBOH",50,-1.,1.,"BTagCVSB OH","");




   BTagCSVSHMatchedToTau =HConfig.GetTH1D(Name+"_BTagCSVSHMatchedToTau","BTagCSVSHMatchedToTau",50,0.,1.,"BTagCSV SH","");
   BTagMVASHMatchedToTau =HConfig.GetTH1D(Name+"_BTagMVASHMatchedToTau","BTagMVASHMatchedToTau",50,-1.,1.,"BTagMVA SH","");
   BTagCVSBSHMatchedToTau=HConfig.GetTH1D(Name+"_BTagCVSBSHMatchedToTau","BTagCVSBSHMatchedToTau",50,-1.,1.,"BTagCVSB SH","");

   BTagCSVOHMatchedToTau =HConfig.GetTH1D(Name+"_BTagCSVOHMatchedToTau","BTagCSVOHMatchedToTau",50,0.,1.,"BTagCSV OH ","");
   BTagMVAOHMatchedToTau =HConfig.GetTH1D(Name+"_BTagMVAOHMatchedToTau","BTagMVAOHMatchedToTau",50,-1.,1.,"BTagMVA OH","");
   BTagCVSBOHMatchedToTau=HConfig.GetTH1D(Name+"_BTagCVSBOHMatchedToTau","BTagCVSBOHMatchedToTau",50,-1.,1.,"BTagCVSB OH","");



   BTagCSVSHVsNJets =HConfig.GetTH2D(Name+"_BTagCSVSHVsNJets","BTagCSVSHVsNJets",50,0.,1.,10,-0.5,9.5,"BTagCSV SH","N");
   BTagMVASHVsNJets =HConfig.GetTH2D(Name+"_BTagMVASHVsNJets","BTagMVASHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagMVA SH","N");
   BTagCVSBSHVsNJets=HConfig.GetTH2D(Name+"_BTagCVSBSHVsNJets","BTagCVSBSHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagCVSB SH","N");

   BTagCSVOHVsNJets =HConfig.GetTH2D(Name+"_BTagCSVOHVsNJets","BTagCSVOHVsNJets",50,0.,1.,10,-0.5,9.5,"BTagCSV OH ","N");
   BTagMVAOHVsNJets =HConfig.GetTH2D(Name+"_BTagMVAOHVsNJets","BTagMVAOHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagMVA OH","N");
   BTagCVSBOHVsNJets=HConfig.GetTH2D(Name+"_BTagCVSBOHVsNJets","BTagCVSBOHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagCVSB OH","N");

   Muon1ImpactAngle =HConfig.GetTH1D(Name+"_Muon1ImpactAngle","Muon1ImpactAngle",30,-1,1,"#mu_{1} impact angle","");
   Muon2ImpactAngle =HConfig.GetTH1D(Name+"_Muon2ImpactAngle","Muon2ImpactAngle",30,-1,1,"#mu_{2} impact angle","");
   Muon3ImpactAngle =HConfig.GetTH1D(Name+"_Muon3ImpactAngle","Muon3ImpactAngle",30,-1,1,"#mu_{3} impact angle","");
   MinMuonImpacAngle =HConfig.GetTH1D(Name+"_MinMuonImpactAngle","MinMuonImpactAngle",30,-1,1,"Min #mu impact angle","");
   MaxMuonImpacAngle =HConfig.GetTH1D(Name+"_MaxMuonImpactAngle","MaxMuonImpactAngle",30,-1,1,"Max #mu impact angle","");


   Muon1StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon1StandardSelectorPass","Muon1StandardSelectorPass",23,-0.5,22.5,"#mu_{1} standard selector(pass); bin 0 - no ID","Events");
   Muon2StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon2StandardSelectorPass","Muon2StandardSelectorPass",23,-0.5,22.5,"#mu_{2} standard selector(pass); bin 0 - no ID","Events");
   Muon3StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon3StandardSelectorPass","Muon3StandardSelectorPass",23,-0.5,22.5,"#mu_{3} standard selector(pass); bin 0 - no ID","Events");

   Muon1StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon1StandardSelectorFail","Muon1StandardSelectorFail",23,-0.5,22.5,"#mu_{1} standard selector(fail); bin 0 - no ID","Events");
   Muon2StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon2StandardSelectorFail","Muon2StandardSelectorFail",23,-0.5,22.5,"#mu_{2} standard selector(fail); bin 0 - no ID","Events");
   Muon3StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon3StandardSelectorFail","Muon3StandardSelectorFail",23,-0.5,22.5,"#mu_{3} standard selector(fail); bin 0 - no ID","Events");




   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove

}

void  FillTMVATrees::Store_ExtraDist(){ 


   Extradist1d.push_back(&Muon1Pt);
   Extradist1d.push_back(&Muon2Pt);
   Extradist1d.push_back(&Muon3Pt);

   Extradist1d.push_back(&TauEta);
   Extradist1d.push_back(&TauPt);
   Extradist1d.push_back(&TauP);


   Extradist1d.push_back(&SVPVTauDirAngle);
   Extradist1d.push_back(&FLSignificance);
   Extradist1d.push_back(&FL);
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

   Extradist1d.push_back(&Isolation_maxdxy);

   Extradist1d.push_back(&EventMassResolution_PtEtaPhi);
   Extradist2d.push_back(&EMR_tau_eta);

   Extradist1d.push_back(&MaxMuonsDca);
   Extradist1d.push_back(&MinMuonsDca);


   Extradist1d.push_back(&deltaMuZ12);
   Extradist1d.push_back(&deltaMuZ13);
   Extradist1d.push_back(&deltaMuZ23);
   Extradist1d.push_back(&MaxdeltaMuZ);
   Extradist1d.push_back(&MindeltaMuZ);

   Extradist1d.push_back(&NSV);
   Extradist1d.push_back(&SVDeltaR);
   Extradist1d.push_back(&SVDistance);
   Extradist1d.push_back(&SV_Mass);
   Extradist1d.push_back(&NtracksClose);

   Extradist1d.push_back(&VertexMu1D0SigPVReco);
   Extradist1d.push_back(&VertexMu2D0SigPVReco);
   Extradist1d.push_back(&VertexMu3D0SigPVReco);
   Extradist1d.push_back(&MaxD0SigPV);
   Extradist1d.push_back(&MinD0SigPV);


   Extradist1d.push_back(&VertexMu1D0SigBSReco);
   Extradist1d.push_back(&VertexMu2D0SigBSReco);
   Extradist1d.push_back(&VertexMu3D0SigBSReco);
   Extradist1d.push_back(&MaxD0SigBS);
   Extradist1d.push_back(&MinD0SigBS);


   Extradist1d.push_back(&VertexMu1D0SigSVReco);
   Extradist1d.push_back(&VertexMu2D0SigSVReco);
   Extradist1d.push_back(&VertexMu3D0SigSVReco);
   Extradist1d.push_back(&MaxD0SigSV);
   Extradist1d.push_back(&MinD0SigSV);

   Extradist1d.push_back(&MinMuon_chi2LocalPosition);
   Extradist1d.push_back(&MaxMuon_chi2LocalPosition);

   Extradist1d.push_back(&Muon1_chi2LocalPosition);
   Extradist1d.push_back(&Muon2_chi2LocalPosition);
   Extradist1d.push_back(&Muon3_chi2LocalPosition);


   Extradist1d.push_back(&MinMuon_chi2LocalMomentum);
   Extradist1d.push_back(&MaxMuon_chi2LocalMomentum);

   Extradist1d.push_back(&Muon1_chi2LocalMomentum);
   Extradist1d.push_back(&Muon2_chi2LocalMomentum);
   Extradist1d.push_back(&Muon3_chi2LocalMomentum);

   Extradist1d.push_back(&MintrkKink);
   Extradist1d.push_back(&MaxtrkKink);
   Extradist1d.push_back(&MinglbKink);
   Extradist1d.push_back(&MaxglbKink);

   Extradist1d.push_back(&MuonglbkinkSum);

   Extradist1d.push_back(&MaxVertexPairQuality);
   Extradist1d.push_back(&MinVertexPairQuality);


   Extradist1d.push_back(&Mu1TrackMass);
   Extradist1d.push_back(&Mu2TrackMass);
   Extradist1d.push_back(&Mu3TrackMass);

   Extradist1d.push_back(&Iso02);
   Extradist1d.push_back(&Iso04);
   Extradist1d.push_back(&Iso06);
   Extradist1d.push_back(&Iso08);
   Extradist1d.push_back(&Iso1);
   Extradist1d.push_back(&Iso12);
   Extradist1d.push_back(&Iso14);
   Extradist1d.push_back(&Iso16);
   Extradist1d.push_back(&Iso18);
   Extradist1d.push_back(&Iso2);

   Extradist1d.push_back(&Iso02Mu1);
   Extradist1d.push_back(&Iso04Mu1);
   Extradist1d.push_back(&Iso06Mu1);
   Extradist1d.push_back(&Iso08Mu1);
   Extradist1d.push_back(&Iso1Mu1);
   Extradist1d.push_back(&Iso12Mu1);
   Extradist1d.push_back(&Iso14Mu1);
   Extradist1d.push_back(&Iso16Mu1);
   Extradist1d.push_back(&Iso18Mu1);
   Extradist1d.push_back(&Iso2Mu1);

   Extradist1d.push_back(&Iso02Mu2);
   Extradist1d.push_back(&Iso04Mu2);
   Extradist1d.push_back(&Iso06Mu2);
   Extradist1d.push_back(&Iso08Mu2);
   Extradist1d.push_back(&Iso1Mu2);
   Extradist1d.push_back(&Iso12Mu2);
   Extradist1d.push_back(&Iso14Mu2);
   Extradist1d.push_back(&Iso16Mu2);
   Extradist1d.push_back(&Iso18Mu2);
   Extradist1d.push_back(&Iso2Mu2);

   Extradist1d.push_back(&Iso02Mu3);
   Extradist1d.push_back(&Iso04Mu3);
   Extradist1d.push_back(&Iso06Mu3);
   Extradist1d.push_back(&Iso08Mu3);
   Extradist1d.push_back(&Iso1Mu3);
   Extradist1d.push_back(&Iso12Mu3);
   Extradist1d.push_back(&Iso14Mu3);
   Extradist1d.push_back(&Iso16Mu3);
   Extradist1d.push_back(&Iso18Mu3);
   Extradist1d.push_back(&Iso2Mu3);

   Extradist1d.push_back(&Iso08MuMax);
   Extradist1d.push_back(&Iso08MuMin);


   Extradist1d.push_back(&MindcaTrackSV);
   Extradist1d.push_back(&dcaTrackPV);

   Extradist1d.push_back(&MinMatchedStations);
   Extradist1d.push_back(&MaxMatchedStations);
   Extradist1d.push_back(&Mu1MatchedStations);
   Extradist1d.push_back(&Mu2MatchedStations);
   Extradist1d.push_back(&Mu3MatchedStations);


   Extradist1d.push_back(&MinMuon_numberOfChambers);
   Extradist1d.push_back(&MaxMuon_numberOfChambers);
   Extradist1d.push_back(&Mu1Muon_numberOfChambers);
   Extradist1d.push_back(&Mu2Muon_numberOfChambers);
   Extradist1d.push_back(&Mu3Muon_numberOfChambers);



   Extradist1d.push_back(&MinMuon_numberOfMatches);
   Extradist1d.push_back(&MaxMuon_numberOfMatches);
   Extradist1d.push_back(&Mu1Muon_numberOfMatches);
   Extradist1d.push_back(&Mu2Muon_numberOfMatches);
   Extradist1d.push_back(&Mu3Muon_numberOfMatches);

   Extradist1d.push_back(&MinMuon_hitPattern_numberOfValidMuonHits);
   Extradist1d.push_back(&MaxMuon_hitPattern_numberOfValidMuonHits);
   Extradist1d.push_back(&Mu1Muon_hitPattern_numberOfValidMuonHits);
   Extradist1d.push_back(&Mu2Muon_hitPattern_numberOfValidMuonHits);
   Extradist1d.push_back(&Mu3Muon_hitPattern_numberOfValidMuonHits);


   //  Extradist1d.push_back(&NBJet4pi);
   //  Extradist1d.push_back(&BTagCSV);
   //  Extradist1d.push_back(&BTagMVA);
   //  Extradist1d.push_back(&BTagCVSB);

   Extradist1d.push_back(&NBJet4piSH);
   Extradist1d.push_back(&NBJet4piOH);

   Extradist1d.push_back(&BTagCSVSH);
   Extradist1d.push_back(&BTagMVASH);
   Extradist1d.push_back(&BTagCVSBSH);

   Extradist1d.push_back(&BTagCSVOH);
   Extradist1d.push_back(&BTagMVAOH);
   Extradist1d.push_back(&BTagCVSBOH);

   Extradist1d.push_back(&BTagCSVSHMatchedToTau);
   Extradist1d.push_back(&BTagMVASHMatchedToTau);
   Extradist1d.push_back(&BTagCVSBSHMatchedToTau);

   Extradist1d.push_back(&BTagCSVOHMatchedToTau);
   Extradist1d.push_back(&BTagMVAOHMatchedToTau);
   Extradist1d.push_back(&BTagCVSBOHMatchedToTau);





   Extradist2d.push_back(&BTagCSVSHVsNJets);
   Extradist2d.push_back(&BTagMVASHVsNJets);
   Extradist2d.push_back(&BTagCVSBSHVsNJets);

   Extradist2d.push_back(&BTagCSVOHVsNJets);
   Extradist2d.push_back(&BTagMVAOHVsNJets);
   Extradist2d.push_back(&BTagCVSBOHVsNJets);





   Extradist1d.push_back(&Muon1ImpactAngle);
   Extradist1d.push_back(&Muon2ImpactAngle);
   Extradist1d.push_back(&Muon3ImpactAngle);
   Extradist1d.push_back(&MinMuonImpacAngle);
   Extradist1d.push_back(&MaxMuonImpacAngle);



   Extradist1d.push_back(&Muon1StandardSelectorPass);
   Extradist1d.push_back(&Muon2StandardSelectorPass);
   Extradist1d.push_back(&Muon3StandardSelectorPass);

   Extradist1d.push_back(&Muon1StandardSelectorFail);
   Extradist1d.push_back(&Muon2StandardSelectorFail);
   Extradist1d.push_back(&Muon3StandardSelectorFail);
}



void  FillTMVATrees::doEvent(){ 


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
         if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                     DoubleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9")) TripleMuFired = Ntp->L1Decision(il1);
         if(L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2"))                       DoubleMuFired = Ntp->L1Decision(il1);
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
            pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS1)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.01); // extend the veto region for two global + tracker category
            pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.01); // extend the veto region for two global + tracker category
         }
         pass.at(OmegaVetoOS1) = (fabs(value.at(OmegaVetoOS1)-PDG_Var::Omega_mass())> omegaVetoSigma);
         pass.at(OmegaVetoOS2) = (fabs(value.at(OmegaVetoOS2)-PDG_Var::Omega_mass())> omegaVetoSigma);
         pass.at(DsGenMatch) = ( value.at(DsGenMatch) < cut.at(DsGenMatch) );
         pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
         if(id!=1) pass.at(TauMassCut) = ( value.at(TauMassCut) > (tauMinMass_) && value.at(TauMassCut) < (tauMaxMass_));
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
         pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS1)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.01); // extend the veto region for two global + tracker category
         pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma + 0.01); // extend the veto region for two global + tracker category
      }
      pass.at(OmegaVetoOS1) = (fabs(value.at(OmegaVetoOS1)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(OmegaVetoOS2) = (fabs(value.at(OmegaVetoOS2)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(DsGenMatch) = ( value.at(DsGenMatch) < cut.at(DsGenMatch) );
      pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );

      if(id!=1) pass.at(TauMassCut) = ( value.at(TauMassCut) > (tauMinMass_) && value.at(TauMassCut) < (tauMaxMass_));
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
   }

   // Weights
   double wobs=1;
   double w;
   double w_normalization=0.76;

   if(!Ntp->isData()){
      w = w_normalization*(puWeights->GetBinContent(Ntp->TruthNumberOfInteraction())); // Weight MC according to truth number of vertices
   }
   //  No weights to data
   else{w=1;}

   bool status=AnalysisCuts(t,w,wobs);

   if(status){

      /*
         int NJets(){return Ntp->Jet_p4->size();}
         TLorentzVector Jet_P4(unsigned int i){return TLorentzVector(Ntp->Jet_p4->at(i)->at(1), Ntp->Jet_p4->at(i)->at(2), Ntp->Jet_p4->at(i)->at(3), Ntp->Jet_p4->at(i)->at(0));}
         double JetBTagCVSB(unsigned int i){return Ntp->Jet_BTagCVSB->at(i);}
         double JetBTagMVA(unsigned int i){return Ntp->Jet_BTagMVA->at(i);}
         double JetBTagCSV(unsigned int i){return Ntp->Jet_BTagCSV->at(i);}
         */

      Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      std::vector<unsigned int> indices;
      indices.push_back(Muon_index_1);
      indices.push_back(Muon_index_2);
      indices.push_back(Muon_index_3);

      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
      TLorentzVector TauLV = Muon1LV + Muon2LV + Muon3LV;
      TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);

      int OH_jet_count(0);
      int SH_jet_count(0);
      int jet_index_closest_to_tau(-1);
      int jet_index_closest_to_tau_opposite_in_phi(-1);
      double TJDeltaR(99.);
      double OHDeltaPhi(TMath::Pi());
      for(int ijet= 0; ijet <  Ntp->NJets(); ijet++){
         if(TauLV.DeltaR(Ntp->Jet_P4(ijet)) < TJDeltaR){

            TJDeltaR = TauLV.DeltaR(Ntp->Jet_P4(ijet));
            jet_index_closest_to_tau = ijet;

         }



         if(fabs(Ntp->Jet_P4(ijet).Phi()  - TauLV.Phi()) < TMath::Pi() ){
            SH_jet_count++;
            BTagCSVSH.at(t).Fill(Ntp->JetBTagCSV(ijet),w);
            BTagMVASH.at(t).Fill(Ntp->JetBTagMVA(ijet),w);
            BTagCVSBSH.at(t).Fill(Ntp->JetBTagCVSB(ijet),w);
         }


         if(fabs(Ntp->Jet_P4(ijet).Phi()  - TauLV.Phi()) > TMath::Pi() ){
            OH_jet_count++;  
            BTagCSVOH.at(t).Fill(Ntp->JetBTagCSV(ijet),w);
            BTagMVAOH.at(t).Fill(Ntp->JetBTagMVA(ijet),w);
            BTagCVSBOH.at(t).Fill(Ntp->JetBTagCVSB(ijet),w);
         }


      }
      for(int ijet= 0; ijet <  Ntp->NJets(); ijet++){
         if(jet_index_closest_to_tau!=-1){
            //    if(fabs(Ntp->Jet_P4(ijet).Phi() - fabs(Ntp->Jet_P4(jet_index_closest_to_tau).Phi() - TMath::Pi())) < OHDeltaPhi  )  {
            if(fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() )) < OHDeltaPhi  )  {
               OHDeltaPhi = fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() ));
               jet_index_closest_to_tau_opposite_in_phi = ijet;
            }

         }
         }
         NBJet4piSH.at(t).Fill(SH_jet_count,w);
         NBJet4piOH.at(t).Fill(OH_jet_count,w);

         if(jet_index_closest_to_tau!=-1){
            BTagCSVSHMatchedToTau.at(t).Fill(Ntp->JetBTagCSV(jet_index_closest_to_tau),w);
            BTagMVASHMatchedToTau.at(t).Fill(Ntp->JetBTagMVA(jet_index_closest_to_tau),w);
            BTagCVSBSHMatchedToTau.at(t).Fill(Ntp->JetBTagCVSB(jet_index_closest_to_tau),w);
         }
         if(jet_index_closest_to_tau_opposite_in_phi!=-1){
            BTagCSVOHMatchedToTau.at(t).Fill(Ntp->JetBTagCSV(jet_index_closest_to_tau_opposite_in_phi),w);
            BTagMVAOHMatchedToTau.at(t).Fill(Ntp->JetBTagMVA(jet_index_closest_to_tau_opposite_in_phi),w);
            BTagCVSBOHMatchedToTau.at(t).Fill(Ntp->JetBTagCVSB(jet_index_closest_to_tau_opposite_in_phi),w);
         }




         //------------- Muon ID

         for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_1).size(); iMuSelector++ ){
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)==1)  Muon1StandardSelectorPass.at(t).Fill(iMuSelector,1);
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)!=1)  Muon1StandardSelectorFail.at(t).Fill(iMuSelector,1);
         }

         for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_2).size(); iMuSelector++ ){
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)==1)  Muon2StandardSelectorPass.at(t).Fill(iMuSelector,1);
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)!=1)  Muon2StandardSelectorFail.at(t).Fill(iMuSelector,1);
         }

         for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_3).size(); iMuSelector++ ){
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)==1)  Muon3StandardSelectorPass.at(t).Fill(iMuSelector,1);
            if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)!=1)  Muon3StandardSelectorFail.at(t).Fill(iMuSelector,1);
         }



         //------------- Muon ID




         unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
         unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
         unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

         std::vector<unsigned int> EtaSortedIndices;


         EtaSortedIndices.push_back(Muon_Eta_index_1);
         EtaSortedIndices.push_back(Muon_Eta_index_2);
         EtaSortedIndices.push_back(Muon_Eta_index_3);
         EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());
         EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);



         float tauMassRes = Ntp->TauMassResolution(EtaSortedIndices,1,false);



         Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),w);
         Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),w);
         Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),w);


         TauEta.at(t).Fill(TauLV.Eta(),w);
         TauPt.at(t).Fill(TauLV.Pt(),w);
         TauP.at(t).Fill(TauLV.P(),w);




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

         TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));

         TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_1),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
         TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_2),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
         TVector3 Mu3ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_3),Ntp->Vertex_MatchedPrimaryVertex(final_idx));



         //---------------------------------------------------------------
         //------------------ calculate var_flightLenSig ---------------------
         TMatrixTSym<double> fls_PVcov = Ntp->Vertex_PrimaryVertex_Covariance(final_idx);
         TMatrixTSym<double> fls_SVcov = Ntp->Vertex_Signal_KF_Covariance(final_idx);



         //      if (Ntp->Vertex_RefitPVisValid(final_idx)==1)
         {

            Muon1ImpactAngle.at(t).Fill(SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),w );
            Muon2ImpactAngle.at(t).Fill(SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),w );
            Muon3ImpactAngle.at(t).Fill(SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag()),w );
            MinMuonImpacAngle.at(t).Fill(std::min({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())}),w);
            MaxMuonImpacAngle.at(t).Fill(std::min({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())}),w);

            // ----- Fill the histograms -----
            SVPVTauDirAngle.at(t).Fill(var_svpvTauAngle,w);
            FLSignificance.at(t).Fill(var_flightLenSig,w);
            FL.at(t).Fill(SVPV.Mag(),w);
            VertexChi2KF.at(t).Fill(var_vertexKFChi2,w);

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




            deltaMuZ12.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),1);
            deltaMuZ13.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),1);
            deltaMuZ23.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),1);


            MaxdeltaMuZ.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),
                     fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),
                     fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z())}),1);

            MindeltaMuZ.at(t).Fill( std::min({fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),
                     fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),
                     fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z())}),1);

            MaxMuonsDca.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)}),1);
            MinMuonsDca.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)}),1);

            var_MaxdeltaMuZ = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
            var_MindeltaMuZ = std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});




            MinMatchedStations.at(t).Fill(std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)}),w);
            MaxMatchedStations.at(t).Fill(std::max({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)}),w);
            Mu1MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_1),w);
            Mu2MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_2),w);
            Mu3MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_3),w);


            MinMuon_numberOfChambers.at(t).Fill(std::min({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)}),w);
            MaxMuon_numberOfChambers.at(t).Fill(std::max({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)}),w);
            Mu1Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_1),w);
            Mu2Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_2),w);
            Mu3Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_3),w);



            MinMuon_numberOfMatches.at(t).Fill(std::min({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)}),w);
            MaxMuon_numberOfMatches.at(t).Fill(std::max({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)}),w);
            Mu1Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_1),w);
            Mu2Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_2),w);
            Mu3Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_3),w);

            MinMuon_hitPattern_numberOfValidMuonHits.at(t).Fill(std::min({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)}),w);
            MaxMuon_hitPattern_numberOfValidMuonHits.at(t).Fill(std::max({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)}),w);
            Mu1Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1),w);
            Mu2Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2),w);
            Mu3Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3),w);



            var_MinMatchedStations = std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});
            var_MaxMatchedStations = std::max({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});
            var_Mu1MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_1);
            var_Mu2MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_2);
            var_Mu3MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_3);


            var_MinMuon_numberOfChambers = std::min({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)});
            var_MaxMuon_numberOfChambers = std::max({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)});
            var_Mu1Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_1);
            var_Mu2Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_2);
            var_Mu3Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_3);



            var_MinMuon_numberOfMatches = std::min({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)});
            var_MaxMuon_numberOfMatches = std::max({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)});
            var_Mu1Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_1);
            var_Mu2Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_2);
            var_Mu3Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_3);

            var_MinMuon_hitPattern_numberOfValidMuonHits = std::min({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)});
            var_MaxMuon_hitPattern_numberOfValidMuonHits = std::max({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)});
            var_Mu1Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1);
            var_Mu2Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2);
            var_Mu3Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3);



            float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0),
                  Ntp->Vertex_d0sig_reco(final_idx,1),
                  Ntp->Vertex_d0sig_reco(final_idx,2)});


            float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0),
                  Ntp->Vertex_d0sig_reco(final_idx,1),
                  Ntp->Vertex_d0sig_reco(final_idx,2)});


            float MaxD0BSSignificance = std::max({Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),
                  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),
                  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2)});


            float MinD0BSSignificance = std::min({Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),
                  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),
                  Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2)});



            float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(final_idx,0),
                  Ntp->Vertex_d0sigSV_reco(final_idx,1),
                  Ntp->Vertex_d0sigSV_reco(final_idx,2)});


            float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(final_idx,0),
                  Ntp->Vertex_d0sigSV_reco(final_idx,1),
                  Ntp->Vertex_d0sigSV_reco(final_idx,2)});

            VertexMu1D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,0),w);
            VertexMu2D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,1),w);
            VertexMu3D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,2),w);

            MaxD0SigPV.at(t).Fill(MaxD0Significance,1);
            MinD0SigPV.at(t).Fill(MinD0Significance,1);



            VertexMu1D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),w);
            VertexMu2D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),w);
            VertexMu3D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2),w);

            MaxD0SigBS.at(t).Fill(MaxD0BSSignificance,1);
            MinD0SigBS.at(t).Fill(MinD0BSSignificance,1);



            VertexMu1D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,0),w);
            VertexMu2D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,1),w);
            VertexMu3D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,2),w);


            MaxD0SigSV.at(t).Fill(MaxD0SVSignificance,1);
            MinD0SigSV.at(t).Fill(MinD0SVSignificance,1);


            MinMuon_chi2LocalPosition.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  }));
            MaxMuon_chi2LocalPosition.at(t).Fill(std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  }));


            Muon1_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),1);
            Muon2_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),1);
            Muon3_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),1);


            MinMuon_chi2LocalMomentum.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  }));
            MaxMuon_chi2LocalMomentum.at(t).Fill(std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  }));


            Muon1_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1), 1);
            Muon2_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2), 1);
            Muon3_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3), 1);



            MintrkKink.at(t).Fill(std::min({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)}),w);
            MaxtrkKink.at(t).Fill(std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)}),w);
            MinglbKink.at(t).Fill(std::min({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)}),w);
            MaxglbKink.at(t).Fill(std::max({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)}),w);

            MuonglbkinkSum.at(t).Fill(var_sumMuTrkKinkChi2,w);

            MaxVertexPairQuality.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)}),w);
            MinVertexPairQuality.at(t).Fill(  std::min({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)}),w);


            //  ----------------------------- secondary vertices ----------------
            int NumberOfPrimaryVertices(0);
            for(unsigned int iVertex=0; iVertex < Ntp->NSecondaryVertices(); iVertex++){
               SV_Mass.at(t).Fill(Ntp->SecondaryVertexMass(iVertex),1);
               TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
               TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(iVertex),Ntp->Vertex_MatchedPrimaryVertex(final_idx));

               SVDeltaR.at(t).Fill(SVfakePV.DeltaR(SVsignalPV),1);
               SVDistance.at(t).Fill((Ntp->Vertex_Signal_KF_pos(final_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag(),1);


               if(SVfakePV.DeltaR(SVsignalPV) < 0.3 && (Ntp->Vertex_Signal_KF_pos(final_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag() > 0.015){ // sv in the tau cone but  displaced
                  NumberOfPrimaryVertices++;

               }
            }
            NSV.at(t).Fill(NumberOfPrimaryVertices,1);


            // -----------------------------------------------------------------




            // --------------------------- Isolation --------------------------------

            int NcloseTracksCount(0);
            double SumPT02(0),SumPT04(0),SumPT06(0),SumPT08(0),SumPT1(0),SumPT12(0),SumPT14(0),SumPT16(0),SumPT18(0),SumPT2(0);
            double Iso08Muon1,Iso08Muon2,Iso08Muon3;
            int TrackIndex(0);
            int TrackIndex_closestToPV(0);
            double TrackPTtreschold(0.8);
            double dca_temp(999.);
            double dcaPV_temp(999.);
            for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){


               if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   
                        pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.03)
               {
                  NcloseTracksCount++;
               }
               if( sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2) ) <  dca_temp){
                  dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2));
                  TrackIndex = i;
               }

               if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05 && 
                     sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.05){


                  if( sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2) ) <  dcaPV_temp){
                     dcaPV_temp = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2));
                     TrackIndex_closestToPV = i;
                  }


                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.2){
                     SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.4){
                     SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.6){
                     SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.8){
                     SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.0){
                     SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.2){
                     SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.4){
                     SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }  
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.6){
                     SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.8){
                     SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 2.0){
                     SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
               }
            }


            NtracksClose.at(t).Fill(NcloseTracksCount,1);


            Iso02.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT02),1);
            Iso04.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT04),1);
            Iso06.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT06),1);
            Iso08.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT08),1);
            Iso1.at(t).Fill(TauRefitLV.Pt()/  (TauRefitLV.Pt() + SumPT1),1);
            Iso12.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT12),1);
            Iso14.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT14),1);
            Iso16.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT16),1);
            Iso18.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT18),1);
            Iso2.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT2),1);


            var_Iso02 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT02);
            var_Iso04 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT04);
            var_Iso06 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT06);
            var_Iso08 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT08);
            var_Iso1 = TauRefitLV.Pt()/  (TauRefitLV.Pt() + SumPT1);
            var_Iso12 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT12);



            if(Ntp->NIsolationTrack(final_idx)!=0){ 
               MindcaTrackSV.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2)),1);
            }

            if(Ntp->NIsolationTrack(final_idx)!=0){ 
               dcaTrackPV.at(t).Fill(sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV),2)),1 );
            }


            // ---------------------- I_mu1
            SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
            for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

               if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
                     && Ntp->IsolationTrack_DocaMu1(final_idx,i) < 0.1){
                  if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(0))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){
                     Mu1TrackMass.at(t).Fill(  (Muon1LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
                     var_Mu1TrackMass = (Muon1LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
                  }else var_Mu1TrackMass = -1;

                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.2){
                     SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.4){
                     SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.6){
                     SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.8){
                     SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.0){
                     SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.2){
                     SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.4){
                     SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.6){
                     SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.8){
                     SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 2.0){
                     SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
               }
            }

            Iso02Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT02),1);
            Iso04Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT04),1);
            Iso06Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT06),1);
            Iso08Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08),1); Iso08Muon1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08);
            Iso1Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1);
            Iso12Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT12),1);
            Iso14Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT14),1);
            Iso16Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT16),1);
            Iso18Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT18),1);
            Iso2Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT2),1);


            var_Iso02Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT02);
            var_Iso04Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT04);
            var_Iso06Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT06);
            var_Iso08Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08);
            var_Iso1Mu1 = Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1);
            var_Iso12Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT12);



            // ---------------------- I_mu2

            SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
            for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

               if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
                     && Ntp->IsolationTrack_DocaMu2(final_idx,i) < 0.1){


                  if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(1))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){

                     Mu2TrackMass.at(t).Fill(  (Muon2LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
                     var_Mu2TrackMass = (Muon2LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
                  }else var_Mu2TrackMass = -1;


                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.2){
                     SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.4){
                     SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.6){
                     SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.8){
                     SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.0){
                     SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.2){
                     SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.4){
                     SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.6){
                     SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.8){
                     SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 2.0){
                     SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
               }
            }

            Iso02Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT02),1);
            Iso04Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT04),1);
            Iso06Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT06),1);
            Iso08Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08),1);Iso08Muon2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08);
            Iso1Mu2.at(t).Fill(Muon2LV.Pt()/  (Muon2LV.Pt() + SumPT1),1);
            Iso12Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT12),1);
            Iso14Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT14),1);
            Iso16Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT16),1);
            Iso18Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT18),1);
            Iso2Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT2),1);

            var_Iso02Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT02);
            var_Iso04Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT04);
            var_Iso06Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT06);
            var_Iso08Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08);
            var_Iso1Mu2 = Muon2LV.Pt()/  (Muon2LV.Pt() + SumPT1);
            var_Iso12Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT12);



            // ---------------------- I_mu3

            SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
            for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

               if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
                     && Ntp->IsolationTrack_DocaMu3(final_idx,i) < 0.1){


                  if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(2))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){
                     Mu3TrackMass.at(t).Fill(  (Muon3LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
                     var_Mu3TrackMass = (Muon3LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
                  }else var_Mu3TrackMass = -1;


                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.2){
                     SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.4){
                     SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.6){
                     SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.8){
                     SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.0){
                     SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.2){
                     SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.4){
                     SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.6){
                     SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.8){
                     SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
                  if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 2.0){
                     SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
                  }
               }
            }

            Iso02Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT02),1);    
            Iso04Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT04),1);
            Iso06Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT06),1);
            Iso08Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08),1);Iso08Muon3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08);
            Iso1Mu3.at(t).Fill(Muon3LV.Pt()/  (Muon3LV.Pt() + SumPT1),1);
            Iso12Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT12),1);
            Iso14Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT14),1);
            Iso16Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT16),1);
            Iso18Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT18),1);
            Iso2Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT2),1);


            var_Iso02Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT02);
            var_Iso04Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT04);
            var_Iso06Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT06);
            var_Iso08Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08);
            var_Iso1Mu3 = Muon3LV.Pt()/  (Muon3LV.Pt() + SumPT1);
            var_Iso12Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT12);
            var_Iso08MuMax = std::max({Iso08Muon1,Iso08Muon2,Iso08Muon3});
            var_Iso08MuMin = std::min({Iso08Muon1,Iso08Muon2,Iso08Muon3});


            Iso08MuMax.at(t).Fill(std::max({Iso08Muon1,Iso08Muon2,Iso08Muon3}),1);
            Iso08MuMin.at(t).Fill(std::min({Iso08Muon1,Iso08Muon2,Iso08Muon3}),1);


            // -----------------------------------------------------------------------

            // -------------------------- Fill MVA mini tree ------------------------

            var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx);
            var_svpvTauAngle = SVPV.Angle(TauLV.Vect());

            var_flightLenSig =  Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx),Ntp->Vertex_PrimaryVertex_Covariance(final_idx),
                  Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx));


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

            var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
            var_minMuonsDca = std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});

            var_nsv = NumberOfPrimaryVertices;


            var_VertexMu1D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,0);
            var_VertexMu2D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,1);
            var_VertexMu3D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,2);

            var_MaxD0SigPV = MaxD0Significance;
            var_MinD0SigPV = MinD0Significance;



            var_VertexMu1D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0);
            var_VertexMu2D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1);
            var_VertexMu3D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2);

            var_MaxD0SigBS = MaxD0BSSignificance;
            var_MinD0SigBS = MinD0BSSignificance;



            var_VertexMu1D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,0);
            var_VertexMu2D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,1);
            var_VertexMu3D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,2);


            var_MaxD0SigSV = MaxD0SVSignificance;
            var_MinD0SigSV = MinD0SVSignificance;


            var_MinMuon_chi2LocalPosition =   std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });
            var_MaxMuon_chi2LocalPosition = std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });


            var_Muon1_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1);
            var_Muon2_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2);
            var_Muon3_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3);


            var_MinMuon_chi2LocalMomentum = std::min({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });
            var_MaxMuon_chi2LocalMomentum = std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });



            var_Muon1_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1);
            var_Muon2_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2);
            var_Muon3_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3);



            var_MintrkKink = std::min({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});
            var_MaxtrkKink = std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});
            var_MinglbKink = std::min({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});
            var_MaxglbKink = std::max({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});

            var_MuonglbkinkSum = var_sumMuTrkKinkChi2;

            var_MaxVertexPairQuality =   std::max({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)});
            var_MinVertexPairQuality =   std::min({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)});

            //---------------------
            var_NtracksClose = NcloseTracksCount;
            var_Muon1Pt = Ntp->Muon_P4(Muon_index_1).Pt();
            var_Muon2Pt = Ntp->Muon_P4(Muon_index_2).Pt();
            var_Muon3Pt = Ntp->Muon_P4(Muon_index_3).Pt();




            if(Ntp->NIsolationTrack(final_idx)!=0){ 
               var_MindcaTrackSV = sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2));
            } else var_MindcaTrackSV = -1;

            if(Ntp->NIsolationTrack(final_idx)!=0){ 
               var_dcaTrackPV = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV),2));
            } else var_dcaTrackPV = -1;



            if (id==1) MC=0;
            else  MC=1;

            if (tauMassRes<tauMassResCutLow) category = 1;
            if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh) category = 2;
            if (tauMassRes>tauMassResCutHigh) category = 3;



            TMVA_Tree->Fill();


         }
      }
   }

   void  FillTMVATrees::Finish(){
      if(mode == RECONSTRUCT){
         double scale(1.);
         double scaleDsTau(0.637);
         double scaleBpTau(0.262);
         double scaleB0Tau(0.099);
         scale = 0.8*Nminus0.at(0).at(0).Integral();
         ScaleAllHistOfType(2,scale*scaleDsTau/Nminus0.at(0).at(2).Integral());
         ScaleAllHistOfType(3,scale*scaleB0Tau/Nminus0.at(0).at(3).Integral());
         ScaleAllHistOfType(4,scale*scaleBpTau/Nminus0.at(0).at(4).Integral());
      }

      file= new TFile("FillTMVATreesInput.root","recreate");
      TMVA_Tree->SetDirectory(file);

      file->Write();
      file->Close();

      Selection::Finish();
   }
