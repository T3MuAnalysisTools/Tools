#include "ZTau3MuTaumu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTaumu::ZTau3MuTaumu(TString Name_, TString id_):
  Selection(Name_,id_),
  AnalysisName(Name_)
{
  // This is a class constructor;
}

ZTau3MuTaumu::~ZTau3MuTaumu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTaumu::Configure(){

  //  Mini tree for limit extraction


  TString treeprefix;
  if(Ntp->GetInputNtuplePath().Contains("z2tautau")) treeprefix="z2tautau";
  if(Ntp->GetInputNtuplePath().Contains("DoubleMuonLowMass")) treeprefix="DoubleMuonLowMass";


  T3MMiniTree= new TTree(treeprefix + "_" + AnalysisName,"Mini Tree Input for mva");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  
  T3MMiniTree->Branch("var_TripletPT",&var_TripletPT);
  T3MMiniTree->Branch("var_TripletEta",&var_TripletEta);
  T3MMiniTree->Branch("var_Tau3MuIsolation",&var_Tau3MuIsolation);
  T3MMiniTree->Branch("var_mu1_pT",&var_mu1_pT);
  T3MMiniTree->Branch("var_mu2_pT",&var_mu2_pT);
  T3MMiniTree->Branch("var_mu3_pT",&var_mu3_pT);
  
  T3MMiniTree->Branch("var_MuonIsolation",&var_MuonIsolation);
  T3MMiniTree->Branch("var_Muon_pT",&var_Muon_pT);
  
  T3MMiniTree->Branch("var_FLSignificance",&var_FLSignificance);
  T3MMiniTree->Branch("var_SVPVTauDirAngle",&var_SVPVTauDirAngle);
  T3MMiniTree->Branch("var_ThreeMuVertexChi2KF",&var_ThreeMuVertexChi2KF);
  T3MMiniTree->Branch("var_MinDistToIsoTrack",&var_MinDistToIsoTrack);
  T3MMiniTree->Branch("var_DeltaPhi",&var_DeltaPhi);
  
  T3MMiniTree->Branch("var_MET_Et",&var_MET_Et);
  T3MMiniTree->Branch("var_MET_Phi",&var_MET_Phi);
  
  T3MMiniTree->Branch("var_VisMass",&var_VisMass);
  T3MMiniTree->Branch("var_DiTauMass_Collinear",&var_DiTauMass_Collinear);
  T3MMiniTree->Branch("var_4Mu_Chi2",&var_4Mu_Chi2);
  
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==PassedFiducialCuts) cut.at(PassedFiducialCuts)=1;
    if(i==L1_TriggerOk)       cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)      cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==TripletPT)          cut.at(TripletPT)=20;
    if(i==OSCharge)           cut.at(OSCharge)=1;
    if(i==nMuons_PF_GL)       cut.at(nMuons_PF_GL)=1;
    if(i==nMuons_pT)          cut.at(nMuons_pT)=3.5;
    if(i==nMuons_eta)         cut.at(nMuons_eta)=2.4;
    if(i==nMuons_dR)          cut.at(nMuons_dR)=1.5;
    if(i==Tau3MuIsolation)    cut.at(Tau3MuIsolation)=0.75;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==MuonIsolation)      cut.at(MuonIsolation)=0.7;
    if(i==VisMass)            cut.at(VisMass)=1;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==L1_TriggerOk){
      title.at(i)="L1 Trigger ";
      hlabel="L1 Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT_TriggerOk){
      title.at(i)="HLT Trigger ";
      hlabel="HLT Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PassedFiducialCuts){
      title.at(i)="Passed fiducial cuts";
      hlabel="Passed fiducial cuts ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PassedFiducialCuts_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PassedFiducialCuts_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nMuons_PF_GL){
      //      title.at(i)=" At least one extra loose(PF+Gl/Tr) $\\mu$, $pT>15 GeV, |\\eta| < 2.4$ ";
      title.at(i)=" At least one extra (PF+GL) $\\mu$";
      hlabel="At least one extra muon (loose)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_PF_GL_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_PF_GL_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }
    else if(i==nMuons_pT){
      //      title.at(i)=" At least one extra loose(PF+Gl/Tr) $\\mu$, $pT>15 GeV, |\\eta| < 2.4$ ";
      title.at(i)=" At least one extra $\\mu$, $pT>3.5 GeV$";
      hlabel="pT, GeV";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_pT_",htitle,80,0.0,40,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_pT_",htitle,80,0.0,40,hlabel,"Events"));
    }
    else if(i==nMuons_eta){
      //      title.at(i)=" At least one extra loose(PF+Gl/Tr) $\\mu$, $pT>15 GeV, |\\eta| < 2.4$ ";
      title.at(i)=" At least one extra $\\mu$, $|\\eta| < 2.4$";
      hlabel="$|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
    }
    else if(i==nMuons_dR){
      //      title.at(i)=" At least one extra loose(PF+Gl/Tr) $\\mu$, $pT>15 GeV, |\\eta| < 2.4$ ";
      title.at(i)=" At least one extra $\\mu$, $\\Delta R (\\mu-3\\mu) >$ 1.5";
      hlabel="$\\Delta R (\\mu-3\\mu) $";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
    }

    else if(i==TripletPT){
      title.at(i)="pT(3$\\mu$)  $>$ 20 GeV";
      htitle=title.at(i);
      hlabel="pT(#tau_{3#mu}) , GeV ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TripletPT_",htitle,50,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TripletPT_",htitle,50,0,150,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate (3,3,2 GeV,  $|\\eta| < 2.4$, dz($\\mu_{i} , \\mu_{j}$)$<$0.5, dR($\\mu_{i} , \\mu_{j}$)$<$0.8, $\\Sigma \\mu_{charge}$ = +-1)";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)="Charge $\\mu$ * $\\tau_{3\\mu}$ = -1 ";
      title.at(i)+=" (at least one)";
      htitle=title.at(i);
      hlabel="Opposite charge ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }


    else if(i==Tau3MuIsolation){
      title.at(i)="$3\\mu $ Relative Isolation  $ > $ 0.8 ";
      //      title.at(i)+= cut.at(Tau3MuIsolation);
      htitle=title.at(i);
      hlabel="I(#tau_{3#mu})= p_{T}(#tau_{3#mu})/(p_{T}(#tau_{3#mu}) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau3MuIsolation_",htitle,50,0,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau3MuIsolation_",htitle,50,0,1.2,hlabel,"Events"));
    }


    else if(i==MuonIsolation){
      title.at(i)="Opposite $\\mu $ Relative Isolation  $ > $ 0.7";
      //      title.at(i)+= cut.at(MuonIsolation);
      htitle=title.at(i);
      hlabel="I(#mu)= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,50,0,1.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,50,0,1.1,hlabel,"Events"));
    }

    else if(i==VisMass){
      title.at(i)="25 GeV $< M(\\tau(\\mu) + \\tau(3\\mu))  < $ 85 GeV";
      htitle=title.at(i);
      hlabel="M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VisMass_",htitle,70,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VisMass_",htitle,70,0,150,hlabel,"Events"));
    }

    else if(i==TriggerMatch){
      title.at(i)="Selected (by $\\chi^2$) 3$\\mu$ matched to trg ";
      hlabel="Trigger Matched ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams
  
  PairMass_OppositeSign_dR12=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR12","PairMass_OppositeSign_dR12",40,0.2,2.,"M_{1}, GeV (OS - SS dR sorted)","Events");
  PairMass_OppositeSign_dR13=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR13","PairMass_OppositeSign_dR13",40,0.2,2.,"M_{2}, GeV (OS - SS dR sorted)","Events");




  matched_pdgId=HConfig.GetTH1D(Name+"_matched_pdgId","matched_pdgId",25,-0.5,24.5,"pdgID MC matched","Events");
  matched_dR=HConfig.GetTH1D(Name+"_matched_dR","matched_dR",50,-0.1,0.5,"#Delta R(MC-RECO) Object opposite to #tau_{3#mu}","Events");


  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{1} ","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{2} ","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{3} ","Events");
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,5,"#Delta R Truth #tau's prods ","Events");
  
  dR_betweenTruth_NeutrinoGuess=HConfig.GetTH1D(Name+"_dR_betweenTruth_NeutrinoGuess","dR_betweenTruth_NeutrinoGuess",20,0,0.5,"#Delta R: Truth to Guessed #nu","Events");
  dR_betweenTruth_Tau=HConfig.GetTH1D(Name+"_dR_betweenTruth_Tau","dR_betweenTruth_Tau",20,0,0.5,"#Delta R: Truth to Guessed #tau","Events");
  Z_Pt=HConfig.GetTH1D(Name+"_Z_Pt","Z_Pt",50,0,70,"Z_{pT}","Events");
  OS_vs_3mu_trigger=HConfig.GetTH2D(Name+"_OS_vs_3mu_trigger","OS_vs_3mu_trigger",2,-0.5,1.5,2,-0.5,1.5,"Whether 3mu Triggered","Whether OS #tau Triggered");
  
  Selection_Cut_3mu_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Pt","Selection_Cut_3mu_Pt",100,0,50.0,"3#mu p_{T}, GeV","Events");
  Selection_Cut_3mu_Rel_Iso=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Rel_Iso","Selection_Cut_3mu_Rel_Iso",50,0,1.1,"3 #mu Relative Isolation, p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","Events");
  Selection_Cut_muon_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_muon_Pt","Selection_Cut_muon_Pt",100,0,50.0,"#mu p_{T}, GeV","Events");
  Selection_Cut_muon_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_muon_Eta","Selection_Cut_muon_Eta",30,0,3.0,"#mu |#eta|","Events");
  Selection_Cut_muon_DeltaR_3mu=HConfig.GetTH1D(Name+"_Selection_Cut_muon_DeltaR_3mu","Selection_Cut_muon_DeltaR_3mu",60,0,1.2,"#Delta R (#mu-3#mu)","Events");
  Selection_Cut_muon_Rel_Iso=HConfig.GetTH1D(Name+"_Selection_Cut_muon_Rel_Iso","Selection_Cut_muon_Rel_Iso",50,0,1.1,"Opposite #mu Relative Isolation, p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T})","Events");
  Selection_Cut_Vis_InvM=HConfig.GetTH1D(Name+"_Selection_Cut_Vis_InvM","Selection_Cut_Vis_InvM",75,0,150.0,"M_{#mu + #tau(3#mu)}, GeV (visible mass)","Events");
  
  Selection_Cut_Mu1_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_before","Selection_Cut_Mu1_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after","Selection_Cut_Mu1_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after_reco","Selection_Cut_Mu1_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu2_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_before","Selection_Cut_Mu2_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after","Selection_Cut_Mu2_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after_reco","Selection_Cut_Mu2_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu3_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_before","Selection_Cut_Mu3_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after","Selection_Cut_Mu3_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after_reco","Selection_Cut_Mu3_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_mu_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_mu_p_eta_before","Selection_Cut_mu_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu p, GeV","#mu |#eta|");
  Selection_Cut_mu_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_mu_p_eta_after","Selection_Cut_mu_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu p, GeV","#mu |#eta|");
  Selection_Cut_mu_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_mu_p_eta_after_reco","Selection_Cut_mu_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu p, GeV","#mu |#eta|");
  Selection_Cut_Mu1_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_dR","Selection_Cut_Mu1_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu2_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_dR","Selection_Cut_Mu2_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu3_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_dR","Selection_Cut_Mu3_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_mu_dR=HConfig.GetTH1D(Name+"_Selection_Cut_mu_dR","Selection_Cut_mu_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu1_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu1_dR_large_scale","Selection_Cut_Mu1_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu2_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu2_dR_large_scale","Selection_Cut_Mu2_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu3_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu3_dR_large_scale","Selection_Cut_Mu3_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_mu_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_mu_dR_large_scale","Selection_Cut_mu_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_RecoMu_P=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_P","Selection_Cut_RecoMu_P",100,0.0,5.0,"#mu p, GeV","Events");
  Selection_Cut_RecoMu_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_Eta","Selection_Cut_RecoMu_Eta",50,2,3.0,"#mu |#eta|","Events");
  
  PostSelection_Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_PostSelection_Tau3MuRelativeIsolation","PostSelection_Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");InputFeatureCollection.push_back(&PostSelection_Tau3MuRelativeIsolation);
  PostSelection_OppositeMuRelativeIsolation=HConfig.GetTH1D(Name+"_PostSelection_OppositeMuRelativeIsolation","PostSelection_OppositeMuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");InputFeatureCollection.push_back(&PostSelection_OppositeMuRelativeIsolation);
  PostSelection_VisibleDiTauMass=HConfig.GetTH1D(Name+"_PostSelection_VisibleDiTauMass","PostSelection_VisibleDiTauMass",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");InputFeatureCollection.push_back(&PostSelection_VisibleDiTauMass);
  PostSelection_MTT=HConfig.GetTH1D(Name+"_PostSelection_MTT","PostSelection_MTT",70,0.,140,"M_{#tau(#mu) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  PostSelection_TripletMass=HConfig.GetTH1D(Name+"_PostSelection_TripletMass","PostSelection_TripletMass",40,1.1,2.2,"M_{3#mu}, GeV","Events");
  
  PostSelection_TripletPt=HConfig.GetTH1D(Name+"_PostSelection_TripletPt","PostSelection_TripletPt",50,2,80,"pT(3#mu), GeV ","Events");InputFeatureCollection.push_back(&PostSelection_TripletPt);
  PostSelection_OppositeMuonPt=HConfig.GetTH1D(Name+"_PostSelection_OppositeMuonPt","PostSelection_OppositeMuonPt",50,2,40,"pT(#mu), GeV ","Events");
  PostSelection_TripletEta=HConfig.GetTH1D(Name+"_PostSelection_TripletEta","PostSelection_TripletEta",50,-2.5,2.5,"#eta(3#mu)","Events");InputFeatureCollection.push_back(&PostSelection_TripletEta);
  PostSelection_OppositeMuonEta=HConfig.GetTH1D(Name+"_PostSelection_OppositeMuonEta","PostSelection_OppositeMuonEta",50,-2.5,2.5,"#eta(#mu)","Events");
  
  PostSelection_MET_Et=HConfig.GetTH1D(Name+"_PostSelection_MET_Et","PostSelection_MET_Et",100,0.0,100.0,"MET Et, GeV","Events");InputFeatureCollection.push_back(&PostSelection_MET_Et);
  PostSelection_MET_Phi=HConfig.GetTH1D(Name+"_PostSelection_MET_Phi","PostSelection_MET_Phi",20,-3.2,3.2,"MET #Phi ","Events");InputFeatureCollection.push_back(&PostSelection_MET_Phi);
  PostSelection_MET_Phi_vs_NeutrinoPhi=HConfig.GetTH2D(Name+"_PostSelection_MET_Phi_vs_NeutrinoPhi","PostSelection_MET_Phi_vs_NeutrinoPhi",40,-3.2,3.2,40,-3.2,3.2,"MET #Phi","#nu #Phi");
  PostSelection_MET_vs_NeutrinoPt=HConfig.GetTH2D(Name+"_PostSelection_MET_vs_NeutrinoPt","PostSelection_MET_vs_NeutrinoPt",50,0,100,50,0,100,"MET Et, GeV","#nu p_{T}, GeV");
  
  PostSelection_Mu1_Pt=HConfig.GetTH1D(Name+"_PostSelection_Mu1_Pt","PostSelection_Mu1_Pt",160,0.0,80.0,"#mu_{1} p, GeV","Events");InputFeatureCollection.push_back(&PostSelection_Mu1_Pt);
  PostSelection_Mu1_Eta=HConfig.GetTH1D(Name+"_PostSelection_Mu1_Eta","PostSelection_Mu1_Eta",30,0,3.14,"#mu_{1} |#eta|","Events");InputFeatureCollection.push_back(&PostSelection_Mu1_Eta);
  PostSelection_Mu2_Pt=HConfig.GetTH1D(Name+"_PostSelection_Mu2_Pt","PostSelection_Mu2_Pt",160,0.0,80.0,"#mu_{2} p, GeV","Events");InputFeatureCollection.push_back(&PostSelection_Mu2_Pt);
  PostSelection_Mu2_Eta=HConfig.GetTH1D(Name+"_PostSelection_Mu2_Eta","PostSelection_Mu2_Eta",30,0,3.14,"#mu_{2} |#eta|","Events");InputFeatureCollection.push_back(&PostSelection_Mu2_Eta);
  PostSelection_Mu3_Pt=HConfig.GetTH1D(Name+"_PostSelection_Mu3_Pt","PostSelection_Mu3_Pt",160,0.0,80.0,"#mu_{3} p, GeV","Events");InputFeatureCollection.push_back(&PostSelection_Mu3_Pt);
  PostSelection_Mu3_Eta=HConfig.GetTH1D(Name+"_PostSelection_Mu3_Eta","PostSelection_Mu3_Eta",30,0,3.14,"#mu_{3} |#eta|","Events");InputFeatureCollection.push_back(&PostSelection_Mu3_Eta);
  PostSelection_mu_Pt=HConfig.GetTH1D(Name+"_PostSelection_mu_Pt","PostSelection_mu_Pt",40,0.0,80.0,"#mu p_{T}, GeV","Events");InputFeatureCollection.push_back(&PostSelection_mu_Pt);
  PostSelection_mu_Eta=HConfig.GetTH1D(Name+"_PostSelection_mu_Eta","PostSelection_mu_Eta",30,0,3.14,"#mu |#eta|","Events");InputFeatureCollection.push_back(&PostSelection_mu_Eta);
  
  PostSelection_FLSignificance=HConfig.GetTH1D(Name+"_PostSelection_FLSignificance","PostSelection_FLSignificance",60,0,30,"PV - SV distance  significance","Events"); InputFeatureCollection.push_back(&PostSelection_FLSignificance);
  PostSelection_SVPVTauDirAngle=HConfig.GetTH1D(Name+"_PostSelection_SVPVTauDirAngle","PostSelection_SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{3#mu}, rad","Events"); InputFeatureCollection.push_back(&PostSelection_SVPVTauDirAngle);
  PostSelection_SVPVTauDirAngle_largescale=HConfig.GetTH1D(Name+"_PostSelection_SVPVTauDirAngle_largescale","PostSelection_SVPVTauDirAngle_largescale",50,-3.2,3.2,"Angle btw #vec{SV}-#vec{PV} and #vec{3#mu}, rad","Events");
  PostSelection_VertexChi2KF=HConfig.GetTH1D(Name+"_PostSelection_VertexChi2KF","PostSelection_VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events"); InputFeatureCollection.push_back(&PostSelection_VertexChi2KF);
  PostSelection_MinDistToIsoTrack=HConfig.GetTH1D(Name+"_PostSelection_MinDistToIsoTrack","PostSelection_MinDistToIsoTrack",100,0,0.1,"Min dR To IsoTrack","Events");InputFeatureCollection.push_back(&PostSelection_MinDistToIsoTrack);
  PostSelection_Kinematics_MissingTrMass=HConfig.GetTH1D(Name+"_PostSelection_Kinematics_MissingTrMass","PostSelection_Kinematics_MissingTrMass",100,0,100.,"M_{T}, GeV","Events");
  PostSelection_Kinematics_MissingTrMass_cos=HConfig.GetTH1D(Name+"_PostSelection_Kinematics_MissingTrMass_cos","PostSelection_Kinematics_MissingTrMass_cos",100,-1,1.,"cos(#theta)","Events");InputFeatureCollection.push_back(&PostSelection_Kinematics_MissingTrMass_cos);
  PostSelection_VisibleDiTauMass_Collinear=HConfig.GetTH1D(Name+"_PostSelection_VisibleDiTauMass_Collinear","PostSelection_VisibleDiTauMass_Collinear",70,30.,180,"M_{#tau(#mu) + #tau(3#mu) + #nu}, GeV","Events");InputFeatureCollection.push_back(&PostSelection_VisibleDiTauMass_Collinear);
  
  PostSelection_prod_size=HConfig.GetTH1D(Name+"PostSelection__prod_size","PostSelection_prod_size",7,-0.5,6.5,"no. of visible products","Events");InputFeatureCollection.push_back(&PostSelection_prod_size);
  
  
  PostSelection_Vertex_Dist=HConfig.GetTH1D(Name+"_PostSelection_Vertex_Dist","PostSelection_Vertex_Dist",100,0.0,1.0,"Vertex shift, cm","Events");InputFeatureCollection.push_back(&PostSelection_Vertex_Dist);
  PostSelection_Vertex_Chi2=HConfig.GetTH1D(Name+"_PostSelection_Vertex_Chi2","PostSelection_Vertex_Chi2",50,0,20,"Vertex #chi^{2}","Events");InputFeatureCollection.push_back(&PostSelection_Vertex_Chi2);
  
  PostSelection_VisibleDiTauMass_peakspectra1=HConfig.GetTH1D(Name+"_PostSelection_VisibleDiTauMass_peakspectra1","PostSelection_VisibleDiTauMass_peakspectra1",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");InputFeatureCollection.push_back(&PostSelection_VisibleDiTauMass_peakspectra1);
  PostSelection_VisibleDiTauMass_peakspectra2=HConfig.GetTH1D(Name+"_PostSelection_VisibleDiTauMass_peakspectra2","PostSelection_VisibleDiTauMass_peakspectra2",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");InputFeatureCollection.push_back(&PostSelection_VisibleDiTauMass_peakspectra2);

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTaumu::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  
  Extradist1d.push_back(&matched_pdgId);
  Extradist1d.push_back(&matched_dR);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);
  Extradist1d.push_back(&dR_betweenTruth_VisibleTaus);

  Extradist1d.push_back(&PairMass_OppositeSign_dR12);
  Extradist1d.push_back(&PairMass_OppositeSign_dR13);
  
  Extradist1d.push_back(&dR_betweenTruth_NeutrinoGuess);
  Extradist1d.push_back(&dR_betweenTruth_Tau);
  Extradist1d.push_back(&Z_Pt);
  Extradist2d.push_back(&OS_vs_3mu_trigger);
  
  Extradist1d.push_back(&Selection_Cut_3mu_Pt);
  Extradist1d.push_back(&Selection_Cut_3mu_Rel_Iso);
  Extradist1d.push_back(&Selection_Cut_muon_Pt);
  Extradist1d.push_back(&Selection_Cut_muon_Eta);
  Extradist1d.push_back(&Selection_Cut_muon_DeltaR_3mu);
  Extradist1d.push_back(&Selection_Cut_muon_Rel_Iso);
  Extradist1d.push_back(&Selection_Cut_Vis_InvM);
  
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_mu_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_mu_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_mu_p_eta_after_reco);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR);
  Extradist1d.push_back(&Selection_Cut_mu_dR);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_mu_dR_large_scale);
  
  Extradist1d.push_back(&Selection_Cut_RecoMu_P);
  Extradist1d.push_back(&Selection_Cut_RecoMu_Eta);
  
  Extradist1d.push_back(&PostSelection_Tau3MuRelativeIsolation);
  Extradist1d.push_back(&PostSelection_OppositeMuRelativeIsolation);
  Extradist1d.push_back(&PostSelection_VisibleDiTauMass);
  Extradist1d.push_back(&PostSelection_MTT);
  Extradist1d.push_back(&PostSelection_TripletMass);
  
  Extradist1d.push_back(&PostSelection_TripletPt);
  Extradist1d.push_back(&PostSelection_OppositeMuonPt);
  Extradist1d.push_back(&PostSelection_TripletEta);
  Extradist1d.push_back(&PostSelection_OppositeMuonEta);
  
  Extradist1d.push_back(&PostSelection_MET_Et);
  Extradist1d.push_back(&PostSelection_MET_Phi);
  Extradist2d.push_back(&PostSelection_MET_Phi_vs_NeutrinoPhi);
  Extradist2d.push_back(&PostSelection_MET_vs_NeutrinoPt);
  
  Extradist1d.push_back(&PostSelection_Mu1_Pt);
  Extradist1d.push_back(&PostSelection_Mu1_Eta);
  Extradist1d.push_back(&PostSelection_Mu2_Pt);
  Extradist1d.push_back(&PostSelection_Mu2_Eta);
  Extradist1d.push_back(&PostSelection_Mu3_Pt);
  Extradist1d.push_back(&PostSelection_Mu3_Eta);
  Extradist1d.push_back(&PostSelection_mu_Pt);
  Extradist1d.push_back(&PostSelection_mu_Eta);
  
  Extradist1d.push_back(&PostSelection_FLSignificance);
  Extradist1d.push_back(&PostSelection_SVPVTauDirAngle);
  Extradist1d.push_back(&PostSelection_SVPVTauDirAngle_largescale);
  Extradist1d.push_back(&PostSelection_VertexChi2KF);
  Extradist1d.push_back(&PostSelection_MinDistToIsoTrack);
  Extradist1d.push_back(&PostSelection_Kinematics_MissingTrMass);
  Extradist1d.push_back(&PostSelection_Kinematics_MissingTrMass_cos);
  Extradist1d.push_back(&PostSelection_VisibleDiTauMass_Collinear);
  
  Extradist1d.push_back(&PostSelection_prod_size);
  
  Extradist1d.push_back(&PostSelection_Vertex_Dist);
  Extradist1d.push_back(&PostSelection_Vertex_Chi2);
  Extradist1d.push_back(&PostSelection_VisibleDiTauMass_peakspectra1);
  Extradist1d.push_back(&PostSelection_VisibleDiTauMass_peakspectra2);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTaumu::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  //  std::cout<<"   id   "<< id << std::endl;

  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu0Fired(false);
  bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);
  bool HLT_OppositeSide(false);
  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //    if(Ntp->HLTDecision(iTrigger))    std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
    if(HLTName.Contains("HLT_IsoMu24") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
    if(HLTName.Contains("HLT_IsoMu27") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
  }
  
  OS_vs_3mu_trigger.at(t).Fill(HLTOk,HLT_OppositeSide,1 );
  
  random_num = rndm.Rndm();
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    // if(Ntp->L1Decision(il1))    std::cout<<" l1 name  "<< Ntp->L1Name(il1) <<  " =   "<<Ntp->L1Decision(il1) <<std::endl;
    //    std::cout<<" random number   "<< random_num << std::endl;
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }
  value.at(L1_TriggerOk)=0;
  value.at(HLT_TriggerOk)=0;
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;

  value.at(L1_TriggerOk)=(L1Ok);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  
  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));



  value.at(SignalCandidate) = Ntp->NThreeMuons();

  int  signal_idx=-1;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  TLorentzVector Muon_LV;
  TLorentzVector MC_NeutrinoSum_LV(0.,0.,0.,0.);
  int Whether_decay_found(0);
  
  bool WhetherSignalMC = id==210||id==210231||id==210232||id==210233;
  if(WhetherSignalMC){
  
  //int Whether_decay_found(0);
  int TausFromZ_Count(0);
  int tau_3mu_idx(-1);
  int tau_mu_idx(-1);
  int tau_e_idx(-1);
  int tau_h_idx(-1);
  for(unsigned int imc =0; imc< Ntp->NMCParticles(); imc++){
    if(abs(Ntp->MCParticle_pdgid(imc)) == 15 && Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(imc) ) == 23){
      TausFromZ_Count++;
          
          //correcting for cases where tau decays multiple times
          int ChildIdx = imc;
          bool Whether_Another_Tau_Decay_In_Chain(false);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15) Whether_Another_Tau_Decay_In_Chain=true;
          }
          while (Whether_Another_Tau_Decay_In_Chain){
            Whether_Another_Tau_Decay_In_Chain=false;
            for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
              if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15){
                ChildIdx = Ntp->MCParticle_childidx(ChildIdx).at(i);
                Whether_Another_Tau_Decay_In_Chain=true;
              }
            }
          }
          
          int nmuons_temp(0);
          int nelectrons_temp(0);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            //std::cout<<"Child "<< i <<" of tau is :  "<< Ntp->MCParticle_childpdgid(ChildIdx).at(i) << std::endl;
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==13){
              nmuons_temp++;
            }
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==11){
              nelectrons_temp++;
            }
          }
          
          if(nmuons_temp==3){
            tau_3mu_idx=ChildIdx;
          }
          else if(nmuons_temp==1){
            tau_mu_idx=ChildIdx;
          }
          else if(nelectrons_temp==1){
            tau_e_idx=ChildIdx;
          }
          else{
            tau_h_idx=ChildIdx;
          }
          
          
      }
  }
    
  if(tau_3mu_idx>-0.5 && tau_mu_idx>-0.5) Whether_decay_found=1;
  int Whether_decay_found_temp(0);
  if(tau_3mu_idx>-0.5) Whether_decay_found_temp=1;
  
  //std::cout<<"  No of taus from Z:  "<< TausFromZ_Count << std::endl;
  //if(Whether_decay_found&&id==210233)      Ntp->printMCDecayChainOfEvent(true,true,true,true);//&&id!=210
  
  TLorentzVector Mu1_LV;
  TLorentzVector Mu2_LV;
  TLorentzVector Mu3_LV;
  //TLorentzVector Muon_LV;
  
  if(Whether_decay_found==1){
    std::vector<int> Sorted_MC_Indices = Ntp->SortedPtMuons_MC(Ntp->MCParticle_childidx(tau_3mu_idx));
    
    Mu1_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(0));
    Mu2_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(1));
    Mu3_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(2));
    
    for(int i = 0; i < Ntp->MCParticle_childpdgid(tau_mu_idx).size(); i++){
      if(abs(Ntp->MCParticle_childpdgid(tau_mu_idx).at(i))==13){
        Muon_LV=Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_mu_idx).at(i));
      }
      if(abs(Ntp->MCParticle_childpdgid(tau_mu_idx).at(i))==12||abs(Ntp->MCParticle_childpdgid(tau_mu_idx).at(i))==14||abs(Ntp->MCParticle_childpdgid(tau_mu_idx).at(i))==16){
        MC_NeutrinoSum_LV=MC_NeutrinoSum_LV+Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_mu_idx).at(i));
      }
    }
  }
  
  double cut_Mu_Candidate_p=2.49;
  double cut_Mu_Candidate_eta=2.41;
  double cut_Tau_mu_Candidate_p=2.49;
  double cut_Tau_mu_Candidate_eta=2.41;
  
  bool var_Whether_decay_found = Whether_decay_found;
  bool var_Mu1_Candidate_p = (Mu1_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu1_Candidate_eta = (abs(Mu1_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Mu2_Candidate_p = (Mu2_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu2_Candidate_eta = (abs(Mu2_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Mu3_Candidate_p = (Mu3_LV.Vect().Mag()) > cut_Mu_Candidate_p;
  bool var_Mu3_Candidate_eta = (abs(Mu3_LV.Eta())) < cut_Mu_Candidate_eta;
  bool var_Tau_mu_Candidate_p = (Muon_LV.Vect().Mag()) > cut_Tau_mu_Candidate_p;
  bool var_Tau_mu_Candidate_eta = (abs(Muon_LV.Eta())) < cut_Tau_mu_Candidate_eta;

  double dR1_max(99.0);
  double dR2_max(99.0);
  double dR3_max(99.0);
  double dR4_max(99.0);
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      TLorentzVector T1=Ntp->Muon_P4(imu);
      
      if(Ntp->Muon_P4(imu).DeltaR(Mu1_LV)<dR1_max){
        dR1_max=Ntp->Muon_P4(imu).DeltaR(Mu1_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu2_LV)<dR2_max){
        dR2_max=Ntp->Muon_P4(imu).DeltaR(Mu2_LV);
      }

      if(Ntp->Muon_P4(imu).DeltaR(Mu3_LV)<dR3_max){
        dR3_max=Ntp->Muon_P4(imu).DeltaR(Mu3_LV);
      }


      if(Ntp->Muon_P4(imu).DeltaR(Muon_LV)<dR4_max){
        dR4_max=Ntp->Muon_P4(imu).DeltaR(Muon_LV);
      }
      
      Selection_Cut_RecoMu_P.at(t).Fill(Ntp->Muon_P4(imu).Vect().Mag(),1 );
      Selection_Cut_RecoMu_Eta.at(t).Fill(abs(Ntp->Muon_P4(imu).Eta()),1 );

    }
  
  bool var_Mu1_Candidate_recod = (dR1_max<0.01);
  bool var_Mu2_Candidate_recod = (dR2_max<0.01);
  bool var_Mu3_Candidate_recod = (dR3_max<0.01);
  bool var_Tau_mu_Candidate_recod = (dR4_max<0.01);
  
  value.at(PassedFiducialCuts)=var_Whether_decay_found&&var_Mu1_Candidate_p&&var_Mu1_Candidate_eta&&var_Mu2_Candidate_p&&var_Mu2_Candidate_eta&&var_Mu3_Candidate_p&&var_Mu3_Candidate_eta&&var_Tau_mu_Candidate_p&&var_Tau_mu_Candidate_eta&&var_Mu1_Candidate_recod&&var_Mu2_Candidate_recod&&var_Mu3_Candidate_recod&&var_Tau_mu_Candidate_recod;
  pass.at(PassedFiducialCuts)=(value.at(PassedFiducialCuts)==cut.at(PassedFiducialCuts));
  
  // This is to print out selected event content
  if(id==210232){
          
          //std::cout<<"Whether tau e decay found: "<< Whether_decay_found << ". Event id: "<<id<< std::endl;
          
          if(Whether_decay_found){
                  /*
                  std::cout<<"------------------------------- "<< std::endl;
                  std::cout<<"Event Content "<< std::endl;
                  std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Mu1_LV) << std::endl;
                  std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Mu2_LV) << std::endl;
                  std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Mu3_LV) << std::endl;
                  std::cout<<" idx4:  "<<Ntp->getMatchTruthIndex(Muon_LV) << std::endl;
                  Ntp->printMCDecayChainOfEvent(true, true, true, true);
                  std::cout<< "\n\n\n\n\n\n";
                  
                  for(unsigned int i=0; i < Ntp->NMCTaus(); i++){
                    for(unsigned int j=0; j < Ntp->NMCTauDecayProducts(i); j++){
                      std::cout<<"Product coming from index: "<<i<<" with product index: "<<j<<" : "<<Ntp->MCTauandProd_pdgid(i, j)<< std::endl;
                    }
                  }
                  */
                  /*
                  Selection_Cut_Mu1_P.at(t).Fill(Mu1_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu1_Eta.at(t).Fill(abs(Mu1_LV.Eta()),1 );
                  Selection_Cut_Mu2_P.at(t).Fill(Mu2_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu2_Eta.at(t).Fill(abs(Mu2_LV.Eta()),1 );
                  Selection_Cut_Mu3_P.at(t).Fill(Mu3_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu3_Eta.at(t).Fill(abs(Mu3_LV.Eta()),1 );
                  Selection_Cut_mu_P.at(t).Fill(Muon_LV.Vect().Mag(),1 );
                  Selection_Cut_mu_Eta.at(t).Fill(abs(Muon_LV.Eta()),1 );
                  */
                  /*
                  Selection_Cut_Mu1_p_eta_before.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                  if(pass.at(Mu1_Candidate_p)&&pass.at(Mu1_Candidate_eta)){
                    Selection_Cut_Mu1_p_eta_after.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    if(!pass.at(Mu1_Candidate_recod)){
                      Selection_Cut_Mu1_p_eta_after_reco.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_Mu2_p_eta_before.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                  if(pass.at(Mu2_Candidate_p)&&pass.at(Mu2_Candidate_eta)){
                    Selection_Cut_Mu2_p_eta_after.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    if(!pass.at(Mu2_Candidate_recod)){
                      Selection_Cut_Mu2_p_eta_after_reco.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu2_Candidate_eta)){
                  Selection_Cut_Mu3_p_eta_before.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                  if(pass.at(Mu3_Candidate_p)&&pass.at(Mu3_Candidate_eta)){
                    Selection_Cut_Mu3_p_eta_after.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    if(!pass.at(Mu3_Candidate_recod)){
                      Selection_Cut_Mu3_p_eta_after_reco.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu3_Candidate_eta)){
                  Selection_Cut_mu_p_eta_before.at(t).Fill(Muon_LV.Vect().Mag(),abs(Muon_LV.Eta()));
                  if(pass.at(Tau_mu_Candidate_p)&&pass.at(Tau_mu_Candidate_eta)){
                    Selection_Cut_mu_p_eta_after.at(t).Fill(Muon_LV.Vect().Mag(),abs(Muon_LV.Eta()));
                    if(!pass.at(Tau_mu_Candidate_recod)){
                      Selection_Cut_mu_p_eta_after_reco.at(t).Fill(Muon_LV.Vect().Mag(),abs(Muon_LV.Eta()));
                    }
                  }
                  }
                  */
                  Selection_Cut_Mu1_dR.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR.at(t).Fill(dR3_max,1 );
                  Selection_Cut_mu_dR.at(t).Fill(dR4_max,1 );
                  
                  Selection_Cut_Mu1_dR_large_scale.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR_large_scale.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR_large_scale.at(t).Fill(dR3_max,1 );
                  Selection_Cut_mu_dR_large_scale.at(t).Fill(dR4_max,1 );
                  
          }
          
  }
  
  }//if(id!=1)
  
  if(!WhetherSignalMC){
    value.at(PassedFiducialCuts)=1;
    pass.at(PassedFiducialCuts)=1;
  }
  

  TLorentzVector Tau3MuLV(0,0,0,0);
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));


  value.at(TripletPT)=0;
  if(signal_idx!=-1)
    {
      Tau3MuLV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
      value.at(TripletPT) = Tau3MuLV.Pt();
      Selection_Cut_3mu_Pt.at(t).Fill(value.at(TripletPT));

    }

  pass.at(TripletPT) = ( value.at(TripletPT) >= cut.at(TripletPT) );
  //  std::cout<<"   index  "<< signal_idx << "  id  "<< id << std::endl;

  std::vector<int> Muons;
  std::vector<int> Muons_OppositeHemisphere;
  std::vector<int> Muons_OppositeHemisphere_PF_GL;
  std::vector<int> Muons_OppositeHemisphere_pT;
  std::vector<int> Muons_OppositeHemisphere_eta;
  std::vector<int> Muons_OppositeHemisphere_dR;
  std::vector<std::vector<double>> Muons_OppositeHemisphere_OppositeCharge;




  value.at(nMuons_PF_GL)  = -1;
  value.at(nMuons_pT)  = -1;
  value.at(nMuons_eta)  = 99.0;
  value.at(nMuons_dR)  = -1;
  
  bool highest_PF_GL(0);//highest value determines whether it passes (when there are multiple muons)
  double highest_pT(-1.0);
  double lowest_eta(10.0);
  double highest_dR(-1.0);
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      if(signal_idx!=-1)
	{
           bool WhetherMuonDifferentFrom3Mu(true);
           for (int i=0; i<3; i++){
             TLorentzVector Muon_Temp_LV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(i));
             double dR_temp = Muon_Temp_LV.DeltaR(Ntp->Muon_P4(imu));
             //std::cout<<"dR: "<< dR_temp << std::endl;
             if(dR_temp<0.001){
               WhetherMuonDifferentFrom3Mu=false;
             }
           }
           
           if(WhetherMuonDifferentFrom3Mu){
                   
                   //value.at(nMuons_PF_GL)  = Ntp->Muon_isPFMuon(imu) && Ntp->Muon_isGlobalMuon(imu);
                   //value.at(nMuons_pT)  = Ntp->Muon_P4(imu).Pt();
                   //value.at(nMuons_eta)  = fabs(Ntp->Muon_P4(imu).Eta());
                   //value.at(nMuons_dR)  = Ntp->Muon_P4(imu).DeltaR(Tau3MuLV);
                   
                   if(Ntp->Muon_isPFMuon(imu) && Ntp->Muon_isGlobalMuon(imu)  ){
                           highest_PF_GL=true;
                           Muons_OppositeHemisphere_PF_GL.push_back(imu);
                           
                           if(Ntp->Muon_P4(imu).DeltaR(Tau3MuLV)>highest_dR){
                             highest_dR=Ntp->Muon_P4(imu).DeltaR(Tau3MuLV);
                           }
                           
                           if(Ntp->Muon_P4(imu).DeltaR(Tau3MuLV) > cut.at(nMuons_dR)  ){
                                   Muons_OppositeHemisphere_dR.push_back(imu);
                                   
                                   if(Ntp->Muon_P4(imu).Pt() > highest_pT){
                                     highest_pT=Ntp->Muon_P4(imu).Pt();
                                   }
                                   
                                   if(Ntp->Muon_P4(imu).Pt() > cut.at(nMuons_pT)  ){
                                           Muons_OppositeHemisphere_pT.push_back(imu);
                                           
                                           if(fabs(Ntp->Muon_P4(imu).Eta()) < lowest_eta){
                                             lowest_eta=fabs(Ntp->Muon_P4(imu).Eta());
                                           }
                                           
                                           if(fabs(Ntp->Muon_P4(imu).Eta()) < cut.at(nMuons_eta)  ){
                                                   Muons_OppositeHemisphere_eta.push_back(imu);
                                           }
                                   }
                           }
                   }
                   
                   
                   if(Ntp->Muon_P4(imu).Pt() > cut.at(nMuons_pT) && fabs(Ntp->Muon_P4(imu).Eta()) < cut.at(nMuons_eta) && Ntp->Muon_isPFMuon(imu) && Ntp->Muon_isGlobalMuon(imu) &&
                   Ntp->Muon_P4(imu).DeltaR(Tau3MuLV) > cut.at(nMuons_dR)  )Muons_OppositeHemisphere.push_back(imu);
                   
                   
           }
                   
           
	}
    }

  
  value.at(nMuons_PF_GL)  = highest_PF_GL;
  value.at(nMuons_pT)  = highest_pT;
  value.at(nMuons_eta)  = lowest_eta;
  value.at(nMuons_dR)  = highest_dR;
  
  
  Selection_Cut_muon_Pt.at(t).Fill(highest_pT);
  Selection_Cut_muon_Eta.at(t).Fill(lowest_eta);
  Selection_Cut_muon_DeltaR_3mu.at(t).Fill(highest_dR);
  
  
  pass.at(nMuons_PF_GL) = ( Muons_OppositeHemisphere_PF_GL.size() > 0 );
  pass.at(nMuons_pT)    = ( Muons_OppositeHemisphere_pT.size() > 0 );
  pass.at(nMuons_eta)   = ( Muons_OppositeHemisphere_eta.size() > 0 );
  pass.at(nMuons_dR)    = ( Muons_OppositeHemisphere_dR.size() > 0 );


  value.at(OSCharge)        =  0;
  value.at(Tau3MuIsolation) = -1;
  value.at(MuonIsolation)   = -1;
  value.at(VisMass)         = -1;

  value.at(TriggerMatch) = 0;
  if(signal_idx!=-1)
    {

      int index_mu_1 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);  
      int index_mu_2 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);  
      int index_mu_3 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);  

      TLorentzVector TripletmuLV = Ntp->Muon_P4(index_mu_1) +  Ntp->Muon_P4(index_mu_2) +  Ntp->Muon_P4(index_mu_3);
      
      value.at(Tau3MuIsolation) = TripletmuLV.Pt()/  (TripletmuLV.Pt()  + Ntp->Muon_RelIso(index_mu_1) +
 						                          Ntp->Muon_RelIso(index_mu_2) +
						                          Ntp->Muon_RelIso(index_mu_3) );
                                                                                                    
      Selection_Cut_3mu_Rel_Iso.at(t).Fill(value.at(Tau3MuIsolation));
    
    
    //Trigger Matching

      int triggerCheck = 0.1;
      if(pass.at(HLT_TriggerOk))
	{
	  vector<TLorentzVector> trigobjTriplet;
	  for (int i=0; i<Ntp->NTriggerObjects(); i++)
	    {
	      TString name = Ntp->TriggerObject_name(i);
	      //	if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
	      TLorentzVector tmp;
	      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
	      trigobjTriplet.push_back(tmp);
	    }
	  std::vector<TLorentzVector> muonTriplet;
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)));
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)));
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)));
	  //	  std::cout<<" trigobjTriplet.size()  "<< trigobjTriplet.size() << std::endl;
	  if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
	}

        value.at(TriggerMatch) = triggerCheck;      

      
	if( Muons_OppositeHemisphere.size()>0)
	  {
	    for(auto i : Muons_OppositeHemisphere)
	      {
		
		int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
		  Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
		  Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));
		
		if(Ntp->Muon_charge(i)*Tau3MuCharge == -1) Muons_OppositeHemisphere_OppositeCharge.push_back({Ntp->Muon_P4(i).DeltaR(Tau3MuLV),i});
                
	      }
	    value.at(OSCharge) = Muons_OppositeHemisphere_OppositeCharge.size();
	  }
	
	
    }
  
    pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));
    //    std::cout<<" (value.at(TriggerMatch)   " << value.at(TriggerMatch)  << " cut.at(TriggerMatch)   " << cut.at(TriggerMatch) << "  ??   " 
    //	     <<(value.at(TriggerMatch)  ==  cut.at(TriggerMatch) )<< std::endl;
    pass.at(OSCharge) = (value.at(OSCharge) >= cut.at(OSCharge)); 



    double highest_MuonIsolation(-1.0);//highest value determines whether it passes (when there are multiple muons)
    double central_VisMass(200.0);
    
    if(pass.at(OSCharge))
      {
      
      for(unsigned int i = 0 ; i< Muons_OppositeHemisphere_OppositeCharge.size(); i++){
           
           int muon_index = Muons_OppositeHemisphere_OppositeCharge[i][1];
           
           double MuonIsolation_val = Ntp->Muon_P4(muon_index).Pt()  /  (Ntp->Muon_P4(muon_index).Pt()  +   Ntp->Muon_RelIso(muon_index) );
           double VisMass_val = (Tau3MuLV + Ntp->Muon_P4(muon_index)).M();
           if(MuonIsolation_val>highest_MuonIsolation){
             highest_MuonIsolation=MuonIsolation_val;
             
             if(fabs(VisMass_val-60.0)   < fabs(central_VisMass-60.0)  ){
               central_VisMass=VisMass_val;
             }
             
           }
      }
      
    }
    
    value.at(MuonIsolation)   = highest_MuonIsolation;
    value.at(VisMass) = central_VisMass;

    pass.at(Tau3MuIsolation) = (value.at(Tau3MuIsolation) > cut.at(Tau3MuIsolation));
    pass.at(MuonIsolation)   = (value.at(MuonIsolation)   > cut.at(MuonIsolation));
    pass.at(VisMass)         = (value.at(VisMass) > 25 && value.at(VisMass) < 85);
    //    std::cout<<"  ----  "<< std::endl;
    /*    for(unsigned int i = 0 ; i< NCuts; i++)
      {

	std::cout<<"  i "<<i<< "    "<<pass.at(i) << std::endl;
	}*/
    

    /*    if(i==WhetherDecayFound)        cut.at(WhetherDecayFound)=1;
    if(i==Mu1_Candidate_p)          cut.at(Mu1_Candidate_p)=2.49;

    if(passAllUntil(Mu1_Candidate_p)) 
      Eta_muon1_beforecuts.at(t).Fill();

    if(passAllUntil(Mu2_Candidate_p)) 
      Eta_muon1_beforecuts.at(t).Fill();
    */
    
    pass.at(TripletPT) = 1;
    pass.at(nMuons_pT) = 1;
    pass.at(MuonIsolation) = 1;
    pass.at(Tau3MuIsolation) = 1;
    pass.at(VisMass) = 1;

    double wobs=1;
    double w;  
             
    if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
    else{w=1;}
    
    std::vector<unsigned int> exclude_cuts;
    exclude_cuts.push_back(VisMass);
    if(passAllBut(exclude_cuts)) Selection_Cut_Vis_InvM.at(t).Fill(value.at(VisMass));
    Selection_Cut_muon_Rel_Iso.at(t).Fill(highest_MuonIsolation);
  

  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 
    
    //    std::cout<<"pass nMuons "<< pass.at(nMuons) <<    "  triplet pT    "<< pass.at(TripletPT)  << std::endl;
    //    std::cout<<"   how many muons i have  " << Muons_OppositeHemisphere_OppositeCharge .size() << std::endl;
    sort( Muons_OppositeHemisphere_OppositeCharge.rbegin(), Muons_OppositeHemisphere_OppositeCharge.rend() ); // sort based on first column, highest first
    unsigned int muon_idx = Muons_OppositeHemisphere_OppositeCharge[0][1];
    TLorentzVector MuLV = Ntp->Muon_P4(muon_idx);
    
    PostSelection_prod_size.at(t).Fill(Muons_OppositeHemisphere_OppositeCharge.size());


    unsigned int muon_1_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int muon_2_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int muon_3_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);
    
    TLorentzVector Tau3muLV = Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));
    ////////////////////////   sort muons by charge and dR and fill pair masses :
    /////
    vector<unsigned int> idx_vec;
    idx_vec.push_back(muon_1_idx);
    idx_vec.push_back(muon_2_idx);
    idx_vec.push_back(muon_3_idx);
    
    PostSelection_Mu1_Pt.at(t).Fill(Ntp->Muon_P4(muon_1_idx).Pt(),1 );
    PostSelection_Mu1_Eta.at(t).Fill(abs(Ntp->Muon_P4(muon_1_idx).Eta()),1 );
    PostSelection_Mu2_Pt.at(t).Fill(Ntp->Muon_P4(muon_2_idx).Pt(),1 );
    PostSelection_Mu2_Eta.at(t).Fill(abs(Ntp->Muon_P4(muon_2_idx).Eta()),1 );
    PostSelection_Mu3_Pt.at(t).Fill(Ntp->Muon_P4(muon_3_idx).Pt(),1 );
    PostSelection_Mu3_Eta.at(t).Fill(abs(Ntp->Muon_P4(muon_3_idx).Eta()),1 );
    
    PostSelection_mu_Pt.at(t).Fill(highest_pT,1 );
    PostSelection_mu_Eta.at(t).Fill(lowest_eta,1 );
    
    //Primary Vertex
    double val_FLSignificance=Ntp->FlightLength_significance(Ntp->Vertex_HighestPt_PrimaryVertex(),Ntp->Vertex_HighestPt_PrimaryVertex_Covariance(),
	   							Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx));
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_HighestPt_PrimaryVertex());
    double val_SVPVTauDirAngle=SVPV.Angle(Tau3muLV.Vect());
    double val_ThreeMuVertexChi2KF=Ntp->Vertex_signal_KF_Chi2(signal_idx);
    double val_MinDistToIsoTrack=Ntp->Isolation_MinDist(signal_idx);
    double val_DeltaPhi=(1-TMath::Cos(Ntp->METPhi()-(MuLV.Vect()).Phi()));
    
    PostSelection_FLSignificance.at(t).Fill(val_FLSignificance);
    PostSelection_VertexChi2KF.at(t).Fill(val_ThreeMuVertexChi2KF);
    PostSelection_SVPVTauDirAngle.at(t).Fill(val_SVPVTauDirAngle);
    PostSelection_SVPVTauDirAngle_largescale.at(t).Fill(val_SVPVTauDirAngle);
    PostSelection_MinDistToIsoTrack.at(t).Fill(val_MinDistToIsoTrack);
    
    // Missing transverse mass
    PostSelection_Kinematics_MissingTrMass.at(t).Fill(sqrt(   2*Ntp->METEt()*TMath::Sqrt(MuLV.Px()*MuLV.Px()+MuLV.Py()*MuLV.Py())*(1-TMath::Cos(Ntp->METPhi()-(MuLV.Vect()).Phi()))   )); //use definition transverse mass for 2 particles
    PostSelection_Kinematics_MissingTrMass_cos.at(t).Fill(  val_DeltaPhi   );
    
    TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(MuLV.Theta()));
    TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());
    double val_DiTauMass_Collinear=(MuLV + Tau3muLV + Neutrino_LV).M();
    PostSelection_VisibleDiTauMass_Collinear.at(t).Fill((MuLV + Tau3muLV + Neutrino_LV).M(), 1);
    
    // Common vertex
    TVector3 FirstGuess(0.1,0.1,0.1);
    std::vector<TrackParticle> Triplet_set;
    Triplet_set.push_back(Ntp->Muon_TrackParticle(muon_1_idx));
    Triplet_set.push_back(Ntp->Muon_TrackParticle(muon_2_idx));
    Triplet_set.push_back(Ntp->Muon_TrackParticle(muon_3_idx));
    Chi2VertexFitter  Triplet_FittedVertex(Triplet_set,Ntp->Vertex_Signal_KF_pos(signal_idx));
    Triplet_FittedVertex.Fit();
    
    std::vector<TrackParticle> Triplet_plus_mu_set;
    Triplet_plus_mu_set=Triplet_set;
    Triplet_plus_mu_set.push_back(Ntp->Muon_TrackParticle(muon_idx));
    Chi2VertexFitter  Four_Particle_FittedVertex(Triplet_plus_mu_set,Triplet_FittedVertex.GetVertex());
    Four_Particle_FittedVertex.Fit();
    
    //std::cout << "Triplet_FittedVertex mag: " << Triplet_FittedVertex.GetVertex().Mag() << std::endl;
    //std::cout << "Four_Particle_FittedVertex mag: " << Four_Particle_FittedVertex.GetVertex().Mag() << std::endl;
    
    double val_4Mu_Chi2=Four_Particle_FittedVertex.ChiSquare();
    PostSelection_Vertex_Dist.at(t).Fill(fabs((Triplet_FittedVertex.GetVertex()-Four_Particle_FittedVertex.GetVertex()).Mag()),1 );
    PostSelection_Vertex_Chi2.at(t).Fill(val_4Mu_Chi2,1 );
    

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    TLorentzVector MuonOS  = Ntp->Muon_P4(os_mu_idx);  
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);



    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){

      PairMass_OppositeSign_dR12.at(t).Fill((MuonOS+MuonSS2).M(),1 );
      PairMass_OppositeSign_dR13.at(t).Fill((MuonOS+MuonSS1).M(),1 );


    }else{
      
      PairMass_OppositeSign_dR12.at(t).Fill((MuonOS+MuonSS1).M(),1 );
      PairMass_OppositeSign_dR13.at(t).Fill((MuonOS+MuonSS2).M(),1 );
      
    }
    //////
    ///////////////////////////
    
    double val_MET_Et=Ntp->METEt();
    double val_MET_Phi=Ntp->METPhi();
    
    PostSelection_MET_Et.at(t).Fill(val_MET_Et );
    PostSelection_MET_Phi.at(t).Fill( val_MET_Phi );
    PostSelection_MET_Phi_vs_NeutrinoPhi.at(t).Fill( Ntp->METPhi(),(MC_NeutrinoSum_LV.Vect()).Phi() );
    PostSelection_MET_vs_NeutrinoPt.at(t).Fill( Ntp->METEt(),MC_NeutrinoSum_LV.Pt() );
    


    TLorentzVector Muon1LV = Ntp->Muon_P4(muon_1_idx);
    TLorentzVector Muon2LV = Ntp->Muon_P4(muon_2_idx);
    TLorentzVector Muon3LV = Ntp->Muon_P4(muon_3_idx);

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);


    



    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    



    float RelativeIsolationMu1 = Ntp->Muon_RelIso(muon_1_idx);
    float RelativeIsolationMu2 = Ntp->Muon_RelIso(muon_2_idx);
    float RelativeIsolationMu3 = Ntp->Muon_RelIso(muon_3_idx);
    //////////// kinematics 

    PostSelection_TripletPt.at(t).Fill(Tau3muLV.Pt(),1);
    PostSelection_OppositeMuonPt.at(t).Fill(MuLV.Pt(),1);
    PostSelection_TripletEta.at(t).Fill(Tau3muLV.Eta(),1);
    PostSelection_OppositeMuonEta.at(t).Fill(MuLV.Eta(),1);
    //////////// kinematics 


    PostSelection_Tau3MuRelativeIsolation.at(t).Fill(    Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt()),1);
    PostSelection_OppositeMuRelativeIsolation.at(t).Fill(    Ntp->Muon_P4(muon_idx).Pt()/(Ntp->Muon_RelIso(muon_idx) +  Ntp->Muon_P4(muon_idx).Pt()),1);

    PostSelection_VisibleDiTauMass.at(t).Fill((MuLV + Tau3muLV).M(), 1);
    PostSelection_MTT.at(t).Fill( (Tau3muLV + MuLV  + Neutrino_LV).M(), 1);
    
    if(Ntp->FlightLength_significance(Ntp->Vertex_HighestPt_PrimaryVertex(),Ntp->Vertex_HighestPt_PrimaryVertex_Covariance(),
	   							Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx))<2.0){
      PostSelection_VisibleDiTauMass_peakspectra1.at(t).Fill((MuLV + Tau3muLV).M(), 1);
    }
    if(fabs((Triplet_FittedVertex.GetVertex()-Four_Particle_FittedVertex.GetVertex()).Mag())<0.02){
      PostSelection_VisibleDiTauMass_peakspectra2.at(t).Fill((MuLV + Tau3muLV).M(), 1);
    }



    bool PlotMCOnly(false);  // and blind for data
    if(id!=1) PlotMCOnly = true;
    if(id==1 && ( (TauRefitLV.M() > 1.1 && TauRefitLV.M() < 1.7233333) or (TauRefitLV.M() > 1.8333333 && TauRefitLV.M() < 2.2)) ) PlotMCOnly=true;
    if(PlotMCOnly)  PostSelection_TripletMass.at(t).Fill(TauRefitLV.M(),1);
    

    TLorentzVector OppositeSideLV = MuLV;
    TLorentzVector Neutrino_LV_Guess_Result;//Guessed Neutrino LV
    if(false)
      {

        
        matched_pdgId.at(t).Fill(Ntp->matchTruth(OppositeSideLV),1);
        matched_dR.at(t).Fill(Ntp->matchToTruthTauDecay(OppositeSideLV).DeltaR(OppositeSideLV),1);

	TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;

	Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
	dR_betweenTruth_VisibleTaus.at(t).Fill(MCTauLV.DeltaR(Ntp->matchToTruthTauDecay(OppositeSideLV)),1);
        
        //////// Getting the Neutrino MC
        int MC_Muon_idx = Ntp->getMatchTruthIndex(MuLV);//Getting the index of Muon from Spectator Tau in MC
        int Tau_MC_idx = Ntp->MCParticle_midx(MC_Muon_idx);//Getting the index of Spectator Tau in MC
        
        TLorentzVector MuonMC;
        TLorentzVector TauMC;
        TLorentzVector NeutrinoMC;
        bool FoundMCTau(false);
        if(Tau_MC_idx>=0){ // Tau_MC_idx is -1 for some events?
                TLorentzVector MuonMC=Ntp->MCParticle_p4(MC_Muon_idx);
                TLorentzVector TauMC=Ntp->MCParticle_p4(Tau_MC_idx);
                TLorentzVector NeutrinoMC=TauMC-MuonMC;
                if(abs(Ntp->MCParticle_pdgid(Tau_MC_idx))==15){
                        FoundMCTau = true;
                }
        }
        if(FoundMCTau){// Neutrinos satisfying tau invariant mass can be created with any direction of neutrino, PROVING THERE ARE MULTIPLE SOLUTIONS (without taking Z into account). We can construct a neutrino LV consistent with tau and Z mass. The same thing will happen with data, so we can't use this method to distinguish between data and MC? (Or, is it?)
                
                //std::cout << "Found mu neutrinos: " << std::endl;
                double Max_Diff(999.0);
                int Division_No_Thet(1);//200
                int Division_No_Phi(1);//1000
                
                double Phi_Visible_Tau = MuonMC.Phi();
                double Theta_Visible_Tau = MuonMC.Theta();
                TVector3 Visible_Tau_Vec = MuonMC.Vect();
                
                for(int Theta_1_idx=0; Theta_1_idx<Division_No_Thet;Theta_1_idx++){
                        for(int Phi_1_idx=0; Phi_1_idx<Division_No_Phi;Phi_1_idx++){
                                
                                double Theta_1 = 0.0 + Theta_1_idx*(0.5)/Division_No_Thet;
                                double Phi_1 = 0.0 + Phi_1_idx*(2*TMath::Pi())/Division_No_Phi;
                                
                                TVector3 Neutrino_Vect_Guess(TMath::Sin(Theta_1)*TMath::Cos(Phi_1),TMath::Sin(Theta_1)*TMath::Sin(Phi_1),TMath::Cos(Theta_1));
                                if(Phi_Visible_Tau >= TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau-2*TMath::Pi();
                                if(Phi_Visible_Tau < -TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau+2*TMath::Pi();
                                Neutrino_Vect_Guess.RotateY(Theta_Visible_Tau);
                                Neutrino_Vect_Guess.RotateZ(Phi_Visible_Tau);
                                
                                double scaling_x = (1.77686*1.77686 - MuonMC.M()*MuonMC.M())/(2*(MuonMC.E()-MuonMC.Vect().Mag()*TMath::Cos(Neutrino_Vect_Guess.Angle(Visible_Tau_Vec))));//Imposing tau mass criteria
                                TLorentzVector Neutrino_LV_Guess(scaling_x*Neutrino_Vect_Guess,scaling_x);
                                
                                // Boosting to COM frame
                                //TLorentzVector MuonMC_Boosted_Guess = MuonMC;
                                //MuonMC_Boosted_Guess.Boost(-1*(MuonMC+Neutrino_LV_Guess).BoostVector());
                                //TLorentzVector Neutrino_LV_Boosted_Guess = Neutrino_LV_Guess;
                                //Neutrino_LV_Boosted_Guess.Boost(-1*(MuonMC+Neutrino_LV_Guess).BoostVector());
                                
                                double Diff_Val = abs((MuonMC+Neutrino_LV_Guess+Tau3muLV).M() - 91.1876);// Imposing Z mass
                                if(Diff_Val < Max_Diff){
                                        Max_Diff = Diff_Val;
                                        Neutrino_LV_Guess_Result = Neutrino_LV_Guess;
                                }
        
                        }
                }
                
                //std::cout << "Angle from guess result to actual: " << Neutrino_LV_Guess_Result.Vect().Angle(NeutrinoMC.Vect()) << " with Max_Diff: " << Max_Diff << std::endl;
                
                dR_betweenTruth_NeutrinoGuess.at(t).Fill(Neutrino_LV_Guess_Result.DeltaR(NeutrinoMC),1);
                dR_betweenTruth_Tau.at(t).Fill((Neutrino_LV_Guess_Result+MuLV).DeltaR(TauMC),1);
                Z_Pt.at(t).Fill((TauMC+Tau3muLV).Pt(),1);
                
                TVector3 MC_MET = (TauMC-NeutrinoMC+Tau3muLV).Vect();
                MC_MET.SetZ(0.);
                MC_MET = -1 * MC_MET;
                
                //std::cout << "MC Delta Phi from actual to MET: " << MC_MET.Phi()-NeutrinoMC.Vect().Phi() << std::endl;
                //std::cout << "Z angle: " << (TauMC+Tau3muLV).Vect().Theta() << " Z (Pt): "<< (TauMC+Tau3muLV).Pt() << " Tau vis to Tau3mu angle: " << (TauMC-NeutrinoMC).Vect().Phi() - Tau3muLV.Vect().Phi() << " Tau to Tau3mu angle: " << (TauMC).Vect().Phi() - Tau3muLV.Vect().Phi() << std::endl;
                
        }//FoundMCNeutrino
        
        
      }
      
      
        //*** fill up the T3MMiniTree.root for statistical analysis
        
        m3m = Tau3muLV.M();
        
        dataMCtype = id;
        event_weight =1; // 1 for data
        if(dataMCtype == 1){event_weight =1;}
        else if(dataMCtype == 210233){event_weight =5.00e-04;}
        else if(dataMCtype == 210232){event_weight =4.93e-04;}
        else if(dataMCtype == 210231){event_weight =4.87e-04;}
        
        m12 = (MuonOS+MuonSS1).M();
        m13 = (MuonOS+MuonSS2).M();
        
        var_TripletPT=Tau3muLV.Pt();
        var_TripletEta=Tau3muLV.Eta();
        var_Tau3MuIsolation=Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt());
        var_mu1_pT=Muon1LV.Pt();
        var_mu2_pT=Muon2LV.Pt();
        var_mu3_pT=Muon3LV.Pt();
        
        var_MuonIsolation=Ntp->Muon_P4(muon_idx).Pt()/(Ntp->Muon_RelIso(muon_idx) +  Ntp->Muon_P4(muon_idx).Pt());
        var_Muon_pT=MuLV.Pt();
        
        var_FLSignificance=val_FLSignificance;
        var_SVPVTauDirAngle=val_SVPVTauDirAngle;
        var_ThreeMuVertexChi2KF=val_ThreeMuVertexChi2KF;
        var_MinDistToIsoTrack=val_MinDistToIsoTrack;
        var_DeltaPhi=val_DeltaPhi;
        
        var_MET_Et=val_MET_Et;
        var_MET_Phi=val_MET_Phi;
        
        var_DiTauMass_Collinear=val_DiTauMass_Collinear;
        var_VisMass=(Tau3muLV + MuLV ).M();
        var_4Mu_Chi2=val_4Mu_Chi2;
        
        T3MMiniTree->Fill();
  
  }
}


void  ZTau3MuTaumu::Finish(){

  if (mode==RECONSTRUCT){
    double lumi_scale_1_taue(0.000293); //need to be entered manually
    double lumi_scale_2_taumu(0.000283);
    double lumi_scale_3_tauh(0.000280);
    for ( unsigned int j=0; j<InputFeatureCollection.size(); ++j){
      double scale(1.0);
      if(InputFeatureCollection.at(j)->size()>=4){
        if(InputFeatureCollection.at(j)->at(2).Integral()!=0) scale = InputFeatureCollection.at(j)->at(0).Integral()/(InputFeatureCollection.at(j)->at(2).Integral()*lumi_scale_2_taumu);
        InputFeatureCollection.at(j)->at(2).Scale(scale);
      }
    }
  }
  
  /*
  if(mode == RECONSTRUCT){
      double scale(1.0);
      if(Nminus0.at(0).at(3).Integral()!=0) scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();//Gives raw event no as Integral. Incorrect
      
      ScaleAllHistOfType(2,scale);
  }
  */

  //*** write down the T3MMiniTree.root for statistical analysis
  TString out_file_name  = "MVA_Mini_Tree_"+ AnalysisName+".root";
  T3MFMiniTree = new TFile(out_file_name,"recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();


  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





