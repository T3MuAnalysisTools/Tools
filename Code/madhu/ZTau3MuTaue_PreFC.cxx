#include "ZTau3MuTaue_PreFC.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTaue_PreFC::ZTau3MuTaue_PreFC(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTaue_PreFC::~ZTau3MuTaue_PreFC(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTaue_PreFC::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==WhetherDecayFound)        cut.at(WhetherDecayFound)=1;
    if(i==Mu1_Candidate_p)          cut.at(Mu1_Candidate_p)=2.49;
    if(i==Mu1_Candidate_eta)        cut.at(Mu1_Candidate_eta)=2.41;
    if(i==Mu2_Candidate_p)          cut.at(Mu2_Candidate_p)=2.49;
    if(i==Mu2_Candidate_eta)        cut.at(Mu2_Candidate_eta)=2.41;
    if(i==Mu3_Candidate_p)          cut.at(Mu3_Candidate_p)=2.49;
    if(i==Mu3_Candidate_eta)        cut.at(Mu3_Candidate_eta)=2.41;
    if(i==Tau_e_Candidate_p)        cut.at(Tau_e_Candidate_p)=4.9;
    if(i==Tau_e_Candidate_eta)      cut.at(Tau_e_Candidate_eta)=2.41;
    if(i==Mu1_Candidate_recod)      cut.at(Mu1_Candidate_recod)=1;
    if(i==Mu2_Candidate_recod)      cut.at(Mu2_Candidate_recod)=1;
    if(i==Mu3_Candidate_recod)      cut.at(Mu3_Candidate_recod)=1;
    if(i==Tau_e_Candidate_recod)    cut.at(Tau_e_Candidate_recod)=1;
    if(i==L1_TriggerOk)       cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)      cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==OSCharge)           cut.at(OSCharge)=1;
    if(i==nElectrons_PF_cut)  cut.at(nElectrons_PF_cut)=1;
    if(i==nElectrons_pT)      cut.at(nElectrons_pT)=5;
    if(i==nElectrons_eta)     cut.at(nElectrons_eta)=2.41;
    if(i==nElectrons_dR)      cut.at(nElectrons_dR)=1.8;
    if(i==ElectronIsolation)  cut.at(ElectronIsolation)=1;
    if(i==Tau3MuIsolation)    cut.at(Tau3MuIsolation)=0.5;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==TripletPT)          cut.at(TripletPT)=20;
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
    else if(i==WhetherDecayFound){
      title.at(i)="$Z \\rightarrow \\tau_{e}, \\tau_{3\\mu}$ decay information found in ntuple";
      hlabel="Decay information found in ntuple ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_WhetherDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_WhetherDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{1}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{1}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{1}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{1}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{2}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{2}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{2}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{1}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{3}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{3}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{3}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{3}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Tau_e_Candidate_p){
      title.at(i)=" Whether GEN level $e$ has $pT>4.9 GeV$ ";
      hlabel="$\\tau_{e}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_e_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_e_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Tau_e_Candidate_eta){
      title.at(i)=" Whether GEN level $e$ has $|\\eta| < 2.41$ ";
      hlabel="$\\tau_{e}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_e_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_e_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_recod){
      title.at(i)=" Whether $\\mu_{1}$ is reconstructed in MC ";
      hlabel="If $\\mu_{1}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_recod){
      title.at(i)=" Whether $\\mu_{2}$ is reconstructed in MC ";
      hlabel="If $\\mu_{2}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_recod){
      title.at(i)=" Whether $\\mu_{3}$ is reconstructed in MC ";
      hlabel="If $\\mu_{3}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau_e_Candidate_recod){
      title.at(i)=" Whether $e$ is reconstructed in MC ";
      hlabel="If $e$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_e_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_e_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nElectrons_PF_cut){
      title.at(i)=" At least 1 e, PF + cutBasedElectronID ( loose WP )";
      hlabel="number of electrons";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nElectrons_PF_cut_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nElectrons_PF_cut_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }
    else if(i==nElectrons_pT){
      title.at(i)=" At least 1 e, $pT> 6 GeV$ ";
      hlabel="pT, GeV";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nElectrons_pT_",htitle,80,0.0,40,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nElectrons_pT_",htitle,80,0.0,40,hlabel,"Events"));
    }
    else if(i==nElectrons_eta){
      title.at(i)=" At least 1 e, $|\\eta| < 2.4$ ";
      hlabel="$|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nElectrons_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nElectrons_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
    }
    else if(i==nElectrons_dR){
      title.at(i)=" At least 1 e, $\\Delta R (e-3\\mu) >$ 1.8 ";
      hlabel=" $\\Delta R (e-3\\mu) $";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nElectrons_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nElectrons_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
    }
    else if(i==TripletPT){
      title.at(i)="pT(3$\\mu$)  $>$ 20 GeV";
      htitle=title.at(i);
      hlabel="pT(#tau_{3#mu}) , GeV ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TripletPT_",htitle,50,5,80,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TripletPT_",htitle,50,5,80,hlabel,"Events"));
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
      title.at(i)="Charge e * $\\tau_{3\\mu}$ =  -1";
      title.at(i)+=" (at least one)";
      htitle=title.at(i);
      hlabel="Opposite charge? ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }
    else if(i==Tau3MuIsolation){
      title.at(i)="$3\\mu $ Relative Isolation  $ > $ 0.5 ";
      //      title.at(i)+= cut.at(Tau3MuIsolation);
      htitle=title.at(i);
      hlabel="I(3#mu)= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
    }
    else if(i==ElectronIsolation){
      title.at(i)="Electron Isolation ";
      //      title.at(i)+= cut.at(Tau3MuIsolation);
      htitle=title.at(i);
      hlabel="Electron Isolation Sum (dummy for now)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ElectronIsolation_",htitle,25,0,20.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ElectronIsolation_",htitle,25,0,20.1,hlabel,"Events"));
    }

    else if(i==VisMass){
      title.at(i)="35 GeV $< M(\\tau(\\mu) + \\tau(3\\mu))  < $ 90 GeV";
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




  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  ElectronSumIsolation=HConfig.GetTH1D(Name+"_ElectronSumIsolation","ElectronSumIsolation",50,0.,10,"I= neutralH + chargedH + photon Iso, GeV","Events ");
  
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",70,0.,140,"M_{#tau(#mu) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",40,1.1,2.2,"M_{3#mu}, GeV","Events");
  PairMass_OppositeSign_dR12=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR12","PairMass_OppositeSign_dR12",40,0.2,2.,"M_{1}, GeV (OS - SS dR sorted)","Events");
  PairMass_OppositeSign_dR13=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR13","PairMass_OppositeSign_dR13",40,0.2,2.,"M_{2}, GeV (OS - SS dR sorted)","Events");





  matched_pdgId=HConfig.GetTH1D(Name+"_matched_pdgId","matched_pdgId",25,-0.5,24.5,"pdgID MC matched","Events");
  matched_dR=HConfig.GetTH1D(Name+"_matched_dR","matched_dR",50,-0.1,0.5,"#Delta R(MC-RECO) Object opposite to #tau_{3#mu}","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{1} ","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{2} ","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{3} ","Events");
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,0.5,"#Delta R Truth #tau's prods ","Events");
  
  dR_betweenTruth_NeutrinoGuess=HConfig.GetTH1D(Name+"_dR_betweenTruth_NeutrinoGuess","dR_betweenTruth_NeutrinoGuess",20,0,0.5,"#Delta R: Truth to Guessed #nu","Events");
  dR_betweenTruth_Tau=HConfig.GetTH1D(Name+"_dR_betweenTruth_Tau","dR_betweenTruth_Tau",20,0,0.5,"#Delta R: Truth to Guessed #tau","Events");
  Z_Pt=HConfig.GetTH1D(Name+"_Z_Pt","Z_Pt",50,0,70,"Z_{pT}","Events");
  OS_vs_3mu_trigger=HConfig.GetTH2D(Name+"_OS_vs_3mu_trigger","OS_vs_3mu_trigger",2,-0.5,1.5,2,-0.5,1.5,"Whether 3mu Triggered","Whether OS #tau Triggered");
  
  Selection_Cut_3mu_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Pt","Selection_Cut_3mu_Pt",100,0,50.0,"3#mu p_{T}, GeV","Events");
  Selection_Cut_3mu_Rel_Iso=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Rel_Iso","Selection_Cut_3mu_Rel_Iso",50,0,1.1,"3 #mu Relative Isolation, p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","Events");
  Selection_Cut_elect_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_elect_Pt","Selection_Cut_elect_Pt",100,0,50.0,"e p_{T}, GeV","Events");
  Selection_Cut_elect_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_elect_Eta","Selection_Cut_elect_Eta",30,0,3.0,"e |#eta|","Events");
  Selection_Cut_elect_DeltaR_3mu=HConfig.GetTH1D(Name+"_Selection_Cut_elect_DeltaR_3mu","Selection_Cut_elect_DeltaR_3mu",60,0,1.2,"#Delta R (e-3#mu)","Events");
  Selection_Cut_Vis_InvM=HConfig.GetTH1D(Name+"_Selection_Cut_Vis_InvM","Selection_Cut_Vis_InvM",75,0,150.0,"M_{e + #tau(3#mu)}, GeV (visible mass)","Events");
  
  Selection_Cut_Mu1_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_P","Selection_Cut_Mu1_P",160,0.0,80.0,"#mu_{1} p, GeV","Events");
  Selection_Cut_Mu1_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_Eta","Selection_Cut_Mu1_Eta",30,0,3.14,"#mu_{1} |#eta|","Events");
  Selection_Cut_Mu1_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_before","Selection_Cut_Mu1_p_eta_before",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after","Selection_Cut_Mu1_p_eta_after",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after_reco","Selection_Cut_Mu1_p_eta_after_reco",40,0.0,80.0,30,0,3.14,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  
  Selection_Cut_Mu2_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_P","Selection_Cut_Mu2_P",160,0.0,80.0,"#mu_{2} p, GeV","Events");
  Selection_Cut_Mu2_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_Eta","Selection_Cut_Mu2_Eta",30,0,3.14,"#mu_{2} |#eta|","Events");
  Selection_Cut_Mu2_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_before","Selection_Cut_Mu2_p_eta_before",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after","Selection_Cut_Mu2_p_eta_after",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after_reco","Selection_Cut_Mu2_p_eta_after_reco",40,0.0,80.0,30,0,3.14,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  
  Selection_Cut_Mu3_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_P","Selection_Cut_Mu3_P",160,0.0,80.0,"#mu_{3} p, GeV","Events");
  Selection_Cut_Mu3_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_Eta","Selection_Cut_Mu3_Eta",30,0,3.14,"#mu_{3} |#eta|","Events");
  Selection_Cut_Mu3_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_before","Selection_Cut_Mu3_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after","Selection_Cut_Mu3_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after_reco","Selection_Cut_Mu3_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  
  Selection_Cut_El_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_El_Pt","Selection_Cut_El_Pt",40,0.0,80.0,"e p_{T}, GeV","Events");
  Selection_Cut_El_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_El_Eta","Selection_Cut_El_Eta",30,0,3.14,"e |#eta|","Events");
  Selection_Cut_El_pt_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_before","Selection_Cut_El_pt_eta_before",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  Selection_Cut_El_pt_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_after","Selection_Cut_El_pt_eta_after",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  Selection_Cut_El_pt_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_El_pt_eta_after_reco","Selection_Cut_El_pt_eta_after_reco",200,0.0,100.0,100,0,5.0,"e pT, GeV","e |#eta|");
  
  Selection_Cut_Mu1_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_dR","Selection_Cut_Mu1_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu2_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_dR","Selection_Cut_Mu2_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu3_dR=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_dR","Selection_Cut_Mu3_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_El_dR=HConfig.GetTH1D(Name+"_Selection_Cut_El_dR","Selection_Cut_El_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu1_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu1_dR_large_scale","Selection_Cut_Mu1_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu2_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu2_dR_large_scale","Selection_Cut_Mu2_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu3_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu3_dR_large_scale","Selection_Cut_Mu3_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_El_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_El_dR_large_scale","Selection_Cut_El_dR_large_scale",100,0,2.0,"#Delta R","Events");
  
  Selection_Cut_RecoMu_P=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_P","Selection_Cut_RecoMu_P",100,0.0,5.0,"#mu p, GeV","Events");
  Selection_Cut_RecoMu_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_Eta","Selection_Cut_RecoMu_Eta",50,2,3.0,"#mu |#eta|","Events");
  Selection_Cut_RecoEl_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_RecoEl_Pt","Selection_Cut_RecoEl_Pt",100,0.0,10.0,"e p_{T}, GeV","Events");
  Selection_Cut_RecoEl_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoEl_Eta","Selection_Cut_RecoEl_Eta",50,2,3.0,"e |#eta|","Events");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTaue_PreFC::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output


  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&ElectronSumIsolation);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);
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
  Extradist1d.push_back(&Selection_Cut_elect_Pt);
  Extradist1d.push_back(&Selection_Cut_elect_Eta);
  Extradist1d.push_back(&Selection_Cut_elect_DeltaR_3mu);
  Extradist1d.push_back(&Selection_Cut_Vis_InvM);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_P);
  Extradist1d.push_back(&Selection_Cut_Mu1_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu2_P);
  Extradist1d.push_back(&Selection_Cut_Mu2_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu3_P);
  Extradist1d.push_back(&Selection_Cut_Mu3_Eta);
  Extradist1d.push_back(&Selection_Cut_El_Pt);
  Extradist1d.push_back(&Selection_Cut_El_Eta);
  
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_before);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_after);
  Extradist2d.push_back(&Selection_Cut_El_pt_eta_after_reco);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR);
  Extradist1d.push_back(&Selection_Cut_El_dR);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_El_dR_large_scale);
  
  Extradist1d.push_back(&Selection_Cut_RecoMu_P);
  Extradist1d.push_back(&Selection_Cut_RecoMu_Eta);
  Extradist1d.push_back(&Selection_Cut_RecoEl_Pt);
  Extradist1d.push_back(&Selection_Cut_RecoEl_Eta);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTaue_PreFC::doEvent(){ 

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
  bool HLT_OppositeSide(false);
  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
    if(HLTName.Contains("HLT_Ele32_WPTight_Gsf") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
    if(HLTName.Contains("HLT_Ele35_WPTight_Gsf") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
  }
  
  OS_vs_3mu_trigger.at(t).Fill(HLTOk,HLT_OppositeSide,1 );
  
  random_num = rndm.Rndm();
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }
  value.at(L1_TriggerOk)=0;  
  value.at(HLT_TriggerOk)= 0 ;
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
  
  bool WhetherSignalMC = id==210||id==210231||id==210232||id==210233;
  if(WhetherSignalMC){
  
  int Whether_decay_found(0);
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
    
  if(tau_3mu_idx>-0.5 && tau_e_idx>-0.5) Whether_decay_found=1;
  int Whether_decay_found_temp(0);
  if(tau_3mu_idx>-0.5) Whether_decay_found_temp=1;
  
  //std::cout<<"  No of taus from Z:  "<< TausFromZ_Count << std::endl;
  //if(Whether_decay_found&&id==210233)      Ntp->printMCDecayChainOfEvent(true,true,true,true);
  
  TLorentzVector Mu1_LV;
  TLorentzVector Mu2_LV;
  TLorentzVector Mu3_LV;
  TLorentzVector Electron_LV;
  if(Whether_decay_found==1){
    std::vector<int> Sorted_MC_Indices = Ntp->SortedPtMuons_MC(Ntp->MCParticle_childidx(tau_3mu_idx));
    
    Mu1_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(0));
    Mu2_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(1));
    Mu3_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(2));
    
    for(int i = 0; i < Ntp->MCParticle_childpdgid(tau_e_idx).size(); i++){
      if(abs(Ntp->MCParticle_childpdgid(tau_e_idx).at(i))==11){
        Electron_LV=Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_e_idx).at(i));
      }
    }
  }
  
  value.at(WhetherDecayFound)=(Whether_decay_found);
  pass.at(WhetherDecayFound)=(value.at(WhetherDecayFound)==cut.at(WhetherDecayFound));
  
  value.at(Mu1_Candidate_p)=(Mu1_LV.Vect().Mag());
  pass.at(Mu1_Candidate_p)=value.at(Mu1_Candidate_p)>cut.at(Mu1_Candidate_p);
  
  value.at(Mu1_Candidate_eta)=(abs(Mu1_LV.Eta()));
  pass.at(Mu1_Candidate_eta)=value.at(Mu1_Candidate_eta)<cut.at(Mu1_Candidate_eta);
  
  value.at(Mu2_Candidate_p)=(Mu2_LV.Vect().Mag());
  pass.at(Mu2_Candidate_p)=value.at(Mu2_Candidate_p)>cut.at(Mu2_Candidate_p);
  
  value.at(Mu2_Candidate_eta)=(abs(Mu2_LV.Eta()));
  pass.at(Mu2_Candidate_eta)=value.at(Mu2_Candidate_eta)<cut.at(Mu2_Candidate_eta);
  
  value.at(Mu3_Candidate_p)=(Mu3_LV.Vect().Mag());
  pass.at(Mu3_Candidate_p)=value.at(Mu3_Candidate_p)>cut.at(Mu3_Candidate_p);
  
  value.at(Mu3_Candidate_eta)=(abs(Mu3_LV.Eta()));
  pass.at(Mu3_Candidate_eta)=value.at(Mu3_Candidate_eta)<cut.at(Mu3_Candidate_eta);
  
  value.at(Tau_e_Candidate_p)=(Electron_LV.Pt());
  pass.at(Tau_e_Candidate_p)=value.at(Tau_e_Candidate_p)>cut.at(Tau_e_Candidate_p);
  
  value.at(Tau_e_Candidate_eta)=(abs(Electron_LV.Eta()));
  pass.at(Tau_e_Candidate_eta)=value.at(Tau_e_Candidate_eta)<cut.at(Tau_e_Candidate_eta);
  
  
  double dR1_max(99.0);
  double dR2_max(99.0);
  double dR3_max(99.0);
  double dR4_max(99.0);
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      if(Ntp->Muon_P4(imu).DeltaR(Mu1_LV)<dR1_max){
        dR1_max=Ntp->Muon_P4(imu).DeltaR(Mu1_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu2_LV)<dR2_max){
        dR2_max=Ntp->Muon_P4(imu).DeltaR(Mu2_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu3_LV)<dR3_max){
        dR3_max=Ntp->Muon_P4(imu).DeltaR(Mu3_LV);
      }
      
      Selection_Cut_RecoMu_P.at(t).Fill(Ntp->Muon_P4(imu).Vect().Mag(),1 );
      Selection_Cut_RecoMu_Eta.at(t).Fill(abs(Ntp->Muon_P4(imu).Eta()),1 );

    }
  
  for(unsigned int ie=0; ie < Ntp->NElectrons(); ie++)
    {
      if(Ntp->Electron_P4(ie).DeltaR(Electron_LV)<dR4_max){
        dR4_max=Ntp->Electron_P4(ie).DeltaR(Electron_LV);
      }
      
      Selection_Cut_RecoEl_Pt.at(t).Fill(Ntp->Electron_P4(ie).Vect().Mag(),1 );
      Selection_Cut_RecoEl_Eta.at(t).Fill(abs(Ntp->Electron_P4(ie).Eta()),1 );
      
    }
  
  value.at(Mu1_Candidate_recod)=(dR1_max<0.01);
  pass.at(Mu1_Candidate_recod)=value.at(Mu1_Candidate_recod);
  value.at(Mu2_Candidate_recod)=(dR2_max<0.01);
  pass.at(Mu2_Candidate_recod)=value.at(Mu2_Candidate_recod);
  value.at(Mu3_Candidate_recod)=(dR3_max<0.01);
  pass.at(Mu3_Candidate_recod)=value.at(Mu3_Candidate_recod);
  
  value.at(Tau_e_Candidate_recod)=(dR4_max<0.01);
  pass.at(Tau_e_Candidate_recod)=value.at(Tau_e_Candidate_recod);
  
  // This is to print out selected event content
  if(id==210231){
          
          //std::cout<<"Whether tau e decay found: "<< Whether_decay_found << ". Event id: "<<id<< std::endl;
          
          if(Whether_decay_found){
                  /*
                  std::cout<<"------------------------------- "<< std::endl;
                  std::cout<<"Event Content "<< std::endl;
                  std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Mu1_LV) << std::endl;
                  std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Mu2_LV) << std::endl;
                  std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Mu3_LV) << std::endl;
                  std::cout<<" idx4:  "<<Ntp->getMatchTruthIndex(Electron_LV) << std::endl;
                  Ntp->printMCDecayChainOfEvent(true, true, true, true);
                  std::cout<< "\n\n\n\n\n\n";
                  
                  for(unsigned int i=0; i < Ntp->NMCTaus(); i++){
                    for(unsigned int j=0; j < Ntp->NMCTauDecayProducts(i); j++){
                      std::cout<<"Product coming from index: "<<i<<" with product index: "<<j<<" : "<<Ntp->MCTauandProd_pdgid(i, j)<< std::endl;
                    }
                  }
                  */
                  Selection_Cut_Mu1_P.at(t).Fill(Mu1_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu1_Eta.at(t).Fill(abs(Mu1_LV.Eta()),1 );
                  Selection_Cut_Mu2_P.at(t).Fill(Mu2_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu2_Eta.at(t).Fill(abs(Mu2_LV.Eta()),1 );
                  Selection_Cut_Mu3_P.at(t).Fill(Mu3_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu3_Eta.at(t).Fill(abs(Mu3_LV.Eta()),1 );
                  Selection_Cut_El_Pt.at(t).Fill(Electron_LV.Pt(),1 );
                  Selection_Cut_El_Eta.at(t).Fill(abs(Electron_LV.Eta()),1 );
                  
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
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_Mu3_p_eta_before.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                  if(pass.at(Mu3_Candidate_p)&&pass.at(Mu3_Candidate_eta)){
                    Selection_Cut_Mu3_p_eta_after.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    if(!pass.at(Mu3_Candidate_recod)){
                      Selection_Cut_Mu3_p_eta_after_reco.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_El_pt_eta_before.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                  if(pass.at(Tau_e_Candidate_p)&&pass.at(Tau_e_Candidate_eta)){
                    Selection_Cut_El_pt_eta_after.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                    if(!pass.at(Tau_e_Candidate_recod)){
                      Selection_Cut_El_pt_eta_after_reco.at(t).Fill(Electron_LV.Pt(),abs(Electron_LV.Eta()));
                    }
                  }
                  }
                  
                  Selection_Cut_Mu1_dR.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR.at(t).Fill(dR3_max,1 );
                  Selection_Cut_El_dR.at(t).Fill(dR4_max,1 );
                  
                  Selection_Cut_Mu1_dR_large_scale.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR_large_scale.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR_large_scale.at(t).Fill(dR3_max,1 );
                  Selection_Cut_El_dR_large_scale.at(t).Fill(dR4_max,1 );
                  
          }
          
  }
  
  }//if(id!=1)
  
  if(!WhetherSignalMC){
    pass.at(WhetherDecayFound)=1;
    pass.at(Mu1_Candidate_p)=1;
    pass.at(Mu1_Candidate_eta)=1;
    pass.at(Mu2_Candidate_p)=1;
    pass.at(Mu2_Candidate_eta)=1;
    pass.at(Mu3_Candidate_p)=1;
    pass.at(Mu3_Candidate_eta)=1;
    pass.at(Tau_e_Candidate_p)=1;
    pass.at(Tau_e_Candidate_eta)=1;
    pass.at(Mu1_Candidate_recod)=1;
    pass.at(Mu2_Candidate_recod)=1;
    pass.at(Mu3_Candidate_recod)=1;
    pass.at(Tau_e_Candidate_recod)=1;
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



  std::vector<int> Electrons;
  std::vector<int> Electrons_OppositeHemisphere;
  std::vector<int> Electrons_OppositeHemisphere_PF_cut;
  std::vector<int> Electrons_OppositeHemisphere_pT;
  std::vector<int> Electrons_OppositeHemisphere_eta;
  std::vector<int> Electrons_OppositeHemisphere_dR;
  std::vector<int> Electrons_OppositeHemisphere_OppositeCharge;




  value.at(nElectrons_PF_cut)  = -1;
  value.at(nElectrons_pT)  = -1;
  value.at(nElectrons_eta)  = 99.0;
  value.at(nElectrons_dR)  = -1;
  for(unsigned int ie=0; ie < Ntp->NElectrons(); ie++)
    {
      
      //value.at(nElectrons_PF_cut)  = Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_loose(ie);
      //value.at(nElectrons_pT)  = Ntp->Electron_P4(ie).Pt();
      //value.at(nElectrons_eta)  = fabs(Ntp->Electron_P4(ie).Eta());
      //value.at(nElectrons_dR)  = Ntp->Electron_P4(ie).DeltaR(Tau3MuLV);
      
      if(Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_loose(ie) ) Electrons_OppositeHemisphere_PF_cut.push_back(ie);
      if(Ntp->Electron_P4(ie).Pt() > cut.at(nElectrons_pT)   ) Electrons_OppositeHemisphere_pT.push_back(ie);
      if(fabs(Ntp->Electron_P4(ie).Eta()) < cut.at(nElectrons_eta)   ) Electrons_OppositeHemisphere_eta.push_back(ie);
      if(Ntp->Electron_P4(ie).DeltaR(Tau3MuLV) > cut.at(nElectrons_dR)   ) Electrons_OppositeHemisphere_dR.push_back(ie);
      
      if(Ntp->Electron_P4(ie).Pt() > cut.at(nElectrons_pT) && fabs(Ntp->Electron_P4(ie).Eta()) < cut.at(nElectrons_eta) && 
	 Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_loose(ie) &&
	 Ntp->Electron_P4(ie).DeltaR(Tau3MuLV) > cut.at(nElectrons_dR)   ) Electrons_OppositeHemisphere.push_back(ie);
   
    Selection_Cut_elect_Pt.at(t).Fill(Ntp->Electron_P4(ie).Pt());
    Selection_Cut_elect_Eta.at(t).Fill(fabs(Ntp->Electron_P4(ie).Eta()));
    Selection_Cut_elect_DeltaR_3mu.at(t).Fill(Ntp->Electron_P4(ie).DeltaR(Tau3MuLV));
   
    }

  pass.at(nElectrons_PF_cut)  = ( Electrons_OppositeHemisphere_PF_cut.size() > 0 );
  pass.at(nElectrons_pT)      = ( Electrons_OppositeHemisphere_pT.size() > 0 );
  pass.at(nElectrons_eta)     = ( Electrons_OppositeHemisphere_eta.size() > 0 );
  pass.at(nElectrons_dR)      = ( Electrons_OppositeHemisphere_dR.size() > 0 );

  if(Electrons_OppositeHemisphere.size()>0)
    {
      unsigned int electron_idx = Electrons_OppositeHemisphere.at(0);
      //value.at(ElectronIsolation)=  Ntp->Electron_puppiPhotonIso(electron_idx)  + Ntp->Electron_trackIso(electron_idx)  + Ntp->Electron_puppiNeutralHadronIso(electron_idx) + Ntp->Electron_puppiChargedHadronIso(electron_idx);
      value.at(ElectronIsolation)=  Ntp->Electron_relativeIsolation(electron_idx);
    }

  pass.at(ElectronIsolation) = true;

  value.at(OSCharge)        =0;
  value.at(Tau3MuIsolation) = -1;
  value.at(VisMass)         = -1;
  value.at(TriggerMatch)    = 0;

  
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
    bool triggerCheck = 0.1;
    if(pass.at(HLT_TriggerOk))
      {
	vector<TLorentzVector> trigobjTriplet;
	for (int i=0; i<Ntp->NTriggerObjects(); i++)
	  {
	    TString name = Ntp->TriggerObject_name(i);
	    //        if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
	    TLorentzVector tmp;
	    tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
	    trigobjTriplet.push_back(tmp);
	  }
	std::vector<TLorentzVector> muonTriplet;
	muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)));
	muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)));
	muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)));

	if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
      }
    value.at(TriggerMatch) = triggerCheck;
    if(Electrons_OppositeHemisphere.size()>0)
      {
      for(auto i : Electrons_OppositeHemisphere)
	{

	  int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
	    Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
	    Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

	  if(Ntp->Electron_charge(i)*Tau3MuCharge == -1) Electrons_OppositeHemisphere_OppositeCharge.push_back(i);
	}
      value.at(OSCharge) = Electrons_OppositeHemisphere_OppositeCharge.size();
      }
    }
  
    pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));
    pass.at(OSCharge) = (value.at(OSCharge) >= cut.at(OSCharge));
    
    if(pass.at(OSCharge))
      {
    	unsigned int electron_index = Electrons_OppositeHemisphere_OppositeCharge.at(0);
    	value.at(VisMass) = (Tau3MuLV + Ntp->Electron_P4(electron_index)).M();
      Selection_Cut_Vis_InvM.at(t).Fill(value.at(VisMass));
      }

    pass.at(Tau3MuIsolation) = (value.at(Tau3MuIsolation) > cut.at(Tau3MuIsolation));
    pass.at(VisMass)         = (value.at(VisMass) > 35  && value.at(VisMass) < 90);
    



    double wobs=1;
    double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 

    //    std::cout<<"   how many electrons i have  " << Electrons_OppositeHemisphere_OppositeCharge .size() << std::endl;
    unsigned int electron_idx = Electrons_OppositeHemisphere_OppositeCharge.at(0);


    unsigned int muon_1_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int muon_2_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int muon_3_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);


    ////////////////////////   sort muons by charge and dR and fill pair masses :
    /////
    vector<unsigned int> idx_vec;
    idx_vec.push_back(muon_1_idx);
    idx_vec.push_back(muon_2_idx);
    idx_vec.push_back(muon_3_idx);

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



    TLorentzVector Muon1LV = Ntp->Muon_P4(muon_1_idx);
    TLorentzVector Muon2LV = Ntp->Muon_P4(muon_2_idx);
    TLorentzVector Muon3LV = Ntp->Muon_P4(muon_3_idx);



    TLorentzVector Tau3muLV = Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) +
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) +
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);



    TLorentzVector ElectronLV = Ntp->Electron_P4(electron_idx);


    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(ElectronLV.Theta()));
    TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());



    float RelativeIsolationMu1 = Ntp->Muon_RelIso(muon_1_idx);
    float RelativeIsolationMu2 = Ntp->Muon_RelIso(muon_2_idx);
    float RelativeIsolationMu3 = Ntp->Muon_RelIso(muon_3_idx);




    Tau3MuRelativeIsolation.at(t).Fill(    Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt()),1);
    ElectronSumIsolation.at(t).Fill( Ntp->Electron_puppiPhotonIso(electron_idx)  + Ntp->Electron_trackIso(electron_idx)  + Ntp->Electron_puppiNeutralHadronIso(electron_idx)   ,1);
    //    ElectronSumIsolation.at(t).Fill( Ntp->Electron_puppiPhotonIso(electron_idx)  + Ntp->Electron_puppiChargedHadronIso(electron_idx)  + Ntp->Electron_puppiNeutralHadronIso(electron_idx)   ,1);

    VisibleDiTauMass.at(t).Fill((ElectronLV + Tau3muLV).M(), 1);
    MTT.at(t).Fill( (Tau3muLV + ElectronLV  + Neutrino_LV).M(), 1);


    bool PlotMCOnly(false);  // and blind for data
    if(id!=1) PlotMCOnly = true;
    if(id==1 && (  (TauRefitLV.M() > 1.1 && TauRefitLV.M() < 1.7233333) or (TauRefitLV.M() > 1.8333333 && TauRefitLV.M() < 2.2)) ) PlotMCOnly=true;


    if(PlotMCOnly)  TripletMass.at(t).Fill(TauRefitLV.M(),1);



    ////////////////////////////////////////
    ///
    TLorentzVector OppositeSideLV = ElectronLV;
    TLorentzVector Neutrino_LV_Guess_Result;//Guessed Neutrino LV
    if(id != WhetherSignalMC)
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
        int MC_Electron_idx = Ntp->getMatchTruthIndex(ElectronLV);//Getting the index of electron from Spectator Tau in MC
        int Tau_MC_idx = Ntp->MCParticle_midx(MC_Electron_idx);//Getting the index of Spectator Tau in MC
        
        TLorentzVector ElectronMC=Ntp->MCParticle_p4(MC_Electron_idx);
        TLorentzVector TauMC=Ntp->MCParticle_p4(Tau_MC_idx);
        TLorentzVector NeutrinoMC=TauMC-ElectronMC;
        bool FoundMCTau(false);
        if(abs(Ntp->MCParticle_pdgid(Tau_MC_idx))==15){
                FoundMCTau = true;
        }
        
        if(FoundMCTau){// Neutrinos satisfying tau invariant mass can be created with any direction of neutrino, PROVING THERE ARE MULTIPLE SOLUTIONS (without taking Z into account). We can construct a neutrino LV consistent with tau and Z mass. The same thing will happen with data, so we can't use this method to distinguish between data and MC? (Or, is it?)
                
                //std::cout << "Found e neutrinos: " << std::endl;
                double Max_Diff(999.0);
                int Division_No_Thet(1);//200
                int Division_No_Phi(1);//1000
                
                double Phi_Visible_Tau = ElectronMC.Phi();
                double Theta_Visible_Tau = ElectronMC.Theta();
                TVector3 Visible_Tau_Vec = ElectronMC.Vect();
                
                for(int Theta_1_idx=0; Theta_1_idx<Division_No_Thet;Theta_1_idx++){
                        for(int Phi_1_idx=0; Phi_1_idx<Division_No_Phi;Phi_1_idx++){
                                
                                double Theta_1 = 0.0 + Theta_1_idx*(0.5)/Division_No_Thet;
                                double Phi_1 = 0.0 + Phi_1_idx*(2*TMath::Pi())/Division_No_Phi;
                                
                                TVector3 Neutrino_Vect_Guess(TMath::Sin(Theta_1)*TMath::Cos(Phi_1),TMath::Sin(Theta_1)*TMath::Sin(Phi_1),TMath::Cos(Theta_1));
                                if(Phi_Visible_Tau >= TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau-2*TMath::Pi();
                                if(Phi_Visible_Tau < -TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau+2*TMath::Pi();
                                Neutrino_Vect_Guess.RotateY(Theta_Visible_Tau);
                                Neutrino_Vect_Guess.RotateZ(Phi_Visible_Tau);
                                
                                double scaling_x = (1.77686*1.77686 - ElectronMC.M()*ElectronMC.M())/(2*(ElectronMC.E()-ElectronMC.Vect().Mag()*TMath::Cos(Neutrino_Vect_Guess.Angle(Visible_Tau_Vec))));//Imposing tau mass criteria
                                TLorentzVector Neutrino_LV_Guess(scaling_x*Neutrino_Vect_Guess,scaling_x);
                                
                                // Boosting to COM frame
                                //TLorentzVector ElectronMC_Boosted_Guess = ElectronMC;
                                //ElectronMC_Boosted_Guess.Boost(-1*(ElectronMC+Neutrino_LV_Guess).BoostVector());
                                //TLorentzVector Neutrino_LV_Boosted_Guess = Neutrino_LV_Guess;
                                //Neutrino_LV_Boosted_Guess.Boost(-1*(ElectronMC+Neutrino_LV_Guess).BoostVector());
                                
                                double Diff_Val = abs((ElectronMC+Neutrino_LV_Guess+Tau3muLV).M() - 91.1876);// Imposing Z mass
                                if(Diff_Val < Max_Diff){
                                        Max_Diff = Diff_Val;
                                        Neutrino_LV_Guess_Result = Neutrino_LV_Guess;
                                }
        
                        }
                }
                
                //std::cout << "Angle from guess result to actual: " << Neutrino_LV_Guess_Result.Vect().Angle(NeutrinoMC.Vect()) << " with Max_Diff: " << Max_Diff << std::endl;
                
                dR_betweenTruth_NeutrinoGuess.at(t).Fill(Neutrino_LV_Guess_Result.DeltaR(NeutrinoMC),1);
                dR_betweenTruth_Tau.at(t).Fill((Neutrino_LV_Guess_Result+ElectronLV).DeltaR(TauMC),1);
                Z_Pt.at(t).Fill((TauMC+Tau3muLV).Pt(),1);
                
                TVector3 MC_MET = (TauMC-NeutrinoMC+Tau3muLV).Vect();
                MC_MET.SetZ(0.);
                MC_MET = -1 * MC_MET;
                
                //std::cout << "MC Delta Phi from actual to MET: " << MC_MET.Phi()-NeutrinoMC.Vect().Phi() << std::endl;
                //std::cout << "Z angle: " << (TauMC+Tau3muLV).Vect().Theta() << " Z (Pt): "<< (TauMC+Tau3muLV).Pt() << " Tau vis to Tau3mu angle: " << (TauMC-NeutrinoMC).Vect().Phi() - Tau3muLV.Vect().Phi() << " Tau to Tau3mu angle: " << (TauMC).Vect().Phi() - Tau3muLV.Vect().Phi() << std::endl;
                
        }//FoundMCNeutrino

      }


  }
}


void  ZTau3MuTaue_PreFC::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





