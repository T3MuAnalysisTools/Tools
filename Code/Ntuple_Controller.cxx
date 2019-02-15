//Ntuple_Controller.cxx IMPLEMENTATION FILE
 

#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"
#include "TF1.h"
#include "Parameters.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include <tuple>

// External code
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
 
///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,cannotObtainHiggsMass(false)
  ,ObjEvent(-1)
  ,isInit(false)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("HTauTauTree");
  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
  for(unsigned int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;
  if(chain==0){
	Logger(Logger::Error) << "chain points to NULL" << std::endl;
  }
  Logger(Logger::Info) << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  Logger(Logger::Info) << "Ntuple Configured" << std::endl;

  gRandom->SetSeed(1234);
  tauCorrection = "";
  muonCorrection = "";
  elecCorrection = "";
  jetCorrection = "";
}

///////////////////////////////////////////////////////////////////////
//
// Function: void InitEvent()
//
// Purpose: Initialize variables etc on event base
//
///////////////////////////////////////////////////////////////////////

void Ntuple_Controller::InitEvent(){
	Muon_corrected_p4.clear();
	//	Muon_corrected_p4.resize(NMuons());
	Muon_isCorrected = false;

	// after everything is initialized
	isInit = true;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
  isInit = false;
  InitEvent();
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  Logger(Logger::Verbose) << "Cleaning up" << std::endl;
  delete Ntp;
  //  delete rmcor;
  Logger(Logger::Verbose) << "Complete." << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
	Logger(Logger::Info) << "Starting D3PD cloning" << std::endl;
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  Logger(Logger::Info) << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  Logger(Logger::Warning) << "ThinTree not implemented." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
Long64_t  Ntuple_Controller::GetMCID(){
	Long64_t  DataMCTypeFromTupel = Ntp->DataMC_Type_idx;
	if (DataMCTypeFromTupel == DataMCType::DY_ll_Signal && HistoC.hasID(DataMCType::DY_ll_Signal)) {
		for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
			if (abs(MCSignalParticle_pdgid(i)) == PDGInfo::Z0) {
				if (fabs(MCSignalParticle_p4(i).M() - PDG_Var::Z_mass()) < 3 * PDG_Var::Z_width()) {
					return DataMCType::Signal;
				}
			}
		}
		return DataMCTypeFromTupel;
	}

	int dmcType = -999;

	 if (HistoC.hasID(DataMCTypeFromTupel % 100)) {
	dmcType = DataMCTypeFromTupel % 100;
	 }
	 return dmcType;
}

// return DataMCType without mass information
int Ntuple_Controller::GetStrippedMCID(){
	return GetMCID() % 100;
}

// return path of input dataset (as given in the line InputNtuples in Input.txt)
TString Ntuple_Controller::GetInputNtuplePath(){
	Parameters Par; // assumes configured in Analysis.cxx
	TString dsPath;
	Par.GetString("InputNtuples:",dsPath);
	return dsPath;
}

// return name of input dataset
TString Ntuple_Controller::GetInputDatasetName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName( gSystem->DirName(dsPath) ) );
}

//
TString Ntuple_Controller::GetInputPublishDataName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName(dsPath) );
}

// determine Higgs mass (from Dataset name or fallback options)
int Ntuple_Controller::getSampleHiggsMass(){
	int mass = -999;

	// default method: analyze dataset name
	mass = readHiggsMassFromString( GetInputNtuplePath() );
	if (mass >= 0) return mass;

	// first fallback: analyze filename (only working when running on GRID)
	mass = readHiggsMassFromString( Get_File_Name() );
	if (mass >= 0) return mass;

	// second fallback: get Higgs mass from MC info
	Logger(Logger::Warning) << "Not able to obtain Higgs mass neither from dataset nor from file name."
	<< "\n\tMake sure the line InputNtuples is set correctly in your Input.txt. "
	<< "\n\tFor now, we will fall back to obtaining the Higgs mass from the generator information."
	<< "\n\tBe aware that SOME EVENTS WILL END UP IN THE WRONG HISTOGRAM!!!" << std::endl;
	return getHiggsSampleMassFromGenInfo();
}

int Ntuple_Controller::readHiggsMassFromString(TString input){
	// loop over possible masses
	for (int m = 100; m < 200; m = m+5){
		if ( input.Contains("M-" + TString::Itoa(m, 10) + "_") ) return m;
	}

	// mass not found
	return -999;
}

double Ntuple_Controller::getResonanceMassFromGenInfo(bool useZ0 /* = true*/, bool useHiggs0 /* = true*/, bool useW /* = true*/){
	for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
		unsigned pdgid = abs(MCSignalParticle_pdgid(i));
		bool Z0 = useZ0 && (pdgid == PDGInfo::Z0);
		bool H0 = useHiggs0 && (pdgid == PDGInfo::Higgs0);
		bool W = useW && (pdgid == PDGInfo::W_plus);
		if (Z0 || H0 || W) {
			return MCSignalParticle_p4(i).M();
		}
	}
	return -999;
}

int Ntuple_Controller::getHiggsSampleMassFromGenInfo(){
	double resMass = getResonanceMassFromGenInfo(false, true, false);
	for (int m = 100; m < 200; m = m+5){
		if (fabs(resMass - m) < 2.5 ) {
			return m;
		}
	}
	return -999;
}



 double Ntuple_Controller::DeltaPhi(double angle1,double angle2)
  {
    double diff=angle1-angle2;
    while(diff>=TMath::Pi())diff=diff-TMath::TwoPi();
    while(diff<=-TMath::Pi())diff=diff+TMath::TwoPi();
    return diff;
  }

float Ntuple_Controller::DeltaRDau(int dau1idx, int dau2idx)
{
  TLorentzVector v1, v2;
  v1 =Daughters_P4(dau1idx);
  v2 =Daughters_P4(dau2idx);
  return v1.DeltaR(v2);
}


bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &tau_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==(int)parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	if((unsigned int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
	  return false;
	}
      }
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	unsigned int tauidx=MCSignalParticle_Tauidx(i).at(j);
	Logger(Logger::Verbose) << "MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << " " << Ntp->MCTau_JAK->size() << " " << Ntp->MCTauandProd_pdgid->size() << std::endl;
	if((int)MCTau_JAK(tauidx)==tau_jak){ tau_idx=tauidx;Boson_idx=i;return true;}
      }
    }
  }
  return false;
}

bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
        if((unsigned int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
          Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
          return false;
        }
      }
      if(MCSignalParticle_Tauidx(i).size()==2){
	tau1_idx=MCSignalParticle_Tauidx(i).at(0);
        tau2_idx=MCSignalParticle_Tauidx(i).at(1);
	Boson_idx=i;
	return true;
      }
    }
  }
  return false;
}


bool Ntuple_Controller::GetTriggerIndex(TString n,  int &i){
  for(i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      if(name.Contains(n))return true;
    } 
	return false;
 }


bool  Ntuple_Controller::CheckIfAnyPassed(  std::vector<int> list){
  for(unsigned int itrig = 0; itrig < list.size(); itrig++){
    if(TriggerAccept(list.at(itrig)))      return true;
    
  }
  return false;
}


std::vector<int> Ntuple_Controller::GetVectorTriggers(TString n){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
	TString name=TriggerName(i);
	if(name.Contains(n)) out.push_back(i) ;
      } 
	  return out;
 }

std::vector<int> Ntuple_Controller::GetVectorTriggers(std::vector<TString>  v){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      for(unsigned int j=0; j<v.size(); j++){
	if(name.Contains(v.at(j))) out.push_back(i) ;
      }
    } 
    return out;
}

std::vector<int> Ntuple_Controller::GetVectorTriggersFullMatch(std::vector<TString>  v){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      bool cpattern(true);
      for(unsigned int j=0; j<v.size(); j++){
	if(!name.Contains(v.at(j))) cpattern =false;
      }
      if(cpattern)  out.push_back(i) ;
    } 
    return out;
}

std::vector<int> Ntuple_Controller::GetVectorCrossTriggers(TString n1,TString n2,TString f1,TString f2){
    std::vector<int> out;



    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      if(name.Contains(n1) &&  name.Contains(n2)  && (!name.Contains(f1)  && !name.Contains(f2))   ) out.push_back(i) ;
    } 
    return out;
}



// // PFTau significance, using the reffited primary and secondary vertices
// double Ntuple_Controller::PFTau_FlightLength_significance(unsigned int i) {
// 	TVector3 PVpos = PFTau_TIP_primaryVertex_pos(i);
// 	TMatrixTSym<double> PVcov = PFTau_TIP_primaryVertex_cov(i);
// 	TVector3 SVpos = PFTau_TIP_secondaryVertex_pos(i);
// 	TMatrixTSym<double> SVcov = PFTau_TIP_secondaryVertex_cov(i);

// 	return PFTau_FlightLength_significance(PVpos, PVcov, SVpos, SVcov);
// }

// // calculate flight length significance from primary and secondary vertex info
 double Ntuple_Controller::PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
   TVector3 SVPV = sv - pv;
   TVectorF FD;
   FD.ResizeTo(3);
   FD(0) = SVPV.X();
   FD(1) = SVPV.Y();
   FD(2) = SVPV.Z();

   TMatrixT<double> PVcv;
   PVcv.ResizeTo(3,3);
   for(int nr =0; nr<PVcov.GetNrows(); nr++){
     for(int nc =0; nc<PVcov.GetNcols(); nc++){
       PVcv(nr,nc) = PVcov(nr,nc);
     }
   }
   TMatrixT<double> SVcv;
   SVcv.ResizeTo(3,3);
   for(int nr =0; nr<SVcov.GetNrows(); nr++){
     for(int nc =0; nc<SVcov.GetNcols(); nc++){
       SVcv(nr,nc) = SVcov(nr,nc);
     }
   }

   TMatrixT<double> SVPVMatrix(3,1);
   for(int i=0; i<SVPVMatrix.GetNrows();i++){
     SVPVMatrix(i,0)=FD(i);
   }

   TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT.T();

   TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   double sigmaabs = sqrt(lambda2(0,0))/SVPV.Mag();
   double sign = SVPV.Mag()/sigmaabs;

   return sign;
 }

//// Generator Information
int Ntuple_Controller::matchTruth(TLorentzVector tvector){
	double testdr=0.3;
	int pdgid = 0;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(tvector.DeltaR(MCParticle_p4(i))<testdr){
				testdr = tvector.DeltaR(MCParticle_p4(i));
				pdgid = MCParticle_pdgid(i);
			}
		}
	}
	return pdgid;
}
bool Ntuple_Controller::matchTruth(TLorentzVector tvector, int pid, double dr){
	if (getMatchTruthIndex(tvector, pid, dr) >= 0) return true;
	return false;
}
int Ntuple_Controller::getMatchTruthIndex(TLorentzVector tvector, int pid, double dr){
	int index = -9;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(fabs(MCParticle_pdgid(i))==pid){
				if(tvector.DeltaR(MCParticle_p4(i))<dr) index = i;
			}
		}
	}
	return index;
}

int Ntuple_Controller::MCTau_getDaughterOfType(unsigned int i_mcTau, int daughter_pdgid, bool ignoreCharge /*= true*/){
	int matchedIndex = -1;
	for(int i_dau=0; i_dau < NMCTauDecayProducts(i_mcTau); i_dau++){
		if( ignoreCharge ){
			if( abs(MCTauandProd_pdgid(i_mcTau, i_dau)) == abs(daughter_pdgid) )
				matchedIndex = i_dau;
		}
		else{
			if( MCTauandProd_pdgid(i_mcTau, i_dau) == daughter_pdgid )
				matchedIndex = i_dau;
		}
	}
	return matchedIndex;
}

TLorentzVector Ntuple_Controller::BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2){
	TVector3 boostvector = TLV1.BoostVector();
	TLorentzVector boosted_TLV2(TLV2);
	boosted_TLV2.Boost(-boostvector);
	return boosted_TLV2;
}

// get visible/invisible part of gen Tau 4-vector
TLorentzVector Ntuple_Controller::MCTau_invisiblePart(unsigned int i){
	TLorentzVector lv(0,0,0,0);
	for(int i_dau=1; i_dau < NMCTauDecayProducts(i); i_dau++){
		if( abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_e ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_mu ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_tau)
				lv += MCTauandProd_p4(i, i_dau);
	}
	if (lv.Pt() == 0.0)
		Logger(Logger::Warning) << "This tau decay has no neutrinos!?" << std::endl;

	return lv;
}
TLorentzVector Ntuple_Controller::MCTau_visiblePart(unsigned int i){
	TLorentzVector lv = MCTau_p4(i);
	lv -= MCTau_invisiblePart(i);
	return lv;
}


//// draw decay chain
void Ntuple_Controller::printMCDecayChainOfEvent(bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "=== Draw MC decay chain of event ===" << std::endl;
	for(unsigned int par = 0; par < NMCParticles(); par++){
		if ( !MCParticle_hasMother(par) )
			printMCDecayChainOfMother(par, printStatus, printPt, printEtaPhi, printQCD);
	}
}
void Ntuple_Controller::printMCDecayChainOfMother(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "Draw decay chain for mother particle at index " << par << " :" << std::endl;
	printMCDecayChain(par, 0, printStatus, printPt, printEtaPhi, printQCD);
}
void Ntuple_Controller::printMCDecayChain(unsigned int par, unsigned int level, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	std::ostringstream out;
	for(unsigned int l = 0; l < level; l++) out << "    ";
	out << MCParticleToString(par, printStatus, printPt, printEtaPhi);
	if (!printQCD && MCParticle_pdgid(par) == 92) {out << "    QCD stuff truncated ...\n"; std::cout << out.str(); return;}
	out << "\n";
	std::cout << out.str();
	for (unsigned int i_dau = 0; i_dau < MCParticle_childidx(par).size(); i_dau++){
		if (MCParticle_childidx(par).at(i_dau) != (int)par) // skip cases where particles are daughters of themselves
			printMCDecayChain(MCParticle_childidx(par).at(i_dau), level+1, printStatus, printPt, printEtaPhi, printQCD);
	}
}

std::string Ntuple_Controller::MCParticleToString(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi){
	std::ostringstream out;
	out << "+-> ";
	out << PDGInfo::pdgIdToName( MCParticle_pdgid(par) );
	if (printStatus) out << " (status = " << MCParticle_status(par) << ", idx = " << par << ")";
	if (printPt || printEtaPhi) out <<  " [";
	if (printPt) out << "pT = " << MCParticle_p4(par).Pt() << "GeV";
	if (printEtaPhi) out << " eta = " << MCParticle_p4(par).Eta() << " phi = " << MCParticle_p4(par).Phi();
	if (printPt || printEtaPhi) out << "]";
	return out.str();
}
bool Ntuple_Controller::CheckDecayID(unsigned  int jak1, unsigned int jak2){
  if(isData()) return false;
  // jak  = 2 - muon
  // jak  = 3 - pion
  // jak  = 4 - rho
  // jak  = 5 - a1
  // if(id!=998) return false;
  bool decayid= false;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
    // std::cout<<"NC: Ntp->MCSignalParticle_p4->size()  "<<Ntp->MCSignalParticle_p4->size()<<std::endl;
    // std::cout<<"NC: Ntp->MCSignalParticle_Tauidx->at(iz).size()  "<<Ntp->MCSignalParticle_Tauidx->at(iz).size()<<std::endl;
    // std::cout<<"  MCSignalParticle_pdgid->at(iz)  "<< Ntp-> MCSignalParticle_pdgid->at(iz) <<std::endl;
    if(fabs(Ntp-> MCSignalParticle_pdgid->at(iz) )!=23) return false;
    if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      //     std::cout<<" Ntp->MCTau_JAK-> "<< Ntp->MCTau_JAK->size()<<std::endl;
      if(Ntp->MCTau_JAK->at(0) == jak1 and Ntp->MCTau_JAK->at(1) ==jak2 ){ decayid = true;}
      else if(Ntp->MCTau_JAK->at(0) ==jak2  and Ntp->MCTau_JAK->at(1) ==jak1){decayid  = true;}
    }
  }
  return decayid;
}
TLorentzVector Ntuple_Controller::GetTruthTauLV(unsigned int jak,unsigned int number){
  //number = 0 or 1 respectively for the first tau or for the second tau;
  TLorentzVector tau(0,0,0,0);
  bool DecayOK = false;
  unsigned int tauIndex;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
    if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      if(Ntp->MCTau_JAK->at(0) == jak){tauIndex=0; DecayOK = true;}
      else if(  Ntp->MCTau_JAK->at(1) ==jak ){ tauIndex=1; DecayOK = true;}
      if(DecayOK){
        tau = TLorentzVector(Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(1),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(2),
                             Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(3),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(0));
      }

      if((Ntp->MCTau_JAK->at(0)==Ntp->MCTau_JAK->at(1)) && (Ntp->MCTau_JAK->at(0)==jak))
	{
	  tau = TLorentzVector(Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(number)).at(0).at(1),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(number)).at(0).at(2),
                             Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(number)).at(0).at(3),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(number)).at(0).at(0));
	}
    }
  }
  return tau;
}


TLorentzVector Ntuple_Controller::GetTruthTauProductLV(unsigned int jak, int pdgID,unsigned int number){
  TLorentzVector tauProd(0,0,0,0);
  bool DecayOK = false;
  unsigned int tauIndex;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
     if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      if(Ntp->MCTau_JAK->at(0) == jak){tauIndex=0; DecayOK = true;}
      else if(  Ntp->MCTau_JAK->at(1) ==jak ){ tauIndex=1; DecayOK = true;}
      
      unsigned int NDec;
      if(DecayOK){
        if(0<=tauIndex && tauIndex<NMCTaus()){ NDec = Ntp->MCTauandProd_p4->at(tauIndex).size();}
        else NDec= 0;
        for(unsigned int iProd =0; iProd < NDec; iProd++ ){//cout<<Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)<<endl;
          if(abs( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd))==pdgID){
            //if( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)==pdgID){
              tauProd = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(2),
                                       Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(0));
              // }
          }
        }
      }
      if((Ntp->MCTau_JAK->at(0)==Ntp->MCTau_JAK->at(1)) && (Ntp->MCTau_JAK->at(0)==jak))
	{
	  if(0<=number && number<NMCTaus()){ NDec = Ntp->MCTauandProd_p4->at(number).size();}
	  else NDec= 0;
	  for(unsigned int iProd =0; iProd < NDec; iProd++ ){
	    if(abs( Ntp->MCTauandProd_pdgid->at(number).at(iProd))==pdgID){
	      //if( Ntp->MCTauandProd_pdgid->at(number).at(iProd)==pdgID){
              tauProd = TLorentzVector(Ntp->MCTauandProd_p4->at(number).at(iProd).at(1),Ntp->MCTauandProd_p4->at(number).at(iProd).at(2),
                                       Ntp->MCTauandProd_p4->at(number).at(iProd).at(3),Ntp->MCTauandProd_p4->at(number).at(iProd).at(0));
              // }
	    }
	  }
	}
    }
  }
    return tauProd;
}


bool Ntuple_Controller::isPVCovAvailable(){ // sometimes returns zero size matrix (rare)
  if(Ntp->pv_cov->size()!=6)  return false; 
  return true;

}
TMatrixTSym<float> Ntuple_Controller::PFTau_TIP_primaryVertex_cov(){
  TMatrixTSym<float> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
      V_cov(j,k)=Ntp->pv_cov->at(l);
      V_cov(k,j)=Ntp->pv_cov->at(l);
      l++;
    }
  }
  //  std::cout<<"  pv is good"<< std::endl; V_cov.Print();
  return  V_cov;
}

TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_secondaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      V_cov(j,k)=Ntp->PFTauSVCov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTauSVCov->at(i).at(l);
      l++;
    }
  }
  //    std::cout<<"  sv is good"<< std::endl; V_cov.Print();
  return  V_cov;
}

 LorentzVectorParticle Ntuple_Controller::PFTau_a1_lvp(unsigned int i){
   TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
   TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
   int l=0;
   if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
     for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
       a1_par(k,0)=Ntp->PFTau_a1_lvp->at(i).at(k);
       for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
 	a1_cov(k,j)=Ntp->PFTau_a1_cov->at(i).at(l);
 	a1_cov(j,k)=Ntp->PFTau_a1_cov->at(i).at(l);
 	l++;
       } 
     }
   }
   return LorentzVectorParticle(a1_par,a1_cov,Ntp->PFTau_a1_pdgid->at(i),Ntp->PFTau_a1_charge->at(i),Ntp->PFTau_a1_B->at(i));
 }


   
double Ntuple_Controller::dxySigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return (-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt();
}
double Ntuple_Controller::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dxySigned(fourvector, poca, vtx));
}


double Ntuple_Controller::dzSigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2);
}
double Ntuple_Controller::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dzSigned(fourvector, poca, vtx));
}


// function to sort any objects by any value in descending order
std::vector<int> Ntuple_Controller::sortObjects(std::vector<int> indices, std::vector<double> values){
	if (indices.size() != values.size()){
		Logger(Logger::Warning) << "Please make sure indices and values have same size for sorting. Abort." << std::endl;
		return std::vector<int>();
	}
	// create vector of pairs to allow for sorting by value
	std::vector< std::pair<int, double> > pairs;
	for(unsigned int i = 0; i<values.size(); i++ ){
		pairs.push_back( std::make_pair(indices.at(i),values.at(i)) );
	}
	// sort vector of pairs
	std::sort(pairs.begin(), pairs.end(), sortIdxByValue());
	// create vector of indices in correct order
	std::vector<int> sortedIndices;
	for(unsigned int i = 0; i<pairs.size(); i++){
		sortedIndices.push_back(pairs.at(i).first);
	}
	return sortedIndices;
}



