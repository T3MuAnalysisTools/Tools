//Ntuple_Controller.cxx IMPLEMENTATION FILE
#include "Ntuple_Controller.h"
#include "Parameters.h"
#include <tuple>
#include "Logger.h"
#include "PDGInfo.h"
using namespace std;

// External code

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,ObjEvent(-1)
  ,isInit(false)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("T3MTree/t3mtree");
  //  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
  for(unsigned int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;
  if(chain==0){
    //	Logger(Logger::Error) << "chain points to NULL" << std::endl;
  }
  std::cout << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  std::cout << "Ntuple Configured" << std::endl;

  gRandom->SetSeed(1234);
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
  std::cout << "Cleaning up" << std::endl;
  delete Ntp;
  //  delete rmcor;
  std::cout << "Complete." << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
    std::cout << "Starting D3PD cloning" << std::endl;
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
  std::cout << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  std::cout << "ThinTree not implemented..." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
  }
}



//Physics get Functions
Long64_t  Ntuple_Controller::GetMCID(){

  Long64_t  DataMCTypeFromTupel =  Ntp->Event_DataMC_Type;
  return DataMCTypeFromTupel;

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

double Ntuple_Controller::DeltaPhi(double angle1,double angle2)
{
  double diff=angle1-angle2;
  while(diff>=TMath::Pi())diff=diff-TMath::TwoPi();
  while(diff<=-TMath::Pi())diff=diff+TMath::TwoPi();
  return diff;
}

//-------------------------- print vector -------------------
template<typename T>
void printVec(int size, T& vec){
    for (int i = 0; i<size; i++) cout<<vec[i]<<" ";
    cout<<endl;
}
//-------------------------- --------------------------------- 

//-------------------------- return deltaR -------------------

double Ntuple_Controller::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}

//-------------------------- --------------------------------- 


//-------------------------- return GEN pdgId -------------------
// Returns the pdgId of the GEN particle matched to a RECO particle

int Ntuple_Controller::GENMatchedPdgId(TLorentzVector vec){
	float dR_min = 999.0;
	int pdgId = 9999999;
	for (unsigned int it=0; it<NMCSignalParticles(); ++it){
		float tmp_dR = MCSignalParticle_p4(it).DeltaR(vec);
		if (tmp_dR<dR_min && tmp_dR<0.05){
			dR_min = tmp_dR;
			pdgId = MCSignalParticle_pdgid(it);
		}
	} 
	return pdgId;
}

//-------------------------- --------------------------------- -------------------

//-------------------------- return GEN LV -------------------
// Returns the lorentz vector of the GEN particle matched to a RECO particle

TLorentzVector Ntuple_Controller::GENMatchedLV(TLorentzVector vec){
	float dR_min = 999.0;
	TLorentzVector tmp_vec(0,0,0,0);
	for (unsigned int it=0; it<NMCSignalParticles(); ++it){
		float tmp_dR = MCSignalParticle_p4(it).DeltaR(vec);
		if (tmp_dR<dR_min && tmp_dR<0.05){
			dR_min = tmp_dR;
			tmp_vec = MCSignalParticle_p4(it);
		}
	} 
	return tmp_vec;
}

//-------------------------- --------------------------------- -------------------

//-------------------------- Match GEN Ds to 2mu+trk candidate -------------------
// Returns the dR of best GEN machted Ds
float Ntuple_Controller::DsGenMatch(unsigned int tmp_idx){

	float ds_dR = 999.0;

	if (tmp_idx>=NTwoMuonsTrack()) {
		cout<<"Index out of rangei"<<endl;
		return 999.0;
		}
   
   for (unsigned int ngen=0; ngen<NMCSignalParticles(); ngen++){
      bool dspi_flag = 0, dsphi_flag = 0;
		float tmp_phi_dR = 999.0, tmp_pi_dR = 999.0;

		int track_idx = -1, mu1_idx = -1, mu2_idx = -1;

		TLorentzVector tmp_pi_p4(0,0,0,0);
		TLorentzVector tmp_phi_p4(0,0,0,0);

      if (fabs(MCSignalParticle_pdgid(ngen))==431){
         dspi_flag = 0;
         dsphi_flag = 0;

			tmp_phi_dR = 999.0;
			tmp_pi_dR = 999.0;

         for (int nchild=0; nchild<MCSignalParticle_Nchilds(ngen); nchild++){

            //Find a pi that is a decay product of ds
            if(fabs(MCSignalParticle_childpdgid(ngen,nchild))==211){
               dspi_flag = 1;
               unsigned int tmp_track_idx = TwoMuonsTrackTrackIndex(tmp_idx).at(0);
               float dR = Track_P4(tmp_track_idx).DeltaR(MCSignalParticle_child_p4(ngen,nchild));
               if (dR<0.05 && dR<tmp_pi_dR) { 
						tmp_pi_dR = dR;
						tmp_pi_p4 = MCSignalParticle_child_p4(ngen,nchild);
						track_idx = tmp_track_idx;
                  }
              }

              //Find a phi that is a decay product of Ds
              if (fabs(MCSignalParticle_childpdgid(ngen,nchild))==333){
                dsphi_flag = 1;
                unsigned int tmp_mu1_idx = TwoMuonsTrackMuonIndices(tmp_idx).at(0);
                unsigned int tmp_mu2_idx = TwoMuonsTrackMuonIndices(tmp_idx).at(1);
                float dR = MCSignalParticle_child_p4(ngen,nchild).DeltaR(Muon_P4(tmp_mu1_idx)+Muon_P4(tmp_mu2_idx));
                if (dR<0.05 && dR<tmp_phi_dR) { 
					 	tmp_phi_dR = dR;
						tmp_phi_p4 = MCSignalParticle_child_p4(ngen,nchild);
						mu1_idx = tmp_mu1_idx;
						mu2_idx = tmp_mu2_idx;
                  }
              }

            }

            if (dspi_flag && dsphi_flag){
					ds_dR = (tmp_phi_p4+tmp_pi_p4).DeltaR(Muon_P4(mu1_idx)+Muon_P4(mu2_idx)+Track_P4(track_idx));
            }
        } 
      }
	return ds_dR;
}

//-------------------------- --------------------------------- -------------------

std::vector<unsigned int> Ntuple_Controller::SortedPtMuons(std::vector<unsigned int> indices){
  
  std::vector<unsigned int> out;
  unsigned int i1,i2,i3;

  double pt1 = Muon_P4(indices.at(0)).Pt();
  double pt2 = Muon_P4(indices.at(1)).Pt();
  double pt3 = Muon_P4(indices.at(2)).Pt();

  if(pt1>pt2)
    {
      if(pt2>pt3)
	{
	  i1=indices.at(0); i2 = indices.at(1); i3 = indices.at(2);
	}
      else if(pt1>pt3)
	{
	  i1=indices.at(0); i2 = indices.at(2); i3 = indices.at(1);
	}
      else
	{
	  i1=indices.at(2); i2 = indices.at(0); i3 = indices.at(1);
	}
    }
  else
    {
      if(pt1>pt3)
	{
	  i1=indices.at(1); i2 = indices.at(0); i3 = indices.at(2);
	}
      else if(pt2>pt3)
	{
	  i1=indices.at(1); i2 = indices.at(2); i3 = indices.at(0);
	}
      else
	{
	  i1=indices.at(2); i2 = indices.at(1); i3 = indices.at(0);
	}
    }
  out.push_back(i1);
  out.push_back(i2);
  out.push_back(i3);

  return out;
}



std::vector<unsigned int> Ntuple_Controller::SortedChargeMuons(std::vector<unsigned int> indices){

  std::vector<unsigned int> out;
  
  double q1 = Muon_charge(indices[0]);
  double q2 = Muon_charge(indices[1]);
  double q3 = Muon_charge(indices[2]);

  double pt1 = Muon_P4(indices[0]).Pt();
  double pt2 = Muon_P4(indices[1]).Pt();
  double pt3 = Muon_P4(indices[2]).Pt();

  unsigned int i1,i2,i3;

  if (q1==q2){
    i1 = indices[2];
    i2 = pt1>pt2?indices[0]:indices[1];
    i3 = pt1>pt2?indices[1]:indices[0];
  }
  else if (q2==q3){
    i1 = indices[1];
    i2 = pt2>pt3?indices[1]:indices[2];
    i3 = pt2>pt3?indices[2]:indices[1];
  }
  else if (q1==q3){
    i1 = indices[2];
    i2 = pt1>pt3?indices[0]:indices[2];
    i3 = pt1>pt3?indices[2]:indices[0];
  }

  if (abs(q1+q2+q3)>1.1) Logger(Logger::Warning)<< "Sum of charges is greater than 1"<<endl;

  out.push_back(i1);
  out.push_back(i2);
  out.push_back(i3);

  return out;
}



TLorentzVector Ntuple_Controller::matchToTruthTauDecay(TLorentzVector vector){
  TLorentzVector out(0,0,0,0);
  double dr(0.3);
  if(NMCTaus()==0)
    {
      Logger(Logger::Warning) << "No truth tau leptons in this event found; return TLorentzVector(0,0,0,0) " << std::endl; return out;
    }
  for(int i=0; i < NMCTaus(); i++ ){
    for(int j =0; j < NMCTauDecayProducts(i); j++){
      if(MCTauandProd_p4(i,j).DeltaR(vector) < dr){
	out=MCTauandProd_p4(i,j);
	dr = MCTauandProd_p4(i,j).DeltaR(vector);
      }
    }
  }
  return out;
}



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

int Ntuple_Controller::getMatchTruthIndex(TLorentzVector tvector){
  int index = -9;
  double dr=1.;
  for(unsigned int i=0;i<NMCParticles();i++){
    if(MCParticle_p4(i).Pt()>0.){
      if(tvector.DeltaR(MCParticle_p4(i))<dr){
	index = i;
	dr = tvector.DeltaR(MCParticle_p4(i));
      }
    }
  }
  return index;
}



std::vector<int> Ntuple_Controller::MuonStandardSelectorBitMask(unsigned int MuonIndex){

  std::vector<int> out;
  out.push_back(1);
 
  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdMedium))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdMediumPrompt))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdGlobalHighPt))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdTrkHighPt))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoVeryLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoMedium))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoVeryTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::TkIsoLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::TkIsoTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::SoftCutBasedId))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::SoftMvaId))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaMedium))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoLoose))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoMedium))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoTight))out.push_back(1);
  else out.push_back(0);

  if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoVeryTight))out.push_back(1);
  else out.push_back(0);

  return out;
}


//// draw decay chain
void Ntuple_Controller::printMCDecayChainOfEvent(bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
  Logger(Logger::Info) << "=== Draw MC decay chain of event ===" << std::endl;
  for(unsigned int par = 0; par < NMCParticles(); par++){
    if ( !MCParticle_hasMother(par) )
      printMCDecayChainOfMother(par, printStatus, printPt, printEtaPhi, printQCD);
  }
}

void Ntuple_Controller::printMCDecayChainOfParticle(unsigned int index, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
  Logger(Logger::Info) << "=== Draw MC decay chain of particle ===" << std::endl;
      printMCDecayChainOfMother(index, printStatus, printPt, printEtaPhi, printQCD);
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
