#include "SignalSelector.h"
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

SignalSelector::SignalSelector(TString Name_, TString id_):
  Selection(Name_,id_),
  tauMinMass_(1.73),
  tauMaxMass_(1.81),
  tauMinSideBand_(1.6),
  tauMaxSideBand_(2.0),
  phiVetoCut1(0.96),
  phiVetoCut2(1.07),
  rmgCutVeto1(0.77),  // rmg = rho&omega
  rmgCutVeto2(0.812), // rmg = rho&omega
  PEMassResolutionCut1_(0.007),
  PEMassResolutionCut2_(0.01),
  mvaA1_(0.125), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2_(0.228),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1_(0.135),
  mvaB2_(0.254),
  mvaC1_(0.125),
  mvaC2_(0.230),
  mvaA1train11_(0.116), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2train11_(0.248),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1train11_(0.151),
  mvaB2train11_(0.256),
  mvaC1train11_(0.127),
  mvaC2train11_(0.2327),
  mvaA1train8_(0.12), // optimal cuts for trainings weights/August_A(BC)_BDT.weights.xml
  mvaA2train8_(0.25),  // obtained by Code/CommonUtils/tmva/Get_BDT_cut.cxx
  mvaB1train8_(0.14),
  mvaB2train8_(0.24),
  mvaC1train8_(0.163),
  mvaC2train8_(0.280),
  mvaBTrainA1_(0.090),
  mvaBTrainA2_(0.182),
  mvaBTrainB1_(0.121),
  mvaBTrainB2_(0.211),
  mvaBTrainC1_(0.127),
  mvaBTrainC2_(0.211),
  mvaDTrainA1_(0.153),
  mvaDTrainA2_(0.212),
  mvaDTrainB1_(0.154),
  mvaDTrainB2_(0.224),
  mvaDTrainC1_(0.129),
  mvaDTrainC2_(0.210)
{
  // This is a class constructor;
}



SignalSelector::~SignalSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SignalSelector::Configure(){


  //  This mini tree is for limit extraction

  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("xv",&xv);
  T3MMiniTree->Branch("phiv",&phiv);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("bdt",&bdt);
  T3MMiniTree->Branch("category",&category);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  T3MMiniTree->Branch("mDr1",&mDr1);
  T3MMiniTree->Branch("mDr2",&mDr2);
  T3MMiniTree->Branch("LumiScale",&LumiScale);
  T3MMiniTree->Branch("A1",&mvaA1);
  T3MMiniTree->Branch("A2",&mvaA2);
  T3MMiniTree->Branch("B1",&mvaB1);
  T3MMiniTree->Branch("B2",&mvaB2);
  T3MMiniTree->Branch("C1",&mvaC1);
  T3MMiniTree->Branch("C2",&mvaC2);


  TString basedir_July = "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVATreeJul_13_2021/Code/CommonUtils/IterativeTrain/";


  std::vector<string> names_train0_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08"};
  std::vector<float> values_train0_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08};
  reader_train0_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train0_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train0_C = new TMVA::Reader( "!Color:!Silent" );
  for(unsigned int i =0 ; i < names_train0_A.size(); i++)
    {

      reader_train0_A->AddVariable(names_train0_A.at(i),  &values_train0_A.at(i));
      reader_train0_B->AddVariable(names_train0_A.at(i),  &values_train0_A.at(i));
      reader_train0_C->AddVariable(names_train0_A.at(i),  &values_train0_A.at(i));

    }
      reader_train0_A->BookMVA( "BDT", basedir_July + "output_0_A/weights/TMVAClassification_BDT.weights.xml" );
      reader_train0_B->BookMVA( "BDT", basedir_July + "output_0_B/weights/TMVAClassification_BDT.weights.xml" );
      reader_train0_C->BookMVA( "BDT", basedir_July + "output_0_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train1_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca"};
  std::vector<float> values_train1_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca};

  reader_train1_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train1_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train1_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train1_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train1_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train1_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train1_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train1_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train1_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train1_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train1_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train1_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train1_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train1_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);



  reader_train1_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train1_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train1_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train1_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train1_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train1_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train1_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train1_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train1_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train1_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train1_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);


  reader_train1_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train1_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train1_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train1_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train1_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train1_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train1_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train1_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train1_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train1_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train1_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);



  reader_train1_A->BookMVA( "BDT", basedir_July + "output_1_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train1_B->BookMVA( "BDT", basedir_July + "output_1_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train1_C->BookMVA( "BDT", basedir_July + "output_1_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train2_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3"};

  std::vector<float>  values_train2_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3};
  reader_train2_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train2_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train2_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train2_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train2_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train2_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train2_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train2_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train2_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train2_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train2_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train2_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train2_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train2_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train2_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train2_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train2_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);



  reader_train2_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train2_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train2_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train2_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train2_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train2_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train2_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train2_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train2_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train2_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train2_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train2_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train2_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train2_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);



  reader_train2_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train2_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train2_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train2_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train2_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train2_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train2_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train2_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train2_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train2_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train2_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train2_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train2_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train2_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);





  reader_train2_A->BookMVA( "BDT", basedir_July + "output_2_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train2_B->BookMVA( "BDT", basedir_July + "output_2_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train2_C->BookMVA( "BDT", basedir_July + "output_2_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train3_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoPhiKKMass_Mu1","var_IsoPhiKKMass_Mu2","var_IsoPhiKKMass_Mu3"};
  std::vector<float>  values_train3_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoPhiKKMass_Mu1,var_IsoPhiKKMass_Mu2,var_IsoPhiKKMass_Mu3};

  reader_train3_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train3_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train3_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train3_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train3_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train3_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train3_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train3_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train3_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train3_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train3_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train3_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train3_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train3_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train3_A->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train3_A->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train3_A->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);


  reader_train3_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train3_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train3_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train3_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train3_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train3_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train3_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train3_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train3_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train3_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train3_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train3_B->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train3_B->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train3_B->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);

  reader_train3_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train3_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train3_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train3_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train3_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train3_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train3_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train3_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train3_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train3_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train3_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train3_C->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train3_C->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train3_C->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);





  reader_train3_A->BookMVA( "BDT", basedir_July + "output_3_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train3_B->BookMVA( "BDT", basedir_July + "output_3_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train3_C->BookMVA( "BDT", basedir_July + "output_3_C/weights/TMVAClassification_BDT.weights.xml" );

  std::vector<string> names_train4_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3"};
  std::vector<float>  values_train4_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3};

  reader_train4_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train4_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train4_C = new TMVA::Reader( "!Color:!Silent" );



  reader_train4_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train4_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train4_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train4_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train4_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train4_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train4_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train4_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train4_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train4_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train4_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train4_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train4_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train4_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);


  reader_train4_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train4_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train4_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train4_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train4_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train4_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train4_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train4_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train4_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train4_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train4_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train4_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train4_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train4_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);


  reader_train4_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train4_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train4_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train4_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train4_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train4_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train4_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train4_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train4_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train4_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train4_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train4_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train4_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train4_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);


  reader_train4_A->BookMVA( "BDT", basedir_July + "output_4_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train4_B->BookMVA( "BDT", basedir_July + "output_4_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train4_C->BookMVA( "BDT", basedir_July + "output_4_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train5_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					"var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3"};
  std::vector<float>  values_train5_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3};

  reader_train5_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train5_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train5_C = new TMVA::Reader( "!Color:!Silent" );




  reader_train5_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train5_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train5_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train5_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train5_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train5_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train5_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train5_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train5_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train5_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train5_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train5_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train5_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train5_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train5_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train5_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train5_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);




  reader_train5_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train5_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train5_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train5_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train5_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train5_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train5_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train5_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train5_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train5_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train5_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train5_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train5_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train5_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train5_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train5_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train5_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);



  reader_train5_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train5_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train5_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train5_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train5_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train5_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train5_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train5_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train5_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train5_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train5_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train5_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train5_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train5_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train5_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train5_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train5_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);


  reader_train5_A->BookMVA( "BDT", basedir_July + "output_5_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train5_B->BookMVA( "BDT", basedir_July + "output_5_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train5_C->BookMVA( "BDT", basedir_July + "output_5_C/weights/TMVAClassification_BDT.weights.xml" );




  std::vector<string> names_train6_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					"var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3",
					"var_IsoPhiKKMass_Mu1","var_IsoPhiKKMass_Mu2","var_IsoPhiKKMass_Mu3"};
  std::vector<float>  values_train6_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,
					var_IsoPhiKKMass_Mu1,var_IsoPhiKKMass_Mu2,var_IsoPhiKKMass_Mu3};

  reader_train6_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train6_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train6_C = new TMVA::Reader( "!Color:!Silent" );




  reader_train6_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train6_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train6_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train6_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train6_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train6_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train6_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train6_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train6_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train6_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train6_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train6_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train6_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train6_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train6_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train6_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train6_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train6_A->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train6_A->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train6_A->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);





  reader_train6_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train6_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train6_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train6_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train6_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train6_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train6_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train6_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train6_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train6_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train6_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train6_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train6_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train6_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train6_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train6_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train6_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train6_B->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train6_B->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train6_B->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);





  reader_train6_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train6_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train6_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train6_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train6_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train6_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train6_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train6_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train6_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train6_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train6_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train6_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train6_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train6_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train6_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train6_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train6_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train6_C->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train6_C->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train6_C->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);


  reader_train6_A->BookMVA( "BDT", basedir_July + "output_6_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train6_B->BookMVA( "BDT", basedir_July + "output_6_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train6_C->BookMVA( "BDT", basedir_July + "output_6_C/weights/TMVAClassification_BDT.weights.xml" );



  std::vector<string> names_train7_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					"var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3",
					"var_IsoPhiKKMass_Mu1","var_IsoPhiKKMass_Mu2","var_IsoPhiKKMass_Mu3"};
  std::vector<float>  values_train7_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,
					var_IsoPhiKKMass_Mu1,var_IsoPhiKKMass_Mu2,var_IsoPhiKKMass_Mu3};
  reader_train7_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train7_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train7_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train7_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train7_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train7_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train7_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train7_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train7_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train7_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train7_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train7_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train7_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train7_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train7_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train7_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train7_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train7_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train7_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train7_A->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train7_A->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train7_A->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);


  reader_train7_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train7_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train7_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train7_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train7_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train7_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train7_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train7_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train7_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train7_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train7_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train7_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train7_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train7_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train7_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train7_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train7_B->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train7_B->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train7_B->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);


  reader_train7_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train7_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train7_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train7_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train7_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train7_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train7_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train7_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train7_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train7_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train7_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train7_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train7_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train7_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train7_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train7_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train7_C->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train7_C->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train7_C->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);




  reader_train7_A->BookMVA( "BDT", basedir_July + "output_7_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train7_B->BookMVA( "BDT", basedir_July + "output_7_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train7_C->BookMVA( "BDT", basedir_July + "output_7_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train8_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					"var_MindcaTrackSV","var_Muon1DetID",
					"var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					"var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					"var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3","var_BvsDSeprator"};
  std::vector<float>  values_train8_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					var_MindcaTrackSV,var_Muon1DetID,
					var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,var_BvsDSeprator};

  reader_train8_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train8_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train8_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train8_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train8_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train8_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train8_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train8_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train8_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train8_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train8_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train8_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train8_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train8_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train8_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train8_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train8_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train8_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train8_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train8_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train8_A->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train8_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train8_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train8_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train8_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train8_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train8_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train8_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train8_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train8_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train8_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train8_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train8_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train8_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train8_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train8_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train8_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train8_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train8_B->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);



  reader_train8_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train8_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train8_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train8_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train8_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train8_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train8_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train8_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train8_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train8_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train8_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train8_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train8_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train8_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train8_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train8_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train8_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train8_C->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);






  reader_train8_A->BookMVA( "BDT", basedir_July + "output_8_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train8_B->BookMVA( "BDT", basedir_July + "output_8_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train8_C->BookMVA( "BDT", basedir_July + "output_8_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train9_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					 "var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3","var_Vertex2muTrkKF"};
  std::vector<float>  values_train9_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					 var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,var_Vertex2muTrkKF};
  reader_train9_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train9_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train9_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train9_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train9_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train9_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train9_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train9_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train9_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train9_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train9_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train9_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train9_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train9_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train9_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train9_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train9_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train9_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train9_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train9_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train9_A->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);


  reader_train9_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train9_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train9_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train9_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train9_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train9_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train9_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train9_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train9_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train9_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train9_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train9_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train9_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train9_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train9_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train9_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train9_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train9_B->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);

  reader_train9_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train9_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train9_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train9_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train9_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train9_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train9_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train9_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train9_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train9_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train9_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train9_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train9_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train9_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train9_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train9_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train9_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train9_C->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);



  reader_train9_A->BookMVA( "BDT", basedir_July + "output_9_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train9_B->BookMVA( "BDT", basedir_July + "output_9_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train9_C->BookMVA( "BDT", basedir_July + "output_9_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train10_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					 "var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3","var_Vertex2muTrkKF","var_BvsDSeprator"};
  std::vector<float>  values_train10_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					 var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,var_Vertex2muTrkKF,var_BvsDSeprator};

  reader_train10_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train10_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train10_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train10_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train10_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train10_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train10_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train10_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train10_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train10_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train10_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train10_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train10_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train10_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train10_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train10_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train10_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train10_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train10_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train10_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train10_A->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train10_A->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train10_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train10_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train10_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train10_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train10_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train10_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train10_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train10_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train10_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train10_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train10_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train10_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train10_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train10_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train10_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train10_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train10_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train10_B->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train10_B->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);

  reader_train10_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train10_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train10_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train10_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train10_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train10_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train10_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train10_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train10_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train10_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train10_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train10_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train10_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train10_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train10_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train10_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train10_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);
  reader_train10_C->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train10_C->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);





  reader_train10_A->BookMVA( "BDT", basedir_July + "output_10_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train10_B->BookMVA( "BDT", basedir_July + "output_10_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train10_C->BookMVA( "BDT", basedir_July + "output_10_C/weights/TMVAClassification_BDT.weights.xml" );

  std::vector<string> names_train11_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					 "var_Vertex2muTrkKF","var_BvsDSeprator"};
  std::vector<float>  values_train11_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					 var_Vertex2muTrkKF,var_BvsDSeprator};

  reader_train11_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train11_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train11_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train11_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train11_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train11_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train11_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train11_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train11_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train11_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train11_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train11_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train11_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train11_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train11_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train11_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train11_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train11_A->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train11_A->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);



  reader_train11_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train11_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train11_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train11_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train11_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train11_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train11_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train11_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train11_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train11_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train11_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train11_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train11_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train11_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train11_B->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train11_B->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train11_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train11_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train11_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train11_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train11_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train11_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train11_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train11_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train11_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train11_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train11_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train11_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train11_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train11_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train11_C->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train11_C->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train11_A->BookMVA( "BDT", basedir_July + "output_11_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train11_B->BookMVA( "BDT", basedir_July + "output_11_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train11_C->BookMVA( "BDT", basedir_July + "output_11_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train12_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					 "var_Vertex2muTrkKF","var_BvsDSeprator","var_segCompMuMin","var_MaxMuon_chi2LocalPosition","var_MaxtrkKink"};
  std::vector<float>  values_train12_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					 var_Vertex2muTrkKF,var_BvsDSeprator,var_segCompMuMin,var_MaxMuon_chi2LocalPosition,var_MaxtrkKink};


  reader_train12_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train12_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train12_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train12_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train12_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train12_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train12_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train12_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train12_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train12_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train12_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train12_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train12_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train12_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train12_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train12_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train12_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train12_A->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train12_A->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);
  reader_train12_A->AddVariable("var_segCompMuMin",  &var_segCompMuMin);
  reader_train12_A->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);
  reader_train12_A->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);


  reader_train12_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train12_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train12_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train12_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train12_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train12_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train12_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train12_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train12_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train12_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train12_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train12_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train12_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train12_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train12_B->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train12_B->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);
  reader_train12_B->AddVariable("var_segCompMuMin",  &var_segCompMuMin);
  reader_train12_B->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);
  reader_train12_B->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);

  reader_train12_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train12_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train12_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train12_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train12_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train12_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train12_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train12_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train12_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train12_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train12_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train12_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train12_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train12_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train12_C->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train12_C->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);
  reader_train12_C->AddVariable("var_segCompMuMin",  &var_segCompMuMin);
  reader_train12_C->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);
  reader_train12_C->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);


      
  reader_train12_A->BookMVA( "BDT", basedir_July + "output_12_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train12_B->BookMVA( "BDT", basedir_July + "output_12_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train12_C->BookMVA( "BDT", basedir_July + "output_12_C/weights/TMVAClassification_BDT.weights.xml" );


  std::vector<string> names_train13_A = {"var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_maxMuonsDca","var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3",
					 "var_Vertex2muTrkKF","var_BvsDSeprator"};
  std::vector<float>  values_train13_A = {var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_maxMuonsDca,var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3,
					 var_Vertex2muTrkKF,var_BvsDSeprator};


  reader_train13_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train13_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train13_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train13_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train13_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train13_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train13_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train13_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train13_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train13_A->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train13_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train13_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train13_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train13_A->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train13_A->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train13_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train13_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train13_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train13_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train13_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train13_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train13_B->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train13_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train13_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train13_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train13_B->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train13_B->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);


  reader_train13_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train13_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train13_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train13_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train13_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train13_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train13_C->AddVariable("var_maxMuonsDca",  &var_maxMuonsDca);
  reader_train13_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train13_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train13_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  reader_train13_C->AddVariable("var_Vertex2muTrkKF",  &var_Vertex2muTrkKF);
  reader_train13_C->AddVariable("var_BvsDSeprator",  &var_BvsDSeprator);



  reader_train13_A->BookMVA( "BDT", basedir_July + "output_13_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train13_B->BookMVA( "BDT", basedir_July + "output_13_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train13_C->BookMVA( "BDT", basedir_July + "output_13_C/weights/TMVAClassification_BDT.weights.xml" );
  std::vector<string> names_train14_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_IsoPhiKKMass_Mu1","var_IsoPhiKKMass_Mu2","var_IsoPhiKKMass_Mu3"};
  std::vector<float>  values_train14_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_IsoPhiKKMass_Mu1,var_IsoPhiKKMass_Mu2,var_IsoPhiKKMass_Mu3};


  reader_train14_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train14_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train14_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train14_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train14_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train14_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train14_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train14_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train14_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train14_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train14_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train14_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train14_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train14_A->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train14_A->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train14_A->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);

  reader_train14_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train14_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train14_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train14_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train14_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train14_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train14_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train14_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train14_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train14_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train14_B->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train14_B->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train14_B->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);

  reader_train14_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train14_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train14_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train14_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train14_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train14_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train14_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train14_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train14_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train14_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train14_C->AddVariable("var_IsoPhiKKMass_Mu1",  &var_IsoPhiKKMass_Mu1);
  reader_train14_C->AddVariable("var_IsoPhiKKMass_Mu2",  &var_IsoPhiKKMass_Mu2);
  reader_train14_C->AddVariable("var_IsoPhiKKMass_Mu3",  &var_IsoPhiKKMass_Mu3);




  reader_train14_A->BookMVA( "BDT", basedir_July + "output_14_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train14_B->BookMVA( "BDT", basedir_July + "output_14_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train14_C->BookMVA( "BDT", basedir_July + "output_14_C/weights/TMVAClassification_BDT.weights.xml" );

  std::vector<string> names_train15_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_IsoKStarMass_Mu1","var_IsoKStarMass_Mu2","var_IsoKStarMass_Mu3"};
  std::vector<float>  values_train15_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_IsoKStarMass_Mu1,var_IsoKStarMass_Mu2,var_IsoKStarMass_Mu3};

  reader_train15_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train15_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train15_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train15_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train15_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train15_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train15_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train15_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train15_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train15_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train15_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train15_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train15_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train15_A->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train15_A->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train15_A->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);


  reader_train15_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train15_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train15_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train15_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train15_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train15_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train15_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train15_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train15_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train15_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train15_B->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train15_B->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train15_B->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);


  reader_train15_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train15_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train15_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train15_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train15_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train15_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train15_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train15_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train15_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train15_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train15_C->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  reader_train15_C->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  reader_train15_C->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);

  reader_train15_A->BookMVA( "BDT", basedir_July + "output_15_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train15_B->BookMVA( "BDT", basedir_July + "output_15_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train15_C->BookMVA( "BDT", basedir_July + "output_15_C/weights/TMVAClassification_BDT.weights.xml" );

  std::vector<string> names_train16_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3"};
  std::vector<float>  values_train16_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3};

  reader_train16_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train16_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train16_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train16_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train16_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train16_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train16_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train16_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train16_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train16_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train16_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train16_A->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train16_A->AddVariable("var_Iso08",  &var_Iso08);
  reader_train16_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train16_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train16_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);




  reader_train16_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train16_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train16_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train16_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train16_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train16_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train16_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train16_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train16_B->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train16_B->AddVariable("var_Iso08",  &var_Iso08);
  reader_train16_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train16_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train16_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);




  reader_train16_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train16_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train16_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train16_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train16_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train16_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train16_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train16_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  reader_train16_C->AddVariable("var_NtracksClose",  &var_NtracksClose);
  reader_train16_C->AddVariable("var_Iso08",  &var_Iso08);
  reader_train16_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1);
  reader_train16_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2);
  reader_train16_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3);






  reader_train16_A->BookMVA( "BDT", basedir_July + "output_16_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train16_B->BookMVA( "BDT", basedir_July + "output_16_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train16_C->BookMVA( "BDT", basedir_July + "output_16_C/weights/TMVAClassification_BDT.weights.xml" );

  std::vector<string> names_train17_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_MindcaTrackSV","var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID","var_MaxVertexPairQuality","var_NtracksClose","var_Iso08",
					 "var_IsoMuMuMass_Mu1","var_IsoMuMuMass_Mu2","var_IsoMuMuMass_Mu3","var_segCompMuMin","var_MaxMuon_chi2LocalPosition","var_MaxtrkKink"};
  std::vector<float>  values_train17_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_MindcaTrackSV,var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID,var_MaxVertexPairQuality,var_NtracksClose,var_Iso08,
					 var_IsoMuMuMass_Mu1,var_IsoMuMuMass_Mu2,var_IsoMuMuMass_Mu3,var_segCompMuMin,var_MaxMuon_chi2LocalPosition,var_MaxtrkKink};

  reader_train17_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train17_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train17_C = new TMVA::Reader( "!Color:!Silent" );

  reader_train17_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2 );
  reader_train17_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle );
  reader_train17_A->AddVariable("var_flightLenSig",  &var_flightLenSig );
  reader_train17_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV );
  reader_train17_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID );
  reader_train17_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID );
  reader_train17_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID );
  reader_train17_A->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality );
  reader_train17_A->AddVariable("var_NtracksClose",  &var_NtracksClose );
  reader_train17_A->AddVariable("var_Iso08",  &var_Iso08 );
  reader_train17_A->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1 );
  reader_train17_A->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2 );
  reader_train17_A->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3 );
  reader_train17_A->AddVariable("var_segCompMuMin",  &var_segCompMuMin );
  reader_train17_A->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition );
  reader_train17_A->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink );



  reader_train17_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2 );
  reader_train17_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle );
  reader_train17_B->AddVariable("var_flightLenSig",  &var_flightLenSig );
  reader_train17_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV );
  reader_train17_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID );
  reader_train17_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID );
  reader_train17_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID );
  reader_train17_B->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality );
  reader_train17_B->AddVariable("var_NtracksClose",  &var_NtracksClose );
  reader_train17_B->AddVariable("var_Iso08",  &var_Iso08 );
  reader_train17_B->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1 );
  reader_train17_B->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2 );
  reader_train17_B->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3 );
  reader_train17_B->AddVariable("var_segCompMuMin",  &var_segCompMuMin );
  reader_train17_B->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition );
  reader_train17_B->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink );



  reader_train17_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2 );
  reader_train17_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle );
  reader_train17_C->AddVariable("var_flightLenSig",  &var_flightLenSig );
  reader_train17_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV );
  reader_train17_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID );
  reader_train17_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID );
  reader_train17_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID );
  reader_train17_C->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality );
  reader_train17_C->AddVariable("var_NtracksClose",  &var_NtracksClose );
  reader_train17_C->AddVariable("var_Iso08",  &var_Iso08 );
  reader_train17_C->AddVariable("var_IsoMuMuMass_Mu1",  &var_IsoMuMuMass_Mu1 );
  reader_train17_C->AddVariable("var_IsoMuMuMass_Mu2",  &var_IsoMuMuMass_Mu2 );
  reader_train17_C->AddVariable("var_IsoMuMuMass_Mu3",  &var_IsoMuMuMass_Mu3 );
  reader_train17_C->AddVariable("var_segCompMuMin",  &var_segCompMuMin );
  reader_train17_C->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition );
  reader_train17_C->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink );




  reader_train17_A->BookMVA( "BDT", basedir_July + "output_17_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train17_B->BookMVA( "BDT", basedir_July + "output_17_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train17_C->BookMVA( "BDT", basedir_July + "output_17_C/weights/TMVAClassification_BDT.weights.xml" );




  std::vector<string> names_train18_A = {"var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig",
					 "var_Muon1DetID",
					 "var_Muon2DetID","var_Muon3DetID"};
  std::vector<float>  values_train18_A = {var_vertexKFChi2, var_svpvTauAngle, var_flightLenSig,
					 var_Muon1DetID,
					 var_Muon2DetID,var_Muon3DetID};

  reader_train18_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train18_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train18_C = new TMVA::Reader( "!Color:!Silent" );




  reader_train18_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train18_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train18_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train18_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train18_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train18_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);


  reader_train18_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train18_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train18_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train18_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train18_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train18_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);


  reader_train18_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train18_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train18_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train18_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train18_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train18_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);

  reader_train18_A->BookMVA( "BDT", basedir_July + "output_18_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train18_B->BookMVA( "BDT", basedir_July + "output_18_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train18_C->BookMVA( "BDT", basedir_July + "output_18_C/weights/TMVAClassification_BDT.weights.xml" );




  reader_train19_A = new TMVA::Reader( "!Color:!Silent" );
  reader_train19_B = new TMVA::Reader( "!Color:!Silent" );
  reader_train19_C = new TMVA::Reader( "!Color:!Silent" );


  reader_train19_A->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train19_A->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train19_A->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train19_A->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train19_A->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train19_A->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);
  reader_train19_A->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train19_A->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train19_A->AddVariable("var_segCompMuMax",  &var_segCompMuMax);
  reader_train19_A->AddVariable("var_MaxD0SigBS",  &var_MaxD0SigBS);
  reader_train19_A->AddVariable("var_MinD0SigBS",  &var_MinD0SigBS);
  reader_train19_A->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);



  reader_train19_B->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train19_B->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train19_B->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train19_B->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train19_B->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train19_B->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);
  reader_train19_B->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train19_B->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train19_B->AddVariable("var_segCompMuMax",  &var_segCompMuMax);
  reader_train19_B->AddVariable("var_MaxD0SigBS",  &var_MaxD0SigBS);
  reader_train19_B->AddVariable("var_MinD0SigBS",  &var_MinD0SigBS);
  reader_train19_B->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);


  reader_train19_C->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  reader_train19_C->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  reader_train19_C->AddVariable("var_flightLenSig",  &var_flightLenSig);
  reader_train19_C->AddVariable("var_Muon1DetID",  &var_Muon1DetID);
  reader_train19_C->AddVariable("var_MindcaTrackSV",  &var_MindcaTrackSV);
  reader_train19_C->AddVariable("var_MaxtrkKink",  &var_MaxtrkKink);
  reader_train19_C->AddVariable("var_Muon2DetID",  &var_Muon2DetID);
  reader_train19_C->AddVariable("var_Muon3DetID",  &var_Muon3DetID);
  reader_train19_C->AddVariable("var_segCompMuMax",  &var_segCompMuMax);
  reader_train19_C->AddVariable("var_MaxD0SigBS",  &var_MaxD0SigBS);
  reader_train19_C->AddVariable("var_MinD0SigBS",  &var_MinD0SigBS);
  reader_train19_C->AddVariable("var_MaxMuon_chi2LocalPosition",  &var_MaxMuon_chi2LocalPosition);



  //  for(unsigned int i =0 ; i < names_train19_A.size(); i++)
  //    { 
      //      reader_train19_A->AddVariable(names_train19_A.at(i),  &values_train19_A.at(i));
  //      reader_train19_B->AddVariable(names_train19_A.at(i),  &values_train19_A.at(i));
  //      reader_train19_C->AddVariable(names_train19_A.at(i),  &values_train19_A.at(i));
  //    }



  reader_train19_A->BookMVA( "BDT", basedir_July + "output_19_A/weights/TMVAClassification_BDT.weights.xml" );
  reader_train19_B->BookMVA( "BDT", basedir_July + "output_19_B/weights/TMVAClassification_BDT.weights.xml" );
  reader_train19_C->BookMVA( "BDT", basedir_July + "output_19_C/weights/TMVAClassification_BDT.weights.xml" );

  TString basedir = "";
  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

  //*** defined the bdt reader for event selection; readerA- category A, readerB - category B ...
  readerA = new TMVA::Reader( "!Color:!Silent" );
  readerA->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerA->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerA->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerA->AddVariable("var_MaxD0SigSV", &var_MaxD0SigSV);
  readerA->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerA->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerA->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerA->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerA->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerA->AddVariable("var_Iso08", &var_Iso08);
  readerA->AddVariable("var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerA->AddVariable("var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerA->AddVariable("var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerA->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_0_A/weights/TMVAClassification_BDT.weights.xml"); 




  readerB = new TMVA::Reader( "!Color:!Silent" );
  readerB->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerB->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerB->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerB->AddVariable("var_MaxD0SigSV", &var_MaxD0SigSV);
  readerB->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerB->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerB->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerB->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerB->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerB->AddVariable("var_Iso08", &var_Iso08);
  readerB->AddVariable("var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerB->AddVariable("var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerB->AddVariable("var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerB->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_0_B/weights/TMVAClassification_BDT.weights.xml"); 


  readerC = new TMVA::Reader( "!Color:!Silent" );
  readerC->AddVariable("var_vertexKFChi2",  &var_vertexKFChi2);
  readerC->AddVariable("var_svpvTauAngle",  &var_svpvTauAngle);
  readerC->AddVariable("var_flightLenSig", &var_flightLenSig);
  readerC->AddVariable("var_MaxD0SigSV", &var_MaxD0SigSV);
  readerC->AddVariable("var_MindcaTrackSV", &var_MindcaTrackSV);
  readerC->AddVariable("var_Muon1DetID", &var_Muon1DetID);
  readerC->AddVariable("var_Muon2DetID", &var_Muon2DetID);
  readerC->AddVariable("var_Muon3DetID", &var_Muon3DetID);
  readerC->AddVariable("var_MaxVertexPairQuality", &var_MaxVertexPairQuality);
  readerC->AddVariable("var_Iso08", &var_Iso08);
  readerC->AddVariable("var_IsoKStarMass_Mu1", &var_IsoKStarMass_Mu1);
  readerC->AddVariable("var_IsoKStarMass_Mu2", &var_IsoKStarMass_Mu2);
  readerC->AddVariable("var_IsoKStarMass_Mu3", &var_IsoKStarMass_Mu3);

  //  readerA->AddSpectator("var_tauMass",&var_tauMass);
  readerC->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_0_C/weights/TMVAClassification_BDT.weights.xml"); 





  readerA_train11= new TMVA::Reader( "!Color:!Silent" );
  readerA_train11->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerA_train11->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerA_train11->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerA_train11->AddVariable("var_MindcaTrackSV",         &var_MindcaTrackSV);
  readerA_train11->AddVariable("var_Muon1DetID",            &var_Muon1DetID);
  readerA_train11->AddVariable("var_Muon2DetID",            &var_Muon2DetID);
  readerA_train11->AddVariable("var_Muon3DetID",            &var_Muon3DetID);
  readerA_train11->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  readerA_train11->AddVariable("var_NtracksClose",      &var_NtracksClose);
  readerA_train11->AddVariable("var_Iso08",             &var_Iso08);
  readerA_train11->AddVariable("var_maxMuonsDca",       &var_maxMuonsDca);
  readerA_train11->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  readerA_train11->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  readerA_train11->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  readerA_train11->AddVariable("var_Vertex2muTrkKF",    &var_Vertex2muTrkKF);
  readerA_train11->AddVariable("var_BvsDSeprator",      &var_BvsDSeprator);
  readerA_train11->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_11_A/weights/TMVAClassification_BDT.weights.xml"); 



  readerB_train11= new TMVA::Reader( "!Color:!Silent" );
  readerB_train11->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerB_train11->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerB_train11->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerB_train11->AddVariable("var_MindcaTrackSV",         &var_MindcaTrackSV);
  readerB_train11->AddVariable("var_Muon1DetID",            &var_Muon1DetID);
  readerB_train11->AddVariable("var_Muon2DetID",            &var_Muon2DetID);
  readerB_train11->AddVariable("var_Muon3DetID",            &var_Muon3DetID);
  readerB_train11->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  readerB_train11->AddVariable("var_NtracksClose",      &var_NtracksClose);
  readerB_train11->AddVariable("var_Iso08",             &var_Iso08);
  readerB_train11->AddVariable("var_maxMuonsDca",       &var_maxMuonsDca);
  readerB_train11->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  readerB_train11->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  readerB_train11->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  readerB_train11->AddVariable("var_Vertex2muTrkKF",    &var_Vertex2muTrkKF);
  readerB_train11->AddVariable("var_BvsDSeprator",      &var_BvsDSeprator);
  readerB_train11->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_11_B/weights/TMVAClassification_BDT.weights.xml"); 



  readerC_train11= new TMVA::Reader( "!Color:!Silent" );
  readerC_train11->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerC_train11->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerC_train11->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerC_train11->AddVariable("var_MindcaTrackSV",         &var_MindcaTrackSV);
  readerC_train11->AddVariable("var_Muon1DetID",            &var_Muon1DetID);
  readerC_train11->AddVariable("var_Muon2DetID",            &var_Muon2DetID);
  readerC_train11->AddVariable("var_Muon3DetID",            &var_Muon3DetID);
  readerC_train11->AddVariable("var_MaxVertexPairQuality",  &var_MaxVertexPairQuality);
  readerC_train11->AddVariable("var_NtracksClose",      &var_NtracksClose);
  readerC_train11->AddVariable("var_Iso08",             &var_Iso08);
  readerC_train11->AddVariable("var_maxMuonsDca",       &var_maxMuonsDca);
  readerC_train11->AddVariable("var_IsoKStarMass_Mu1",  &var_IsoKStarMass_Mu1);
  readerC_train11->AddVariable("var_IsoKStarMass_Mu2",  &var_IsoKStarMass_Mu2);
  readerC_train11->AddVariable("var_IsoKStarMass_Mu3",  &var_IsoKStarMass_Mu3);
  readerC_train11->AddVariable("var_Vertex2muTrkKF",    &var_Vertex2muTrkKF);
  readerC_train11->AddVariable("var_BvsDSeprator",      &var_BvsDSeprator);

  readerC_train11->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_11_C/weights/TMVAClassification_BDT.weights.xml"); 




  readerA_train8= new TMVA::Reader( "!Color:!Silent" );
  readerA_train8->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerA_train8->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerA_train8->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerA_train8->AddVariable("var_MindcaTrackSV",          &var_MindcaTrackSV);
  readerA_train8->AddVariable("var_Muon1DetID",          &var_Muon1DetID);
  readerA_train8->AddVariable("var_Muon2DetID",          &var_Muon2DetID);
  readerA_train8->AddVariable("var_Muon3DetID",          &var_Muon3DetID);
  readerA_train8->AddVariable("var_MaxVertexPairQuality",          &var_MaxVertexPairQuality);
  readerA_train8->AddVariable("var_NtracksClose",          &var_NtracksClose);
  readerA_train8->AddVariable("var_Iso08",          &var_Iso08);
  readerA_train8->AddVariable("var_maxMuonsDca",          &var_maxMuonsDca);
  readerA_train8->AddVariable("var_IsoKStarMass_Mu1",          &var_IsoKStarMass_Mu1);
  readerA_train8->AddVariable("var_IsoKStarMass_Mu2",          &var_IsoKStarMass_Mu2);
  readerA_train8->AddVariable("var_IsoKStarMass_Mu3",          &var_IsoKStarMass_Mu3);
  readerA_train8->AddVariable("var_IsoMuMuMass_Mu1",          &var_IsoMuMuMass_Mu1);
  readerA_train8->AddVariable("var_IsoMuMuMass_Mu2",          &var_IsoMuMuMass_Mu2);
  readerA_train8->AddVariable("var_IsoMuMuMass_Mu3",          &var_IsoMuMuMass_Mu3);
  readerA_train8->AddVariable("var_BvsDSeprator",          &var_BvsDSeprator);
  readerA_train8->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_8_A/weights/TMVAClassification_BDT.weights.xml"); 


  readerB_train8= new TMVA::Reader( "!Color:!Silent" );
  readerB_train8->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerB_train8->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerB_train8->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerB_train8->AddVariable("var_MindcaTrackSV",          &var_MindcaTrackSV);
  readerB_train8->AddVariable("var_Muon1DetID",          &var_Muon1DetID);
  readerB_train8->AddVariable("var_Muon2DetID",          &var_Muon2DetID);
  readerB_train8->AddVariable("var_Muon3DetID",          &var_Muon3DetID);
  readerB_train8->AddVariable("var_MaxVertexPairQuality",          &var_MaxVertexPairQuality);
  readerB_train8->AddVariable("var_NtracksClose",          &var_NtracksClose);
  readerB_train8->AddVariable("var_Iso08",          &var_Iso08);
  readerB_train8->AddVariable("var_maxMuonsDca",          &var_maxMuonsDca);
  readerB_train8->AddVariable("var_IsoKStarMass_Mu1",          &var_IsoKStarMass_Mu1);
  readerB_train8->AddVariable("var_IsoKStarMass_Mu2",          &var_IsoKStarMass_Mu2);
  readerB_train8->AddVariable("var_IsoKStarMass_Mu3",          &var_IsoKStarMass_Mu3);
  readerB_train8->AddVariable("var_IsoMuMuMass_Mu1",          &var_IsoMuMuMass_Mu1);
  readerB_train8->AddVariable("var_IsoMuMuMass_Mu2",          &var_IsoMuMuMass_Mu2);
  readerB_train8->AddVariable("var_IsoMuMuMass_Mu3",          &var_IsoMuMuMass_Mu3);
  readerB_train8->AddVariable("var_BvsDSeprator",          &var_BvsDSeprator);
  readerB_train8->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_8_B/weights/TMVAClassification_BDT.weights.xml"); 

  readerC_train8= new TMVA::Reader( "!Color:!Silent" );
  readerC_train8->AddVariable("var_vertexKFChi2",          &var_vertexKFChi2);
  readerC_train8->AddVariable("var_svpvTauAngle",          &var_svpvTauAngle);
  readerC_train8->AddVariable("var_flightLenSig",          &var_flightLenSig);
  readerC_train8->AddVariable("var_MindcaTrackSV",          &var_MindcaTrackSV);
  readerC_train8->AddVariable("var_Muon1DetID",          &var_Muon1DetID);
  readerC_train8->AddVariable("var_Muon2DetID",          &var_Muon2DetID);
  readerC_train8->AddVariable("var_Muon3DetID",          &var_Muon3DetID);
  readerC_train8->AddVariable("var_MaxVertexPairQuality",          &var_MaxVertexPairQuality);
  readerC_train8->AddVariable("var_NtracksClose",          &var_NtracksClose);
  readerC_train8->AddVariable("var_Iso08",          &var_Iso08);
  readerC_train8->AddVariable("var_maxMuonsDca",          &var_maxMuonsDca);
  readerC_train8->AddVariable("var_IsoKStarMass_Mu1",          &var_IsoKStarMass_Mu1);
  readerC_train8->AddVariable("var_IsoKStarMass_Mu2",          &var_IsoKStarMass_Mu2);
  readerC_train8->AddVariable("var_IsoKStarMass_Mu3",          &var_IsoKStarMass_Mu3);
  readerC_train8->AddVariable("var_IsoMuMuMass_Mu1",          &var_IsoMuMuMass_Mu1);
  readerC_train8->AddVariable("var_IsoMuMuMass_Mu2",          &var_IsoMuMuMass_Mu2);
  readerC_train8->AddVariable("var_IsoMuMuMass_Mu3",          &var_IsoMuMuMass_Mu3);
  readerC_train8->AddVariable("var_BvsDSeprator",          &var_BvsDSeprator);
  readerC_train8->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVAApr_19_2021/Code/CommonUtils/IterativeTrain/output_8_C/weights/TMVAClassification_BDT.weights.xml"); 








  //*** Muon MVA ID  readers

  readerMuIDBarrel= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  readerMuIDBarrel->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDBarrel->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDBarrel->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDBarrel->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDBarrel->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDBarrel->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDBarrel->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDBarrel->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDBarrel->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDBarrel->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDBarrel->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDBarrel->BookMVA( "BDT", basedir+"MuonMVA_02may_barrel/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles


  readerMuIDEndcap= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDEndcap->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDEndcap->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDEndcap->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDEndcap->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDEndcap->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDEndcap->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDEndcap->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDEndcap->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDEndcap->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDEndcap->BookMVA( "BDT", basedir+"MuonMVA_02may_endcap/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles


  readerBvsD= new TMVA::Reader( "!Color:!Silent" );
  readerBvsD->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBvsD->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBvsD->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBvsD->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBvsD->AddVariable("var_nsv",&var_nsv);
  readerBvsD->AddVariable("var_MaxD0SigBS",&var_MaxD0SigBS);
  readerBvsD->AddVariable("var_MinD0SigBS",&var_MinD0SigBS);
  readerBvsD->AddVariable("var_Iso08",&var_Iso08);
  readerBvsD->AddVariable("var_dcaTrackPV",&var_dcaTrackPV);
  readerBvsD->AddVariable("var_MinMuonImpactAngle",&var_MinMuonImpactAngle);
  readerBvsD->AddVariable("var_flightLenDist",&var_flightLenDist);
  readerBvsD->BookMVA( "BDTG", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_0_MCTrainA/weights/TMVAClassification_BDTG.weights.xml" );



  readerBTrainA= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainA->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainA->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainA->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainA->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerBTrainA->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainA->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainA->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainA->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainA->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainA->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainA->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainA->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainA->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainA->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainA->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_A_B/weights/TMVAClassification_BDT.weights.xml");


  readerBTrainB= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainB->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainB->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainB->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainB->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerBTrainB->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainB->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainB->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainB->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainB->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainB->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainB->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainB->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainB->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainB->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainB->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_B_B/weights/TMVAClassification_BDT.weights.xml");



  readerBTrainC= new TMVA::Reader( "!Color:!Silent" );
  readerBTrainC->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerBTrainC->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBTrainC->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBTrainC->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerBTrainC->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBTrainC->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerBTrainC->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerBTrainC->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerBTrainC->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerBTrainC->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBTrainC->AddVariable("var_Iso08",&var_Iso08);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerBTrainC->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerBTrainC->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerBTrainC->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerBTrainC->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_C_B/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainA= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainA->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainA->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainA->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainA->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerDTrainA->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainA->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainA->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainA->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainA->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainA->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainA->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainA->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainA->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainA->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainA->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_A_DS/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainB= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainB->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainB->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainB->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainB->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerDTrainB->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainB->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainB->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainB->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainB->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainB->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainB->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainB->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainB->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainB->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainB->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_B_DS/weights/TMVAClassification_BDT.weights.xml");




  readerDTrainC= new TMVA::Reader( "!Color:!Silent" );
  readerDTrainC->AddVariable("var_vertexKFChi2",&var_vertexKFChi2);
  readerDTrainC->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerDTrainC->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerDTrainC->AddVariable("var_MaxD0SigSV",&var_MaxD0SigSV);
  readerDTrainC->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerDTrainC->AddVariable("var_Muon1DetID",&var_Muon1DetID);
  readerDTrainC->AddVariable("var_Muon2DetID",&var_Muon2DetID);
  readerDTrainC->AddVariable("var_Muon3DetID",&var_Muon3DetID);
  readerDTrainC->AddVariable("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  readerDTrainC->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerDTrainC->AddVariable("var_Iso08",&var_Iso08);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2);
  readerDTrainC->AddVariable("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2);
  readerDTrainC->AddVariable("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2);
  readerDTrainC->AddVariable("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3);
  readerDTrainC->BookMVA( "BDT", "/afs/cern.ch/work/c/cherepan/Analysis/workdirMakeMVADec_14_2020/Code/CommonUtils/IterativeTrain/output_7_C_DS/weights/TMVAClassification_BDT.weights.xml");



  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1T)                cut.at(L1T)=1;
    if(i==HLT)                cut.at(HLT)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==PhiVeto1)           cut.at(PhiVeto1)=0; // defined below
    if(i==OmegaVeto1)         cut.at(OmegaVeto1)=0; // defined below
    if(i==PhiVeto2)           cut.at(PhiVeto2)=0; // defined below
    if(i==OmegaVeto2)         cut.at(OmegaVeto2)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==TauMassCut)         cut.at(TauMassCut)=1;// true for MC and mass side band for data
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;

    if(i==L1T){
      title.at(i)="L1T trigger ";
      hlabel="Level 1 Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1T_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1T_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
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

    else if(i==MuonID){
      title.at(i)="Muons GL and PF";
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

    else if(i==TauMassCut){
      title.at(i)="$\\tau$ mass 1.6 - 2 GeV ";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,60,2.1,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,60,1.4,2.2,hlabel,"Events"));
    }

  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms


  EtaGenSource=HConfig.GetTH1D(Name+"_EtaGenSource","EtaGenSource",2,-0.5,1.5," 1 - D, 2 - D_{s} ","Events");

  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassResolutionHelixRefit=HConfig.GetTH1D(Name+"_TauMassResolutionHelixRefit","TauMassResolutionHelixRefit",50,-0.2,0.2,"Helix refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  TauMass_all_nophiVeto =HConfig.GetTH2D(Name+"_TauMass_all_nophiVeto","3#mu mass vs phimass ",60,1.5,2.1,50,0.8,1.2,"3#mu mass, GeV","#phi mass, GeV");
  TauMass_all =HConfig.GetTH1D(Name+"_TauMass_all","3#mu  mass",60,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMass_allVsBDTA=HConfig.GetTH2D(Name+"_TauMass_allVsBDTA","3#mu mass vs BDTa",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTB=HConfig.GetTH2D(Name+"_TauMass_allVsBDTB","3#mu mass vs BDTb",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");
  TauMass_allVsBDTC=HConfig.GetTH2D(Name+"_TauMass_allVsBDTC","3#mu mass vs BDTc",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","BDT");

  TauMass_vs_BDT_train0_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train0_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train0_A");
  TauMass_vs_BDT_train1_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train1_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train1_A");
  TauMass_vs_BDT_train2_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train2_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train2_A");
  TauMass_vs_BDT_train3_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train3_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train3_A");
  TauMass_vs_BDT_train4_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train4_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train4_A");
  TauMass_vs_BDT_train5_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train5_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train5_A");
  TauMass_vs_BDT_train6_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train6_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train6_A");
  TauMass_vs_BDT_train7_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train7_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train7_A");
  TauMass_vs_BDT_train8_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train8_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train8_A");
  TauMass_vs_BDT_train9_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train9_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train9_A");
  TauMass_vs_BDT_train10_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train10_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train10_A");
  TauMass_vs_BDT_train11_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train11_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train11_A");
  TauMass_vs_BDT_train12_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train12_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train12_A");
  TauMass_vs_BDT_train13_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train13_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train13_A");
  TauMass_vs_BDT_train14_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train14_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train14_A");
  TauMass_vs_BDT_train15_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train15_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train15_A");
  TauMass_vs_BDT_train16_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train16_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train16_A");
  TauMass_vs_BDT_train17_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train17_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train17_A");
  TauMass_vs_BDT_train18_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train18_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train18_A");
  TauMass_vs_BDT_train19_A=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train19_A","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train19_A");



  TauMass_vs_BDT_train0_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train0_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train0_B");
  TauMass_vs_BDT_train1_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train1_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train1_B");
  TauMass_vs_BDT_train2_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train2_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train2_B");
  TauMass_vs_BDT_train3_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train3_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train3_B");
  TauMass_vs_BDT_train4_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train4_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train4_B");
  TauMass_vs_BDT_train5_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train5_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train5_B");
  TauMass_vs_BDT_train6_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train6_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train6_B");
  TauMass_vs_BDT_train7_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train7_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train7_B");
  TauMass_vs_BDT_train8_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train8_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train8_B");
  TauMass_vs_BDT_train9_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train9_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train9_B");
  TauMass_vs_BDT_train10_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train10_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train10_B");
  TauMass_vs_BDT_train11_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train11_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train11_B");
  TauMass_vs_BDT_train12_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train12_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train12_B");
  TauMass_vs_BDT_train13_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train13_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train13_B");
  TauMass_vs_BDT_train14_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train14_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train14_B");
  TauMass_vs_BDT_train15_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train15_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train15_B");
  TauMass_vs_BDT_train16_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train16_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train16_B");
  TauMass_vs_BDT_train17_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train17_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train17_B");
  TauMass_vs_BDT_train18_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train18_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train18_B");
  TauMass_vs_BDT_train19_B=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train19_B","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train19_B");


  TauMass_vs_BDT_train0_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train0_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train0_C");
  TauMass_vs_BDT_train1_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train1_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train1_C");
  TauMass_vs_BDT_train2_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train2_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train2_C");
  TauMass_vs_BDT_train3_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train3_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train3_C");
  TauMass_vs_BDT_train4_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train4_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train4_C");
  TauMass_vs_BDT_train5_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train5_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train5_C");
  TauMass_vs_BDT_train6_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train6_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train6_C");
  TauMass_vs_BDT_train7_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train7_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train7_C");
  TauMass_vs_BDT_train8_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train8_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train8_C");
  TauMass_vs_BDT_train9_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train9_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train9_C");
  TauMass_vs_BDT_train10_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train10_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train10_C");
  TauMass_vs_BDT_train11_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train11_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train11_C");
  TauMass_vs_BDT_train12_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train12_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train12_C");
  TauMass_vs_BDT_train13_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train13_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train13_C");
  TauMass_vs_BDT_train14_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train14_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train14_C");
  TauMass_vs_BDT_train15_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train15_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train15_C");
  TauMass_vs_BDT_train16_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train16_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train16_C");
  TauMass_vs_BDT_train17_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train17_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train17_C");
  TauMass_vs_BDT_train18_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train18_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train18_C");
  TauMass_vs_BDT_train19_C=HConfig.GetTH2D(Name+"_TauMass_vs_BDT_train19_C","3#mu mass vs BDT",60,1.5,2.1,50,-0.3,0.3,"3#mu mass, GeV","TauMass_vs_BDT_train19_C");


  BDTOutput_train0_A= HConfig.GetTH1D(Name+"_BDTOutput_train0_A","BDTOutput_train0_A",40,-0.4,0.4,"BDTOutput_train0_A","Events");
  BDTOutput_train1_A= HConfig.GetTH1D(Name+"_BDTOutput_train1_A","BDTOutput_train1_A",40,-0.4,0.4,"BDTOutput_train1_A","Events");
  BDTOutput_train2_A= HConfig.GetTH1D(Name+"_BDTOutput_train2_A","BDTOutput_train2_A",40,-0.4,0.4,"BDTOutput_train2_A","Events");
  BDTOutput_train3_A= HConfig.GetTH1D(Name+"_BDTOutput_train3_A","BDTOutput_train3_A",40,-0.4,0.4,"BDTOutput_train3_A","Events");
  BDTOutput_train4_A= HConfig.GetTH1D(Name+"_BDTOutput_train4_A","BDTOutput_train4_A",40,-0.4,0.4,"BDTOutput_train4_A","Events");
  BDTOutput_train5_A= HConfig.GetTH1D(Name+"_BDTOutput_train5_A","BDTOutput_train5_A",40,-0.4,0.4,"BDTOutput_train5_A","Events");
  BDTOutput_train6_A= HConfig.GetTH1D(Name+"_BDTOutput_train6_A","BDTOutput_train6_A",40,-0.4,0.4,"BDTOutput_train6_A","Events");
  BDTOutput_train7_A= HConfig.GetTH1D(Name+"_BDTOutput_train7_A","BDTOutput_train7_A",40,-0.4,0.4,"BDTOutput_train7_A","Events");
  BDTOutput_train8_A= HConfig.GetTH1D(Name+"_BDTOutput_train8_A","BDTOutput_train8_A",40,-0.4,0.4,"BDTOutput_train8_A","Events");
  BDTOutput_train9_A= HConfig.GetTH1D(Name+"_BDTOutput_train9_A","BDTOutput_train9_A",40,-0.4,0.4,"BDTOutput_train9_A","Events");
  BDTOutput_train10_A= HConfig.GetTH1D(Name+"_BDTOutput_train10_A","BDTOutput_train10_A",40,-0.4,0.4,"BDTOutput_train10_A","Events");
  BDTOutput_train11_A= HConfig.GetTH1D(Name+"_BDTOutput_train11_A","BDTOutput_train11_A",40,-0.4,0.4,"BDTOutput_train11_A","Events");
  BDTOutput_train12_A= HConfig.GetTH1D(Name+"_BDTOutput_train12_A","BDTOutput_train12_A",40,-0.4,0.4,"BDTOutput_train12_A","Events");
  BDTOutput_train13_A= HConfig.GetTH1D(Name+"_BDTOutput_train13_A","BDTOutput_train13_A",40,-0.4,0.4,"BDTOutput_train13_A","Events");
  BDTOutput_train14_A= HConfig.GetTH1D(Name+"_BDTOutput_train14_A","BDTOutput_train14_A",40,-0.4,0.4,"BDTOutput_train14_A","Events");
  BDTOutput_train15_A= HConfig.GetTH1D(Name+"_BDTOutput_train15_A","BDTOutput_train15_A",40,-0.4,0.4,"BDTOutput_train15_A","Events");
  BDTOutput_train16_A= HConfig.GetTH1D(Name+"_BDTOutput_train16_A","BDTOutput_train16_A",40,-0.4,0.4,"BDTOutput_train16_A","Events");
  BDTOutput_train17_A= HConfig.GetTH1D(Name+"_BDTOutput_train17_A","BDTOutput_train17_A",40,-0.4,0.4,"BDTOutput_train17_A","Events");
  BDTOutput_train18_A= HConfig.GetTH1D(Name+"_BDTOutput_train18_A","BDTOutput_train18_A",40,-0.4,0.4,"BDTOutput_train18_A","Events");
  BDTOutput_train19_A= HConfig.GetTH1D(Name+"_BDTOutput_train19_A","BDTOutput_train19_A",40,-0.4,0.4,"BDTOutput_train19_A","Events");



  BDTOutput_train0_B= HConfig.GetTH1D(Name+"_BDTOutput_train0_B","BDTOutput_train0_B",40,-0.4,0.4,"BDTOutput_train0_B","Events");
  BDTOutput_train1_B= HConfig.GetTH1D(Name+"_BDTOutput_train1_B","BDTOutput_train1_B",40,-0.4,0.4,"BDTOutput_train1_B","Events");
  BDTOutput_train2_B= HConfig.GetTH1D(Name+"_BDTOutput_train2_B","BDTOutput_train2_B",40,-0.4,0.4,"BDTOutput_train2_B","Events");
  BDTOutput_train3_B= HConfig.GetTH1D(Name+"_BDTOutput_train3_B","BDTOutput_train3_B",40,-0.4,0.4,"BDTOutput_train3_B","Events");
  BDTOutput_train4_B= HConfig.GetTH1D(Name+"_BDTOutput_train4_B","BDTOutput_train4_B",40,-0.4,0.4,"BDTOutput_train4_B","Events");
  BDTOutput_train5_B= HConfig.GetTH1D(Name+"_BDTOutput_train5_B","BDTOutput_train5_B",40,-0.4,0.4,"BDTOutput_train5_B","Events");
  BDTOutput_train6_B= HConfig.GetTH1D(Name+"_BDTOutput_train6_B","BDTOutput_train6_B",40,-0.4,0.4,"BDTOutput_train6_B","Events");
  BDTOutput_train7_B= HConfig.GetTH1D(Name+"_BDTOutput_train7_B","BDTOutput_train7_B",40,-0.4,0.4,"BDTOutput_train7_B","Events");
  BDTOutput_train8_B= HConfig.GetTH1D(Name+"_BDTOutput_train8_B","BDTOutput_train8_B",40,-0.4,0.4,"BDTOutput_train8_B","Events");
  BDTOutput_train9_B= HConfig.GetTH1D(Name+"_BDTOutput_train9_B","BDTOutput_train9_B",40,-0.4,0.4,"BDTOutput_train9_B","Events");
  BDTOutput_train10_B= HConfig.GetTH1D(Name+"_BDTOutput_train10_B","BDTOutput_train10_B",40,-0.4,0.4,"BDTOutput_train10_B","Events");
  BDTOutput_train11_B= HConfig.GetTH1D(Name+"_BDTOutput_train11_B","BDTOutput_train11_B",40,-0.4,0.4,"BDTOutput_train11_B","Events");
  BDTOutput_train12_B= HConfig.GetTH1D(Name+"_BDTOutput_train12_B","BDTOutput_train12_B",40,-0.4,0.4,"BDTOutput_train12_B","Events");
  BDTOutput_train13_B= HConfig.GetTH1D(Name+"_BDTOutput_train13_B","BDTOutput_train13_B",40,-0.4,0.4,"BDTOutput_train13_B","Events");
  BDTOutput_train14_B= HConfig.GetTH1D(Name+"_BDTOutput_train14_B","BDTOutput_train14_B",40,-0.4,0.4,"BDTOutput_train14_B","Events");
  BDTOutput_train15_B= HConfig.GetTH1D(Name+"_BDTOutput_train15_B","BDTOutput_train15_B",40,-0.4,0.4,"BDTOutput_train15_B","Events");
  BDTOutput_train16_B= HConfig.GetTH1D(Name+"_BDTOutput_train16_B","BDTOutput_train16_B",40,-0.4,0.4,"BDTOutput_train16_B","Events");
  BDTOutput_train17_B= HConfig.GetTH1D(Name+"_BDTOutput_train17_B","BDTOutput_train17_B",40,-0.4,0.4,"BDTOutput_train17_B","Events");
  BDTOutput_train18_B= HConfig.GetTH1D(Name+"_BDTOutput_train18_B","BDTOutput_train18_B",40,-0.4,0.4,"BDTOutput_train18_B","Events");
  BDTOutput_train19_B= HConfig.GetTH1D(Name+"_BDTOutput_train19_B","BDTOutput_train19_B",40,-0.4,0.4,"BDTOutput_train19_B","Events");


  BDTOutput_train0_C= HConfig.GetTH1D(Name+"_BDTOutput_train0_C","BDTOutput_train0_C",40,-0.4,0.4,"BDTOutput_train0_C","Events");
  BDTOutput_train1_C= HConfig.GetTH1D(Name+"_BDTOutput_train1_C","BDTOutput_train1_C",40,-0.4,0.4,"BDTOutput_train1_C","Events");
  BDTOutput_train2_C= HConfig.GetTH1D(Name+"_BDTOutput_train2_C","BDTOutput_train2_C",40,-0.4,0.4,"BDTOutput_train2_C","Events");
  BDTOutput_train3_C= HConfig.GetTH1D(Name+"_BDTOutput_train3_C","BDTOutput_train3_C",40,-0.4,0.4,"BDTOutput_train3_C","Events");
  BDTOutput_train4_C= HConfig.GetTH1D(Name+"_BDTOutput_train4_C","BDTOutput_train4_C",40,-0.4,0.4,"BDTOutput_train4_C","Events");
  BDTOutput_train5_C= HConfig.GetTH1D(Name+"_BDTOutput_train5_C","BDTOutput_train5_C",40,-0.4,0.4,"BDTOutput_train5_C","Events");
  BDTOutput_train6_C= HConfig.GetTH1D(Name+"_BDTOutput_train6_C","BDTOutput_train6_C",40,-0.4,0.4,"BDTOutput_train6_C","Events");
  BDTOutput_train7_C= HConfig.GetTH1D(Name+"_BDTOutput_train7_C","BDTOutput_train7_C",40,-0.4,0.4,"BDTOutput_train7_C","Events");
  BDTOutput_train8_C= HConfig.GetTH1D(Name+"_BDTOutput_train8_C","BDTOutput_train8_C",40,-0.4,0.4,"BDTOutput_train8_C","Events");
  BDTOutput_train9_C= HConfig.GetTH1D(Name+"_BDTOutput_train9_C","BDTOutput_train9_C",40,-0.4,0.4,"BDTOutput_train9_C","Events");
  BDTOutput_train10_C= HConfig.GetTH1D(Name+"_BDTOutput_train10_C","BDTOutput_train10_C",40,-0.4,0.4,"BDTOutput_train10_C","Events");
  BDTOutput_train11_C= HConfig.GetTH1D(Name+"_BDTOutput_train11_C","BDTOutput_train11_C",40,-0.4,0.4,"BDTOutput_train11_C","Events");
  BDTOutput_train12_C= HConfig.GetTH1D(Name+"_BDTOutput_train12_C","BDTOutput_train12_C",40,-0.4,0.4,"BDTOutput_train12_C","Events");
  BDTOutput_train13_C= HConfig.GetTH1D(Name+"_BDTOutput_train13_C","BDTOutput_train13_C",40,-0.4,0.4,"BDTOutput_train13_C","Events");
  BDTOutput_train14_C= HConfig.GetTH1D(Name+"_BDTOutput_train14_C","BDTOutput_train14_C",40,-0.4,0.4,"BDTOutput_train14_C","Events");
  BDTOutput_train15_C= HConfig.GetTH1D(Name+"_BDTOutput_train15_C","BDTOutput_train15_C",40,-0.4,0.4,"BDTOutput_train15_C","Events");
  BDTOutput_train16_C= HConfig.GetTH1D(Name+"_BDTOutput_train16_C","BDTOutput_train16_C",40,-0.4,0.4,"BDTOutput_train16_C","Events");
  BDTOutput_train17_C= HConfig.GetTH1D(Name+"_BDTOutput_train17_C","BDTOutput_train17_C",40,-0.4,0.4,"BDTOutput_train17_C","Events");
  BDTOutput_train18_C= HConfig.GetTH1D(Name+"_BDTOutput_train18_C","BDTOutput_train18_C",40,-0.4,0.4,"BDTOutput_train18_C","Events");
  BDTOutput_train19_C= HConfig.GetTH1D(Name+"_BDTOutput_train19_C","BDTOutput_train19_C",40,-0.4,0.4,"BDTOutput_train19_C","Events");




  EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

  TauMassA1 =HConfig.GetTH1D(Name+"_TauMassA1","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitA1 =HConfig.GetTH1D(Name+"_TauMassRefitA1","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A1)","Events");






  TauMassRefitABC1 =HConfig.GetTH1D(Name+"_TauMassRefitABC1","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2 =HConfig.GetTH1D(Name+"_TauMassRefitABC2","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2)","Events");




  TauMassRefitABC1_train11 =HConfig.GetTH1D(Name+"_TauMassRefitABC1_train11","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1) (train 11)","Events");
  TauMassRefitABC2_train11 =HConfig.GetTH1D(Name+"_TauMassRefitABC2_train11","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2) (train 11)","Events");
  TauMassRefitA1_train11 =HConfig.GetTH1D(Name+"_TauMassRefitA1_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A1) (train 11)","Events");
  TauMassRefitB1_train11 =HConfig.GetTH1D(Name+"_TauMassRefitB1_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B1) (train 11)","Events");
  TauMassRefitC1_train11 =HConfig.GetTH1D(Name+"_TauMassRefitC1_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C1) (train 11)","Events");
  TauMassRefitA2_train11 =HConfig.GetTH1D(Name+"_TauMassRefitA2_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A2) (train 11)","Events");
  TauMassRefitB2_train11 =HConfig.GetTH1D(Name+"_TauMassRefitB2_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B2) (train 11)","Events");
  TauMassRefitC2_train11 =HConfig.GetTH1D(Name+"_TauMassRefitC2_train11","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C2) (train 11)","Events");
  AllignSortMass1_train11=HConfig.GetTH1D(Name+"_AllignSortMass1_train11","AllignSortMass1",40,0.2,1.777,"M_{1} (#Delta R OS sorted) (train 11), GeV","");
  AllignSortMass2_train11=HConfig.GetTH1D(Name+"_AllignSortMass2_train11","AllignSortMass2",40,0.2,1.777,"M_{2} (#Delta R OS sorted) (train 11), GeV","");

  BDTOutputA_train11 = HConfig.GetTH1D(Name+"_BDTOutputA_train11","BDTOutputA",40,-0.4,0.4,"BDT Output (train 11) cat A","Events");
  BDTOutputB_train11 = HConfig.GetTH1D(Name+"_BDTOutputB_train11","BDTOutputB",40,-0.4,0.4,"BDT Output (train 11) cat B","Events");
  BDTOutputC_train11 = HConfig.GetTH1D(Name+"_BDTOutputC_train11","BDTOutputC",40,-0.4,0.4,"BDT Output (train 11) cat C","Events");

  TauMassRefitABC1_train8 =HConfig.GetTH1D(Name+"_TauMassRefitABC1_train8","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1) (train 8)","Events");
  TauMassRefitABC2_train8 =HConfig.GetTH1D(Name+"_TauMassRefitABC2_train8","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2) (train 8)","Events");
  TauMassRefitA1_train8 =HConfig.GetTH1D(Name+"_TauMassRefitA1_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A1) (train 8)","Events");
  TauMassRefitB1_train8 =HConfig.GetTH1D(Name+"_TauMassRefitB1_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B1) (train 8)","Events");
  TauMassRefitC1_train8 =HConfig.GetTH1D(Name+"_TauMassRefitC1_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C1) (train 8)","Events");
  TauMassRefitA2_train8 =HConfig.GetTH1D(Name+"_TauMassRefitA2_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A2) (train 8)","Events");
  TauMassRefitB2_train8 =HConfig.GetTH1D(Name+"_TauMassRefitB2_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B2) (train 8)","Events");
  TauMassRefitC2_train8 =HConfig.GetTH1D(Name+"_TauMassRefitC2_train8","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C2) (train 8)","Events");
  AllignSortMass1_train8=HConfig.GetTH1D(Name+"_AllignSortMass1_train8","AllignSortMass1",40,0.2,1.777,"M_{1} (#Delta R OS sorted) (train 8), GeV","");
  AllignSortMass2_train8=HConfig.GetTH1D(Name+"_AllignSortMass2_train8","AllignSortMass2",40,0.2,1.777,"M_{2} (#Delta R OS sorted) (train 8), GeV","");

  BDTOutputA_train8 = HConfig.GetTH1D(Name+"_BDTOutputA_train8","BDTOutputA",40,-0.4,0.4,"BDT Output (train 8) cat A","Events");
  BDTOutputB_train8 = HConfig.GetTH1D(Name+"_BDTOutputB_train8","BDTOutputB",40,-0.4,0.4,"BDT Output (train 8) cat B","Events");
  BDTOutputC_train8 = HConfig.GetTH1D(Name+"_BDTOutputC_train8","BDTOutputC",40,-0.4,0.4,"BDT Output (train 8) cat C","Events");









  TauMassRefitABC1_BDSeparateTrain=HConfig.GetTH1D(Name+"_TauMassRefitABC1_BDSeparateTrain","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC1)","Events");
  TauMassRefitABC2_BDSeparateTrain=HConfig.GetTH1D(Name+"_TauMassRefitABC2_BDSeparateTrain","Refit #tau lepton mass",25,1.5,2.1,"M_{3#mu} , GeV (inclusive ABC2)","Events");



  TauMassRefitABC1_eta =HConfig.GetTH2D(Name+"_TauMassRefitABC1_eta","Refit #tau lepton mass vs eta",25,1.5,2.1,30,0,2.5,"M_{3#mu} , GeV (inclusive ABC1)","#eta_{#tau}");
  TauMassRefitABC2_eta =HConfig.GetTH2D(Name+"_TauMassRefitABC2_eta","Refit #tau lepton mass vs eta",25,1.5,2.1,30,0,2.5,"M_{3#mu} , GeV (inclusive ABC2)","#eta_{#tau}");


  TauMassB1 =HConfig.GetTH1D(Name+"_TauMassB1","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitB1 =HConfig.GetTH1D(Name+"_TauMassRefitB1","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B1)","Events");









  TauMassRefitABC1FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitABC1FullEtaVetoCut","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (ABC1)","Events");
  TauMassRefitABC2FullEtaVetoCut =HConfig.GetTH1D(Name+"_TauMassRefitABC2FullEtaVetoCut","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau}, GeV (#eta veto) (ABC2)","Events");





  TauMassC1 =HConfig.GetTH1D(Name+"_TauMassC1","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV","Events");
  TauMassRefitC1 =HConfig.GetTH1D(Name+"_TauMassRefitC1","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C1)","Events");



  TauMassA2 =HConfig.GetTH1D(Name+"_TauMassA2","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV (A2)","Events");
  TauMassRefitA2 =HConfig.GetTH1D(Name+"_TauMassRefitA2","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (A2)","Events");


  TauMassB2 =HConfig.GetTH1D(Name+"_TauMassB2","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV (B2)","Events");
  TauMassRefitB2 =HConfig.GetTH1D(Name+"_TauMassRefitB2","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (B2)","Events");


  TauMassC2 =HConfig.GetTH1D(Name+"_TauMassC2","#tau lepton mass",25,1.5,2.1,"  M_{#tau} , GeV (C2)","Events");
  TauMassRefitC2 =HConfig.GetTH1D(Name+"_TauMassRefitC2","Refit #tau lepton mass",25,1.5,2.1,"KF refit  M_{#tau} , GeV (C2)","Events");
  //TauMassRefitC2FullEtaVetoCut

  EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
  EventMassResolution_PtEtaPhi_TauEta1p2 = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi_TauEta1p2","EventMassResolution_PtEtaPhi_TauEta1p2",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi) |#eta_{#tau}| < 1.2","Events");

  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
  VertexChi2KF_vs_HelixFit=HConfig.GetTH2D(Name+"_VertexChi2KF_vs_HelixFit","VertexChi2KF_vs_HelixFit",100,0,100,100,0,100,"Kalman Vertex #chi^{2}","Helix Vertex  Fitter #chi^{2}");

  KF_Helix_deltaX=HConfig.GetTH1D(Name+"_KF_Helix_deltaX","KF_Helix_deltaX",50,-0.05,0.05,"#Delta X, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaY=HConfig.GetTH1D(Name+"_KF_Helix_deltaY","KF_Helix_deltaY",50,-0.05,0.05,"#Delta Y, cm (Helix Fitter - Kalman Fitter)","Events");
  KF_Helix_deltaZ=HConfig.GetTH1D(Name+"_KF_Helix_deltaZ","KF_Helix_deltaZ",50,-0.05,0.05,"#Delta Z, cm (Helix Fitter - Kalman Fitter)","Events");

  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",50,0,15,"PV - SV distance  significance","Events");

  SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",50,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.1,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.1,"reco - mc #mu_{2} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.1,"reco - mc #mu_{3} #Delta R","Events");

  TriggerMatchdR1 =HConfig.GetTH1D(Name+"_TriggerMatchdR1","TriggerMatchdR1",50,0,1,"trigger match #Delta R 1","Events");
  TriggerMatchdR2 =HConfig.GetTH1D(Name+"_TriggerMatchdR2","TriggerMatchdR2",50,0,1,"trigger match #Delta R 2","Events");
  TriggerMatchdR3 =HConfig.GetTH1D(Name+"_TriggerMatchdR3","TriggerMatchdR3",50,0,1,"trigger match #Delta R 3","Events");

  BDTOutputA = HConfig.GetTH1D(Name+"_BDTOutputA","BDTOutputA",40,-0.4,0.4,"BDT Output (train 2) cat A","Events");
  BDTOutputB = HConfig.GetTH1D(Name+"_BDTOutputB","BDTOutputB",40,-0.4,0.4,"BDT Output (train 2) cat B","Events");
  BDTOutputC = HConfig.GetTH1D(Name+"_BDTOutputC","BDTOutputC",40,-0.4,0.4,"BDT Output (train 2) cat C","Events");

  BvsDBDTG  = HConfig.GetTH1D(Name+"_BvsDBDTG","BvsDBDTG",50,-1.0,1.0," B vs D BDTG","Events");

  BvsDBDTG_ABC1  = HConfig.GetTH1D(Name+"_BvsDBDTG_ABC1","BvsDBDTG_ABC1",40,-1.0,1.0," B vs D BDTG (ABC1 inclusive)","Events");
  BvsDBDTG_ABC2  = HConfig.GetTH1D(Name+"_BvsDBDTG_ABC2","BvsDBDTG_ABC2",40,-1.0,1.0," B vs D BDTG (ABC2 inclusive)","Events");

  NSignalCandidates =HConfig.GetTH1D(Name+"_NSignalCandidates","NSignalCandidates",5,-0.5,4.5,"Number of signal candidates","Events");
  PairMass=HConfig.GetTH2D(Name+"_PairMass","PairMass",100,0.2,1.8,100,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");


 
  KpiIsolationMass_OS=HConfig.GetTH1D(Name+"_KpiIsolationMass_OS","KpiIsolationMass_OS",100,0.6,1.8,"M_{1}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS1=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS1","KpiIsolationMass_SS1",100,0.6,1.8,"M_{2}(K#pi), GeV (comb. iso #pi)","");
  KpiIsolationMass_SS2=HConfig.GetTH1D(Name+"_KpiIsolationMass_SS2","KpiIsolationMass_SS2",100,0.6,1.8,"M_{3}(K#pi), GeV (comb. iso #pi)","");
 

  PairMass1NoSorting=HConfig.GetTH1D(Name+"_PairMass1NoSorting","PairMass1NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","Events");
  PairMass2NoSorting=HConfig.GetTH1D(Name+"_PairMass2NoSorting","PairMass2NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2 no sorting), GeV","Events");
  MuMuMassNoSorting=HConfig.GetTH2D(Name+"_MuMuMassNoSorting","MuMuMassNoSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","M_{#mu#mu} (OS-SS, 2 no sorting sorting), GeV");


  PairMass1PTSorting=HConfig.GetTH1D(Name+"_PairMass1PTSorting","PairMass1PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","Events");
  PairMass2PTSorting=HConfig.GetTH1D(Name+"_PairMass2PTSorting","PairMass2PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV","Events");
  MuMuMassPTSorting=HConfig.GetTH2D(Name+"_MuMuMassPTSorting","MuMuMassPTSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV");


  PairMass1AllignedSorting=HConfig.GetTH1D(Name+"_PairMass1AllignedSorting","PairMass1AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV","Events");
  PairMass2AllignedSorting=HConfig.GetTH1D(Name+"_PairMass2AllignedSorting","PairMass2AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","Events");
  MuMuMassAllignedSorting=HConfig.GetTH2D(Name+"_MuMuMassAllignedSorting","MuMuMassAllignedSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","M_{#mu#mu} (OS-SS, 1st collimated pair) GeV");



  PairMassDRSorted1A=HConfig.GetTH1D(Name+"_PairMassDRSorted1A","PairMassDRSorted1A",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 1st pair), GeV; cat A","Events");
  PairMassDRSorted2A=HConfig.GetTH1D(Name+"_PairMassDRSorted2A","PairMassDRSorted2A",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 2nd pair), GeV; cat A","Events");
  PairMassDRSorted1B=HConfig.GetTH1D(Name+"_PairMassDRSorted1B","PairMassDRSorted1B",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 1st pair), GeV; cat B","Events");
  PairMassDRSorted2B=HConfig.GetTH1D(Name+"_PairMassDRSorted2B","PairMassDRSorted2B",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 2nd pair), GeV; cat B","Events");
  PairMassDRSorted1C=HConfig.GetTH1D(Name+"_PairMassDRSorted1C","PairMassDRSorted1C",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 1st pair), GeV; cat C","Events");
  PairMassDRSorted2C=HConfig.GetTH1D(Name+"_PairMassDRSorted2C","PairMassDRSorted2C",60,0.25,1.7,"M_{#mu#mu} (OS-SS, 2nd pair), GeV; cat C","Events");



  PairMassdRSorted=HConfig.GetTH2D(Name+"_PairMassdRSorted","PairMassdRSorted",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassVertexSorted=HConfig.GetTH2D(Name+"_PairMassVertexSorted","PairMassVertexSorted",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1VertexSorting=HConfig.GetTH1D(Name+"_PairMass1VertexSorting","PairMass1VertexSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st vertex sorting pair), GeV","Events");
  PairMass2VertexSorting=HConfig.GetTH1D(Name+"_PairMass2VertexSorting","PairMass2VertexSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd vertex sorting pair), GeV","Events");




  PairMassPhiMassSorting=HConfig.GetTH2D(Name+"_PairMassPhiMassSorting","PairMassPhiMassSorting",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1PhiMassSorting=HConfig.GetTH1D(Name+"_PairMass1PhiMassSorting","PairMass1PhiMassSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st #phi mass sorting pair), GeV","Events");
  PairMass2PhiMassSorting=HConfig.GetTH1D(Name+"_PairMass2PhiMassSorting","PairMass2PhiMassSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd #phi mass sorting pair), GeV","Events");



  PairMass1TauPhiMassSorting=HConfig.GetTH2D(Name+"_PairMass1TauPhiMassSorting","PairMass1TauPhiMassSorting",200,0.2,1.8,100,1.6,2.0,"M_{OS}, GeV","M_{#tau}, GeV");
  PairMass2TauPhiMassSorting=HConfig.GetTH2D(Name+"_PairMass2TauPhiMassSorting","PairMass2TauPhiMassSorting",200,0.2,1.8,100,1.6,2.0,"M_{OS}, GeV","M_{#tau}, GeV");


  PairMassdRSortedXVeto=HConfig.GetTH2D(Name+"_PairMassdRSortedXVeto","PairMassdRSortedXVeto",200,0.2,1.8,200,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");



  PairMassFinalSel=HConfig.GetTH2D(Name+"_PairMassFinalSel","PairMassFinalSel",60,0.2,1.8,60,0.2,1.8,"M_{OS}, GeV","M_{OS}, GeV");
  PairMass1=HConfig.GetTH1D(Name+"_PairMass1","PairMass1",40,0.2,1.777,"M_{1}, GeV","");
  PairMass2=HConfig.GetTH1D(Name+"_PairMass2","PairMass2",40,0.2,1.777,"M_{2}, GeV","");



  AllignSortMass1=HConfig.GetTH1D(Name+"_AllignSortMass1","AllignSortMass1",40,0.2,1.777,"M_{1} (#Delta R OS sorted), GeV","");
  AllignSortMass2=HConfig.GetTH1D(Name+"_AllignSortMass2","AllignSortMass2",40,0.2,1.777,"M_{2} (#Delta R OS sorted), GeV","");





  //  AllignSortMass1XVeto=HConfig.GetTH1D(Name+"_AllignSortMass1XVeto","AllignSortMass1XVeto",80,0.2,1.777,"M_{1} (#Delta R OS sorted), GeV","");
  //  AllignSortMass2XVeto=HConfig.GetTH1D(Name+"_AllignSortMass2XVeto","AllignSortMass2XVeto",80,0.2,1.777,"M_{2} (#Delta R OS sorted), GeV","");




  PairMassWithCut=HConfig.GetTH2D(Name+"_PairMassWithCut","PairMassWithCut",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEta=HConfig.GetTH2D(Name+"_PairMassEta","PairMassEta",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");
  PairMassEtaPrime=HConfig.GetTH2D(Name+"_PairMassEtaPrime","PairMassEtaPrime",40,0,2,40,0,2,"M_{OS}, GeV","M_{OS}, GeV");



  Muon1MVAID=HConfig.GetTH1D(Name+"_Muon1MVAID","Muon1MVAID",50,0.0,1.0,"#mu_{1} MVA","Events");
  Muon2MVAID=HConfig.GetTH1D(Name+"_Muon2MVAID","Muon2MVAID",50,0.0,1.0,"#mu_{2} MVA","Events");
  Muon3MVAID=HConfig.GetTH1D(Name+"_Muon3MVAID","Muon3MVAID",50,0.0,1.0,"#mu_{3} MVA","Events");



  BetterMuMuVertex=HConfig.GetTH1D(Name+"_BetterMuMuVertex","BetterMuMuVertex",30,0,5,"vertex pair quality (close)","");
  WorseMuMuVertex=HConfig.GetTH1D(Name+"_WorseMuMuVertex","WorseMuMuVertex",30,0,5,"vertex pair quality (far)","");

  OSSS1Angle_TRF=HConfig.GetTH1D(Name+"_OSSS1Angle_TRF","OSSS1Angle_TRF",50,-1.1,1.1,"cos#alpha_{1} (OS-SS1), 3#mu RF","");
  OSSS2Angle_TRF=HConfig.GetTH1D(Name+"_OSSS2Angle_TRF","OSSS2Angle_TRF",50,-1.1,1.1,"cos#alpha_{2} (OS-SS2), 3#mu RF","");


  OSSS1Angle_RRF=HConfig.GetTH1D(Name+"_OSSS1Angle_RRF","OSSS1Angle_RRF",50,-1.1,1.1,"cos#theta (n_{#rho} - n_{#mu})","");
  OSSS2Angle_RRF=HConfig.GetTH1D(Name+"_OSSS2Angle_RRF","OSSS2Angle_RRF",50,-1.1,1.1,"cos#theta (n_{#rho} - n_{#mu})","");

  cTheta_TRF_SSSS=HConfig.GetTH1D(Name+"_cTheta_TRF_SSSS","cTheta_TRF_SSSS",50,-1.1,1.1,"cos#theta (SS1-SS2), 3#mu RF","");
  cTheta_TRF_OSSS=HConfig.GetTH1D(Name+"_cTheta_TRF_OSSS","cTheta_TRF_OSSS",50,-1.1,1.1,"cos#theta (OS -SS1), 3#mu RF","");
  cTheta_MuonOS_TauPol_TRF=HConfig.GetTH1D(Name+"_cTheta_MuonOS_TauPol_TRF","cTheta_MuonOS_TauPol_TRF",50,-1.1,1.1,"cos#Theta (OS - #tau), 3#mu RF","");

  pTMu1OverMass_TRF=HConfig.GetTH1D(Name+"_pTMu1OverMass_TRF","pTMu1OverMass_TRF",50,0,0.6,"P_{#mu_{OS}}/M_{3#mu}  ","");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}



void  SignalSelector::Store_ExtraDist(){ 



  Extradist1d.push_back(&pTMu1OverMass_TRF);
  Extradist1d.push_back(&OSSS1Angle_TRF);
  Extradist1d.push_back(&OSSS2Angle_TRF);

  Extradist1d.push_back(&OSSS1Angle_RRF);
  Extradist1d.push_back(&OSSS2Angle_RRF);


  Extradist1d.push_back(&cTheta_MuonOS_TauPol_TRF);


  Extradist1d.push_back(&cTheta_TRF_SSSS);
  Extradist1d.push_back(&cTheta_TRF_OSSS);





  Extradist1d.push_back(&EtaGenSource);
  Extradist1d.push_back(&PairMassDRSorted1A);
  Extradist1d.push_back(&PairMassDRSorted2A);


  Extradist1d.push_back(&PairMassDRSorted1B);
  Extradist1d.push_back(&PairMassDRSorted2B);

  Extradist1d.push_back(&PairMassDRSorted1C);
  Extradist1d.push_back(&PairMassDRSorted2C);



  Extradist1d.push_back(&TauMassRefitABC1);  
  Extradist1d.push_back(&TauMassRefitABC2);

  Extradist1d.push_back(&TauMassRefitABC1_BDSeparateTrain);
  Extradist1d.push_back(&TauMassRefitABC2_BDSeparateTrain);


  Extradist2d.push_back(&TauMassRefitABC1_eta);  
  Extradist2d.push_back(&TauMassRefitABC2_eta);

  Extradist1d.push_back(&TauMassRefitA1);
  Extradist1d.push_back(&TauMassRefitB1);
  Extradist1d.push_back(&TauMassRefitC1);
  Extradist1d.push_back(&TauMassRefitA2);
  Extradist1d.push_back(&TauMassRefitB2);
  Extradist1d.push_back(&TauMassRefitC2);




  Extradist1d.push_back(&TauMassRefitABC1FullEtaVetoCut);
  Extradist1d.push_back(&TauMassRefitABC2FullEtaVetoCut);



  Extradist1d.push_back(&AllignSortMass1);
  Extradist1d.push_back(&AllignSortMass2);

  //  Extradist1d.push_back(&AllignSortMass1XVeto);
  //  Extradist1d.push_back(&AllignSortMass2XVeto);

  Extradist1d.push_back(&BetterMuMuVertex);
  Extradist1d.push_back(&WorseMuMuVertex);




  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);
  Extradist1d.push_back(&TauMassResolutionHelixRefit);

  Extradist2d.push_back(&TauMass_all_nophiVeto);
  Extradist1d.push_back(&TauMass_all);
  Extradist2d.push_back(&TauMass_allVsBDTA);
  Extradist2d.push_back(&TauMass_allVsBDTB);
  Extradist2d.push_back(&TauMass_allVsBDTC);

  Extradist2d.push_back(&EMR_tau_eta);

  Extradist1d.push_back(&VertexChi2KF);
  Extradist2d.push_back(&VertexChi2KF_vs_HelixFit);
  Extradist1d.push_back(&KF_Helix_deltaX);
  Extradist1d.push_back(&KF_Helix_deltaY);
  Extradist1d.push_back(&KF_Helix_deltaZ);

  Extradist1d.push_back(&SVPVTauDirAngle);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);


  Extradist1d.push_back(&TriggerMatchdR1);
  Extradist1d.push_back(&TriggerMatchdR2);
  Extradist1d.push_back(&TriggerMatchdR3);


  // Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&EventMassResolution_PtEtaPhi);
  Extradist1d.push_back(&EventMassResolution_PtEtaPhi_TauEta1p2);
  Extradist1d.push_back(&NSignalCandidates);


  Extradist1d.push_back(&BDTOutputA);
  Extradist1d.push_back(&BDTOutputB);
  Extradist1d.push_back(&BDTOutputC);
  Extradist1d.push_back(&BvsDBDTG);


  Extradist1d.push_back(&BDTOutput_train0_A);
  Extradist1d.push_back(&BDTOutput_train1_A);
  Extradist1d.push_back(&BDTOutput_train2_A);
  Extradist1d.push_back(&BDTOutput_train3_A);
  Extradist1d.push_back(&BDTOutput_train4_A);
  Extradist1d.push_back(&BDTOutput_train5_A);
  Extradist1d.push_back(&BDTOutput_train6_A);
  Extradist1d.push_back(&BDTOutput_train7_A);
  Extradist1d.push_back(&BDTOutput_train8_A);
  Extradist1d.push_back(&BDTOutput_train9_A);
  Extradist1d.push_back(&BDTOutput_train10_A);
  Extradist1d.push_back(&BDTOutput_train11_A);
  Extradist1d.push_back(&BDTOutput_train12_A);
  Extradist1d.push_back(&BDTOutput_train13_A);
  Extradist1d.push_back(&BDTOutput_train14_A);
  Extradist1d.push_back(&BDTOutput_train15_A);
  Extradist1d.push_back(&BDTOutput_train16_A);
  Extradist1d.push_back(&BDTOutput_train17_A);
  Extradist1d.push_back(&BDTOutput_train18_A);
  Extradist1d.push_back(&BDTOutput_train19_A);

  Extradist1d.push_back(&BDTOutput_train0_B);
  Extradist1d.push_back(&BDTOutput_train1_B);
  Extradist1d.push_back(&BDTOutput_train2_B);
  Extradist1d.push_back(&BDTOutput_train3_B);
  Extradist1d.push_back(&BDTOutput_train4_B);
  Extradist1d.push_back(&BDTOutput_train5_B);
  Extradist1d.push_back(&BDTOutput_train6_B);
  Extradist1d.push_back(&BDTOutput_train7_B);
  Extradist1d.push_back(&BDTOutput_train8_B);
  Extradist1d.push_back(&BDTOutput_train9_B);
  Extradist1d.push_back(&BDTOutput_train10_B);
  Extradist1d.push_back(&BDTOutput_train11_B);
  Extradist1d.push_back(&BDTOutput_train12_B);
  Extradist1d.push_back(&BDTOutput_train13_B);
  Extradist1d.push_back(&BDTOutput_train14_B);
  Extradist1d.push_back(&BDTOutput_train15_B);
  Extradist1d.push_back(&BDTOutput_train16_B);
  Extradist1d.push_back(&BDTOutput_train17_B);
  Extradist1d.push_back(&BDTOutput_train18_B);
  Extradist1d.push_back(&BDTOutput_train19_B);

  Extradist1d.push_back(&BDTOutput_train0_C);
  Extradist1d.push_back(&BDTOutput_train1_C);
  Extradist1d.push_back(&BDTOutput_train2_C);
  Extradist1d.push_back(&BDTOutput_train3_C);
  Extradist1d.push_back(&BDTOutput_train4_C);
  Extradist1d.push_back(&BDTOutput_train5_C);
  Extradist1d.push_back(&BDTOutput_train6_C);
  Extradist1d.push_back(&BDTOutput_train7_C);
  Extradist1d.push_back(&BDTOutput_train8_C);
  Extradist1d.push_back(&BDTOutput_train9_C);
  Extradist1d.push_back(&BDTOutput_train10_C);
  Extradist1d.push_back(&BDTOutput_train11_C);
  Extradist1d.push_back(&BDTOutput_train12_C);
  Extradist1d.push_back(&BDTOutput_train13_C);
  Extradist1d.push_back(&BDTOutput_train14_C);
  Extradist1d.push_back(&BDTOutput_train15_C);
  Extradist1d.push_back(&BDTOutput_train16_C);
  Extradist1d.push_back(&BDTOutput_train17_C);
  Extradist1d.push_back(&BDTOutput_train18_C);
  Extradist1d.push_back(&BDTOutput_train19_C);


  Extradist2d.push_back(&TauMass_vs_BDT_train0_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train1_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train2_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train3_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train4_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train5_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train6_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train7_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train8_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train9_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train10_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train11_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train12_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train13_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train14_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train15_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train16_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train17_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train18_A);
  Extradist2d.push_back(&TauMass_vs_BDT_train19_A);




  Extradist2d.push_back(&TauMass_vs_BDT_train0_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train1_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train2_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train3_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train4_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train5_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train6_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train7_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train8_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train9_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train10_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train11_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train12_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train13_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train14_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train15_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train16_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train17_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train18_B);
  Extradist2d.push_back(&TauMass_vs_BDT_train19_B);



  Extradist2d.push_back(&TauMass_vs_BDT_train0_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train1_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train2_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train3_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train4_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train5_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train6_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train7_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train8_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train9_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train10_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train11_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train12_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train13_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train14_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train15_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train16_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train17_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train18_C);
  Extradist2d.push_back(&TauMass_vs_BDT_train19_C);





  Extradist1d.push_back(&BvsDBDTG_ABC1);
  Extradist1d.push_back(&BvsDBDTG_ABC2);

  Extradist2d.push_back(&PairMass);


  Extradist1d.push_back(&KpiIsolationMass_OS);
  Extradist1d.push_back(&KpiIsolationMass_SS1);
  Extradist1d.push_back(&KpiIsolationMass_SS2);


  Extradist2d.push_back(&PairMassdRSorted);
  Extradist2d.push_back(&PairMassVertexSorted);

  Extradist1d.push_back(&PairMass1VertexSorting);
  Extradist1d.push_back(&PairMass2VertexSorting);


  Extradist2d.push_back(&PairMassdRSortedXVeto);
  Extradist2d.push_back(&PairMassPhiMassSorting);

  Extradist2d.push_back(&PairMassFinalSel);
  Extradist1d.push_back(&PairMass1);
  Extradist1d.push_back(&PairMass2);
  Extradist2d.push_back(&PairMassWithCut);
  Extradist2d.push_back(&PairMassEta);
  Extradist2d.push_back(&PairMassEtaPrime);



  Extradist1d.push_back(&Muon1MVAID);
  Extradist1d.push_back(&Muon2MVAID);
  Extradist1d.push_back(&Muon3MVAID);


  Extradist1d.push_back(&PairMass1NoSorting);
  Extradist1d.push_back(&PairMass2NoSorting);
  Extradist2d.push_back(&MuMuMassNoSorting);

  Extradist1d.push_back(&PairMass1PTSorting);
  Extradist1d.push_back(&PairMass2PTSorting);
  Extradist2d.push_back(&MuMuMassPTSorting);


  Extradist1d.push_back(&PairMass1AllignedSorting);
  Extradist1d.push_back(&PairMass2AllignedSorting);
  Extradist2d.push_back(&MuMuMassAllignedSorting);


  Extradist1d.push_back(&PairMass1PhiMassSorting);
  Extradist1d.push_back(&PairMass2PhiMassSorting);

  Extradist2d.push_back(&PairMass1TauPhiMassSorting);
  Extradist2d.push_back(&PairMass2TauPhiMassSorting);


  Extradist1d.push_back(&TauMassRefitABC1_train11);
  Extradist1d.push_back(&TauMassRefitABC2_train11);
  Extradist1d.push_back(&TauMassRefitA1_train11);
  Extradist1d.push_back(&TauMassRefitB1_train11);
  Extradist1d.push_back(&TauMassRefitC1_train11);
  Extradist1d.push_back(&TauMassRefitA2_train11);
  Extradist1d.push_back(&TauMassRefitB2_train11);
  Extradist1d.push_back(&TauMassRefitC2_train11);
  Extradist1d.push_back(&AllignSortMass1_train11);
  Extradist1d.push_back(&AllignSortMass2_train11);
  
  Extradist1d.push_back(&BDTOutputA_train11);
  Extradist1d.push_back(&BDTOutputB_train11);
  Extradist1d.push_back(&BDTOutputC_train11);

  Extradist1d.push_back(&TauMassRefitABC1_train8);
  Extradist1d.push_back(&TauMassRefitABC2_train8);
  Extradist1d.push_back(&TauMassRefitA1_train8);
  Extradist1d.push_back(&TauMassRefitB1_train8);
  Extradist1d.push_back(&TauMassRefitC1_train8);
  Extradist1d.push_back(&TauMassRefitA2_train8);
  Extradist1d.push_back(&TauMassRefitB2_train8);
  Extradist1d.push_back(&TauMassRefitC2_train8);
  Extradist1d.push_back(&AllignSortMass1_train8);
  Extradist1d.push_back(&AllignSortMass2_train8);
  
  Extradist1d.push_back(&BDTOutputA_train8);
  Extradist1d.push_back(&BDTOutputB_train8);
  Extradist1d.push_back(&BDTOutputC_train8);







}


void  SignalSelector::doEvent(){ 

  
  unsigned int t;
  int id(Ntp->GetMCID());
  //  std::cout<<" id   "<< id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}


  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu0Fired(false);
  bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);

  random_num = rndm.Rndm();


  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //    std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("DoubleMu3_TkMu_DsTau3Mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}
  }

  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //    std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;

  if (HLTOk) value.at(HLT) = true;
  else value.at(HLT) = false;

  if (L1Ok) value.at(L1T) = true;
  else value.at(L1T) = false;


  //  if(DoubleMuFired) value.at(L1T)=1;

  //  std::cout<<"  "<< value.at(L1T) << "  "<<value.at(HLT)  << std::endl;


  pass.at(L1T)= (value.at(L1T)==cut.at(L1T));
  pass.at(HLT)= (value.at(HLT)==cut.at(HLT));



  value.at(SignalCandidate)=0;
  unsigned int  signal_idx=0;
  value.at(TriggerMatch)=0;

  double min_chi2(99.);
  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }


  NSignalCandidates.at(t).Fill(Ntp->NThreeMuons(),1);
  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();

    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);

    unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);
    //*** Muon ID cut
    value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) &&  Ntp->Muon_isPFMuon(mu1_pt_idx) &&
    			Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&  Ntp->Muon_isPFMuon(mu2_pt_idx) &&
			Ntp->Muon_isGlobalMuon(mu3_pt_idx) &&  Ntp->Muon_isPFMuon(mu3_pt_idx));



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
    TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);


    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();




    value.at(PhiVeto1) =  M_osss1;//fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(PhiVeto2) =  M_osss2;//fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 
    value.at(OmegaVeto1) = M_osss1;//fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;
    value.at(OmegaVeto2) = M_osss2;//fabs(M_osss1-PDG_Var::Omega_mass())< fabs(M_osss2-PDG_Var::Omega_mass()) ? M_osss1 : M_osss2;


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
    if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
    value.at(TriggerMatch) = triggerCheck;

    value.at(TauMassCut) = TauRefittedLV.M();
  }

  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut));
  pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut));
  pass.at(Mu3PtCut) = (value.at(Mu3PtCut) >= cut.at(Mu3PtCut));
  pass.at(MuonID)   =(value.at(MuonID)  == cut.at(MuonID));
  pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));
  pass.at(PhiVeto1) = true;//(value.at(PhiVeto1) < 0.98 || value.at(PhiVeto1) > 1.06 );
  pass.at(OmegaVeto1) = true;//(value.at(OmegaVeto1) < 0.742 || value.at(OmegaVeto1) > 0.822 );
  pass.at(PhiVeto2) = true;//(value.at(PhiVeto2) < 0.98 || value.at(PhiVeto2) > 1.06 );
  pass.at(OmegaVeto2) = true;//(value.at(OmegaVeto2) < 0.742 || value.at(OmegaVeto2) > 0.822 );
  pass.at(TauMassCut) =( (value.at(TauMassCut) > tauMinSideBand_)  &&   (value.at(TauMassCut) < tauMaxSideBand_ ));

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(PhiVeto1);
  exclude_cuts.push_back(OmegaVeto1);


  exclude_cuts.push_back(PhiVeto2);
  exclude_cuts.push_back(OmegaVeto2);

  if(passAllBut(exclude_cuts)){


    TLorentzVector TauRefittedLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);

    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
    vector<unsigned int> idx_vec;
    
    idx_vec.push_back(mu1_idx);
    idx_vec.push_back(mu2_idx);
    idx_vec.push_back(mu3_idx);

    unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
    double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();

    double pmass  = fabs(M_osss1-PDG_Var::Phi_mass())  < fabs(M_osss2-PDG_Var::Phi_mass()) ? M_osss1 : M_osss2; 

    TauMass_all_nophiVeto.at(t).Fill(TauRefittedLV.M(),pmass,1);
  }


  double wobs=1;
  double w; 

  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}


    if(id == 128){
      TLorentzVector EtaMu1(0,0,0,0),EtaMu2(0,0,0,0),Eta(0,0,0,0),EtaGamma(0,0,0,0);
      int EtaGenIndex(-1);
      for(unsigned int imc = 0; imc < Ntp->NMCParticles(); imc++){
	
	if(Ntp->MCParticle_pdgid(imc) == 331){//   or Ntp->MCParticle_pdgid(imc) == 333  or Ntp->MCParticle_pdgid(imc) == 331){

	  Eta = Ntp->MCParticle_p4(imc);
	  for(unsigned int imc_childs=0; imc_childs < Ntp->MCParticle_childpdgid(imc).size(); imc_childs++)
	    {
	      
	      if( abs(Ntp->MCParticle_childpdgid(imc).at(imc_childs)) == 13)
		{
		  EtaMu1 = Ntp->MCParticle_p4(Ntp->MCParticle_childidx(imc).at(0));
		  EtaMu2 = Ntp->MCParticle_p4(Ntp->MCParticle_childidx(imc).at(1));
		  EtaGenIndex = imc;	  
		}
	    }
	}
      }
    

    
      //      std::cout<<"  mother index  "<< Ntp->MCParticle_midx(EtaGenIndex) <<"   "  <<Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(EtaGenIndex)) <<std::endl;

      if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(EtaGenIndex)))== 411) EtaGenSource.at(t).Fill(0.,w);
      if(abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(EtaGenIndex)))== 431) EtaGenSource.at(t).Fill(1.,w);
    }


  bool status=AnalysisCuts(t,w,wobs);


  //  if(id!=1)  std::cout<<" id:   "<< id << "  NMCSignalParticles  "<< Ntp->NMCSignalParticles() << "  NMCTaus   "<< Ntp->NMCTaus() << std::endl;


  if(status){


    var_Vertex2muTrkKF = Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(signal_idx);


    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);



    std::vector<unsigned int> EtaSortedIndices;
    
    EtaSortedIndices.push_back(Muon_Eta_index_1);
    EtaSortedIndices.push_back(Muon_Eta_index_2);
    EtaSortedIndices.push_back(Muon_Eta_index_3);

    EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);



    //*** Rapidity sorted muons
    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);  
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);


    if(fabs(   (Muon1LV + Muon2LV + Muon3LV).Eta() ) < 1.2 )    EventMassResolution_PtEtaPhi_TauEta1p2.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);



    //    std::cout <<"  " << Muon1LV.M() <<"   " <<Muon2LV.M() <<"   "  <<  Muon3LV.M() <<std::endl;

    vector<unsigned int> idx_vec;
    idx_vec.push_back(Muon_index_1);
    idx_vec.push_back(Muon_index_2);
    idx_vec.push_back(Muon_index_3);

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    //*** With such sorting pT(ss1) > pT(ss2)
    TLorentzVector MuonOS  = Ntp->Muon_P4(os_mu_idx);  
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);

    TLorentzVector MuonOS_Rotated = MuonOS; MuonOS_Rotated.SetVect(Ntp->Rotate(MuonOS_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
    TLorentzVector MuonSS1_Rotated = MuonSS1; MuonSS1_Rotated.SetVect(Ntp->Rotate(MuonSS1_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
    TLorentzVector MuonSS2_Rotated = MuonSS2; MuonSS2_Rotated.SetVect(Ntp->Rotate(MuonSS2_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
    TLorentzVector ThreeMuon_Rotated = MuonOS + MuonSS1 +MuonSS2; ThreeMuon_Rotated.SetVect(Ntp->Rotate(ThreeMuon_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));


    TLorentzVector MuonOS_TRF  = Ntp->Boost(MuonOS_Rotated, ThreeMuon_Rotated);    // Muon in Tau REst Frame with Z axis alligned on tau directions
    TLorentzVector MuonSS1_TRF = Ntp->Boost(MuonSS1_Rotated, ThreeMuon_Rotated);
    TLorentzVector MuonSS2_TRF = Ntp->Boost(MuonSS2_Rotated, ThreeMuon_Rotated);



    TLorentzVector MuonOS_R1_RF_rotated = MuonOS_TRF; MuonOS_R1_RF_rotated.SetVect(Ntp->Rotate(MuonOS_R1_RF_rotated.Vect(), (MuonOS_TRF + MuonSS1_TRF).Vect())   );
    TLorentzVector MuonSS1_R1_RF_rotated = MuonSS1_TRF; MuonSS1_R1_RF_rotated.SetVect(Ntp->Rotate(MuonSS1_R1_RF_rotated.Vect(), (MuonOS_TRF + MuonSS1_TRF).Vect())   );



    TLorentzVector MuonOS_R2_RF_rotated = MuonOS_TRF; MuonOS_R2_RF_rotated.SetVect(Ntp->Rotate(MuonOS_R2_RF_rotated.Vect(), (MuonOS_TRF + MuonSS2_TRF).Vect())   );
    TLorentzVector MuonSS2_R2_RF_rotated = MuonSS2_TRF; MuonSS2_R2_RF_rotated.SetVect(Ntp->Rotate(MuonSS2_R2_RF_rotated.Vect(), (MuonOS_TRF + MuonSS2_TRF).Vect()  )   );

    TLorentzVector MuonOS_R1_RF = Ntp->Boost(MuonOS_R1_RF_rotated, MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated);
    TLorentzVector MuonSS1_R1_RF=Ntp->Boost(MuonSS1_R1_RF_rotated, MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated);

    TLorentzVector MuonOS_R2_RF = Ntp->Boost(MuonOS_R2_RF_rotated, MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated);
    TLorentzVector MuonSS2_R2_RF=Ntp->Boost(MuonSS2_R2_RF_rotated, MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated);

    double R1_RF_Theta =     MuonOS_R1_RF.Vect().Dot( (MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated).Vect() ) * (1./MuonOS_R1_RF.Vect().Mag()/(MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated).Vect().Mag());
    double R2_RF_Theta =     MuonOS_R2_RF.Vect().Dot( (MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated).Vect() ) * (1./MuonOS_R2_RF.Vect().Mag()/(MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated).Vect().Mag());

    OSSS1Angle_RRF.at(t).Fill(R1_RF_Theta,1.);
    OSSS2Angle_RRF.at(t).Fill(R2_RF_Theta,1.);


    OSSS1Angle_TRF.at(t).Fill(( MuonOS_TRF.Vect() * MuonSS1_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS1_TRF.Vect().Mag()),1.);
    OSSS2Angle_TRF.at(t).Fill(( MuonOS_TRF.Vect() * MuonSS2_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS2_TRF.Vect().Mag()),1.);

    double costheta_TRF_SSSS;
    double costheta_TRF_OSSS;

    TVector3 n_tau(0,0,1);




    //    (Ntp->Rotate(MuonOS.Vect(),MuonOS.Vect())).Print();
      
    if(MuonSS1_TRF.Vect().Mag() > MuonSS2_TRF.Vect().Mag()){

      costheta_TRF_SSSS = MuonSS1_TRF.Vect().Cross(MuonSS2_TRF.Vect()).Dot(n_tau) *(1/ (MuonSS1_TRF.Vect().Cross(MuonSS2_TRF.Vect())).Mag());
      costheta_TRF_OSSS = MuonOS_TRF.Vect().Cross(MuonSS1_TRF.Vect()).Dot(n_tau) *(1/ ( MuonOS_TRF.Vect().Cross(MuonSS1_TRF.Vect())).Mag());

    } else {
      
      costheta_TRF_SSSS = MuonSS2_TRF.Vect().Cross(MuonSS1_TRF.Vect()).Dot(n_tau) *(1/ (MuonSS2_TRF.Vect().Cross(MuonSS1_TRF.Vect())).Mag());
      costheta_TRF_OSSS = MuonOS_TRF.Vect().Cross(MuonSS2_TRF.Vect()).Dot(n_tau) *(1/  (MuonOS_TRF.Vect().Cross(MuonSS2_TRF.Vect())).Mag());

    }
      

    cTheta_TRF_SSSS.at(t).Fill(costheta_TRF_SSSS,1.);
    cTheta_TRF_OSSS.at(t).Fill(costheta_TRF_OSSS,1.);

    cTheta_MuonOS_TauPol_TRF.at(t).Fill(MuonOS_TRF.Vect().Dot(n_tau)*(1./MuonOS_TRF.Vect().Mag()));
    pTMu1OverMass_TRF.at(t).Fill(MuonOS_TRF.P()/(MuonOS_TRF + MuonSS1_TRF + MuonSS2_TRF).M(),1.);



    //    std::cout<<" angle OS-SS1, OS - SS2  "<<
      // ( MuonOS_TRF.Vect() * MuonSS1_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS1_TRF.Vect().Mag())<< "   " <<
      // ( MuonOS_TRF.Vect() * MuonSS2_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS2_TRF.Vect().Mag())<< std::endl;

    double VertexQuality_OS_SS1;    
    double VertexQuality_OS_SS2;    


    //    std::cout<<"lets check "<< std::endl;
    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,0);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,2);
	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,2);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,0);
	  }
      }


    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,0);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,1);
	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,1);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,0);
	  }
      }


    if(MuonOS.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2))) ==0 )
      {
	if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)))==0)
	  {
	    VertexQuality_OS_SS1 = Ntp->Vertex_pair_quality(signal_idx,2);
	    VertexQuality_OS_SS2 = Ntp->Vertex_pair_quality(signal_idx,1);

	  }
	else if(MuonSS1.DeltaR(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)))==0)
	  {
	    VertexQuality_OS_SS1  = Ntp->Vertex_pair_quality(signal_idx,1);
	    VertexQuality_OS_SS2  = Ntp->Vertex_pair_quality(signal_idx,2);
	  }
      }


    PairMass.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1); 


    float VertexQualitySortedMass1,VertexQualitySortedMass2;
    float BetterPhiVertex, WorsePhiVertex;
    if(VertexQuality_OS_SS1 < VertexQuality_OS_SS2){
      VertexQualitySortedMass1 = (MuonOS+MuonSS1).M();
      VertexQualitySortedMass2 = (MuonOS+MuonSS2).M();
      BetterPhiVertex=VertexQuality_OS_SS1;
      WorsePhiVertex=VertexQuality_OS_SS2;
    }else{
      VertexQualitySortedMass1 = (MuonOS+MuonSS2).M();
      VertexQualitySortedMass2 = (MuonOS+MuonSS1).M();
      BetterPhiVertex=VertexQuality_OS_SS2;
      WorsePhiVertex=VertexQuality_OS_SS1;
    }

    BetterMuMuVertex.at(t).Fill(BetterPhiVertex,1);
    WorseMuMuVertex.at(t).Fill(WorsePhiVertex,1);
    PairMassVertexSorted.at(t).Fill(VertexQualitySortedMass1,VertexQualitySortedMass2 ,1);
    PairMass1VertexSorting.at(t).Fill(VertexQualitySortedMass1,w);
    PairMass2VertexSorting.at(t).Fill(VertexQualitySortedMass2,w);



    std::vector<unsigned int> Indices;
    Indices.push_back(ss1_mu_idx);
    Indices.push_back(ss2_mu_idx);



    unsigned int SS1RandomIndex(0);
    unsigned int SS2RandomIndex(0);


    float random_muon_index = rndm.Uniform();
    if(random_muon_index >= 0.5 ){SS1RandomIndex =  Indices.at(0); SS2RandomIndex = Indices.at(1) ; }
    if(random_muon_index <  0.5 ){SS1RandomIndex =  Indices.at(1); SS2RandomIndex = Indices.at(0) ; }

    TLorentzVector MuonLV_RandomSS1 = Ntp->Muon_P4(SS1RandomIndex);
    TLorentzVector MuonLV_RandomSS2 = Ntp->Muon_P4(SS2RandomIndex);
    
    
    PairMass1NoSorting.at(t).Fill((MuonOS+MuonLV_RandomSS1).M(),w);
    PairMass2NoSorting.at(t).Fill((MuonOS+MuonLV_RandomSS2).M(),w);
    MuMuMassNoSorting.at(t).Fill((MuonOS+MuonLV_RandomSS2).M(),(MuonOS+MuonLV_RandomSS1).M());


    PairMass1PTSorting.at(t).Fill((MuonOS+MuonSS1).M(),w);
    PairMass2PTSorting.at(t).Fill((MuonOS+MuonSS2).M(),w);
    MuMuMassPTSorting.at(t).Fill((MuonOS+MuonSS1).M(),(MuonOS+MuonSS2).M());


    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      PairMass1AllignedSorting.at(t).Fill((MuonOS+MuonSS2).M(),w);
      PairMass2AllignedSorting.at(t).Fill((MuonOS+MuonSS1).M(),w);
      MuMuMassAllignedSorting.at(t).Fill((MuonOS+MuonSS1).M(),(MuonOS+MuonSS2).M());
    }else{
      PairMass1AllignedSorting.at(t).Fill((MuonOS+MuonSS1).M(),w);
      PairMass2AllignedSorting.at(t).Fill((MuonOS+MuonSS2).M(),w);
      MuMuMassAllignedSorting.at(t).Fill((MuonOS+MuonSS2).M(),(MuonOS+MuonSS1).M());
    }
  


      //      (EtaMu1+EtaMu2).Print();
            if(id == 120){
	      std::cout<<"-------------- All categoris -------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
	      std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	      std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	      std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	      Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
	      Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
	      Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
	      Ntp->printMCDecayChainOfEvent(true, true, true, true);
	      std::cout<< "\n\n\n\n\n\n";
	    }
      
  
    float dRSortedMassPair1,dRSortedMassPair2;
    unsigned int ss1_mu_idx_dr, ss2_mu_idx_dr;

    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      ss1_mu_idx_dr=ss2_mu_idx;
      ss2_mu_idx_dr=ss1_mu_idx;
      dRSortedMassPair1 = (MuonOS+MuonSS2).M();
      dRSortedMassPair2 = (MuonOS+MuonSS1).M();
    }else{
      ss1_mu_idx_dr=ss1_mu_idx;
      ss2_mu_idx_dr=ss2_mu_idx;
      dRSortedMassPair1 = (MuonOS+MuonSS1).M();
      dRSortedMassPair2 = (MuonOS+MuonSS2).M();
    }

    TLorentzVector KaonOS(0,0,0,0);  KaonOS.SetXYZM(MuonOS.Px(),MuonOS.Py(),MuonOS.Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS1(0,0,0,0); KaonSS1.SetXYZM(Ntp->Muon_P4(ss1_mu_idx_dr).Px(),Ntp->Muon_P4(ss1_mu_idx_dr).Py(),Ntp->Muon_P4(ss1_mu_idx_dr).Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS2(0,0,0,0); KaonSS2.SetXYZM(Ntp->Muon_P4(ss2_mu_idx_dr).Px(),Ntp->Muon_P4(ss2_mu_idx_dr).Py(),Ntp->Muon_P4(ss2_mu_idx_dr).Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS1pTsort(0,0,0,0); KaonSS1pTsort.SetXYZM(MuonSS1.Px(),MuonSS1.Py(),MuonSS1.Pz(),PDG_Var::Kp_mass());
    TLorentzVector KaonSS2pTsort(0,0,0,0); KaonSS2pTsort.SetXYZM(MuonSS2.Px(),MuonSS2.Py(),MuonSS2.Pz(),PDG_Var::Kp_mass());

    







    for(unsigned int iIsoTrack=0; iIsoTrack < Ntp->NIsolationTrack(signal_idx); iIsoTrack++){


      if(Ntp->Muon_charge(os_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_OS.at(t).Fill((KaonOS + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
      if(Ntp->Muon_charge(ss1_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_SS1.at(t).Fill((KaonSS1pTsort + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
      if(Ntp->Muon_charge(ss2_mu_idx)*Ntp->IsolationTrack_charge(signal_idx,iIsoTrack)==-1)
	{
	  KpiIsolationMass_SS2.at(t).Fill((KaonSS2pTsort + Ntp->IsolationTrack_p4(signal_idx,iIsoTrack)).M(),1);
	}
    }



    
    PairMassdRSorted.at(t).Fill(dRSortedMassPair2,dRSortedMassPair1 ,1); 
    
    
    float Mass_osss1 = (Ntp->Muon_P4(os_mu_idx)+Ntp->Muon_P4(ss1_mu_idx)).M();
    float Mass_osss2 = (Ntp->Muon_P4(os_mu_idx)+Ntp->Muon_P4(ss2_mu_idx)).M();
    
    float CloserToPhiMassPair  = fabs(Mass_osss1-PDG_Var::Phi_mass())  < fabs(Mass_osss2-PDG_Var::Phi_mass()) ? Mass_osss1 : Mass_osss2;
    
    if(CloserToPhiMassPair==Mass_osss1)
      {

	PairMass1PhiMassSorting.at(t).Fill(Mass_osss1,1);
	PairMass2PhiMassSorting.at(t).Fill(Mass_osss2,1);

	PairMass1TauPhiMassSorting.at(t).Fill(Mass_osss1,(MuonOS+MuonSS1+MuonSS2).M());
	PairMass2TauPhiMassSorting.at(t).Fill(Mass_osss2,(MuonOS+MuonSS1+MuonSS2).M());
 
	PairMassPhiMassSorting.at(t).Fill(Mass_osss1,Mass_osss2);
 

      }
    if(CloserToPhiMassPair==Mass_osss2)
      {

	PairMass1PhiMassSorting.at(t).Fill(Mass_osss2,1);
	PairMass2PhiMassSorting.at(t).Fill(Mass_osss1,1);

	PairMass1TauPhiMassSorting.at(t).Fill(Mass_osss2,(MuonOS+MuonSS1+MuonSS2).M());
	PairMass2TauPhiMassSorting.at(t).Fill(Mass_osss1,(MuonOS+MuonSS1+MuonSS2).M());
 

	PairMassPhiMassSorting.at(t).Fill(Mass_osss2,Mass_osss1);

      }
    ////////////////////////////////////////////////////////  var_Iso08MuMin







    //////////////////////////////////////////////////////////////////////////
    

    bool RemoveEta(false);
    bool RemoveHalfEta(false);
    
    
    bool phiVeto(false);
    bool rmgVeto(false);
    
    bool CrossVeto1(true);
    bool CrossVeto2(true);
    bool CrossVeto3(true);
    bool CrossVeto(true);
    double    m12v = (MuonOS+MuonSS1).M();
    double    m13v = (MuonOS+MuonSS2).M();
    
    
    //    if(( m12v < phiVetoCut1  || m12v > phiVetoCut2 )  && (m13v < phiVetoCut1 || m13v > phiVetoCut2)  )  phiVeto=true;
    if(( m12v < rmgCutVeto1 || m12v > rmgCutVeto2 )  && (m13v < rmgCutVeto1 || m13v > rmgCutVeto2))  rmgVeto=true;
    
    
    //    if((dRSortedMassPair1 < phiVetoCut1 || dRSortedMassPair1 > phiVetoCut2 ) && (dRSortedMassPair2 < 0.65 || dRSortedMassPair2 > 1.6) )CrossVeto=true;
    if((dRSortedMassPair1 > phiVetoCut1 &&  dRSortedMassPair2 > 0.65)  &&  (dRSortedMassPair1 < phiVetoCut2 &&  dRSortedMassPair2 < 1.6))CrossVeto1=false;
    if((dRSortedMassPair2 > phiVetoCut1 &&  dRSortedMassPair1 > 0.2)  &&  (dRSortedMassPair2 < phiVetoCut2 &&  dRSortedMassPair1 < 1.4))CrossVeto2=false;
    if((dRSortedMassPair1 > rmgCutVeto1 &&  dRSortedMassPair2 > 0.95)  &&  (dRSortedMassPair1 < rmgCutVeto2 &&  dRSortedMassPair2 < 1.5))CrossVeto3=false;


    CrossVeto = CrossVeto1*CrossVeto2*CrossVeto3;
    if(CrossVeto)   {PairMassdRSortedXVeto.at(t).Fill(dRSortedMassPair2,dRSortedMassPair1 ,1);
    }


    if((MuonOS+MuonSS1).M() > 0.549 && (MuonOS+MuonSS2).M() > 0.549) RemoveEta = true;
    if((MuonOS+MuonSS2).M() > 0.549) RemoveHalfEta = true;
    if(RemoveEta && phiVeto && rmgVeto)    PairMassWithCut.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);

    

    //*** 3-mu mass after the KF vertex constrain

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) + 
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);




    //    std::cout<<" M1   "<<  Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0).M() <<" M2 "  << Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1).M() <<"  M3   "  <<Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2).M() <<std::endl;
    //*** uncomment if you want to have print outs
    /*
    if(id ==120 ){// or id == 40){
      std::cout<<"-------------- All categoris ----------------"<< std::endl;
      std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Muon1LV.Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Muon1LV) << std::endl;
      Muon2LV.Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Muon2LV) << std::endl;
      Muon3LV.Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Muon3LV) << std::endl;
      Ntp->printMCDecayChainOfEvent(true, true, true, true);
      std::cout<< "\n\n\n\n\n\n";
    }
    */



    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);

    EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());  // Event Mass resolution


    TriggerMatchdR1.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(0),1);
    TriggerMatchdR2.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(1),1);
    TriggerMatchdR3.at(t).Fill(Ntp->ThreeMuons_TriggerMatch_dR(signal_idx).at(2),1);

    VertexChi2KF.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(signal_idx),w);
    FLSignificance.at(t).Fill(sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(signal_idx),Ntp->Vertex_PrimaryVertex_Covariance(signal_idx),
								   Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx))),w);
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    SVPVTauDirAngle.at(t).Fill(SVPV.Angle(TauLV.Vect()),w);

    
    //***  define the mva varables used for evaluation of BDT weights for selection

    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      var_mass12_dRsorting = (MuonOS+MuonSS2).M();
      var_mass13_drSorting = (MuonOS+MuonSS1).M();
    }else{
      var_mass12_dRsorting = (MuonOS+MuonSS1).M();
      var_mass13_drSorting = (MuonOS+MuonSS2).M();
    }
    
    var_vertexKFChi2 =Ntp->Vertex_signal_KF_Chi2(signal_idx);
    var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
    var_flightLenSig = sqrt( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(signal_idx),Ntp->Vertex_PrimaryVertex_Covariance(signal_idx),
							    Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_Signal_KF_Covariance(signal_idx)));
    var_flightLenDist =  (Ntp->Vertex_MatchedPrimaryVertex(signal_idx) - Ntp->Vertex_Signal_KF_pos(signal_idx)).Mag();
    var_sumMuTrkKinkChi2= (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));
    var_segCompMuMin  = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
    var_segCompMuMax  = std::max({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
    var_MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});


    var_MuMu_mindR = std::min({Muon1LV.DeltaR(Muon2LV),Muon2LV.DeltaR(Muon3LV),Muon1LV.DeltaR(Muon3LV)});
    var_RelPt_Mu1Tau = Ntp->Muon_P4(Muon_index_1).Pt()/TauLV.Pt();
    var_Eta_au = TauLV.Eta();
    var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(signal_idx,0), Ntp->Vertex_pair_quality(signal_idx,1), Ntp->Vertex_pair_quality(signal_idx,2)});
    var_maxdca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});
    var_MuTau_maxdR = std::max({Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV),Muon1LV.DeltaR(TauLV)});

    var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});

    var_MaxMuon_chi2LocalPosition = std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });
    var_MaxtrkKink = std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});

    var_MaxD0SigSV=    std::max({Ntp->Vertex_d0sigSV_reco(signal_idx,0),
	  Ntp->Vertex_d0sigSV_reco(signal_idx,1),
	  Ntp->Vertex_d0sigSV_reco(signal_idx,2)});
  


    var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(signal_idx),Ntp->Vertex_DCA23(signal_idx),Ntp->Vertex_DCA31(signal_idx)});



    var_MaxVertexPairQuality =   std::max({Ntp->Vertex_pair_quality(signal_idx,0),Ntp->Vertex_pair_quality(signal_idx,1),Ntp->Vertex_pair_quality(signal_idx,2)});
    var_MaxMuon_chi2LocalMomentum = std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });


    var_MuonglbkinkSum    = (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));



    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(signal_idx,0),
	  Ntp->Vertex_d0sig_reco(signal_idx,1),
	  Ntp->Vertex_d0sig_reco(signal_idx,2)});

    var_MaxD0Significance = MaxD0Significance;
    var_IsolationMinDist = Ntp->Isolation_MinDist(signal_idx);



    float MaxD0BSSignificance = std::max({Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,0),
	  Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,1),
	  Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,2)});
    
    
    float MinD0BSSignificance = std::min({Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,0),
	  Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,1),
	  Ntp->Vertex_d0BeamSpot_reco_sig(signal_idx,2)});
    

    var_MaxD0SigBS = MaxD0BSSignificance;
    var_MinD0SigBS = MinD0BSSignificance;




    TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_1),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_2),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
    TVector3 Mu3ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_3),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));


    
    var_MinMuonImpactAngle = std::min({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())});
    var_MaxMuonImpactAngle = std::max({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())});
    //  ----------------------------- secondary vertices ----------------
    int NumberOfPrimaryVertices(0);
    for(unsigned int iVertex=0; iVertex < Ntp->NSecondaryVertices(); iVertex++){
      TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(signal_idx),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
      TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(iVertex),Ntp->Vertex_MatchedPrimaryVertex(signal_idx));
        
        
      if(SVfakePV.DeltaR(SVsignalPV) < 0.3 && (Ntp->Vertex_Signal_KF_pos(signal_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag() > 0.015){ // sv in the tau cone but  displaced
	NumberOfPrimaryVertices++;
	    
      }
    }
    var_nsv = NumberOfPrimaryVertices;


    // -----------------------------------------------------------------



    var_tauMass=TauRefitLV.M();
    TauMass_all.at(t).Fill(TauRefitLV.M(),1);

    // ** Isolation mass spectra
    double Chi2IsoTrackVertexToMuon3(999.);
    TLorentzVector IsoTrack_P4_closestToMu3(0,0,0,0);

    double Chi2IsoTrackVertexToMuon2(999.);
    TLorentzVector IsoTrack_P4_closestToMu2(0,0,0,0);

    double Chi2IsoTrackVertexToMuon1(999.);
    TLorentzVector IsoTrack_P4_closestToMu1(0,0,0,0);
    int NcloseTracksCount(0);
    int TrackIndex_closestToPV(0);
    int TrackIndex(0);
    double dca_temp(999.);
    double dcaPV_temp(999.);
    double TrackPTtreschold(0.4);
    double SumPT08;
    for(int i =0; i< Ntp->NIsolationTrack(signal_idx); i++){

      if(Ntp->IsolationTrack_p4(signal_idx,i).DeltaR(TauRefitLV) < 0.8){
	SumPT08 += Ntp->IsolationTrack_p4(signal_idx,i).Pt();
      }
        
      if(Ntp->IsolationTrack_p4(signal_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(signal_idx,i),2)   +   
								  pow(Ntp->IsolationTrack_dxySV(signal_idx,i),2)) < 0.03)
	{
	  NcloseTracksCount++;
	}


      if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(2)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon3)
	       {
		 Chi2IsoTrackVertexToMuon3= Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu3 = Ntp->IsolationTrack_p4(signal_idx,i);
	       }
	   }


	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(1)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon2)
	       {
		 Chi2IsoTrackVertexToMuon2= Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu2 = Ntp->IsolationTrack_p4(signal_idx,i);
	       }
	   }

	   
	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(signal_idx).at(0)) * Ntp->IsolationTrack_charge(signal_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,i) < Chi2IsoTrackVertexToMuon1)
	       {
		 Chi2IsoTrackVertexToMuon1= Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(signal_idx,i);
		 IsoTrack_P4_closestToMu1 = Ntp->IsolationTrack_p4(signal_idx,i);
	       }
	   }



	   if( sqrt(  pow(Ntp->IsolationTrack_dzSV(signal_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(signal_idx,i),2) ) <  dca_temp){
	     dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(signal_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(signal_idx,i),2));
	     TrackIndex = i;
	   }

	   if(Ntp->IsolationTrack_p4(signal_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(signal_idx,i)) < 0.05 && 
	      sqrt(  pow(Ntp->IsolationTrack_dzSV(signal_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(signal_idx,i),2)) < 0.05){
	     if( sqrt(  pow(Ntp->IsolationTrack_dzPV(signal_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(signal_idx,i),2) ) <  dcaPV_temp){
	       dcaPV_temp = sqrt(  pow(Ntp->IsolationTrack_dzPV(signal_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(signal_idx,i),2));
	       TrackIndex_closestToPV = i;
	     }
	   }
    }

    var_Iso08 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT08);

    if(Ntp->NIsolationTrack(signal_idx)!=0){ 
      var_MindcaTrackSV = sqrt( pow(Ntp->IsolationTrack_dzSV(signal_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(signal_idx,TrackIndex),2));
    } else var_MindcaTrackSV = -1;

    if(Ntp->NIsolationTrack(signal_idx)!=0){ 
      var_dcaTrackPV = sqrt(  pow(Ntp->IsolationTrack_dzPV(signal_idx,TrackIndex_closestToPV),2)   +   pow(Ntp->IsolationTrack_dxyPV(signal_idx,TrackIndex_closestToPV),2));
    } else var_dcaTrackPV = -1;



    IsoTrack_P4_closestToMu3.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
				       IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
				       IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.493677*0.493677));


    TLorentzVector IsoTrack_P4_closestToMu3_PiMass = IsoTrack_P4_closestToMu3;
    IsoTrack_P4_closestToMu3_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
					      IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
					      IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.135*0.135));


    TLorentzVector IsoTrack_P4_closestToMu3_MuMass = IsoTrack_P4_closestToMu3;
    IsoTrack_P4_closestToMu3_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
					      IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
					      IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.106*0.106));


    TLorentzVector Mu3_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
    
    Mu3_WithKMass.SetE(sqrt(Mu3_WithKMass.Px()*Mu3_WithKMass.Px() + 
			    Mu3_WithKMass.Py()*Mu3_WithKMass.Py() + 
			    Mu3_WithKMass.Pz()*Mu3_WithKMass.Pz() + 0.493677*0.493677));
    
  
    
    var_IsoPhiKKMass_Mu3 = (IsoTrack_P4_closestToMu3+Mu3_WithKMass).M();
    var_IsoKStarMass_Mu3 = (IsoTrack_P4_closestToMu3_PiMass+Mu3_WithKMass).M();
    var_IsoMuMuMass_Mu3 = (IsoTrack_P4_closestToMu3_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2))).M();
    
    
    
    
    IsoTrack_P4_closestToMu2.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
				       IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
				       IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.493677*0.493677));
    TLorentzVector Mu2_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1));
    
    Mu2_WithKMass.SetE(sqrt(Mu2_WithKMass.Px()*Mu2_WithKMass.Px() + 
			    Mu2_WithKMass.Py()*Mu2_WithKMass.Py() + 
			    Mu2_WithKMass.Pz()*Mu2_WithKMass.Pz() + 0.493677*0.493677));
    
    TLorentzVector IsoTrack_P4_closestToMu2_PiMass = IsoTrack_P4_closestToMu2;
    IsoTrack_P4_closestToMu2_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
					      IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
					      IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.135*0.135));
    
    TLorentzVector IsoTrack_P4_closestToMu2_MuMass = IsoTrack_P4_closestToMu2;
    IsoTrack_P4_closestToMu2_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
					      IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
					      IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.106*0.106));
    
   
    
    var_IsoPhiKKMass_Mu2 = (IsoTrack_P4_closestToMu2+Mu2_WithKMass).M();
    var_IsoKStarMass_Mu2 = (IsoTrack_P4_closestToMu2_PiMass+Mu2_WithKMass).M();
    var_IsoMuMuMass_Mu2 = (IsoTrack_P4_closestToMu2_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))).M();
    


    IsoTrack_P4_closestToMu1.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
				       IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
				       IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.493677*0.493677));
    
    TLorentzVector IsoTrack_P4_closestToMu1_PiMass = IsoTrack_P4_closestToMu1;
    IsoTrack_P4_closestToMu1_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
					      IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
					      IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.135*0.135));
    
    TLorentzVector IsoTrack_P4_closestToMu1_MuMass = IsoTrack_P4_closestToMu1;
    IsoTrack_P4_closestToMu1_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
					      IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
					      IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.106*0.106));
    
    
    
    TLorentzVector Mu1_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0));
    
    Mu1_WithKMass.SetE(sqrt(Mu1_WithKMass.Px()*Mu1_WithKMass.Px() + 
			    Mu1_WithKMass.Py()*Mu1_WithKMass.Py() + 
			    Mu1_WithKMass.Pz()*Mu1_WithKMass.Pz() + 0.493677*0.493677));
    
   
    
    var_IsoPhiKKMass_Mu1 = (IsoTrack_P4_closestToMu1+Mu1_WithKMass).M();
    var_IsoKStarMass_Mu1 = (IsoTrack_P4_closestToMu1_PiMass+Mu1_WithKMass).M();
    var_IsoMuMuMass_Mu1 = (IsoTrack_P4_closestToMu1_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))).M();
    
    var_NtracksClose = NcloseTracksCount;


    //*** define variables for Mu ID and evaluate the BDT

    for(unsigned int imu=0; imu<3;imu++){

      mu_combinedQuality_chi2LocalMomentum=Ntp->Muon_combinedQuality_chi2LocalMomentum(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_chi2LocalPosition=Ntp->Muon_combinedQuality_chi2LocalPosition(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_staRelChi2=Ntp->Muon_combinedQuality_staRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_trkRelChi2=Ntp->Muon_combinedQuality_trkRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_globalDeltaEtaPhi=Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_combinedQuality_trkKink=log(Ntp->Muon_combinedQuality_glbKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)));
      mu_combinedQuality_glbKink=log(Ntp->Muon_combinedQuality_trkKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)));
      mu_combinedQuality_glbTrackProbability=Ntp->Muon_combinedQuality_glbTrackProbability(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_Numberofvalidtrackerhits=Ntp->Muon_numberofValidPixelHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_Numberofvalidpixelhits=Ntp->Muon_innerTrack_numberOfValidTrackerHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_validMuonHitComb=Ntp->Muon_hitPattern_numberOfValidMuonHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_numberOfMatchedStations=Ntp->Muon_numberOfMatchedStations(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_segmentCompatibility=Ntp->Muon_segmentCompatibility(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_timeAtIpInOutErr=Ntp->Muon_timeAtIpInOutErr(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_GLnormChi2=Ntp->Muon_normChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      
      mu_innerTrack_normalizedChi2=Ntp->Muon_innerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_outerTrack_normalizedChi2= Ntp->Muon_outerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));
      mu_innerTrack_validFraction=Ntp->Muon_innerTrack_validFraction(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu));

      if(fabs(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(imu)).Eta()) < 1.2    )
	{
	  if(imu==0)
	    {
	      Muon1DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }
	  if(imu==1)
	    {
	      Muon2DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }

	  if(imu==2)
	    {
	      Muon3DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	    }
	}
      else
	{
	  if(imu==0)
	    {
	      Muon1DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	  if(imu==1)
	    {
	      Muon2DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	  if(imu==2)
	    {
	      Muon3DetID= readerMuIDEndcap->EvaluateMVA("BDT");
	    }
	}
    }

    Muon1MVAID.at(t).Fill(Muon1DetID);
    Muon2MVAID.at(t).Fill(Muon2DetID);
    Muon3MVAID.at(t).Fill(Muon3DetID);

    var_Muon1DetID = Muon1DetID; 
    var_Muon2DetID = Muon2DetID;
    var_Muon3DetID = Muon3DetID;




    // for this particluar binnin Lower

    bool KeepSignalRegionForMC(false);
    if(id!=1) KeepSignalRegionForMC = true;
    if(id==1 && (TauRefitLV.M() > tauMinSideBand_ && TauRefitLV.M() < 1.74) or (TauRefitLV.M() > 1.812 && TauRefitLV.M() < tauMaxSideBand_) ) KeepSignalRegionForMC=true;





  


    //  define also phi veto here

    if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_){


      //      if((dRSortedMassPair1<0.994 ||  dRSortedMassPair1 > 1.044) && (dRSortedMassPair2<0.994 || dRSortedMassPair2> 1.044) ) phiVeto = true;
      if((dRSortedMassPair1<0.993 ||  dRSortedMassPair1 > 1.047) && (dRSortedMassPair2<0.994 || dRSortedMassPair2> 1.044) ) phiVeto = true;
      PairMassDRSorted1A.at(t).Fill(var_mass12_dRsorting,1);
      PairMassDRSorted2A.at(t).Fill(var_mass13_drSorting,1);





    } else if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_){

      //      if((dRSortedMassPair1<0.988 ||  dRSortedMassPair1 > 1.053) && (dRSortedMassPair2<0.985 || dRSortedMassPair2> 1.053) ) phiVeto = true;
      if((dRSortedMassPair1<0.986 ||  dRSortedMassPair1 > 1.054) && (dRSortedMassPair2<0.986 || dRSortedMassPair2> 1.054) ) phiVeto = true;
      PairMassDRSorted1B.at(t).Fill(var_mass12_dRsorting,1);
      PairMassDRSorted2B.at(t).Fill(var_mass13_drSorting,1);




    }else if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_){

      //      if((dRSortedMassPair1<0.974 ||  dRSortedMassPair1 > 1.064) && (dRSortedMassPair2<0.974 || dRSortedMassPair2> 1.064) ) phiVeto = true;
      if((dRSortedMassPair1<0.974 ||  dRSortedMassPair1 > 1.064) && (dRSortedMassPair2<0.974 || dRSortedMassPair2> 1.064) ) phiVeto = true;
      PairMassDRSorted1C.at(t).Fill(var_mass12_dRsorting,1);
      PairMassDRSorted2C.at(t).Fill(var_mass13_drSorting,1);


    }



    var_BvsDSeprator = readerBvsD->EvaluateMVA("BDTG");
    BvsDBDTG.at(t).Fill(var_BvsDSeprator,1);


    




    //    std::cout<<" Btrain   "<< readerBTrainA->EvaluateMVA("BDT")<< "   "<< readerBTrainB->EvaluateMVA("BDT") << "   "<< readerBTrainC->EvaluateMVA("BDT") << std::endl;
    //    std::cout<<" Dtrain   "<< readerDTrainA->EvaluateMVA("BDT")<< "   "<< readerDTrainB->EvaluateMVA("BDT") << "   "<< readerDTrainC->EvaluateMVA("BDT") << std::endl;

    bool BLikeEvent(false);
    bool DLikeEvent(false);

    if(var_BvsDSeprator  > 0) DLikeEvent = true;
    if(var_BvsDSeprator <= 0) BLikeEvent = true;

    //****************************************************   B-D separate plots ***************************
    if(DLikeEvent)
      {
	//category A
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {
	    if(phiVeto)
	      {
		if(readerDTrainA->EvaluateMVA("BDT") > mvaDTrainA2_)
		  {
		    //here to be filled DLike Category A1
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);



		      }
		  }
	      }
	  }



	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		if(readerDTrainB->EvaluateMVA("BDT") > mvaDTrainB2_)
		  {
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);
		      }
		    //here to be filled DLike Category B1
		  }
	      }
	  }
	


        if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
          {
	    if(phiVeto)
	      {
		if(readerDTrainC->EvaluateMVA("BDT") > mvaDTrainB2_)
		  {
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);
		      }
                    //here to be filled DLike Category C1
                  }
              }
          }

	// ==================== D like subcategory 2


	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {
	    if(phiVeto)
	      {
		
		if(readerDTrainA->EvaluateMVA("BDT") > mvaDTrainA1_ && readerDTrainA->EvaluateMVA("BDT") < mvaDTrainA2_ )
		  {
		    //here to be filled DLike Category A2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }



	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		
		if(readerDTrainB->EvaluateMVA("BDT") > mvaDTrainB1_ && readerDTrainB->EvaluateMVA("BDT") < mvaDTrainB2_ )
		  {
		    //here to be filled DLike Category B2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }


	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		
		if(readerDTrainC->EvaluateMVA("BDT") > mvaDTrainC1_ && readerDTrainC->EvaluateMVA("BDT") < mvaDTrainC2_ )
		  {
		    //here to be filled DLike Category C2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }


      }

    ////////////////////////////////////////////  Now B only


    if(BLikeEvent)
      {
	//category A
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {
	    if(phiVeto)
	      {
		
		if(readerDTrainB->EvaluateMVA("BDT") > mvaBTrainA2_)
		  {
		    //here to be filled BLike Category A1
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		
		if(readerBTrainB->EvaluateMVA("BDT") > mvaBTrainB2_)
		  {
		    //here to be filled BLike Category B1
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }



        if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
          {
	    if(phiVeto)
	      {
		
		if(readerBTrainC->EvaluateMVA("BDT") > mvaBTrainB2_)
		  {
                    //here to be filled BLike Category C1
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC1_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC1.at(t).Fill(var_BvsDSeprator,1);
		      }
                  }
              }
          }

	// ==================== B like subcategory 2


	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {
	    if(phiVeto)
	      {
		
		if(readerBTrainA->EvaluateMVA("BDT") > mvaBTrainA1_ && readerBTrainA->EvaluateMVA("BDT") < mvaBTrainA2_ )
		  {
		    //here to be filled BLike Category A2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }



	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		if(readerBTrainB->EvaluateMVA("BDT") > mvaBTrainB1_ && readerBTrainB->EvaluateMVA("BDT") < mvaBTrainB2_ )
		  {
		    //here to be filled BLike Category B2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }
	
	
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {
	    if(phiVeto)
	      {
		if(readerBTrainC->EvaluateMVA("BDT") > mvaBTrainC1_ && readerBTrainC->EvaluateMVA("BDT") < mvaBTrainC2_ )
		  {
		    //here to be filled DLike Category C2
		    if(KeepSignalRegionForMC)
		      {
			TauMassRefitABC2_BDSeparateTrain.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
			BvsDBDTG_ABC2.at(t).Fill(var_BvsDSeprator,1);
		      }
		  }
	      }
	  }


      }






    //****************************************************   B-D separate plots ***************************


    if(phiVeto)
      {
	//A
        if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
          {




	    TauMass_vs_BDT_train0_A.at(t).Fill(TauRefitLV.M(),   reader_train0_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train1_A.at(t).Fill(TauRefitLV.M(),   reader_train1_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train2_A.at(t).Fill(TauRefitLV.M(),   reader_train2_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train3_A.at(t).Fill(TauRefitLV.M(),   reader_train3_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train4_A.at(t).Fill(TauRefitLV.M(),   reader_train4_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train5_A.at(t).Fill(TauRefitLV.M(),   reader_train5_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train6_A.at(t).Fill(TauRefitLV.M(),   reader_train6_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train7_A.at(t).Fill(TauRefitLV.M(),   reader_train7_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train8_A.at(t).Fill(TauRefitLV.M(),   reader_train8_A ->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train9_A.at(t).Fill(TauRefitLV.M(),   reader_train9_A ->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train10_A.at(t).Fill(TauRefitLV.M(),  reader_train10_A->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train11_A.at(t).Fill(TauRefitLV.M(),  reader_train11_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train12_A.at(t).Fill(TauRefitLV.M(),  reader_train12_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train13_A.at(t).Fill(TauRefitLV.M(),  reader_train13_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train14_A.at(t).Fill(TauRefitLV.M(),  reader_train14_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train15_A.at(t).Fill(TauRefitLV.M(),  reader_train15_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train16_A.at(t).Fill(TauRefitLV.M(),  reader_train16_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train17_A.at(t).Fill(TauRefitLV.M(),  reader_train17_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train18_A.at(t).Fill(TauRefitLV.M(),  reader_train18_A->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train19_A.at(t).Fill(TauRefitLV.M(),  reader_train19_A->EvaluateMVA("BDT")  );

	    BDTOutput_train0_A.at(t).Fill(     reader_train0_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train1_A.at(t).Fill(     reader_train1_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train2_A.at(t).Fill(     reader_train2_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train3_A.at(t).Fill(     reader_train3_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train4_A.at(t).Fill(     reader_train4_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train5_A.at(t).Fill(     reader_train5_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train6_A.at(t).Fill(     reader_train6_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train7_A.at(t).Fill(     reader_train7_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train8_A.at(t).Fill(     reader_train8_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train9_A.at(t).Fill(     reader_train9_A ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train10_A.at(t).Fill(    reader_train10_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train11_A.at(t).Fill(    reader_train11_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train12_A.at(t).Fill(    reader_train12_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train13_A.at(t).Fill(    reader_train13_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train14_A.at(t).Fill(    reader_train14_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train15_A.at(t).Fill(    reader_train15_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train16_A.at(t).Fill(    reader_train16_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train17_A.at(t).Fill(    reader_train17_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train18_A.at(t).Fill(    reader_train18_A->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train19_A.at(t).Fill(    reader_train19_A->EvaluateMVA("BDT"), 1 );

	    BDTOutput_train19_A.at(t).Fill(    reader_train19_A->EvaluateMVA("BDT"), 1 );



	  }

	//B
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {

	    TauMass_vs_BDT_train0_B.at(t).Fill(TauRefitLV.M(),  reader_train0_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train1_B.at(t).Fill(TauRefitLV.M(),  reader_train1_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train2_B.at(t).Fill(TauRefitLV.M(),  reader_train2_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train3_B.at(t).Fill(TauRefitLV.M(),  reader_train3_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train4_B.at(t).Fill(TauRefitLV.M(),  reader_train4_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train5_B.at(t).Fill(TauRefitLV.M(),  reader_train5_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train6_B.at(t).Fill(TauRefitLV.M(),  reader_train6_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train7_B.at(t).Fill(TauRefitLV.M(),  reader_train7_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train8_B.at(t).Fill(TauRefitLV.M(),  reader_train8_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train9_B.at(t).Fill(TauRefitLV.M(),  reader_train9_B->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train10_B.at(t).Fill(TauRefitLV.M(),  reader_train10_B->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train11_B.at(t).Fill(TauRefitLV.M(),  reader_train11_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train12_B.at(t).Fill(TauRefitLV.M(),  reader_train12_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train13_B.at(t).Fill(TauRefitLV.M(),  reader_train13_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train14_B.at(t).Fill(TauRefitLV.M(),  reader_train14_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train15_B.at(t).Fill(TauRefitLV.M(),  reader_train15_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train16_B.at(t).Fill(TauRefitLV.M(),  reader_train16_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train17_B.at(t).Fill(TauRefitLV.M(),  reader_train17_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train18_B.at(t).Fill(TauRefitLV.M(),  reader_train18_B->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train19_B.at(t).Fill(TauRefitLV.M(),  reader_train19_B->EvaluateMVA("BDT")  );

	    BDTOutput_train0_B.at(t).Fill(     reader_train0_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train1_B.at(t).Fill(     reader_train1_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train2_B.at(t).Fill(     reader_train2_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train3_B.at(t).Fill(     reader_train3_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train4_B.at(t).Fill(     reader_train4_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train5_B.at(t).Fill(     reader_train5_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train6_B.at(t).Fill(     reader_train6_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train7_B.at(t).Fill(     reader_train7_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train8_B.at(t).Fill(     reader_train8_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train9_B.at(t).Fill(     reader_train9_B ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train10_B.at(t).Fill(    reader_train10_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train11_B.at(t).Fill(    reader_train11_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train12_B.at(t).Fill(    reader_train12_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train13_B.at(t).Fill(    reader_train13_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train14_B.at(t).Fill(    reader_train14_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train15_B.at(t).Fill(    reader_train15_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train16_B.at(t).Fill(    reader_train16_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train17_B.at(t).Fill(    reader_train17_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train18_B.at(t).Fill(    reader_train18_B->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train19_B.at(t).Fill(    reader_train19_B->EvaluateMVA("BDT"), 1 );
	  }



	//C

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {
	    TauMass_vs_BDT_train0_C.at(t).Fill(TauRefitLV.M(),  reader_train0_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train1_C.at(t).Fill(TauRefitLV.M(),  reader_train1_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train2_C.at(t).Fill(TauRefitLV.M(),  reader_train2_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train3_C.at(t).Fill(TauRefitLV.M(),  reader_train3_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train4_C.at(t).Fill(TauRefitLV.M(),  reader_train4_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train5_C.at(t).Fill(TauRefitLV.M(),  reader_train5_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train6_C.at(t).Fill(TauRefitLV.M(),  reader_train6_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train7_C.at(t).Fill(TauRefitLV.M(),  reader_train7_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train8_C.at(t).Fill(TauRefitLV.M(),  reader_train8_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train9_C.at(t).Fill(TauRefitLV.M(),  reader_train9_C->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train10_C.at(t).Fill(TauRefitLV.M(),  reader_train10_C->EvaluateMVA("BDT")  );
	    TauMass_vs_BDT_train11_C.at(t).Fill(TauRefitLV.M(),  reader_train11_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train12_C.at(t).Fill(TauRefitLV.M(),  reader_train12_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train13_C.at(t).Fill(TauRefitLV.M(),  reader_train13_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train14_C.at(t).Fill(TauRefitLV.M(),  reader_train14_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train15_C.at(t).Fill(TauRefitLV.M(),  reader_train15_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train16_C.at(t).Fill(TauRefitLV.M(),  reader_train16_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train17_C.at(t).Fill(TauRefitLV.M(),  reader_train17_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train18_C.at(t).Fill(TauRefitLV.M(),  reader_train18_C->EvaluateMVA("BDT")  );
            TauMass_vs_BDT_train19_C.at(t).Fill(TauRefitLV.M(),  reader_train19_C->EvaluateMVA("BDT")  );

	    BDTOutput_train0_C.at(t).Fill(     reader_train0_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train1_C.at(t).Fill(     reader_train1_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train2_C.at(t).Fill(     reader_train2_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train3_C.at(t).Fill(     reader_train3_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train4_C.at(t).Fill(     reader_train4_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train5_C.at(t).Fill(     reader_train5_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train6_C.at(t).Fill(     reader_train6_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train7_C.at(t).Fill(     reader_train7_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train8_C.at(t).Fill(     reader_train8_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train9_C.at(t).Fill(     reader_train9_C ->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train10_C.at(t).Fill(    reader_train10_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train11_C.at(t).Fill(    reader_train11_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train12_C.at(t).Fill(    reader_train12_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train13_C.at(t).Fill(    reader_train13_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train14_C.at(t).Fill(    reader_train14_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train15_C.at(t).Fill(    reader_train15_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train16_C.at(t).Fill(    reader_train16_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train17_C.at(t).Fill(    reader_train17_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train18_C.at(t).Fill(    reader_train18_C->EvaluateMVA("BDT"), 1 );
	    BDTOutput_train19_C.at(t).Fill(    reader_train19_C->EvaluateMVA("BDT"), 1 ); 
	  }

      }




    double dRSortedMass;

    //*** define per event resolution categroies 
    //*** Category A1
    //    if(phiVeto && rhoVeto)




    //*** defined the pair with SS best alligned to OS
    double dRSortedPairMass1(0.);
    double dRSortedPairMass2(0.);

    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
      dRSortedPairMass1 = (MuonOS+MuonSS2).M();
      dRSortedPairMass2 = (MuonOS+MuonSS1).M();
    }else{
      dRSortedPairMass2 = (MuonOS+MuonSS2).M();
      dRSortedPairMass1 = (MuonOS+MuonSS1).M();
    }

    if(phiVeto)
      {
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {

	    TauMass_allVsBDTA.at(t).Fill(TauRefitLV.M(),readerA->EvaluateMVA("BDT"));
	    BDTOutputA.at(t).Fill(    readerA->EvaluateMVA("BDT"),1 );
	    
	    if(readerA->EvaluateMVA("BDT") > mvaA2_){

		PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
		
		AllignSortMass1.at(t).Fill(dRSortedPairMass1,1);


		std::cout<<" ================================= pair mass "<< dRSortedPairMass1 << std::endl;
		if(id == 120){
		  if(dRSortedPairMass1 < 0.550){

		    std::cout<<"-------------- A1---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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


		AllignSortMass2.at(t).Fill(dRSortedPairMass2,1);

		TauMassA1.at(t).Fill(TauLV.M(),1);                 // three mu mass 
		TauMassRefitA1.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
		if(KeepSignalRegionForMC)
		  {
		    TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		    TauMassRefitABC1_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));

		    if(id == 120){
		      std::cout<<"-------------- A1 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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


		if(dRSortedMass > 0.549) 
		  {

		    TauMassRefitABC1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
		  }
		
	    }
	  }




	    //-----------


	    BDTOutputA_train11.at(t).Fill(    readerA_train11->EvaluateMVA("BDT") );
	    if(readerA_train11->EvaluateMVA("BDT") > mvaA2train11_){
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitA1_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }
      
	    
	    BDTOutputA_train8.at(t).Fill(    readerA_train8->EvaluateMVA("BDT") );
	    if(readerA_train8->EvaluateMVA("BDT") > mvaA2train8_){
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitA1_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }
      
	    

  

      

	//Category B1
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {

	    TauMass_allVsBDTB.at(t).Fill(TauRefitLV.M(),readerB->EvaluateMVA("BDT"));
	    BDTOutputB.at(t).Fill(readerB->EvaluateMVA("BDT"), 1);

	    if(readerB->EvaluateMVA("BDT") > mvaB2_){
	      
	      PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
	      PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
	      PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	      
	      //*** defined the pair with SS best alligned to OS
	      if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		dRSortedMass = (MuonOS+MuonSS2).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS2).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS1).M(),1);
		
	      }else{
		dRSortedMass = (MuonOS+MuonSS1).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS1).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS2).M(),1);
		
		
	      }
	      //***
	      
	      if(id == 120){
		if(dRSortedPairMass1 < 0.550){

		  std::cout<<"-------------- B1---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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
	      
	      
	      TauMassB1.at(t).Fill(TauLV.M(),1);                  // three mu mass 
	      TauMassRefitB1.at(t).Fill(TauRefitLV.M(),1);        // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);      // fill up all categories inclusive
		  TauMassRefitABC1_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));
		  		  if(id == 120){
		    std::cout<<"-------------- B1 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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


	      if(dRSortedMass > 0.549) {

		TauMassRefitABC1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
	      }
	    }




  
	    BDTOutputB_train11.at(t).Fill(    readerB_train11->EvaluateMVA("BDT") );
	    if(readerB_train11->EvaluateMVA("BDT") > mvaB2train11_){
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitB1_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		    {
		      TauMassRefitABC1_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		    }
	    }

	    BDTOutputB_train8.at(t).Fill(    readerB_train8->EvaluateMVA("BDT") );
	    if(readerB_train8 ->EvaluateMVA("BDT") > mvaB2train8_){
	      AllignSortMass1_train8.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train8.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitB1_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		    {
		      TauMassRefitABC1_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		    }
	    }



	  }
	

	
	
	
	//Category C1
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_)
	  {

	    TauMass_allVsBDTC.at(t).Fill(TauRefitLV.M(),readerC->EvaluateMVA("BDT"));
	    BDTOutputC.at(t).Fill(    readerC->EvaluateMVA("BDT") );

	    if(readerC->EvaluateMVA("BDT") > mvaC2_){
	      
	      PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
	      PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
	      PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	      
	      //*** defined the pair with SS best alligned to OS
	      if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		dRSortedMass = (MuonOS+MuonSS2).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS2).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS1).M(),1);
		
		
	      }else{
		dRSortedMass = (MuonOS+MuonSS1).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS1).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS2).M(),1);
		
	      }
	      //***
	      
              if(id == 120){
                if(dRSortedPairMass1 < 0.550){

		  std::cout<<"-------------- C1---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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

	      
	      TauMassC1.at(t).Fill(TauLV.M(),1);	          // three mu mass 
	      TauMassRefitC1.at(t).Fill(TauRefitLV.M(),1);      // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1.at(t).Fill(TauRefitLV.M(),1);    // fill up all categories inclusive
		  TauMassRefitABC1_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));
		  		  if(id == 120){
		    std::cout<<"-------------- C1 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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


	      if(dRSortedMass > 0.549){

		TauMassRefitABC1FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
	      }
	    }
	    


  
	    BDTOutputC_train11.at(t).Fill(    readerC_train11->EvaluateMVA("BDT") );
	    if(readerC_train11->EvaluateMVA("BDT") > mvaC2train11_){
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitC1_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }

  
	    BDTOutputC_train8.at(t).Fill(    readerC_train8->EvaluateMVA("BDT") );
	    if(readerC_train8->EvaluateMVA("BDT") > mvaC2train8_){
	      AllignSortMass1_train8.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train8.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitC1_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC1_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }


	  }
	    
      

	//Category A2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_)
	  {

	    if(readerA->EvaluateMVA("BDT") > mvaA1_ && readerA->EvaluateMVA("BDT") < mvaA2_){

	      PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
	      PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
	      PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	      
	      //*** defined the pair with SS best alligned to OS
	      if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		dRSortedMass = (MuonOS+MuonSS2).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS2).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS1).M(),1);
		
	      }else{
		dRSortedMass = (MuonOS+MuonSS1).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS1).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS2).M(),1);
		
	      }
	      //***
	      
	      
	      if(id == 120){
		if(dRSortedPairMass1 < 0.550){

		  std::cout<<"-------------- A2---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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


	      
	      TauMassA2.at(t).Fill(TauLV.M(),1);
	      TauMassRefitA2.at(t).Fill(TauRefitLV.M(),1);    
	      if(KeepSignalRegionForMC){
		TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    
		TauMassRefitABC2_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));
				if(id == 120){
		  std::cout<<"-------------- A2 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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


	      if(dRSortedMass > 0.549){

		TauMassRefitABC2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
	      }
	      
	    }
	  
	    


	    if(readerA_train11->EvaluateMVA("BDT") > mvaA1train11_ && readerA_train11->EvaluateMVA("BDT") < mvaA2train11_ ){
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
		  
	      TauMassRefitA2_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }


	    if(readerA_train8->EvaluateMVA("BDT") > mvaA1train8_ && readerA_train8->EvaluateMVA("BDT") < mvaA2train8_ ){
	      AllignSortMass1_train8.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train8.at(t).Fill(dRSortedPairMass2,1);
		  
	      TauMassRefitA2_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }


	  }
	    



      

	//Category B2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_)
	  {
	    if(readerB->EvaluateMVA("BDT") > mvaB1_ && readerB->EvaluateMVA("BDT") < mvaB2_){

		  PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
		  PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
		  PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	  
                  //*** defined the pair with SS best alligned to OS
                  if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
                    dRSortedMass = (MuonOS+MuonSS2).M();
                    AllignSortMass1.at(t).Fill((MuonOS+MuonSS2).M(),1);
                    AllignSortMass2.at(t).Fill((MuonOS+MuonSS1).M(),1);

                  }else{
                    dRSortedMass = (MuonOS+MuonSS1).M();
                    AllignSortMass1.at(t).Fill((MuonOS+MuonSS1).M(),1);
                    AllignSortMass2.at(t).Fill((MuonOS+MuonSS2).M(),1);

                  }


		  if(id == 120){
		    if(dRSortedPairMass1 < 0.550){

		      std::cout<<"-------------- B2---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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


                  //***

		  //		  std::cout<<"dRSortedMass  "<<dRSortedMass << std::endl;

		  TauMassB2.at(t).Fill(TauLV.M(),1);
		  TauMassRefitB2.at(t).Fill(TauRefitLV.M(),1);    
		  if(KeepSignalRegionForMC)
		    {
		      TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    
		      TauMassRefitABC2_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));
		      		      if(id == 120){
			std::cout<<"-------------- B2 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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


		  if(dRSortedMass > 0.549){

		    TauMassRefitABC2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
		  }

	    }
	  



	    if(readerB_train11->EvaluateMVA("BDT") > mvaB1train11_ && readerB_train11->EvaluateMVA("BDT") < mvaB2train11_ ){
	      
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitB2_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }



	    if(readerB_train8->EvaluateMVA("BDT") > mvaB1train8_ && readerB_train8->EvaluateMVA("BDT") < mvaB2train8_ ){
	      
	      AllignSortMass1_train8.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train8.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitB2_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }


	  }
	    

    
	//Category C2
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_){

	    if(readerC->EvaluateMVA("BDT") > mvaC1_ && readerC->EvaluateMVA("BDT")< mvaC2_){

	      PairMass1.at(t).Fill((MuonOS+MuonSS1).M() ,1);
	      PairMass2.at(t).Fill((MuonOS+MuonSS2).M() ,1);
	      PairMassFinalSel.at(t).Fill((MuonOS+MuonSS1).M(), (MuonOS+MuonSS2).M() ,1);
	      
	      //*** defined the pair with SS best alligned to OS
	      if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
		dRSortedMass = (MuonOS+MuonSS2).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS2).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS1).M(),1);
		
	      }else{
		dRSortedMass = (MuonOS+MuonSS1).M();
		AllignSortMass1.at(t).Fill((MuonOS+MuonSS1).M(),1);
		AllignSortMass2.at(t).Fill((MuonOS+MuonSS2).M(),1);
		
	      }
	      //***
	      


              if(id == 120){
                if(dRSortedPairMass1 < 0.550){

		  std::cout<<"-------------- C2---------------- pair mass:  "<< dRSortedPairMass1 <<std::endl;
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


	      
	      TauMassC2.at(t).Fill(TauLV.M(),1);	      
	      TauMassRefitC2.at(t).Fill(TauRefitLV.M(),1);    
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2.at(t).Fill(TauRefitLV.M(),1);    
		  TauMassRefitABC2_eta.at(t).Fill(TauRefitLV.M(),fabs(TauRefitLV.Eta()));
		  		  if(id == 120){
		    std::cout<<"-------------- C2 ----------------mass:  "<< (Muon1LV + Muon2LV + Muon3LV).M() <<std::endl;
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
	      

	      if(dRSortedMass > 0.549){

		TauMassRefitABC2FullEtaVetoCut.at(t).Fill(TauRefitLV.M(),1);
	      }
	    }
	    
	
	
	    if(readerC_train11->EvaluateMVA("BDT") > mvaC1train11_ && readerC_train11->EvaluateMVA("BDT") < mvaC2train11_ ){
	      
	      AllignSortMass1_train11.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train11.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitC2_train11.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train11.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }


	    if(readerC_train8->EvaluateMVA("BDT") > mvaC1train8_ && readerC_train8 ->EvaluateMVA("BDT") < mvaC2train8_ ){
	      
	      AllignSortMass1_train8.at(t).Fill(dRSortedPairMass1,1);
	      AllignSortMass2_train8.at(t).Fill(dRSortedPairMass2,1);
	      
	      TauMassRefitC2_train8.at(t).Fill(TauRefitLV.M(),1);       // three mu KF reffited mass
	      if(KeepSignalRegionForMC)
		{
		  TauMassRefitABC2_train8.at(t).Fill(TauRefitLV.M(),1);     // fill up all categories inclusive
		}
	    }




	}
      }
      
  

	//*** below are basic  purity and resolution plots
	if(id==40 || id == 60 || id ==90){
	  if(Ntp->MCEventIsReconstructed()){
	    TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
	    TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
	    TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
	    TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
	
	    TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	    TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);
	    //	    TauMassResolutionHelixRefit.at(t).Fill((MotherParticle.LV().M() - MCTauLV.M())/MCTauLV.M(),1);


	    Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
	    Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
	    Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);
	  }
	}
    
    
	//*** fill up the T3MMiniTree.root for statistical analysis
    
	m3m = TauRefitLV.M();
	if(CrossVeto)	xv = 1;
	if(!CrossVeto)	xv = 0;
	
	if(phiVeto) phiv=1;
	if(!phiVeto) phiv=0;

	dataMCtype = id;
	event_weight =1; // 1 for data
	if(dataMCtype == 1){event_weight =1;}
	else if(dataMCtype == 40){event_weight =9.2e-04;} // event_weight is a value Lumi Scale 
	else if(dataMCtype == 60){event_weight =4.61e-04;}
	else if(dataMCtype == 90){event_weight =6.25e-04;}
    
    
	mvaA1 = mvaA1_;
	mvaA2 = mvaA2_;
	mvaB1 = mvaB1_;
	mvaB2 = mvaB2_;
	mvaC1 = mvaC1_;
	mvaC2 = mvaC2_;
    

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut1_ ){
	  category = 1;
	  bdt = readerA->EvaluateMVA("BDT");
	}
    
	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut1_ && Ntp->TauMassResolution(EtaSortedIndices,1,false) < PEMassResolutionCut2_){
	  category = 2;
	  bdt = readerB->EvaluateMVA("BDT");
	}

	if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > PEMassResolutionCut2_){
	  category = 3;
	  bdt = readerC->EvaluateMVA("BDT");
	}


	if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
	  mDr1 = (MuonOS+MuonSS2).M();
	  mDr2 = (MuonOS+MuonSS1).M();
	}else{
	  mDr1 = (MuonOS+MuonSS1).M();
	  mDr2 = (MuonOS+MuonSS2).M();
	}


	m12 = (MuonOS+MuonSS1).M();
	m13 = (MuonOS+MuonSS2).M();
	LumiScale = 1.;
	T3MMiniTree->Fill();
  }
}



void  SignalSelector::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis
  T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();


  //*** extra actions
  if(mode == RECONSTRUCT){

    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    //    int id(Ntp->GetMCID());
    //    double scale(1.);
    //    double scaleDsTau(0.637);
    //    double scaleBpTau(0.262);
    //    double scaleB0Tau(0.099);
    //total xsection of producing taus is 12.848 ub 
    // if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    // ScaleAllHistOfType(2,scale*scaleDsTau);
    // if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    // ScaleAllHistOfType(3,scale*scaleB0Tau);
    // if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    // ScaleAllHistOfType(4,scale*scaleBpTau);
    //    }
  }
  Selection::Finish();
}





