// T3MSelectionTree
#include "T3MSelectionTree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

T3MSelectionTree::T3MSelectionTree(TString Name_, TString id_):
   Selection(Name_,id_),
   tauMinMass_(1.75),
   tauMaxMass_(1.80),
   tauMinSideBand_(1.62),
   tauMaxSideBand_(2.0),
   tauMassResCutLow(0.007),
   tauMassResCutHigh(0.0105),
   phiVetoSigma(0.022),
   omegaVetoSigma(0.017)
{
   // This is a class constructor;
   TString basedir = "";
   basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2018";
   PUWeightFile = new TFile(basedir+"/PUWeights_Run2018.root");
   puWeights = (TH1D*)PUWeightFile->Get("h1_weights");
   file= new TFile("T3MSelectionTreeInput_v2.root","recreate");
}

T3MSelectionTree::~T3MSelectionTree(){
   for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
         << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
         << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
   }
   Logger(Logger::Info) << "complete." << std::endl;
}
void  T3MSelectionTree::Configure(){
   // Set tree branches

   TMVA_Tree= new TTree("tree","tree");

   TMVA_Tree->Branch("MC",&MC);
   TMVA_Tree->Branch("category",&category);
   TMVA_Tree->Branch("DataMCType",&DataMCType);
   TMVA_Tree->Branch("threeGlobal",&threeGlobal);
   TMVA_Tree->Branch("EventWeight",&EventWeight);
   TMVA_Tree->Branch("l1seed",&l1seed);
   TMVA_Tree->Branch("var_vertexKFChi2",&var_vertexKFChi2);
   TMVA_Tree->Branch("var_svpvTauAngle",&var_svpvTauAngle);
   TMVA_Tree->Branch("l1seed",&l1seed);
   TMVA_Tree->Branch("run",&run);
   TMVA_Tree->Branch("lumi",&lumi);
   TMVA_Tree->Branch("eventNumber", &eventNumber);
   TMVA_Tree->Branch("var_Mosss1", &var_Mosss1);
   TMVA_Tree->Branch("var_Mosss2", &var_Mosss2);

   //commmon variables
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
   TMVA_Tree->Branch("var_nMatches_mu3",&var_nMatches_mu3); //Number of matches (mu3)

   TMVA_Tree->Branch("var_tauMass",&var_tauMass);
   TMVA_Tree->Branch("var_ntracks",&var_ntracks);

   TMVA_Tree->Branch("var_nsv",&var_nsv);
   TMVA_Tree->Branch("var_VertexMu1D0SigPVReco",&var_VertexMu1D0SigPVReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigPVReco",&var_VertexMu2D0SigPVReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigPVReco",&var_VertexMu3D0SigPVReco);

   TMVA_Tree->Branch("var_VertexMu1D0SigBSReco",&var_VertexMu1D0SigBSReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigBSReco",&var_VertexMu2D0SigBSReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigBSReco",&var_VertexMu3D0SigBSReco);

   TMVA_Tree->Branch("var_VertexMu1D0SigSVReco",&var_VertexMu1D0SigSVReco);
   TMVA_Tree->Branch("var_VertexMu2D0SigSVReco",&var_VertexMu2D0SigSVReco);
   TMVA_Tree->Branch("var_VertexMu3D0SigSVReco",&var_VertexMu3D0SigSVReco);

   // Muon station information (Muon1)
   TMVA_Tree->Branch("Muon1_station1_status", &Muon1_station_vars[0][0]);
   TMVA_Tree->Branch("Muon1_station1_TrackX", &Muon1_station_vars[0][1]);
   TMVA_Tree->Branch("Muon1_station1_TrackY", &Muon1_station_vars[0][2]);
   TMVA_Tree->Branch("Muon1_station1_pullX", &Muon1_station_vars[0][3]);
   TMVA_Tree->Branch("Muon1_station1_pullY", &Muon1_station_vars[0][4]);
   TMVA_Tree->Branch("Muon1_station1_pullDxDz", &Muon1_station_vars[0][5]);
   TMVA_Tree->Branch("Muon1_station1_pullDyDz", &Muon1_station_vars[0][6]);

   TMVA_Tree->Branch("Muon1_station2_status", &Muon1_station_vars[1][0]);
   TMVA_Tree->Branch("Muon1_station2_TrackX", &Muon1_station_vars[1][1]);
   TMVA_Tree->Branch("Muon1_station2_TrackY", &Muon1_station_vars[1][2]);
   TMVA_Tree->Branch("Muon1_station2_pullX", &Muon1_station_vars[1][3]);
   TMVA_Tree->Branch("Muon1_station2_pullY", &Muon1_station_vars[1][4]);
   TMVA_Tree->Branch("Muon1_station2_pullDxDz", &Muon1_station_vars[1][5]);
   TMVA_Tree->Branch("Muon1_station2_pullDyDz", &Muon1_station_vars[1][6]);

   TMVA_Tree->Branch("Muon1_station3_status", &Muon1_station_vars[2][0]);
   TMVA_Tree->Branch("Muon1_station3_TrackX", &Muon1_station_vars[2][1]);
   TMVA_Tree->Branch("Muon1_station3_TrackY", &Muon1_station_vars[2][2]);
   TMVA_Tree->Branch("Muon1_station3_pullX", &Muon1_station_vars[2][3]);
   TMVA_Tree->Branch("Muon1_station3_pullY", &Muon1_station_vars[2][4]);
   TMVA_Tree->Branch("Muon1_station3_pullDxDz", &Muon1_station_vars[2][5]);
   TMVA_Tree->Branch("Muon1_station3_pullDyDz", &Muon1_station_vars[2][6]);

   TMVA_Tree->Branch("Muon1_station4_status", &Muon1_station_vars[3][0]);
   TMVA_Tree->Branch("Muon1_station4_TrackX", &Muon1_station_vars[3][1]);
   TMVA_Tree->Branch("Muon1_station4_TrackY", &Muon1_station_vars[3][2]);
   TMVA_Tree->Branch("Muon1_station4_pullX", &Muon1_station_vars[3][3]);
   TMVA_Tree->Branch("Muon1_station4_pullY", &Muon1_station_vars[3][4]);
   TMVA_Tree->Branch("Muon1_station4_pullDxDz", &Muon1_station_vars[3][5]);
   TMVA_Tree->Branch("Muon1_station4_pullDyDz", &Muon1_station_vars[3][6]);

   TMVA_Tree->Branch("Muon1_station5_status", &Muon1_station_vars[4][0]);
   TMVA_Tree->Branch("Muon1_station5_TrackX", &Muon1_station_vars[4][1]);
   TMVA_Tree->Branch("Muon1_station5_TrackY", &Muon1_station_vars[4][2]);
   TMVA_Tree->Branch("Muon1_station5_pullX", &Muon1_station_vars[4][3]);
   TMVA_Tree->Branch("Muon1_station5_pullY", &Muon1_station_vars[4][4]);
   TMVA_Tree->Branch("Muon1_station5_pullDxDz", &Muon1_station_vars[4][5]);
   TMVA_Tree->Branch("Muon1_station5_pullDyDz", &Muon1_station_vars[4][6]);

   TMVA_Tree->Branch("Muon1_station6_status", &Muon1_station_vars[5][0]);
   TMVA_Tree->Branch("Muon1_station6_TrackX", &Muon1_station_vars[5][1]);
   TMVA_Tree->Branch("Muon1_station6_TrackY", &Muon1_station_vars[5][2]);
   TMVA_Tree->Branch("Muon1_station6_pullX", &Muon1_station_vars[5][3]);
   TMVA_Tree->Branch("Muon1_station6_pullY", &Muon1_station_vars[5][4]);
   TMVA_Tree->Branch("Muon1_station6_pullDxDz", &Muon1_station_vars[5][5]);
   TMVA_Tree->Branch("Muon1_station6_pullDyDz", &Muon1_station_vars[5][6]);

   TMVA_Tree->Branch("Muon1_station7_status", &Muon1_station_vars[6][0]);
   TMVA_Tree->Branch("Muon1_station7_TrackX", &Muon1_station_vars[6][1]);
   TMVA_Tree->Branch("Muon1_station7_TrackY", &Muon1_station_vars[6][2]);
   TMVA_Tree->Branch("Muon1_station7_pullX", &Muon1_station_vars[6][3]);
   TMVA_Tree->Branch("Muon1_station7_pullY", &Muon1_station_vars[6][4]);
   TMVA_Tree->Branch("Muon1_station7_pullDxDz", &Muon1_station_vars[6][5]);
   TMVA_Tree->Branch("Muon1_station7_pullDyDz", &Muon1_station_vars[6][6]);

   TMVA_Tree->Branch("Muon1_station8_status", &Muon1_station_vars[7][0]);
   TMVA_Tree->Branch("Muon1_station8_TrackX", &Muon1_station_vars[7][1]);
   TMVA_Tree->Branch("Muon1_station8_TrackY", &Muon1_station_vars[7][2]);
   TMVA_Tree->Branch("Muon1_station8_pullX", &Muon1_station_vars[7][3]);
   TMVA_Tree->Branch("Muon1_station8_pullY", &Muon1_station_vars[7][4]);
   TMVA_Tree->Branch("Muon1_station8_pullDxDz", &Muon1_station_vars[7][5]);
   TMVA_Tree->Branch("Muon1_station8_pullDyDz", &Muon1_station_vars[7][6]);

   // Muon station information (Muon2)
   TMVA_Tree->Branch("Muon2_station1_status", &Muon2_station_vars[0][0]);
   TMVA_Tree->Branch("Muon2_station1_TrackX", &Muon2_station_vars[0][1]);
   TMVA_Tree->Branch("Muon2_station1_TrackY", &Muon2_station_vars[0][2]);
   TMVA_Tree->Branch("Muon2_station1_pullX", &Muon2_station_vars[0][3]);
   TMVA_Tree->Branch("Muon2_station1_pullY", &Muon2_station_vars[0][4]);
   TMVA_Tree->Branch("Muon2_station1_pullDxDz", &Muon2_station_vars[0][5]);
   TMVA_Tree->Branch("Muon2_station1_pullDyDz", &Muon2_station_vars[0][6]);

   TMVA_Tree->Branch("Muon2_station2_status", &Muon2_station_vars[1][0]);
   TMVA_Tree->Branch("Muon2_station2_TrackX", &Muon2_station_vars[1][1]);
   TMVA_Tree->Branch("Muon2_station2_TrackY", &Muon2_station_vars[1][2]);
   TMVA_Tree->Branch("Muon2_station2_pullX", &Muon2_station_vars[1][3]);
   TMVA_Tree->Branch("Muon2_station2_pullY", &Muon2_station_vars[1][4]);
   TMVA_Tree->Branch("Muon2_station2_pullDxDz", &Muon2_station_vars[1][5]);
   TMVA_Tree->Branch("Muon2_station2_pullDyDz", &Muon2_station_vars[1][6]);

   TMVA_Tree->Branch("Muon2_station3_status", &Muon2_station_vars[2][0]);
   TMVA_Tree->Branch("Muon2_station3_TrackX", &Muon2_station_vars[2][1]);
   TMVA_Tree->Branch("Muon2_station3_TrackY", &Muon2_station_vars[2][2]);
   TMVA_Tree->Branch("Muon2_station3_pullX", &Muon2_station_vars[2][3]);
   TMVA_Tree->Branch("Muon2_station3_pullY", &Muon2_station_vars[2][4]);
   TMVA_Tree->Branch("Muon2_station3_pullDxDz", &Muon2_station_vars[2][5]);
   TMVA_Tree->Branch("Muon2_station3_pullDyDz", &Muon2_station_vars[2][6]);

   TMVA_Tree->Branch("Muon2_station4_status", &Muon2_station_vars[3][0]);
   TMVA_Tree->Branch("Muon2_station4_TrackX", &Muon2_station_vars[3][1]);
   TMVA_Tree->Branch("Muon2_station4_TrackY", &Muon2_station_vars[3][2]);
   TMVA_Tree->Branch("Muon2_station4_pullX", &Muon2_station_vars[3][3]);
   TMVA_Tree->Branch("Muon2_station4_pullY", &Muon2_station_vars[3][4]);
   TMVA_Tree->Branch("Muon2_station4_pullDxDz", &Muon2_station_vars[3][5]);
   TMVA_Tree->Branch("Muon2_station4_pullDyDz", &Muon2_station_vars[3][6]);

   TMVA_Tree->Branch("Muon2_station5_status", &Muon2_station_vars[4][0]);
   TMVA_Tree->Branch("Muon2_station5_TrackX", &Muon2_station_vars[4][1]);
   TMVA_Tree->Branch("Muon2_station5_TrackY", &Muon2_station_vars[4][2]);
   TMVA_Tree->Branch("Muon2_station5_pullX", &Muon2_station_vars[4][3]);
   TMVA_Tree->Branch("Muon2_station5_pullY", &Muon2_station_vars[4][4]);
   TMVA_Tree->Branch("Muon2_station5_pullDxDz", &Muon2_station_vars[4][5]);
   TMVA_Tree->Branch("Muon2_station5_pullDyDz", &Muon2_station_vars[4][6]);

   TMVA_Tree->Branch("Muon2_station6_status", &Muon2_station_vars[5][0]);
   TMVA_Tree->Branch("Muon2_station6_TrackX", &Muon2_station_vars[5][1]);
   TMVA_Tree->Branch("Muon2_station6_TrackY", &Muon2_station_vars[5][2]);
   TMVA_Tree->Branch("Muon2_station6_pullX", &Muon2_station_vars[5][3]);
   TMVA_Tree->Branch("Muon2_station6_pullY", &Muon2_station_vars[5][4]);
   TMVA_Tree->Branch("Muon2_station6_pullDxDz", &Muon2_station_vars[5][5]);
   TMVA_Tree->Branch("Muon2_station6_pullDyDz", &Muon2_station_vars[5][6]);

   TMVA_Tree->Branch("Muon2_station7_status", &Muon2_station_vars[6][0]);
   TMVA_Tree->Branch("Muon2_station7_TrackX", &Muon2_station_vars[6][1]);
   TMVA_Tree->Branch("Muon2_station7_TrackY", &Muon2_station_vars[6][2]);
   TMVA_Tree->Branch("Muon2_station7_pullX", &Muon2_station_vars[6][3]);
   TMVA_Tree->Branch("Muon2_station7_pullY", &Muon2_station_vars[6][4]);
   TMVA_Tree->Branch("Muon2_station7_pullDxDz", &Muon2_station_vars[6][5]);
   TMVA_Tree->Branch("Muon2_station7_pullDyDz", &Muon2_station_vars[6][6]);

   TMVA_Tree->Branch("Muon2_station8_status", &Muon2_station_vars[7][0]);
   TMVA_Tree->Branch("Muon2_station8_TrackX", &Muon2_station_vars[7][1]);
   TMVA_Tree->Branch("Muon2_station8_TrackY", &Muon2_station_vars[7][2]);
   TMVA_Tree->Branch("Muon2_station8_pullX", &Muon2_station_vars[7][3]);
   TMVA_Tree->Branch("Muon2_station8_pullY", &Muon2_station_vars[7][4]);
   TMVA_Tree->Branch("Muon2_station8_pullDxDz", &Muon2_station_vars[7][5]);
   TMVA_Tree->Branch("Muon2_station8_pullDyDz", &Muon2_station_vars[7][6]);

   // Muon station information (Muon3)
   TMVA_Tree->Branch("Muon3_station1_status", &Muon3_station_vars[0][0]);
   TMVA_Tree->Branch("Muon3_station1_TrackX", &Muon3_station_vars[0][1]);
   TMVA_Tree->Branch("Muon3_station1_TrackY", &Muon3_station_vars[0][2]);
   TMVA_Tree->Branch("Muon3_station1_pullX", &Muon3_station_vars[0][3]);
   TMVA_Tree->Branch("Muon3_station1_pullY", &Muon3_station_vars[0][4]);
   TMVA_Tree->Branch("Muon3_station1_pullDxDz", &Muon3_station_vars[0][5]);
   TMVA_Tree->Branch("Muon3_station1_pullDyDz", &Muon3_station_vars[0][6]);

   TMVA_Tree->Branch("Muon3_station2_status", &Muon3_station_vars[1][0]);
   TMVA_Tree->Branch("Muon3_station2_TrackX", &Muon3_station_vars[1][1]);
   TMVA_Tree->Branch("Muon3_station2_TrackY", &Muon3_station_vars[1][2]);
   TMVA_Tree->Branch("Muon3_station2_pullX", &Muon3_station_vars[1][3]);
   TMVA_Tree->Branch("Muon3_station2_pullY", &Muon3_station_vars[1][4]);
   TMVA_Tree->Branch("Muon3_station2_pullDxDz", &Muon3_station_vars[1][5]);
   TMVA_Tree->Branch("Muon3_station2_pullDyDz", &Muon3_station_vars[1][6]);

   TMVA_Tree->Branch("Muon3_station3_status", &Muon3_station_vars[2][0]);
   TMVA_Tree->Branch("Muon3_station3_TrackX", &Muon3_station_vars[2][1]);
   TMVA_Tree->Branch("Muon3_station3_TrackY", &Muon3_station_vars[2][2]);
   TMVA_Tree->Branch("Muon3_station3_pullX", &Muon3_station_vars[2][3]);
   TMVA_Tree->Branch("Muon3_station3_pullY", &Muon3_station_vars[2][4]);
   TMVA_Tree->Branch("Muon3_station3_pullDxDz", &Muon3_station_vars[2][5]);
   TMVA_Tree->Branch("Muon3_station3_pullDyDz", &Muon3_station_vars[2][6]);

   TMVA_Tree->Branch("Muon3_station4_status", &Muon3_station_vars[3][0]);
   TMVA_Tree->Branch("Muon3_station4_TrackX", &Muon3_station_vars[3][1]);
   TMVA_Tree->Branch("Muon3_station4_TrackY", &Muon3_station_vars[3][2]);
   TMVA_Tree->Branch("Muon3_station4_pullX", &Muon3_station_vars[3][3]);
   TMVA_Tree->Branch("Muon3_station4_pullY", &Muon3_station_vars[3][4]);
   TMVA_Tree->Branch("Muon3_station4_pullDxDz", &Muon3_station_vars[3][5]);
   TMVA_Tree->Branch("Muon3_station4_pullDyDz", &Muon3_station_vars[3][6]);

   TMVA_Tree->Branch("Muon3_station5_status", &Muon3_station_vars[4][0]);
   TMVA_Tree->Branch("Muon3_station5_TrackX", &Muon3_station_vars[4][1]);
   TMVA_Tree->Branch("Muon3_station5_TrackY", &Muon3_station_vars[4][2]);
   TMVA_Tree->Branch("Muon3_station5_pullX", &Muon3_station_vars[4][3]);
   TMVA_Tree->Branch("Muon3_station5_pullY", &Muon3_station_vars[4][4]);
   TMVA_Tree->Branch("Muon3_station5_pullDxDz", &Muon3_station_vars[4][5]);
   TMVA_Tree->Branch("Muon3_station5_pullDyDz", &Muon3_station_vars[4][6]);

   TMVA_Tree->Branch("Muon3_station6_status", &Muon3_station_vars[5][0]);
   TMVA_Tree->Branch("Muon3_station6_TrackX", &Muon3_station_vars[5][1]);
   TMVA_Tree->Branch("Muon3_station6_TrackY", &Muon3_station_vars[5][2]);
   TMVA_Tree->Branch("Muon3_station6_pullX", &Muon3_station_vars[5][3]);
   TMVA_Tree->Branch("Muon3_station6_pullY", &Muon3_station_vars[5][4]);
   TMVA_Tree->Branch("Muon3_station6_pullDxDz", &Muon3_station_vars[5][5]);
   TMVA_Tree->Branch("Muon3_station6_pullDyDz", &Muon3_station_vars[5][6]);

   TMVA_Tree->Branch("Muon3_station7_status", &Muon3_station_vars[6][0]);
   TMVA_Tree->Branch("Muon3_station7_TrackX", &Muon3_station_vars[6][1]);
   TMVA_Tree->Branch("Muon3_station7_TrackY", &Muon3_station_vars[6][2]);
   TMVA_Tree->Branch("Muon3_station7_pullX", &Muon3_station_vars[6][3]);
   TMVA_Tree->Branch("Muon3_station7_pullY", &Muon3_station_vars[6][4]);
   TMVA_Tree->Branch("Muon3_station7_pullDxDz", &Muon3_station_vars[6][5]);
   TMVA_Tree->Branch("Muon3_station7_pullDyDz", &Muon3_station_vars[6][6]);

   TMVA_Tree->Branch("Muon3_station8_status", &Muon3_station_vars[7][0]);
   TMVA_Tree->Branch("Muon3_station8_TrackX", &Muon3_station_vars[7][1]);
   TMVA_Tree->Branch("Muon3_station8_TrackY", &Muon3_station_vars[7][2]);
   TMVA_Tree->Branch("Muon3_station8_pullX", &Muon3_station_vars[7][3]);
   TMVA_Tree->Branch("Muon3_station8_pullY", &Muon3_station_vars[7][4]);
   TMVA_Tree->Branch("Muon3_station8_pullDxDz", &Muon3_station_vars[7][5]);
   TMVA_Tree->Branch("Muon3_station8_pullDyDz", &Muon3_station_vars[7][6]);

   // Muon variables
   TMVA_Tree->Branch("var_Muon1_Pt", &var_Muon1_Pt);
   TMVA_Tree->Branch("var_Muon1_Eta", &var_Muon1_Eta);
   TMVA_Tree->Branch("var_Muon1_Phi", &var_Muon1_Phi);
   TMVA_Tree->Branch("var_Muon2_Pt", &var_Muon2_Pt);
   TMVA_Tree->Branch("var_Muon2_Eta", &var_Muon2_Eta);
   TMVA_Tree->Branch("var_Muon2_Phi", &var_Muon2_Phi);
   TMVA_Tree->Branch("var_Muon3_Pt", &var_Muon3_Pt);
   TMVA_Tree->Branch("var_Muon3_Eta", &var_Muon3_Eta);
   TMVA_Tree->Branch("var_Muon3_Phi", &var_Muon3_Phi);

   TMVA_Tree->Branch("Muon1_vx",&var_Muon1_vx);
   TMVA_Tree->Branch("Muon1_vy",&var_Muon1_vy);
   TMVA_Tree->Branch("Muon1_vz",&var_Muon1_vz);
   TMVA_Tree->Branch("Muon1_IsGlobalMuon",&var_Muon1_IsGlobalMuon);
   TMVA_Tree->Branch("Muon1_IsStandAloneMuon",&var_Muon1_IsStandAloneMuon);
   TMVA_Tree->Branch("Muon1_IsTrackerMuon",&var_Muon1_IsTrackerMuon);
   TMVA_Tree->Branch("Muon1_IsCaloMuon",&var_Muon1_IsCaloMuon);
   TMVA_Tree->Branch("Muon1_IsIsolationValid",&var_Muon1_IsIsolationValid);
   TMVA_Tree->Branch("Muon1_IsQualityValid",&var_Muon1_IsQualityValid);
   TMVA_Tree->Branch("Muon1_IsTimeValid",&var_Muon1_IsTimeValid);
   TMVA_Tree->Branch("Muon1_IsPFMuon",&var_Muon1_IsPFMuon);
   TMVA_Tree->Branch("Muon1_IsRPCMuon",&var_Muon1_IsRPCMuon);
   TMVA_Tree->Branch("Muon1_emEt03",&var_Muon1_emEt03);
   TMVA_Tree->Branch("Muon1_emVetoEt03",&var_Muon1_emVetoEt03);
   TMVA_Tree->Branch("Muon1_hadEt03",&var_Muon1_hadEt03);
   TMVA_Tree->Branch("Muon1_hadVetoEt03",&var_Muon1_hadVetoEt03);
   TMVA_Tree->Branch("Muon1_nJets03",&var_Muon1_nJets03);
   TMVA_Tree->Branch("Muon1_nTracks03",&var_Muon1_nTracks03);
   TMVA_Tree->Branch("Muon1_StandardSelection",&var_Muon1_StandardSelection);
   TMVA_Tree->Branch("Muon1_trackerVetoPt03",&var_Muon1_trackerVetoPt03);
   TMVA_Tree->Branch("Muon1_sumChargedHadronPt03",&var_Muon1_sumChargedHadronPt03);
   TMVA_Tree->Branch("Muon1_sumChargedParticlePt03",&var_Muon1_sumChargedParticlePt03);
   TMVA_Tree->Branch("Muon1_sumNeutralHadronEt03",&var_Muon1_sumNeutralHadronEt03);
   TMVA_Tree->Branch("Muon1_sumNeutralHadronEtHighThreshold03",&var_Muon1_sumNeutralHadronEtHighThreshold03);
   TMVA_Tree->Branch("Muon1_sumPhotonEt03",&var_Muon1_sumPhotonEt03);
   TMVA_Tree->Branch("Muon1_sumPhotonEtHighThreshold03",&var_Muon1_sumPhotonEtHighThreshold03);
   TMVA_Tree->Branch("Muon1_sumPUPt03",&var_Muon1_sumPUPt03);
   TMVA_Tree->Branch("Muon1_numberOfChambers",&var_Muon1_numberOfChambers);
   TMVA_Tree->Branch("Muon1_Track_idx",&var_Muon1_Track_idx);
   TMVA_Tree->Branch("Muon1_combinedQuality_updatedSta",&var_Muon1_combinedQuality_updatedSta);
   TMVA_Tree->Branch("Muon1_combinedQuality_trkKink",&var_Muon1_combinedQuality_trkKink);
   TMVA_Tree->Branch("Muon1_combinedQuality_glbKink",&var_Muon1_combinedQuality_glbKink);
   TMVA_Tree->Branch("Muon1_combinedQuality_trkRelChi2",&var_Muon1_combinedQuality_trkRelChi2);
   TMVA_Tree->Branch("Muon1_combinedQuality_staRelChi2",&var_Muon1_combinedQuality_staRelChi2);
   TMVA_Tree->Branch("Muon1_combinedQuality_chi2LocalPosition",&var_Muon1_combinedQuality_chi2LocalPosition);
   TMVA_Tree->Branch("Muon1_combinedQuality_chi2LocalMomentum",&var_Muon1_combinedQuality_chi2LocalMomentum);
   TMVA_Tree->Branch("Muon1_combinedQuality_localDistance",&var_Muon1_combinedQuality_localDistance);
   TMVA_Tree->Branch("Muon1_combinedQuality_globalDeltaEtaPhi",&var_Muon1_combinedQuality_globalDeltaEtaPhi);
   TMVA_Tree->Branch("Muon1_combinedQuality_tightMatch",&var_Muon1_combinedQuality_tightMatch);
   TMVA_Tree->Branch("Muon1_combinedQuality_glbTrackProbability",&var_Muon1_combinedQuality_glbTrackProbability);
   TMVA_Tree->Branch("Muon1_prod_inner_outer_charge",&var_Muon1_prod_inner_outer_charge);
   TMVA_Tree->Branch("Muon1_innerTrack_quality",&var_Muon1_innerTrack_quality);
   TMVA_Tree->Branch("Muon1_ptErrOverPt",&var_Muon1_ptErrOverPt);
   TMVA_Tree->Branch("Muon1_calEnergy_hadS9",&var_Muon1_calEnergy_hadS9);
   TMVA_Tree->Branch("Muon1_calEnergy_had",&var_Muon1_calEnergy_had);
   TMVA_Tree->Branch("Muon1_calEnergy_emS25",&var_Muon1_calEnergy_emS25);
   TMVA_Tree->Branch("Muon1_calEnergy_emS9",&var_Muon1_calEnergy_emS9);
   TMVA_Tree->Branch("Muon1_calEnergy_em",&var_Muon1_calEnergy_em);
   TMVA_Tree->Branch("Muon1_charge",&var_Muon1_charge);
   TMVA_Tree->Branch("Muon1_trackCharge",&var_Muon1_trackCharge);
   TMVA_Tree->Branch("Muon1_hitPattern_pixelLayerwithMeas",&var_Muon1_hitPattern_pixelLayerwithMeas);
   TMVA_Tree->Branch("Muon1_numberOfMatchedStations",&var_Muon1_numberOfMatchedStations);
   TMVA_Tree->Branch("Muon1_normChi2",&var_Muon1_normChi2);
   TMVA_Tree->Branch("Muon1_hitPattern_numberOfValidMuonHits",&var_Muon1_hitPattern_numberOfValidMuonHits);
   TMVA_Tree->Branch("Muon1_innerTrack_numberofValidHits",&var_Muon1_innerTrack_numberofValidHits);
   TMVA_Tree->Branch("Muon1_numberofValidPixelHits",&var_Muon1_numberofValidPixelHits);
   TMVA_Tree->Branch("Muon1_numberOfMatches",&var_Muon1_numberOfMatches);
   TMVA_Tree->Branch("Muon1_trackerLayersWithMeasurement",&var_Muon1_trackerLayersWithMeasurement);
   TMVA_Tree->Branch("Muon1_segmentCompatibility",&var_Muon1_segmentCompatibility);
   TMVA_Tree->Branch("Muon1_caloCompatibility",&var_Muon1_caloCompatibility);
   TMVA_Tree->Branch("Muon1_innerTrack_validFraction",&var_Muon1_innerTrack_validFraction);
   TMVA_Tree->Branch("Muon1_innerTrack_pixelLayersWithMeasurement",&var_Muon1_innerTrack_pixelLayersWithMeasurement);
   TMVA_Tree->Branch("Muon1_innerTrack_numberOfValidTrackerHits",&var_Muon1_innerTrack_numberOfValidTrackerHits);
   TMVA_Tree->Branch("Muon1_innerTrack_numberOfLostTrackerHits",&var_Muon1_innerTrack_numberOfLostTrackerHits);
   TMVA_Tree->Branch("Muon1_innerTrack_numberOfLostTrackerInnerHits",&var_Muon1_innerTrack_numberOfLostTrackerInnerHits);
   TMVA_Tree->Branch("Muon1_innerTrack_numberOfLostTrackerOuterHits",&var_Muon1_innerTrack_numberOfLostTrackerOuterHits);
   TMVA_Tree->Branch("Muon1_innerTrack_normalizedChi2",&var_Muon1_innerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon1_vmuonhitcomb_reco",&var_Muon1_vmuonhitcomb_reco);
   TMVA_Tree->Branch("Muon1_rpchits_reco",&var_Muon1_rpchits_reco);
   TMVA_Tree->Branch("Muon1_outerTrack_normalizedChi2",&var_Muon1_outerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon1_outerTrack_muonStationsWithValidHits",&var_Muon1_outerTrack_muonStationsWithValidHits);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TM2DCompatibility",&var_Muon1_isGoodMuon_TM2DCompatibility);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TrackerMuonArbitrated",&var_Muon1_isGoodMuon_TrackerMuonArbitrated);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMOneStationTight",&var_Muon1_isGoodMuon_TMOneStationTight);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMOneStationAngTight",&var_Muon1_isGoodMuon_TMOneStationAngTight);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMLastStationTight",&var_Muon1_isGoodMuon_TMLastStationTight);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMLastStationAngTight",&var_Muon1_isGoodMuon_TMLastStationAngTight);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMLastStationOptimizedLowPtTight",&var_Muon1_isGoodMuon_TMLastStationOptimizedLowPtTight);
   TMVA_Tree->Branch("Muon1_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight",&var_Muon1_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);

   TMVA_Tree->Branch("Muon2_vx",&var_Muon2_vx);
   TMVA_Tree->Branch("Muon2_vy",&var_Muon2_vy);
   TMVA_Tree->Branch("Muon2_vz",&var_Muon2_vz);
   TMVA_Tree->Branch("Muon2_IsGlobalMuon",&var_Muon2_IsGlobalMuon);
   TMVA_Tree->Branch("Muon2_IsStandAloneMuon",&var_Muon2_IsStandAloneMuon);
   TMVA_Tree->Branch("Muon2_IsTrackerMuon",&var_Muon2_IsTrackerMuon);
   TMVA_Tree->Branch("Muon2_IsCaloMuon",&var_Muon2_IsCaloMuon);
   TMVA_Tree->Branch("Muon2_IsIsolationValid",&var_Muon2_IsIsolationValid);
   TMVA_Tree->Branch("Muon2_IsQualityValid",&var_Muon2_IsQualityValid);
   TMVA_Tree->Branch("Muon2_IsTimeValid",&var_Muon2_IsTimeValid);
   TMVA_Tree->Branch("Muon2_IsPFMuon",&var_Muon2_IsPFMuon);
   TMVA_Tree->Branch("Muon2_IsRPCMuon",&var_Muon2_IsRPCMuon);
   TMVA_Tree->Branch("Muon2_emEt03",&var_Muon2_emEt03);
   TMVA_Tree->Branch("Muon2_emVetoEt03",&var_Muon2_emVetoEt03);
   TMVA_Tree->Branch("Muon2_hadEt03",&var_Muon2_hadEt03);
   TMVA_Tree->Branch("Muon2_hadVetoEt03",&var_Muon2_hadVetoEt03);
   TMVA_Tree->Branch("Muon2_nJets03",&var_Muon2_nJets03);
   TMVA_Tree->Branch("Muon2_nTracks03",&var_Muon2_nTracks03);
   TMVA_Tree->Branch("Muon2_StandardSelection",&var_Muon2_StandardSelection);
   TMVA_Tree->Branch("Muon2_trackerVetoPt03",&var_Muon2_trackerVetoPt03);
   TMVA_Tree->Branch("Muon2_sumChargedHadronPt03",&var_Muon2_sumChargedHadronPt03);
   TMVA_Tree->Branch("Muon2_sumChargedParticlePt03",&var_Muon2_sumChargedParticlePt03);
   TMVA_Tree->Branch("Muon2_sumNeutralHadronEt03",&var_Muon2_sumNeutralHadronEt03);
   TMVA_Tree->Branch("Muon2_sumNeutralHadronEtHighThreshold03",&var_Muon2_sumNeutralHadronEtHighThreshold03);
   TMVA_Tree->Branch("Muon2_sumPhotonEt03",&var_Muon2_sumPhotonEt03);
   TMVA_Tree->Branch("Muon2_sumPhotonEtHighThreshold03",&var_Muon2_sumPhotonEtHighThreshold03);
   TMVA_Tree->Branch("Muon2_sumPUPt03",&var_Muon2_sumPUPt03);
   TMVA_Tree->Branch("Muon2_numberOfChambers",&var_Muon2_numberOfChambers);
   TMVA_Tree->Branch("Muon2_Track_idx",&var_Muon2_Track_idx);
   TMVA_Tree->Branch("Muon2_combinedQuality_updatedSta",&var_Muon2_combinedQuality_updatedSta);
   TMVA_Tree->Branch("Muon2_combinedQuality_trkKink",&var_Muon2_combinedQuality_trkKink);
   TMVA_Tree->Branch("Muon2_combinedQuality_glbKink",&var_Muon2_combinedQuality_glbKink);
   TMVA_Tree->Branch("Muon2_combinedQuality_trkRelChi2",&var_Muon2_combinedQuality_trkRelChi2);
   TMVA_Tree->Branch("Muon2_combinedQuality_staRelChi2",&var_Muon2_combinedQuality_staRelChi2);
   TMVA_Tree->Branch("Muon2_combinedQuality_chi2LocalPosition",&var_Muon2_combinedQuality_chi2LocalPosition);
   TMVA_Tree->Branch("Muon2_combinedQuality_chi2LocalMomentum",&var_Muon2_combinedQuality_chi2LocalMomentum);
   TMVA_Tree->Branch("Muon2_combinedQuality_localDistance",&var_Muon2_combinedQuality_localDistance);
   TMVA_Tree->Branch("Muon2_combinedQuality_globalDeltaEtaPhi",&var_Muon2_combinedQuality_globalDeltaEtaPhi);
   TMVA_Tree->Branch("Muon2_combinedQuality_tightMatch",&var_Muon2_combinedQuality_tightMatch);
   TMVA_Tree->Branch("Muon2_combinedQuality_glbTrackProbability",&var_Muon2_combinedQuality_glbTrackProbability);
   TMVA_Tree->Branch("Muon2_prod_inner_outer_charge",&var_Muon2_prod_inner_outer_charge);
   TMVA_Tree->Branch("Muon2_innerTrack_quality",&var_Muon2_innerTrack_quality);
   TMVA_Tree->Branch("Muon2_ptErrOverPt",&var_Muon2_ptErrOverPt);
   TMVA_Tree->Branch("Muon2_calEnergy_hadS9",&var_Muon2_calEnergy_hadS9);
   TMVA_Tree->Branch("Muon2_calEnergy_had",&var_Muon2_calEnergy_had);
   TMVA_Tree->Branch("Muon2_calEnergy_emS25",&var_Muon2_calEnergy_emS25);
   TMVA_Tree->Branch("Muon2_calEnergy_emS9",&var_Muon2_calEnergy_emS9);
   TMVA_Tree->Branch("Muon2_calEnergy_em",&var_Muon2_calEnergy_em);
   TMVA_Tree->Branch("Muon2_charge",&var_Muon2_charge);
   TMVA_Tree->Branch("Muon2_trackCharge",&var_Muon2_trackCharge);
   TMVA_Tree->Branch("Muon2_hitPattern_pixelLayerwithMeas",&var_Muon2_hitPattern_pixelLayerwithMeas);
   TMVA_Tree->Branch("Muon2_numberOfMatchedStations",&var_Muon2_numberOfMatchedStations);
   TMVA_Tree->Branch("Muon2_normChi2",&var_Muon2_normChi2);
   TMVA_Tree->Branch("Muon2_hitPattern_numberOfValidMuonHits",&var_Muon2_hitPattern_numberOfValidMuonHits);
   TMVA_Tree->Branch("Muon2_innerTrack_numberofValidHits",&var_Muon2_innerTrack_numberofValidHits);
   TMVA_Tree->Branch("Muon2_numberofValidPixelHits",&var_Muon2_numberofValidPixelHits);
   TMVA_Tree->Branch("Muon2_numberOfMatches",&var_Muon2_numberOfMatches);
   TMVA_Tree->Branch("Muon2_trackerLayersWithMeasurement",&var_Muon2_trackerLayersWithMeasurement);
   TMVA_Tree->Branch("Muon2_segmentCompatibility",&var_Muon2_segmentCompatibility);
   TMVA_Tree->Branch("Muon2_caloCompatibility",&var_Muon2_caloCompatibility);
   TMVA_Tree->Branch("Muon2_innerTrack_validFraction",&var_Muon2_innerTrack_validFraction);
   TMVA_Tree->Branch("Muon2_innerTrack_pixelLayersWithMeasurement",&var_Muon2_innerTrack_pixelLayersWithMeasurement);
   TMVA_Tree->Branch("Muon2_innerTrack_numberOfValidTrackerHits",&var_Muon2_innerTrack_numberOfValidTrackerHits);
   TMVA_Tree->Branch("Muon2_innerTrack_numberOfLostTrackerHits",&var_Muon2_innerTrack_numberOfLostTrackerHits);
   TMVA_Tree->Branch("Muon2_innerTrack_numberOfLostTrackerInnerHits",&var_Muon2_innerTrack_numberOfLostTrackerInnerHits);
   TMVA_Tree->Branch("Muon2_innerTrack_numberOfLostTrackerOuterHits",&var_Muon2_innerTrack_numberOfLostTrackerOuterHits);
   TMVA_Tree->Branch("Muon2_innerTrack_normalizedChi2",&var_Muon2_innerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon2_vmuonhitcomb_reco",&var_Muon2_vmuonhitcomb_reco);
   TMVA_Tree->Branch("Muon2_rpchits_reco",&var_Muon2_rpchits_reco);
   TMVA_Tree->Branch("Muon2_outerTrack_normalizedChi2",&var_Muon2_outerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon2_outerTrack_muonStationsWithValidHits",&var_Muon2_outerTrack_muonStationsWithValidHits);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TM2DCompatibility",&var_Muon2_isGoodMuon_TM2DCompatibility);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TrackerMuonArbitrated",&var_Muon2_isGoodMuon_TrackerMuonArbitrated);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMOneStationTight",&var_Muon2_isGoodMuon_TMOneStationTight);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMOneStationAngTight",&var_Muon2_isGoodMuon_TMOneStationAngTight);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMLastStationTight",&var_Muon2_isGoodMuon_TMLastStationTight);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMLastStationAngTight",&var_Muon2_isGoodMuon_TMLastStationAngTight);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMLastStationOptimizedLowPtTight",&var_Muon2_isGoodMuon_TMLastStationOptimizedLowPtTight);
   TMVA_Tree->Branch("Muon2_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight",&var_Muon2_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);

   TMVA_Tree->Branch("Muon3_vx",&var_Muon3_vx);
   TMVA_Tree->Branch("Muon3_vy",&var_Muon3_vy);
   TMVA_Tree->Branch("Muon3_vz",&var_Muon3_vz);
   TMVA_Tree->Branch("Muon3_IsGlobalMuon",&var_Muon3_IsGlobalMuon);
   TMVA_Tree->Branch("Muon3_IsStandAloneMuon",&var_Muon3_IsStandAloneMuon);
   TMVA_Tree->Branch("Muon3_IsTrackerMuon",&var_Muon3_IsTrackerMuon);
   TMVA_Tree->Branch("Muon3_IsCaloMuon",&var_Muon3_IsCaloMuon);
   TMVA_Tree->Branch("Muon3_IsIsolationValid",&var_Muon3_IsIsolationValid);
   TMVA_Tree->Branch("Muon3_IsQualityValid",&var_Muon3_IsQualityValid);
   TMVA_Tree->Branch("Muon3_IsTimeValid",&var_Muon3_IsTimeValid);
   TMVA_Tree->Branch("Muon3_IsPFMuon",&var_Muon3_IsPFMuon);
   TMVA_Tree->Branch("Muon3_IsRPCMuon",&var_Muon3_IsRPCMuon);
   TMVA_Tree->Branch("Muon3_emEt03",&var_Muon3_emEt03);
   TMVA_Tree->Branch("Muon3_emVetoEt03",&var_Muon3_emVetoEt03);
   TMVA_Tree->Branch("Muon3_hadEt03",&var_Muon3_hadEt03);
   TMVA_Tree->Branch("Muon3_hadVetoEt03",&var_Muon3_hadVetoEt03);
   TMVA_Tree->Branch("Muon3_nJets03",&var_Muon3_nJets03);
   TMVA_Tree->Branch("Muon3_nTracks03",&var_Muon3_nTracks03);
   TMVA_Tree->Branch("Muon3_StandardSelection",&var_Muon3_StandardSelection);
   TMVA_Tree->Branch("Muon3_trackerVetoPt03",&var_Muon3_trackerVetoPt03);
   TMVA_Tree->Branch("Muon3_sumChargedHadronPt03",&var_Muon3_sumChargedHadronPt03);
   TMVA_Tree->Branch("Muon3_sumChargedParticlePt03",&var_Muon3_sumChargedParticlePt03);
   TMVA_Tree->Branch("Muon3_sumNeutralHadronEt03",&var_Muon3_sumNeutralHadronEt03);
   TMVA_Tree->Branch("Muon3_sumNeutralHadronEtHighThreshold03",&var_Muon3_sumNeutralHadronEtHighThreshold03);
   TMVA_Tree->Branch("Muon3_sumPhotonEt03",&var_Muon3_sumPhotonEt03);
   TMVA_Tree->Branch("Muon3_sumPhotonEtHighThreshold03",&var_Muon3_sumPhotonEtHighThreshold03);
   TMVA_Tree->Branch("Muon3_sumPUPt03",&var_Muon3_sumPUPt03);
   TMVA_Tree->Branch("Muon3_numberOfChambers",&var_Muon3_numberOfChambers);
   TMVA_Tree->Branch("Muon3_Track_idx",&var_Muon3_Track_idx);
   TMVA_Tree->Branch("Muon3_combinedQuality_updatedSta",&var_Muon3_combinedQuality_updatedSta);
   TMVA_Tree->Branch("Muon3_combinedQuality_trkKink",&var_Muon3_combinedQuality_trkKink);
   TMVA_Tree->Branch("Muon3_combinedQuality_glbKink",&var_Muon3_combinedQuality_glbKink);
   TMVA_Tree->Branch("Muon3_combinedQuality_trkRelChi2",&var_Muon3_combinedQuality_trkRelChi2);
   TMVA_Tree->Branch("Muon3_combinedQuality_staRelChi2",&var_Muon3_combinedQuality_staRelChi2);
   TMVA_Tree->Branch("Muon3_combinedQuality_chi2LocalPosition",&var_Muon3_combinedQuality_chi2LocalPosition);
   TMVA_Tree->Branch("Muon3_combinedQuality_chi2LocalMomentum",&var_Muon3_combinedQuality_chi2LocalMomentum);
   TMVA_Tree->Branch("Muon3_combinedQuality_localDistance",&var_Muon3_combinedQuality_localDistance);
   TMVA_Tree->Branch("Muon3_combinedQuality_globalDeltaEtaPhi",&var_Muon3_combinedQuality_globalDeltaEtaPhi);
   TMVA_Tree->Branch("Muon3_combinedQuality_tightMatch",&var_Muon3_combinedQuality_tightMatch);
   TMVA_Tree->Branch("Muon3_combinedQuality_glbTrackProbability",&var_Muon3_combinedQuality_glbTrackProbability);
   TMVA_Tree->Branch("Muon3_prod_inner_outer_charge",&var_Muon3_prod_inner_outer_charge);
   TMVA_Tree->Branch("Muon3_innerTrack_quality",&var_Muon3_innerTrack_quality);
   TMVA_Tree->Branch("Muon3_ptErrOverPt",&var_Muon3_ptErrOverPt);
   TMVA_Tree->Branch("Muon3_calEnergy_hadS9",&var_Muon3_calEnergy_hadS9);
   TMVA_Tree->Branch("Muon3_calEnergy_had",&var_Muon3_calEnergy_had);
   TMVA_Tree->Branch("Muon3_calEnergy_emS25",&var_Muon3_calEnergy_emS25);
   TMVA_Tree->Branch("Muon3_calEnergy_emS9",&var_Muon3_calEnergy_emS9);
   TMVA_Tree->Branch("Muon3_calEnergy_em",&var_Muon3_calEnergy_em);
   TMVA_Tree->Branch("Muon3_charge",&var_Muon3_charge);
   TMVA_Tree->Branch("Muon3_trackCharge",&var_Muon3_trackCharge);
   TMVA_Tree->Branch("Muon3_hitPattern_pixelLayerwithMeas",&var_Muon3_hitPattern_pixelLayerwithMeas);
   TMVA_Tree->Branch("Muon3_numberOfMatchedStations",&var_Muon3_numberOfMatchedStations);
   TMVA_Tree->Branch("Muon3_normChi2",&var_Muon3_normChi2);
   TMVA_Tree->Branch("Muon3_hitPattern_numberOfValidMuonHits",&var_Muon3_hitPattern_numberOfValidMuonHits);
   TMVA_Tree->Branch("Muon3_innerTrack_numberofValidHits",&var_Muon3_innerTrack_numberofValidHits);
   TMVA_Tree->Branch("Muon3_numberofValidPixelHits",&var_Muon3_numberofValidPixelHits);
   TMVA_Tree->Branch("Muon3_numberOfMatches",&var_Muon3_numberOfMatches);
   TMVA_Tree->Branch("Muon3_trackerLayersWithMeasurement",&var_Muon3_trackerLayersWithMeasurement);
   TMVA_Tree->Branch("Muon3_segmentCompatibility",&var_Muon3_segmentCompatibility);
   TMVA_Tree->Branch("Muon3_caloCompatibility",&var_Muon3_caloCompatibility);
   TMVA_Tree->Branch("Muon3_innerTrack_validFraction",&var_Muon3_innerTrack_validFraction);
   TMVA_Tree->Branch("Muon3_innerTrack_pixelLayersWithMeasurement",&var_Muon3_innerTrack_pixelLayersWithMeasurement);
   TMVA_Tree->Branch("Muon3_innerTrack_numberOfValidTrackerHits",&var_Muon3_innerTrack_numberOfValidTrackerHits);
   TMVA_Tree->Branch("Muon3_innerTrack_numberOfLostTrackerHits",&var_Muon3_innerTrack_numberOfLostTrackerHits);
   TMVA_Tree->Branch("Muon3_innerTrack_numberOfLostTrackerInnerHits",&var_Muon3_innerTrack_numberOfLostTrackerInnerHits);
   TMVA_Tree->Branch("Muon3_innerTrack_numberOfLostTrackerOuterHits",&var_Muon3_innerTrack_numberOfLostTrackerOuterHits);
   TMVA_Tree->Branch("Muon3_innerTrack_normalizedChi2",&var_Muon3_innerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon3_vmuonhitcomb_reco",&var_Muon3_vmuonhitcomb_reco);
   TMVA_Tree->Branch("Muon3_rpchits_reco",&var_Muon3_rpchits_reco);
   TMVA_Tree->Branch("Muon3_outerTrack_normalizedChi2",&var_Muon3_outerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon3_outerTrack_muonStationsWithValidHits",&var_Muon3_outerTrack_muonStationsWithValidHits);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TM2DCompatibility",&var_Muon3_isGoodMuon_TM2DCompatibility);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TrackerMuonArbitrated",&var_Muon3_isGoodMuon_TrackerMuonArbitrated);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMOneStationTight",&var_Muon3_isGoodMuon_TMOneStationTight);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMOneStationAngTight",&var_Muon3_isGoodMuon_TMOneStationAngTight);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMLastStationTight",&var_Muon3_isGoodMuon_TMLastStationTight);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMLastStationAngTight",&var_Muon3_isGoodMuon_TMLastStationAngTight);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMLastStationOptimizedLowPtTight",&var_Muon3_isGoodMuon_TMLastStationOptimizedLowPtTight);
   TMVA_Tree->Branch("Muon3_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight",&var_Muon3_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);

   // Time variables
   TMVA_Tree->Branch("var_Muon1_timeAtIpInOutErr", &var_Muon1_timeAtIpInOutErr);
   TMVA_Tree->Branch("var_Muon1_timeAtIpInOut", &var_Muon1_timeAtIpInOut);
   TMVA_Tree->Branch("var_Muon2_timeAtIpInOutErr", &var_Muon2_timeAtIpInOutErr);
   TMVA_Tree->Branch("var_Muon2_timeAtIpInOut", &var_Muon2_timeAtIpInOut);
   TMVA_Tree->Branch("var_Muon3_timeAtIpInOutErr", &var_Muon3_timeAtIpInOutErr);
   TMVA_Tree->Branch("var_Muon3_timeAtIpInOut", &var_Muon3_timeAtIpInOut);

   // Isolation variables
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

   TMVA_Tree->Branch("var_NtracksClose",&var_NtracksClose);
   TMVA_Tree->Branch("var_MindcaTrackSV",&var_MindcaTrackSV);
   TMVA_Tree->Branch("var_Mu1TrackMass",&var_Mu1TrackMass);
   TMVA_Tree->Branch("var_Mu2TrackMass",&var_Mu2TrackMass);
   TMVA_Tree->Branch("var_Mu3TrackMass",&var_Mu3TrackMass);
   TMVA_Tree->Branch("var_dcaTrackPV",&var_dcaTrackPV);

   // Track
   TMVA_Tree->Branch("var_NMuTrkOS1", &var_NMuTrkOS1);
   TMVA_Tree->Branch("var_NMuTrkOS2", &var_NMuTrkOS2);
   TMVA_Tree->Branch("var_NMuTrkOS3", &var_NMuTrkOS3);

   TMVA_Tree->Branch("var_MuTrkInvMassOS1Pair1", &var_MuTrkInvMassOS1Pair1);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS1Pair1", &var_MuKTrkInvMassOS1Pair1);
   TMVA_Tree->Branch("var_MuTrkOS1Pair1dR", &var_MuTrkOS1Pair1dR);
   TMVA_Tree->Branch("var_MuTrkInvMassOS2Pair1", &var_MuTrkInvMassOS2Pair1);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS2Pair1", &var_MuKTrkInvMassOS2Pair1);
   TMVA_Tree->Branch("var_MuTrkOS2Pair1dR", &var_MuTrkOS2Pair1dR);
   TMVA_Tree->Branch("var_MuTrkInvMassOS1Pair2", &var_MuTrkInvMassOS1Pair2);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS1Pair2", &var_MuKTrkInvMassOS1Pair2);
   TMVA_Tree->Branch("var_MuTrkOS1Pair2dR", &var_MuTrkOS1Pair2dR);
   TMVA_Tree->Branch("var_MuTrkInvMassOS2Pair2", &var_MuTrkInvMassOS2Pair2);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS2Pair2 ", &var_MuKTrkInvMassOS2Pair2 );
   TMVA_Tree->Branch("var_MuTrkOS2Pair2dR", &var_MuTrkOS2Pair2dR);
   TMVA_Tree->Branch("var_MuTrkInvMassOS1Pair3", &var_MuTrkInvMassOS1Pair3);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS1Pair3", &var_MuKTrkInvMassOS1Pair3);
   TMVA_Tree->Branch("var_MuTrkOS1Pair3dR", &var_MuTrkOS1Pair3dR);
   TMVA_Tree->Branch("var_MuTrkInvMassOS2Pair3", &var_MuTrkInvMassOS2Pair3);
   TMVA_Tree->Branch("var_MuKTrkInvMassOS2Pair3", &var_MuKTrkInvMassOS2Pair3);
   TMVA_Tree->Branch("var_MuTrkOS2Pair3dR", &var_MuTrkOS2Pair3dR);

   TMVA_Tree->Branch("mu3_trk_InvMass", &mu3_trk_InvMass);

   for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==SignalCandidate)     cut.at(SignalCandidate)=1;
      if(i==L1Fired)             cut.at(L1Fired)=1;
      if(i==HLTFired)            cut.at(HLTFired)=1;
      if(i==KFChi2)              cut.at(KFChi2)=999;
      if(i==BSSVSig)             cut.at(BSSVSig)=2;
      if(i==PFMuons)             cut.at(PFMuons)=1;
      if(i==Mu1PtCut)            cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)            cut.at(Mu2PtCut)=3.0;
      if(i==Mu3PtCut)            cut.at(Mu3PtCut)=1.2;
      if(i==TauMassCut)          cut.at(TauMassCut)=1; // true for MC and mass side band for data
      if(i==TriggerMatch)        cut.at(TriggerMatch)=1;
      if(i==MuonID)              cut.at(MuonID)=1;
      if(i==PhiVetoOS1)          cut.at(PhiVetoOS1)=0; // defined below
      if(i==PhiVetoOS2)          cut.at(PhiVetoOS2)=0; // defined below
      if(i==OmegaVetoOS1)        cut.at(OmegaVetoOS1)=0; // defined below
      if(i==OmegaVetoOS2)        cut.at(OmegaVetoOS2)=0; // defined below
      if(i==TrackerLayers)       cut.at(TrackerLayers)=1;
   }

   TString hlabel;
   TString htitle;
   for(int i=0; i<NCuts; i++){
      title.push_back("");
      distindx.push_back(false);
      dist.push_back(std::vector<float>());
      TString c="_Cut_";c+=i;
      if(i==SignalCandidate){
         title.at(i)="signal candidates";
         hlabel="3mu candidates";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidates_",htitle,19,1.0,20.0,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidates_",htitle,19,1.0,20.0,hlabel,"Entries"));
      }
      else if(i==L1Fired){
         title.at(i)="L1 Fired";
         hlabel="DoubleMu/TripleMu fired";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1Fired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1Fired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==HLTFired){
         title.at(i)="HLT Fired";
         hlabel="DoubleMu3_Trk_Tau3mu";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTFired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTFired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==KFChi2){
         title.at(i)="chi square of the vertex";
         hlabel="chi square of the vertex";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_KFChi2_",htitle,50,1.0,500.0,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_KFChi2_",htitle,50,1.0,500.0,hlabel,"Entries"));
      }
      else if(i==BSSVSig){
         title.at(i)="displacement significance of SV with beam spot";
         hlabel="BS-SV displacement significance";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_BSSVSig_",htitle,50,1.0,500.0,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_BSSVSig_",htitle,50,1.0,500.0,hlabel,"Entries"));
      }
      if(i==PFMuons){
         title.at(i)="3 PF Muons";
         hlabel="Particle Flow Muons";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PFMuons_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PFMuons_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==Mu1PtCut){
         title.at(i)="$p_{T}(\\mu_{1}) >$";
         title.at(i)+=cut.at(Mu1PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Entries"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="$p_{T}(\\mu_{2}) >$";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,0,20,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,0,20,hlabel,"Entries"));
      }
      else if(i==Mu3PtCut){
         title.at(i)="$p_{T}(\\mu_{3}) >$";
         title.at(i)+="(2.0/1.2) GeV";
         htitle=title.at(i);
         htitle.ReplaceAll("$","");
         htitle.ReplaceAll("\\","#");
         hlabel="Muon3 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,30,0,15,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,30,0,15,hlabel,"Entries"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="gl,gl,gl";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==PhiVetoOS1){
         title.at(i)="phi mass veto (OS1)";
         hlabel="Phi mass Veto (OS1), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVetoOS1_",htitle,60,0.8,1.2,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVetoOS1_",htitle,60,0.8,1.2,hlabel,"Entries"));
      }
      else if(i==PhiVetoOS2){
         title.at(i)="phi mass veto (OS2)";
         hlabel="Phi mass Veto (OS2), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVetoOS2_",htitle,60,0.8,1.2,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVetoOS2_",htitle,60,0.8,1.2,hlabel,"Entries"));
      }
      else if(i==OmegaVetoOS1){
         title.at(i)="omega mass veto (OS1)";
         hlabel="Omega mass veto (OS1), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVetoOS1_",htitle,50,0.4,0.9,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVetoOS1_",htitle,50,0.4,0.9,hlabel,"Entries"));
      }
      else if(i==OmegaVetoOS2){
         title.at(i)="omega mass veto (OS2)";
         hlabel="Omega mass veto (OS2), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVetoOS2_",htitle,50,0.4,0.9,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVetoOS2_",htitle,50,0.4,0.9,hlabel,"Entries"));
      }
      else if(i==TriggerMatch){
         title.at(i)="Trigger Matching";
         hlabel="trigger matching";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==TauMassCut){
         title.at(i)="Tau Mass";
         hlabel="three mu mass, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,50,1.60,2.1,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,50,1.60,2.1,hlabel,"Entries"));
      }
      else if(i==TrackerLayers){
         title.at(i)="Tracker Layers with measurements";
         hlabel="tracker layers";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMassCut_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMassCut_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
   } 

   Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Entries"); // Do not remove
   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove

}

void  T3MSelectionTree::Store_ExtraDist(){ 
   // Store additional histograms
}



void  T3MSelectionTree::doEvent(){ 
   value.at(KFChi2)=999.0;
   value.at(BSSVSig)=0;
   value.at(PFMuons)=0;
   value.at(HLTFired)=0;
   value.at(L1Fired)=0;
   value.at(SignalCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(Mu3PtCut)=0;
   value.at(MuonID)=0;
   value.at(TriggerMatch)=0;
   value.at(PhiVetoOS1)=99.0;
   value.at(PhiVetoOS2)=99.0;
   value.at(OmegaVetoOS1)=99.0;
   value.at(OmegaVetoOS2)=99.0;
   value.at(TauMassCut)=0;
   value.at(TrackerLayers)=0;

   eventNumber = Ntp->EventNumber();
   run = Ntp->RunNumber();
   lumi = Ntp->LuminosityBlock();

   random_num = rndm.Uniform();

   unsigned int t;
   int id(Ntp->GetMCID());

   DataMCType = id;

   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
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
   if (DoubleMuFired || TripleMuFired) L1Ok = true;

   if (HLTOk) value.at(HLTFired) = true;
   else value.at(HLTFired) = false;

   if (L1Ok) value.at(L1Fired) = true;
   else value.at(L1Fired) = false;

   if (HLTFired && !L1Fired && !randomFailed) cout<<"wrong hlt"<<endl;

   pass.at(HLTFired) = (value.at(HLTFired) == cut.at(HLTFired));
   pass.at(L1Fired) = (value.at(L1Fired) == cut.at(L1Fired));

   l1seed = DoubleMu0Fired+2*DoubleMu4Fired+4*TripleMuFired;

   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_tau = 99.0;

   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_tau = 0.;
   double sumPtTracks_iso05 = 0.;

   double sumTracksIsoTracks_mu1 = 0;
   double sumTracksIsoTracks_mu3 = 0;
   double sumTracksIsoTracks_mu2 = 0;

   double sumTracksIsoTracks_tau = 0.;
   double sumTracksIsoTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_tau = 0;

   // Initialize invariant mass variables

   var_MuTrkInvMassOS1Pair1 = -1;
   var_MuKTrkInvMassOS1Pair1 = -1;
   var_MuTrkOS1Pair1dR = -1;
   var_MuTrkInvMassOS2Pair1 = -1;
   var_MuKTrkInvMassOS2Pair1 = -1;
   var_MuTrkOS2Pair1dR = -1;
   var_MuTrkInvMassOS1Pair2 = -1;
   var_MuKTrkInvMassOS1Pair2 = -1;
   var_MuTrkOS1Pair2dR = -1;
   var_MuTrkInvMassOS2Pair2 = -1;
   var_MuKTrkInvMassOS2Pair2  = -1;
   var_MuTrkOS2Pair2dR = -1;
   var_MuTrkInvMassOS1Pair3 = -1;
   var_MuKTrkInvMassOS1Pair3 = -1;
   var_MuTrkOS1Pair3dR = -1;
   var_MuTrkInvMassOS2Pair3 = -1;
   var_MuKTrkInvMassOS2Pair3 = -1;
   var_MuTrkOS2Pair3dR = -1;

   // Selection of the best candidate

   vector<unsigned int> selectedIndices_threeGlobal;
   vector<unsigned int> selectedIndices_twoGlobalTracker;
   vector<unsigned int> candidateRank_threeGlobal;
   vector<unsigned int> candidateRank_twoGlobalTracker;

   unsigned int final_idx = 0;
   double minChiSq = 999.0;
   value.at(SignalCandidate) = Ntp->NThreeMuons();

   if (Ntp->NThreeMuons()>0){
      for (size_t j=0; j<Ntp->NThreeMuons(); ++j){
         Muon_index_1 = Ntp->ThreeMuonIndices(j).at(0);
         Muon_index_2 = Ntp->ThreeMuonIndices(j).at(1);
         Muon_index_3 = Ntp->ThreeMuonIndices(j).at(2);

         value.at(KFChi2) = Ntp->Vertex_signal_KF_Chi2(j, false);
         value.at(BSSVSig) = Ntp->Vertex_signal_KF_BS_significance(j, false);

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

         value.at(PFMuons) = (Ntp->Muon_isPFMuon(Muon_index_1) && Ntp->Muon_isPFMuon(Muon_index_2) && Ntp->Muon_isPFMuon(Muon_index_3));
         //
         value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) &&
               Ntp->Muon_isGlobalMuon(Muon_index_2) &&
               ( Ntp->Muon_isGlobalMuon(Muon_index_3) || Ntp->Muon_isTrackerMuon(Muon_index_3) ));
         //------------------------------------------------------------------------------------------------------

         if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal==true;
         else threeGlobal = false;
         value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
         value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
         value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();

         TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)+Ntp->Muon_P4(Muon_index_2)+Ntp->Muon_P4(Muon_index_3);

         value.at(PhiVetoOS1) = M_osss1;
         value.at(PhiVetoOS2) = M_osss2;
         value.at(OmegaVetoOS1) = M_osss1;
         value.at(OmegaVetoOS2) = M_osss2;
         vector<TLorentzVector> trigobjTriplet;
         for (int i=0; i<Ntp->NTriggerObjects(); i++){
            TString name = Ntp->TriggerObject_name(i);
            if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
            TLorentzVector tmp;
            tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
            trigobjTriplet.push_back(tmp);
         }

         std::vector<TLorentzVector> muonTriplet;
         muonTriplet.push_back(Ntp->Muon_P4(Muon_index_1));
         muonTriplet.push_back(Ntp->Muon_P4(Muon_index_2));
         muonTriplet.push_back(Ntp->Muon_P4(Muon_index_3));

         bool triggerCheck = false;
         if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet);
         value.at(TriggerMatch) = triggerCheck;
         value.at(TauMassCut) = TauLV.M();

         value.at(TrackerLayers) = ( (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7) && (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7) && (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7));
         pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
         pass.at(KFChi2) = (value.at(KFChi2)<cut.at(KFChi2));
         pass.at(BSSVSig) = (value.at(BSSVSig)>=cut.at(BSSVSig));
         pass.at(PFMuons) = (value.at(PFMuons) == cut.at(PFMuons));
         pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
         pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
         if (threeGlobal) pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > 2.0);
         else pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > 1.2);
         pass.at(MuonID) = (value.at(MuonID)  == cut.at(MuonID));
         pass.at(TauMassCut) = ( value.at(TauMassCut) > tauMinSideBand_ && value.at(TauMassCut) < tauMaxSideBand_);
         pass.at(TriggerMatch) = true;
         pass.at(PhiVetoOS1) = true;
         pass.at(PhiVetoOS2) = true;
         pass.at(OmegaVetoOS1) = true;
         pass.at(OmegaVetoOS2) = true;
         pass.at(TrackerLayers) = true;
         unsigned int score = 0;         
         bool done=false;
         for (unsigned int k=0; k<NCuts; ++k) {
            if (done) continue;
            if (pass.at(k)) score++;
            else done=true;
         }

         if (score==NCuts){
            if(threeGlobal) selectedIndices_threeGlobal.push_back(j);
            else selectedIndices_twoGlobalTracker.push_back(j);
         }

         if (threeGlobal) candidateRank_threeGlobal.push_back(score);
         else candidateRank_twoGlobalTracker.push_back(score);

      }

      if (selectedIndices_threeGlobal.size()==0) {
         if (selectedIndices_twoGlobalTracker.size()==0) {
            if (candidateRank_threeGlobal.size()!=0) final_idx=std::distance(candidateRank_threeGlobal.begin(), std::max_element(candidateRank_threeGlobal.begin(), candidateRank_threeGlobal.end()));
            else final_idx = std::distance(candidateRank_twoGlobalTracker.begin(), std::max_element(candidateRank_twoGlobalTracker.begin(), candidateRank_twoGlobalTracker.end()));
         }
         else {
            for (size_t j=0; j<selectedIndices_twoGlobalTracker.size(); j++){
               int _ = selectedIndices_twoGlobalTracker.at(j);
               double tmpchi = Ntp->ThreeMuons_SV_Chi2(_);
               if (tmpchi<minChiSq) {
                  minChiSq = tmpchi;
                  final_idx = _;
               }
            }
         }
      }

      else{
         for (size_t j=0; j<selectedIndices_threeGlobal.size(); ++j){
            int _ = selectedIndices_threeGlobal.at(j);
            double tmpchi = Ntp->ThreeMuons_SV_Chi2(_);
            if (tmpchi<minChiSq) {
               minChiSq = tmpchi;
               final_idx = _;
            }
         }
      }

      Muon_index_1 = Ntp->ThreeMuonIndices(final_idx).at(0);
      Muon_index_2 = Ntp->ThreeMuonIndices(final_idx).at(1);
      Muon_index_3 = Ntp->ThreeMuonIndices(final_idx).at(2);

      value.at(KFChi2) = Ntp->Vertex_signal_KF_Chi2(final_idx, false);
      value.at(BSSVSig) = Ntp->Vertex_signal_KF_BS_significance(final_idx, false);

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
      value.at(PFMuons) = (Ntp->Muon_isPFMuon(Muon_index_1) && Ntp->Muon_isPFMuon(Muon_index_2) && Ntp->Muon_isPFMuon(Muon_index_3));
      //
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Muon_index_1) &&
            Ntp->Muon_isGlobalMuon(Muon_index_2) &&
            (Ntp->Muon_isGlobalMuon(Muon_index_3) || Ntp->Muon_isTrackerMuon(Muon_index_3)));
      //------------------------------------------------------------------------------------------------------

      if (Ntp->Muon_isGlobalMuon(Muon_index_3)) threeGlobal = true;
      else threeGlobal = false;
      value.at(Mu1PtCut) = Ntp->Muon_P4(Muon_index_1).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(Muon_index_2).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(Muon_index_3).Pt();
      vector<TLorentzVector> trigobjTriplet;

      for (int i=0; i<Ntp->NTriggerObjects(); i++){
         TString name = Ntp->TriggerObject_name(i);
         if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
         TLorentzVector tmp;
         tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
         trigobjTriplet.push_back(tmp);
      }

      std::vector<TLorentzVector> muonTriplet;
      muonTriplet.push_back(Ntp->Muon_P4(Muon_index_1));
      muonTriplet.push_back(Ntp->Muon_P4(Muon_index_2));
      muonTriplet.push_back(Ntp->Muon_P4(Muon_index_3));

      bool triggerCheck = false;
      if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet);
      value.at(TriggerMatch) = triggerCheck;

      value.at(TauMassCut) = TauLV.M();
      value.at(TrackerLayers) = ( (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7) && (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7) && (Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1)>=7));

      pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
      pass.at(KFChi2) = (value.at(KFChi2)<cut.at(KFChi2));
      pass.at(BSSVSig) = (value.at(BSSVSig)>=cut.at(BSSVSig));
      pass.at(PFMuons) = ( value.at(PFMuons) == cut.at(PFMuons) );
      pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
      pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
      if (threeGlobal) pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > 2.0);
      else pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > 1.2);
      pass.at(MuonID) = (value.at(MuonID)  == cut.at(MuonID));
      pass.at(TriggerMatch) = (value.at(TriggerMatch) == cut.at(TriggerMatch));
      pass.at(PhiVetoOS1) = (fabs(value.at(PhiVetoOS1)-PDG_Var::Phi_mass()) > phiVetoSigma);
      pass.at(PhiVetoOS2) = (fabs(value.at(PhiVetoOS2)-PDG_Var::Phi_mass()) > phiVetoSigma);
      pass.at(OmegaVetoOS1) = (fabs(value.at(OmegaVetoOS1)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(OmegaVetoOS2) = (fabs(value.at(OmegaVetoOS2)-PDG_Var::Omega_mass())> omegaVetoSigma);
      pass.at(TauMassCut) = ( value.at(TauMassCut) > tauMinSideBand_ && value.at(TauMassCut) <  tauMaxSideBand_);
      pass.at(TrackerLayers) = ( value.at(TrackerLayers) == cut.at(TrackerLayers) );
   }


   // Weights
   double wobs=1;
   double w;

   if(!Ntp->isData()){
      w = puWeights->GetBinContent(Ntp->TruthNumberOfInteraction()); // Weight MC according to truth number of vertices
   }
   //  No weights to data
   else{w=1;}

   EventWeight = w;

   bool status = AnalysisCuts(t, w, wobs);    
   if(status){
      
      /*
         int NJets(){return Ntp->Jet_p4->size();}
         TLorentzVector Jet_P4(unsigned int i){return TLorentzVector(Ntp->Jet_p4->at(i)->at(1), Ntp->Jet_p4->at(i)->at(2), Ntp->Jet_p4->at(i)->at(3), Ntp->Jet_p4->at(i)->at(0));}
         double JetBTagCVSB(unsigned int i){return Ntp->Jet_BTagCVSB->at(i);}
         double JetBTagMVA(unsigned int i){return Ntp->Jet_BTagMVA->at(i);}
         double JetBTagCSV(unsigned int i){return Ntp->Jet_BTagCSV->at(i);}
         */

      double mumutrk_chi1(999.), mumutrk_chi2(999.), mumutrk_chi3(999.);
      double mutrk_chi12(999.), mutrk_chi13(999.);
      double mutrk_chi21(999.), mutrk_chi23(999.);
      double mutrk_chi31(999.), mutrk_chi32(999.);
      int NMuTrkOS1(0), NMuTrkOS2(0), NMuTrkOS3(0);

      int mumutrk_trk1_idx = -1, mumutrk_trk2_idx = -1, mumutrk_trk3_idx = -1;
      int mutrk_trk32_idx = -1, mutrk_trk31_idx = -1;
      int mutrk_trk13_idx = -1, mutrk_trk12_idx = -1;
      int mutrk_trk23_idx = -1, mutrk_trk21_idx = -1;
      int addTrackCounter1=0, addTrackCounter2=0, addTrackCounter3=0;
      std::vector<unsigned int> indices;

      indices.push_back(Ntp->ThreeMuonIndices(final_idx).at(0));
      indices.push_back(Ntp->ThreeMuonIndices(final_idx).at(1));
      indices.push_back(Ntp->ThreeMuonIndices(final_idx).at(2));

      Muon_index_1 = Ntp->SortedPtMuons(indices).at(0);
      Muon_index_2 = Ntp->SortedPtMuons(indices).at(1);
      Muon_index_3 = Ntp->SortedPtMuons(indices).at(2);

      int closest_track_index = -1;
      float min_track_dr = 999.;

      // Find if a macthing isolation track is found
      for (int itrk=0; itrk<Ntp->NIsolationTrack(final_idx); itrk++){
         if (Ntp->IsolationTrack_charge(final_idx, itrk, false)==Ntp->Muon_charge(Muon_index_3)) continue;
         TLorentzVector trk_p4 = Ntp->IsolationTrack_p4(final_idx, itrk, false);
         if (trk_p4.Pt()<1.0) continue;
         float tmp_dr = trk_p4.DeltaR(Ntp->Muon_P4(Muon_index_3));
         if (tmp_dr<min_track_dr){
            min_track_dr = tmp_dr;
            closest_track_index = itrk;
         }
      }

      mu3_trk_InvMass = -1;
      if (closest_track_index!=-1){
         mu3_trk_InvMass = (Ntp->IsolationTrack_p4(final_idx, closest_track_index, false)+Ntp->Muon_P4(Muon_index_3)).M();
      }

      for (unsigned int iTwoMu=0; iTwoMu<Ntp->NTwoMuonsTrack(); iTwoMu++){
         auto mu1_idx = Ntp->TwoMuonsTrackMuonIndices(iTwoMu).at(0);
         auto mu2_idx = Ntp->TwoMuonsTrackMuonIndices(iTwoMu).at(1);
         auto trk_idx = Ntp->TwoMuonsTrackTrackIndex(iTwoMu).at(0);
         if (Ntp->Track_P4(trk_idx).Pt()<1) continue;
         if ( (Muon_index_1==mu1_idx && Muon_index_2==mu2_idx) || (Muon_index_2==mu1_idx && Muon_index_1==mu2_idx) ) {
            int comb1 = ((Muon_index_1==mu1_idx && Muon_index_2==mu2_idx) ? 1:2);
            int comb2 = ((Muon_index_1==mu1_idx && Muon_index_2==mu2_idx) ? 2:1);
            float dr = Ntp->Track_P4(trk_idx).DeltaR(Ntp->Muon_P4(Muon_index_3));
            if (dr>0.01) {
               // find the best fitting track (mutrk)
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb1,true) &&
                     mutrk_chi32>Ntp->Vertex_pair_quality(iTwoMu,comb1,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_2)) {
                  mutrk_chi32 = Ntp->Vertex_pair_quality(iTwoMu,comb1,true);
                  mutrk_trk32_idx = trk_idx;
               }
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb2,true) &&
                     mutrk_chi31>Ntp->Vertex_pair_quality(iTwoMu,comb2,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_1)) {
                  mutrk_chi31 = Ntp->Vertex_pair_quality(iTwoMu,comb2,true);
                  mutrk_trk31_idx = trk_idx;
               }
               // find the best fitting track (mumutrk)
               double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(iTwoMu);
               addTrackCounter3++;
               if(tmpchi<mumutrk_chi3) {
                  mumutrk_chi3 = tmpchi;
                  mumutrk_trk3_idx = trk_idx;
               }
            }
         }
         else if ( (Muon_index_3==mu1_idx && Muon_index_1==mu2_idx) || (Muon_index_1==mu1_idx && Muon_index_3==mu2_idx) ) {
            float dr = Ntp->Track_P4(trk_idx).DeltaR(Ntp->Muon_P4(Muon_index_2));
            if (dr>0.01) {
               int comb1 = ((Muon_index_3==mu1_idx && Muon_index_1==mu2_idx) ? 0:1);
               int comb2 = ((Muon_index_3==mu1_idx && Muon_index_1==mu2_idx) ? 1:0);
               // find the best fitting track (mutrk)
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb1,true) &&
                     mutrk_chi21>Ntp->Vertex_pair_quality(iTwoMu,comb1,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_1)) {
                  mutrk_chi21 = Ntp->Vertex_pair_quality(iTwoMu,comb1,true);
                  mutrk_trk21_idx = trk_idx;
               }
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb2,true) &&
                     mutrk_chi23>Ntp->Vertex_pair_quality(iTwoMu,comb2,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_3)) {
                  mutrk_chi23 = Ntp->Vertex_pair_quality(iTwoMu,comb2,true);
                  mutrk_trk23_idx = trk_idx;
               }

               // find the best fitting track
               double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(iTwoMu);
               addTrackCounter2++;
               if(tmpchi<mumutrk_chi2) {
                  mumutrk_chi2 = tmpchi;
                  mumutrk_trk2_idx = trk_idx;
               }
            }
         }
         else if ( (Muon_index_2==mu1_idx && Muon_index_3==mu2_idx) || (Muon_index_3==mu1_idx && Muon_index_2==mu2_idx) ) {
            float dr = Ntp->Track_P4(trk_idx).DeltaR(Ntp->Muon_P4(Muon_index_1));
            if (dr>0.01) {
               int comb1 = ((Muon_index_2==mu1_idx && Muon_index_3==mu2_idx) ? 0:2);
               int comb2 = ((Muon_index_2==mu1_idx && Muon_index_3==mu2_idx) ? 2:0);
               // find the best fitting track (mutrk)
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb1,true) &&
                     mutrk_chi12>Ntp->Vertex_pair_quality(iTwoMu,comb1,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_2)) {
                  mutrk_chi12 = Ntp->Vertex_pair_quality(iTwoMu,comb1,true);
                  mutrk_trk12_idx = trk_idx;
               }
               if ( Ntp->Vertex_pairfit_status(iTwoMu,comb2,true) &&
                     mutrk_chi13>Ntp->Vertex_pair_quality(iTwoMu,comb2,true) &&
                     Ntp->Track_charge(trk_idx)!=Ntp->Muon_charge(Muon_index_3)) {
                  mutrk_chi13 = Ntp->Vertex_pair_quality(iTwoMu,comb2,true);
                  mutrk_trk13_idx = trk_idx;
               }
               // find the best fitting track
               double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(iTwoMu);
               addTrackCounter1++;
               if(tmpchi<mumutrk_chi1) {
                  mumutrk_chi1 = tmpchi;
                  mumutrk_trk1_idx = trk_idx;
               }
            }
         }

         // Get track indices not matched to muons
         TLorentzVector mu1_p4 = Ntp->Muon_P4(mu1_idx);
         TLorentzVector mu2_p4 = Ntp->Muon_P4(mu2_idx);
         TLorentzVector trk_p4 = Ntp->Track_P4(trk_idx);
      }

      TLorentzVector vecMu1 = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector vecMu2 = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector vecMu3 = Ntp->Muon_P4(Muon_index_3);
      auto vecTau = vecMu1+vecMu2+vecMu3;
      double tauMass = vecTau.M();

      if (mutrk_chi31>=15 && mutrk_chi32>=15) NMuTrkOS3 = 0;
      else if (mutrk_chi31<15 && mutrk_chi32>=15) NMuTrkOS3 = 1;
      else if (mutrk_chi31>=15 && mutrk_chi32<15) NMuTrkOS3 = 2;
      else NMuTrkOS3 = 3;

      if (mutrk_chi21>=15 && mutrk_chi23>=15) NMuTrkOS2 = 0;
      else if (mutrk_chi21<15 && mutrk_chi23>=15) NMuTrkOS2 = 1;
      else if (mutrk_chi21>=15 && mutrk_chi23<15) NMuTrkOS2 = 2;
      else NMuTrkOS2 = 3;

      if (mutrk_chi12>=15 && mutrk_chi13>=15) NMuTrkOS1 = 0;
      else if (mutrk_chi12<15 && mutrk_chi13>=15) NMuTrkOS1 = 1;
      else if (mutrk_chi12>=15 && mutrk_chi13<15) NMuTrkOS1 = 2;
      else NMuTrkOS1 = 3;
      // First pair (muon1+muon2+track)
      if (mutrk_chi31<=mutrk_chi32){
         if (mutrk_chi31<15) {
            var_MuTrkInvMassOS1Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->Track_P4(mutrk_trk31_idx)).M();
            var_MuKTrkInvMassOS1Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->KaonTrack_P4(mutrk_trk31_idx)).M();
            var_MuTrkOS1Pair1dR = Ntp->Muon_P4(Muon_index_3).DeltaR(Ntp->Track_P4(mutrk_trk31_idx));
         }
         if (mutrk_chi32<15) {
            var_MuTrkInvMassOS2Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->Track_P4(mutrk_trk32_idx)).M();
            var_MuKTrkInvMassOS2Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->KaonTrack_P4(mutrk_trk32_idx)).M();
            var_MuTrkOS2Pair1dR = Ntp->Muon_P4(Muon_index_3).DeltaR(Ntp->Track_P4(mutrk_trk32_idx));
         }
      } else {
         if (mutrk_chi32<15) {
            var_MuTrkInvMassOS1Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->Track_P4(mutrk_trk32_idx)).M();
            var_MuKTrkInvMassOS1Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->KaonTrack_P4(mutrk_trk32_idx)).M();
            var_MuTrkOS1Pair1dR = Ntp->Muon_P4(Muon_index_3).DeltaR(Ntp->Track_P4(mutrk_trk32_idx));
         }
         if (mutrk_chi31<15) {
            var_MuTrkInvMassOS2Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->Track_P4(mutrk_trk31_idx)).M();
            var_MuKTrkInvMassOS2Pair1 = (Ntp->Muon_P4(Muon_index_3)+Ntp->KaonTrack_P4(mutrk_trk31_idx)).M();
            var_MuTrkOS2Pair1dR = Ntp->Muon_P4(Muon_index_3).DeltaR(Ntp->Track_P4(mutrk_trk31_idx));
         }
      }

      // Second pair (muon2+muon3+track)
      if (mutrk_chi12<=mutrk_chi13){
         if (mutrk_chi12<15) {
            var_MuTrkInvMassOS1Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->Track_P4(mutrk_trk12_idx)).M();
            var_MuKTrkInvMassOS1Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->KaonTrack_P4(mutrk_trk12_idx)).M();
            var_MuTrkOS1Pair2dR = Ntp->Muon_P4(Muon_index_1).DeltaR(Ntp->Track_P4(mutrk_trk12_idx));
         }
         if (mutrk_chi13<15) {
            var_MuTrkInvMassOS2Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->Track_P4(mutrk_trk13_idx)).M();
            var_MuKTrkInvMassOS2Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->KaonTrack_P4(mutrk_trk13_idx)).M();
            var_MuTrkOS2Pair2dR = Ntp->Muon_P4(Muon_index_1).DeltaR(Ntp->Track_P4(mutrk_trk13_idx));
         }
      } else {
         if (mutrk_chi13<15) {
            var_MuTrkInvMassOS1Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->Track_P4(mutrk_trk13_idx)).M();
            var_MuKTrkInvMassOS1Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->KaonTrack_P4(mutrk_trk13_idx)).M();
            var_MuTrkOS1Pair2dR = Ntp->Muon_P4(Muon_index_1).DeltaR(Ntp->Track_P4(mutrk_trk13_idx));
         }
         if (mutrk_chi12<15) {
            var_MuTrkInvMassOS2Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->Track_P4(mutrk_trk12_idx)).M();
            var_MuKTrkInvMassOS2Pair2 = (Ntp->Muon_P4(Muon_index_1)+Ntp->KaonTrack_P4(mutrk_trk12_idx)).M();
            var_MuTrkOS2Pair2dR = Ntp->Muon_P4(Muon_index_1).DeltaR(Ntp->Track_P4(mutrk_trk12_idx));
         }
      }

      // Third pair (muon1+muon3+track)
      if (mutrk_chi21<=mutrk_chi23){
         if (mutrk_chi21<15) {
            var_MuTrkInvMassOS1Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->Track_P4(mutrk_trk21_idx)).M();
            var_MuKTrkInvMassOS1Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->KaonTrack_P4(mutrk_trk21_idx)).M();
            var_MuTrkOS1Pair3dR = Ntp->Muon_P4(Muon_index_2).DeltaR(Ntp->Track_P4(mutrk_trk21_idx));
         }
         if (mutrk_chi23<15) {
            var_MuTrkInvMassOS2Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->Track_P4(mutrk_trk23_idx)).M();
            var_MuKTrkInvMassOS2Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->KaonTrack_P4(mutrk_trk23_idx)).M();
            var_MuTrkOS2Pair3dR = Ntp->Muon_P4(Muon_index_2).DeltaR(Ntp->Track_P4(mutrk_trk23_idx));
         }
      } else {
         if (mutrk_chi23<15) {
            var_MuTrkInvMassOS1Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->Track_P4(mutrk_trk23_idx)).M();
            var_MuKTrkInvMassOS1Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->KaonTrack_P4(mutrk_trk23_idx)).M();
            var_MuTrkOS1Pair3dR = Ntp->Muon_P4(Muon_index_2).DeltaR(Ntp->Track_P4(mutrk_trk23_idx));
         }
         if (mutrk_chi21<15) {
            var_MuTrkInvMassOS2Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->Track_P4(mutrk_trk21_idx)).M();
            var_MuKTrkInvMassOS2Pair3 = (Ntp->Muon_P4(Muon_index_2)+Ntp->KaonTrack_P4(mutrk_trk21_idx)).M();
            var_MuTrkOS2Pair3dR = Ntp->Muon_P4(Muon_index_2).DeltaR(Ntp->Track_P4(mutrk_trk21_idx));
         }
      }

      Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
      TLorentzVector TauLV = Muon1LV + Muon2LV + Muon3LV;
      TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,false)+
         Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,false)+
         Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,false);

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
            BTagCSVSH = Ntp->JetBTagCSV(ijet);
            BTagMVASH = Ntp->JetBTagMVA(ijet);
            BTagCVSBSH = Ntp->JetBTagCVSB(ijet);
         }


         if(fabs(Ntp->Jet_P4(ijet).Phi()  - TauLV.Phi()) > TMath::Pi() ){
            OH_jet_count++;  
            BTagCSVOH = Ntp->JetBTagCSV(ijet);
            BTagMVAOH = Ntp->JetBTagMVA(ijet);
            BTagCVSBOH = Ntp->JetBTagCVSB(ijet);
         }


         if(jet_index_closest_to_tau!=-1){
            if(fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() )) < OHDeltaPhi  )  {
               OHDeltaPhi = fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() ));
               jet_index_closest_to_tau_opposite_in_phi = ijet;
            }
         }
      }


      if(jet_index_closest_to_tau!=-1){
         BTagCSVSHMatchedToTau = Ntp->JetBTagCSV(jet_index_closest_to_tau);
         BTagMVASHMatchedToTau = Ntp->JetBTagMVA(jet_index_closest_to_tau);
         BTagCVSBSHMatchedToTau = Ntp->JetBTagCVSB(jet_index_closest_to_tau);
      }
      if(jet_index_closest_to_tau_opposite_in_phi!=-1){
         BTagCSVOHMatchedToTau = Ntp->JetBTagCSV(jet_index_closest_to_tau_opposite_in_phi);
         BTagMVAOHMatchedToTau = Ntp->JetBTagMVA(jet_index_closest_to_tau_opposite_in_phi);
         BTagCVSBOHMatchedToTau = Ntp->JetBTagCVSB(jet_index_closest_to_tau_opposite_in_phi);
      }

      //------------- Muon ID

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_1).size(); iMuSelector++ ){
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)==1)  Muon1StandardSelectorPass = iMuSelector;
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)!=1)  Muon1StandardSelectorFail = iMuSelector;
      }

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_2).size(); iMuSelector++ ){
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)==1)  Muon2StandardSelectorPass = iMuSelector;
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)!=1)  Muon2StandardSelectorFail = iMuSelector;
      }

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_3).size(); iMuSelector++ ){
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)==1)  Muon3StandardSelectorPass = iMuSelector;
         if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)!=1)  Muon3StandardSelectorFail = iMuSelector;
      }


      //------------- Muon ID

      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      std::vector<unsigned int> EtaSortedIndices;

      EtaSortedIndices.push_back(Muon_Eta_index_1);
      EtaSortedIndices.push_back(Muon_Eta_index_2);
      EtaSortedIndices.push_back(Muon_Eta_index_3);
      EMR_tau_eta = Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta();
      EventMassResolution_PtEtaPhi = Ntp->TauMassResolution(EtaSortedIndices,1,false);

      float tauMassRes = Ntp->TauMassResolution(EtaSortedIndices,1,false);

      var_Muon1_Pt = Ntp->Muon_P4(Muon_index_1).Pt();
      var_Muon2_Pt = Ntp->Muon_P4(Muon_index_2).Pt();
      var_Muon3_Pt = Ntp->Muon_P4(Muon_index_3).Pt();

      var_Tau_Eta = TauLV.Eta();
      var_Tau_Pt = TauLV.Pt();
      var_Tau_Phi = TauLV.Phi();

      //------------------ calculate var_svpvTauAngle ---------------------
      TVector3 vec_sv = Ntp->Vertex_Signal_KF_pos(final_idx,false);
      TVector3 vec_pv(0,0,0);
      if (Ntp->Vertex_RefitPVisValid(final_idx,false)==1){
         vec_pv = Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx,false);
      }
      TVector3 vec_tau = TauLV.Vect();
      TVector3 d_pvsv = vec_sv - vec_pv;

      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,false),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));
      TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_1),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));
      TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_2),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));
      TVector3 Mu3ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_3),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));


      //---------------------------------------------------------------
      //------------------ calculate var_flightLenSig ---------------------
      TMatrixTSym<double> fls_PVcov = Ntp->Vertex_PrimaryVertex_Covariance(final_idx,false);
      TMatrixTSym<double> fls_SVcov = Ntp->Vertex_Signal_KF_Covariance(final_idx,false);



      //      if (Ntp->Vertex_RefitPVisValid(final_idx)==1)

      var_Muon1ImpactAngle = SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag());
      var_Muon2ImpactAngle = SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag());
      var_Muon3ImpactAngle = SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag());

      // ----- Fill the histograms -----

      Isolation_NTracks = Ntp->Isolation_NTracks(final_idx,false);
      Isolation_RelPt = Ntp->Isolation_RelPt(final_idx,false);
      Isolation_maxdxy = Ntp->Isolation_maxdy(final_idx,false);

      deltaMuZ12 = fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z());
      deltaMuZ13 = fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z());
      deltaMuZ23 = fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z());

      VertexMu1D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,0,false);
      VertexMu2D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,1,false);
      VertexMu3D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,2,false);

      VertexMu1D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0,false);
      VertexMu2D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1,false);
      VertexMu3D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2,false);

      VertexMu1D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,0,false);
      VertexMu2D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,1,false);
      VertexMu3D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,2,false);


      //  ----------------------------- secondary vertices ----------------
      int NumberOfPrimaryVertices(0);
      for(unsigned int iVertex=0; iVertex < Ntp->NSecondaryVertices(); iVertex++){
         SV_Mass = Ntp->SecondaryVertexMass(iVertex);
         TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,false),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));
         TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(iVertex),Ntp->Vertex_MatchedPrimaryVertex(final_idx,false));

         SVDeltaR = SVfakePV.DeltaR(SVsignalPV);
         SVDistance = (Ntp->Vertex_Signal_KF_pos(final_idx,false) - Ntp->SecondaryVertexPosition(iVertex)).Mag();

         if(SVfakePV.DeltaR(SVsignalPV) < 0.3 && (Ntp->Vertex_Signal_KF_pos(final_idx,false) - Ntp->SecondaryVertexPosition(iVertex)).Mag() > 0.015){ // sv in the tau cone but  displaced
            NumberOfPrimaryVertices++;

         }
      }
      NSV = NumberOfPrimaryVertices;


      // ----------------------------------------------------------------------
      // --------------------------- Isolation --------------------------------
      // ----------------------------------------------------------------------

      int NcloseTracksCount(0);
      int TrackIndex(0);
      int TrackIndex_closestToPV(0);
      double TrackPTtreschold(0.8);
      double dca_temp(999.);
      double dcaPV_temp(999.);

      for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

         if(Ntp->IsolationTrack_p4(final_idx,i,false).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i,false),2)   +   
                  pow(Ntp->IsolationTrack_dxySV(final_idx,i,false),2)) < 0.03)
         {
            NcloseTracksCount++;
         }
         if( sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i,false),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i,false),2) ) <  dca_temp){
            dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i,false),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i,false),2));
            TrackIndex = i;
         }

         if(Ntp->IsolationTrack_p4(final_idx,i,false).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i,false)) < 0.05 && 
               sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i,false),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i,false),2)) < 0.05){


            if( sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i,false),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i,false),2) ) <  dcaPV_temp){
               dcaPV_temp = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i,false),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i,false),2));
               TrackIndex_closestToPV = i;
            }


            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 0.2){
               sumTracksIso02 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 0.4){
               sumTracksIso04 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 0.6){
               sumTracksIso06 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 0.8){
               sumTracksIso08 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 1.0){
               sumTracksIso1 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 1.2){
               sumTracksIso12 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 1.4){
               sumTracksIso14 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }  
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 1.6){
               sumTracksIso16 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 1.8){
               sumTracksIso18 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(TauRefitLV) < 2.0){
               sumTracksIso2 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
         }
      }


      var_NtracksClose = NcloseTracksCount;

      var_Iso02 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso02);
      var_Iso04 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso04);
      var_Iso06 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso06);
      var_Iso08 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso08);
      var_Iso1 = Muon2LV.Pt()/  (Muon2LV.Pt() + sumTracksIso1);
      var_Iso12 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso12);
      var_Iso14 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso14);
      var_Iso16 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso16);
      var_Iso18 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso18);
      var_Iso2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso2);

      if(Ntp->NIsolationTrack(final_idx,false)!=0){ 
         var_MindcaTrackSV = sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex,false),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex,false),2));
      }

      if(Ntp->NIsolationTrack(final_idx)!=0){ 
         var_dcaTrackPV = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV,false),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV,false),2));
      }


      // ---------------------- I_mu1
      sumTracksIso02=0;sumTracksIso04=0;sumTracksIso06=0;sumTracksIso08=0;sumTracksIso1=0;
      sumTracksIso12=0;sumTracksIso14=0;sumTracksIso16=0;sumTracksIso18=0;sumTracksIso2=0;

      for(int i =0; i< Ntp->NIsolationTrack(final_idx,false); i++){

         if(Ntp->IsolationTrack_p4(final_idx,i,false).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i,false)) < 0.05
               && Ntp->IsolationTrack_DocaMu1(final_idx,i,false) < 0.1){
            if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(0))* Ntp->IsolationTrack_charge(final_idx,i,false)==-1 ){
               var_Mu1TrackMass = (Muon1LV +Ntp->IsolationTrack_p4(final_idx,i,false)).M();
            }else var_Mu1TrackMass = -1;

            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 0.2){
               sumTracksIso02 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 0.4){
               sumTracksIso04 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 0.6){
               sumTracksIso06 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 0.8){
               sumTracksIso08 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 1.0){
               sumTracksIso1 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 1.2){
               sumTracksIso12 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 1.4){
               sumTracksIso14 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 1.6){
               sumTracksIso16 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 1.8){
               sumTracksIso18 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon1LV) < 2.0){
               sumTracksIso2 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
         }
      }

      var_Iso02Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso02);
      var_Iso04Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso04);
      var_Iso06Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso06);
      var_Iso08Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso08);
      var_Iso1Mu1 = Muon2LV.Pt()/  (Muon2LV.Pt() + sumTracksIso1);
      var_Iso12Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso12);
      var_Iso14Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso14);
      var_Iso16Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso16);
      var_Iso18Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso18);
      var_Iso2Mu1 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso2);

      // ---------------------- I_mu2

      sumTracksIso02=0; sumTracksIso04=0; sumTracksIso06=0; sumTracksIso08=0; sumTracksIso1=0;
      sumTracksIso12=0; sumTracksIso14=0; sumTracksIso16=0; sumTracksIso18=0; sumTracksIso2=0;

      for(int i =0; i< Ntp->NIsolationTrack(final_idx,false); i++){

         if(Ntp->IsolationTrack_p4(final_idx,i,false).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i,false)) < 0.05
               && Ntp->IsolationTrack_DocaMu2(final_idx,i,false) < 0.1){


            if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(1))* Ntp->IsolationTrack_charge(final_idx,i,false)==-1 ){

               var_Mu2TrackMass = (Muon2LV +Ntp->IsolationTrack_p4(final_idx,i,false)).M();
            }else var_Mu2TrackMass = -1;


            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 0.2){
               sumTracksIso02 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 0.4){
               sumTracksIso04 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 0.6){
               sumTracksIso06 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 0.8){
               sumTracksIso08 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 1.0){
               sumTracksIso1 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 1.2){
               sumTracksIso12 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 1.4){
               sumTracksIso14 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 1.6){
               sumTracksIso16 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 1.8){
               sumTracksIso18 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon2LV) < 2.0){
               sumTracksIso2 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
         }
      }

      var_Iso02Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso02);
      var_Iso04Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso04);
      var_Iso06Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso06);
      var_Iso08Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso08);
      var_Iso1Mu2 = Muon2LV.Pt()/  (Muon2LV.Pt() + sumTracksIso1);
      var_Iso12Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso12);
      var_Iso14Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso14);
      var_Iso16Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso16);
      var_Iso18Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso18);
      var_Iso2Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso2);

      // ---------------------- I_mu3

      sumTracksIso02=0;sumTracksIso04=0;sumTracksIso06=0;sumTracksIso08=0;sumTracksIso1=0;
      sumTracksIso12=0;sumTracksIso14=0;sumTracksIso16=0;sumTracksIso18=0;sumTracksIso2=0;

      for(int i =0; i< Ntp->NIsolationTrack(final_idx,false); i++){

         if(Ntp->IsolationTrack_p4(final_idx,i,false).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i,false)) < 0.05
               && Ntp->IsolationTrack_DocaMu3(final_idx,i,false) < 0.1){


            if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(2))* Ntp->IsolationTrack_charge(final_idx,i,false)==-1 ){
               var_Mu3TrackMass = (Muon3LV +Ntp->IsolationTrack_p4(final_idx,i,false)).M();
            }else var_Mu3TrackMass = -1;


            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 0.2){
               sumTracksIso02 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 0.4){
               sumTracksIso04 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 0.6){
               sumTracksIso06 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 0.8){
               sumTracksIso08 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 1.0){
               sumTracksIso1 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 1.2){
               sumTracksIso12 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 1.4){
               sumTracksIso14 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 1.6){
               sumTracksIso16 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 1.8){
               sumTracksIso18 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
            if(Ntp->IsolationTrack_p4(final_idx,i,false).DeltaR(Muon3LV) < 2.0){
               sumTracksIso2 += Ntp->IsolationTrack_p4(final_idx,i,false).Pt();
            }
         }
      }

      var_Iso02Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso02);
      var_Iso04Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso04);
      var_Iso06Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso06);
      var_Iso08Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso08);
      var_Iso1Mu3 = Muon2LV.Pt()/  (Muon2LV.Pt() + sumTracksIso1);
      var_Iso12Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso12);
      var_Iso14Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso14);
      var_Iso16Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso16);
      var_Iso18Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso18);
      var_Iso2Mu3 = Muon2LV.Pt()/ (Muon2LV.Pt() + sumTracksIso2);

      int nTracksIso0p3 = 0;
      int nTracksIso0p5 = 0;
      int nTracksIso0p8 = 0;
      int nTracksIso1p2 = 0;
      int nTracksIso1p4 = 0;

      TLorentzVector sumTracksIso0p3 = vecTau;
      TLorentzVector sumTracksIso0p5 = vecTau;
      TLorentzVector sumTracksIso0p8 = vecTau;
      TLorentzVector sumTracksIso1p2 = vecTau;
      TLorentzVector sumTracksIso1p4 = vecTau;

      for (int it=0; it<Ntp->NIsolationTrack(final_idx); it++){
         TLorentzVector TrackLV = Ntp->IsolationTrack_p4(final_idx, it);
         if (TrackLV.Pt()<1.0) continue;

         if (TrackLV.DeltaR(vecTau)<0.3){
            sumTracksIso0p3 += TrackLV;
            nTracksIso0p3 += 1;
         }
         if (TrackLV.DeltaR(vecTau)<0.5){
            sumTracksIso0p5 += TrackLV;
            nTracksIso0p5 += 1;
         }
         if (TrackLV.DeltaR(vecTau)<0.8){
            sumTracksIso0p8 += TrackLV;
            nTracksIso0p8 += 1;
         }
         if (TrackLV.DeltaR(vecTau)<1.2){
            sumTracksIso1p2 += TrackLV;
            nTracksIso1p2 += 1;
         }
         if (TrackLV.DeltaR(vecTau)<1.4){
            sumTracksIso1p4 += TrackLV;
            nTracksIso1p4 += 1;
         }
      }

      //----------------------------------------------------------------------
      // count number of stations crossed by muon
      // Segment compatibility variables
      for (int istation=0; istation<8; istation++){
         int station_status = 0;

         if (Ntp->Muon_TrackX(Muon_index_1, istation)<999999.0 && Ntp->Muon_TrackY(Muon_index_1, istation)<999999.0) station_status++;
         if (Ntp->Muon_dX(Muon_index_1, istation)>=999999.0 && station_status) station_status = -1;
         Muon1_station_vars[istation][0] = station_status;
         Muon1_station_vars[istation][1] = Ntp->Muon_TrackX(Muon_index_1, istation);
         Muon1_station_vars[istation][2] = Ntp->Muon_TrackY(Muon_index_1, istation);
         Muon1_station_vars[istation][3] = Ntp->Muon_pullX(Muon_index_1, istation);
         Muon1_station_vars[istation][4] = Ntp->Muon_pullY(Muon_index_1, istation);
         Muon1_station_vars[istation][5] = Ntp->Muon_pullDxDz(Muon_index_1, istation);
         Muon1_station_vars[istation][6] = Ntp->Muon_pullDyDz(Muon_index_1, istation);

         station_status = 0;
         if (Ntp->Muon_TrackX(Muon_index_2, istation)<999999.0 && Ntp->Muon_TrackY(Muon_index_2, istation)<999999.0) station_status++;
         if (Ntp->Muon_dX(Muon_index_2, istation)>=999999.0 && station_status) station_status = -1;
         Muon2_station_vars[istation][0] = station_status;
         Muon2_station_vars[istation][1] = Ntp->Muon_TrackX(Muon_index_2, istation);
         Muon2_station_vars[istation][2] = Ntp->Muon_TrackY(Muon_index_2, istation);
         Muon2_station_vars[istation][3] = Ntp->Muon_pullX(Muon_index_2, istation);
         Muon2_station_vars[istation][4] = Ntp->Muon_pullY(Muon_index_2, istation);
         Muon2_station_vars[istation][5] = Ntp->Muon_pullDxDz(Muon_index_2, istation);
         Muon2_station_vars[istation][6] = Ntp->Muon_pullDyDz(Muon_index_2, istation);

         station_status = 0;
         if (Ntp->Muon_TrackX(Muon_index_3, istation)<999999.0 && Ntp->Muon_TrackY(Muon_index_3, istation)<999999.0) station_status++;
         if (Ntp->Muon_dX(Muon_index_3, istation)>=999999.0 && station_status) station_status = -1;
         Muon3_station_vars[istation][0] = station_status;
         Muon3_station_vars[istation][1] = Ntp->Muon_TrackX(Muon_index_3, istation);
         Muon3_station_vars[istation][2] = Ntp->Muon_TrackY(Muon_index_3, istation);
         Muon3_station_vars[istation][3] = Ntp->Muon_pullX(Muon_index_3, istation);
         Muon3_station_vars[istation][4] = Ntp->Muon_pullY(Muon_index_3, istation);
         Muon3_station_vars[istation][5] = Ntp->Muon_pullDxDz(Muon_index_3, istation);
         Muon3_station_vars[istation][6] = Ntp->Muon_pullDyDz(Muon_index_3, istation);

   }
   // -------------------------- Fill MVA mini tree ------------------------
   // Muon variables
   var_Muon1_Pt = Ntp->Muon_P4(Muon_index_1).Pt();
   var_Muon2_Pt = Ntp->Muon_P4(Muon_index_2).Pt();
   var_Muon3_Pt = Ntp->Muon_P4(Muon_index_3).Pt();

   var_Muon1_Eta = Ntp->Muon_P4(Muon_index_1).Eta();
   var_Muon2_Eta = Ntp->Muon_P4(Muon_index_2).Eta();
   var_Muon3_Eta = Ntp->Muon_P4(Muon_index_3).Eta();

   var_Muon1_Phi = Ntp->Muon_P4(Muon_index_1).Phi();
   var_Muon2_Phi = Ntp->Muon_P4(Muon_index_2).Phi();
   var_Muon3_Phi = Ntp->Muon_P4(Muon_index_3).Phi();

   var_Muon1_vx = Ntp->Muon_Poca(Muon_index_1).X();
   var_Muon1_vy = Ntp->Muon_Poca(Muon_index_1).Y();
   var_Muon1_vz = Ntp->Muon_Poca(Muon_index_1).Z();
   var_Muon1_IsGlobalMuon = Ntp->Muon_isGlobalMuon(Muon_index_1);
   var_Muon1_IsStandAloneMuon = Ntp->Muon_isStandAloneMuon(Muon_index_1);
   var_Muon1_IsTrackerMuon = Ntp->Muon_isTrackerMuon(Muon_index_1);
   var_Muon1_IsCaloMuon = Ntp->Muon_isCaloMuon(Muon_index_1);
   var_Muon1_IsIsolationValid = Ntp->Muon_isIsolationValid(Muon_index_1);
   var_Muon1_IsQualityValid = Ntp->Muon_isQualityValid(Muon_index_1);
   var_Muon1_IsTimeValid = Ntp->Muon_isTimeValid(Muon_index_1);
   var_Muon1_IsPFMuon = Ntp->Muon_isPFMuon(Muon_index_1);
   var_Muon1_IsRPCMuon = Ntp->Muon_isRPCMuon(Muon_index_1);
   var_Muon1_emEt03 = Ntp->Muon_emEt03(Muon_index_1);
   var_Muon1_emVetoEt03 = Ntp->Muon_emVetoEt03(Muon_index_1);
   var_Muon1_hadEt03 = Ntp->Muon_hadEt03(Muon_index_1);
   var_Muon1_hadVetoEt03 = Ntp->Muon_hadVetoEt03(Muon_index_1);
   var_Muon1_nJets03 = Ntp->Muon_nJets03(Muon_index_1);
   var_Muon1_nTracks03 = Ntp->Muon_nTracks03(Muon_index_1);
   var_Muon1_StandardSelection = Ntp->Muon_StandardSelection(Muon_index_1);
   var_Muon1_sumPt03 = Ntp->Muon_sumPt03(Muon_index_1);
   var_Muon1_trackerVetoPt03 = Ntp->Muon_trackerVetoPt03(Muon_index_1);
   var_Muon1_sumChargedHadronPt03 = Ntp->Muon_sumChargedHadronPt03(Muon_index_1);
   var_Muon1_sumChargedParticlePt03 = Ntp->Muon_sumChargedParticlePt03(Muon_index_1);
   var_Muon1_sumNeutralHadronEt03 = Ntp->Muon_sumNeutralHadronEt03(Muon_index_1);
   var_Muon1_sumNeutralHadronEtHighThreshold03 = Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_1);
   var_Muon1_sumPhotonEt03 = Ntp->Muon_sumPhotonEt03(Muon_index_1);
   var_Muon1_sumPhotonEtHighThreshold03 = Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_1);
   var_Muon1_sumPUPt03 = Ntp->Muon_sumPUPt03(Muon_index_1);
   var_Muon1_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_1);
   var_Muon1_Track_idx = Ntp->Muon_Track_idx(Muon_index_1);
   var_Muon1_combinedQuality_updatedSta = Ntp->Muon_combinedQuality_updatedSta(Muon_index_1);
   var_Muon1_combinedQuality_trkKink = Ntp->Muon_combinedQuality_trkKink(Muon_index_1);
   var_Muon1_combinedQuality_glbKink = Ntp->Muon_combinedQuality_glbKink(Muon_index_1);
   var_Muon1_combinedQuality_trkRelChi2 = Ntp->Muon_combinedQuality_trkRelChi2(Muon_index_1);
   var_Muon1_combinedQuality_staRelChi2 = Ntp->Muon_combinedQuality_staRelChi2(Muon_index_1);
   var_Muon1_combinedQuality_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1);
   var_Muon1_combinedQuality_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1);
   var_Muon1_combinedQuality_localDistance = Ntp->Muon_combinedQuality_localDistance(Muon_index_1);
   var_Muon1_combinedQuality_globalDeltaEtaPhi = Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_1);
   var_Muon1_combinedQuality_tightMatch = Ntp->Muon_combinedQuality_tightMatch(Muon_index_1);
   var_Muon1_combinedQuality_glbTrackProbability = Ntp->Muon_combinedQuality_glbTrackProbability(Muon_index_1);
   var_Muon1_prod_inner_outer_charge = Ntp->Muon_prod_inner_outer_charge(Muon_index_1);
   var_Muon1_innerTrack_quality = Ntp->Muon_innerTrack_quality(Muon_index_1);
   var_Muon1_ptErrOverPt = Ntp->Muon_ptErrOverPt(Muon_index_1);
   var_Muon1_calEnergy_hadS9 = Ntp->Muon_calEnergy_hadS9(Muon_index_1);
   var_Muon1_calEnergy_had = Ntp->Muon_calEnergy_had(Muon_index_1);
   var_Muon1_calEnergy_emS25 = Ntp->Muon_calEnergy_emS25(Muon_index_1);
   var_Muon1_calEnergy_emS9 = Ntp->Muon_calEnergy_emS9(Muon_index_1);
   var_Muon1_calEnergy_em = Ntp->Muon_calEnergy_em(Muon_index_1);
   var_Muon1_charge = Ntp->Muon_charge(Muon_index_1);
   var_Muon1_trackCharge = Ntp->Muon_trackCharge(Muon_index_1);
   var_Muon1_hitPattern_pixelLayerwithMeas = Ntp->Muon_hitPattern_pixelLayerwithMeas(Muon_index_1);
   var_Muon1_numberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_1);
   var_Muon1_normChi2 = Ntp->Muon_normChi2(Muon_index_1);
   var_Muon1_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1);
   var_Muon1_innerTrack_numberofValidHits = Ntp->Muon_innerTrack_numberofValidHits(Muon_index_1);
   var_Muon1_numberofValidPixelHits = Ntp->Muon_numberofValidPixelHits(Muon_index_1);
   var_Muon1_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_1);
   var_Muon1_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(Muon_index_1);
   var_Muon1_segmentCompatibility = Ntp->Muon_segmentCompatibility(Muon_index_1);
   var_Muon1_caloCompatibility = Ntp->Muon_caloCompatibility(Muon_index_1);
   var_Muon1_innerTrack_validFraction = Ntp->Muon_innerTrack_validFraction(Muon_index_1);
   var_Muon1_innerTrack_pixelLayersWithMeasurement = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(Muon_index_1);
   var_Muon1_innerTrack_numberOfValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(Muon_index_1);
   var_Muon1_innerTrack_numberOfLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(Muon_index_1);
   var_Muon1_innerTrack_numberOfLostTrackerInnerHits = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(Muon_index_1);
   var_Muon1_innerTrack_numberOfLostTrackerOuterHits = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(Muon_index_1);
   var_Muon1_innerTrack_normalizedChi2 = Ntp->Muon_innerTrack_normalizedChi2(Muon_index_1);
   var_Muon1_vmuonhitcomb_reco = Ntp->Muon_vmuonhitcomb_reco(Muon_index_1);
   var_Muon1_rpchits_reco = Ntp->Muon_rpchits_reco(Muon_index_1);
   var_Muon1_outerTrack_normalizedChi2 = Ntp->Muon_outerTrack_normalizedChi2(Muon_index_1);
   var_Muon1_outerTrack_muonStationsWithValidHits = Ntp->Muon_outerTrack_muonStationsWithValidHits(Muon_index_1);
   var_Muon1_isGoodMuon_TM2DCompatibility = Ntp->Muon_isGoodMuon_TM2DCompatibility(Muon_index_1);
   var_Muon1_isGoodMuon_TrackerMuonArbitrated = Ntp->Muon_isGoodMuon_TrackerMuonArbitrated(Muon_index_1);
   var_Muon1_isGoodMuon_TMOneStationTight = Ntp->Muon_isGoodMuon_TMOneStationTight(Muon_index_1);
   var_Muon1_isGoodMuon_TMOneStationAngTight = Ntp->Muon_isGoodMuon_TMOneStationAngTight(Muon_index_1);
   var_Muon1_isGoodMuon_TMLastStationTight = Ntp->Muon_isGoodMuon_TMLastStationTight(Muon_index_1);
   var_Muon1_isGoodMuon_TMLastStationAngTight = Ntp->Muon_isGoodMuon_TMLastStationAngTight(Muon_index_1);
   var_Muon1_isGoodMuon_TMLastStationOptimizedLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(Muon_index_1);
   var_Muon1_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(Muon_index_1);

   var_Muon2_vx = Ntp->Muon_Poca(Muon_index_2).X();
   var_Muon2_vy = Ntp->Muon_Poca(Muon_index_2).Y();
   var_Muon2_vz = Ntp->Muon_Poca(Muon_index_2).Z();
   var_Muon2_IsGlobalMuon = Ntp->Muon_isGlobalMuon(Muon_index_2);
   var_Muon2_IsStandAloneMuon = Ntp->Muon_isStandAloneMuon(Muon_index_2);
   var_Muon2_IsTrackerMuon = Ntp->Muon_isTrackerMuon(Muon_index_2);
   var_Muon2_IsCaloMuon = Ntp->Muon_isCaloMuon(Muon_index_2);
   var_Muon2_IsIsolationValid = Ntp->Muon_isIsolationValid(Muon_index_2);
   var_Muon2_IsQualityValid = Ntp->Muon_isQualityValid(Muon_index_2);
   var_Muon2_IsTimeValid = Ntp->Muon_isTimeValid(Muon_index_2);
   var_Muon2_IsPFMuon = Ntp->Muon_isPFMuon(Muon_index_2);
   var_Muon2_IsRPCMuon = Ntp->Muon_isRPCMuon(Muon_index_2);
   var_Muon2_emEt03 = Ntp->Muon_emEt03(Muon_index_2);
   var_Muon2_emVetoEt03 = Ntp->Muon_emVetoEt03(Muon_index_2);
   var_Muon2_hadEt03 = Ntp->Muon_hadEt03(Muon_index_2);
   var_Muon2_hadVetoEt03 = Ntp->Muon_hadVetoEt03(Muon_index_2);
   var_Muon2_nJets03 = Ntp->Muon_nJets03(Muon_index_2);
   var_Muon2_nTracks03 = Ntp->Muon_nTracks03(Muon_index_2);
   var_Muon2_StandardSelection = Ntp->Muon_StandardSelection(Muon_index_2);
   var_Muon2_sumPt03 = Ntp->Muon_sumPt03(Muon_index_2);
   var_Muon2_trackerVetoPt03 = Ntp->Muon_trackerVetoPt03(Muon_index_2);
   var_Muon2_sumChargedHadronPt03 = Ntp->Muon_sumChargedHadronPt03(Muon_index_2);
   var_Muon2_sumChargedParticlePt03 = Ntp->Muon_sumChargedParticlePt03(Muon_index_2);
   var_Muon2_sumNeutralHadronEt03 = Ntp->Muon_sumNeutralHadronEt03(Muon_index_2);
   var_Muon2_sumNeutralHadronEtHighThreshold03 = Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_2);
   var_Muon2_sumPhotonEt03 = Ntp->Muon_sumPhotonEt03(Muon_index_2);
   var_Muon2_sumPhotonEtHighThreshold03 = Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_2);
   var_Muon2_sumPUPt03 = Ntp->Muon_sumPUPt03(Muon_index_2);
   var_Muon2_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_2);
   var_Muon2_Track_idx = Ntp->Muon_Track_idx(Muon_index_2);
   var_Muon2_combinedQuality_updatedSta = Ntp->Muon_combinedQuality_updatedSta(Muon_index_2);
   var_Muon2_combinedQuality_trkKink = Ntp->Muon_combinedQuality_trkKink(Muon_index_2);
   var_Muon2_combinedQuality_glbKink = Ntp->Muon_combinedQuality_glbKink(Muon_index_2);
   var_Muon2_combinedQuality_trkRelChi2 = Ntp->Muon_combinedQuality_trkRelChi2(Muon_index_2);
   var_Muon2_combinedQuality_staRelChi2 = Ntp->Muon_combinedQuality_staRelChi2(Muon_index_2);
   var_Muon2_combinedQuality_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2);
   var_Muon2_combinedQuality_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2);
   var_Muon2_combinedQuality_localDistance = Ntp->Muon_combinedQuality_localDistance(Muon_index_2);
   var_Muon2_combinedQuality_globalDeltaEtaPhi = Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_2);
   var_Muon2_combinedQuality_tightMatch = Ntp->Muon_combinedQuality_tightMatch(Muon_index_2);
   var_Muon2_combinedQuality_glbTrackProbability = Ntp->Muon_combinedQuality_glbTrackProbability(Muon_index_2);
   var_Muon2_prod_inner_outer_charge = Ntp->Muon_prod_inner_outer_charge(Muon_index_2);
   var_Muon2_innerTrack_quality = Ntp->Muon_innerTrack_quality(Muon_index_2);
   var_Muon2_ptErrOverPt = Ntp->Muon_ptErrOverPt(Muon_index_2);
   var_Muon2_calEnergy_hadS9 = Ntp->Muon_calEnergy_hadS9(Muon_index_2);
   var_Muon2_calEnergy_had = Ntp->Muon_calEnergy_had(Muon_index_2);
   var_Muon2_calEnergy_emS25 = Ntp->Muon_calEnergy_emS25(Muon_index_2);
   var_Muon2_calEnergy_emS9 = Ntp->Muon_calEnergy_emS9(Muon_index_2);
   var_Muon2_calEnergy_em = Ntp->Muon_calEnergy_em(Muon_index_2);
   var_Muon2_charge = Ntp->Muon_charge(Muon_index_2);
   var_Muon2_trackCharge = Ntp->Muon_trackCharge(Muon_index_2);
   var_Muon2_hitPattern_pixelLayerwithMeas = Ntp->Muon_hitPattern_pixelLayerwithMeas(Muon_index_2);
   var_Muon2_numberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_2);
   var_Muon2_normChi2 = Ntp->Muon_normChi2(Muon_index_2);
   var_Muon2_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2);
   var_Muon2_innerTrack_numberofValidHits = Ntp->Muon_innerTrack_numberofValidHits(Muon_index_2);
   var_Muon2_numberofValidPixelHits = Ntp->Muon_numberofValidPixelHits(Muon_index_2);
   var_Muon2_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_2);
   var_Muon2_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(Muon_index_2);
   var_Muon2_segmentCompatibility = Ntp->Muon_segmentCompatibility(Muon_index_2);
   var_Muon2_caloCompatibility = Ntp->Muon_caloCompatibility(Muon_index_2);
   var_Muon2_innerTrack_validFraction = Ntp->Muon_innerTrack_validFraction(Muon_index_2);
   var_Muon2_innerTrack_pixelLayersWithMeasurement = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(Muon_index_2);
   var_Muon2_innerTrack_numberOfValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(Muon_index_2);
   var_Muon2_innerTrack_numberOfLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(Muon_index_2);
   var_Muon2_innerTrack_numberOfLostTrackerInnerHits = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(Muon_index_2);
   var_Muon2_innerTrack_numberOfLostTrackerOuterHits = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(Muon_index_2);
   var_Muon2_innerTrack_normalizedChi2 = Ntp->Muon_innerTrack_normalizedChi2(Muon_index_2);
   var_Muon2_vmuonhitcomb_reco = Ntp->Muon_vmuonhitcomb_reco(Muon_index_2);
   var_Muon2_rpchits_reco = Ntp->Muon_rpchits_reco(Muon_index_2);
   var_Muon2_outerTrack_normalizedChi2 = Ntp->Muon_outerTrack_normalizedChi2(Muon_index_2);
   var_Muon2_outerTrack_muonStationsWithValidHits = Ntp->Muon_outerTrack_muonStationsWithValidHits(Muon_index_2);
   var_Muon2_isGoodMuon_TM2DCompatibility = Ntp->Muon_isGoodMuon_TM2DCompatibility(Muon_index_2);
   var_Muon2_isGoodMuon_TrackerMuonArbitrated = Ntp->Muon_isGoodMuon_TrackerMuonArbitrated(Muon_index_2);
   var_Muon2_isGoodMuon_TMOneStationTight = Ntp->Muon_isGoodMuon_TMOneStationTight(Muon_index_2);
   var_Muon2_isGoodMuon_TMOneStationAngTight = Ntp->Muon_isGoodMuon_TMOneStationAngTight(Muon_index_2);
   var_Muon2_isGoodMuon_TMLastStationTight = Ntp->Muon_isGoodMuon_TMLastStationTight(Muon_index_2);
   var_Muon2_isGoodMuon_TMLastStationAngTight = Ntp->Muon_isGoodMuon_TMLastStationAngTight(Muon_index_2);
   var_Muon2_isGoodMuon_TMLastStationOptimizedLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(Muon_index_2);
   var_Muon2_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(Muon_index_2);

   var_Muon3_vx = Ntp->Muon_Poca(Muon_index_3).X();
   var_Muon3_vy = Ntp->Muon_Poca(Muon_index_3).Y();
   var_Muon3_vz = Ntp->Muon_Poca(Muon_index_3).Z();
   var_Muon3_IsGlobalMuon = Ntp->Muon_isGlobalMuon(Muon_index_3);
   var_Muon3_IsStandAloneMuon = Ntp->Muon_isStandAloneMuon(Muon_index_3);
   var_Muon3_IsTrackerMuon = Ntp->Muon_isTrackerMuon(Muon_index_3);
   var_Muon3_IsCaloMuon = Ntp->Muon_isCaloMuon(Muon_index_3);
   var_Muon3_IsIsolationValid = Ntp->Muon_isIsolationValid(Muon_index_3);
   var_Muon3_IsQualityValid = Ntp->Muon_isQualityValid(Muon_index_3);
   var_Muon3_IsTimeValid = Ntp->Muon_isTimeValid(Muon_index_3);
   var_Muon3_IsPFMuon = Ntp->Muon_isPFMuon(Muon_index_3);
   var_Muon3_IsRPCMuon = Ntp->Muon_isRPCMuon(Muon_index_3);
   var_Muon3_emEt03 = Ntp->Muon_emEt03(Muon_index_3);
   var_Muon3_emVetoEt03 = Ntp->Muon_emVetoEt03(Muon_index_3);
   var_Muon3_hadEt03 = Ntp->Muon_hadEt03(Muon_index_3);
   var_Muon3_hadVetoEt03 = Ntp->Muon_hadVetoEt03(Muon_index_3);
   var_Muon3_nJets03 = Ntp->Muon_nJets03(Muon_index_3);
   var_Muon3_nTracks03 = Ntp->Muon_nTracks03(Muon_index_3);
   var_Muon3_StandardSelection = Ntp->Muon_StandardSelection(Muon_index_3);
   var_Muon3_sumPt03 = Ntp->Muon_sumPt03(Muon_index_3);
   var_Muon3_trackerVetoPt03 = Ntp->Muon_trackerVetoPt03(Muon_index_3);
   var_Muon3_sumChargedHadronPt03 = Ntp->Muon_sumChargedHadronPt03(Muon_index_3);
   var_Muon3_sumChargedParticlePt03 = Ntp->Muon_sumChargedParticlePt03(Muon_index_3);
   var_Muon3_sumNeutralHadronEt03 = Ntp->Muon_sumNeutralHadronEt03(Muon_index_3);
   var_Muon3_sumNeutralHadronEtHighThreshold03 = Ntp->Muon_sumNeutralHadronEtHighThreshold03(Muon_index_3);
   var_Muon3_sumPhotonEt03 = Ntp->Muon_sumPhotonEt03(Muon_index_3);
   var_Muon3_sumPhotonEtHighThreshold03 = Ntp->Muon_sumPhotonEtHighThreshold03(Muon_index_3);
   var_Muon3_sumPUPt03 = Ntp->Muon_sumPUPt03(Muon_index_3);
   var_Muon3_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_3);
   var_Muon3_Track_idx = Ntp->Muon_Track_idx(Muon_index_3);
   var_Muon3_combinedQuality_updatedSta = Ntp->Muon_combinedQuality_updatedSta(Muon_index_3);
   var_Muon3_combinedQuality_trkKink = Ntp->Muon_combinedQuality_trkKink(Muon_index_3);
   var_Muon3_combinedQuality_glbKink = Ntp->Muon_combinedQuality_glbKink(Muon_index_3);
   var_Muon3_combinedQuality_trkRelChi2 = Ntp->Muon_combinedQuality_trkRelChi2(Muon_index_3);
   var_Muon3_combinedQuality_staRelChi2 = Ntp->Muon_combinedQuality_staRelChi2(Muon_index_3);
   var_Muon3_combinedQuality_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3);
   var_Muon3_combinedQuality_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3);
   var_Muon3_combinedQuality_localDistance = Ntp->Muon_combinedQuality_localDistance(Muon_index_3);
   var_Muon3_combinedQuality_globalDeltaEtaPhi = Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Muon_index_3);
   var_Muon3_combinedQuality_tightMatch = Ntp->Muon_combinedQuality_tightMatch(Muon_index_3);
   var_Muon3_combinedQuality_glbTrackProbability = Ntp->Muon_combinedQuality_glbTrackProbability(Muon_index_3);
   var_Muon3_prod_inner_outer_charge = Ntp->Muon_prod_inner_outer_charge(Muon_index_3);
   var_Muon3_innerTrack_quality = Ntp->Muon_innerTrack_quality(Muon_index_3);
   var_Muon3_ptErrOverPt = Ntp->Muon_ptErrOverPt(Muon_index_3);
   var_Muon3_calEnergy_hadS9 = Ntp->Muon_calEnergy_hadS9(Muon_index_3);
   var_Muon3_calEnergy_had = Ntp->Muon_calEnergy_had(Muon_index_3);
   var_Muon3_calEnergy_emS25 = Ntp->Muon_calEnergy_emS25(Muon_index_3);
   var_Muon3_calEnergy_emS9 = Ntp->Muon_calEnergy_emS9(Muon_index_3);
   var_Muon3_calEnergy_em = Ntp->Muon_calEnergy_em(Muon_index_3);
   var_Muon3_charge = Ntp->Muon_charge(Muon_index_3);
   var_Muon3_trackCharge = Ntp->Muon_trackCharge(Muon_index_3);
   var_Muon3_hitPattern_pixelLayerwithMeas = Ntp->Muon_hitPattern_pixelLayerwithMeas(Muon_index_3);
   var_Muon3_numberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_3);
   var_Muon3_normChi2 = Ntp->Muon_normChi2(Muon_index_3);
   var_Muon3_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3);
   var_Muon3_innerTrack_numberofValidHits = Ntp->Muon_innerTrack_numberofValidHits(Muon_index_3);
   var_Muon3_numberofValidPixelHits = Ntp->Muon_numberofValidPixelHits(Muon_index_3);
   var_Muon3_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_3);
   var_Muon3_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(Muon_index_3);
   var_Muon3_segmentCompatibility = Ntp->Muon_segmentCompatibility(Muon_index_3);
   var_Muon3_caloCompatibility = Ntp->Muon_caloCompatibility(Muon_index_3);
   var_Muon3_innerTrack_validFraction = Ntp->Muon_innerTrack_validFraction(Muon_index_3);
   var_Muon3_innerTrack_pixelLayersWithMeasurement = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(Muon_index_3);
   var_Muon3_innerTrack_numberOfValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(Muon_index_3);
   var_Muon3_innerTrack_numberOfLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(Muon_index_3);
   var_Muon3_innerTrack_numberOfLostTrackerInnerHits = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(Muon_index_3);
   var_Muon3_innerTrack_numberOfLostTrackerOuterHits = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(Muon_index_3);
   var_Muon3_innerTrack_normalizedChi2 = Ntp->Muon_innerTrack_normalizedChi2(Muon_index_3);
   var_Muon3_vmuonhitcomb_reco = Ntp->Muon_vmuonhitcomb_reco(Muon_index_3);
   var_Muon3_rpchits_reco = Ntp->Muon_rpchits_reco(Muon_index_3);
   var_Muon3_outerTrack_normalizedChi2 = Ntp->Muon_outerTrack_normalizedChi2(Muon_index_3);
   var_Muon3_outerTrack_muonStationsWithValidHits = Ntp->Muon_outerTrack_muonStationsWithValidHits(Muon_index_3);
   var_Muon3_isGoodMuon_TM2DCompatibility = Ntp->Muon_isGoodMuon_TM2DCompatibility(Muon_index_3);
   var_Muon3_isGoodMuon_TrackerMuonArbitrated = Ntp->Muon_isGoodMuon_TrackerMuonArbitrated(Muon_index_3);
   var_Muon3_isGoodMuon_TMOneStationTight = Ntp->Muon_isGoodMuon_TMOneStationTight(Muon_index_3);
   var_Muon3_isGoodMuon_TMOneStationAngTight = Ntp->Muon_isGoodMuon_TMOneStationAngTight(Muon_index_3);
   var_Muon3_isGoodMuon_TMLastStationTight = Ntp->Muon_isGoodMuon_TMLastStationTight(Muon_index_3);
   var_Muon3_isGoodMuon_TMLastStationAngTight = Ntp->Muon_isGoodMuon_TMLastStationAngTight(Muon_index_3);
   var_Muon3_isGoodMuon_TMLastStationOptimizedLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(Muon_index_3);
   var_Muon3_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(Muon_index_3);

   // muon time variabels
   var_Muon1_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(Muon_index_1);
   var_Muon1_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(Muon_index_1);
   var_Muon2_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(Muon_index_2);
   var_Muon2_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(Muon_index_2);
   var_Muon3_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(Muon_index_3);
   var_Muon3_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(Muon_index_3);


   // Dimuon variables
   var_NMuTrkOS1 = NMuTrkOS1;
   var_NMuTrkOS2 = NMuTrkOS2;
   var_NMuTrkOS3 = NMuTrkOS3;

   // Dimuon variables
   var_Mosss1 = M_osss1;
   var_Mosss2 = M_osss2;

   // standard variables
   var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx,false);
   var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
   var_flightLenSig =  Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,false),
         Ntp->Vertex_PrimaryVertex_Covariance(final_idx,false),
         Ntp->Vertex_Signal_KF_pos(final_idx,false),
         Ntp->Vertex_Signal_KF_Covariance(final_idx,false));

   var_tauMass = TauLV.M();
   var_ntracks = Ntp->Isolation_NTracks(final_idx,false);
   var_relPt = Ntp->Isolation_RelPt(final_idx,false);
   var_isoMax = Ntp->Isolation_maxdy(final_idx,false);

   double maxMuondR = std::max({Muon1LV.DeltaR(TauLV), Muon2LV.DeltaR(TauLV), Muon3LV.DeltaR(TauLV)});
   double minMuonPt = std::min({Muon1LV.Pt(), Muon2LV.Pt(), Muon3LV.Pt()});

   float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,false),
         Ntp->Vertex_d0sig_reco(final_idx,1,false),
         Ntp->Vertex_d0sig_reco(final_idx,2,false)});
   float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0,false),
         Ntp->Vertex_d0sig_reco(final_idx,1,false),
         Ntp->Vertex_d0sig_reco(final_idx,2,false)});

   var_nsv = NumberOfPrimaryVertices;

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

   var_minMatchedStations = std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});
   var_pmin = std::min(Ntp->Muon_P4(Muon_index_1).P(), std::min(Ntp->Muon_P4(Muon_index_2).P(), Ntp->Muon_P4(Muon_index_3).P()));

   // Muon ID variables
   var_max_cLP = std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1), std::max(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2), Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)));
   var_max_tKink = std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_1), std::max(Ntp->Muon_combinedQuality_trkKink(Muon_index_2), Ntp->Muon_combinedQuality_trkKink(Muon_index_3)));
   var_segCompMuMin  = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
   var_MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
   var_sumMuTrkKinkChi2= (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));

   // vertex variables
   var_maxdca = std::max({Ntp->Vertex_DCA12(final_idx, false),Ntp->Vertex_DCA23(final_idx, false),Ntp->Vertex_DCA31(final_idx, false)});
   var_MuMu_minKFChi2 = std::min({Ntp->Vertex_pair_quality(final_idx, 0, false), Ntp->Vertex_pair_quality(final_idx, 1, false), Ntp->Vertex_pair_quality(final_idx, 2, false)});

   // Isolation variables
   var_MinD0Significance = MinD0Significance;
   var_MaxD0Significance = MaxD0Significance;
   var_mindca_iso = mindca_iso;
   var_trk_relPt = std::max({mu1_relPt, mu2_relPt, mu3_relPt});
   var_nMatches_mu3 = Ntp->Muon_numberOfMatches(Muon_index_3);

   var_VertexMu1D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,0,false);
   var_VertexMu2D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,1,false);
   var_VertexMu3D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,2,false);

   var_VertexMu1D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0,false);
   var_VertexMu2D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1,false);
   var_VertexMu3D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2,false);

   var_VertexMu1D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,0,false);
   var_VertexMu2D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,1,false);
   var_VertexMu3D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,2,false);

   // Isolation
   var_InvMassConeIso0p3 = sumTracksIso0p3.M();
   var_InvMassConeIso0p5 = sumTracksIso0p5.M();
   var_InvMassConeIso0p8 = sumTracksIso0p8.M();
   var_InvMassConeIso1p2 = sumTracksIso1p2.M();
   var_InvMassConeIso1p4 = sumTracksIso1p4.M();

   var_NTrksIso0p3 = nTracksIso0p3;
   var_NTrksIso0p5 = nTracksIso0p5;
   var_NTrksIso0p8 = nTracksIso0p8;
   var_NTrksIso1p2 = nTracksIso1p2;
   var_NTrksIso1p4 = nTracksIso1p4;


   //---------------------
   var_NtracksClose = NcloseTracksCount;

   if(Ntp->NIsolationTrack(final_idx,false)!=0){ 
      var_MindcaTrackSV = sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex,false),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex,false),2));
   } else var_MindcaTrackSV = -1;

   if(Ntp->NIsolationTrack(final_idx,false)!=0){ 
      var_dcaTrackPV = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV,false),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV,false),2));
   } else var_dcaTrackPV = -1;

   if (id==1) MC=0;
   else  MC=1;

   if (tauMassRes<tauMassResCutLow) category = 1;
   if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh) category = 2;
   if (tauMassRes>tauMassResCutHigh) category = 3;

   TMVA_Tree->Fill();
}
}

void  T3MSelectionTree::Finish(){


   TMVA_Tree->Write();
   file->Close();

   Selection::Finish();
}

double T3MSelectionTree::GetReferenceFrameAngle(TLorentzVector vec1, TLorentzVector vec2, TLorentzVector refVector){
   TVector3 boostVector = -(vec1+vec2).BoostVector();
   vec1.Boost(boostVector);
   vec2.Boost(boostVector);
   refVector.Boost(boostVector);
   return vec1.Vect().Angle(refVector.Vect());
}

double T3MSelectionTree::GetReferenceFrameCosine(TLorentzVector vec1, TLorentzVector vec2, TLorentzVector refVector){
   TVector3 boostVector = -(vec1+vec2).BoostVector();
   vec1.Boost(boostVector);
   vec2.Boost(boostVector);
   refVector.Boost(boostVector);
   return TMath::Cos(vec1.Vect().Angle(refVector.Vect()));
}

template <typename T>
int T3MSelectionTree::minQuantityIndex(std::vector<T>& vec){
   if (vec.at(0)<=vec.at(1) && vec.at(0)<=vec.at(2)) return 0;
   if (vec.at(1)<=vec.at(2) && vec.at(1)<=vec.at(0)) return 1;
   if (vec.at(2)<=vec.at(0) && vec.at(2)<=vec.at(1)) return 2;
   return -1;
}

template <typename T>
int T3MSelectionTree::maxQuantityIndex(std::vector<T>& vec){
   if (vec.at(0)>=vec.at(1) && vec.at(0)>=vec.at(2)) return 0;
   if (vec.at(1)>=vec.at(2) && vec.at(1)>=vec.at(0)) return 1;
   if (vec.at(2)>=vec.at(0) && vec.at(2)>=vec.at(1)) return 2;
   return -1;
}

