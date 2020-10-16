void makeWeighedTree(TString fname, TString path_to_weight_file){

	  TFile* weight_file = new TFile(path_to_weight_file, "READ");
	  TH2F* weights = (TH2F*)weight_file->Get("weights");

	  TFile* file = new TFile(fname, "update");
	  TTree* tree = (TTree*)file->Get("tree");

	  float pt, eta, phi;
	  int ptBin, etaBin, phiBin;

	  float eventWeight;
	  int fake;
	  float muonPt ;
	  float muonEta ;
	  float muonPhi ;
	  float muonInnerNC2 ;
	  float muonValidFraction ;
	  int muonInnerNValidHits ;
	  float muonInnerTrackQuality ;
	  int muonNValidPixelHits ;
	  int muonNValidTrackerHits ;
	  int muonNLostTrackerHits ;
	  int muonNLostTrackerHitsInner ;
	  int muonNLostTrackerHitsOuter ;
	  int muonPixelLayers ;
	  int muonTrackerLayers ;
	  int muonNMatchedStations ;
	  int muonNMatches ;
	  bool muonRPC ;
	  float muonPtErrPt ;
	  float muonSegComp ;
	  float muonCaloComp ;
	  float muonHadS9 ;
	  float muonHad ;
	  float muonEM ;
	  float muonEMS9 ;
	  float muonEMS25 ;
	  float muonKink;

	  tree->SetBranchAddress("fake",&fake);
	  tree->SetBranchAddress("muonPt",&muonPt);
	  tree->SetBranchAddress("muonEta",&muonEta);
	  tree->SetBranchAddress("muonPhi",&muonPhi);
	  tree->SetBranchAddress("muonInnerNC2",&muonInnerNC2);
	  tree->SetBranchAddress("muonValidFraction",&muonValidFraction);
	  tree->SetBranchAddress("muonInnerNValidHits",&muonInnerNValidHits);
	  tree->SetBranchAddress("muonInnerTrackQuality",&muonInnerTrackQuality);
	  tree->SetBranchAddress("muonNValidPixelHits",&muonNValidPixelHits);
	  tree->SetBranchAddress("muonNValidTrackerHits",&muonNValidTrackerHits);
	  tree->SetBranchAddress("muonNLostTrackerHits",&muonNLostTrackerHits);
	  tree->SetBranchAddress("muonNLostTrackerHitsInner",&muonNLostTrackerHitsInner);
	  tree->SetBranchAddress("muonNLostTrackerHitsOuter",&muonNLostTrackerHitsOuter);
	  tree->SetBranchAddress("muonPixelLayers",&muonPixelLayers);
	  tree->SetBranchAddress("muonTrackerLayers",&muonTrackerLayers);
	  tree->SetBranchAddress("muonNMatchedStations",&muonNMatchedStations);
	  tree->SetBranchAddress("muonNMatches",&muonNMatches);
	  tree->SetBranchAddress("muonRPC",&muonRPC);
	  tree->SetBranchAddress("muonPtErrPt",&muonPtErrPt);
	  tree->SetBranchAddress("muonSegComp",&muonSegComp);
	  tree->SetBranchAddress("muonCaloComp",&muonCaloComp);
	  tree->SetBranchAddress("muonHadS9",&muonHadS9);
	  tree->SetBranchAddress("muonHad",&muonHad);
	  tree->SetBranchAddress("muonEM",&muonEM);
	  tree->SetBranchAddress("muonEMS9",&muonEMS9);
	  tree->SetBranchAddress("muonEMS25",&muonEMS25);
	  tree->SetBranchAddress("muonKink",&muonKink);

	  TBranch *br = (TBranch*)tree->Branch("eventWeight",&eventWeight);

	  for (int i=0; i<tree->GetEntries(); i++){
				 tree->GetEntry(i);
				 if (!fake) eventWeight = 1.0;
				 else{
							pt = muonPt;
							eta = muonEta;
							phi = muonPhi;
							TAxis* xaxis = weights->GetXaxis();
							TAxis* yaxis = weights->GetYaxis();
							Int_t xbin = xaxis->FindBin(pt);
							Int_t ybin = yaxis->FindBin(eta);
							eventWeight = weights->GetBinContent(xbin, ybin);
							//cout<<"Entry: "<<i<<", weight: "<<eventWeight<<endl;
				 }
				 br->Fill();
	  }

	  tree->Write();
	  delete file;
}
