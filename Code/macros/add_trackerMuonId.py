from array import array
import ROOT
from ROOT import TMVA

Muon3_Pt = array('f', [0.])
Muon3_Phi = array('f', [0.])
Muon3_Eta = array('f', [0.])

Muon3_caloCompatibility = array('f', [0.])
Muon3_segmentCompatibility = array('f', [0.])
Muon3_innerTrack_numberOfLostTrackerHits = array('f', [0.])
Muon3_innerTrack_numberOfLostTrackerInnerHits = array('f', [0.])
Muon3_numberOfMatchedStations = array('f', [0.])
Muon3_calEnergy_had = array('f', [0.])
Muon3_innerTrack_numberofValidHits = array('f', [0.])
Muon3_calEnergy_em = array('f', [0.])
Muon3_innerTrack_numberOfLostTrackerOuterHits = array('f', [0.])
Muon3_innerTrack_normalizedChi2 = array('f', [0.])
Muon3_innerTrack_validFraction = array('f', [0.])
Muon3_innerTrack_pixelLayersWithMeasurement = array('f', [0.])
Muon3_ptErrOverPt = array('f', [0.])

fake = array('f', [0.]);

MuonId = array('f', [-99.0]);

readerMuonId = TMVA.Reader("!Color:!Silent");
readerMuonId.AddVariable("Muon_caloCompatibility", Muon3_caloCompatibility);
readerMuonId.AddVariable("Muon_segmentCompatibility", Muon3_segmentCompatibility);
readerMuonId.AddVariable("Muon_innerTrack_numberOfLostTrackerHits", Muon3_innerTrack_numberOfLostTrackerHits);
readerMuonId.AddVariable("Muon_innerTrack_numberOfLostTrackerInnerHits", Muon3_innerTrack_numberOfLostTrackerInnerHits);
readerMuonId.AddVariable("Muon_numberOfMatchedStations", Muon3_numberOfMatchedStations);
readerMuonId.AddVariable("Muon_calEnergy_had", Muon3_calEnergy_had);
readerMuonId.AddVariable("Muon_innerTrack_numberofValidHits", Muon3_innerTrack_numberofValidHits);
readerMuonId.AddVariable("Muon_calEnergy_em", Muon3_calEnergy_em);
readerMuonId.AddVariable("Muon_innerTrack_numberOfLostTrackerOuterHits", Muon3_innerTrack_numberOfLostTrackerOuterHits);
readerMuonId.AddVariable("Muon_innerTrack_normalizedChi2", Muon3_innerTrack_normalizedChi2);
readerMuonId.AddVariable("Muon_innerTrack_validFraction", Muon3_innerTrack_validFraction);
readerMuonId.AddVariable("Muon_innerTrack_pixelLayersWithMeasurement", Muon3_innerTrack_pixelLayersWithMeasurement);
readerMuonId.AddVariable("Muon_ptErrOverPt", Muon3_ptErrOverPt);
readerMuonId.AddSpectator("Muon_Eta", Muon3_Eta);
readerMuonId.AddSpectator("Muon_Phi", Muon3_Phi);
readerMuonId.AddSpectator("Muon_Pt", Muon3_Pt);
readerMuonId.AddSpectator("fake", fake);
readerMuonId.BookMVA( "BDT", "../CommonFiles/weights/TrackerMuonClassification_2018_standard_variables/\
weights/MuPiTMVA_mod_2018_BDT.weights.xml"); #weights weights.xml file after training, place it to CommonFiles


file_ = ROOT.TFile("/eos/user/b/bjoshi/RunIITau23Mu/AnalysisTrees/T3MSelectionTree_veto_optimized_dnn_globalMuonId.root","READ");
old_tree_list = [ file_.Get("TreeB"), file_.Get('TreeS_Ds'), file_.Get('TreeS_Bu'), file_.Get('TreeS_Bd')];
new_file = ROOT.TFile("/eos/user/b/bjoshi/RunIITau23Mu/AnalysisTrees/T3MSelectionTree_veto_optimized_dnn_globalMuonId_trackerMuonId_v2.root","RECREATE");
new_tree_list = []
for i in xrange(4): new_tree_list.append(old_tree_list[i].CloneTree(0));

for i in xrange(4):
   new_tree_list[i].Branch("var_trackerMuonId_v2", MuonId, "var_trackerMuonId/F");


for i in xrange(4):
   tree = old_tree_list[i] 
   print tree.GetName()
   nentries = tree.GetEntriesFast();
   for j in xrange(nentries):
      old_tree_list[i].GetEntry(j);
      if (j%100000==0): print "Processing ",j,"/",nentries," ..."
      Muon3_caloCompatibility[0] = float(getattr(old_tree_list[i], 'Muon3_caloCompatibility'));
      Muon3_segmentCompatibility[0] = float(getattr(old_tree_list[i], 'Muon3_segmentCompatibility'));
      Muon3_innerTrack_numberOfLostTrackerHits[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_numberOfLostTrackerHits'));
      Muon3_innerTrack_numberOfLostTrackerInnerHits[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_numberOfLostTrackerInnerHits'));
      Muon3_numberOfMatchedStations[0] = float(getattr(old_tree_list[i], 'Muon3_numberOfMatchedStations'));
      Muon3_calEnergy_had[0] = float(getattr(old_tree_list[i], 'Muon3_calEnergy_had'));
      Muon3_innerTrack_numberofValidHits[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_numberofValidHits'));
      Muon3_calEnergy_em[0] = float(getattr(old_tree_list[i], 'Muon3_calEnergy_em'));
      Muon3_innerTrack_numberOfLostTrackerOuterHits[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_numberOfLostTrackerOuterHits'));
      Muon3_innerTrack_normalizedChi2[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_normalizedChi2'));
      Muon3_innerTrack_validFraction[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_validFraction'));
      Muon3_innerTrack_pixelLayersWithMeasurement[0] = float(getattr(old_tree_list[i], 'Muon3_innerTrack_pixelLayersWithMeasurement'));
      Muon3_ptErrOverPt[0] = float(getattr(old_tree_list[i], 'Muon3_ptErrOverPt'));
      muon_score = readerMuonId.EvaluateMVA("BDT");
      #print muon_score
    
      MuonId[0] = muon_score;
    
      new_tree_list[i].Fill();

   new_file.cd();
   new_tree_list[i].Write();

new_file.Close();
file_.Close();
