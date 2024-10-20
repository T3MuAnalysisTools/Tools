import ROOT

hList_reweighed = {
"muonPt":{"branch":"muonPt", "xlabel":"muonPt", "ylabel":"Entries","hist_pi":ROOT.TH1F("h_muonPt_pi_reweighed","",30,0,15)},
"muonEta":{"branch":"muonEta", "xlabel":"muonEta", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonEta_pi_reweighed","",25,-2.5,2.5)},
"muonPhi":{"branch":"muonPhi", "xlabel":"muonPhi", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonPhi_pi_reweighed","",36,-3.6,3.6)},
"muonInnerNC2":{"branch":"muonInnerNC2", "xlabel":"muonInnerNC2", "ylabel":"Entries","hist_pi":ROOT.TH1F("h_muonInnerNC2_pi_reweighed","",25,0,5)},
"muonInnerNValidHits":{"branch":"muonInnerNValidHits", "xlabel":"muonInnerNValidHits", "ylabel":"Entries","hist_pi":ROOT.TH1F("h_muonInnerNValidHits_pi_reweighed","",30,0,30)},
"muonNValidTrackerHits":{"branch":"muonNValidTrackerHits", "xlabel":"muonNValidTrackerHits", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNValidTrackerHits_pi_reweighed","",30,0,30)},
"muonNLostTrackerHits":{"branch":"muonNLostTrackerHits", "xlabel":"muonNLostTrackerHits", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNLostTrackerHits_pi_reweighed","",7,0,7)},
"muonNLostTrackerHitsInner":{"branch":"muonNLostTrackerHitsInner", "xlabel":"muonNLostTrackerHitsInner", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNLostTrackerHitsInner_pi_reweighed","",7,0,7)},
"muonNLostTrackerHitsOuter":{"branch":"muonNLostTrackerHitsOuter", "xlabel":"muonNLostTrackerHitsOuter", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNLostTrackerHitsOuter_pi_reweighed","",7,0,7)},
"muonPixelLayers":{"branch":"muonPixelLayers", "xlabel":"muonPixelLayers", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonPixelLayers_pi_reweighed","",6,0,6)},
"muonNMatchedStations":{"branch":"muonNMatchedStations", "xlabel":"muonNMatchedStations", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNMatchedStations_pi_reweighed","",6,0,6)},
"muonNMatches":{"branch":"muonNMatches", "xlabel":"muonNMatches", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonNMatches_pi_reweighed","",10,0,10)},
"muonRPC":{"branch":"muonRPC", "xlabel":"muonRPC", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonRPC_pi_reweighed","",2,0,2)},
"muonPtErrPt":{"branch":"muonPtErrPt", "xlabel":"muonPtErrPt", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonPtErrPt_pi_reweighed","",20,0,0.4)},
"muonSegComp":{"branch":"muonSegComp", "xlabel":"muonSegComp", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonSegComp_pi_reweighed","",25,0,1)},
"muonCaloComp":{"branch":"muonCaloComp", "xlabel":"muonCaloComp", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonCaloComp_pi_reweighed","",25,0,1)},
"muonHadS9":{"branch":"muonHadS9", "xlabel":"muonHadS9", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonHadS9_pi_reweighed","",30,0,30)},
"muonHad":{"branch":"muonHad", "xlabel":"muonHad", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonHad_pi_reweighed","",30,0,30)},
"muonEM":{"branch":"muonEM", "xlabel":"muonEM", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonEM_pi_reweighed","",25,0,2)},
"muonEMS9":{"branch":"muonEMS9", "xlabel":"muonEMS9", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonEMS9_pi_reweighed","",25,0,2)},
"muonEMS25":{"branch":"muonEMS25", "xlabel":"muonEMS25", "ylabel":"Entries", "hist_pi":ROOT.TH1F("h_muonEMS25_pi_reweighed","",25,0,2)},
"muonKink":{"branch":"muonKink", "xlabel":"muonKink", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonKink_mu_reweighed","",100,0,100),"hist_pi":ROOT.TH1F("h_muonKink_pi_reweighed","",100,0,100)}
}

hList = {
"muonPt":{"branch":"muonPt", "xlabel":"muonPt", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonPt_mu","",30,0,15),"hist_pi":ROOT.TH1F("h_muonPt_pi","",30,0,15)},
"muonEta":{"branch":"muonEta", "xlabel":"muonEta", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonEta_mu","",25,-2.5,2.5),"hist_pi":ROOT.TH1F("h_muonEta_pi","",25,-2.5,2.5)},
"muonPhi":{"branch":"muonPhi", "xlabel":"muonPhi", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonPhi_mu","",36,-3.6,3.6),"hist_pi":ROOT.TH1F("h_muonPhi_pi","",36,-3.6,3.6)},
"muonInnerNC2":{"branch":"muonInnerNC2", "xlabel":"muonInnerNC2", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonInnerNC2_mu","",25,0,5),"hist_pi":ROOT.TH1F("h_muonInnerNC2_pi","",25,0,5)},
"muonInnerNValidHits":{"branch":"muonInnerNValidHits", "xlabel":"muonInnerNValidHits", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonInnerNValidHits_mu","",30,0,30),"hist_pi":ROOT.TH1F("h_muonInnerNValidHits_pi","",30,0,30)},
"muonNValidTrackerHits":{"branch":"muonNValidTrackerHits", "xlabel":"muonNValidTrackerHits", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNValidTrackerHits_mu","",30,0,30),"hist_pi":ROOT.TH1F("h_muonNValidTrackerHits_pi","",30,0,30)},
"muonNLostTrackerHits":{"branch":"muonNLostTrackerHits", "xlabel":"muonNLostTrackerHits", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNLostTrackerHits_mu","",7,0,7),"hist_pi":ROOT.TH1F("h_muonNLostTrackerHits_pi","",7,0,7)},
"muonNLostTrackerHitsInner":{"branch":"muonNLostTrackerHitsInner", "xlabel":"muonNLostTrackerHitsInner", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNLostTrackerHitsInner_mu","",7,0,7),"hist_pi":ROOT.TH1F("h_muonNLostTrackerHitsInner_pi","",7,0,7)},
"muonNLostTrackerHitsOuter":{"branch":"muonNLostTrackerHitsOuter", "xlabel":"muonNLostTrackerHitsOuter", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNLostTrackerHitsOuter_mu","",7,0,7),"hist_pi":ROOT.TH1F("h_muonNLostTrackerHitsOuter_pi","",7,0,7)},
"muonPixelLayers":{"branch":"muonPixelLayers", "xlabel":"muonPixelLayers", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonPixelLayers_mu","",6,0,6),"hist_pi":ROOT.TH1F("h_muonPixelLayers_pi","",6,0,6)},
"muonTrackerLayers":{"branch":"muonTrackerLayers", "xlabel":"muonTrackerLayers", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonTrackerLayers_mu","",10,0,10),"hist_pi":ROOT.TH1F("h_muonTrackerLayers_pi","",10,0,10)},
"muonNMatchedStations":{"branch":"muonNMatchedStations", "xlabel":"muonNMatchedStations", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNMatchedStations_mu","",6,0,6),"hist_pi":ROOT.TH1F("h_muonNMatchedStations_pi","",6,0,6)},
"muonNMatches":{"branch":"muonNMatches", "xlabel":"muonNMatches", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonNMatches_mu","",10,0,10),"hist_pi":ROOT.TH1F("h_muonNMatches_pi","",10,0,10)},
"muonRPC":{"branch":"muonRPC", "xlabel":"muonRPC", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonRPC_mu","",2,0,2),"hist_pi":ROOT.TH1F("h_muonRPC_pi","",2,0,2)},
"muonPtErrPt":{"branch":"muonPtErrPt", "xlabel":"muonPtErrPt", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonPtErrPt_mu","",20,0,0.4),"hist_pi":ROOT.TH1F("h_muonPtErrPt_pi","",20,0,0.4)},
"muonSegComp":{"branch":"muonSegComp", "xlabel":"muonSegComp", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonSegComp_mu","",25,0,1),"hist_pi":ROOT.TH1F("h_muonSegComp_pi","",25,0,1)},
"muonCaloComp":{"branch":"muonCaloComp", "xlabel":"muonCaloComp", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonCaloComp_mu","",25,0,1),"hist_pi":ROOT.TH1F("h_muonCaloComp_pi","",25,0,1)},
"muonHadS9":{"branch":"muonHadS9", "xlabel":"muonHadS9", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonHadS9_mu","",30,0,30),"hist_pi":ROOT.TH1F("h_muonHadS9_pi","",30,0,30)},
"muonHad":{"branch":"muonHad", "xlabel":"muonHad", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonHad_mu","",30,0,30),"hist_pi":ROOT.TH1F("h_muonHad_pi","",30,0,30)},
"muonEM":{"branch":"muonEM", "xlabel":"muonEM", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonEM_mu","",25,0,2),"hist_pi":ROOT.TH1F("h_muonEM_pi","",25,0,2)},
"muonEMS9":{"branch":"muonEMS9", "xlabel":"muonEMS9", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonEMS9_mu","",25,0,2),"hist_pi":ROOT.TH1F("h_muonEMS9_pi","",25,0,2)},
"muonEMS25":{"branch":"muonEMS25", "xlabel":"muonEMS25", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonEMS25_mu","",25,0,2),"hist_pi":ROOT.TH1F("h_muonEMS25_pi","",25,0,2)},
"muonKink":{"branch":"muonKink", "xlabel":"muonKink", "ylabel":"Entries", "hist_mu":ROOT.TH1F("h_muonKink_mu","",100,0,100),"hist_pi":ROOT.TH1F("h_muonKink_pi","",100,0,100)}
}
