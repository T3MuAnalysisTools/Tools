from utils.event_weights import *
_ = find_event_weights("abs(tree.muonEta)<=2.4 and tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==0", "abs(tree.muonEta)<=2.4 and tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==1", "../Ntuples/MuonPionTree_combined_2017.root","../Weights/event_weights_2017.root")
_.fill_histogram()
