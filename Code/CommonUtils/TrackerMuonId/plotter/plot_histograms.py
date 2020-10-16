from plotting_tools import *
from histogram_list import *
from utils.event_weights import *

import ROOT

rootfile = ROOT.TFile("../Ntuples/MuonPionTree_combined_2017.root","READ")
tree = rootfile.Get("tree")

fillHistograms(tree, hList)
#for h_ in hList: makePlot(c, getOverflowHist(hList[h_]["hist_mu"]), getOverflowHist(hList[h_]["hist_pi"]),h_,hList[h_]["xlabel"],hList[h_]["ylabel"], "../Plots/2017/Unweighed/")

reweighHistograms(tree, hList_reweighed, "../Weights/event_weights_2017.root")
for h_ in hList_reweighed: makePlot(c, getOverflowHist(hList[h_]["hist_mu"]), getOverflowHist(hList_reweighed[h_]["hist_pi"]),h_,hList[h_]["xlabel"],hList[h_]["ylabel"], "../Plots/2017/Reweighed/")
