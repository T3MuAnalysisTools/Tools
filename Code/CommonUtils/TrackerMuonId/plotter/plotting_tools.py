import ROOT
from tdrstyle import setTDRStyle

#data visualization
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from root_pandas import read_root
import math
import rootplot.root2matplotlib as r2m

from utils.event_weights import *

setTDRStyle()
c = ROOT.TCanvas("c","",990,660)
#correlation heatmap of dataset
def correlation_heatmap(df, path_to_ouput_dir, filename):
    _ , ax = plt.subplots(figsize =(18, 16))
    colormap = sns.diverging_palette(220, 10, as_cmap = True)

    _ = sns.heatmap(
       df.corr(), 
       cmap = colormap,
       square=True, 
       cbar_kws={'shrink':.9 }, 
       ax=ax,
       annot=True, 
       linewidths=0.1,vmax=1.0, linecolor='white',
       annot_kws={'fontsize':10 }
    )
    plt.title('Pearson Correlation ('+filename+')', y=1.05, size=20)
    plt.savefig(path_to_output_dir)
    plt.close()

def fillHistograms(tree, hList):
    Nevents = tree.GetEntriesFast()
    for i in xrange(Nevents):
        tree.GetEntry(i)
        if (i%10000==0): print "Event %d/%d" % (i,Nevents)
        for h_ in hList: 
            if (tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==0): hList[h_]["hist_mu"].Fill(getattr(tree, hList[h_]["branch"]))
            if (tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==1): hList[h_]["hist_pi"].Fill(getattr(tree, hList[h_]["branch"]))    

def makePlot(c, h_mu, h_pi, variable, xlabel, ylabel, path_to_output_dir):
    leg = ROOT.TLegend(0.8, 0.8, 0.9, 0.9)
    h_mu.Draw("hist E")
    h_pi.Draw("same hist E")
    h_mu.GetXaxis().SetTitle(xlabel)
    h_mu.GetYaxis().SetTitle(ylabel)
    h_mu.SetLineColor(ROOT.kBlue)
    h_pi.SetLineColor(ROOT.kRed)
    h_mu.SetMarkerColor(ROOT.kBlue)
    h_pi.SetMarkerColor(ROOT.kRed)
    h_mu.SetLineWidth(2)
    h_mu.SetFillColorAlpha(ROOT.kBlue, 0.3)
    h_pi.SetFillColorAlpha(ROOT.kRed, 0.3)
    leg.AddEntry(h_pi,"#pi (D_{s}#rightarrow#phi(#mu#mu)#pi)","f")
    leg.AddEntry(h_mu,"#mu (D_{s}#rightarrow#tau#rightarrow3#mu)","f")
    leg.Draw("same")
    h_mu.GetYaxis().SetRangeUser(0,1.1*max(h_pi.GetMaximum(), h_mu.GetMaximum()))
    c.SaveAs(path_to_output_dir+variable+".pdf")

def getOverflowHist(h):
    nbins = h.GetNbinsX()
    binWidth = (h.GetBinLowEdge(h.GetNbinsX())-h.GetBinLowEdge(1))/(nbins-1)
    h_new = ROOT.TH1F(h.GetName()+"_overflow","",nbins+1, h.GetBinLowEdge(1), h.GetBinLowEdge(nbins)+2*binWidth)
    for i in range(1,nbins+2):
        h_new.Fill(h.GetBinCenter(i), h.GetBinContent(i))
        h_new.SetBinError(i, h.GetBinError(i))
    if not (h_new.Integral()==0): h_new.Scale(1.0/h_new.Integral())
    return h_new

def reweighHistograms(tree, hList, path_to_weightfile):
    mc_weights = event_weights()
    mc_weights.read_weights(path_to_weightfile)
    weight_ = 0
    Nevents = tree.GetEntriesFast()
    print "[plotter]: plotting reweighed histograms ..."
    for i in xrange(Nevents):
        tree.GetEntry(i)
        if (tree.fake==1 and tree.isPF==1 and tree.isTracker==1 and tree.isGlobal==0):
            pt_ = tree.muonPt
            eta_ = tree.muonEta
            weight_ = mc_weights.get_weight(pt_, eta_)
            if (i%1000==0): print pt_, eta_, weight_
            for h_ in hList: 
                hList[h_]["hist_pi"].Fill(getattr(tree, hList[h_]["branch"]), weight_)
