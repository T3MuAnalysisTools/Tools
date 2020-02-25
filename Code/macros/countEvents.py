import sys
import ROOT
from tdrstyle import setTDRStyle
from CMSStyle import CMS_lumi

#------------------
# define histograms
#------------------

setTDRStyle()

h1_3glb_A1_signal=ROOT.TH1D("h1_3glb_A1_signal"," 3 global muons (A1)",25,1.65,1.90)
h1_3glb_A2_signal=ROOT.TH1D("h1_3glb_A2_signal"," 3 global muons (A2)",25,1.65,1.90)
h1_3glb_B1_signal=ROOT.TH1D("h1_3glb_B1_signal"," 3 global muons (B1)",25,1.65,1.90)
h1_3glb_B2_signal=ROOT.TH1D("h1_3glb_B2_signal"," 3 global muons (B2)",25,1.65,1.90)
h1_3glb_C1_signal=ROOT.TH1D("h1_3glb_C1_signal"," 3 global muons (C1)",25,1.65,1.90)
h1_3glb_C2_signal=ROOT.TH1D("h1_3glb_C2_signal"," 3 global muons (C2)",25,1.65,1.90)
h1_2glbTrk_A1_signal=ROOT.TH1D("h1_2glbTrk_A1_signal"," 2 global muons and tracker muon (A1)",25,1.65,1.90)
h1_2glbTrk_A2_signal=ROOT.TH1D("h1_2glbTrk_A2_signal"," 2 global muons and tracker muon (A2)",25,1.65,1.90)
h1_2glbTrk_B1_signal=ROOT.TH1D("h1_2glbTrk_B1_signal"," 2 global muons and tracker muon (B1)",25,1.65,1.90)
h1_2glbTrk_B2_signal=ROOT.TH1D("h1_2glbTrk_B2_signal"," 2 global muons and tracker muon (B2)",25,1.65,1.90)
h1_2glbTrk_C1_signal=ROOT.TH1D("h1_2glbTrk_C1_signal"," 2 global muons and tracker muon (C1)",25,1.65,1.90)
h1_2glbTrk_C2_signal=ROOT.TH1D("h1_2glbTrk_C2_signal"," 2 global muons and tracker muon (C2)",25,1.65,1.90)

h1_3glb_A1_bkg=ROOT.TH1D("h1_3glb_A1_bkg"," 3 global muons (A1)",25,1.65,1.90)
h1_3glb_A2_bkg=ROOT.TH1D("h1_3glb_A2_bkg"," 3 global muons (A2)",25,1.65,1.90)
h1_3glb_B1_bkg=ROOT.TH1D("h1_3glb_B1_bkg"," 3 global muons (B1)",25,1.65,1.90)
h1_3glb_B2_bkg=ROOT.TH1D("h1_3glb_B2_bkg"," 3 global muons (B2)",25,1.65,1.90)
h1_3glb_C1_bkg=ROOT.TH1D("h1_3glb_C1_bkg"," 3 global muons (C1)",25,1.65,1.90)
h1_3glb_C2_bkg=ROOT.TH1D("h1_3glb_C2_bkg"," 3 global muons (C2)",25,1.65,1.90)
h1_2glbTrk_A1_bkg=ROOT.TH1D("h1_2glbTrk_A1_bkg"," 2 global muons and tracker muon (A1)",25,1.65,1.90)
h1_2glbTrk_A2_bkg=ROOT.TH1D("h1_2glbTrk_A2_bkg"," 2 global muons and tracker muon (A2)",25,1.65,1.90)
h1_2glbTrk_B1_bkg=ROOT.TH1D("h1_2glbTrk_B1_bkg"," 2 global muons and tracker muon (B1)",25,1.65,1.90)
h1_2glbTrk_B2_bkg=ROOT.TH1D("h1_2glbTrk_B2_bkg"," 2 global muons and tracker muon (B2)",25,1.65,1.90)
h1_2glbTrk_C1_bkg=ROOT.TH1D("h1_2glbTrk_C1_bkg"," 2 global muons and tracker muon (C1)",25,1.65,1.90)
h1_2glbTrk_C2_bkg=ROOT.TH1D("h1_2glbTrk_C2_bkg"," 2 global muons and tracker muon (C2)",25,1.65,1.90)

bkg_category = {1: h1_3glb_A1_bkg, 2: h1_3glb_A2_bkg, 3: h1_3glb_B1_bkg, 4: h1_3glb_B2_bkg, 5: h1_3glb_C1_bkg, 6: h1_3glb_C2_bkg, 7: h1_2glbTrk_A1_bkg, 8: h1_2glbTrk_A2_bkg, 9: h1_2glbTrk_B1_bkg, 10: h1_2glbTrk_B2_bkg, 11: h1_2glbTrk_C1_bkg, 12: h1_2glbTrk_C2_bkg}
signal_category = {1: h1_3glb_A1_signal, 2: h1_3glb_A2_signal, 3: h1_3glb_B1_signal, 4: h1_3glb_B2_signal, 5: h1_3glb_C1_signal, 6: h1_3glb_C2_signal, 7: h1_2glbTrk_A1_signal, 8: h1_2glbTrk_A2_signal, 9: h1_2glbTrk_B1_signal, 10: h1_2glbTrk_B2_signal, 11: h1_2glbTrk_C1_signal, 12: h1_2glbTrk_C2_signal}


c = ROOT.TCanvas("c","Tau Mass",990,660)

file = ROOT.TFile("T3MMiniTree.root","READ")
tree = file.Get("T3MMiniTree")

for i in xrange(tree.GetEntriesFast()):
    tree.GetEntry(i)
    if (tree.dataMCType==1):
        bkg_category[tree.category].Fill(tree.m3m)
    else:
        signal_category[tree.category].Fill(tree.m3m,(tree.LumiScale*tree.event_weight))

for j in bkg_category:
    c.cd()
    signal_category[j].SetLineColor(ROOT.kRed)
    signal_category[j].SetLineWidth(2)
    bkg_category[j].Draw("pe")
    signal_category[j].Draw("hist same")
    CMS_lumi(c, 4, 0, lumi_13TeV = 'Run 2017 ( 41.22 fb^{-1})', relPosX = 0.130)
    c.SaveAs(bkg_category[j].GetName()+".png")
    c.Clear()

sys.stdout.write("-------------------------------\n")
for j in bkg_category:
    sys.stdout.write(str(signal_category[j].GetName()).replace("_signal",""))
    sys.stdout.write("\t")
    sys.stdout.write(str(signal_category[j].Integral()))
    sys.stdout.write("\t")
    sys.stdout.write(str(bkg_category[j].Integral()))
    sys.stdout.write("\n")
sys.stdout.write("-------------------------------\n")
