import ROOT
import math
import numpy as np
import os
from tdrstyle import setTDRStyle
from varlist import var_limits, varsets

setTDRStyle()

# define histogram related parameters

Color = [ROOT.kBlack, ROOT.kRed, ROOT.kOrange+6, ROOT.kGreen+3, ROOT.kBlue, ROOT.kViolet-1]

# significance optimization

def TH1_integral (h, xmin, xmax):
    axis = h.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    integral = h.Integral(bmin,bmax)
    integral -= h.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin))/axis.GetBinWidth(bmin)
    integral -= h.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax)/ axis.GetBinWidth(bmax)
    if (integral<0): return 0
    return integral


def BDT_optimal_cut(filename, tag, norm_factor, path_to_output='./Plots/bdt_cuts/'):
    a = 0.0
    b = 0.0
    N_s_1=0
    N_b_1=0
    N_s_2=0
    N_b_2=0
    S1=0
    S2=0
    S=0
    S1_list = [] ; S2_list = [] ;S_list = [] ;a_list = [] ;b_list = [];
    sig_norm = norm_factor # normalization
    #Double_t sig_norm = 1. #average normalization factor for the three signal samples
    f = ROOT.TFile(filename,"READ")
    if (f.IsZombie()):
        print "Unable to open file: "+filename
        return

    h_test_signal = ROOT.TH1F()
    h_test_bkg = ROOT.TH1F()
    h_test_signal = f.Get("datasets/Method_BDT/BDT/MVA_BDT_S_high")
    h_test_bkg = f.Get("datasets/Method_BDT/BDT/MVA_BDT_B_high")

    #Signal is normalized to "sig_norm" factor
    h_test_signal.Scale(1.0/h_test_signal.Integral())
    h_test_bkg.Scale(1.0/h_test_bkg.Integral())
    #Make up on plots
    h_test_signal.GetXaxis().SetRangeUser(-0.5,0.5)
    h_test_bkg.SetLineColor(ROOT.kBlue)
    h_test_signal.SetLineColor(ROOT.kRed)
    h_test_bkg.SetFillColorAlpha(ROOT.kBlue,0.35)
    h_test_signal.SetFillColorAlpha(ROOT.kRed,0.35)
    h_test_bkg.SetMarkerColor(ROOT.kBlue)
    h_test_signal.SetMarkerColor(ROOT.kRed)
    h_test_bkg.SetMarkerStyle(21)
    h_test_signal.SetMarkerStyle(21)
    h_test_bkg.SetFillStyle(3004)
    h_test_signal.SetFillStyle(3004)

    #Loop on both cuts in [-0.50.5]
    X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_bkg.GetXaxis().GetXmin())
    X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_bkg.GetXaxis().GetXmax())
    N = 50
    step = (X_max - X_min)/N
    dim = 0

    for i in xrange(N-1):
       b = X_min + i*step
       for j in range(i+1,N):
           a = X_min + j * step
           N_s_1 = TH1_integral(h_test_signal,a,X_max)*sig_norm
           N_b_1 = TH1_integral(h_test_bkg,a,X_max)
           N_s_2 = TH1_integral(h_test_signal,b,a)*sig_norm
           N_b_2 = TH1_integral(h_test_bkg,b,a)
           if(N_s_1 < TH1_integral(h_test_signal,X_min,X_max)*sig_norm*0.0005): continue
           if(N_b_1 < TH1_integral(h_test_bkg,X_min,X_max)*0.0005): continue
           if(N_s_2 < TH1_integral(h_test_signal,X_min, X_max)*sig_norm*0.0001): continue
           if(N_b_2 < TH1_integral(h_test_signal,X_min, X_max)*0.0001): continue
           if (N_s_1<0 or N_s_2<0 or N_b_1<0 or N_b_2<0): print "values:",N_s_1, N_s_2, N_b_1, N_b_2
           if ( (N_b_1)>0 and (N_b_2)>0 ):
              #S1 = N_s_1 / np.sqrt(N_s_1 + N_b_1)
              #S2 = N_s_2 / np.sqrt(N_s_2 + N_b_2)
              S1 = np.sqrt( 2*( (N_s_1+N_b_1)*np.log(1+ (N_s_1/N_b_1)) - N_s_1 ) )
              S2 = np.sqrt( 2*( (N_s_2+N_b_2)*np.log(1+ (N_s_2/N_b_2)) - N_s_2 ) )
              #Combined significance
              S = np.sqrt(S1*S1 + S2*S2)
              a_list.append(a)
              b_list.append(b)
              S_list.append(S)
              dim += 1

    c3 = ROOT.TCanvas("c3","c3",150,10,990,660)
    g2d = ROOT.TGraph2D(dim, np.array(a_list), np.array(b_list), np.array(S_list))
    #g2d.SetTitle("SignalSignificance")
    g2d.GetHistogram().GetXaxis().SetTitle("a")
    g2d.GetHistogram().GetYaxis().SetTitle("b")
    g2d.GetHistogram().GetXaxis().SetLimits(min(a_list),max(a_list))
    g2d.GetHistogram().GetYaxis().SetLimits(min(b_list),max(b_list))
    g2d.Draw("colz1")

    #Taking absolute maximum of the combined significance
    S_max = max(S_list)
    S_maxIndex = S_list.index(S_max)

    a_max = a_list[S_maxIndex]
    b_max = b_list[S_maxIndex]

    #Computing cut efficiency on signal
    N_S_12 = TH1_integral(h_test_signal,b_max,0.5)
    N_S_tot = TH1_integral(h_test_signal,-0.5,0.5)

    l = ROOT.TLine()
    l.DrawLine(a_max,min(b_list),a_max,max(b_list))
    l.DrawLine(min(a_list),b_max,max(a_list),b_max)

    c3.Update()
    c3.SaveAs(path_to_output+tag+"_2dplot.png")

    #Drawing BDT score from scratch without signal normalization
    c2 = ROOT.TCanvas("c2","c2",150,10,990,660)
    c2.cd()
    h_test_signal2 = ROOT.TH1F()
    h_test_bkg2 = ROOT.TH1F()
    h_test_signal2 = f.Get("datasets/Method_BDT/BDT/MVA_BDT_S")
    h_test_bkg2 = f.Get("datasets/Method_BDT/BDT/MVA_BDT_B")
    h_test_signal2.GetXaxis().SetRangeUser(-0.5,0.5)
    h_test_bkg2.SetLineColor(ROOT.kBlue)
    h_test_signal2.SetLineColor(ROOT.kRed)
    h_test_bkg2.SetMarkerColor(ROOT.kBlue)
    h_test_signal2.SetMarkerColor(ROOT.kRed)
    h_test_bkg2.SetMarkerStyle(21)
    h_test_signal2.SetMarkerStyle(21)

    h_test_bkg2.Draw()
    h_test_signal2.Draw("same")
    l.DrawLine(a_max,0,a_max,max(h_test_signal2.GetMaximum(),h_test_bkg2.GetMaximum()))
    l.DrawLine(b_max,0,b_max,max(h_test_signal2.GetMaximum(),h_test_bkg2.GetMaximum()))
    h_test_bkg2.GetYaxis().SetRangeUser(0,1.1*max(h_test_signal2.GetMaximum(),h_test_bkg2.GetMaximum()))
    h_test_signal2.GetYaxis().SetRangeUser(0,1.1*max(h_test_signal2.GetMaximum(),h_test_bkg2.GetMaximum()))

    leg2 = ROOT.TLegend(0.2,0.8,0.4,0.9)
    leg2.AddEntry(h_test_signal2,"signal","p")
    leg2.AddEntry(h_test_bkg2,"bkg","p")
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    #leg2.SetNColumns(2)
    leg2.Draw("same")
    c2.Update()
    c2.Modified()
    c2.SaveAs(path_to_output+"optimal_cuts/"+tag+".png")
    f.Close()
    return [a_list[S_maxIndex], b_list[S_maxIndex]] 


def CompareROC(inputFiles, titles, plotName):
    rootFiles = []
    for iFile in inputFiles: rootFiles.append(ROOT.TFile(iFile, "READ"))

    rocCurves = []
    for rFile in rootFiles: rocCurves.append(rFile.Get("datasets/Method_BDT/BDT/MVA_BDT_rejBvsS"))

    c1 = ROOT.TCanvas("c1","c1",150,10,990,660)
    leg = ROOT.TLegend(0.15,0.2,0.45,0.4)

    for i, roc in enumerate(rocCurves,0):
       if (i==0): roc.Draw()
       else: roc.Draw("same")
       roc.GetXaxis().SetTitle("Bkg rejection")
       roc.GetYaxis().SetTitle("Signal efficiency")
       roc.SetLineColor(Color[i])
       roc.SetMarkerColor(Color[i])
       roc.SetMarkerStyle(21)
       leg.AddEntry(roc,titles[i],"p")

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("same")

    c1.SaveAs(plotName)
    for rFile in rootFiles: rFile.Close()


def plotVariables(filename, var, plotName):

    c1 = ROOT.TCanvas("c1","c1",150,10,990,660)
    c1.cd()
    rootFile = ROOT.TFile(filename,"READ")

    # trees
    _treeB = rootFile.Get("TreeB")
    _treeDs = rootFile.Get("TreeS_Ds")
    _treeBu = rootFile.Get("TreeS_Bu")
    _treeBd = rootFile.Get("TreeS_Bd")

    # histograms
    h1_bkg = ROOT.TH1F("h1_bkg","Background",100,var_limits[var][0],var_limits[var][1])
    h1_ds = ROOT.TH1F("h1_ds","Signal (Ds)",100,var_limits[var][0],var_limits[var][1])
    h1_bu = ROOT.TH1F("h1_bu","Signal (Bu)",100,var_limits[var][0],var_limits[var][1])
    h1_bd = ROOT.TH1F("h1_bd","Signal (Bd)",100,var_limits[var][0],var_limits[var][1])

    #get histograms
    _treeB.Draw(var+">>h1_bkg","","hist")
    _treeDs.Draw(var+">>h1_ds","","hist")
    _treeBu.Draw(var+">>h1_bu","","hist")
    _treeBd.Draw(var+">>h1_bd","","hist")

    # Add colors
    h1_bkg.SetMarkerColor(255)
    h1_bkg.SetLineColor(255)
    h1_bkg.SetMarkerStyle(21)
    h1_ds.SetMarkerColor(ROOT.kBlue)
    h1_ds.SetLineColor(ROOT.kBlue)
    h1_ds.SetFillColorAlpha(ROOT.kBlue,0.35)
    h1_ds.SetMarkerStyle(21)
    h1_ds.SetFillStyle(3004)
    h1_bu.SetMarkerColor(ROOT.kGreen+3)
    h1_bu.SetLineColor(ROOT.kGreen+3)
    h1_bu.SetFillColorAlpha(ROOT.kGreen+3,0.35)
    h1_bu.SetMarkerStyle(21)
    h1_bu.SetFillStyle(3004)
    h1_bd.SetMarkerColor(ROOT.kOrange+7)
    h1_bd.SetLineColor(ROOT.kOrange+7)
    h1_bd.SetFillColorAlpha(ROOT.kOrange+7,0.35)
    h1_bd.SetMarkerStyle(21)
    h1_bd.SetFillStyle(3004)

    # legend
    leg = ROOT.TLegend(0.6,0.7,0.8,0.9)
    leg.AddEntry(h1_bkg,'Background','p')
    leg.AddEntry(h1_ds,'D_{s}#rightarrow#tau#rightarrow3#mu','p')
    leg.AddEntry(h1_bu,'B_{u}#rightarrow#tau#rightarrow3#mu','p')
    leg.AddEntry(h1_bd,'B_{d}#rightarrow#tau#rightarrow3#mu','p')

    #normalize
    h1_bkg.Scale(1/h1_bkg.Integral())
    h1_ds.Scale(1/h1_ds.Integral())
    h1_bu.Scale(1/h1_bu.Integral())
    h1_bd.Scale(1/h1_bd.Integral())

    #Draw
    h1_bkg.Draw("")
    h1_ds.Draw("same hist")
    h1_bd.Draw("same hist")
    h1_bu.Draw("same hist")

    leg.Draw("same")
    c1.SaveAs(plotName)

    rootFile.Close()


def compareVariables(inputFiles, treeName, titles, variable, plotName):
    rootFiles = []
    for iFile in inputFiles: 
       rootFiles.append(ROOT.TFile(iFile, "READ"))

    plotList = []
    for i, rFile in enumerate(rootFiles):
       plotList.append(ROOT.TH1F("h_"+titles[i],variable,100,var_limits[variable][0],var_limits[variable][1]))
       _tree = rFile.Get(treeName)
       _tree.Draw(variable+">>h_"+titles[i].replace('.root',''),"","hist")

    c1 = ROOT.TCanvas("c1","c1",150,10,990,660)
    leg = ROOT.TLegend(0.75,0.7,0.85,0.9)

    for i, plot in enumerate(plotList,0):
       if (i==0): plot.Draw("hist")
       else: plot.Draw("same hist")
       plot.GetXaxis().SetTitle(variable)
       plot.GetYaxis().SetTitle("a.u.")
       plot.SetLineColor(Color[i])
       plot.SetMarkerColor(Color[i])
       plot.SetMarkerStyle(21)
       plot.Scale(1.0/(plot.GetEntries()))
       plot.SetFillStyle(3004)
       leg.AddEntry(plot,titles[i],"p")

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("same")

    c1.SaveAs(plotName)
    for rFile in rootFiles: rFile.Close()


def plotTMVAInputVars(filename, varlist, output):

   c1 = ROOT.TCanvas("c1","",1980,1980)
   c2 = ROOT.TCanvas("c2","",1980,1980)

   c1.Divide(3,2)

   rootfile = ROOT.TFile(filename, "READ")
   tree_bkg = rootfile.Get("TreeB")
   tree_ds = rootfile.Get("TreeS_Ds")
   tree_bu = rootfile.Get("TreeS_Bu")
   tree_bd = rootfile.Get("TreeS_Bd")

   h1_bkg = []
   h1_ds = []
   h1_bu = []
   h1_bd = []

   for num, var in enumerate(varlist):
      
      if(num<6): c1.cd((num/3)+(num%3))
      else: c2.cd((num-6)/3+(num-6)%3)

      h1_bkg.append(ROOT.TH1F("h1_bkg_"+var,"Background "+var,100,var_limits[var][0],var_limits[var][1]))
      h1_ds.append(ROOT.TH1F("h1_ds_"+var,"Signal (Ds) "+var,100,var_limits[var][0],var_limits[var][1]))
      h1_bu.append(ROOT.TH1F("h1_bu_"+var,"Signal (Bu) "+var,100,var_limits[var][0],var_limits[var][1]))
      h1_bd.append(ROOT.TH1F("h1_bd_"+var,"Signal (Bd) "+var,100,var_limits[var][0],var_limits[var][1]))

      tree_bkg.Draw(var+">>h1_bkg_"+var,"","hist")
      tree_ds.Draw(var+">>h1_ds_"+var,"","hist")
      tree_bu.Draw(var+">>h1_bu_"+var,"","hist")
      tree_bd.Draw(var+">>h1_bd_"+var,"","hist")

      # Add colors
      h1_bkg[num].SetTitle(var)
      h1_bkg[num].SetMarkerColor(255)
      h1_bkg[num].SetLineColor(255)
      h1_bkg[num].SetLineWidth(2)
      h1_bkg[num].SetMarkerStyle(21)
      h1_ds[num].SetMarkerColor(ROOT.kBlue)
      h1_ds[num].SetLineColor(ROOT.kBlue)
      h1_ds[num].SetLineWidth(2)
      h1_ds[num].SetFillColorAlpha(ROOT.kBlue,0.35)
      h1_ds[num].SetMarkerStyle(21)
      h1_ds[num].SetFillStyle(3004)
      h1_bu[num].SetMarkerColor(ROOT.kGreen+3)
      h1_bu[num].SetLineColor(ROOT.kGreen+3)
      h1_bu[num].SetLineWidth(2)
      h1_bu[num].SetFillColorAlpha(ROOT.kGreen+3,0.35)
      h1_bu[num].SetMarkerStyle(21)
      h1_bu[num].SetFillStyle(3004)
      h1_bd[num].SetMarkerColor(ROOT.kOrange+7)
      h1_bd[num].SetLineColor(ROOT.kOrange+7)
      h1_bd[num].SetLineWidth(2)
      h1_bd[num].SetFillColorAlpha(ROOT.kOrange+7,0.35)
      h1_bd[num].SetMarkerStyle(21)
      h1_bd[num].SetFillStyle(3004)
      
      #normalize
      h1_bkg[num].Scale(1/h1_bkg[num].Integral())
      h1_ds[num].Scale(1/h1_ds[num].Integral())
      h1_bu[num].Scale(1/h1_bu[num].Integral())
      h1_bd[num].Scale(1/h1_bd[num].Integral())

      #Draw
      DrawOverflowBin(h1_bkg[num]).Draw("same")
      DrawOverflowBin(h1_ds[num]).Draw("same hist")
      DrawOverflowBin(h1_bu[num]).Draw("same hist")
      DrawOverflowBin(h1_bd[num]).Draw("same hist")

   c1.SaveAs(output+"_c1.png")
   c2.SaveAs(output+"_c2.png")


def DrawOverflowBin(h):

    #function to paint the histogram h with an extra bin for overflows
    nx = h.GetNbinsX()+1
    xbins = np.zeros(nx+1)
    for i in xrange(nx):
       xbins[i]=h.GetBinLowEdge(i+1)
    
    xbins[nx]=xbins[nx-1]+h.GetBinWidth(nx)
    
    # book a temporary histogram having extra bins for overflows
    htmp = ROOT.TH1F(h.GetName(), h.GetTitle(), nx, xbins)
    htmp.Sumw2()

    # fill the new histogram including the overflows
    for i in range(1,nx+1):
       htmp.SetBinContent(htmp.FindBin(htmp.GetBinCenter(i)),h.GetBinContent(i))
       htmp.SetBinError(htmp.FindBin(htmp.GetBinCenter(i)),h.GetBinError(i))

       htmp.SetBinContent(htmp.FindBin(h.GetBinLowEdge(1)-1), h.GetBinContent(0))
       htmp.SetBinError(htmp.FindBin(h.GetBinLowEdge(1)-1), h.GetBinError(0))

    # Restore the number of entries
    htmp.SetEntries(h.GetEffectiveEntries())
    return htmp
