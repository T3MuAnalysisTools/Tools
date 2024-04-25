#!/usr/bin/env python


import os
import sys
from ROOT import *
import argparse
import ROOT



parser = argparse.ArgumentParser()
parser.add_argument("-f","--files", help="List of files (coma separated); [Default: %(default)s]",  nargs="+", action="store", default = ["a.root","b.root"])
parser.add_argument("-m","--method", help="name of method; [Default: %(default)s]",  type = str, action="store", default = "BDT")
parser.add_argument("-o","--output", help="output files prefix; [Default: %(default)s]",  type = str, action="store", default = "Out")
args = parser.parse_args()
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# use your args
#print("connecting to {}".format(args.host))



def get_file_dict(rfilename):
    my_dict = {}
    rfile = ROOT.TFile(rfilename, "OPEN")
    for k in rfile.GetListOfKeys():
        cat = k.ReadObj()
        if isinstance(cat, ROOT.TDirectoryFile):
            my_dict[k.GetName()] = {}
            for kk in cat.GetListOfKeys():
                samp = kk.ReadObj()
                if isinstance(samp, ROOT.TDirectoryFile):
                    my_dict[k.GetName()][kk.GetName()] = {}
                    for kkk in samp.GetListOfKeys():
                        subsamp = kkk.ReadObj()
                        if isinstance(subsamp, ROOT.TDirectoryFile):
                            my_dict[k.GetName()][kk.GetName()][kkk.GetName()] = {}
                            for kkkk in subsamp.GetListOfKeys():
                                hist=kkkk.ReadObj()
                                if isinstance(hist, ROOT.TH1D):
                                    my_dict[k.GetName()][kk.GetName()][kkk.GetName()][kkkk.GetName()] = hist
                                    hist.SetDirectory(0)
    return my_dict





if __name__== "__main__":
 
    cat=0 #tauha: 0, #tauhb: 1, #taumu: 2, #taue: 3
    no_of_vars_tested = [15,18,19,16] # no of trainings tested - 1. tauha: 15, #tauhb: 18, #taumu: 19, #taue: 16
    title_type = ["{h,A}","{h,B}","{#mu}","{e}"]
    
    w = 1100 
    h =  700
    can  = TCanvas("can", "", w, h)
    can.SetTitle("")
    can.SetFrameLineWidth(3);
    can.SetTickx();
    can.SetTicky();

    legend = ROOT.TLegend(0.15,0.20,0.45,0.50)
    #legend = ROOT.TLegend(0.15,0.20,0.45,0.50)
    legend.SetHeader("Training")

    color = 1
    for ifile in args.files:
        print ifile
        datasets_prefix = ifile.split("/")[-1][9:-5]
        print datasets_prefix
        
        
        
        cat_no = ifile.split("/")[-1][9:-5].split("_")[0]
        no_of_vars = no_of_vars_tested[cat]+5-int(cat_no); 
        label = str(no_of_vars)+" variables, MaxDepth "+ifile.split("/")[0][-1]
        print label
        
        color +=2
        rdict = get_file_dict(ifile)
#        hist = rdict["output_"+datasets_prefix]["Method_"+args.method][args.method]["MVA_"+args.method+"_trainingRejBvsS"]
        hist = rdict["output_"+datasets_prefix]["Method_"+args.method][args.method]["MVA_"+args.method+"_trainingRejBvsS"]
#output_4_A

        hist.SetTitle("")
        hist.SetLineColor(color)
        hist.SetStats(0)
        hist.SetTitle("ROC Curve, #tau_"+title_type[cat]);
        hist.GetXaxis().SetTitle("Signal Efficiency");
        hist.GetYaxis().SetTitle("Background Rejection");
        #hist.GetYaxis().SetRangeUser(0.8,1.05);
        #hist.GetXaxis().SetRangeUser(0.9,1.05);
        hist.GetYaxis().SetRangeUser(0.80,1.05);
        hist.GetXaxis().SetRangeUser(0.90,1.0);
        
        hist.SetLineWidth(2)
        hist.Draw("same")
        
        legend.AddEntry(hist,str(label),"LL")
        can.Update()

    #Only because some part of histogram is drawn over the axes
    #axis_y = ROOT.TLine(0.90, 0.80, 0.90, 1.05)
    axis_y = ROOT.TLine(0.90, 0.80, 0.90, 1.05)
    axis_y.SetLineWidth(3)
    axis_y.Draw()
    
    cmd = 'mkdir plots'
    os.system(cmd)
    legend.Draw()


    can.SaveAs("plots/Compare_"+args.output+".png")
    can.SaveAs("plots/Compare_"+args.output+".root")


