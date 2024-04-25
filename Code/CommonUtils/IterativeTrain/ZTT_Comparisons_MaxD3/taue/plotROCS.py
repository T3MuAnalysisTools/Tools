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
 
    w = 1400 
    h =  700
    can  = TCanvas("can", "", w, h)
    can.SetTitle("")
    can.SetFrameLineWidth(3);
    can.SetTickx();
    can.SetTicky();

    legend = ROOT.TLegend(0.12,0.50,0.43,0.66)
    legend.SetHeader("Trains")

    color = 1
    for ifile in args.files:
        print ifile
        datasets_prefix = ifile[:-5]
        print datasets_prefix
        color +=2
        rdict = get_file_dict(ifile)
#        hist = rdict["output_"+datasets_prefix]["Method_"+args.method][args.method]["MVA_"+args.method+"_trainingRejBvsS"]
        hist = rdict[datasets_prefix]["Method_"+args.method][args.method]["MVA_"+args.method+"_trainingRejBvsS"]
#output_4_A

        hist.SetTitle("")
        hist.SetLineColor(color)
        hist.SetStats(0)
        hist.GetXaxis().SetTitle("Efficiency");
        hist.GetYaxis().SetTitle("Rejection");
        hist.GetYaxis().SetRangeUser(0.8,1.05);
        hist.SetLineWidth(2)
        hist.Draw("same,C")
        legend.AddEntry(hist,str(ifile),"L")
        can.Update()

    cmd = 'mkdir plots'
    os.system(cmd)
    legend.Draw()


    can.SaveAs("plots/Compare_"+args.output+".png")
    can.SaveAs("plots/Compare_"+args.output+".root")


