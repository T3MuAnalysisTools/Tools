#!/usr/bin/env python

import argparse
import ROOT
import sys
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-data", "--data",help="Input root file with data tree", type = str, action="store", default = 'FillMVATree_2017_merged.root')
    parser.add_argument("-ds",   "--ds",help="Input root file with  MC (DsTau) tree", type = str, action="store", default = 'FillMVATree_DsToTau_merged.root')
    parser.add_argument("-bp",   "--bp",help="Input root file with  MC (BpTau) tree", type = str, action="store", default = 'FillMVATree_BpToTau_merged.root')
    parser.add_argument("-b0",   "--b0",help="Input root file with  MC (B0Tau) tree", type = str, action="store", default = 'FillMVATree_B0ToTau_merged.root')
    parser.add_argument("-o",    "--output-file", help="Output root file", type = str, action="store", default='FillMVATreeInput_combined.root')
    args = parser.parse_args()
    cwd = os.getcwd()

    input_file_data = ROOT.TFile(args.data, "READ")
    tree_data = input_file_data.Get("tree")

    input_file_ds = ROOT.TFile(args.ds, "READ")
    tree_ds = input_file_ds.Get("tree")
    
    input_file_bp = ROOT.TFile(args.bp, "READ")
    tree_bp = input_file_bp.Get("tree")

    input_file_b0 = ROOT.TFile(args.b0, "READ")
    tree_b0 = input_file_b0.Get("tree")   
    
    output_file = ROOT.TFile(args.output_file, "RECREATE")

    new_tree_data=tree_data.CloneTree()
    new_tree_data.SetName("TreeB")

    new_tree_ds=tree_ds.CloneTree()
    new_tree_ds.SetName("TreeS_Ds")
    output_file.Write()

    new_tree_bp=tree_bp.CloneTree()
    new_tree_bp.SetName("TreeS_Bu")
    output_file.Write()
    
    new_tree_b0=tree_b0.CloneTree()
    new_tree_b0.SetName("TreeS_Bd")
    output_file.Write()
