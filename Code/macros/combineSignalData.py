#!/usr/bin/env python

import argparse
import ROOT
import sys
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", "--first-file",help="Input root file with data tree", type = str, action="store", default = '')
    parser.add_argument("-f2", "--second-file",help="Input root file with  MC (DsTau) tree", type = str, action="store", default = '')
    parser.add_argument("-f3", "--third-file",help="Input root file with  MC (BuTau) tree", type = str, action="store", default = '')
    parser.add_argument("-f4", "--fourth-file",help="Input root file with  MC (BdTau) tree", type = str, action="store", default = '')
    parser.add_argument("-o",  "--output-file", help="Output root file", type = str, action="store", default='TMVATress.root')
    args = parser.parse_args()
    cwd = os.getcwd()

    input_file_1 = ROOT.TFile(args.first_file, "READ")
    tree1 = input_file_1.Get("tree")

    input_file_2 = ROOT.TFile(args.second_file, "READ")
    tree2 = input_file_2.Get("tree")
    
    input_file_3 = ROOT.TFile(args.third_file, "READ")
    tree3 = input_file_3.Get("tree")

    input_file_4 = ROOT.TFile(args.fourth_file, "READ")
    tree4 = input_file_4.Get("tree")   
    
    output_file = ROOT.TFile(args.output_file, "RECREATE")

    new_tree1=tree1.CloneTree()
    new_tree1.SetName("TreeB")

    new_tree2=tree2.CloneTree()
    new_tree2.SetName("TreeS_Ds")
    output_file.Write()

    new_tree3=tree3.CloneTree()
    new_tree3.SetName("TreeS_Bu")
    output_file.Write()
    
    new_tree4=tree4.CloneTree()
    new_tree4.SetName("TreeS_Bd")
    output_file.Write()
