#!/usr/bin/env python

import argparse
import ROOT
import sys
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", "--first-file",help="Input root file with data tree", type = str, action="store", default = '')
    parser.add_argument("-f2", "--second-file",help="Input root file with  MC tree", type = str, action="store", default = '')
    args = parser.parse_args()
    cwd = os.getcwd()

    input_file_1 = ROOT.TFile(args.first_file, "READ")
    tree1 = input_file_1.Get("tree")


    input_file_2 = ROOT.TFile(args.second_file, "READ")
    tree2 = input_file_2.Get("tree")


    output_file = ROOT.TFile("TMVATrees.root", "RECREATE")

    new_tree1=tree1.CloneTree()
    new_tree1.SetName("TreeB")

    new_tree2=tree2.CloneTree()
    new_tree2.SetName("TreeS") 
    output_file.Write()


