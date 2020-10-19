#! /usr/bin/env python
from utils.event_weights import *

# ---------- Argument parser -------------
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input-file", help="Name of the input file for TMVA classification; [Default: %(default)s]", default="../Ntuples/MuonPionTree_combined_2017.root")
parser.add_argument("-o","--output-file", help="Name of the output file; [Default: %(default)s]", default="../Weights/event_weights_2017.root")
args = parser.parse_args()
# ----------------------------------------

cuts = "abs(tree.muonEta)<=2.4 and tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==0"
cutb = "abs(tree.muonEta)<=2.4 and tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1 and tree.fake==1"

_ = find_event_weights(cuts, cutb, args.input_file, args.output_file)
_.fill_histogram()
