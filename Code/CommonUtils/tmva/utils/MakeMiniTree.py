#!/usr/bin/env python
 
from ROOT import TMVA, TFile, TString, TTree
from array import array
#from subprocess import call
#from os.path import isfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--input-file",help="input file root file with MVA tree; [Default: %(default)s] ", action="store", default = 'FillMVATreeInput_combined.root')
parser.add_argument("-w", "--weights-prefix",help="The BDT weight prefix; [Default: %(default)s] ", action="store", default = 'MyTMVAClassification_Test_vars_')
parser.add_argument("-ls", "--lumi-scale",help="The BDT weight prefix; [Default: %(default)s] ", action="store", default = 'MyTMVAClassification_Test_vars_')
args = parser.parse_args()

BASEDIR = 'datasets/weights'
filename = args.input_file
print "Opening file: "+filename

def fileStr(str):
    category = str.split('_')[0].replace('reader','')
    return category



DsScale = 0.001241
B0Scale = 0.000388
BpScale = 0.000056



# list of variables
varlist_spec = ['var_tauMass']
#varlist_train = ['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_pmin','var_max_cLP','var_max_tKink','var_MinD0Significance','var_mindca_iso','var_trk_relPt']
#varlist_train = ['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_pmin','var_max_cLP','var_max_tKink','var_MinD0Significance','var_mindca_iso','var_trk_relPt']
#varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
#                 'var_MaxD0SigBS','var_maxMuonsDca','var_minMuonsDca',
#                 'var_segCompMuMin','var_Iso08MuMin',
#                 'var_Iso08MuMax','var_MaxVertexPairQuality','var_MinMuon_chi2LocalPosition','var_NtracksClose','var_nsv']

#varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_segCompMuMin', 'var_Muon3P',
#                 'var_MinMuon_chi2LocalPosition', 'var_MaxtrkKink', 'var_MinD0SigPV', 'var_minMuonsDca', 'var_relPt']

#varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
#                   'var_MaxD0SigBS','var_MaxD0SigSV','var_MindcaTrackSV','var_maxMuonsDca','var_nsv',
#                   'var_segCompMuMin','var_MaxMuon_chi2LocalPosition','var_NtracksClose','var_MaxVertexPairQuality','var_Iso08','var_Iso08MuMin',
#                   'var_MaxMuonImpactAngle']

#varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_segCompMuMin', 'var_Muon3P',
#                   'var_MinMuon_chi2LocalPosition', 'var_MaxtrkKink', 'var_MinD0SigPV', 'var_minMuonsDca', 'var_Muon1PFIsoTight','var_Muon2PFIsoTight','var_Muon3PFIsoTight']

#varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_segCompMuMin', 'var_Muon3P',
#                   'var_MinMuon_chi2LocalPosition', 'var_MaxtrkKink', 'var_MinD0SigPV', 'var_minMuonsDca', 'var_Muon1TightId','var_Muon2TightId','var_Muon3TightId']

varlist_train = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                   'var_MaxD0SigBS','var_MaxD0SigSV','var_MindcaTrackSV','var_maxMuonsDca','var_nsv',
                   'var_segCompMuMin','var_MaxMuon_chi2LocalPosition','var_NtracksClose','var_MaxVertexPairQuality','var_Iso08','var_Iso08MuMin',
                   'var_MaxMuonImpactAngle','var_MaxtrkKink','var_MuonglbkinkSum','var_MinMatchedStations','var_MaxdeltaMuZ']







varlist_help = ['var_mass12','var_mass13']




print "Initializing TMVA..."
# Setup TMVA
TMVA.Tools.Instance()
mvaReaderList = ['readerA','readerB','readerC','readerBC']
mvaReaders = {}
mvaWeights = {}

print "Setting up reader..."
for reader in mvaReaderList:
    mvaReaders[reader] = TMVA.Reader('Color:!Silent')
#    tmpstr = 'MyTMVAClassification_Test_vars_'+fileStr(reader)
    tmpstr = args.weights_prefix+'_vars_'+fileStr(reader)
#--weights-prefix

    mvaWeights[reader] = TString(BASEDIR+'/'+tmpstr+'_BDT.weights.xml') # weights weights.xml file after training, place it to CommonFiles

#Vars2016_vars_



# Open input file and set branch addresses
rootfile = TFile(filename, 'READ')
rootfile.ls()
tree_bkg = rootfile.Get('TreeB')
tree_ds  = rootfile.Get('TreeS_Ds')
tree_bu  = rootfile.Get('TreeS_Bu')
tree_bd  = rootfile.Get('TreeS_Bd')


# List of branches
branches = {}
branchList_miniTree = ['m3m','dataMCType','event_weight','bdt','category','LumiScale','m12','m13']
branches_miniTree = {}

for branchName in varlist_train:
    branches[branchName] = array('f', [-999])
    for reader in mvaReaders: mvaReaders[reader].AddVariable(branchName, branches[branchName])
    tree_ds.SetBranchAddress(branchName, branches[branchName])
    tree_bu.SetBranchAddress(branchName, branches[branchName])
    tree_bd.SetBranchAddress(branchName, branches[branchName])
    tree_bkg.SetBranchAddress(branchName, branches[branchName])


for branchName in varlist_spec:
    branches[branchName] = array('f', [-999])
    for reader in mvaReaders: mvaReaders[reader].AddSpectator(branchName, branches[branchName])
    tree_ds.SetBranchAddress(branchName, branches[branchName])
    tree_bu.SetBranchAddress(branchName, branches[branchName])
    tree_bd.SetBranchAddress(branchName, branches[branchName])
    tree_bkg.SetBranchAddress(branchName, branches[branchName])


for branchName in varlist_help:
    branches[branchName] = array('f', [-999])
    tree_ds.SetBranchAddress(branchName, branches[branchName])
    tree_bu.SetBranchAddress(branchName, branches[branchName])
    tree_bd.SetBranchAddress(branchName, branches[branchName])
    tree_bkg.SetBranchAddress(branchName, branches[branchName])



branches['category'] = array('i', [-999])
tree_ds.SetBranchAddress('category', branches['category'])
tree_bu.SetBranchAddress('category', branches['category'])
tree_bd.SetBranchAddress('category', branches['category'])
tree_bkg.SetBranchAddress('category', branches['category'])




#Setup a new tree
T3MFMiniTree = TFile("T3MMiniTree"+"_"+args.weights_prefix+".root","recreate");
T3MMiniTree = TTree('T3MMiniTree','T3MMiniTree');
T3MMiniTree.SetDirectory(T3MFMiniTree);

for branchName in branchList_miniTree:
    branches_miniTree[branchName] = array('f', [-999])
    T3MMiniTree.Branch(branchName, branches_miniTree[branchName],branchName+'/F')


# Book methods
for reader in mvaReaderList:
     mvaReaders[reader].BookMVA('BDT', mvaWeights[reader]) 


print branches_miniTree

# print some example classifications
print('Classifying signal events (Ds)')
for i in range(tree_ds.GetEntriesFast()):
    tree_ds.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 40
    if ['category'][0]==1: 
        reader = 'readerA'

    if ['category'][0]==2: 
        reader = 'readerB'

    if ['category'][0]==3: 
        reader = 'readerC'

    if ['category'][0]==2 or ['category'][0]==3: 
        reader = 'readerBC'

    score = mvaReaders[reader].EvaluateMVA('BDT')
    branches_miniTree['bdt'][0] = score
    branches_miniTree['category'][0] = branches['category'][0]
    branches_miniTree['m12'][0] = branches['var_mass12'][0]
    branches_miniTree['m13'][0] = branches['var_mass13'][0]
    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
    branches_miniTree['event_weight'][0] = DsScale
    branches_miniTree['LumiScale'][0] = 1.0
    T3MMiniTree.Fill()



print('Classifying signal events (Bu)')
for i in range(tree_bu.GetEntriesFast()):
    tree_bu.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 60
    if ['category'][0]==1: 
        reader = 'readerA'

    if ['category'][0]==2: 
        reader = 'readerB'

    if ['category'][0]==3: 
        reader = 'readerC'

    if ['category'][0]==2 or ['category'][0]==3: 
        reader = 'readerBC'
    score = mvaReaders[reader].EvaluateMVA('BDT')
    branches_miniTree['bdt'][0] = score
    branches_miniTree['category'][0] = branches['category'][0]
    branches_miniTree['m12'][0] = branches['var_mass12'][0]
    branches_miniTree['m13'][0] = branches['var_mass13'][0]
    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
    branches_miniTree['event_weight'][0] = B0Scale
    branches_miniTree['LumiScale'][0] = 1.0
    T3MMiniTree.Fill()



print('Classifying signal events (Bd)')
for i in range(tree_bd.GetEntriesFast()):
    tree_bd.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 90
    if ['category'][0]==1: 
        reader = 'readerA'

    if ['category'][0]==2: 
        reader = 'readerB'

    if ['category'][0]==3: 
        reader = 'readerC'

    if ['category'][0]==2 or ['category'][0]==3: 
        reader = 'readerBC'
    score = mvaReaders[reader].EvaluateMVA('BDT')
    branches_miniTree['bdt'][0] = score
    branches_miniTree['category'][0] = branches['category'][0]
    branches_miniTree['m12'][0] = branches['var_mass12'][0]
    branches_miniTree['m13'][0] = branches['var_mass13'][0]
    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
    branches_miniTree['event_weight'][0] = BpScale
    branches_miniTree['LumiScale'][0] = 1.0
    T3MMiniTree.Fill()


 
print('Classifying background events')
for i in range(tree_bkg.GetEntriesFast()):
    tree_bkg.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 1
    if ['category'][0]==1: 
        reader = 'readerA'

    if ['category'][0]==2: 
        reader = 'readerB'

    if ['category'][0]==3: 
        reader = 'readerC'

    if ['category'][0]==2 or ['category'][0]==3: 
        reader = 'readerBC'
    score = mvaReaders[reader].EvaluateMVA('BDT')
    branches_miniTree['bdt'][0] = score
    branches_miniTree['category'][0] = branches['category'][0]
    branches_miniTree['m12'][0] = branches['var_mass12'][0]
    branches_miniTree['m13'][0] = branches['var_mass13'][0]
    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
    branches_miniTree['event_weight'][0] = 1.0
    branches_miniTree['LumiScale'][0] = 1.0
    T3MMiniTree.Fill()



T3MFMiniTree.Write();
T3MFMiniTree.Close();
rootfile.Close();
