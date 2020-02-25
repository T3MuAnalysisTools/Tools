#!/usr/bin/env python
 
from ROOT import TMVA, TFile, TString, TTree
from array import array
#from subprocess import call
#from os.path import isfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="Input file used for MVA training", default="FillMVATreeInput_combined.root", type=str)
args = parser.parse_args()

BASEDIR = '../CommonFiles/weights/'
filename = args.filename
print "Opening file: "+filename
# define and initialize cuts
bdt_cuts = {}
for cat in ['A','B','C']:
   for Type in ['3glb', '2glbTrk']:
        for rank in ['2','1']: bdt_cuts['bdt_cut'+cat+rank+'_'+Type] = -99.0;

# import cuts from a card
_ = open('bdt_cuts.dat')
cutValues = _.readlines()
_.close()

print("Reading bdt cuts...")
for line in cutValues: bdt_cuts[line.split(' ')[0]] = float(line.split(' ')[1].strip('\n'))

print bdt_cuts
def fileStr(str):
    category = str.split('_')[0].replace('reader','')
    Type = str.split('_')[1]
    if (Type=='3glb'): return category+'_threeGlobal'
    elif (Type=='2glbTrk'): return category+'_twoGlobaltrk'
    return ''

def getRank(reader, score):
    Type = reader.split('_')[1]
    cut_min = 'bdt_cut'+reader.split('_')[0].replace('reader','')+'2_'+Type
    cut_max = 'bdt_cut'+reader.split('_')[0].replace('reader','')+'1_'+Type
    if ( score > bdt_cuts[cut_min] ) and ( score <= bdt_cuts[cut_max] ): return 2
    elif ( score > bdt_cuts[cut_max] ): return 1
    else: return 3

# list of variables
varlist_spec = ['var_tauMass', 'var_tauMassRes', 'var_Eta_Tau','threeGlobal']
varlist_train = ['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_pmin','var_max_cLP','var_max_tKink','var_MinD0Significance','var_mindca_iso','var_trk_relPt']

# List of branches
branches = {}
branchList_miniTree = ['m3m','dataMCType','event_weight','bdt','category','eta','LumiScale']
branches_miniTree = {}

print "Initializing TMVA..."
# Setup TMVA
TMVA.Tools.Instance()
#TMVA.PyMethodBase.PyInitialize()
mvaReaderList = ['readerA_3glb','readerB_3glb','readerC_3glb','readerA_2glbTrk','readerB_2glbTrk','readerC_2glbTrk']
mvaReaders = {}
mvaWeights = {}

print "Setting up reader..."
for reader in mvaReaderList:
    mvaReaders[reader] = TMVA.Reader('Color:!Silent')
    tmpstr = 'TMVAClassification_2016vars_'+fileStr(reader)
    mvaWeights[reader] = TString(BASEDIR+tmpstr+'/weights/'+tmpstr+'_BDT.weights.xml') # weights weights.xml file after training, place it to CommonFiles


# Open input file and set branch addresses
rootfile = TFile(filename, 'READ')
rootfile.ls()
tree_bkg = rootfile.Get('TreeB')
tree_ds = rootfile.Get('TreeS_Ds')
tree_bu = rootfile.Get('TreeS_Bu')
tree_bd = rootfile.Get('TreeS_Bd')

for branchName in varlist_train:
    branches[branchName] = array('f', [-999])
    for reader in mvaReaders: mvaReaders[reader].AddVariable(branchName, branches[branchName])
    tree_ds.SetBranchAddress(branchName, branches[branchName])
    tree_bu.SetBranchAddress(branchName, branches[branchName])
    tree_bd.SetBranchAddress(branchName, branches[branchName])
    tree_bkg.SetBranchAddress(branchName, branches[branchName])

for branchName in varlist_spec:
    if 'threeGlobal' in branchName: branches[branchName] = array('i', [0])
    else: branches[branchName] = array('f', [-999])
    for reader in mvaReaders: mvaReaders[reader].AddSpectator(branchName, branches[branchName])
    tree_ds.SetBranchAddress(branchName, branches[branchName])
    tree_bu.SetBranchAddress(branchName, branches[branchName])
    tree_bd.SetBranchAddress(branchName, branches[branchName])
    tree_bkg.SetBranchAddress(branchName, branches[branchName])

#Setup a new tree
T3MFMiniTree = TFile("T3MMiniTree.root","recreate");
T3MMiniTree = TTree('T3MMiniTree','T3MMiniTree');
T3MMiniTree.SetDirectory(T3MFMiniTree);
for branchName in branchList_miniTree:
    branches_miniTree[branchName] = array('f', [-999])
    T3MMiniTree.Branch(branchName, branches_miniTree[branchName],branchName+'/F')

# Book methods
for reader in mvaReaderList:
     mvaReaders[reader].BookMVA('BDT', mvaWeights[reader]) 

# Function return (1, reader) of event passes mass cuts (0, reader) otherwise
def invMassCut(branches, isData):
    category = ''
    Type = ''
    passed = False
    if (branches['var_tauMassRes'][0] <= 0.007): category = 'A'
    elif (branches['var_tauMassRes'][0] > 0.007 and branches['var_tauMassRes'][0] <= 0.01): category = 'B'
    elif (branches['var_tauMassRes'][0] > 0.01): category = 'C'
    if (branches['threeGlobal'][0]==1): Type = '3glb'
    else: Type = '2glbTrk'
    reader = 'reader'+category+'_'+Type
    tauMass = branches['var_tauMass'][0]
    if ( isData and ( (tauMass > 1.65 and tauMass < 1.75) or (tauMass > 1.80 and tauMass<1.90) ) ): passed = True
    elif (( not isData ) and ( tauMass < 1.80 and tauMass > 1.75 ) ): passed = True
    return (passed, reader)

print branches_miniTree

# print some example classifications
print('Classifying signal events (Ds)')
print("Total number of events in Ds: ", tree_ds.GetEntriesFast())
for i in range(tree_ds.GetEntriesFast()):
    tree_ds.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 40
    (passed, reader) = invMassCut(branches, 0)
    if (passed):
        score = mvaReaders[reader].EvaluateMVA('BDT')
        branches_miniTree['bdt'][0] = score
        rank = getRank(reader, score)
        if (rank < 3):
            if (i%1000==0): print('BDT Score: ',score)
            branches_miniTree['category'][0] = 2*mvaReaderList.index(reader)+rank
	    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
	    branches_miniTree['event_weight'][0] = 0.11970
            branches_miniTree['eta'][0] = branches['var_Eta_Tau'][0]
            branches_miniTree['LumiScale'][0] = 0.0225
            T3MMiniTree.Fill()

print('Classifying signal events (Bu)')
print("Total number of events in Bu: ", tree_bu.GetEntriesFast())
for i in range(tree_bu.GetEntriesFast()):
    tree_bu.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 60
    (passed, reader) = invMassCut(branches, 0)
    if (passed):
        score = mvaReaders[reader].EvaluateMVA('BDT')
        branches_miniTree['bdt'][0] = score
        rank = getRank(reader, score)
        if (rank < 3):
            if (i%1000==0): print('BDT Score: ',score)
            type(branches_miniTree['category'][0])
            branches_miniTree['category'][0] = 2*mvaReaderList.index(reader)+rank
	    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
	    branches_miniTree['event_weight'][0] = 0.03067
	    branches_miniTree['eta'][0] = branches['var_Eta_Tau'][0]
	    branches_miniTree['LumiScale'][0] = 0.01372
	    T3MMiniTree.Fill()

print('Classifying signal events (Bd)')
print("Total number of events in Bd: ", tree_bd.GetEntriesFast())
for i in range(tree_bd.GetEntriesFast()):
    tree_bd.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 90
    (passed, reader) = invMassCut(branches, 0)
    if (passed):
        score = mvaReaders[reader].EvaluateMVA('BDT')
        branches_miniTree['bdt'][0] = score
        rank = getRank(reader, score)
        if (rank < 3):
            if (i%1000==0): print('BDT Score: ',score)
            branches_miniTree['category'][0] = 2*mvaReaderList.index(reader)+rank
	    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
	    branches_miniTree['event_weight'][0] = 0.02963
	    branches_miniTree['eta'][0] = branches['var_Eta_Tau'][0]
	    branches_miniTree['LumiScale'][0] = 0.020555
	    T3MMiniTree.Fill()
print('Classifying background events')
print("Total number of background events: ",tree_bkg.GetEntriesFast())
for i in range(tree_bkg.GetEntriesFast()):
    tree_bkg.GetEntry(i)
    branches_miniTree['dataMCType'][0] = 1
    (passed, reader) = invMassCut(branches, 1)
    if (passed):
        score = mvaReaders[reader].EvaluateMVA('BDT')
        branches_miniTree['bdt'][0] = score
        rank = getRank(reader, score)
        if (rank < 3):
            if (i%1000==0): print('Reader, BDT score, rank: ', reader, score, rank)
            branches_miniTree['category'][0] = 2*mvaReaderList.index(reader)+rank
	    branches_miniTree['m3m'][0] = branches['var_tauMass'][0]
	    branches_miniTree['event_weight'][0] = 1.0
	    branches_miniTree['eta'][0] = branches['var_Eta_Tau'][0]
	    branches_miniTree['LumiScale'][0] = 1.0
	    T3MMiniTree.Fill()

T3MFMiniTree.Write();
T3MFMiniTree.Close();
rootfile.Close();
