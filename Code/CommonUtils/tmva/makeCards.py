#! /usr/bin/env python
# script to create yaml file containing all the config parameters

import ROOT
from utils.MVAPlottingTools import *
from utils.varlist import *
import argparse
import yaml
from collections import OrderedDict
from utils.yaml_utils import * 

TAU_MASS_MIN = 1.60
TAU_MASS_MAX = 2.00

BASEDIR = './weights/' 

def get_dump(key, item):
    _ = {}
    _[key] = item
    return _

def makeCard(input_file, card_name, varset_):

   weights_to_dump = {}
   cuts_to_dump = {}
   varset_to_dump = {}

   tmva_filename = ''
   path_to_weight_file = ''
   
   SCALE_DS = 0.00132 # Normalization factor for Ds->Tau channel
   SCALE_BU = 0.000441 # Normalization factor for Bu->Tau channel
   SCALE_BD = 0.00122 # Normalization factor for Bd->Tau channel
   DATA_NORM = 0.93 # Normalization factor from DsToPhiPi channel
   
   rootfile = ROOT.TFile(input_file, 'READ')
   treeB = rootfile.Get("TreeB")
   treeDs = rootfile.Get("TreeS_Ds")
   treeBu = rootfile.Get("TreeS_Bu")
   treeBd = rootfile.Get("TreeS_Bd")
    
   for i, cat in enumerate(['A','B','C'], 1):
      for type_ in ['3glb', '2glbTrk']:
         category_cuts = {}
         cuts = 'category=='+str(i)
         weight_name = ''
         category_cuts['category'] = [i, i]
         if (type_=='3glb'):
            weight_name = 'TMVAClassification_'+varset_+'_'+cat+'_threeGlobal'
            tmva_filename = weight_name+'_Output.root'
            path_to_weight_file = BASEDIR+weight_name+'/weights/'+weight_name+'.weights.xml'
            cuts += " && threeGlobal==1"
            category_cuts['threeGlobal'] = [1,1]
         elif (type_=='2glbTrk'):
            weight_name = 'TMVAClassification_'+varset_+'_'+cat+'_twoGlobalTracker'
            tmva_filename = weight_name+'_Output.root' 
            cuts += " && threeGlobal==0"
            category_cuts['threeGlobal'] = [0,0]
            path_to_weight_file = BASEDIR+weight_name+'/weights/'+weight_name+'_BDT.weights.xml'
         
         nevents_data = treeB.Draw("category", cuts)
         nevents_ds = treeDs.Draw("category", cuts)
         nevents_bu = treeBu.Draw("category", cuts) 
         nevents_bd = treeBd.Draw("category", cuts)

         reader = 'reader'+cat+'_'+type_

         path_to_tmva_file = './tmva_output/'+tmva_filename
         print path_to_tmva_file
         norm = DATA_NORM*(SCALE_DS*nevents_ds+SCALE_BU*nevents_bu+SCALE_BD*nevents_bd)/(nevents_ds+nevents_bu+nevents_bd)
         if (os.path.exists(path_to_tmva_file)): cuts_to_dump[reader] = BDT_optimal_cut(path_to_tmva_file, weight_name, norm)
         if (os.path.exists(path_to_weight_file)):
             weights_to_dump[reader] = { 'file':  path_to_weight_file, 'cuts': category_cuts }
   
   rootfile.Close()

   cardfile = open(card_name, 'w+')
   cardfile.write('INPUTFILE: '+input_file+'\n\n')
   cardfile.write('VARSET: '+varset_+'\n\n')
   cardfile.write(yaml.dump(get_dump('WEIGHTS', weights_to_dump), \
           Dumper=indent_dumper, default_flow_style=False))
   cardfile.write('\n')
   cardfile.write(yaml.dump(get_dump('BDT_CUTS', cuts_to_dump), \
           Dumper=indent_dumper, default_flow_style=False))
   cardfile.close()

def main():

   #--------------------- --------------------------- Argument Parser ------------------------ ------------------------
   parser = argparse.ArgumentParser()
   parser.add_argument('-i','--input-file', type=str, help='Input file; Default:(%default)s', default='')
   parser.add_argument('-c','--card_name',type=str,help='Output yaml card; Default:(%default)s',default='example.yml')
   args = parser.parse_args()
   #--------------------- --------------------------- -------------- ------------------------- ------------------------

   from utils.varlist import varsets

   for varset in varsets:
      card_name = './cards/card_'+varset+'.yml'
      makeCard(args.input_file, card_name, varset)


if __name__=='__main__':
   main()
