import os
from array import array

import ROOT

from collections import OrderedDict
from utils.yaml_utils import * 
from utils.varlist import *


def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

class DatacardHist():
   
   def __init__(self):

      self.input_file = ROOT.TFile()

      self.bdt_varlist = OrderedDict()
      self.bdt_speclist = OrderedDict()
      self.branches = {}

      self.trees = {
         'TreeB': { 'tree': None, 'channel': 'Bkg', 'SF': 0.0 },
         'TreeS_Ds': { 'tree': None, 'channel': 'DsToTau', 'SF': 0.0},
         'TreeS_Bu': { 'tree': None, 'channel': 'BuToTau', 'SF': 0.0},
         'TreeS_Bd': { 'tree':None, 'channel': 'BdToTau', 'SF': 0.0}
      }

      self.bdt_weightfiles = OrderedDict()
      self.classifiers = OrderedDict()
      self.bdt_cuts = {}
      self.histogram_dict = {}
      
      self.config = None
   
   def __enter__(self):
       return self

   def read_config_file(self, config_filename):

      self.config = ordered_load(open(config_filename), Loader=yaml.SafeLoader)

      if not ('INPUTFILE' in self.config):
         print "Specify the input file!"
         return

      else:
         if not os.path.exists(self.config['INPUTFILE']): 
            print "Input file does not exist!"
            return
         
         self.input_file = ROOT.TFile(self.config['INPUTFILE'], 'READ')
         for tree_ in self.trees: self.trees[tree_]['tree'] = self.input_file.Get(tree_)

      if not ('WEIGHTS' in self.config):
         print "Specify the paths to XML files containing BDT weights!"
         return
      for reader_ in self.config['WEIGHTS']:
         path_to_weights = self.config['WEIGHTS'][reader_]['file']
         if not os.path.exists(path_to_weights): 
            print path_to_weights+" does not exist!"
            return
         self.bdt_weightfiles[reader_] = path_to_weights
      
      
      if not ('BDT_CUTS' in self.config):
         print "Provide the BDT cuts!"
         return
      for cut_ in self.config['BDT_CUTS']:
         self.bdt_cuts[cut_] = self.config['BDT_CUTS'][cut_]

      if not ('VARSET' in self.config):
         print "Specify the set of variables used for training the BDT!"
         return
      else:
         for var_ in varsets[self.config['VARSET']]['training']:
            if var_ not in var_limits: self.bdt_varlist[var_] = var_ 
            else:   
               tmpvar_limits = var_limits[var_]
               min_expression = '('+var_+'>'+str(tmpvar_limits[0])+'?'+var_+':'+str(tmpvar_limits[0])+')'
               final_expression = min_expression+'<'+str(tmpvar_limits[1])+'?'+\
                       min_expression+':'+str(tmpvar_limits[1])
               self.bdt_varlist[var_] = final_expression

         if 'spectator' in varsets[self.config['VARSET']]:
            for var_ in varsets[self.config['VARSET']]['spectator']:
               if var_ not in var_limits: self.bdt_speclist[var_] = var_
               else:  
                  tmpvar_limits = var_limits[var_]
                  min_expression = '('+var_+'>'+str(tmpvar_limits[0])+'?'+var_+':'+str(tmpvar_limits[0])+')'
                  final_expression = min_expression+'<'+str(tmpvar_limits[1])+'?'+\
                          min_expression+':'+str(tmpvar[1])
                  self.bdt_speclist[var_] = final_expression
   
   def set_bdt_branch_address(self):

      for branchName in merge_two_dicts(self.bdt_varlist, self.bdt_speclist):
         if 'threeGlobal' in branchName: self.branches[branchName] = array('i', [0])
         else: self.branches[branchName] = array('f', [-999.])
         for tree_ in self.trees:
            self.trees[tree_]['tree'].SetBranchAddress(branchName, self.branches[branchName])
      
      for reader_ in self.bdt_weightfiles:
         cut_dict = self.config['WEIGHTS'][reader_]['cuts']
         for branchName in cut_dict:
            if branchName in self.branches: continue
            else: self.branches[branchName] = array('f', [-999.])
            for tree_ in self.trees:
               self.trees[tree_]['tree'].SetBranchAddress(branchName, self.branches[branchName])

   def initialize_classifiers(self):

      for reader_ in self.bdt_weightfiles:
         self.classifiers[reader_] = ROOT.TMVA.Reader('Color:!Silent')
         for var_ in self.bdt_varlist:
             self.classifiers[reader_].AddVariable(self.bdt_varlist[var_], self.branches[var_])
         for var_ in self.bdt_speclist:
             self.classifiers[reader_].AddSpectator(self.bdt_speclist[var_], self.branches[var_])
         self.classifiers[reader_].BookMVA('BDT', self.bdt_weightfiles[reader_])
   
   def get_bdt_score(self, reader_name):
       return self.classifiers[reader_name].EvaluateMVA('BDT')
   
   def get_lumi_scale(self, year=2018):
       sf_file = open('scale_factors.yml')
       sf_data = yaml.load(sf_file, Loader=yaml.FullLoader) 
       for tree_ in self.trees:
          if (self.trees[tree_]['channel']=='Bkg'):
              self.trees[tree_]['SF'] = 1.0
              continue
          self.trees[tree_]['SF'] = sf_data['SignalSF'][self.trees[tree_]['channel']][year] *\
                                    sf_data['SignalNormSF'][year]
       sf_file.close()
   
   def check_cuts(self, tree, cut_dict):
       passed = True
       for cut in cut_dict:
          var = self.branches[cut][0] 
          if ( (cut_dict[cut][0] == cut_dict[cut][1] ) and\
                  var==cut_dict[cut][0]): passed &= True
          elif ( var>=cut_dict[cut][0] and \
                  var<cut_dict[cut][1] ): passed &= True
          else: passed &= False
       return passed

   def fill_histogram(self):
       for tree_ in self.trees:
           nEntries = self.trees[tree_]['tree'].GetEntriesFast()
           print "Reading ", nEntries,"from ",tree_
           for i in xrange(nEntries):
               #if (i>10): break 
               if (i%1000==0): print "processing ",i,"/",nEntries,"..."
               self.trees[tree_]['tree'].GetEntry(i)
               for ireader, reader_ in enumerate(self.classifiers): 
                  
                  if not self.check_cuts(self.trees[tree_]['tree'], \
                          self.config['WEIGHTS'][reader_]['cuts']): continue
   
                  bdt_ = self.get_bdt_score(reader_)
                  ncuts = len(self.bdt_cuts[reader_])
                  for icut in xrange(ncuts):
                     cut_low = self.bdt_cuts[reader_][icut]
                     if (icut==0): cut_high = 999.0 
                     else: cut_high = self.bdt_cuts[reader_][icut-1]
   
                     if not ( bdt_ >  cut_low and bdt_ < cut_high ): continue
                     data_str = 'signal_'
                     if self.trees[tree_]['channel']=='Bkg': data_str = 'bkg_'
                     #print 'pass'
                     self.histogram_dict[reader_][data_str+chr(ireader+65)+str(icut+1)]\
                             .Fill(getattr(self.trees[tree_]['tree'], 'var_tauMass'),\
                                           self.trees[tree_]['SF'])

   def make_input_histograms(self, outputfile):
      
      self.get_lumi_scale()
      self.set_bdt_branch_address()
      self.initialize_classifiers()

      file_ = ROOT.TFile(outputfile, 'RECREATE')

      c = ROOT.TCanvas("c","Tau Mass",990,660)

      # Initialize the histograms
      for ireader, reader_ in enumerate(self.classifiers):
          self.histogram_dict[reader_] = {}
          ncuts = len(self.bdt_cuts[reader_])
          for icut in xrange(ncuts):
             
             hname_signal = 'h_signal_'+reader_+str(icut+1)
             hname_bkg = 'h_bkg_'+reader_+str(icut+1) 
             self.histogram_dict[reader_]['signal_'+chr(ireader+65)+str(icut+1)] = ROOT.TH1D(hname_signal,"",38,1.62,2.0)
             self.histogram_dict[reader_]['bkg_'+chr(ireader+65)+str(icut+1)] = ROOT.TH1D(hname_bkg,"",38,1.62,2.0)
             
      self.fill_histogram()
      
      for ireader, reader_ in enumerate(self.classifiers):     
         ncuts = len(self.bdt_cuts[reader_])
         for icut in xrange(ncuts):
            self.histogram_dict[reader_]['bkg_'+chr(ireader+65)+str(icut+1)].Write("data_obs"+chr(ireader+65)+str(icut+1))
            self.histogram_dict[reader_]['bkg_'+chr(ireader+65)+str(icut+1)].Write("background"+chr(ireader+65)+str(icut+1))
            self.histogram_dict[reader_]['signal_'+chr(ireader+65)+str(icut+1)].Write("signal"+chr(ireader+65)+str(icut+1))

      file_.Close()
   
   def makeMiniTree(self, outputfile):

       self.get_lumi_scale()
       self.set_bdt_branch_address()
       self.initialize_classifiers()

       file_ = ROOT.TFile(outputfile, 'RECREATE')

       for ireader, reader_ in enumerate(self.classifiers):
           ncuts = len(self.bdt_cuts[reader_])
           # to be filled later

   def __exit__(self, exc_type, exc_value, traceback):
       self.input_file.Close()

def main():

    #for file_ in ['card_2016vars.yml', 'card_BDTSegmentComp.yml', 'card_Muon3TimeAtIp.yml', 'card_Muon3TimeAtIp_MinBDTSegmComp.yml', 'card_TrackerMuonId.yml']:
    #for file_ in ['card_TrackerMuonId_globalMuonId.yml', 'card_BDTSegmentComp_globalMuonId.yml']:
    for file_ in ['cards_TrackerMuonId.yml']
       with DatacardHist() as _:
          _.read_config_file('./cards/'+file_)
          filename = 'inputs_for_combine/input_histograms_2glbTrk_'+file_.replace('card_','').split('.')[0]+'.root'
          print filename
          _.make_input_histograms(filename)


if __name__=='__main__':
    main()
