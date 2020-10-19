import ROOT

class event_weights():
   
   def __init__(self):
      self.weights = ROOT.TH2F('weights','',30,0,15,25,-2.5,2.5) # pt bins 50 [0,25], eta bins 50 [-2.5,2.5]

      # Initialize weights
      # rows for pt, cols for eta
      for row_ in xrange(1,31):
         for col_ in xrange(1,26):
            bin_ = self.weights.GetBin(row_,col_)
            self.weights.SetBinContent(bin_,0.0)

   def set_weights(self, h2_signal_pt_eta, h2_bkg_pt_eta):
      h2_signal_pt_eta.Scale(1/h2_signal_pt_eta.GetEntries())
      h2_bkg_pt_eta.Scale(1/h2_bkg_pt_eta.GetEntries())
      h2_signal_pt_eta.Divide(h2_bkg_pt_eta)

      for row_ in xrange(0,h2_signal_pt_eta.GetNbinsX()):
         for col_ in xrange(0,h2_signal_pt_eta.GetNbinsY()):
            bin_ = h2_signal_pt_eta.GetBin(row_+1, col_+1)
            self.weights.SetBinContent(bin_, h2_signal_pt_eta.GetBinContent(bin_))
            self.weights.SetBinError(bin_, h2_signal_pt_eta.GetBinError(bin_))

   def get_weight(self, pt, eta):
      xaxis = self.weights.GetXaxis()
      yaxis = self.weights.GetYaxis()
      binx = xaxis.FindBin(pt)
      biny = yaxis.FindBin(eta)
      bin_ = self.weights.GetBin(binx, biny)
      return self.weights.GetBinContent(bin_)
   
   def read_weights(self, filename):
      weight_file = ROOT.TFile(filename, "READ")
      h = weight_file.Get('weights')
      for row_ in xrange(0,h.GetNbinsX()):
         for col_ in xrange(0,h.GetNbinsY()):
            bin_ = h.GetBin(row_+1, col_+1)
            self.weights.SetBinContent(bin_, h.GetBinContent(bin_))
            self.weights.SetBinError(bin_, h.GetBinError(bin_))



class find_event_weights():
   
   def __init__(self, signal_cuts, bkg_cuts, filename, outputfilename):
      self.filename = filename
      self.outputfilename = outputfilename
      self.weights = event_weights()
      self.h2_signal_pt_eta = ROOT.TH2F("h2_signal_pt_eta", "", 30, 0, 15, 25, -2.5, 2.5)
      self.h2_bkg_pt_eta = ROOT.TH2F("h2_bkg_pt_eta", "", 30, 0, 15, 25, -2.5, 2.5)
      self.signal_cuts = signal_cuts
      self.bkg_cuts = bkg_cuts
      
   def fill_histogram(self):
      
      file_ = ROOT.TFile(self.filename, 'READ')

      weights = event_weights()
      tree = file_.Get('tree')
      nevents = tree.GetEntriesFast()
      
      print '[add_event_weights]: filling histograms...'
      for i in xrange(nevents):
         _ = tree.GetEntry(i)
         pt_ = -99.
         eta_ = -99.
         if (i%1000==0): print "processing %d/%d" % (i, nevents)
         #if (abs(tree.muonEta)<=2.4 and tree.isGlobal==0 and tree.isTracker==1 and tree.isPF==1) # replace by user defined string
         if (eval(self.signal_cuts)): 
            pt_ = tree.muonPt
            eta_ = tree.muonEta
            if (tree.muonPt>15): pt_ = 15
            self.h2_signal_pt_eta.Fill(pt_, eta_)
         elif (eval(self.bkg_cuts)): 
            pt_ = tree.muonPt
            eta_ = tree.muonEta
            if (tree.muonPt>15): pt_ = 15
            self.h2_bkg_pt_eta.Fill(pt_, eta_)
      
      output_file = ROOT.TFile(self.outputfilename,"RECREATE")
      self.h2_signal_pt_eta.Write()
      self.h2_bkg_pt_eta.Write()
      weights.set_weights(self.h2_signal_pt_eta, self.h2_bkg_pt_eta)

      weights.weights.Write()

      output_file.Close()
      file_.Close()
