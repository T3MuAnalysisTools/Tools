#------------------------------------------------------------------
#Script to combine the plots from all the Sets and fit the tau mass
#------------------------------------------------------------------

import ROOT
import os, sys
import csv
import argparse
from t3mUtils import *
from ROOT import RooFit
from ROOT import TFile, TH1F, TLine, TCanvas, TPad
from ROOT import RooRealVar, RooCBShape, RooGaussian, RooAddPdf, RooDataHist, RooArgSet, RooArgList, RooAbsArg
from tdrstyle import setTDRStyle

setTDRStyle() # use cms tdrstyle

#----- Argument parser --------
parser = argparse.ArgumentParser()

parser.add_argument("--parameterSet",
                    help=PRINTER.prepend("file containing list of parameters\n \
                    Format for crystal ball: CONSTANT\tCONSTANT_MIN\tCONSTANT_MAX\tMEAN\tMEAN_MIN\t \
                    MEAN_MAX\tSIGMA\tSIGMA_MIN\tSIGMA_MAX\tALPHA\tALPHA_MIN\tALPHA_MAX\tN\tN_MIN\tN_MAX\n \ 
                    Format for gaussian: CONSTANT\tMEAN\tSIGMA"))
parser.add_argument("--fitType", help=PRINTER.prepend("type of fit: (1: gaussian, 2: crystal ball, 3: double gaussian, 4: double crystal ball 5: gaussian + crystalball)"))
parser.add_argument("--inputFile", help=PRINTER.prepend("specify the root file"))

args = parser.parse_args()
#------------------------------

PRINTER = t3m_printer('[Tau23MuFitter]: ') # printer
#------------------------------

FIT_TYPES = ['1','2','3','4']

if not (args.fitType in FIT_TYPES):
    PRINTER.print("Invalid fit type!\n")

if not (os.path.exists(args.inputFile)):
    PRINTER.print(args.inputFile+" root file doesn't exist!")
    parser.print_help()
    sys.exit()
else: rootFile = TFile(args.inputFile)

if not (os.path.exists(args.parameterSet)):
    PRINTER.print(args.parameterSet+" does not exist!")
    parser.print_help()
    sys.exit()
else:
        if ()
#----------------------
# Setup canvas and pads
#----------------------

c = TCanvas("c","Tau Mass Fit",1600,1200)
pad1 = TPad("pad1", "Pad for histogram", 0,0.3,1,1.0)
pad2 = TPad("pad2", "Pad for the ratio plot", 0,0.0,1,0.3)

#---------------
# Get histograms
#---------------

rootFile = TFile("LOCAL_ANALYSIS_threemu_MC_merged.root")

h1 = rootFile.Get("threemu_default_TauMassMC2")
h2 = rootFile.Get("threemu_default_TauMassMC3")
h3 = rootFile.Get("threemu_default_TauMassMC4")

h0 = TH1F("h0","Tau Mass",40,1.5,1.9)
h0.Add(h1)
h0.Add(h2)
h0.Add(h3)
#-------------
# Setup model
#-------------

# Declare variables x,mean,sigma with associated name, title, initial value and allowed range
x = RooRealVar("x","x",1.5,1.9)

# Parameters of crystal ball fit
crystal_mean = RooRealVar("crystal_mean","mean of crystal ball",1.77,1.5,1.9)
crystal_sigma = RooRealVar("crystal_sigma","width of crystal ball",0.02,0.0,0.1)
n = RooRealVar("n","power parameter in the crystal ball",1,0,20)
alpha = RooRealVar("alpha","boundry in the crystal ball",1.75,1.6,1.9)
crystal_constant = RooRealVar("crystal_constant", "constant", 7417.0, 1769.0, 8011.0)

#Parameters of gaussian fit
gaus_mean = RooRealVar("gaus_mean","mean of the gaussian",1.77,1.5,1.9)
gaus_sigma = RooRealVar("gaus_sigma","width of the gaussian",0.01,0.0,0.1)
gaus_constant = RooRealVar("gau_constant","constant of the gaussian",100,0,1000)
# Build crystall p.d.f in terms of x,mean and sigma
crystalball = RooCBShape("CBShape", "Cystal Ball Function", x, crystal_mean, crystal_sigma, alpha, n)
gaus = RooGaussian("GaussShape", "Gaussian Function", x, gaus_mean, gaus_sigma)

model = RooAddPdf("model", "model", RooArgList(crystalball,gaus), RooArgList(crystal_constant,gaus_constant));

# Construct plot frame in 'x'
#xframe = x.frame(Title("Crystalball p.d.f."))
xframe = x.frame()


# ---------------------------------------
# Plot model and change parameter values
# ---------------------------------------

# Change the values of the parameters
crystal_mean.setVal(1.776)
crystal_sigma.setVal(0.02)
alpha.setVal(1.75)
n.setVal(3)

gaus_mean.setVal(1.776)
gaus_sigma.setVal(0.01)

# ------------------------------
# Get points from the histogram
# ------------------------------

tauMassMC = RooDataHist("tauMassMC", "Tau Mass from MC", RooArgList(x), h0)

# ------------------------------
# Fit model on MC (signal)
# ------------------------------

crystalBallFitMass = model.fitTo(tauMassMC)

#------
# Draw
#------
tauMassMC.plotOn(xframe)
model.plotOn(xframe)
#crystalball.plotOn(xframe,LineColor(kRed))

pad1.cd()
xframe.Draw()

pad2.cd()
rp0 = ROOT.TRatioPlot(h0)
rp0.Draw("e")

c.SaveAs("TauMassFit.pdf")
