#------------------------------------------------------------------
#Script to combine the plots from all the Sets and fit the tau mass
#------------------------------------------------------------------

import ROOT
from ROOT import RooFit
from ROOT import TFile, TH1F, TLine, TCanvas
from ROOT import RooRealVar, RooCBShape, RooGaussian, RooAddPdf, RooDataHist, RooArgSet, RooArgList, RooAbsArg
from tdrstyle import setTDRStyle

setTDRStyle() # use cms tdrstyle

c = TCanvas("c","Tau Mass Fit",1600,1200)

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
tauMassMC.plotOn(xframe)
model.plotOn(xframe)
#crystalball.plotOn(xframe,LineColor(kRed))
xframe.Draw()

c.SaveAs("TauMassFit.pdf")
