#include "SimpleFits/FitSoftware/interface/TrackParticle.h"

TrackParticle::TrackParticle():
  Particle(TMatrixT<double>(NHelixPar,1),TMatrixTSym<double>(NHelixPar),0,0,0),
  mass(0)
{
}

TrackParticle::TrackParticle(const TrackParticle& another):
  Particle(another.getParMatrix(), another.getCovMatrix(), another.PDGID(), another.Charge(), another.BField()),
  mass(another.Mass())
{
}

TrackParticle::TrackParticle(TMatrixT<double> par_, TMatrixTSym<double> cov_, int pdgid_,double mass_, double charge_, double b_):
  Particle(par_,cov_,pdgid_,charge_,b_),
  mass(mass_)
{
}

TString TrackParticle::Name(int i){
  if(i==kappa)  return "kappa";
  if(i==lambda) return "lambda";
  if(i==phi)    return "phi";
  if(i==dz)     return "dz";
  if(i==dxy)    return "dxy";
  return "invalid";
}
