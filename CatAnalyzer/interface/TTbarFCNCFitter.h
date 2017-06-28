#ifndef CATTools_CatAnalyzer_FCNCFitter 
#define CATTools_CatAnalyzer_FCNCFitter 

#include "TLorentzVector.h"
#include "TMinuit.h"
#include <iostream>

#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

namespace fcnc{

  double JetEResolution         (double energy);
  double METPhiResolution       (double met);
  double METResolution          (double met);
  double TwoObjectMassResolution(TLorentzVector &j1, double releres1, TLorentzVector &j2, double releres2);
  double TwoJetMassResolution   (TLorentzVector &j1, TLorentzVector &j2);// relative mass resolution
  double SolvettbarLepJets(double &nupz, Double_t &metscale, Double_t &blscale, Double_t &bjscale, Double_t &j1scale, Double_t &j2scale);
  void   FindHadronicTop        (TLorentzVector &lepton, std::vector<cat::Jet>& jets, const cat::MET & met, bool, int * csvid, std::vector<int> &bestindices, float &bestchi2, TLorentzVector &nusol, TLorentzVector &blrefit, TLorentzVector &bjrefit, TLorentzVector &j1refit, TLorentzVector &j2refit);

  void fcnfull(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  void fcn    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  void InitMinuit();

  TMinuit *tm = 0;

  TLorentzVector tmplep, tmpnu, tmpbl, tmpbj, tmpj1, tmpj2;
  float blres, bjres, j1res, j2res, metres;

}

#endif
