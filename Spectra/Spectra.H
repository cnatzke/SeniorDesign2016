
#include "TH1F.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

class Spectra{
  
public:
  
  Spectra();
  void GetData();
  void InitHisto();
  void DrawHisto();
  void FillHisto();

private:

	TTree* tspecData;
	Double_t nwwlln;
	Double_t energy;
	TH1D * enSpec;
};
