#include "TH1F.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "THnSparse"

class angPrep{
  
public:
  
  angPrep();
  void GetData();
  void InitHisto();
  void DrawHisto();
  void FillHisto();

private:

	THnSparse * dataFull;
	TH3D * xyzProj;
	TH2D * ggMatrix;
};
