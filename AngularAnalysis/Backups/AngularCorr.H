#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"

class AngularCorr{

public:
	AngularCorr();
	TH2D* GetData(Double_t, Double_t, Double_t, Double_t);
	vector<Double_t> cutHisto(TH2D*, Double_t, Double_t);
	void printCounts(vector<Double_t>);
	vector<Double_t> WeightAdjust(vector<Double_t>);
	void AngleFold(vector<Double_t>);



private:
	
};

