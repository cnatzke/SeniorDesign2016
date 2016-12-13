#include "TH1F.h"
#include <iostream>
#include "angPrep.H"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

angPrep::angPrep(){
  GetData();
  InitialHisto();
  GGMatrixProjection();
  DrawHisto();

}

void angPrep::GetData(){
  
	TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Converted.root");
	isoData->cd("GriffinND");
	dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_addback_sparse");
}

void angPrep::InitialHisto(){
	xyzProj=(TH3D*) dataFull->Projection(0,1,2);

}


void angPrep::GGMatrixProjection(){
	ggMatrix=(TH2D*) dataFull->Projection(0,1);

}

void angPrep::DrawHisto(){
	TFile histfile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Macros/Step1.root","RECREATE");
   	TCanvas c("c","c",600,600);
  	c.SetLogz();
  	gStyle->SetOptStat(10);
  	ggMatrix->SetTitle("#\gamma #\gamma Correlation Matrix;Energy (keV); Energy (keV)");
  	ggMatrix->GetXaxis()->CenterTitle();
  	ggMatrix->GetYaxis()->CenterTitle();
	  
	ggMatrix->Draw("colz");
  	c.Update();
  	c.SaveAs("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Spectra/ggMatrix.pdf");
  	histfile.Write();


}
