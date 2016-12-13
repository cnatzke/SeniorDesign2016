#include "TH1F.h"
#include <iostream>
#include "FirstAnalysis.H"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

FirstAnalysis::FirstAnalysis(){
  
  GetData();
  InitHisto();
  FillHisto();
  DrawHisto();

}

void FirstAnalysis::GetData(){
  
  TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/g4out.root");
  isoData->cd("ntuple");
  tspecData = (TTree*)gDirectory->Get("ntuple");
  tspecData->SetBranchAddress("depEnergy",&energy);

}

void FirstAnalysis::InitHisto(){
		Int_t nbins=1000;
		Int_t EHigh=2000;
//	TFile histfile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/spectra.root","RECREATE");
	enSpec= new TH1D("enSpec","Energy Spectrum ^{32}Ne;Energy (keV);Counts",nbins,0,EHigh);
//	TCut c_attempt("");
//	TString v_attempt("depEnergy");
//	tspecData->Project("enSpec",v_attempt,c_attempt);
//	histfile.Write();

}


void FirstAnalysis::FillHisto(){
	TCut c_attempt("");
	TString v_attempt("depEnergy");

	tspecData->Project("enSpec",v_attempt,c_attempt);

 //stay away from for loops as much as possible. 
 // for(Double_t i=0; i<nwwlln; i++){
 //   twwlln->GetEntry(i);
 //   hlatitude_wwlln->Fill(latitude_wwlln);
 // }

}

void FirstAnalysis::DrawHisto(){
//	TTree* tree=new TTree("Spectra","Spectra");
//	tree->Branch();


	TCanvas c("c","c",600,600);
	c.SetLogy();
	gStyle->SetOptStat(10);
	enSpec->GetXaxis()->CenterTitle();
	enSpec->GetYaxis()->CenterTitle();
  
  enSpec->Draw();
  c.Update();
  c.SaveAs("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Ne32Spectra.pdf");


}


