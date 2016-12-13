
#include "TH1F.h"
#include <iostream>
#include "Spectra.H"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

Spectra::Spectra(){
  
  GetData();
//  InitHisto();
//  FillHisto();
  DrawHisto();

}

void Spectra::GetData(){
		Int_t nbins=2000;
		Int_t EHigh=2000;
  
	TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Converted.root");

	enSpec= new TH1D("enSpec","Energy Spectrum ^{60}Co;Energy (keV);Counts",nbins,0,EHigh);
	enSpec=(TH1D*)isoData->Get("Griffin1D/griffin_crystal_unsup_edep_cry");

}

void Spectra::InitHisto(){
		Int_t nbins=2000;
		Int_t EHigh=2000;
//	TFile histfile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/spectra.root","RECREATE");
	enSpec= new TH1D("enSpec","Energy Spectrum ^{32}Ne;Energy (keV);Counts",nbins,0,EHigh);
//	TCut c_attempt("");
//	TString v_attempt("depEnergy");
//	tspecData->Project("enSpec",v_attempt,c_attempt);
//	histfile.Write();

}


void Spectra::FillHisto(){
	TCut c_attempt("");
	TString v_attempt("griffin_crystal_unsup_edep_cry");

	enSpec=(TH1D*)isoData->Get("Griffin1D/griffin_crystal_unsup_edep_cry");
 //stay away from for loops as much as possible. 
 // for(Double_t i=0; i<nwwlln; i++){
 //   twwlln->GetEntry(i);
 //   hlatitude_wwlln->Fill(latitude_wwlln);
 // }

}

void Spectra::DrawHisto(){
//	TTree* tree=new TTree("Spectra","Spectra");
//	tree->Branch();


	TCanvas c("c","c",600,600);
	c.SetLogy();
	gStyle->SetOptStat(0);
	enSpec->GetXaxis()->CenterTitle();
	enSpec->GetYaxis()->CenterTitle();
	enSpec->GetXaxis()->SetRangeUser(0,2000);
	enSpec->SetTitle(";Energy (keV);Counts");
	enSpec->GetYaxis()->SetTitleOffset(1.4);

  enSpec->Draw();
  c.Update();
  c.SaveAs("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Ne32Spectra.pdf");


	TFile histfile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/spectra.root","RECREATE");
	histfile.WriteTObject(enSpec);
}


