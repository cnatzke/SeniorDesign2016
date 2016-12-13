#include "TH1F.h"
#include <iostream>
#include "SpectraCompare.H"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

using namespace std;

SpectraCompare::SpectraCompare(){
  
  GetData();
  InitHisto();
  FillHisto();
  DrawHisto();

}

void SpectraCompare::GetData(){
  
  TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/g4out.root");
  isoData->cd("ntuple");
  tspecData = (TTree*)gDirectory->Get("ntuple");
  tspecData->SetBranchAddress("depEnergy",&energy);

}

void SpectraCompare::InitHisto(){
		Int_t nbins=1000;
		Int_t EHigh=2000;
//	TFile histfile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/spectra.root","RECREATE");
	enSpec= new TH1D("enSpec","Energy Spectrum ^{32}Ne;Energy (keV);Counts",nbins,0,EHigh);
//	TCut c_attempt("");
//	TString v_attempt("depEnergy");
//	tspecData->Project("enSpec",v_attempt,c_attempt);
//	histfile.Write();

}


void SpectraCompare::FillHisto(){
	TCut c_attempt("");
	TString v_attempt("depEnergy");

	tspecData->Project("enSpec",v_attempt,c_attempt);

 //stay away from for loops as much as possible. 
 // for(Double_t i=0; i<nwwlln; i++){
 //   twwlln->GetEntry(i);
 //   hlatitude_wwlln->Fill(latitude_wwlln);
 // }

}

void SpectraCompare::DrawHisto(){
//	TTree* tree=new TTree("Spectra","Spectra");
//	tree->Branch();


//	TCanvas c("c","c",600,600);
//	c.SetLogy();
//	gStyle->SetOptStat(10);
//	enSpec->GetXaxis()->CenterTitle();
//	enSpec->GetYaxis()->CenterTitle();
  
//  enSpec->Draw();
//  c.Update();
//  c.SaveAs("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Ne32Spectra.pdf");

	TFile* compData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/spectra.root");

		TCanvas* c1=new TCanvas("c1","c1",800,800);
		c1->SetLogy();
		gStyle->SetOptStat(0);

		enSpec->Draw();
		griffin_crystal_unsup_edep_cry->Draw("SAME");
		enSpec->SetTitle("Beta Filted Spectrum Comparison for ^{32}Ne;Energy (keV);Counts");
		enSpec->SetLineColor(4);
		griffin_crystal_unsup_edep_cry->SetLineColor(8);
		enSpec->GetXaxis()->CenterTitle();
		enSpec->GetYaxis()->CenterTitle();
		enSpec->GetYaxis()->SetTitleOffset(1.4);



		leg=new TLegend(0.6,0.75,0.85,0.85);
		leg->SetHeader("Energy Spectra");
		leg->AddEntry("enSpec","Unfiltered","l");
		leg->AddEntry("griffin_crystal_unsup_edep_cry","Filtered","l");
		leg->Draw("SAME");


		c1->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/SpectraComp.pdf");
		c1->SaveAs(path);

}

