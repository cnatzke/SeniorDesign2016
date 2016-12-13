//This Macro SHOULD compare various ground state angular correlations via other root files.

#include <iostream>
#include <vector>
#include "GndStateCompAlt.H"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TMultiGraph.h"

using namespace std;

GndStateCompAlt::GndStateCompAlt(){
	GetData();
}

void GndStateCompAlt::GetData(){
  
	TFile* ggac0mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac0mm.root");
	TFile* ggac10mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac10mm.root");
	TFile* ggac20mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac20mm.root");
	TFile* ggac30mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac30mm.root");
	TFile* ggac40mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac40mm.root");
	TFile* ggac50mm = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac50mm.root");

	hist0mm=(TGraphErrors*)ggac0mm->Get("Graph");
	hist10mm=(TGraphErrors*)ggac10mm->Get("Graph");
	hist20mm=(TGraphErrors*)ggac20mm->Get("Graph");
	hist30mm=(TGraphErrors*)ggac30mm->Get("Graph");
	hist40mm=(TGraphErrors*)ggac40mm->Get("Graph");
	hist50mm=(TGraphErrors*)ggac50mm->Get("Graph");
		
		TCanvas* c1=new TCanvas("c1","c1",800,800);
		TMultiGraph* mg=new TMultiGraph();
		mg->Add(hist0mm);
		mg->Add(hist10mm);
		mg->Add(hist20mm);
		mg->Add(hist30mm);
		mg->Add(hist40mm);
		mg->Add(hist50mm);

		mg->Draw("ap");
		mg->SetMinimum(110);
		mg->SetMaximum(150);
		mg->SetTitle("Comparison Particle Beam Diameter for ^{60}Co;Cos(#theta);Normalized Counts");
		hist0mm->SetLineColor(1);
		hist10mm->SetLineColor(2);
		hist20mm->SetLineColor(3);
		hist30mm->SetLineColor(4);
		hist40mm->SetLineColor(9);
		hist50mm->SetLineColor(6);
		mg->GetXaxis()->CenterTitle();
		mg->GetYaxis()->CenterTitle();
		mg->GetYaxis()->SetTitleOffset(1.4);
		mg->GetXaxis()->SetNdivisions(10);
		mg->GetYaxis()->SetNdivisions(5);



		leg=new TLegend(0.6,0.75,0.85,0.85);
		leg->SetHeader("Beam Radius");
        leg->SetNColumns(2);
		leg->AddEntry(hist0mm,"Point Source","ep");
		leg->AddEntry(hist10mm,"10 mm","ep");
		leg->AddEntry(hist20mm,"20 mm","ep");
		leg->AddEntry(hist30mm,"30 mm","ep");
		leg->AddEntry(hist40mm,"40 mm","ep");
		leg->AddEntry(hist50mm,"50 mm","ep");
		leg->Draw("SAME");


		c1->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/");
		path+="beamcomp";
		path+=".pdf";
		cout << path.Data() << endl;
		c1->SaveAs(path);
}

