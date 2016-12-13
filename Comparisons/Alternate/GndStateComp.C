#include "TF1.h"
#include <iostream>
#include "GndStateComp.H"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

GndStateComp::GndStateComp(){

	Double_t a12=0.01131308;	//-3 Ground State
	Double_t a14=0;
	Double_t scale1=5.42116;
	Double_t a22=0.0331156;		//-4 Ground State
	Double_t a24=0;
	Double_t scale2=5.4595;

	PlotFull(a12,a14,scale1,a22,a24,scale2);
}

//**********************************************************************************
//Generates plot of the two angular correlations
//**********************************************************************************
void GndStateComp::PlotFull(Double_t a12, Double_t a14, Double_t scale1, Double_t a22, Double_t a24, Double_t scale2){

		TF1* state1=new TF1("state1","(1+1/2*[0]*(3*x*x-1)+1/8*[1]*(35*x*x*x*x-30*x*x+3))*[2]",-1,1);
		state1->SetParameter(0,a12);
		state1->SetParameter(1,a14);
		state1->SetParameter(2,scale1);
		state1->SetParName(0,"a2");
		state1->SetParName(1,"a4");
		state1->SetParName(2,"scale");

		TF1* state2=new TF1("state2","(1+1/2*[0]*(3*x*x-1)+1/8*[1]*(35*x*x*x*x-30*x*x+3))*[2]",-1,1);
		state2->SetParameter(0,a22);
		state2->SetParameter(1,a24);
		state2->SetParameter(2,scale2);
		state2->SetParName(0,"a2");
		state2->SetParName(1,"a4");
		state2->SetParName(2,"scale");



		TCanvas* c1=new TCanvas("c1","c1",800,800);
		state1->Draw();
		state1->SetMinimum(5);
		state1->SetMaximum(6);
		state1->SetTitle("Comparison of Fitted GGAC Curves for Various Ground States of ^{32}Ne;Cos(#theta);Normalized Counts");
		state1->SetLineColor(4);
		state2->SetLineColor(8);
		state1->SetLineStyle(7);
		state2->SetLineStyle(7);
		state1->GetXaxis()->CenterTitle();
		state1->GetYaxis()->CenterTitle();
		state1->GetYaxis()->SetTitleOffset(1.4);
		state1->GetXaxis()->SetNdivisions(10);
		state1->GetYaxis()->SetNdivisions(10);


		state2->Draw("SAME");

		leg=new TLegend(0.6,0.75,0.85,0.85);
		leg->SetHeader("Ground State");
		leg->AddEntry("state1","-3","l");
		leg->AddEntry("state2","-4","l");
		leg->Draw("SAME");


		c1->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/");
		path+="gndcomp";
		path+=".pdf";
		cout << path.Data() << endl;
		c1->SaveAs(path);

}
