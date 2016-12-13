// Connor Natzke 
// July 19, 2016
//
// To Do:
// Plot resultant angular correlations
// subtract backgrounds

#include <iostream>
#include "AngularCorr.H"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLatex.h"

using namespace std;

AngularCorr::AngularCorr(){
	Double_t gateLow=1329;	//First Gate
	Double_t gateHigh=1336;	//First Gate
	Double_t ELow=1160;		//Low Energy Gate
	Double_t EHigh=1180;	//High Energy Gate

	TH2D * yzProj=GetData(gateLow, gateHigh, ELow, EHigh);
	vector<Double_t> counts=cutHisto(yzProj, ELow, EHigh);
	vector<Double_t> wcounts=WeightAdjust(counts);
//	printCounts(wcounts);
	AngularCorrHisto(wcounts);
}

//**********************************************************************************
TH2D* AngularCorr::GetData(Double_t gateLow, Double_t gateHigh, Double_t ELow, Double_t EHigh){
  
	TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Converted.root");
	isoData->cd("GriffinND");
	THnSparse * dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH3D * xyzProj=(TH3D*) dataFull->Projection(0,1,2);	//3D Histogram
	xyzProj->GetXaxis()->SetRangeUser(gateLow,gateHigh); //First Gate

	TH2D * yzProj=(TH2D*) xyzProj->Project3D("zy");

	return yzProj;
}

//********************************************************
vector<Double_t> AngularCorr::cutHisto(TH2D* angHist,Double_t ELow, Double_t EHigh){
	
	Int_t angSplits=52; //number of angular indices

	TFile HISTFILE("Angular_Splits.root","RECREATE");
	vector<TH1D*> histVec(angSplits); //Initialize Vector of Histograms
	vector<Double_t> counts(angSplits); //Initalize Vector for Counts

	//Initialize Canvas
	TCanvas* c = new TCanvas("c","c",800,800);  
	gStyle->SetOptStat(10);
	c->SetLogy();

	for(Int_t i=0;i<angSplits;i++){ //loop through angular indices

		Double_t indLow=i; //Starts at bin0
		Double_t indHigh=indLow+1; //cut at 1 bin higher

       //Histogram Names
		TString h_name("h_angIndex_");
 		h_name+=indLow;
		cout << h_name.Data() << endl;
		
       //histogram titles
		TString h_title("Gated Counts for Index ");
		h_title+=indLow;
		h_title+=";Energy (keV);Counts";
		cout << h_title.Data() << endl;


		delete histVec[i]; //clear memory space
		histVec[i]=angHist->ProjectionX(h_name,indLow,indHigh);
		histVec[i]->SetTitle(h_title);


	   //Draw Histogram
		histVec[i]->Draw();
		TF1 * gFit=new TF1("gFit","gaus",ELow, EHigh);			//Defining Fit
		histVec[i]->GetXaxis()->SetRangeUser(ELow-5,EHigh+5);	//Setting Range of Histogram to allow data fitting
		histVec[i]->Fit("gFit","Q");							//Fitting
		counts[i]=gFit->Integral(ELow,EHigh);					//Integrating Peak to Find Counts
//		histVec[i]->GetXaxis()->SetRangeUser(0,2500);			//Resetting Range for Cosmetic Purposes
		gFit->Draw("SAME");

	   //Save Histogram
		c->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Angular_Histograms/");
		path+=h_name;
		path+=".png";
		cout << path.Data() << endl;
		c->SaveAs(path);

	}

	HISTFILE.Write();
	return counts;
}

//*****************************************************************************
void AngularCorr::printCounts(vector<Double_t> counts){
	for(Int_t i=0;i<52;i++){
			cout << counts[i] << " ";
	}
	cout << endl;
}
//*****************************************************************************
vector<Double_t> AngularCorr::WeightAdjust(vector<Double_t> counts){

	vector<Double_t> wcounts(52);

	for(Int_t i=0;i<52;i++){
		if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50){
				wcounts[i]=counts[i]/(128);
		}
		else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46){
				wcounts[i]=counts[i]/(48);
		}
		else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44){
				wcounts[i]=counts[i]/(96);
		}
		else{
				wcounts[i]=counts[i]/(64);
		}
//		cout << wcounts[i] << " ";
	}
//	cout << endl;
	return wcounts;
}
//*****************************************************************************
void AngularCorr::AngularCorrHisto(vector<Double_t> wcounts){
	//Angular Indices
		Double_t index[]={0.0,18.79097,25.60153, 26.69036,31.94623,33.65414, 44.36426, 46.79372, 48.57554, 49.79788, 53.83362, 60.15106, 62.70487, 63.08604, 65.01569, 66.46082, 67.45617, 69.86404, 70.86009, 73.08384, 76.38138, 78.66898, 83.04252, 86.22840, 86.23761, 88.47356, 91.52644, 93.76239, 93.77160, 96.95749, 101.33102, 103.61822, 106.91616, 109.13991, 110.13596, 112.54383, 113.53918, 114.98431, 116.91396, 117.29513, 119.84894, 126.16638, 130.20212, 131.42446, 133.20628, 135.63574, 146.34586, 148.05377, 153.30964, 154.39847, 161.21315,180.0};

		vector<Double_t> indexVec(52);

		//Fold Angles and convert to vector
		for(Int_t i=0;i<=52;i++){
		if(index[i]>90){
				index[i]=180-index[i];
				indexVec[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
			}
		else{
				indexVec[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
			}
		}

	vector<Double_t> indexAdd(26);
	vector<Double_t> wcountsAdd(26);
	for(Int_t j=0;j<=25;j++){
			indexAdd[j]=indexVec[j];
			wcountsAdd[j]=wcounts[j]+wcounts[52-j];
//			cout << wcountsAdd[j] << " ";
	}
//	cout << endl;

	for(Int_t k=0;k<=25;k++){
			cout << indexAdd[k] << " ";
	}
	cout << endl;

	TCanvas * c1=new TCanvas("c1","c1",600,600);

	TGraph * g=new TGraph(26,&indexAdd[0],&wcountsAdd[0]);
	g->SetTitle("Angular Correlation;Cos(#theta);Counts");
	g->SetMarkerStyle(7);
	g->Draw("AP");


	   //Save Histogram
		c1->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/");
		path+="corr";
		path+=".pdf";
		cout << path.Data() << endl;
		c1->SaveAs(path);
}
