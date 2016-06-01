#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ostream>
#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TFile.h"
#include "TObjString.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

 float calculateerror(float slope, float error){
   float life = (10.0/slope)*log(2);
   float lifeup = (10.0/(slope-error))*log(2);
   return abs(lifeup-life);
 }




void energy(string Finput="BipoEnergySytematics_FC.root",  int lwth=72){
// This function produces a distribution of the measured lifetime from samples generated with different energy cuts.
//Inputs: Finput = contains the distributions with different BiPo samples, lwth = lower threshold of the fitting range. 
//Outputs: A root file with a histogram of the lifetime distribution

   float upth = 1000;

   TFile *f1 = new TFile(Finput.c_str());
   TH1F* lifes = new TH1F("slopes", "slopes", 800, 290, 300);
  
   TIter next(f1->GetListOfKeys()); 
   TKey *key; 
   while ((key=(TKey*)next())) { 
      
      string namehist =key->GetName();
      TH1F* h1 = (TH1F*)f1->Get(namehist.c_str());
      h1->Draw();

   TF1 *fex = new TF1("fex","TMath::Exp([1]*x + [0]) + [2]", lwth, upth);
   h1->Fit("fex", "Q", "", lwth, upth);
   h1->Fit("fex", "LQ", "", lwth, upth);

   float slope = fex->GetParameter(1);
   float life = -10.0/(slope)*log(2);
    
   lifes->Fill(life);

   }

   lifes->Draw();
   TFile *fenergydist = new TFile("systematics_sampleSelection.root ", "RECREATE");
   lifes->Write("slifes"); 
   fenergydist->Close();
}


void sysrange(string Finput="bipo_plot250bins_FC.root", string outfilename="lowthreshold_analysis.root"){
//Produces a TGraph of lifetime measurements as a funciton of lower threshold limit selection.
//inputs: 1)Finput: it's the distribution (TH1F) of the time separation of the two largest S1 peaks of the event. 2) outfilename file to store the TGraph produced.



    TFile* f1 = new TFile(Finput.c_str());
    TH1F* h1 = (TH1F*)f1->Get("htau");
    TGraphErrors *gslopes = new TGraphErrors();
    h1->Rebin(2);
    //a_edgelw[] = array with lower threshold limits 
    float a_edgelw[] = {32, 40, 48, 56, 64,72,80, 88, 96, 104, 112, 120, 128, 132,140, 148, 156, 164, 170, 178};
    int points = sizeof(a_edgelw)/sizeof(int);

    for(int i=0; i<points; i++){
        float edgelw  = a_edgelw[i];
        float edgehi = 1000;  
        TF1 *fex = new TF1("fex","TMath::Exp([1]*x + [0]) + [2]", edgelw, edgehi);
        h1->Fit("fex", "Q", "", edgelw, edgehi);
        h1->Fit("fex", "LQ", "", edgelw, edgehi);
        float slope = fex->GetParameter(1);
        float err = fex->GetParError(1);
        float life = -10.0/(slope)*log(2);
           
        gslopes->SetPoint(i, edgelw, life);
        gslopes->SetMarkerStyle(8);
        gslopes->SetMarkerSize(0.8);
        gslopes->SetPointError(i,0, calculateerror(slope,err));
    }
  
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    gslopes->Draw("a l p");
    gslopes->GetXaxis()->SetTitle("Low threshold of the fitting range (samples)");
    gslopes->GetYaxis()->SetTitle("Fitted half-life (ns)");
    gslopes->Write("g1");
    fout->Close();
}
 




