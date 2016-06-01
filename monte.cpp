#include <TCanvas.h>
#include <string>
#include "TMath.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "TCut.h"
#include <sstream>
#include "TROOT.h"
#include <sstream>
#include <fstream>
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFormula.h"
#include "TEventList.h"
#include "TFile.h"
#include "TObjString.h"
#include <ostream>
#include <utility>
#include <algorithm>
#include "TFriendElement.h"
#include "TKey.h"
#include "TCollection.h"
#include "TRandom3.h"

using namespace std;

//Decay Model: TMath::Exp([1]*x + [0]) + [2]"
//Input values for simulation extracted from best fit. 
float po212par0 =  1.13848e+01;
float po212par1 = -2.35868e-02;
float po212par2 = 4.24153e+00;


float live(int ev, int bins=125, int lwth=72){
//This function generates a distribution according to the radiactive decay model and returns the lifetime from the distribution
//Input parameters are: ev = number of events to be drawn, binning = bins in (0,1000) samples, lwth = lower threshold of the range
//Events are being randomly drawn in the range (lwth-1000) samples

   
  // po212L po212L1;
   TF1 *fex = new TF1("fex","TMath::Exp([1]*x + [0]) + [2]", lwth, 1000); 
   TF1 *fit = new TF1("fit","TMath::Exp([1]*x + [0]) + [2]", lwth, 1000);
   fex->SetParameter(0,po212par0);
   fex->SetParameter(1,po212par1);
   fex->SetParameter(2,po212par2);

   TH1F * h1 = new TH1F("h1", "h1",bins , 0, 1000);
   double event;
   for(int i=0; i< ev; i++){
       event = fex->GetRandom();
       h1->Fill(event);
    }
    h1->Fit("fit", "Q", "", lwth, 1000);
    h1->Fit("fit", "Q l", "", lwth, 1000);
    float slope = fit->GetParameter(1);

    delete fex;
    delete fit;
    delete h1;
    return slope;
}


void slives(int mciter, int  evs, int bins=125, int lwth=72){
//Simulates m distributions of po212 decay, every distribution with n events. Extract lifetime from the distribution and produces a distribution of the lifetimes.
//Input varibales: m = mciter, n=evs, bins = number of bins in the range (0,1000) samples, lwth = Lower threshold limit of the range (lwth,1000).

   TH1F * hlives = new TH1F("hlives", "hlives", 400, 270, 320);
   for(int i=0; i< mciter; i++){
      float half = (-10.0/(live(evs, bins, lwth)))*TMath::Log(2);
      hlives->Fill(half);
   }

   TFile *fout = new TFile(Form("montecarlo_bins%i.root", bins), "RECREATE");
   hlives->Write();
   hlives->Draw();
   fout->Close();
   }


    TH1F* histo(int evs, int bins){
    //po212L po212L1;
    TF1 *fex = new TF1("fex","TMath::Exp([1]*x + [0]) + [2]", 0, 1000); 
    fex->SetParameter(0, po212par0);
    fex->SetParameter(1, po212par1);
    fex->SetParameter(2, po212par2);
    TH1F * h1 = new TH1F("h1", "h1",bins , 0, 1000);

    double event;
    for(int i=0; i< evs; i++){
        event = fex->GetRandom();
        h1->Fill(event);
    }
 
    return h1;
}
 

