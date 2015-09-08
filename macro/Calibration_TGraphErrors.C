//Usage
//Help from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookHowToFit
//===============================================================================
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"

string IntToString (int number1)
{
  std::ostringstream oss1;
  oss1<< number1;
  return oss1.str();
}

//===============================================================================
Double_t langaufun(Double_t *x, Double_t *par) {
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  Double_t mpshift = 0;//-0.22278298; // Landau maximum location
  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  return (par[2] * step * sum * invsq2pi / par[3]);
};

//===============================================================================
Double_t langaufungaus(Double_t *x, Double_t *par) {
return langaufun(x,par)+par[4]*TMath::Gaus(x[0],par[5],par[6],true);
};

//===============================================================================
void Calibration(){

  //===============================================================================
  gStyle->SetOptTitle(0);

  //The root files with the individual plots 
  const int numberoffiles = 4;
  const int numberofchambers = 6;

  TString filename[numberoffiles];
  //mu-
  filename[0]="testcalo_mu-_50.root";
  filename[1]="testcalo_mu-_100.root";
  filename[2]="testcalo_mu-_150.root";
  filename[3]="testcalo_mu-_200.root";

  //Energy distributions we want to fit 
  TString energy[numberofchambers];
  //nosmear
  // energy[0] = "histo/6";
  // energy[1] = "histo/7";
  // energy[2] = "histo/8";
  // energy[3] = "histo/9";
  // energy[4] = "histo/10";
  // energy[5] = "histo/11";

  //smeared
  energy[0] = "histo/18";
  energy[1] = "histo/19";
  energy[2] = "histo/20";
  energy[3] = "histo/21";
  energy[4] = "histo/22";
  energy[5] = "histo/23";

  double mpshift = 0;//-0.22278298;

  TGraphErrors *grMPV[numberoffiles];
  TString mpvname[numberoffiles];
  
  TGraphErrors *grSigma[numberoffiles];
  TString sigmaname[numberoffiles];

  //For the mpv and sigma plot 
  std::vector<double> chamb;
  std::vector<double> mpv;
  std::vector<double> thesigma;
  std::vector<double> chamr;
  std::vector<double> mpvr;
  std::vector<double> thesigmaerr;
  chamb.push_back(1.);chamb.push_back(2.);chamb.push_back(3.);chamb.push_back(4.);chamb.push_back(5.);chamb.push_back(6.);
  
  //======================================================================================
  //Canvas
  TCanvas *can[numberoffiles][numberofchambers];
  TString canname;

  TCanvas *canmpv[numberoffiles];
  //======================================================================================
  //Read the files and histograms here
  TFile* files[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    files[k]= new TFile(filename[k]);
  }
  TH1F * histo[numberoffiles][numberofchambers];

  for (int i=0; i<numberoffiles; i++){
    gROOT->Reset();
    for (int k=0; k<numberofchambers; k++){
      histo[i][k]= (TH1F*) files[i]->Get(energy[k]);

      if (!histo[i][k]){std::cerr << "Could not get histogram " << energy[k] << "in file "
				    <<  files[i]->GetName() << std::endl;}
      std::cout << "=========================================================================" << std::endl;
      std::cout << files[i]->GetName() << std::endl;

      canname = filename[i] + histo[i][k]->GetName();
      can[i][k] = new TCanvas(canname, canname,800,600);
      can[i][k]->cd();

      histo[i][k]->Draw();

      TF1 *fit = 0;
      double mpc = 0;
      double minfit = 0.;
      double maxfit = 0.05;

      histo[i][k]->Fit("landau","R+","same",minfit,maxfit);
      fit = (TF1*) histo[i][k]->GetFunction("landau");

      mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 

      mpv.push_back(mpc);
      mpvr.push_back(fit->GetParError(1));
      chamr.push_back(0.);

      thesigma.push_back(fit->GetParameter(2));
      thesigmaerr.push_back(fit->GetParError(2));

      // if(i==0) {mpvname = "grMPV_50GeV"; sigmaname = "grSigma_50GeV";}
      // if(i==1) {mpvname = "grMPV_100GeV";sigmaname = "grSigma_100GeV";}
      // if(i==2) {mpvname = "grMPV_150GeV";sigmaname = "grSigma_150GeV";}
      // if(i==3) {mpvname = "grMPV_200GeV";sigmaname = "grSigma_200GeV";}
      mpvname[i] = "grMPV_" + filename[i]; 
      sigmaname[i] = "grSigma" + filename[i];
      
      // std::cout << "Chamber " << k << " MPV "<< mpc << std::endl;

      can[i][k]->Print(canname + ".root",".root");
      can[i][k]->Close();
  
      
    }// end of loop over histos   

    TGraphErrors *gr = new TGraphErrors(mpv.size(),&chamb[0],&mpv[0],&mpvr[0],&chamr[0]);
    gr->SetMarkerStyle(20); 
    gr->SetMarkerSize(1.0); 
    // gr->GetYaxis()->SetRangeUser(0.,100.);
    //gr->GetXaxis()->SetRangeUser(0.,300.);
    gr->SetTitle(mpvname[i]);
    gr->GetHistogram()->SetXTitle(" Chamber  ");
    gr->GetHistogram()->SetYTitle(" Landau MPV (MIPs)  ");
    gr->GetYaxis()->SetTitleOffset(1.6);
    TCanvas *c1 = new TCanvas("c1","mpv",200,10,700,500);
    gr->Draw("APL");
    c1->Update();

    c1->Print(mpvname[i]+".png",".png");
    c1->Close();

    mpv.clear();
    mpvr.clear();
    
    TGraphErrors *grs = new TGraphErrors(thesigma.size(),&chamb[0],&thesigma[0],&thesigmaerr[0],&chamr[0]);
    grs->SetMarkerStyle(20); 
    grs->SetMarkerSize(1.0); 
    // grs->GetYaxis()->SetRangeUser(0.,100.);
    //grs->GetXaxis()->SetRangeUser(0.,300.);
    grs->SetTitle(sigmaname[i]);
    grs->GetHistogram()->SetXTitle(" Chamber  ");
    grs->GetHistogram()->SetYTitle(" Sigma ");
    grs->GetYaxis()->SetTitleOffset(1.6);
    TCanvas *c2 = new TCanvas("c2","sigma",200,10,700,500);
    grs->Draw("APL");
    c2->Update();

    c2->Print(sigmaname[i]+".png",".png");
    c2->Close();

    thesigma.clear();
    thesigmaerr.clear();
    chamr.clear();


    gr->Delete();
    c1->Close();
    c1->Close();


  }// end of loop over files
  
  
  


}

