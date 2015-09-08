//===============================================================================
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"

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
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   // his->Fit(FunName,"RB0","same");   // fit within specified range, use ParLimits, do not plot
   his->Fit(FunName,"R+","same");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}




//===============================================================================
void AbsorberStudy(){

  //===============================================================================
  gStyle->SetOptTitle(0);

  //The root files with the individual plots 
  const int numberoffiles = 1;
  const int numberofabsorbers = 1;

  TString filename[numberoffiles];
  //mu-
  filename[0]="Scream_50cmFeAbsorber.root";
  // filename[1]="testcalo_mu-_100.root";
  // filename[2]="testcalo_mu-_150.root";
  // filename[3]="testcalo_mu-_200.root";

  //Energy distributions we want to fit 
  TString energy[numberofabsorbers];
  //nosmear
  // energy[0] = "histo/6";
  // energy[1] = "histo/7";
  // energy[2] = "histo/8";
  // energy[3] = "histo/9";
  // energy[4] = "histo/10";
  // energy[5] = "histo/11";

  //smeared
  energy[0] = "histo/33";

  double mpshift = 0;//-0.22278298;

  // TGraphErrors *grMPV[numberoffiles];
  // TGraphErrors *grSigma[numberoffiles];
  
  // //For the mpv and sigma plot 
  // std::vector<double> chamb;
  // std::vector<double> mpv;
  // std::vector<double> thesigma;
  // std::vector<double> chamr;
  // std::vector<double> mpvr;
  // std::vector<double> thesigmaerr;
  // chamb.push_back(1.);chamb.push_back(2.);chamb.push_back(3.);chamb.push_back(4.);chamb.push_back(5.);chamb.push_back(6.);
  
  //======================================================================================
  //Canvas
  TCanvas *can[numberoffiles][numberofabsorbers];
  TString canname;

  TCanvas *canmpv[numberoffiles];
  
  //======================================================================================
  //Read the files and histograms here
  TFile* files[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    files[k]= new TFile(filename[k]);
  }
  TH1F * histo[numberoffiles][numberofabsorbers];

  for (int i=0; i<numberoffiles; i++){
    gROOT->Reset();
    for (int k=0; k<numberofabsorbers; k++){
      histo[i][k]= (TH1F*) files[i]->Get(energy[k]);

      if (!histo[i][k]){std::cerr << "Could not get histogram " << energy[k] << "in file "
				    <<  files[i]->GetName() << std::endl;}
      std::cout << "=========================================================================" << std::endl;
      std::cout << files[i]->GetName() << std::endl;

      canname = filename[i] + histo[i][k]->GetName();
      histo[i][k]->Scale(0.1);//This is for the 10 events I run
      can[i][k] = new TCanvas(canname, canname,800,600);
      can[i][k]->cd();

      histo[i][k]->Draw();

      // TF1 *fit = 0;
      double mpc = 0;
      double minfit = 0.;
      double maxfit = 500.;

      // histo[i][k]->Fit("landau","R+","same",minfit,maxfit);
      // fit = (TF1*) histo[i][k]->GetFunction("landau");
 
      // TF1 *fit = new TF1("langaufun",langaufun,minfit,maxfit,4);
      // fit->SetParNames("Width","MP","Area","GSigma");  
      // Double_t *startvalues;
      // sv[0]=1.0e+05; sv[1]=7.0e-04; sv[2]=2.0e-04; // sv[3]=0.04;
      // histo[i][k]->Fit(langaufun,"R+","same",minfit,maxfit);
      // fit = (TF1*) histo[i][k]->GetFunction(langaufun);

      //Setting fit range and start values
      // Double_t fr[2];
      // Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      // fr[0]=minfit; fr[1]=maxfit;
      // sv[0]=1.0e-04; sv[1]=1.0e-03; sv[2]=2.0e-04; sv[3]=0.04;
      // pllo[0]=0.000001; pllo[1]=0.0005; pllo[2]=0.000001; pllo[3]=0.000001;
      // plhi[0]=0.05; plhi[1]=0.002; plhi[2]=50; plhi[3]=0.05;
      // Double_t chisqr;
      // Int_t    ndf;

      // fit = langaufit(histo[i][k],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

      // mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 

      // mpv.push_back(mpc * 1000.);
      // mpvr.push_back(fit->GetParError(1) * 1000.);
      // chamr.push_back(0.);

      // thesigma.push_back(fit->GetParameter(2) * 1000.);
      // thesigmaerr.push_back(fit->GetParError(2) * 1000.);
   
      // std::cout << "File " << filename[i] << " Chamber " << k << " MPV "<< mpc << std::endl;

      can[i][k]->Print(canname + ".root",".root");
      //can[i][k]->Close();
  
      
    }// end of loop over histos   
    
    //For the multi plot
    // grMPV[i] = new TGraphErrors(mpv.size(),&chamb[0],&mpv[0],&mpvr[0],&chamr[0]);
    // grSigma[i] = new TGraphErrors(thesigma.size(),&chamb[0],&thesigma[0],&thesigmaerr[0],&chamr[0]);
    
    // mpv.clear();
    // mpvr.clear();
    // thesigma.clear();
    // thesigmaerr.clear();
    // chamr.clear();
    
  }// end of loop over files
  
  // const Int_t nch = 6;
  // char *chamb_labels[nch] = {"Std1","Star","Mirror","Snake","Spider","Std2"};
  // TFile f("mpvandsigma.root","recreate");

  // //======================================================================================
  // //For the multi mpv plot
  // TCanvas *c1 = new TCanvas("c1","mpv",200,10,1200,968); 
  // //Change the bin labels
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(1.0) , chamb_labels[0]); 
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(2.0) , chamb_labels[1]); 
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(3.0) , chamb_labels[2]); 
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(4.0) , chamb_labels[3]); 
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(5.0) , chamb_labels[4]); 
  // grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(6.0) , chamb_labels[5]); 
  // grMPV[0]->GetXaxis()->SetTitleOffset(2.0);
  // grMPV[0]->GetYaxis()->SetTitleOffset(1.2);
  // grMPV[0]->GetXaxis()->SetLabelSize(0.06);
  
  // grMPV[0]->SetMarkerStyle(20); 
  // grMPV[0]->SetMarkerSize(1.2); 
  // grMPV[0]->SetMarkerColor(1);  
  // grMPV[0]->SetLineColor(1);
  // grMPV[0]->Draw("APL");
  // c1->Update();
    
  // grMPV[1]->SetMarkerStyle(21);
  // grMPV[1]->SetMarkerSize(1.2); 
  // grMPV[1]->SetMarkerColor(2);
  // grMPV[1]->SetLineColor(2);
  // grMPV[1]->Draw("PSL");
  // c1->Update();

  // grMPV[2]->SetMarkerStyle(22);
  // grMPV[2]->SetMarkerSize(1.2); 
  // grMPV[2]->SetMarkerColor(3);
  // grMPV[2]->SetLineColor(3);
  // grMPV[2]->Draw("PSL");
  // c1->Update();

  // grMPV[3]->SetMarkerStyle(23);
  // grMPV[3]->SetMarkerSize(1.2); 
  // grMPV[3]->SetMarkerColor(4);
  // grMPV[3]->SetLineColor(4);
  // grMPV[3]->Draw("PSL");
  // c1->Update();
  
  // TLegend *leg = new TLegend(0.8,0.6,0.89,0.89);  //coordinates are fractions of pad dimensions
  // leg->AddEntry(grMPV[0],"50 GeV","LP");  
  // leg->AddEntry(grMPV[1],"100 GeV","LP");  
  // leg->AddEntry(grMPV[2],"150 GeV","LP");  
  // leg->AddEntry(grMPV[3],"200 GeV","LP");  
  // leg->SetFillColor(18);
  // leg->SetHeader(" #mu^{-} energies ");
                                            
  // leg->Draw("PS");
  // c1->Update();
  
  // grMPV[0]-> SetTitle( " Landau Most Probable Value for different #mu^{-} energies  " );
  // grMPV[0]->GetHistogram()->SetXTitle(" Chamber ");
  // grMPV[0]->GetHistogram()->SetYTitle(" Landau MPV (MIPs) (keV) ");
  // grMPV[0]->GetYaxis()->SetRangeUser(0.65,0.95);
  // //grMPV[0]->GetXaxis()->SetRangeUser(0.,5000.);
  // c1->Update();

  // //c1->Print("mpv.png",".png");

  // c1->Write();
  // c1->Close();

 
  // //======================================================================================
  // //For the multi sigma plot
  // TCanvas *c2 = new TCanvas("c2","sigma",200,10,1200,968); 
  // //Change the bin labels
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(1.0) , chamb_labels[0]); 
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(2.0) , chamb_labels[1]); 
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(3.0) , chamb_labels[2]); 
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(4.0) , chamb_labels[3]); 
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(5.0) , chamb_labels[4]); 
  // grSigma[0]->GetXaxis()->SetBinLabel( grSigma[0]->GetXaxis()->FindBin(6.0) , chamb_labels[5]); 
  // grSigma[0]->GetXaxis()->SetTitleOffset(2.0);
  // grSigma[0]->GetYaxis()->SetTitleOffset(1.2);
  // grSigma[0]->GetXaxis()->SetLabelSize(0.06);
  
  // grSigma[0]->SetMarkerStyle(20); 
  // grSigma[0]->SetMarkerSize(1.2); 
  // grSigma[0]->SetMarkerColor(1);  
  // grSigma[0]->SetLineColor(1);
  // grSigma[0]->Draw("APL");
  // c2->Update();
    
  // grSigma[1]->SetMarkerStyle(21);
  // grSigma[1]->SetMarkerSize(1.2); 
  // grSigma[1]->SetMarkerColor(2);
  // grSigma[1]->SetLineColor(2);
  // grSigma[1]->Draw("PSL");
  // c2->Update();

  // grSigma[2]->SetMarkerStyle(22);
  // grSigma[2]->SetMarkerSize(1.2); 
  // grSigma[2]->SetMarkerColor(3);
  // grSigma[2]->SetLineColor(3);
  // grSigma[2]->Draw("PSL");
  // c2->Update();

  // grSigma[3]->SetMarkerStyle(23);
  // grSigma[3]->SetMarkerSize(1.2); 
  // grSigma[3]->SetMarkerColor(4);
  // grSigma[3]->SetLineColor(4);
  // grSigma[3]->Draw("PSL");
  // c2->Update();
  
  // TLegend *leg = new TLegend(0.8,0.6,0.89,0.89);  //coordinates are fractions of pad dimensions
  // leg->AddEntry(grSigma[0],"50 GeV","LP");  
  // leg->AddEntry(grSigma[1],"100 GeV","LP");  
  // leg->AddEntry(grSigma[2],"150 GeV","LP");  
  // leg->AddEntry(grSigma[3],"200 GeV","LP");  
  // leg->SetFillColor(18);
  // leg->SetHeader(" #mu^{-} energies ");
                                            
  // leg->Draw("PS");
  // c2->Update();
  
  // grSigma[0]-> SetTitle( " Landau Sigma Value for different #mu^{-} energies  " );
  // grSigma[0]->GetHistogram()->SetXTitle(" Chamber ");
  // grSigma[0]->GetHistogram()->SetYTitle(" Landau Sigma (keV) ");
  // grSigma[0]->GetYaxis()->SetRangeUser(0.20,0.40);
  // //grSigma[0]->GetXaxis()->SetRangeUser(0.,5000.);
  // c2->Update();

  // //c2->Print("sigma.png",".png");

  // c2->Write();
  // f.Close();


  
  


}

