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
void AbsorberStudy_RooFit(){

  //===============================================================================
  gStyle->SetOptTitle(0);
  
  using namespace RooFit;
  
  //The root files with the individual plots 
  const int numberoffiles = 20;
  const int numberofabsorbers = 1;

  TString filename[numberoffiles];
  //e-
  // filename[0]="Scream_50cmFeAbsorber_ele10GeV.root";
  // filename[1]="Scream_50cmFeAbsorber_ele20GeV.root";
  // filename[2]="Scream_50cmFeAbsorber_ele30GeV.root";
  // filename[3]="Scream_50cmFeAbsorber_ele40GeV.root";
  // filename[4]="Scream_50cmFeAbsorber_ele50GeV.root";
  // filename[5]="Scream_50cmFeAbsorber_ele60GeV.root";
  // filename[6]="Scream_50cmFeAbsorber_ele70GeV.root";
  // filename[7]="Scream_50cmFeAbsorber_ele80GeV.root";
  // filename[8]="Scream_50cmFeAbsorber_ele90GeV.root";
  // filename[9]="Scream_50cmFeAbsorber_ele100GeV.root";
  // filename[10]="Scream_50cmFeAbsorber_ele110GeV.root";
  // filename[11]="Scream_50cmFeAbsorber_ele120GeV.root";
  // filename[12]="Scream_50cmFeAbsorber_ele130GeV.root";
  // filename[13]="Scream_50cmFeAbsorber_ele140GeV.root";
  // filename[14]="Scream_50cmFeAbsorber_ele150GeV.root";
  // filename[15]="Scream_50cmFeAbsorber_ele160GeV.root";
  // filename[16]="Scream_50cmFeAbsorber_ele170GeV.root";
  // filename[17]="Scream_50cmFeAbsorber_ele180GeV.root";
  // filename[18]="Scream_50cmFeAbsorber_ele190GeV.root";
  // filename[19]="Scream_50cmFeAbsorber_ele200GeV.root";

  filename[0]="testcalo_e-_10.root";  
  filename[1]="testcalo_e-_20.root";  
  filename[2]="testcalo_e-_30.root";  
  filename[3]="testcalo_e-_40.root";  
  filename[4]="testcalo_e-_50.root";  
  filename[5]="testcalo_e-_60.root";  
  filename[6]="testcalo_e-_70.root";  
  filename[7]="testcalo_e-_80.root";  
  filename[8]="testcalo_e-_90.root";  
  filename[9]="testcalo_e-_100.root"; 
  filename[10]="testcalo_e-_110.root"; 
  filename[11]="testcalo_e-_120.root"; 
  filename[12]="testcalo_e-_130.root"; 
  filename[13]="testcalo_e-_140.root"; 
  filename[14]="testcalo_e-_150.root"; 
  filename[15]="testcalo_e-_160.root"; 
  filename[16]="testcalo_e-_170.root"; 
  filename[17]="testcalo_e-_180.root"; 
  filename[18]="testcalo_e-_190.root"; 
  filename[19]="testcalo_e-_200.root"; 





  //Energy distributions we want to fit 
  TString energy[numberofabsorbers];
  energy[0] = "histo/33";
 
  //For the maxv and sigma plot 
  std::vector<double> ene;
  std::vector<double> maxv;
  std::vector<double> thesigma;
  std::vector<double> eneerr;
  std::vector<double> maxvr;
  std::vector<double> thesigmaerr;
  ene.push_back(10.);ene.push_back(20.);ene.push_back(30.);ene.push_back(40.);ene.push_back(50.);ene.push_back(60.);ene.push_back(70.);ene.push_back(80.);
  ene.push_back(90.);ene.push_back(100.);ene.push_back(110.);ene.push_back(120.);ene.push_back(130.);ene.push_back(140.);ene.push_back(150.);ene.push_back(160.);
  ene.push_back(170.);ene.push_back(180.);ene.push_back(190.);ene.push_back(200.);
  
  //Useful parameters for the fitting
  double I = 0.;
  double histo_mean = 0.;
  double RMS = 0.;
  int   division = 0;
  float eneMIN = 0.;
  float eneMAX = 0.;
  float BIN_SIZE = 0.;
  float histomaximum = 0.;

  //======================================================================================
  //Canvas
  TCanvas *can[numberoffiles][numberofabsorbers];
  TString canname;

  TCanvas *canmaxv[numberoffiles];
  
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
      can[i][k] = new TCanvas(canname, canname,800,600);
      can[i][k]->cd();

      histo[i][k]->Draw();
      // histo[i][k]->Scale(0.1);//This is for the 10 events I run.
      //Set the initial values for the fitting
      I = histo[i][k]->Integral();
      histo_mean = histo[i][k]->GetMean();
      RMS = histo[i][k]->GetRMS();
      division = histo[i][k]->GetNbinsX();
      eneMIN = histo[i][k]->GetBinLowEdge(1);
      eneMAX = histo[i][k]->GetBinLowEdge(division+1);
      BIN_SIZE = histo[i][k]->GetBinWidth(1);
      histomaximum = histo[i][k]->GetBinContent(histo[i][k]->GetMaximumBin());
      
      // Declare observable x
      RooRealVar x("x","x",eneMIN,eneMAX) ;
      eneMIN = histo[i][k]->GetMaximumBin() - RMS; eneMAX = histo[i][k]->GetMaximumBin() + RMS;
      // std::cout << histo[i][k]->GetMaximumBin() << std::endl;
      // eneMIN = 50; eneMAX = 200.;
      x.setRange("fitRegion1",eneMIN  ,  eneMAX );

      RooDataHist dh("dh","dh",x,Import(*histo[i][k])) ;
      
      RooPlot* frame = x.frame(Title(filename[i]), Bins(1000)) ;
      dh.plotOn(frame,MarkerColor(2),MarkerSize(0.9),MarkerStyle(21));  //this will show histogram data points on canvas 
      dh.statOn(frame);  //this will display hist stat on canvas
      
      RooRealVar mean("mean","mean",histo_mean, eneMIN, eneMAX);
      RooRealVar width("width","width",RMS, eneMIN, eneMAX);
      RooRealVar sigma("sigma","sigma",RMS, eneMIN, eneMAX);
      
       RooGaussian gauss("gauss","gauss",x,mean,sigma);
      // RooExponential expo("expo", "exponential PDF", x, lambda);
      
      // Construct landau(t,ml,sl) ;
      //RooLandau landau("lx","lx",x,mean,sigma) ;
  
      //RooBreitWigner BW1("BW1","Breit Wigner theory",x,mean,sigma);
      RooVoigtian bwandgaus("bwandgaus", "Breit Wigner convoluted with gauss",x,mean,width,sigma);
      //RooCBShape cryBall1("cryBall1","Crystal Ball resolution model", x, CBmean, CBsigma, CBalpha, CBn) ;
      //RooFFTConvPdf bwxCryBall1("bwxCryBall1", "FFT Conv CryBall and BW", x, BW1, cryBall1); 
      // RooBreitWigner gauss("gauss","gauss",x,mean,sigma);
      // RooVoigtian gauss("gauss","gauss",x,mean,width,sigma);

      // RooFitResult* filters = gauss.fitTo(dh, NumCPU(4), Range(("fitRegion"+IntToString(i)).c_str()), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
       RooFitResult* filters = gauss.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      //RooFitResult* filters = bwxCryBall1.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      // RooFitResult* filters = BW1.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      //RooFitResult* filters = bwandgaus.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      // RooFitResult* filters = sum.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      // RooFitResult* filters = gauss.fitTo(dh, NumCPU(4), Save(true), SumW2Error(kTRUE), PrintLevel(-1));

      //RooFitResult* filters = landau.fitTo(dh, NumCPU(4), Range("fitRegion1"), Save(true), SumW2Error(kTRUE), PrintLevel(-1));
      gauss.plotOn(frame,LineColor(4));//this will show fit overlay on canvas 
      gauss.paramOn(frame); //this will display the fit parameters on canvas
      
      filters->Print();
      
      // // Access list of final fit parameter values
      // std::cout << "final value of floating parameters" << std::endl ;
      // filters->floatParsFinal().Print("s") ;
      
      //mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 
      mpc = mean.getVal();

      maxv.push_back(mpc);
      maxvr.push_back( (mean.getError()) );
      eneerr.push_back(0.);

      thesigma.push_back( (sigma.getVal()) );
      thesigmaerr.push_back( (sigma.getError()) );
   
      // std::cout << "File " << filename[i] << " Eneer " << k << " MAXV "<< mpc << std::endl;

      // can[i][k]->Print(canname + ".root",".root");
      // can[i][k]->Close();
  
  	// // Draw all frames on a canvas
	// can[k] = new TCanvas(filename[i], filename[i],800,600);
	// can[k]->cd();
	// gPad->SetLeftMargin(0.15);
	// frame->GetXaxis()->SetTitle("E (GeV)");  
	// frame->GetXaxis()->SetTitleOffset(1.2);
	// frame->GetXaxis()->SetRangeUser(0.,0.05);
	// // float binsize = histo[i][k]->GetBinWidth(1); 
	// // char Bsize[50]; 
	// //sprintf(Bsize,"Events per %2.2f",binsize);
	// frame->GetYaxis()->SetTitle("Events");  
	//frame->GetYaxis()->SetTitleOffset(1.2);
	frame->Draw() ;
	// can[k]->Print(filename[i]+".pdf",".pdf");
	// can[k]->Print(histogramname[k]+".png",".png");
	// can[k]->Close();

	can[i][k]->Print(canname + ".root",".root");
	//can[i][k]->Close();


      

      
    } // end of loop over histos  
    //For the multi plot
 
  } // end of loop over files

   TFile f("maxvandsigma.root","recreate");

  TGraphErrors *grMAXV = new TGraphErrors(maxv.size(),&ene[0],&maxv[0],&maxvr[0],&eneerr[0]);
  TGraphErrors *grSigma = new TGraphErrors(thesigma.size(),&ene[0],&thesigma[0],&thesigmaerr[0],&eneerr[0]);
     
  //======================================================================================
  //For the maxv plot
  TCanvas *c1 = new TCanvas("c1","maxv",200,10,1200,968); 
  grMAXV->SetMarkerStyle(20); 
  grMAXV->SetMarkerSize(1.2); 
  grMAXV->Draw("APL");
  c1->Update();
  grMAXV->SetTitle( " Maximum energy loss for different e^{-} energies in Fe absorber " );
  grMAXV->GetHistogram()->SetXTitle(" Energy (GeV) ");
  grMAXV->GetHistogram()->SetYTitle(" Maximum energy loss location (mm) ");
  grMAXV->GetYaxis()->SetRangeUser(50.,150.);
  c1->Update();

  //c1->Print("maxv.png",".png");

  c1->Write();
  c1->Close();


  //======================================================================================
  //For the sigma plot
  TCanvas *c2 = new TCanvas("c2","sigma",200,10,1200,968); 
  grSigma->SetMarkerStyle(20); 
  grSigma->SetMarkerSize(1.2); 
  grSigma->Draw("APL");
  c2->Update();

  grSigma->SetTitle( " Sigma Value for different e^{-} energies  " );
  grSigma->GetHistogram()->SetXTitle(" Energy (GeV) ");
  grSigma->GetHistogram()->SetYTitle(" Sigma (mm) ");
  //grSigma->GetYaxis()->SetRangeUser(0.20,0.40);
  c2->Update();
  //c2->Print("sigma.png",".png");
  c2->Write();
  c2->Close();
  f.Close();


}

