{
  gROOT->Reset();
  //gROOT->ForceStyle(); 
  gStyle->SetOptStat("ksiourmen");

  // Draw histos filled by Geant4 simulation 
  //   
  // TFile f = TFile("testcalo_e-_150_absorberinmaxofshower.root");  
  //TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_allenergy.root");  
  // TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_4mmradius.root");  
  // TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_20mmradius.root");  
  // TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_allenergy_nogaptodetector.root");  
  // TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_25mmradius_nogaptodetector.root"); 
  TFile f = TFile("testcalo_e-_100_10cmFeAbsorber_30mmradius_nogaptodetector.root"); 
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  //TH1D* hist1 = (TH1D*)f.Get("1");
  //hist1->Draw("HIST");
  
  //TH1D* hist2 = (TH1D*)f.Get("2");
  //hist2->Draw("HIST");
  
  TH1D* hist33 = (TH1D*)f.Get("histo/6");
  hist33->GetXaxis()->SetTitle("Energy in chamber 1 (MeV)"); 
  //hist33->GetYaxis()->SetTitle("Energy deposited (MeV/mm)");
  hist33->SetTitle( "100 GeV e^{-} with 100.0 mm Fe absorber and 3mm Argon layer in 30mm radius  " ); //
  hist33->Rebin(32);

  //c1->SetLogy(1);
  c1->cd();
  c1->Update(); 
  hist33->Draw("HIST");
  //c1->Print("e-_50GeV_absorberinmaxofshower.pdf",".pdf");
}  
