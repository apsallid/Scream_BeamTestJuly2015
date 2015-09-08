//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
//#include "OutfileManager.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
//#include <fstream>
#include "G4Event.hh"
#include "TRandom3.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



EventAction::EventAction(RunAction* run, HistoManager* histo)//, OutfileManager* out)
 :G4UserEventAction(),
  fRunAct(run),fHistoManager(histo),//,fOutfileManager(out),
  fEnergyAbs1(0.), fEnergyAbs2(0.), fEnergyAl1(0.), fEnergyAl2(0.), fEnergyPCB1(0.), fEnergyPCB2(0.),
  //, fEnergyGap1(0.), //fEnergyGap(0.),
  // fTrackLAbs(0.), fTrackLGap(0.),
  fPrintModulo(0)                             
{
  fPrintModulo = 100; 
  outF_=TFile::Open("ScreamTree.root","RECREATE");
  outF_->cd();
  tree_=new TTree("SCREAMTree","SCREAM simulation tree");
  tree_->Branch("fEnergyAbs1",&fEnergyAbs1,"fEnergyAbs1/D");
  tree_->Branch("fEnergyAbs2",&fEnergyAbs2,"fEnergyAbs2/D");
  tree_->Branch("eneAr1",&eneAr1,"eneAr1/D");
  tree_->Branch("eneAr2",&eneAr2,"eneAr2/D");
  tree_->Branch("eneinallcha",&eneinallcha,"eneinallcha/D");
  tree_->Branch("eneAr1sm",&eneAr1,"eneAr1sm/D");
  tree_->Branch("eneAr2sm",&eneAr2,"eneAr2sm/D");
  tree_->Branch("eneinallchasm",&eneinallcha,"eneinallchasm/D");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ 
  outF_->cd();
  tree_->Write();
  outF_->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%fPrintModulo == 0) 
   G4cout << "\n---> Begin of event   : " << evtNb << G4endl;
 counter=evtNb+1;
 // initialisation per event
 eneAr1 = 0.; eneAr2 = 0.; eneinallcha = 0.;
 eneAr1sm = 0.; eneAr2sm = 0.; eneinallchasm = 0.;
 //if (evtNb==1){

  
 for(G4int j=0;j<10;j++){
   for(G4int i=0;i<10;i++){

     fEnergyAbs1 = 0.; fEnergyAbs2 = 0.; fEnergyAl1 = 0.; fEnergyAl2 = 0.; fEnergyPCB1 = 0.; fEnergyPCB2 = 0.;
     fEnergyAr1[i][j]=0.; fEnergyAr2[i][j]=0.; AdcValsAr1[i][j]=0; AdcValsAr2[i][j]=0;

   }
 }


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//std::ofstream outfile("/home/Kalamaris/geant4/Test-build");

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  G4int evtNb = evt->GetEventID();
  G4int  Mapping[96]={89,91,93,95,49,51,53,55,
                      70,72,74,76,78,80,81,83,85,87,
                      77,84,79,82,42,44,46,48,66,68,
                      67,94,69,92,71,90,73,88,75,86,
                      41,56,43,54,45,52,47,50,65,96,
                      18,15,64,33,62,35,60,37,58,39,
                      28,5,26,07,24,9,22,11,20,13,
                      14,16,34,36,38,40,32,1,30,3,
                      25,27,29,31,2,4,6,8,10,12,
                      57,59,61,63,17,19,21,23};
  

  fRunAct->fillPerEvent(fEnergyAbs1, fEnergyAbs2, fEnergyAl1, fEnergyAl2, fEnergyPCB1, fEnergyPCB2, fEnergyAr1[0][0], fEnergyAr2[0][0]);

  G4int padcounter=0;

  std::ofstream outfile("Run.data",std::ios::out);

  //Smearing
  TRandom3 lrndm;
  //If the seed is zero the seed is set to a random value
  lrndm.SetSeed(0);
  
  //Loop over 2 Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber1
      if ( fEnergyAr1[i][j] != 0 ){
	//std::cout<< lrndm.Gaus( 1., (0.237/sqrt(fEnergyGap[k][i][j] * 1000.) ) + 0.08  ) << std::endl;
	eneAr1 = eneAr1 + fEnergyAr1[i][j]; 
	eneAr1sm = eneAr1sm + (fEnergyAr1[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[i][j] * 1000.) ) + 0.08  ) );
      }
      // Chamber2
      if ( fEnergyAr2[i][j] != 0 ){
	  eneAr2 = eneAr2 + fEnergyAr2[i][j]; 
	  eneAr2sm = eneAr2sm + (fEnergyAr2[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr2[i][j] * 1000.) ) + 0.08   ) );
      }


      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);
      if((j==0) && (i==0) )
	{outfile<<"START_RUN"<<G4endl;
	  outfile<<evtNb<<G4endl;
	  outfile<<evtNb<<G4endl;
	  outfile<<"START_EVENT"<<G4endl;
	  outfile<<evtNb+1<<G4endl;
	  outfile<<evtNb+1<<G4endl;
	}

      if ((i==0) && (j==0))  outfile<<96<<G4endl;

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr1[i][j]=G4int( (fEnergyAr1[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 
      AdcValsAr2[i][j]=G4int( (fEnergyAr2[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      if (AdcValsAr1[i][j]>1023) {ovf=true;AdcValsAr1[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr1[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }

      if (AdcValsAr2[i][j]>1023) {ovf=true;AdcValsAr2[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr2[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
  
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
  //fEnergyGap[k][i][j]
  fHistoManager->FillHisto(1, fEnergyAr1[0][4]);//First Detector pad (1,5) 
  fHistoManager->FillHisto(2, fEnergyAr1[1][4]);//First Detector pad (2,5)
  fHistoManager->FillHisto(3, fEnergyAr1[2][4]);//First Detector pad (3,5)
  fHistoManager->FillHisto(4, fEnergyAr1[3][4]);//First Detector pad (4,5)
  fHistoManager->FillHisto(5, fEnergyAr1[4][4]);//First Detector pad (5,5)
  fHistoManager->FillHisto(6, eneAr1);//Energy in chamber1
  fHistoManager->FillHisto(7, eneAr2);//Energy in chamber2

  eneinallcha = eneAr1 + eneAr2;//Energy in all chamber
  fHistoManager->FillHisto(8, eneinallcha);//Energy in all chambers

  //Smeared distributions
  if ( fEnergyAr1[0][4] != 0 ){
    fHistoManager->FillHisto(9, fEnergyAr1[0][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[0][4] * 1000.) )  + 0.08  )  );//First Detector pad (1,5) 
  }
  if ( fEnergyAr1[1][4] != 0 ){
    fHistoManager->FillHisto(10, fEnergyAr1[1][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[1][4] * 1000.) )  + 0.08  )  );//First Detector pad (2,5)
  }
  if ( fEnergyAr1[2][4] != 0 ){
    fHistoManager->FillHisto(11, fEnergyAr1[2][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[2][4] * 1000.) )  + 0.08 )  );//First Detector pad (3,5)
  }
  if ( fEnergyAr1[3][4] != 0 ){
    fHistoManager->FillHisto(12, fEnergyAr1[3][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[3][4] * 1000.) )  + 0.08 )  );//First Detector pad (4,5)
  }
  if ( fEnergyAr1[4][4] != 0 ){
    fHistoManager->FillHisto(13, fEnergyAr1[4][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[4][4] * 1000.) )  + 0.08 )  );//First Detector pad (5,5)
  }
  fHistoManager->FillHisto(14, eneAr1sm );//Smear Energy in chamber1
  fHistoManager->FillHisto(15, eneAr2sm );//Smear Energy in chamber2

  eneinallchasm = eneAr1sm + eneAr2sm;//Energy in all chamber

  fHistoManager->FillHisto(16, eneinallchasm);//Energy in all chamber
  if ( fEnergyAr1[0][4] != 0 ){
    fHistoManager->FillHisto(17, lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[0][4] * 1000.) ) + 0.08  ) );//Testing the smear factor only for the 1st chamber pad. This should be gaussian.
  }
  fHistoManager->FillHisto(18, eneAr1sm * 30.95062);//Energy in chamber1 in adc
  fHistoManager->FillHisto(19, eneAr2sm * 30.95062);//Energy in chamber2 in adc

  fHistoManager->FillHisto(20, eneinallcha * 30.95062);//Energy in all chamber in adc

  tree_->Fill();
  //G4cout<< "First Detector pad (1,5) " << fEnergyAr1[0][4] <<G4endl;
  //G4cout<< "First Detector pad (2,5) " << fEnergyAr1[1][4] <<G4endl;

  //G4cout<<"Padcounter  "<<padcounter<<G4endl;
  //G4cout<<"End of event # "<<counter<<G4endl;
}

  /* for(G4int j=0;j<10;j++){
  for(G4int i=0;i<10;i++){

  if (AdcVals[i][j]=1023) {ovfcounter++;}

}


}
}
  G4cout<<"Overflow # per event "<<ovfcounter<<G4endl;
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
