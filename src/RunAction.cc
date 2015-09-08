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
/// \file analysis/shared/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// $Id: RunAction.cc 74272 2013-10-02 14:48:50Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"
#include "EventAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include <fstream>
//#include "OutfileManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//std::ofstream outfile("/home/Kalamaris/geant4/Test-build");

RunAction::RunAction(HistoManager* histo)//, OutfileManager* out)
: G4UserRunAction(),
  fHistoManager(histo),//fOutfileManager(out),
  fSumEAbs1(0.),fSumEAbs2(0.), //fSum2EAbs(0.),
  fSumEAl1(0.), fSumEAl2(0.),
  fSumEPCB1(0.), fSumEPCB2(0.),
  fSumEAr1(0.), fSumEAr2(0.) //fSum2EGap(0.),
 // fSumEGap1(0.)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    
  //initialize cumulative quantities
  //
  fSumEAbs1 = fSumEAbs2 = fSumEAl1  = fSumEAl2  = fSumEPCB1  = fSumEPCB2 = fSumEAr1  = fSumEAr2  =0; // fSum2EGap1 = 0.;
  //fSumLAbs = fSum2LAbs =fSumLGap = fSum2LGap=0; //=fSumLGap1 = fSum2LGap1= 0.;

  //histograms
  //
  fHistoManager->book(); 
  //fOutfileManager->Outbook();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::fillPerEvent(G4double EAbs1, G4double EAbs2, G4double EAl1, G4double EAl2, G4double EPCB1, G4double EPCB2, G4double EAr1, G4double EAr2)
{
  //accumulate statistic
  //
  fSumEAbs1  += EAbs1;  //fSum2EAbs1 += EAbs1*EAbs1;
  fSumEAbs2  += EAbs2;  //fSum2EAbs2 += EAbs2*EAbs2;
  fSumEAl1  += EAl1;  //fSum2EAr1 += EAr1*EAr1;
  fSumEAl2  += EAl2;  //fSum2EAr2 += EAr2*EAr2;
  fSumEPCB1  += EPCB1;  //fSum2EAr1 += EAr1*EAr1;
  fSumEPCB2  += EPCB2;  //fSum2EAr2 += EAr2*EAr2;
  fSumEAr1  += EAr1;  //fSum2EAr1 += EAr1*EAr1;
  fSumEAr2  += EAr2;  //fSum2EAr2 += EAr2*EAr2;
  //fSumEAr11 += EAr11;
  
 // fSumLAbs += LAbs;  fSum2LAbs += LAbs*LAbs;
 // fSumLAr1 += LAr1;  fSum2LAr1 += LAr1*LAr1;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
//fSumEAbs /= NbOfEvents; //fSum2EAbs /= NbOfEvents;
  //G4double rmsEAbs = fSum2EAbs - fSumEAbs*fSumEAbs;
  //if (rmsEAbs >0.) rmsEAbs = std::sqrt(rmsEAbs); else rmsEAbs = 0.;
  
//  fSumEGap /= NbOfEvents; //fSum2EGap /= NbOfEvents;
  //G4double rmsEGap = fSum2EGap - fSumEGap*fSumEGap;
  //if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;
  
// fSumEGap1 /=NbOfEvents;
  //fSumLAbs /= NbOfEvents; fSum2LAbs /= NbOfEvents;
  // G4double rmsLAbs = fSum2LAbs - fSumLAbs*fSumLAbs;
  // if (rmsLAbs >0.) rmsLAbs = std::sqrt(rmsLAbs); else rmsLAbs = 0.;
  
  //fSumLGap /= NbOfEvents; fSum2LGap /= NbOfEvents;
  //G4double rmsLGap = fSum2LGap - fSumLGap*fSumLGap;
  //if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  //print
  //


  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n Energy in Absorber1 : " << G4BestUnit(fSumEAbs1,"Energy")
     //<< " +- "                          << G4BestUnit(rmsEAbs1,"Energy")  
     << "\n Energy in Absorber2 : " << G4BestUnit(fSumEAbs2,"Energy")
     << "\n Energy in Chamber1      : " << G4BestUnit(fSumEAr1,"Energy")
     << "\n Energy in Chamber2      : " << G4BestUnit(fSumEAr2,"Energy")
     //<< " +- "                          << G4BestUnit(rmsEGap,"Energy")
  
     << G4endl;
     
  /*G4cout
     << "\n mean trackLength in Absorber : " << G4BestUnit(fSumLAbs,"Length")
     << " +- "                               << G4BestUnit(rmsLAbs,"Length")  
     << "\n mean trackLength in Gap      : " << G4BestUnit(fSumLGap,"Length")
     << " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl; */
     
  //save histograms
  //
  //fHistoManager->PrintStatistic();
  fHistoManager->save();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
