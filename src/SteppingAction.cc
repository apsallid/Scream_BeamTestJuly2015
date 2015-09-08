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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
//#include "OutfileManager.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
			       EventAction* evt)
: G4UserSteppingAction(), 
  fDetector(det), fEventAction(evt)                                   
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
 // G4double Earray[10][10]; 
  /*G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();*/

  //Bragg curve
  // G4StepPoint* prePoint  = aStep->GetPreStepPoint();
  // G4StepPoint* postPoint = aStep->GetPostStepPoint();
   
  // G4double x1 = prePoint->GetPosition().x();
  // G4double x2 = postPoint->GetPosition().x();  
  // G4double x  = x1 + fabs(G4UniformRand()*(x2-x1)) + 0.5*(fDetector->GetAbsorber1Thickness());

  // std::cout << "The absorber size is " << fDetector->GetAbsorberThickness() << std::endl;

  // G4double y1 = prePoint->GetPosition().y();
  // G4double y2 = postPoint->GetPosition().y();  
  // G4double y  = y1 + fabs(G4UniformRand()*(y2-y1));

  // G4double z1 = prePoint->GetPosition().z();
  // G4double z2 = postPoint->GetPosition().z();  
  // G4double z  = z1 + fabs(G4UniformRand()*(z2-z1));
  
  // std::cout << "y position is  " <<  y << std::endl;
  // std::cout << "z position is  " <<  z << std::endl;
  // std::cout << "The radius is  " <<  sqrt( pow(y,2.) + pow(z,2.) ) << std::endl;

  //HistoManager* histo;
  //fHistoManager->FillHisto(33, x, edep);//Energy loss 
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        
  // if (volume == fDetector->GetAbsorber()) {
  //   analysisManager->FillH1(33, x, edep);
  //   fEventAction->AddAbs(edep);
  // }

  if (volume == fDetector->GetAbsorber1()) {fEventAction->AddAbs1(edep);}

  if (volume == fDetector->GetAbsorber2()) {fEventAction->AddAbs2(edep);}

  if (volume == fDetector->GetAl1()) {fEventAction->AddAl1(edep);}

  if (volume == fDetector->GetAl2()) {fEventAction->AddAl2(edep);}

  if (volume == fDetector->GetPCB1()) {fEventAction->AddPCB1(edep);}

  if (volume == fDetector->GetPCB2()) {fEventAction->AddPCB2(edep);}

  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      //Disk of radius 2 mm
      //bool radius = ( sqrt( pow(y,2.) + pow(z,2.) ) < 2. ) ; // in mm

      // if ((volume == fDetector->GetGap(i,j))  && (k==0) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap1(i,j)) && (k==1) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap2(i,j)) && (k==2) && radius ) fEventAction->AddGap(edep,k,i,j);
 
      // if ((volume == fDetector->GetGap3(i,j)) && (k==3) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap4(i,j)) && (k==4) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap5(i,j)) && (k==5) && radius ) fEventAction->AddGap(edep,k,i,j);
 
      if ( (volume == fDetector->GetAr1(i,j)) ) fEventAction->AddAr1(edep,i,j);

      if ( (volume == fDetector->GetAr2(i,j)) ) fEventAction->AddAr2(edep,i,j);

    }//for i
  }//for j

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
