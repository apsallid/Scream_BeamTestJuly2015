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
/// \file analysis/shared/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <fstream>

#include "TTree.h"
#include "TFile.h"

class RunAction;
class HistoManager;
//class OutfileManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int chambers =2;

class EventAction : public G4UserEventAction
{
public:
  //G4double fEnergyGap[10][10];
  EventAction(RunAction*, HistoManager*);
  virtual ~EventAction();
 
  virtual void  BeginOfEventAction(const G4Event*);
  virtual void    EndOfEventAction(const G4Event*);
  

  
  void AddAbs1(G4double de) {fEnergyAbs1 += de;};
  void AddAbs2(G4double de) {fEnergyAbs2 += de;};
  void AddAl1(G4double de) {fEnergyAl1 += de;};
  void AddAl2(G4double de) {fEnergyAl2 += de;};
  void AddAr1(G4double de, int i , int j)  {fEnergyAr1[i][j] += de;};
  void AddAr2(G4double de, int i , int j)  {fEnergyAr2[i][j] += de;};
  void AddPCB1(G4double de) {fEnergyPCB1 += de;};
  void AddPCB2(G4double de) {fEnergyPCB2 += de;};



  G4int ovfcounter;
  G4int counter;
  G4bool ovf;
  G4int  AdcValsAr1[10][10], AdcValsAr2[10][10];

 
    
private:
   RunAction*    fRunAct;
   HistoManager* fHistoManager;
 
   
  G4double  fEnergyAbs1, fEnergyAbs2, fEnergyAl1, fEnergyAl2, fEnergyAr1[10][10], fEnergyAr2[10][10], fEnergyPCB1, fEnergyPCB2;

  //Energies
  G4double eneAr1, eneAr2, eneinallcha; 
  G4double eneAr1sm, eneAr2sm, eneinallchasm; 
  
  TFile *outF_;
  TTree *tree_;

  
                    
   G4int     fPrintModulo;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
