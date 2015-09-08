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
/// \file analysis/shared/include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
//
// $Id: DetectorMessenger.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
    virtual ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DetectorConstruction* fDetector;
    
    G4UIdirectory*             fN03Dir;
    G4UIdirectory*             fDetDir;
    G4UIcmdWithAString*        fAbs1MaterCmd;
    G4UIcmdWithAString*        fAbs2MaterCmd;
    G4UIcmdWithAString*        fAr1MaterCmd;
    G4UIcmdWithAString*        fAr2MaterCmd;
    G4UIcmdWithAString*        fAl1MaterCmd;
    G4UIcmdWithAString*        fAl2MaterCmd;
    G4UIcmdWithAString*        fPCB1MaterCmd;
    G4UIcmdWithAString*        fPCB2MaterCmd;
    G4UIcmdWithADoubleAndUnit* fAbs1ThickCmd;
    G4UIcmdWithADoubleAndUnit* fAbs2ThickCmd;
    G4UIcmdWithADoubleAndUnit* fAr1ThickCmd;
    G4UIcmdWithADoubleAndUnit* fAr2ThickCmd;
    G4UIcmdWithADoubleAndUnit* fAl1ThickCmd;
    G4UIcmdWithADoubleAndUnit* fAl2ThickCmd;
    G4UIcmdWithADoubleAndUnit* fPCB1ThickCmd;
    G4UIcmdWithADoubleAndUnit* fPCB2ThickCmd;
    G4UIcmdWithADoubleAndUnit* fSizeYZCmd;
    G4UIcmdWithAnInteger*      fEnvThickCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

