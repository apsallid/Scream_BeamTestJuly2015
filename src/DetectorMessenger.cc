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
/// \file analysis/shared/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
// $Id: DetectorMessenger.cc 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger( DetectorConstruction* Det)
: G4UImessenger(),
 fDetector(Det),
 fN03Dir(0),
 fDetDir(0),
 fAbs1MaterCmd(0),
 fAbs2MaterCmd(0),
 fAr1MaterCmd(0),
 fAr2MaterCmd(0),
 fAl1MaterCmd(0),
 fAl2MaterCmd(0),
 fPCB1MaterCmd(0),
 fPCB2MaterCmd(0),
 fAbs1ThickCmd(0),
 fAbs2ThickCmd(0),
 fAr1ThickCmd(0),
 fAr2ThickCmd(0),
 fAl1ThickCmd(0),
 fAl2ThickCmd(0),
 fPCB1ThickCmd(0),
 fPCB2ThickCmd(0),
 fSizeYZCmd(0),
 fEnvThickCmd(0)    
{ 
  fN03Dir = new G4UIdirectory("/N03/");
  fN03Dir->SetGuidance("UI commands of this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/N03/det/",broadcast);
  fDetDir->SetGuidance("detector control");
       
  fAbs1MaterCmd = new G4UIcmdWithAString("/N03/det/setAbs1Mat",this);
  fAbs1MaterCmd->SetGuidance("Select Material of the Absorber1.");
  fAbs1MaterCmd->SetParameterName("choice",false);
  fAbs1MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAbs2MaterCmd = new G4UIcmdWithAString("/N03/det/setAbs2Mat",this);
  fAbs2MaterCmd->SetGuidance("Select Material of the Absorber2.");
  fAbs2MaterCmd->SetParameterName("choice",false);
  fAbs2MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAr1MaterCmd = new G4UIcmdWithAString("/N03/det/setAr1Mat",this);
  fAr1MaterCmd->SetGuidance("Select Material of the Ar1.");
  fAr1MaterCmd->SetParameterName("choice",false);
  fAr1MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fAr2MaterCmd = new G4UIcmdWithAString("/N03/det/setAr2Mat",this);
  fAr2MaterCmd->SetGuidance("Select Material of the Ar2.");
  fAr2MaterCmd->SetParameterName("choice",false);
  fAr2MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl1MaterCmd = new G4UIcmdWithAString("/N03/det/setAl1Mat",this);
  fAl1MaterCmd->SetGuidance("Select Material of the Al1.");
  fAl1MaterCmd->SetParameterName("choice",false);
  fAl1MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fAl2MaterCmd = new G4UIcmdWithAString("/N03/det/setAl2Mat",this);
  fAl2MaterCmd->SetGuidance("Select Material of the Al2.");
  fAl2MaterCmd->SetParameterName("choice",false);
  fAl2MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPCB1MaterCmd = new G4UIcmdWithAString("/N03/det/setPCB1Mat",this);
  fPCB1MaterCmd->SetGuidance("Select Material of the PCB1.");
  fPCB1MaterCmd->SetParameterName("choice",false);
  fPCB1MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fPCB2MaterCmd = new G4UIcmdWithAString("/N03/det/setPCB2Mat",this);
  fPCB2MaterCmd->SetGuidance("Select Material of the PCB2.");
  fPCB2MaterCmd->SetParameterName("choice",false);
  fPCB2MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbs1ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAbs1Thick",this);
  fAbs1ThickCmd->SetGuidance("Set Thickness of the Absorber1");
  fAbs1ThickCmd->SetParameterName("Size",false);
  fAbs1ThickCmd->SetRange("Size>=0.");
  fAbs1ThickCmd->SetUnitCategory("Length");
  fAbs1ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAbs2ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAbs2Thick",this);
  fAbs2ThickCmd->SetGuidance("Set Thickness of the Absorber2");
  fAbs2ThickCmd->SetParameterName("Size",false);
  fAbs2ThickCmd->SetRange("Size>=0.");
  fAbs2ThickCmd->SetUnitCategory("Length");
  fAbs2ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAr1ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAr1Thick",this);
  fAr1ThickCmd->SetGuidance("Set Thickness of the Ar1");
  fAr1ThickCmd->SetParameterName("Size",false);
  fAr1ThickCmd->SetRange("Size>=0.");
  fAr1ThickCmd->SetUnitCategory("Length");  
  fAr1ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAr2ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAr2Thick",this);
  fAr2ThickCmd->SetGuidance("Set Thickness of the Ar2");
  fAr2ThickCmd->SetParameterName("Size",false);
  fAr2ThickCmd->SetRange("Size>=0.");
  fAr2ThickCmd->SetUnitCategory("Length");  
  fAr2ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl1ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAl1Thick",this);
  fAl1ThickCmd->SetGuidance("Set Thickness of the Al1");
  fAl1ThickCmd->SetParameterName("Size",false);
  fAl1ThickCmd->SetRange("Size>=0.");
  fAl1ThickCmd->SetUnitCategory("Length");  
  fAl1ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl2ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAl2Thick",this);
  fAl2ThickCmd->SetGuidance("Set Thickness of the Al2");
  fAl2ThickCmd->SetParameterName("Size",false);
  fAl2ThickCmd->SetRange("Size>=0.");
  fAl2ThickCmd->SetUnitCategory("Length");  
  fAl2ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPCB1ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setPCB1Thick",this);
  fPCB1ThickCmd->SetGuidance("Set Thickness of the PCB1");
  fPCB1ThickCmd->SetParameterName("Size",false);
  fPCB1ThickCmd->SetRange("Size>=0.");
  fPCB1ThickCmd->SetUnitCategory("Length");  
  fPCB1ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPCB2ThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setPCB2Thick",this);
  fPCB2ThickCmd->SetGuidance("Set Thickness of the PCB2");
  fPCB2ThickCmd->SetParameterName("Size",false);
  fPCB2ThickCmd->SetRange("Size>=0.");
  fPCB2ThickCmd->SetUnitCategory("Length");  
  fPCB2ThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("Size",false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");    
  fSizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fEnvThickCmd = new G4UIcmdWithAnInteger("/N03/det/setEnvThickness",this);
  fEnvThickCmd->SetGuidance("Set EnvThickness");
  fEnvThickCmd->SetParameterName("EnvThick",false);
  fEnvThickCmd->SetRange("EnvThick>0 && EnvThick<500");
  fEnvThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fEnvThickCmd;
  delete fAbs1MaterCmd;delete fAbs2MaterCmd; delete fAr1MaterCmd; delete fAr2MaterCmd; delete fAl1MaterCmd; delete fAl2MaterCmd;delete fPCB1MaterCmd; delete fPCB2MaterCmd;
  delete fAbs1ThickCmd;delete fAbs2ThickCmd; delete fAr1ThickCmd; delete fAr2ThickCmd; delete fAl1ThickCmd; delete fAl2ThickCmd;delete fPCB1ThickCmd; delete fPCB2ThickCmd;
  delete fSizeYZCmd;   
  delete fDetDir;
  delete fN03Dir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fAbs1MaterCmd )
   { fDetector->SetAbsorber1Material(newValue);}
   
  if( command == fAbs2MaterCmd )
   { fDetector->SetAbsorber2Material(newValue);}
   
  if( command == fAr1MaterCmd )
   { fDetector->SetAr1Material(newValue);}
  
  if( command == fAr2MaterCmd )
   { fDetector->SetAr2Material(newValue);}

  if( command == fAl1MaterCmd )
   { fDetector->SetAl1Material(newValue);}
  
  if( command == fAl2MaterCmd )
   { fDetector->SetAl2Material(newValue);}

  if( command == fPCB1MaterCmd )
   { fDetector->SetPCB1Material();}
  
  if( command == fPCB2MaterCmd )
   { fDetector->SetPCB2Material();}

  if( command == fAbs1ThickCmd )
   { fDetector->SetAbsorber1Thickness(fAbs1ThickCmd
                                               ->GetNewDoubleValue(newValue));}
   
  if( command == fAbs2ThickCmd )
   { fDetector->SetAbsorber2Thickness(fAbs2ThickCmd
                                               ->GetNewDoubleValue(newValue));}
 
  if( command == fAr1ThickCmd )
   { fDetector->SetAr1Thickness(fAr1ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fAr2ThickCmd )
   { fDetector->SetAr2Thickness(fAr2ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fAl1ThickCmd )
   { fDetector->SetAl1Thickness(fAl1ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fAl2ThickCmd )
   { fDetector->SetAl2Thickness(fAl2ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fPCB1ThickCmd )
   { fDetector->SetPCB1Thickness(fPCB1ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fPCB2ThickCmd )
   { fDetector->SetPCB2Thickness(fPCB2ThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == fSizeYZCmd )
   { fDetector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == fEnvThickCmd )
   { fDetector->SetEnvThickness(fEnvThickCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
