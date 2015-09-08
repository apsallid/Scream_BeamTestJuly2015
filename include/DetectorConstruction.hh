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
/// \file analysis/shared/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//#include "G4VSensitiveDetector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
     
     void SetAbsorber1Material (G4String);     
     void SetAbsorber1Thickness(G4double);     

     void SetAbsorber2Material (G4String);     
     void SetAbsorber2Thickness(G4double);     

     void SetAl1Material (G4String);     
     void SetAl1Thickness(G4double);
     
     void SetAl2Material (G4String);     
     void SetAl2Thickness(G4double);

     void SetAr1Material (G4String);     
     void SetAr1Thickness(G4double);
     
     void SetAr2Material (G4String);     
     void SetAr2Thickness(G4double);

     void SetPCB1Material ();     
     void SetPCB1Thickness(G4double);
     
     void SetPCB2Material ();     
     void SetPCB2Thickness(G4double);

     void SetCalorSizeYZ(G4double);          
     void SetEnvThickness (G4int);   
     
     virtual G4VPhysicalVolume* Construct();

  public:
  
     void PrintCalorParameters(); 
                    
     G4double GetWorldSizeX()           {return fWorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return fWorldSizeYZ;};
     
     G4double GetCalorThickness()       {return fCalorThickness;}; 
     G4double GetCalorSizeYZ()          {return fCalorSizeYZ;};
      
     G4int GetEnvThickness()              {return fEnvThickness;}; 
     
     G4Material* GetAbsorber1Material()  {return fAbsorber1Material;};
     G4double    GetAbsorber1Thickness() {return fAbsorber1Thickness;};      
     
     G4Material* GetAbsorber2Material()  {return fAbsorber2Material;};
     G4double    GetAbsorber2Thickness() {return fAbsorber2Thickness;};      

     G4Material* GetAl1Material()       {return fAl1Material;};
     G4double    GetAl1Thickness()      {return fAl1Thickness;};

     G4Material* GetAl2Material()       {return fAl2Material;};
     G4double    GetAl2Thickness()      {return fAl2Thickness;};
         
     G4Material* GetAr1Material()       {return fAr1Material;};
     G4double    GetAr1Thickness()      {return fAr1Thickness;};

     G4Material* GetAr2Material()       {return fAr2Material;};
     G4double    GetAr2Thickness()      {return fAr2Thickness;};
     
     G4Material* GetPCB1Material()       {return fPCB1Material;};
     G4double    GetPCB1Thickness()      {return fPCB1Thickness;};

     G4Material* GetPCB2Material()       {return fPCB2Material;};
     G4double    GetPCB2Thickness()      {return fPCB2Thickness;};
 
     G4double    GetAirThickness()      {return fAirThickness;};

     const G4VPhysicalVolume* GetphysiWorld()              {return fPhysiWorld;};           
     const G4VPhysicalVolume* GetAbsorber1()                {return fPhysiAbsorber1;};
     const G4VPhysicalVolume* GetAbsorber2()                {return fPhysiAbsorber2;};
     const G4VPhysicalVolume* GetAl1()                       {return fPhysiAl1;};
     const G4VPhysicalVolume* GetAl2()                       {return fPhysiAl2;};
     const G4VPhysicalVolume* GetAr1(int i, int j)        {return fPhysiAr1[i][j];};
     const G4VPhysicalVolume* GetAr2(int i, int j)        {return fPhysiAr2[i][j];};
     const G4VPhysicalVolume* GetPCB1()                       {return fPhysiPCB1;};
     const G4VPhysicalVolume* GetPCB2()                       {return fPhysiPCB2;};

            
  private:
     
     G4Material*        fAbsorber1Material;
     G4double           fAbsorber1Thickness;
     
     G4Material*        fAbsorber2Material;
     G4double           fAbsorber2Thickness;

     G4Material*        fAl1Material;
     G4double           fAl1Thickness;
     
     G4Material*        fAl2Material;
     G4double           fAl2Thickness;

     G4Material*        fAr1Material;
     G4double           fAr1Thickness;
     
     G4Material*        fAr2Material;
     G4double           fAr2Thickness;

     G4Material*        fPCB1Material;
     G4double           fPCB1Thickness;
     
     G4Material*        fPCB2Material;
     G4double           fPCB2Thickness;

     G4double           fAirThickness;

     G4int              fEnvThickness;
     G4double           fLayerThickness;
          
     G4double           fCalorSizeYZ;
     G4double           fCalorThickness;
     
     G4Material*        fDefaultMaterial;
     G4double           fWorldSizeYZ;
     G4double           fWorldSizeX;
            
     G4Box*             fSolidWorld;    //pointer to the solid World 
     G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

     G4Box*             fSolidCalor;    //pointer to the solid Calor 
     G4LogicalVolume*   fLogicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume* fPhysiCalor;    //pointer to the physical Calor
     
   //  G4Box*             fSolidLayer;    //pointer to the solid Layer 
   //  G4LogicalVolume*   fLogicLayer;    //pointer to the logical Layer
   //  G4VPhysicalVolume* fPhysiLayer;    //pointer to the physical Layer
         
     G4Box*             fSolidAbsorber1; //pointer to the solid Absorber1
     G4LogicalVolume*   fLogicAbsorber1; //pointer to the logical Absorber1
     G4VPhysicalVolume* fPhysiAbsorber1; //pointer to the physical Absorber1
     
     G4Box*             fSolidAbsorber2; //pointer to the solid Absorber2
     G4LogicalVolume*   fLogicAbsorber2; //pointer to the logical Absorber2
     G4VPhysicalVolume* fPhysiAbsorber2; //pointer to the physical Absorber2

     G4Box*             fSolidAl1;      //pointer to the solid Al1
     G4LogicalVolume*   fLogicAl1;      //pointer to the logical Al1
     G4VPhysicalVolume* fPhysiAl1;      //pointer to the physical Al1

     G4Box*             fSolidAl2;      //pointer to the solid Al2
     G4LogicalVolume*   fLogicAl2;      //pointer to the logical Al2
     G4VPhysicalVolume* fPhysiAl2;      //pointer to the physical Al2

     G4Box*             fSolidAr1[10][10];      //pointer to the solid Ar1
     G4LogicalVolume*   fLogicAr1;      //pointer to the logical Ar1
     G4VPhysicalVolume* fPhysiAr1[10][10];      //pointer to the physical Ar1

     G4Box*             fSolidAr2[10][10];      //pointer to the solid Ar2
     G4LogicalVolume*   fLogicAr2;      //pointer to the logical Ar2
     G4VPhysicalVolume* fPhysiAr2[10][10];      //pointer to the physical Ar2

     G4Box*             fSolidPCB1;      //pointer to the solid PCB1
     G4LogicalVolume*   fLogicPCB1;      //pointer to the logical PCB1
     G4VPhysicalVolume* fPhysiPCB1;      //pointer to the physical PCB1
   
     G4Box*             fSolidPCB2;      //pointer to the solid PCB2
     G4LogicalVolume*   fLogicPCB2;      //pointer to the logical PCB2
     G4VPhysicalVolume* fPhysiPCB2;      //pointer to the physical PCB2
     
     DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fCalorThickness = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness;

  fWorldSizeX = 1.2*fCalorThickness;
  fWorldSizeYZ = 1.2*fCalorSizeYZ;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

