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
/// \file analysis/shared/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "G4SDManager.hh"
//#include "G4VSensitiveDetector.hh"
//#include "SDGap.hh"
#include "G4PVReplica.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 :G4VUserDetectorConstruction(),
  fAbsorber1Material(0),
  fAbsorber2Material(0),
  fAl1Material(0),
  fAl2Material(0),
  fAr1Material(0),
  fAr2Material(0),
  fPCB1Material(0),
  fPCB2Material(0),
  fDefaultMaterial(0),
  fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
  fSolidCalor(0),fLogicCalor(0),fPhysiCalor(0),
  fSolidAbsorber1(0),fLogicAbsorber1(0),fPhysiAbsorber1(0),
  fSolidAbsorber2(0),fLogicAbsorber2(0),fPhysiAbsorber2(0),
  fSolidAl1(0),fLogicAl1(0),fPhysiAl1(0),
  fSolidAl2(0),fLogicAl2(0),fPhysiAl2(0),
  fLogicAr1(0),fLogicAr2(0),
  //fSolidAr1(0),fLogicAr1(0),fPhysiAr1(0),
  //fSolidAr2(0),fLogicAr2(0),fPhysiAr2(0),
  fSolidPCB1(0),fLogicPCB1(0),fPhysiPCB1(0),
  fSolidPCB2(0),fLogicPCB2(0),fPhysiPCB2(0),
  fDetectorMessenger(0)
{
  // default parameter values of the calorimeter
  fAbsorber1Thickness = 4.0*cm; 
  fAbsorber2Thickness = 6.0*cm; 
  fAl1Thickness      = 2.*mm;
  fAl2Thickness      = 2.*mm;
  fAr1Thickness      = 3.*mm;
  fAr2Thickness      = 3.*mm;
  fPCB1Thickness     = 2.*mm;
  fPCB2Thickness     = 2.*mm;
  fAirThickness      = 3.*mm;

  fEnvThickness      = 5.*cm;
  fCalorSizeYZ       = 10.*cm;
   
  // materials
  DefineMaterials();
  SetAbsorber1Material("G4_Fe");
  SetAbsorber2Material("G4_Fe");
  SetAl1Material("G4_Al");
  SetAl2Material("G4_Al");
  SetAr1Material("G4_Ar");
  SetAr2Material("G4_Ar");
  SetPCB1Material();
  SetPCB2Material();
  
  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
// use G4-NIST materials data base
//
G4NistManager* man = G4NistManager::Instance();
fDefaultMaterial = man->FindOrBuildMaterial("G4_AIR");
man->FindOrBuildMaterial("G4_Fe");
man->FindOrBuildMaterial("G4_Ar");

// print table
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  G4double posx,posy,posz;
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();
  
  //     
  // World
  //
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);  //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,            //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                 G4ThreeVector(),           //at (0,0,0)
                                 fLogicWorld,             //its logical volume
                                 "World",                   //its name
                                 0,                         //its mother  volume
                                 false,                  //no boolean operation
                                 0);                        //copy number
  
  //                               
  // Calorimeter
  //  
  fSolidCalor=0; fLogicCalor=0; fPhysiCalor=0;
 // fSolidLayer=0; fLogicLayer=0; fPhysiLayer=0;
  
  if (fCalorThickness > 0.)  
    { fSolidCalor = new G4Box("Calorimeter",                //its name
                    4*fCalorThickness/2,4*fCalorSizeYZ,4*fCalorSizeYZ);//size
                                 
      fLogicCalor = new G4LogicalVolume(fSolidCalor,        //its solid
                                        fDefaultMaterial,   //its material
                                        "Calorimeter");     //its name
                                           
      fPhysiCalor = new G4PVPlacement(0,                    //no rotation
                                     G4ThreeVector(),       //at (0,0,0)
                                     fLogicCalor,           //its logical volume
                                     "Calorimeter",         //its name
                                     fLogicWorld,           //its mother  volume
                                     false,              //no boolean operation
                                     0);                    //copy number
  
   }
  //                               
  // Absorber1
  //
  fSolidAbsorber1=0; fLogicAbsorber1=0; fPhysiAbsorber1=0;  
  
  if (fAbsorber1Thickness > 0.) 
    { fSolidAbsorber1 = new G4Box("Absorber1",                //its name
                          fAbsorber1Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      fLogicAbsorber1 = new G4LogicalVolume(fSolidAbsorber1,    //its solid
                                            fAbsorber1Material, //its material
                                            fAbsorber1Material->GetName());//name
                                                
      fPhysiAbsorber1 = new G4PVPlacement(0,                   //no rotation
					  G4ThreeVector(fAbsorber1Thickness/2,0.,0.),  //its position
					  fLogicAbsorber1,     //its logical volume
					  fAbsorber1Material->GetName(), //its name
					  fLogicCalor,          //its mother
					  false,               //no boulean operat
					  0);                   //copy number
                                        
    }
  
  // Al1
  fSolidAl1 = 0; fLogicAl1 = 0; fPhysiAl1 = 0;
  if (fAl1Thickness > 0.){
    fSolidAl1 = new G4Box("Al1",
			 fAl1Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl1 = new G4LogicalVolume(fSolidAl1,
				   fAl1Material,
				   fAl1Material->GetName());
   
     
    posx = fAbsorber1Thickness + (fAl1Thickness/2.);

    fPhysiAl1= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl1,           //its logical volume
				 "Al1",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    // G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl1->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar1--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr1[i][j] = new G4Box("Ar1",
				      fAr1Thickness/2,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr1 = new G4LogicalVolume(fSolidAr1[i][j],
					  fAr1Material,
					  fAr1Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + (fAr1Thickness/2.);
	 
	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr1[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr1,           //its logical volume
					    "Ar1",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr1->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }


   // PCB1
  fSolidPCB1 = 0; fLogicPCB1 = 0; fPhysiPCB1 = 0;
  if (fPCB1Thickness > 0.){
    fSolidPCB1 = new G4Box("PCB1",
			 fPCB1Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB1 = new G4LogicalVolume(fSolidPCB1,
				     fPCB1Material,
				     fPCB1Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + (fPCB1Thickness/2.);
    
    fPhysiPCB1= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicPCB1,           //its logical volume
				 "PCB1",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    // G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB1->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }
 
  // Absorber2
  fSolidAbsorber2=0; fLogicAbsorber2=0; fPhysiAbsorber2=0;  
  
  if (fAbsorber2Thickness > 0.) 
    { fSolidAbsorber2 = new G4Box("Absorber2",                //its name
                          fAbsorber2Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      fLogicAbsorber2 = new G4LogicalVolume(fSolidAbsorber2,    //its solid
                                            fAbsorber2Material, //its material
                                            fAbsorber2Material->GetName());//name

      posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + (fAbsorber2Thickness/2.);


                                                
      fPhysiAbsorber2 = new G4PVPlacement(0,                   //no rotation
					  G4ThreeVector(posx,0.,0.),  //its position
					  fLogicAbsorber2,     //its logical volume
					  fAbsorber2Material->GetName(), //its name
					  fLogicCalor,          //its mother
					  false,               //no boulean operat
					  0);                   //copy number
                                        
    }
  

  // Al2
  fSolidAl2 = 0; fLogicAl2 = 0; fPhysiAl2 = 0;
  if (fAl2Thickness > 0.){
    fSolidAl2 = new G4Box("Al2",
			 fAl2Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl2 = new G4LogicalVolume(fSolidAl2,
				   fAl2Material,
				   fAl2Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + (fAl2Thickness/2.);
    
    fPhysiAl2= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl2,           //its logical volume
				 "Al2",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl2->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar2--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr2[i][j] = new G4Box("Ar2",
				      fAr2Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr2 = new G4LogicalVolume(fSolidAr2[i][j],
					  fAr2Material,
					  fAr2Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + (fAr2Thickness/2.);
	 	 
	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr2[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr2,           //its logical volume
					    "Ar2",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr2->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);
	 fLogicAbsorber2->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }

   // PCB2
  fSolidPCB2 = 0; fLogicPCB2 = 0; fPhysiPCB2 = 0;
  if (fPCB2Thickness > 0.){
    fSolidPCB2 = new G4Box("PCB2",
			 fPCB2Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB2 = new G4LogicalVolume(fSolidPCB2,
				     fPCB2Material,
				     fPCB2Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + (fPCB2Thickness/2.);
    
    fPhysiPCB2= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicPCB2,           //its logical volume
				 "PCB2",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB2->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

  //
  //always return the physical World
  //
  return fPhysiWorld;

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
 /* G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << fEnvThickness << " layers of: [ "
         << fAbsorber1Thickness/mm << "mm of " << fAbsorber1Material->GetName() 
         << " + "
         << fGapThickness/mm << "mm of " << fGapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
*/}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber1Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);  
  if (pttoMaterial)
  {
      fAbsorber1Material = pttoMaterial;
      if ( fLogicAbsorber1 )
      {
          fLogicAbsorber1->SetMaterial(fAbsorber1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber2Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber2Material = pttoMaterial;
      if ( fLogicAbsorber2 )
      {
          fLogicAbsorber2->SetMaterial(fAbsorber2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl1Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl1Material = pttoMaterial;
      if ( fLogicAl1 )
	{
          fLogicAl1->SetMaterial(fAl1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl2Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl2Material = pttoMaterial;
      if ( fLogicAl2 )
	{
          fLogicAl2->SetMaterial(fAl2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAr1Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);  
  if (pttoMaterial)
    {
      fAr1Material = pttoMaterial;
      if ( fLogicAr1 )
	{
	  fLogicAr1->SetMaterial(fAr1Material);
	  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAr2Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr2Material = pttoMaterial;
      if ( fLogicAr2 )
	{
          fLogicAr2->SetMaterial(fAr2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPCB1Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB1Material = pttoMaterial;
      if ( fLogicPCB1 )
	{
          fLogicPCB1->SetMaterial(fPCB1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPCB2Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_0",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB2Material = pttoMaterial;
      if ( fLogicPCB2 )
	{
          fLogicPCB2->SetMaterial(fPCB2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber1Thickness(G4double val)
{
  // change Absorber1 thickness and recompute the calorimeter parameters
  fAbsorber1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber2Thickness(G4double val)
{
  // change Absorber2 thickness and recompute the calorimeter parameters
  fAbsorber2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl1Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl2Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAr1Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAr2Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPCB1Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPCB2Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fCalorSizeYZ = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEnvThickness(G4int val)
{
  fEnvThickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}











//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
