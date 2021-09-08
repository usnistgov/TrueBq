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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4Region.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    :G4VUserDetectorConstruction(),
    fAbsorberMater(0), fLogicAbsorber(0),
    fChipMater(0), fLogicChip(0),
    fWorldMater(0), fPhysiWorld(0),
    detectorMessenger(0)
{
    // Solid Square absorber sitting on Solid Square Chip
    // Solid Square "activity" used by primary generator

    // Activity in the center of absorber by default. Can change in macro
    fActivitySide = 0.0;
    fActivityThickness = 0.0;

    //Absorber side and thickness (small square that will be on top)
    fAbsorberSide = 1.5 * mm; //make it the side
    fAbsorberThickness = 20 * um; //make it thickness of the small square

    //Detector side and thickness measures(Big square that will be down
    fChipSide = 5 * mm;
    fChipThickness = 0.8 * mm;



    //A way to make the word bigger
    fWorldSide = std::max(fAbsorberSide, fChipSide);
    fWorldThickness = fAbsorberThickness + fChipThickness; //Should change the name of this variable


    DefineMaterials();

    detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    // build materials
    //



    G4Element* N = new G4Element("Nitrogen", "N", 7, 14.01 * g / mole);
    G4Element* O = new G4Element("Oxygen", "O", 8, 16.00 * g / mole);
    //
    G4int ncomponents; G4double fractionmass;
    G4Material* Air20 = new G4Material("Air", 1.205 * mg / cm3, ncomponents = 2,
        kStateGas, 293. * kelvin, 1. * atmosphere);
    Air20->AddElement(N, fractionmass = 0.7);
    Air20->AddElement(O, fractionmass = 0.3);
    // or use G4 materials data base
    //
    G4NistManager* man = G4NistManager::Instance();

    //
    // fWorldMater = Air20;
    fWorldMater = man->FindOrBuildMaterial("G4_Galactic"); // RPF - changed from Air20. G4_Galactic is low density 10^-25 g/cm^3. 
    // To do: Add thin layer of air for backscatter
    //

    fAbsorberMater = man->FindOrBuildMaterial("G4_Au"); // changed from "G4_CESIUM_IODIDE" by Ryan 24Feb2021


    fChipMater = man->FindOrBuildMaterial("G4_Si"); // changed from Germanium by Ryan 24Feb2021

    // new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);

   // G4cout << *(G4Material::GetMaterialTable()) << G4endl; //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // ============  WORLD ================= //
    //
    // (re) compute World dimensions if necessary
    fWorldSide = std::max(fAbsorberSide, fChipSide);
    fWorldThickness = fAbsorberThickness + 2 * fChipThickness;


    G4Box*
        sWorld = new G4Box("World",                                     //name
            0.5 * fWorldSide, 0.5 * fWorldSide, 0.5 * fWorldThickness); //dimensions  

    G4LogicalVolume*
        lWorld = new G4LogicalVolume(sWorld,                            //shape
            fWorldMater,                                                //material
            "World");                                                   //name

    fPhysiWorld = new G4PVPlacement(0,      //no rotation
        G4ThreeVector(),                    //at (0,0,0)
        lWorld,                             //logical volume
        "World",                            //name
        0,                                  //mother volume
        false,                              //no boolean operation
        0);                                 //copy number

 // ============  ABSORBER ================= //
//
    G4Box*
        sAbsorber = new G4Box("Absorber",                                   //name
            0.5 * fAbsorberSide, 0.5 * fAbsorberSide, 0.5 * fAbsorberThickness); //dimensions


    fLogicAbsorber = new G4LogicalVolume(sAbsorber,         //shape
        fAbsorberMater,                                     //material
        "Absorber");                                        //name

//Placement of the absorber in the word

    new G4PVPlacement(0,                //no rotation
        G4ThreeVector(0, 0, 0),         //at (0,0,0)
        fLogicAbsorber,                 //logical volume
        "Absorber",                     //name
        lWorld,                         //mother  volume
        false,                          //no boolean operation
        0);                             //copy number

    G4Region* emAbsorber = new G4Region("Absorber"); // create region that is the same as logical volume. Used for setting special production cuts
    fLogicAbsorber->SetRegion(emAbsorber);
    emAbsorber->AddRootLogicalVolume(fLogicAbsorber);
   
 // ============  DETECTOR ================= //
    
    G4Box*
        sDetector = new G4Box("Chip",
            0.5 * fChipSide, 0.5 * fChipSide, 0.5 * fChipThickness);


    fLogicChip = new G4LogicalVolume(sDetector,       //shape
        fChipMater,            //material
        "Chip");               //name

    new G4PVPlacement(0,                         //no rotation
        G4ThreeVector(0, 0, ((-0.5 * fAbsorberThickness) - 0.5 * fChipThickness)),             //at (0,0,0)
        fLogicChip,              //logical volume
        "Chip",                  //name
        lWorld,                  //mother  volume
        false,                   //no boolean operation
        0);                      //copy number


    PrintParameters();

    //always return the root volume
    //
    return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    G4cout << "\n Absorber : Side = " << G4BestUnit(fAbsorberSide, "Length")
        << " Thickness = " << G4BestUnit(fAbsorberThickness, "Length")
        << " Material = " << fAbsorberMater->GetName();
    G4cout << "\n Chip : Length = " << G4BestUnit(fChipSide, "Length")
        << " Tickness = " << G4BestUnit(fChipThickness, "Length")
        << " Material = " << fChipMater->GetName() << G4endl;
    G4cout << "\n" << fAbsorberMater << "\n" << fChipMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial) {
        fAbsorberMater = pttoMaterial;
        if (fLogicAbsorber) { fLogicAbsorber->SetMaterial(fAbsorberMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
    else {
        G4cout << "\n--> warning from DetectorConstruction::SetAbsorberMaterial : "
            << materialChoice << " not found" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial) {
        fChipMater = pttoMaterial;
        if (fLogicChip) { fLogicChip->SetMaterial(fChipMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
    else {
        G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
            << materialChoice << " not found" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double value)
{
    fAbsorberThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}


G4double DetectorConstruction::GetActivitySide()
{
    return fActivitySide;
}
G4double DetectorConstruction::GetActivityThickness()
{
    return fActivityThickness;
}



void DetectorConstruction::SetActivitySide(G4double value)
{
    fActivitySide = value;
}


void DetectorConstruction::SetActivityThickness(G4double value)
{
    fActivityThickness = value;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberSide(G4double value)
{
    fAbsorberSide = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
    fChipThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorSide(G4double value)
{
    fChipSide = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetAbsorberSide()
{
    return fAbsorberSide;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetAbsorberThickness()
{
    return fAbsorberThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetAbsorberMaterial()
{
    return fAbsorberMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicAbsorber()
{
    return fLogicAbsorber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorSide()
{
    return fChipSide;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
    return fChipThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
    return fChipMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
    return fLogicChip;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......