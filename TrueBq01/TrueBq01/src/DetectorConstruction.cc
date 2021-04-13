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
    fTargetMater(0), fLogicTarget(0),
    fDetectorMater(0), fLogicDetector(0),
    fWorldMater(0), fPhysiWorld(0),
    fDetectorMessenger(0)
{
    //Change this to side for a square
    //replace g4tubs to g4boxes
    //make side and thickness 

    //Target side and thickness (small square that will be on top)
    fTargetSide = 1.5 * mm; //make it the side
    fTargetThickness = 20 * um; //make it thickness of the small square

    //Detector side and thickness measures(Big square that will be down
    fDetectorSide = 5 * mm;
    fDetectorThickness = 0.8 * mm;



    //A way to make the word bigger
    fWorldSide = std::max(fTargetSide, fDetectorSide);
    fWorldThickness = fTargetThickness + fDetectorThickness; //Should change the name of this variable


    DefineMaterials();

    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fDetectorMessenger;
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

    fTargetMater = man->FindOrBuildMaterial("G4_Au"); // changed from "G4_CESIUM_IODIDE" by Ryan 24Feb2021


    fDetectorMater = man->FindOrBuildMaterial("G4_Si"); // changed from Germanium by Ryan 24Feb2021

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

    // World
    //
    // (re) compute World dimensions if necessary
    fWorldSide = std::max(fTargetSide, fDetectorSide);
    fWorldThickness = fTargetThickness + 2 * fDetectorThickness;

    //Replace by G4box for the square remove
    G4Box*
        sWorld = new G4Box("World",                                 //name
            0.5 * fWorldSide, 0.5 * fWorldSide, 0.5 * fWorldThickness); //dimensions  

    G4LogicalVolume*
        lWorld = new G4LogicalVolume(sWorld,                  //shape
            fWorldMater,               //material
            "World");                  //name

    fPhysiWorld = new G4PVPlacement(0,                    //no rotation
        G4ThreeVector(),            //at (0,0,0)
        lWorld,                     //logical volume
        "World",                    //name
        0,                          //mother volume
        false,                      //no boolean operation
        0);                         //copy number

// Target
//
    G4Box*
        sTarget = new G4Box("Target",                                   //name
            0.5 * fTargetSide, 0.5 * fTargetSide, 0.5 * fTargetThickness); //dimensions


    fLogicTarget = new G4LogicalVolume(sTarget,           //shape
        fTargetMater,              //material
        "Target");                 //name

//Placement of the target in the word
//Changing the position of the small square
    new G4PVPlacement(0,                         //no rotation
        G4ThreeVector(0, 0, 0),             //at (0,0, ( 1/2small square thickness+ 1/2big square thickness) 
        fLogicTarget,                //logical volume
        "Target",                    //name
        lWorld,                      //mother  volume
        false,                       //no boolean operation
        0);                          //copy number

    G4Region* emTarget = new G4Region("Target"); // create region that is the same as logical volume. Used for setting special production cuts
    fLogicTarget->SetRegion(emTarget);
    emTarget->AddRootLogicalVolume(fLogicTarget);
    // Detector
    //
    G4Box*
        sDetector = new G4Box("Detector",
            0.5 * fDetectorSide, 0.5 * fDetectorSide, 0.5 * fDetectorThickness);


    fLogicDetector = new G4LogicalVolume(sDetector,       //shape
        fDetectorMater,            //material
        "Detector");               //name

    new G4PVPlacement(0,                         //no rotation
        G4ThreeVector(0, 0, ((-0.5 * fTargetThickness) - 0.5 * fDetectorThickness)),             //at (0,0,0)
        fLogicDetector,              //logical volume
        "Detector",                  //name
        lWorld,                      //mother  volume
        false,                       //no boolean operation
        0);                          //copy number


    PrintParameters();

    //always return the root volume
    //
    return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    G4cout << "\n Target : Side = " << G4BestUnit(fTargetSide, "Length")
        << " Thickness = " << G4BestUnit(fTargetThickness, "Length")
        << " Material = " << fTargetMater->GetName();
    G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorSide, "Length")
        << " Tickness = " << G4BestUnit(fDetectorThickness, "Length")
        << " Material = " << fDetectorMater->GetName() << G4endl;
    G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial) {
        fTargetMater = pttoMaterial;
        if (fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
    else {
        G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
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
        fDetectorMater = pttoMaterial;
        if (fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
    else {
        G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
            << materialChoice << " not found" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetThickness(G4double value)
{
    fTargetThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetSide(G4double value)
{
    fTargetSide = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
    fDetectorThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorSide(G4double value)
{
    fDetectorSide = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetSide()
{
    return fTargetSide;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetThickness()
{
    return fTargetThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
    return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
    return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorSide()
{
    return fDetectorSide;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
    return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
    return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
    return fLogicDetector;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......