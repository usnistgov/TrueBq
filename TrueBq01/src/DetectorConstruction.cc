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

#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    :G4VUserDetectorConstruction(),
    fAbsorberMater(0), fLogicAbsorber(0),
    fChipMater(0), fLogicChip(0),fLogicTES(),
    fWorldMater(0), fPhysiWorld(0),
    detectorMessenger(0), fSMater()
{
    // Solid Square absorber sitting on Solid Square Chip
    // Solid Square "activity" used by primary generator

    fEres = 0.0 * keV; // energy resolution. Set by user.

    // Activity in the center of absorber by default. Can change in macro
    fActivitySide = 0.0;
    fActivityThickness = 0.0;
    fActivityZOffset = 0.0; 

    //Absorber side and thickness (small square that will be on top)
    fAbsorberSide = 1.5 * mm; //make it the side
    fAbsorberThickness = 30 * um; //make it thickness of the small square

    // Main body of chip. Excludes legs
    fChipWidth = 7 * mm; // 10.55 mm total
    fChipLength = 7 * mm; // 10.55 mm total

    fChipThickness = 0.275 * mm;

    fAbsorberXOffset = -2.0 * mm; // offset from center of chip (actually moves chip not absorber)
    AbsorberCenter = G4ThreeVector(fAbsorberXOffset, 0.0, fAbsorberThickness / 2.0); // recomputed below, in case parameters changed in macro

    fTESside = 2.0 * mm;
    fTESThickness = 0.200 * um;
    fTESXOffset = 2. * mm; // ofset of  from center of chip

    flegThickness = 0.1 * mm; // legs
    fChipBorderThickness = 0.2 * mm; // border outside legs


    //A way to make the word bigger
    fWorldSide = std::max(fAbsorberSide, fChipLength);
    fWorldThickness = fAbsorberThickness + fChipThickness+ fActivityZOffset; // recomputed below after user input; Should change the name of this variable


    DefineMaterials();

    detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete detectorMessenger;
    delete fStepLimit;
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

    G4Material* PuNitrate = new G4Material("PuNitrate", 1.0 * g / cm3, ncomponents = 3);
    PuNitrate->AddElement(man->FindOrBuildElement("Pu"), 1);
    PuNitrate->AddElement(man->FindOrBuildElement("N"), 4);
    PuNitrate->AddElement(man->FindOrBuildElement("O"), 12);

    G4Material* thorium = new G4Material("thorium", 5.0 * g / cm3, ncomponents = 1);
    thorium->AddElement(man->FindOrBuildElement("Th"), 1);

    //
    // fWorldMater = Air20;
    fWorldMater = man->FindOrBuildMaterial("G4_Galactic"); // RPF - changed from Air20. G4_Galactic is low density 10^-25 g/cm^3. 
    // To do: Add thin layer of air for backscatter
    //

    fAbsorberMater = man->FindOrBuildMaterial("G4_Au"); // G4_Au, 
   

    fChipMater = man->FindOrBuildMaterial("G4_Si"); // G4_Si

    fTESMater = man->FindOrBuildMaterial("G4_Mo"); // G4_Mo

    fSMater = man->FindOrBuildMaterial("G4_Cu"); // G4_Ba for external source


    // Pu(NO3)4



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

 
 //========== RECOMPUTE DIMENSIONS AS NEEDED BAED ON USER INPUT =====================//
    
 // re-compute world to contain everything

    fWorldSide = 2.0 * (fChipLength);
    fWorldThickness = fAbsorberThickness + fTESThickness + 2 * fChipThickness + 2 * fActivityZOffset;

 // re-compute Absorber and TES positions based on their dimensions

    fAbsorberXOffset = -fAbsorberSide/2.0 - 0.1*mm; // leave small gap
    AbsorberCenter =  G4ThreeVector(0.,0.,0.);

    fTESXOffset = fTESside/2.0 + 0.1 * mm; // leave small gap
    G4ThreeVector TESCenter = G4ThreeVector(fTESXOffset-fAbsorberXOffset, 0.0, fAbsorberThickness / 2.0 + fTESThickness / 2.0);


 // ==============  WORLD =================== //
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

// ============  ABSORBER ====== center of absorber is origin of the world (0,0,0) ================= //
//
    G4Box* sAbsorber = new G4Box("Absorber", 0.5 * fAbsorberSide, 0.5 * fAbsorberSide, 0.5 * fAbsorberThickness); 

    fLogicAbsorber = new G4LogicalVolume(sAbsorber, fAbsorberMater, "Absorber");                               

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), fLogicAbsorber, "Absorber", lWorld,false, 0);

    G4Region* emAbsorber = new G4Region("Absorber"); // create region that is the same as logical volume. Used for setting special production cuts
    fLogicAbsorber->SetRegion(emAbsorber);
    emAbsorber->AddRootLogicalVolume(fLogicAbsorber);

    // STEP LIMITER IN ABSORBER //

    G4double maxStep = 1000. * micrometer; // To do, add this to macro
    fStepLimit = new G4UserLimits(maxStep);
    fLogicAbsorber->SetUserLimits(fStepLimit);

    // END STEP LIMITER        //

   
 // ============  CHIP  ================= //
    
    G4Box* sChip = new G4Box("Chip", 0.5 * fChipLength, 0.5 * fChipWidth, 0.5 * fChipThickness);
    fLogicChip = new G4LogicalVolume(sChip, fChipMater, "Chip"); 
    new G4PVPlacement(0, G4ThreeVector(-fAbsorberXOffset, 0.0, -0.5 * fAbsorberThickness - 0.5 * fChipThickness), fLogicChip, "Chip", lWorld, false, 0);        
  
 
    // ============  TES ================= //

    G4Box* sTES = new G4Box("TES", 0.5 * fTESside, 0.5 * fTESside, 0.5 * fTESThickness);
    fLogicTES = new G4LogicalVolume(sTES, fTESMater, "TES");
    new G4PVPlacement(0, G4ThreeVector(-fAbsorberXOffset+fTESXOffset, 0.0, - 0.5 * fAbsorberThickness + 0.5 * fTESThickness), fLogicTES, "TES", lWorld, false, 0);



    // ============  Ext Source Material offset to source location =====To do: Macro inputs for side, material; also add can etc.============ //

    G4double sMat_side = 4.0 * mm;
    G4Box* sSmat = new G4Box("Smat", 0.5 * sMat_side, 0.5 * sMat_side, 0.5 * sMat_side);

    G4LogicalVolume* fLogicSmat = new G4LogicalVolume(sSmat, fSMater, "Smat"); // make out of fSMater
    if (fActivityZOffset > sMat_side) // only make source material if source is far enough away from sensor
    {
        new G4PVPlacement(0, G4ThreeVector(0, 0, fActivityZOffset), fLogicSmat, "Smat", lWorld, false, 0);
    }

// --------------------------------------------------------------------------------------------- //

    
    PrintParameters();
    return fPhysiWorld;   //always return the root volume
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    G4cout << "\n Absorber : Side = " << G4BestUnit(fAbsorberSide, "Length")
        << " Thickness = " << G4BestUnit(fAbsorberThickness, "Length")
        << " Material = " << fAbsorberMater->GetName();
    G4cout << "\n Chip : Length = " << G4BestUnit(fChipLength, "Length");
        G4cout << "\n Chip : Width = " << G4BestUnit(fChipWidth, "Length")
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

void DetectorConstruction::SetChipMaterial(G4String materialChoice)
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
        G4cout << "\n--> warning from DetectorConstruction::SetChipMaterial : "
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

void DetectorConstruction::SetChipThickness(G4double value)
{
    fChipThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetChipLength(G4double value)
{
    fChipLength = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetChipWidth(G4double value)
{
    fChipWidth = value;
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

G4double DetectorConstruction::GetAbsorberXOffset()
{
    return fAbsorberXOffset;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetActivityZOffset(G4double value)
{
    fActivityZOffset = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

G4double DetectorConstruction::GetActivityZOffset()
{
    return fActivityZOffset;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetThetaMax(G4double value)
{
    fThetaMax = value;
 
}

G4double DetectorConstruction::GetThetaMax()
{
    return fThetaMax;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetThetaMin(G4double value)
{
    fThetaMin = value;

}

G4double DetectorConstruction::GetThetaMin()
{
    return fThetaMin;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector DetectorConstruction::GetAbsorberCenter()
{
    return AbsorberCenter;
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

G4double DetectorConstruction::GetChipWidth()
{
    return fChipWidth;
}

G4double DetectorConstruction::GetChipLength()
{
    return fChipLength;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetChipThickness()
{
    return fChipThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetChipMaterial()
{
    return fChipMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicChip()
{
    return fLogicChip;
}
G4LogicalVolume* DetectorConstruction::GetLogicTES()
{
    return fLogicTES;
}
G4double DetectorConstruction::GetEres()
{
    return fEres;
}
void DetectorConstruction::SetEres(G4double myEres)
{
    fEres = myEres;
   
}

G4double DetectorConstruction::GetEtail()
{
    return fEtail;
}
void DetectorConstruction::SetEtail(G4double myEtail)
{
    fEtail = myEtail;

}

G4double DetectorConstruction::GetPtail()
{
    return fPtail;
}
void DetectorConstruction::SetPtail(G4double myPtail)
{
    fPtail = myPtail;

}

G4double DetectorConstruction::GetEtail2()
{
    return fEtail2;
}
void DetectorConstruction::SetEtail2(G4double myEtail2)
{
    fEtail2 = myEtail2;

}

G4double DetectorConstruction::GetPtail2()
{
    return fPtail2;
}
void DetectorConstruction::SetPtail2(G4double myPtail2)
{
    fPtail2 = myPtail2;

}

G4double DetectorConstruction::GetEtailH()
{
    return fEtailH;
}
void DetectorConstruction::SetEtailH(G4double myEtailH)
{
    fEtailH = myEtailH;

}

G4double DetectorConstruction::GetPtailH()
{
    return fPtailH;
}
void DetectorConstruction::SetPtailH(G4double myPtailH)
{
    fPtailH = myPtailH;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetParentOnly(G4bool val)
{
    ParentOnly = val;
    G4cout << "IN DETECTOR ParentONLY = " << ParentOnly << G4endl;
}

G4bool DetectorConstruction::GetParentOnly()
{
    return ParentOnly;
}
