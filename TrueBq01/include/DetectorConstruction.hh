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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class G4UserLimits;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();
    

    void SetActivitySide(G4double value);
    G4double GetActivitySide();
    void SetActivityThickness(G4double value);
    G4double GetActivityThickness();

    void SetActivityZOffset(G4double value); // e.g. for external source
    void SetThetaMax(G4double value);
    void SetThetaMin(G4double value);

    void SetAbsorberSide (G4double value);
    void SetAbsorberThickness (G4double value);
    void SetAbsorberMaterial (G4String);
   

    void SetChipWidth(G4double value);
    void SetChipLength(G4double value);           
    void SetChipThickness(G4double value);  
    void SetChipMaterial(G4String);               
                   
    void PrintParameters();
    
  public:
      
    G4double GetAbsorberSide();
    G4double GetAbsorberThickness();
    G4Material* GetAbsorberMaterial();       
    G4LogicalVolume* GetLogicAbsorber();
    G4double GetAbsorberXOffset();
    G4ThreeVector GetAbsorberCenter();

    G4double GetActivityZOffset(); // e.g. for external source
    G4double GetThetaMax(); 
    G4double GetThetaMin();

    G4double GetChipWidth();
    G4double GetChipLength();
    G4double GetChipThickness();
    G4Material* GetChipMaterial();                 
    G4LogicalVolume* GetLogicChip(); 
    G4LogicalVolume* GetLogicTES();

    G4double GetEres();
    void SetEres(G4double);

    G4double GetEtail();
    void SetEtail(G4double);
    G4double GetPtail();
    void SetPtail(G4double);

    G4double GetEtail2();
    void SetEtail2(G4double);
    G4double GetPtail2();
    void SetPtail2(G4double);

    G4double GetEtailH();
    void SetEtailH(G4double);
    G4double GetPtailH();
    void SetPtailH(G4double);

    void SetParentOnly(G4bool val);
    G4bool GetParentOnly();

                       
  private:
  
    G4String myParticleName; // for histogramming recoil atoms
    G4double            fActivityThickness; // activity region
    G4double            fActivitySide;      
    G4double           fAbsorberSide; 
    G4double           fAbsorberThickness;
    G4Material*        fAbsorberMater;
    G4Material*         fSMater;
    G4LogicalVolume*   fLogicAbsorber;

    G4double fEres; // energy resolution
    G4double fEtail; // low energy tail res
    G4double fPtail; // low energy tail prob
    G4double fEtail2; // low energy tail res
    G4double fPtail2; // low energy tail prob
    G4double fEtailH; // low energy tail res
    G4double fPtailH; // low energy tail prob

    G4double fAbsorberXOffset;
    G4double fTESside;
    G4double fTESThickness;
    G4double fTESXOffset;
    G4ThreeVector AbsorberCenter;
    G4double flegThickness;
    G4double fChipBorderThickness;

    G4double fActivityZOffset; 
    G4double fThetaMax; // for primary generator trajectories
    G4double fThetaMin; // for primary generator trajectories
    G4Material* fTESMater;

    G4double           fChipWidth;
    G4double           fChipLength;
    G4double           fChipThickness;
    G4Material*        fChipMater;
    G4LogicalVolume*   fLogicChip;
    G4LogicalVolume*   fLogicTES;
               
    G4double           fWorldSide;
    G4double           fWorldThickness;
    G4Material*        fWorldMater;     
    G4VPhysicalVolume* fPhysiWorld;
                
    DetectorMessenger* detectorMessenger;

    G4bool  ParentOnly; // limit rad decay to parent

    G4UserLimits* fStepLimit = nullptr; // pointer to user step limits

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

