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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:G4UImessenger(), 

//2 Dir field and 6commands line fields
 fChip(Det), fRdecayDir(0), fDetDir(0),
 fTargMatCmd(0), fChipMatCmd(0), fTargThicknessCmd(0),
 fChipThicknessCmd(0), fTargSideCmd(0), fChipSideCmd(0) 
{ 
  fRdecayDir = new G4UIdirectory("/TrueBq01/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/TrueBq01/det/",broadcast);
  fDetDir->SetGuidance("Detector construction commands");
        
  fTargMatCmd = new G4UIcmdWithAString("/TrueBq01/det/setAbsorberMate",this);
  fTargMatCmd->SetGuidance("Select material of the absorber");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fTargThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setAbsorberThickness", this); //switched from setAbsorberRadius to setAbsorberThickness
  fTargThicknessCmd->SetGuidance("Set the Absorber Thickness."); //switched radius to thickness
  fTargThicknessCmd->SetUnitCategory("Length");
  fTargThicknessCmd->SetParameterName("choice",false);
  fTargThicknessCmd->AvailableForStates(G4State_PreInit);  

  
  fTargSideCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setAbsorberSide", this); //switched to setAbsorberSide
  fTargSideCmd->SetGuidance("Set the Absorber Side.");
  fTargSideCmd->SetUnitCategory("Length");
  fTargSideCmd->SetParameterName("choice",false);
  fTargSideCmd->AvailableForStates(G4State_PreInit);
  

  fChipMatCmd = new G4UIcmdWithAString("/TrueBq01/det/setChipMate",this);
  fChipMatCmd->SetGuidance("Select Material of the Chip.");
  fChipMatCmd->SetParameterName("choice",false);
  fChipMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fChipThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setChipThickness",this);
  fChipThicknessCmd->SetGuidance("Set the Chip Thickness.");
  fChipThicknessCmd->SetUnitCategory("Length");
  fChipThicknessCmd->SetParameterName("choice",false);
  fChipThicknessCmd->AvailableForStates(G4State_PreInit);

  fChipSideCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setChipSide",this);//switch it to side
  fChipSideCmd->SetGuidance("Set the Chip Side.");
  fChipSideCmd->SetUnitCategory("Length");
  fChipSideCmd->SetParameterName("choice",false);
  fChipSideCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTargMatCmd;
  delete fChipMatCmd;
  delete fTargThicknessCmd;
  delete fChipThicknessCmd;
  delete fTargSideCmd;
  delete fChipSideCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fTargMatCmd )
   { fChip->SetAbsorberMaterial(newValue);}
   
  if (command == fTargSideCmd ) 
    { fChip->SetAbsorberSide(fTargSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fTargThicknessCmd ) 
    {fChip->SetAbsorberThickness(fTargSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fChipMatCmd )
    { fChip->SetDetectorMaterial(newValue);}
    
  if (command == fChipSideCmd ) 
    {fChip->SetDetectorSide(
                     fChipSideCmd->GetNewDoubleValue(newValue));}

  if (command == fChipThicknessCmd ) 
    {fChip->SetDetectorThickness(
                     fChipThicknessCmd->GetNewDoubleValue(newValue));}      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
