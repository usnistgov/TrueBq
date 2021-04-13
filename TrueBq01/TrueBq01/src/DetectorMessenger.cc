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
 fDetector(Det), fRdecayDir(0), fDetDir(0),
 fTargMatCmd(0), fDetectMatCmd(0), fTargThicknessCmd(0),
 fDetectThicknessCmd(0), fTargSideCmd(0), fDetectSideCmd(0) 
{ 
  fRdecayDir = new G4UIdirectory("/TrueBq01/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/TrueBq01/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fTargMatCmd = new G4UIcmdWithAString("/TrueBq01/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select material of the target");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fTargThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setTargetThickness", this); //switched from setTargetRadius to setTargetThickness
  fTargThicknessCmd->SetGuidance("Set the Target Thickness."); //switched radius to thickness
  fTargThicknessCmd->SetUnitCategory("Length");
  fTargThicknessCmd->SetParameterName("choice",false);
  fTargThicknessCmd->AvailableForStates(G4State_PreInit);  

  
  fTargSideCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setTargetSide", this); //switched to setTargetSide
  fTargSideCmd->SetGuidance("Set the Target Side.");
  fTargSideCmd->SetUnitCategory("Length");
  fTargSideCmd->SetParameterName("choice",false);
  fTargSideCmd->AvailableForStates(G4State_PreInit);
  

  fDetectMatCmd = new G4UIcmdWithAString("/TrueBq01/det/setDetectorMate",this);
  fDetectMatCmd->SetGuidance("Select Material of the Detector.");
  fDetectMatCmd->SetParameterName("choice",false);
  fDetectMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fDetectThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setDetectorThickness",this);
  fDetectThicknessCmd->SetGuidance("Set the Detector Thickness.");
  fDetectThicknessCmd->SetUnitCategory("Length");
  fDetectThicknessCmd->SetParameterName("choice",false);
  fDetectThicknessCmd->AvailableForStates(G4State_PreInit);

  fDetectSideCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setDetectorSide",this);//switch it to side
  fDetectSideCmd->SetGuidance("Set the Detector Side.");
  fDetectSideCmd->SetUnitCategory("Length");
  fDetectSideCmd->SetParameterName("choice",false);
  fDetectSideCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTargMatCmd;
  delete fDetectMatCmd;
  delete fTargThicknessCmd;
  delete fDetectThicknessCmd;
  delete fTargSideCmd;
  delete fDetectSideCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fTargMatCmd )
   { fDetector->SetTargetMaterial(newValue);}
   
  if (command == fTargSideCmd ) 
    { fDetector->SetTargetSide(fTargSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fTargThicknessCmd ) 
    {fDetector->SetTargetThickness(fTargSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fDetectMatCmd )
    { fDetector->SetDetectorMaterial(newValue);}
    
  if (command == fDetectSideCmd ) 
    {fDetector->SetDetectorSide(
                     fDetectSideCmd->GetNewDoubleValue(newValue));}

  if (command == fDetectThicknessCmd ) 
    {fDetector->SetDetectorThickness(
                     fDetectThicknessCmd->GetNewDoubleValue(newValue));}      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
