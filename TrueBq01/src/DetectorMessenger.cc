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
 detector(Det), fRdecayDir(0), fDetDir(0),
 fAbsorberMatCmd(0), fChipMatCmd(0), fAbsorberThicknessCmd(0),
 fChipThicknessCmd(0), fAbsorberSideCmd(0), fChipSideCmd(0), fActivitySideCmd(0), fActivityThicknessCmd(0)
{ 
  fRdecayDir = new G4UIdirectory("/TrueBq01/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/TrueBq01/det/",broadcast);
  fDetDir->SetGuidance("Detector construction commands");
        
  fAbsorberMatCmd = new G4UIcmdWithAString("/TrueBq01/det/setAbsorberMate",this);
  fAbsorberMatCmd->SetGuidance("Select material of the absorber");
  fAbsorberMatCmd->SetParameterName("choice",false);
  fAbsorberMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAbsorberThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setAbsorberThickness", this); //switched from setAbsorberRadius to setAbsorberThickness
  fAbsorberThicknessCmd->SetGuidance("Set the Absorber Thickness."); //switched radius to thickness
  fAbsorberThicknessCmd->SetUnitCategory("Length");
  fAbsorberThicknessCmd->SetParameterName("choice",false);
  fAbsorberThicknessCmd->AvailableForStates(G4State_PreInit);  

  
  fAbsorberSideCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setAbsorberSide", this); //switched to setAbsorberSide
  fAbsorberSideCmd->SetGuidance("Set the Absorber Side.");
  fAbsorberSideCmd->SetUnitCategory("Length");
  fAbsorberSideCmd->SetParameterName("choice",false);
  fAbsorberSideCmd->AvailableForStates(G4State_PreInit);
  

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

  fActivityThicknessCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setActivityThickness", this);
  fActivityThicknessCmd->SetGuidance("Set the Activity Thickness.");
  fActivityThicknessCmd->SetUnitCategory("Length");
  fActivityThicknessCmd->SetParameterName("choice", false);
  fActivityThicknessCmd->AvailableForStates(G4State_PreInit);

  fActivitySideCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setActivitySide", this);//switch it to side
  fActivitySideCmd->SetGuidance("Set the Activity Side.");
  fActivitySideCmd->SetUnitCategory("Length");
  fActivitySideCmd->SetParameterName("choice", false);
  fActivitySideCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
 
  delete fActivitySideCmd;
  delete fActivityThicknessCmd;
  delete fAbsorberMatCmd;
  delete fChipMatCmd;
  delete fAbsorberThicknessCmd;
  delete fChipThicknessCmd;
  delete fAbsorberSideCmd;
  delete fChipSideCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fAbsorberMatCmd )
   { detector->SetAbsorberMaterial(newValue);}
   
  if (command == fAbsorberSideCmd ) 
    {
      detector->SetAbsorberSide(fAbsorberSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fAbsorberThicknessCmd ) 
    {
      detector->SetAbsorberThickness(fAbsorberSideCmd->GetNewDoubleValue(newValue));}
    
  if (command == fChipMatCmd )
    {
      detector->SetDetectorMaterial(newValue);}
    
  if (command == fChipSideCmd ) 
    {
      detector->SetDetectorSide(
                     fChipSideCmd->GetNewDoubleValue(newValue));}

  if (command == fChipThicknessCmd ) 
    {
      detector->SetDetectorThickness(
                     fChipThicknessCmd->GetNewDoubleValue(newValue));}      


  if (command == fActivitySideCmd)
  {
      detector->SetActivitySide(fActivitySideCmd->GetNewDoubleValue(newValue));
  }

  if (command == fActivityThicknessCmd)
  {
      detector->SetActivityThickness(fAbsorberSideCmd->GetNewDoubleValue(newValue));
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
