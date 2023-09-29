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
#include "Run.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:G4UImessenger(), 

//2 Dir field and 6commands line fields
 detector(Det), fRdecayDir(0), fDetDir(0),
 fAbsorberMatCmd(0), fChipMatCmd(0), fAbsorberThicknessCmd(0),
 fChipThicknessCmd(0), fAbsorberSideCmd(0), fChipLengthCmd(0), fChipWidthCmd(0), fActivitySideCmd(0), fActivityThicknessCmd(0), 
 fEresCmd(0), fActivityZOffsetCmd(),fThetaMaxCmd(), fThetaMinCmd(), fEtailCmd(), fPtailCmd()
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

  fChipWidthCmd =
       new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setChipWidth",this);//width
  fChipWidthCmd->SetGuidance("Set the Chip width.");
  fChipWidthCmd->SetUnitCategory("Length");
  fChipWidthCmd->SetParameterName("choice",false);
  fChipWidthCmd->AvailableForStates(G4State_PreInit);

  fChipLengthCmd =
      new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setChipLength", this);//length
  fChipLengthCmd->SetGuidance("Set the Chip length.");
  fChipLengthCmd->SetUnitCategory("Length");
  fChipLengthCmd->SetParameterName("choice", false);
  fChipLengthCmd->AvailableForStates(G4State_PreInit);

  fActivityThicknessCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setActivityThickness", this);
  fActivityThicknessCmd->SetGuidance("Set the Activity Thickness.");
  fActivityThicknessCmd->SetUnitCategory("Length");
  fActivityThicknessCmd->SetParameterName("choice", false);
  fActivityThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fActivitySideCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setActivitySide", this);//switch it to side
  fActivitySideCmd->SetGuidance("Set the Activity Side.");
  fActivitySideCmd->SetUnitCategory("Length");
  fActivitySideCmd->SetParameterName("choice", false);
  fActivitySideCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fActivityZOffsetCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setActivityZOffset", this);// Z offset of activity center
  fActivityZOffsetCmd->SetGuidance("Set the Activity Z offset.");
  fActivityZOffsetCmd->SetUnitCategory("Length");
  fActivityZOffsetCmd->SetParameterName("choice", false);
  fActivityZOffsetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fThetaMaxCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setThetaMax", this);// max trajectory angle
  fThetaMaxCmd->SetGuidance("Set the maximum angle");
  fThetaMaxCmd->SetUnitCategory("Angle");
  fThetaMaxCmd->SetParameterName("choice", false);
  fThetaMaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fThetaMinCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setThetaMin", this);// min trajectory angle
  fThetaMinCmd->SetGuidance("Set the minimum angle");
  fThetaMinCmd->SetUnitCategory("Angle");
  fThetaMinCmd->SetParameterName("choice", false);
  fThetaMinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEresCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setEres", this);//energy resolution
  fEresCmd->SetGuidance("Set the energy resolution");
  fEresCmd->SetUnitCategory("Energy");
  fEresCmd->SetParameterName("choice", false);
  fEresCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEtailCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setEtail", this);//tail resolution
  fEtailCmd->SetGuidance("Set the low energy tail");
  fEtailCmd->SetUnitCategory("Energy");
  fEtailCmd->SetParameterName("choice", false);
  fEtailCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPtailCmd = new G4UIcmdWithADouble("/TrueBq01/det/setPtail", this);//tail prob
  fPtailCmd->SetGuidance("Set the low energy tail prob");
  fPtailCmd->SetParameterName("choice", false);
  fPtailCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEtail2Cmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setEtail2", this);//tail 2 resolution
  fEtail2Cmd->SetGuidance("Set the low energy tail 2");
  fEtail2Cmd->SetUnitCategory("Energy");
  fEtail2Cmd->SetParameterName("choice", false);
  fEtail2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPtail2Cmd = new G4UIcmdWithADouble("/TrueBq01/det/setPtail2", this);//tail 2 prob
  fPtail2Cmd->SetGuidance("Set the low energy tail 2 prob");
  fPtail2Cmd->SetParameterName("choice", false);
  fPtail2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEtailHCmd = new G4UIcmdWithADoubleAndUnit("/TrueBq01/det/setEtailH", this);//tail 2 resolution
  fEtailHCmd->SetGuidance("Set the low energy tail H");
  fEtailHCmd->SetUnitCategory("Energy");
  fEtailHCmd->SetParameterName("choice", false);
  fEtailHCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPtailHCmd = new G4UIcmdWithADouble("/TrueBq01/det/setPtailH", this);//tail 2 prob
  fPtailHCmd->SetGuidance("Set the low energy tail H prob");
  fPtailHCmd->SetParameterName("choice", false);
  fPtailHCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


  fParentOnlyCmd = new G4UIcmdWithABool("/TrueBq01/det/ParentOnly", this);
  fParentOnlyCmd->SetGuidance("TRUE for parent only, no chain");
  fParentOnlyCmd->SetParameterName("ParentOnly", false);
  fParentOnlyCmd->AvailableForStates(G4State_PreInit, G4State_Idle); // check this
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
  delete fChipLengthCmd;
  delete fChipWidthCmd;
  delete fDetDir;
  delete fRdecayDir;  
  delete fEresCmd;
  delete fEtailCmd;
  delete fPtailCmd;
  delete fEtail2Cmd;
  delete fPtail2Cmd;
  delete fEtailHCmd;
  delete fPtailHCmd;
  delete  fActivityZOffsetCmd;
  delete fThetaMaxCmd;
  delete fParentOnlyCmd;
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
      detector->SetChipMaterial(newValue);}
    
  if (command == fChipWidthCmd ) 
    {
      detector->SetChipWidth(
                     fChipWidthCmd->GetNewDoubleValue(newValue));}

  if (command == fChipLengthCmd)
  {
      detector->SetChipLength(
          fChipLengthCmd->GetNewDoubleValue(newValue));
  }
  if (command == fChipThicknessCmd ) 
    {
      detector->SetChipThickness(
                     fChipThicknessCmd->GetNewDoubleValue(newValue));}      


  if (command == fActivitySideCmd)
  {
      detector->SetActivitySide(fActivitySideCmd->GetNewDoubleValue(newValue));
  }

  if (command == fActivityThicknessCmd)
  {
      detector->SetActivityThickness(fAbsorberSideCmd->GetNewDoubleValue(newValue));
  }

  if (command == fActivityZOffsetCmd)
  {
      detector->SetActivityZOffset(fActivityZOffsetCmd->GetNewDoubleValue(newValue));

  }
  if (command == fThetaMaxCmd)
  {
      detector->SetThetaMax(fThetaMaxCmd->GetNewDoubleValue(newValue));

  }
  if (command == fThetaMinCmd)
  {
      detector->SetThetaMin(fThetaMinCmd->GetNewDoubleValue(newValue));

  }
  if (command == fEresCmd) // pass energy resolution to the run
  {
      detector->SetEres(fEresCmd->GetNewDoubleValue(newValue));
  }
  if (command == fEtailCmd) // pass tail resolution to the run
  {
      detector->SetEtail(fEtailCmd->GetNewDoubleValue(newValue));
  }
  if (command == fPtailCmd) // pass tail prob
  {
      detector->SetPtail(fPtailCmd->GetNewDoubleValue(newValue));
  }
  if (command == fEtail2Cmd) // pass tail 2 resolution to the run
  {
      detector->SetEtail2(fEtail2Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fPtail2Cmd) // pass tail 2 prob
  {
      detector->SetPtail2(fPtail2Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fEtailHCmd) // pass tail H resolution to the run
  {
      detector->SetEtailH(fEtailHCmd->GetNewDoubleValue(newValue));
  }
  if (command == fPtailHCmd) // pass tail 2 prob
  {
      detector->SetPtailH(fPtailHCmd->GetNewDoubleValue(newValue));
  }


  if (command = fParentOnlyCmd)
  {
      detector->SetParentOnly(fParentOnlyCmd->GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
