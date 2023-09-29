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
/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:
  
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
      //radius=thickness
      //length=side
    DetectorConstruction*      detector;
    
    G4UIdirectory*             fRdecayDir;
    G4UIdirectory*             fDetDir;
    G4UIcmdWithAString*        fAbsorberMatCmd;
    G4UIcmdWithAString*        fChipMatCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fChipThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberSideCmd;
    G4UIcmdWithADoubleAndUnit* fChipWidthCmd;  
    G4UIcmdWithADoubleAndUnit* fChipLengthCmd;

    G4UIcmdWithADoubleAndUnit* fActivityZOffsetCmd;
    G4UIcmdWithADoubleAndUnit* fThetaMaxCmd;
    G4UIcmdWithADoubleAndUnit* fThetaMinCmd;

    G4UIcmdWithADoubleAndUnit* fActivitySideCmd;
    G4UIcmdWithADoubleAndUnit* fActivityThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fEresCmd; // energy resolution
    G4UIcmdWithADoubleAndUnit* fEtailCmd; // energy resolution
    G4UIcmdWithADouble* fPtailCmd; // energy resolution
    G4UIcmdWithADoubleAndUnit* fEtail2Cmd; // energy resolution
    G4UIcmdWithADouble* fPtail2Cmd; // energy resolution
    G4UIcmdWithADoubleAndUnit* fEtailHCmd; // energy resolution
    G4UIcmdWithADouble* fPtailHCmd; // energy resolution

    G4UIcmdWithABool* fParentOnlyCmd; // limit rad decay to parent (main nuclide)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

