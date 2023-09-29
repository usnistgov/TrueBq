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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4timer.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void WriteAnAscii();
    
    void AddEdep (G4int iVol, G4double Edep, G4double time, G4double weight);
    void SetParticleCount(G4int n);
    void IncrementParticleCount();
    G4int GetParticleCount();
    G4int GeteventID() { return eventID; };
  private:
    DetectorConstruction* mydet;
    G4double dtReal; // real clock time elapsed for progress bar
    G4double fEdep1,   fEdep2;
    G4double fWeight1, fWeight2;
 
    G4double fEparticle; // energy of created particle
    G4String sParticleName; // chosen particle name to histogram
    G4double fTime0; 
    G4Timer * myTimer;
    G4double timeSoFar;
    G4int numberOfBeams; // number of beams to run
    G4String niceTime(G4int tsec);
    G4int eventID;
    G4String sPrimary; // name of primary particle
    G4double  myRes; // resolution set by user
    G4double myTail; // tail set by user
    G4double myPTail; // tail set by user
    G4double myTail2; // tail set by user
    G4double myPTail2; // tail set by user
    G4double myTailH; // tail set by user
    G4double myPTailH; // tail set by user

    G4int iParticleCount;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
