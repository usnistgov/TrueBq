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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{

    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);


    fParticleGun->SetParticleEnergy(0 * eV);

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));

    fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
        G4int Z = 10, A = 24;
        G4double ionCharge = 0. * eplus;
        G4double excitEnergy = 0. * keV;

        G4ParticleDefinition* ion
            = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }

    // initialize ranges for random theta DIRECTION (0 is along z axis)

    G4double dirThetaMin = 0.0;				// Min angle (0 degrees)
    G4double dirThetaMax = CLHEP::pi;		// Max angle (pi = 180 degrees)

  // initialize ranges for random phi direction (2*pi)
    G4double dirPhiMin = 0.0;				// 0
    G4double dirPhiMax = CLHEP::twopi;		// 2*pi

    G4double b = cos(dirThetaMax);						// for Random angles
    G4double a = cos(dirThetaMin) - b;					// for Random angles
    G4double cosTheta = a * G4UniformRand() + b;		// for 0 to 180 deg, cos form 1 to -1, b = 2, a = -1)

    G4double sinTheta2 = 1. - cosTheta * cosTheta;		// (sintheta)^2
    if (sinTheta2 < 0.)  sinTheta2 = 0.;				// fix rounding error before sqrt!
    G4double sinTheta = std::sqrt(sinTheta2);

    G4double phi = G4RandFlat::shoot(dirPhiMin, dirPhiMax);  // random  phi from from 0 to 2*pi

    G4ThreeVector direction = G4ThreeVector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta).unit(); // create direction vector


    fParticleGun->SetParticleMomentumDirection(direction);
    //create vertex
    //   
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
