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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//


#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh" // second-most accurate
#include "G4EmStandardPhysics_option4.hh" // most accurate

#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4NuclideTable.hh"
#include "BiasedRDPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4IonCoulombScatteringModel.hh"


#include "G4EmStandardPhysicsSS.hh"
#include "G4EmPenelopePhysics.hh"

#include "PhysListEmStandardNR.hh" // local source code copied from example TestEM7
#include "G4EmProcessOptions.hh"
// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4RegionStore.hh"

#include "G4StepLimiterPhysics.hh"


PhysicsList::PhysicsList()
	:G4VModularPhysicsList()
{
	G4int verb = 1;
	SetVerboseLevel(verb);

	//add new units for radioActive decays
	//
	new G4UnitDefinition("millielectronVolt", "meV", "Energy", 1.e-3 * eV);
	// 
	const G4double minute = 60 * second;
	const G4double hour = 60 * minute;
	const G4double day = 24 * hour;
	const G4double year = 365 * day;
	new G4UnitDefinition("minute", "min", "Time", minute);
	new G4UnitDefinition("hour", "h", "Time", hour);
	new G4UnitDefinition("day", "d", "Time", day);
	new G4UnitDefinition("year", "y", "Time", year);

	// Mandatory for G4NuclideTable
	// Half-life threshold must be set small or many short-lived isomers 
	// will not be assigned life times (default to 0) 
	G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond); // 0.1*picosecond; 100*minutes to avoid Pu-239 extra gammas
//	 G4NuclideTable::GetInstance()->SetLevelTolerance(1.0 * eV);

// EM physics: CHOOSE YOUR PHYSICS LIST
//	RegisterPhysics(new G4EmStandardPhysics()); // Standard Physics
	RegisterPhysics(new G4EmStandardPhysics_option4()); // Standard Physics option 4 = most accurate
//	RegisterPhysics(new G4EmStandardPhysicsSS()); // Single Scattering, includes G4IonCoulombScatteringModel
//	RegisterPhysics(new G4EmPenelopePhysics); // Penelope model
//	RegisterPhysics(new PhysListEmStandardNR()); // G4NuclearRecoil

// Step limiter
	RegisterPhysics(new G4StepLimiterPhysics()); // RPF

	G4EmParameters* param = G4EmParameters::Instance();
	param->SetAugerCascade(true);
	param->SetStepFunction(.1, 1 * CLHEP::um);
	param->SetStepFunctionMuHad(.1, 1 * CLHEP::um);
	param->SetMinEnergy(0.4 * eV);
	param->SetLowestElectronEnergy(100 * eV);
	
	 G4EmProcessOptions emOptions;
	
	//physics tables
	//
	emOptions.SetMinEnergy(.04 * eV);        // I changed this to 1 eV (from 10 eV) to get events down to 200 eV
	
	//emOptions.SetMaxEnergy(10 * MeV);      // was 10*TeV
	//.SetDEDXBinning(12 * 20);
	//emOptions.SetLambdaBinning(12 * 20);


	// Decay
	RegisterPhysics(new G4DecayPhysics());

	// Radioactive decay
	RegisterPhysics(new BiasedRDPhysics());

	// Hadron Elastic scattering
	RegisterPhysics(new G4HadronElasticPhysics(verb));

	// Hadron Inelastic physics
	RegisterPhysics(new G4HadronPhysicsFTFP_BERT(verb));
	////RegisterPhysics( new G4HadronInelasticQBBC(verb));        
	////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));

	// Ion Elastic scattering
	RegisterPhysics(new G4IonElasticPhysics(verb));

	// Ion Inelastic physics
	RegisterPhysics(new G4IonPhysics(verb));
	////RegisterPhysics( new G4IonINCLXXPhysics(verb));

	// Gamma-Nuclear Physics
	G4EmExtraPhysics* gnuc = new G4EmExtraPhysics(verb);
	gnuc->ElectroNuclear(false); 
	gnuc->MuonNuclear(false);
	
	RegisterPhysics(gnuc);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	G4BosonConstructor  pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts() // set general cut value (range of 2ndary particle to produce). Then set shorter cut (lower E) for Absorber
// in Gold, 5 um -> 53 keV e' and 5.2 keV gamma, 500 eV 197Au in Single Scattering mode (EmStandardSS)
// in Si, 1 mm -> 548 keV e' and 7 keV gamma
// To do: Add these options to messenger
{
	SetCutValue(100 * um, "proton");
	SetCutValue(1 * um, "e-");
	SetCutValue(1 * um, "e+");
	SetCutValue(10 * um, "gamma");

	// to set cut below 990 eV, in the Marco use /cuts/setLowEdge 100 eV, or use this code below //
	// to do, work this into a custom cuts macro that does regular, regional, and lowedge cuts!
	// can use /run/setCut or run/setCutForAGivenParticle, run/setCutForARegion as well, down to 990 eV
	// 
	// G4double HighEdge = G4ProductionCutsTable::GetProductionCutsTable()->GetHighEdgeEnergy();
	// ::GetProductionCutsTable()->SetEnergyRange(100 * eV, HighEdge);


	// Production thresholds for ABSORBER region
	G4Region* region;
	G4String regName;
	G4ProductionCuts* cuts;
	regName = "Absorber";
	region = G4RegionStore::GetInstance()->GetRegion(regName);
	cuts = new G4ProductionCuts;
	
	cuts->SetProductionCut(1 * um); // general. Includes recoil ions, though can't seem to set them per se. (0.1 um)
	//cuts->SetProductionCut(.25 * um,"e-"); // below 1 um, limit is still 990 eV. Why? (0.1 um)
	//cuts->SetProductionCut(.25 * um, "e+"); // 
	//cuts->SetProductionCut(.25 * um, "gamma"); // 
	region->SetProductionCuts(cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......