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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "G4Timer.hh" 
#include "Run.hh"
#include "HistoManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"  // define pi, kBotzmann etc.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
    :G4UserEventAction(),
    fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
    fTime0(-1 * s),//time in radioactivity world
    timer()
    
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
    fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
    fTime0 = -1 * s;
    timer.Start(); //initializing the G4Timer variable

    
      //  G4cout<< G4Timer::GetClockTime
   
    //print the total number of beams that being run at the beginning.
    //print the start time. 
    //Ex: we are print 10000 beams at 14.23s
    //start time at each 10% 
    //print the radionuclide that being run at the beginning

    
 


   
    //to hold the number of beams that the user decided to run
    G4int numberOfBeams = G4RunManager::GetRunManager()->GetNonConstCurrentRun()->GetNumberOfEventToBeProcessed();

    
   
    G4int previousProgress = (event->GetEventID() - 1) * 100 / (numberOfBeams - 1); //previous progress -1 because getEvenID goes from 0-9 and not from 1-10
    G4int currentProgress = event->GetEventID() * 100 / (numberOfBeams - 1); //current progress
   
    //first case when there is no previous event- print line can be removed if user does not want to print 0%
    if (event->GetEventID() == 0) {
        G4cout << "Number of Beams: " << numberOfBeams << " beams" << G4endl; //printing the number of beams that are being run
        G4cout << "Radionuclide: " << event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName() << G4endl;
        
        G4cout << currentProgress << "% -- " ; //line printing the 0% progress bar

        G4cout << "Run Started on :" << timer.GetClockTime() << G4endl; // print the time when it started to run those beams
        
    } 
   // second case when there is a repetion due to the rounding
    else if (currentProgress % 10 == 0 && currentProgress!=previousProgress) {
       G4cout << currentProgress << "% --" ; //line printing the progress bar on the command line
       G4cout << " Run Started on :" << timer.GetClockTime()<<G4endl;  //Bprint the time when it started to run those beams

      // if (currentProgress == 100) {
         //  timer.Stop();
         //  G4cout << timer.GetRealElapsed() << G4endl; 
      // }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4int iVol, G4double edep,
    G4double time, G4double weight)
{
    // initialize t0
    if (fTime0 < 0.) fTime0 = time;

    // out of time window ?
    const G4double TimeWindow(1 * microsecond);
    if (std::fabs(time - fTime0) > TimeWindow) return;

    if (iVol == 1) { fEdep1 += edep; fWeight1 += edep * weight; }
    if (iVol == 2) { fEdep2 += edep; fWeight2 += edep * weight; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Apply finite energy resulution - RPF 2021_03_01

    G4double Temp = 0.12 * kelvin; // Temperature 
    G4double Cv = 4.3E-4 * joule / kg / kelvin; // Specific heat capacity for gold (Sn is 300 time smaller). Could set for material and temperature
    G4double detrho = 19300. * kg / m3; // density of gold
    G4double detvol = pi * pow(2, 2 * mm) * 0.2 * mm; // 2 mm radius x 0.2 mm height gold
    detvol = 1.8 * mm * 3.6 * mm * 0.015 * mm;; // Hoover 2015 geometry
    G4double detmass = detrho * detvol; // mass of target for Resolution formula

    G4double Resolution = sqrt(k_Boltzmann * Cv * detmass * Temp * Temp) * 7; // 

    // G4cout << "vol/m3: " << detvol / m3 << ", mass/kg: " << detmass / kg << ", Res/eV: " << Resolution / eV << G4endl;

    G4double E1res;
    //
    G4double E2res;

    // Generate random gaussian noise. Distribution is centered on original energy, with distribution set by thermal noise
    E1res = G4RandGauss::shoot(fEdep1, Resolution);
    E2res = G4RandGauss::shoot(fEdep2, Resolution);

    // G4cout << fEdep1 << " " << E1res << G4endl;

    fEdep1 = E1res; // replace event energy with the modified value
    fEdep2 = E2res; // replace event energy with the modified value

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

    // sum energies of the target and detector
    G4double Etot = fEdep1 + fEdep2;
    G4double Wtot = (fWeight1 + fWeight2) / Etot;

    // pulse height in target
    //
    if (fEdep1 > 0.) {
        fWeight1 /= fEdep1;
        analysisManager->FillH1(0, fEdep1, fWeight1);
        analysisManager->FillH1(9, fEdep1, fWeight1);
    }

    // pulse height in detector
    //   
    if (fEdep2 > 0.) {
        fWeight2 /= fEdep2;
        analysisManager->FillH1(1, fEdep2, fWeight2);
    }

    // total
    //
    analysisManager->FillH1(2, Etot, Wtot);

    // threshold in target and detector        
    const G4double Threshold1(10 * keV), Threshold2(10 * keV);

    //coincidence, anti-coincidences 
    //  
    G4bool coincidence = ((fEdep1 >= Threshold1) && (fEdep2 >= Threshold2));
    G4bool anti_coincidence1 = ((fEdep1 >= Threshold1) && (fEdep2 < Threshold2));
    G4bool anti_coincidence2 = ((fEdep1 < Threshold1) && (fEdep2 >= Threshold2));

    if (coincidence)       analysisManager->FillH1(3, fEdep2, fWeight2);
    if (anti_coincidence1) analysisManager->FillH1(4, fEdep1, fWeight1);
    if (anti_coincidence2) analysisManager->FillH1(5, fEdep2, fWeight2);

    // pass energies to Run
    //  
    Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    run->AddEdep(fEdep1, fEdep2);

    //number of event to run
    //minimum amount from printing
}


        

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

