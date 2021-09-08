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

#include "Run.hh"
#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4Timer.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"  // define pi, kBotzmann etc.
#include <time.h>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
    :G4UserEventAction(),
    fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.), myTimer(0),timeSoFar(0.0),numberOfBeams(2),
    eventID(0), fTime0(-1 * s),mydet(det) //time in radioactivity world

{
    myTimer = new G4Timer(); // create a timer to track wall clock time for the program
    myTimer->Start(); // note, timer must be stopped before reading it. restarting rezeroes.
    dtReal = 0.0; // for elapsed time
 
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
    delete myTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String EventAction::niceTime(G4int tsec)
{ // returns time as hh:mm:ss //
    G4String temp;
    div_t result;
    temp = "Hi";
    result = div(tsec, 3600);
    G4int hr = result.quot;
    result= div(result.rem, 60);
    G4int min = result.quot;
    G4int sec = result.rem;
    G4String s_hr = std::to_string(hr);
    if (s_hr.length() < 2) { s_hr = "0" + s_hr; }
    G4String s_min = std::to_string(min);
    if (s_min.length() < 2) { s_min = "0" + s_min; }
    G4String s_sec = std::to_string(sec);
    if (s_sec.length() < 2) { s_sec = "0" + s_sec; }
    return s_hr + ":" + s_min + ":" + s_sec;
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
    fTime0 = -1 * s;
    
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

void EventAction::EndOfEventAction(const G4Event* event)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Apply finite energy resulution - RPF 2021_03_01

    G4double Temp = 0.12 * kelvin; // Temperature. Could set in material definition.
    G4double Cv = 4.3E-4 * joule / kg / kelvin; // Specific heat capacity for gold (Sn is 300 time smaller). Could set for material and temperature
    G4double detrho = mydet->GetAbsorberMaterial()->GetDensity(); // Get absorber density from detector
    G4double detvol = mydet->GetAbsorberSide() * mydet->GetAbsorberSide() * mydet->GetAbsorberThickness(); // Get absorber volume from Detector
    G4double detmass = detrho * detvol; // mass of absorber for Resolution formula


    G4double Hoover_res = 0.53 * keV; // osberved resolution in LANL/Boulder paper for alpha emitter
    // Hoover resolution around 0.6 keV. That is 7*sqrt(k_Boltzmann * Cv * detmass * Temp * Temp) (at 1.8 x 3.6 x 0.015 mm^3)
    // hardwire that in then add in quadrature expected thermal mass. (Assume most of observed resolution was due to crystals etc.)
    G4double Resolution = sqrt(k_Boltzmann * Cv * detmass * Temp * Temp); // 
    Resolution = sqrt(Hoover_res*Hoover_res + Resolution*Resolution);

    if (eventID == (numberOfBeams-1)) // print energy resolution of final event at the end of the run
    {
        G4cout << "vol/m3: " << detvol / m3 << ", mass/kg: " << detmass / kg << ", Res/eV: " << Resolution / eV << G4endl;
    }

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

    // sum energies of the absorber and chip
    G4double Etot = fEdep1 + fEdep2;
    G4double Wtot = (fWeight1 + fWeight2) / Etot;

    // pulse height in absorber
    //
    if (fEdep1 > 0.) {
        fWeight1 /= fEdep1;
        analysisManager->FillH1(0, fEdep1, fWeight1);
        analysisManager->FillH1(9, fEdep1, fWeight1);
    }

    // pulse height in chip
    //   
    if (fEdep2 > 0.) {
        fWeight2 /= fEdep2;
        analysisManager->FillH1(1, fEdep2, fWeight2);
    }

    // total
    //
    analysisManager->FillH1(2, Etot, Wtot);

    // threshold in absorber and chip        
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

    // =========== PERIODIC PRINTING AND SAVING =============== //

    eventID = event->GetEventID();
    if (numberOfBeams < 2) { numberOfBeams = 2; } // avoid div(0) error
    G4int previousProgress = (eventID) * 100 / (numberOfBeams); //previous progress -1 because getEvenID goes from 0-9 and not from 1-10
    G4int currentProgress = (eventID+1) * 100 / (numberOfBeams); //current progress

    //first case when there is no previous event- print line can be removed if user does not want to print 0%
    if (eventID == 0) {
        numberOfBeams = G4RunManager::GetRunManager()->GetNonConstCurrentRun()->GetNumberOfEventToBeProcessed();
        
        sPrimary = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName();
        G4cout << "Primary: " << sPrimary<< G4endl;
        

        myTimer->Stop();
        dtReal = myTimer->GetRealElapsed();
        G4String timestr;
        timestr = myTimer->GetClockTime();
        if (!timestr.empty()) {
            timestr.resize(timestr.size() - 1); // remove end-of-line char
        }
        G4cout << "Begin beamOn: " << numberOfBeams << "\tat " << timestr << " UTC" << G4endl;
        G4cout << "Progress\tTime (UTC)\t\t\tdt\t\tt so far\tt remaining" << G4endl;

        myTimer->Start();

    }

    // periodic status printing. The "&&" takes care of rounding issues
    else if ((currentProgress % 10 == 0 || currentProgress == 1 || currentProgress == 5) && currentProgress != previousProgress)
    {
        G4cout << currentProgress << "%\t"; //line printing the progress bar on the command line
        myTimer->Stop();
        G4double dt = myTimer->GetRealElapsed();
        G4String timestr;
        timestr = myTimer->GetClockTime();
        if (!timestr.empty()) {
            timestr.resize(timestr.size() - 1); // remove end-of-line char
        }
        timeSoFar += dt;
        G4double dt_remaining = dt * (100 - 1.0 * currentProgress) / (currentProgress * 1.0 - previousProgress * 1.0);

        G4cout << currentProgress << "\t" << timestr << "\t" << niceTime(dt) << "\t" << niceTime(timeSoFar) << "\t" << niceTime(dt_remaining) << G4endl;
        myTimer->Start();

        WriteAnAscii(); // write output so far

    }
}


void EventAction::WriteAnAscii()
{
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    for (G4int i = 0; i < analysis->GetNofH1s(); ++i)
    {
        if (analysis->GetH1Activation(i) == 1) // write output for each active hisogram
        {
            G4String myName = analysis->GetFileName() + "_h1_" + std::to_string(i) + ".out"; // set filename based on Analysis filename. Can be set in Macro.
            tools::histo::h1d* myh1d; 
            myh1d = analysis->GetH1(i); // get ith histogram
  
            std::ofstream output; // write headder
            output.open(myName);
            output << "#title\t" << myh1d->get_title() << G4endl;
            output << "#partle\t" << sPrimary<< G4endl;
            output << "#absorb\t"<< mydet->GetAbsorberMaterial()->GetName() << G4endl;
            output << "#Abs s\t" << mydet->GetAbsorberSide() << G4endl;
            output << "#Abs t\t" << mydet->GetAbsorberThickness() << G4endl;
            output << "#Act s\t" << mydet->GetActivitySide() << G4endl;
            output << "#Act t\t" << mydet->GetActivityThickness() << G4endl;
    

            // Write the number of events so far and the total number of beams set. Use scientific notation
            
            std::ios_base::fmtflags oldflags = output.flags(); // original formatting
            output << "#n_done\t" << std::scientific << 1.0*(myh1d->all_entries()) << G4endl;
            output << "#n_set\t" << 1.0*numberOfBeams << G4endl;
    
           // WRITE SPECTRUM 

            output << "#nbins\t" << std::defaultfloat << myh1d->get_bins() << G4endl;
            output << G4endl << "#E_l (MeV)\tcounts" << G4endl;

            G4double Eunit = analysis->GetH1Unit(i); // get unit (1 = MeV; 0.001 = keV, etc.)
            
            for (G4int j = 0; j < myh1d->get_bins(); ++j) // write spectrum
            {
                
                G4double energy =  myh1d->get_axis(0).bin_lower_edge(j);
                G4int counts = myh1d->bins_entries().at(j);

                if (j == 0) // 0 energy for counts below lowest channel
                {
                    energy = 0.0; // catchall for energies less than minimum bin energy

                }
                else if (j < (myh1d->get_bins() - 1)) // spctrum
                {
                    energy = myh1d->get_axis(0).bin_lower_edge(j - 1);
                }
                else // overflow channel
                {
                    energy =  myh1d->get_axis(0).bin_upper_edge(j - 2);
                }
                output << Eunit * energy << "\t" << counts << G4endl;
            }

            output.close();
       

        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
