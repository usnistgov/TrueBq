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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
    : fFileName("TrueBq01")
{
    Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
    // Create or get analysis manager
    // The choice of analysis technology is done via selection of a namespace
    // in HistoManager.hh
    //
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    analysis->SetFileName(fFileName);
    analysis->SetVerboseLevel(1);
    analysis->SetActivation(true);     //enable inactivation of histos, nTuples

    // Default values (to be reset via /analysis/h1/set command)               
    G4int nbins = 100;
    G4double vmin = 0.;
    G4double vmax = 100.;

    // Create all histograms as inactivated 
    // as we have not yet set nbins, vmin, vmax
    //
    ////analysis->SetHistoDirectoryName("histo");  
    ////analysis->SetFirstHistoId(1);

    G4int id = analysis->CreateH1("H10", "Energy_deposit_(MeV)_in_the_absorber",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);


    id = analysis->CreateH1("H11", "Energy_deposit_(MeV) in the Detector",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H12", "Energy_deposit_(MeV) in absorber and Detector",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H13",
        "Coincidence_spectrum_(MeV)_between the absorber and Detector",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H14",
        "Anti-coincidence_spectrum_(MeV) in the traget",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H15",
        "Anti-coincidence_spectrum_(MeV) in the Detector",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H16", "Decay_emission_spectrum_(0 - 10 MeV)",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H17", "Decay_emission_spectrum_(0 - 1 MeV)",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H18", "Decay_emission_spectrum_(0 - 0.1 MeV)",
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    id = analysis->CreateH1("H19", "Energy_deposit_(MeV) in the absorber2", // same as 10 thing but can have different # of channels, range
        nbins, vmin, vmax);
    analysis->SetH1Activation(id, false);

    // nTuples
    //
    ////analysis->SetNtupleDirectoryName("ntuple");
    ////analysis->SetFirstNtupleId(1);
    //       
    analysis->CreateNtuple("T1", "Emitted Particles");
    analysis->CreateNtupleDColumn("PID");       //column 0
    analysis->CreateNtupleDColumn("Energy");    //column 1
    analysis->CreateNtupleDColumn("Time");      //column 2
    analysis->CreateNtupleDColumn("Weight");    //column 3
    analysis->FinishNtuple();

    analysis->CreateNtuple("T2", "RadioIsotopes");
    analysis->CreateNtupleDColumn("PID");       //column 0
    analysis->CreateNtupleDColumn("Time");      //column 1
    analysis->CreateNtupleDColumn("Weight");    //column 2
    analysis->FinishNtuple();

    analysis->CreateNtuple("T3", "Energy depositions");
    analysis->CreateNtupleDColumn("Energy");    //column 0
    analysis->CreateNtupleDColumn("Time");      //column 1
    analysis->CreateNtupleDColumn("Weight");    //column 2
    analysis->FinishNtuple();

    analysis->CreateNtuple("RDecayProducts", "All Products of RDecay");
    analysis->CreateNtupleDColumn("PID");       //column 0
    analysis->CreateNtupleDColumn("Z");         //column 1
    analysis->CreateNtupleDColumn("A");         //column 2    
    analysis->CreateNtupleDColumn("Energy");    //column 3
    analysis->CreateNtupleDColumn("Time");      //column 4
    analysis->CreateNtupleDColumn("Weight");    //column 5
    analysis->FinishNtuple();

    analysis->SetNtupleActivation(false);
}

void HistoManager::WriteAnAscii()
{
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    analysis->GetFirstH1Id();
    analysis->GetFileName();
    analysis->GetActivation();
    analysis->GetH1Activation(0);
    analysis->GetH1Ascii(0);
    analysis->GetNofH1s();
    G4cout << "analysis->GetFileName()" << analysis->GetFileName() << G4endl;
    G4cout << "analysis->GetActivation() " << analysis->GetActivation() << G4endl;
    G4cout << "analysis->GetNofH1s() " << analysis->GetNofH1s() << G4endl;
    G4cout << "analysis->GetNofH1s() " << analysis->GetNofH1s() << G4endl;
    G4cout << "analysis->GetNofH1s() " << analysis->GetNofH1s() << G4endl;

    G4cout << "analysis->GetNofH1s() " << analysis->GetNofH1s() << G4endl;
    G4cout << "analysis->GetH1Ascii(0) " << analysis->GetH1Ascii(0) << G4endl;
    G4cout << "analysis->GetFirstH1Id() " << analysis->GetFirstH1Id() << G4endl;
    G4cout << "analysis->GetH1Activation(0) " << analysis->GetH1Activation(0) << G4endl;
    G4cout << "analysis->GetH1Activation(3) " << analysis->GetH1Activation(3) << G4endl;


    tools::histo::h1d* ryan;


    for (G4int i = 0; i < analysis->GetNofH1s(); ++i)
    {
        if (analysis->GetH1Activation(i) == 1) // write output for each active hisogram
        {
            G4String myName = analysis->GetFileName() + "_h1_" + std::to_string(i)+".dat";
            ryan = analysis->GetH1(i);
            ryan->get_title();
            ryan->get_bins();
            ryan->get_axis(i).bin_center(0);
            ryan->bin_height(0);
            G4cout << "myName" << myName << G4endl;
            G4cout << "ryan->get_title();" << ryan->get_title() << G4endl;
            G4cout << "ryan->get_bins()" << ryan->get_bins() << G4endl;

            
            std::ofstream output; // write headder
            output.open(myName);
            output << "title\t" <<ryan->get_title() << G4endl;
            output << "nbins\t" << ryan->get_bins() << G4endl << G4endl;
            output << "E_min (MeV)\tcounts" << G4endl;

            for (G4int j = 0; j < ryan->get_bins(); ++j) // write spectrum
            {
                G4double energy = ryan->get_axis(0).bin_lower_edge(j);
                G4int counts = ryan->bins_entries().at(j);

                if (j == 0) // 0 energy for counts below lowest channel
                {
                    energy = 0.0; // catchall for energies less than minimum bin energy

                }
                else if (j< (ryan->get_bins()-1)) // spctrum
                {
                    energy = ryan->get_axis(0).bin_lower_edge(j-1);
                }
                else // overflow channel
                {
                    energy =  ryan->get_axis(0).bin_upper_edge(j-2);
                }
                output<< energy << "\t"<< counts << G4endl;

            }

            output.close();
            G4cout << ryan->sum_all_bin_heights() << "\t" << ryan->sum_bin_heights() << "\t"<<ryan->extra_entries() << G4endl;
            
        }


    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......