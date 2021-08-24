# TrueBq01
TrueBq01 is a simple Geant4 simulation for cryogenic decay-energy spectrometery.
It has been compiled using Geant4.10.6 on Windows using Microsoft Visual Studio Community 2019
Version 16.8.3. 

TrueBq01 is based on the the rdecay02 example included with Geant4 distribtuions. See rdecay02 README file within file structure

## Installation Process
1. Install Geant4
2. Download source code folder (TrueBq01)
3. Follow Geant4 user code compilation instructions (I used Cmake version 3.19.2)
4. Copy addtional User files from source code folder to installation folder

## Simulation information

### Geometry:
 1. Rectangular solid absorber (gold as default material) sitting on a 
 2. Rectangualr solid chip (silicon as default)
 
Dimensions of both solids can be set by user (e.g. in batch.mac)

### Source:
Source is selected by the user; can be radionuclide or particle of a given energy.
Source is either centered or uniform in the absorber. Other options can be implemented in the code.

### Physics:
BiasedRDPhysics with Radioactivation

### Event processing:
Decay Energy Spectrum of total energy deposited in the absorber and chip are tallied in separate histograms.
Energy resoultion added as Gaussian noise. Rough resolution is based on Hoover (2015) results from LANL & NIST, scaled for mass and temperature

### Output:
Decay Energy Spectrum of total energy deposited in the absorber is saved in an ascii (csv) file. Can be changed to root in histomanager header file.
For ascii file, header gives number of bins and energy range. The column 1 is energy, column 2 counts, then other moments (anyone see where this is documenented?)

## Contact Information
Ryan Fitzgerald
ryan.fitzgerald@nist.gov

## Acknowledgement for reused code if any
See LICENCE.md

Also see this from rdecay02 source file:
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

