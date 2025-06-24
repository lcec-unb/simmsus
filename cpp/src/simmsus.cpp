/*
 * This program is licensed granted by the University of Brasília - UnB (the “University”)
 * for use of SIMULATION OF MAGNETIC SUSPENSIONS - SIMMSUS (“the Software”) through this website
 * https://github.com/rafaelgabler/simmsus (the ”Website”).
 *
 * By downloading the Software through the Website, you (the “Licensee”) are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about SIMMSUS please contact: rafael.gabler@unb.br (Rafael Gabler Gontijo)
 * or leandro.zanotto@gmail.com (Leandro Zanotto).
 *
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <header/config.hpp> //Migrar para TomL
#include <header/constants.hpp>
#include <header/globals.hpp>
#include <header/particleDistribution.hpp>
#include <header/boxSize.hpp>
#include <header/greenTable.hpp>
#include <header/periodicStructure.hpp>
#include <header/initialCondition.hpp>
#include <header/dipoleDistribute.hpp>

// #include <brownian.hpp>
// #include <repulsion.hpp>
// #include <periodicStructure.hpp>
using namespace std;

int main(int argc, char **argv){

if (argc < 2){
    cout << "File Not Found - Please specify an input configuration file as paramenter\n";
}else{
//Read the Input File and Create the Object

Configuration configuration(argv[1]);
numParticles = configuration.getNumpart();
numRealizations = configuration.getNumreal();

//Globals Initializations
X0 = new double[numRealizations * numParticles]{};
X1 = new double[numRealizations * numParticles]{};
X2 = new double[numRealizations * numParticles]{};
U0 = new double[numRealizations * numParticles]{};
U1 = new double[numRealizations * numParticles]{};
U2 = new double[numRealizations * numParticles]{};
W0 = new double[numRealizations * numParticles]{};
W1 = new double[numRealizations * numParticles]{};
W2 = new double[numRealizations * numParticles]{};
XI0 = new double[nb * numRealizations * numParticles]{};
XI1 = new double[nb * numRealizations * numParticles]{};
XI2 = new double[nb * numRealizations * numParticles]{};
ILF0 = new double[125]{};
ILF1 = new double[125]{};
ILF2 = new double[125]{};
ILR0 = new double[27]{};
ILR1 = new double[27]{};
ILR2 = new double[27]{};
DI0 = new double[numRealizations * numParticles]{};
DI1 = new double[numRealizations * numParticles]{};
DI2 = new double[numRealizations * numParticles]{};
FORCAS10 = new double[numRealizations * numParticles]{};
FORCAS11 = new double[numRealizations * numParticles]{};
FORCAS12 = new double[numRealizations * numParticles]{};
FORCAS20 = new double[numRealizations * numParticles]{};
FORCAS21 = new double[numRealizations * numParticles]{};
FORCAS22 = new double[numRealizations * numParticles]{};
FORCAS30 = new double[numRealizations * numParticles]{};
FORCAS31 = new double[numRealizations * numParticles]{};
FORCAS32 = new double[numRealizations * numParticles]{};
FORCAS40 = new double[numRealizations * numParticles]{};
FORCAS41 = new double[numRealizations * numParticles]{};
FORCAS42 = new double[numRealizations * numParticles]{};
FORCAS50 = new double[numRealizations * numParticles]{};
FORCAS51 = new double[numRealizations * numParticles]{};
FORCAS52 = new double[numRealizations * numParticles]{};
FORCAS60 = new double[numRealizations * numParticles]{};
FORCAS61 = new double[numRealizations * numParticles]{};
FORCAS62 = new double[numRealizations * numParticles]{};
TORQUES10 = new double[numRealizations * numParticles]{};
TORQUES11 = new double[numRealizations * numParticles]{};
TORQUES12 = new double[numRealizations * numParticles]{};
TORQUES20 = new double[numRealizations * numParticles]{};
TORQUES21 = new double[numRealizations * numParticles]{};
TORQUES22 = new double[numRealizations * numParticles]{};
TORQUES30 = new double[numRealizations * numParticles]{};
TORQUES31 = new double[numRealizations * numParticles]{};
TORQUES32 = new double[numRealizations * numParticles]{};
FT0 = new double[numRealizations * numParticles]{};
FT1 = new double[numRealizations * numParticles]{};
FT2 = new double[numRealizations * numParticles]{};
TT0 = new double[numRealizations * numParticles]{};
TT1 = new double[numRealizations * numParticles]{};
TT2 = new double[numRealizations * numParticles]{};
AUX0 = new double[numRealizations * numParticles]{};
AUX1 = new double[numRealizations * numParticles]{};
AUX2 = new double[numRealizations * numParticles]{};
AUX3 = new double[numRealizations * numParticles]{};
HAUX10 = new double[numParticles]{};
HAUX11 = new double[numParticles]{};
HAUX12 = new double[numParticles]{};
HAUX20 = new double[numParticles]{};
HAUX21 = new double[numParticles]{};
HAUX22 = new double[numParticles]{};
C0 = new double[numRealizations * numParticles]{};
C1 = new double[numRealizations * numParticles]{};
C2 = new double[numRealizations * numParticles]{}; 
H10 = new double[numBoxes]{}; 
H11 = new double[numBoxes]{};
H12 = new double[numBoxes]{};
H20 = new double[numRecBoxes]{};
H21 = new double[numRecBoxes]{};
H22 = new double[numRecBoxes]{};
AUXT0 = new double[numParticles]{};
AUXT1 = new double[numParticles]{};
AUXT2 = new double[numParticles]{};
T10 = new double[numBoxes]{}; 
T11 = new double[numBoxes]{};
T12 = new double[numBoxes]{};
T20 = new double[numRecBoxes]{};
T21 = new double[numRecBoxes]{};
T22 = new double[numRecBoxes]{};
F10 = new double[numBoxes]{}; 
F11 = new double[numBoxes]{};
F12 = new double[numBoxes]{};
F20 = new double[numRecBoxes]{};
F21 = new double[numRecBoxes]{};
F22 = new double[numRecBoxes]{};
AUXF0 = new double[numParticles]{};
AUXF1 = new double[numParticles]{};
AUXF2 = new double[numParticles]{};

shearratei = configuration.getDynincrshrate(); //! shear-rate -> Arquivo de configuracao
per = (4.0 / 3.0) * configuration.getBrownianpecletnum(); //! rotational Peclet number -> Arquivo de configuracao

phi = configuration.getVolumefracpart();

totalRealParticle = numParticles * numRealizations;

npast = 2;
if(!configuration.getStatanalysis()){
    npast = configuration.getSimulationtime()  / configuration.getNumtimestep();
}

// Number of random numbers used in each time-step
nnr = 3 * totalRealParticle;

if(configuration.getAccounthi() || configuration.getPmf() || configuration.getPmt()){
    periodicity = true;
}

diam = new double[totalRealParticle]{};
betaVec = new double[totalRealParticle]{};

// Defining the size of the particles
particleDistribution(configuration.getMonopolidisp(), diam, betaVec);

// Calculating the size of the simulation box based on the 
// number of particles and on the volume fraction defined
// by the user in the simconfig.dat

boxSize(numParticles, configuration.getVolumefracpart(), configuration.getBoxaspectratio(), configuration.getInitialspheraggr());

if(configuration.getMp()){ 

    distributeDipole(configuration.getOrdereddipoles(), configuration.getPercentnonmagpart(), configuration.getMixmagnonmagpart());    
}

initialCondition(configuration.getInitialspheraggr());

double qsi = pow(M_PI,0.5)  / (pow(l * l * h,(1.0/3.0)));

// Building a table with all the Green-functions
// And building the periodic structure to compute
// the Ewald summations

if(periodicity){
    double wave = pow(nbr,(1.0/3.0));
    if(wave == 5.0){
        if(configuration.getAccounthi()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            greenTableLigaihWave5(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13);
        }else if(configuration.getPmt()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof05 = new double[nGreen];
            cof15 = new double[nGreen];
            cof07 = new double[nGreen];
            cof17 = new double[nGreen];
            greenTableTmagper5(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13,
            cof04, cof14, cof05, cof15, cof07, cof17);
        }else if(configuration.getPmf()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof06 = new double[nGreen];
            cof16 = new double[nGreen];
            cof08 = new double[nGreen];
            cof18 = new double[nGreen];
            greenTableFmagper3(nGreen, qsi, l, nb, wave, h, cof01,  cof11, cof02, cof12, cof03, cof13,
            cof06, cof16, cof08, cof18);
        }
    }else if(wave == 3.0){
        if(configuration.getAccounthi()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            greenTableLigaihWave3(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13);
        }else if(configuration.getPmt()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof05 = new double[nGreen];
            cof15 = new double[nGreen];
            cof07 = new double[nGreen];
            cof17 = new double[nGreen];
            greenTableTmagper3(nGreen, qsi, l, nb, wave, h, cof01,  cof11, cof02, cof12, cof03, cof13,
            cof04, cof14, cof05, cof15, cof07, cof17);
        }else if(configuration.getPmf()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof06 = new double[nGreen];
            cof16 = new double[nGreen];
            cof08 = new double[nGreen];
            cof18 = new double[nGreen];
            greenTableFmagper3(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13,
            cof06, cof16, cof08, cof18);            

        }
    }
}
periodicStructure();
cout << numParticles << "\n";
for(int j = 0; j < 1; j++){    
    for(int i = 0; i < numParticles; i++){        
        cout << X0[j * numParticles + i] << "\t" << X1[j * numParticles + i] << "\t" << X2[j * numParticles + i] << "\t";
        cout << DI0[j * numParticles + i] << "\t" << DI1[j * numParticles + i] << "\t" << DI2[j * numParticles + i] << "\t" << diam[j * numParticles + i];
        cout << "\n";
    }    
}

//free memory

delete[] X0;
delete[] X1;
delete[] X2;
delete[] U0;
delete[] U1;
delete[] U2;
delete[] W0;
delete[] W1;
delete[] W2;
delete[] XI0;
delete[] XI1;
delete[] XI2;
delete[] ILF0;
delete[] ILF1;
delete[] ILF2;
delete[] ILR0;
delete[] ILR1;
delete[] ILR2;
delete[] DI0;
delete[] DI1;
delete[] DI2;
delete[] FORCAS10;
delete[] FORCAS11;
delete[] FORCAS12;
delete[] FORCAS20;
delete[] FORCAS21;
delete[] FORCAS22;
delete[] FORCAS30;
delete[] FORCAS31;
delete[] FORCAS32;
delete[] FORCAS40;
delete[] FORCAS41;
delete[] FORCAS42;
delete[] FORCAS50;
delete[] FORCAS51;
delete[] FORCAS52;
delete[] FORCAS60;
delete[] FORCAS61;
delete[] FORCAS62;
delete[] TORQUES10;
delete[] TORQUES11;
delete[] TORQUES12;
delete[] TORQUES20;
delete[] TORQUES21;
delete[] TORQUES22;
delete[] TORQUES30;
delete[] TORQUES31;
delete[] TORQUES32;
delete[] FT0;
delete[] FT1;
delete[] FT2;
delete[] TT0;
delete[] TT1;
delete[] TT2;
delete[] AUX0;
delete[] AUX1;
delete[] AUX2;
delete[] AUX3;
delete[] HAUX10;
delete[] HAUX11;
delete[] HAUX12;
delete[] HAUX20;
delete[] HAUX21;
delete[] HAUX22;
delete[] C0;
delete[] C1;
delete[] C2;
delete[] H10; 
delete[] H11;
delete[] H12;
delete[] H20;
delete[] H21;
delete[] H22;
delete[] AUXT0;
delete[] AUXT1;
delete[] AUXT2;
delete[] F10;
delete[] F11;
delete[] F12;
delete[] F20;
delete[] F21;
delete[] F22;
delete[] AUXF0;
delete[] AUXF1;
delete[] AUXF2;

if(cof01 != NULL)
    delete[] cof01;
if(cof11 != NULL)
    delete[] cof11;
if(cof02 != NULL)
    delete[] cof02;
if(cof12 != NULL)
    delete[] cof12;
if(cof03 != NULL)
    delete[] cof03;
if(cof13 != NULL)
    delete[] cof13;
if(cof04 != NULL)
    delete[] cof04;
if(cof14 != NULL)
    delete[] cof14;
if(cof05 != NULL)
    delete[] cof05;
if(cof15 != NULL)
    delete[] cof15;
if(cof06 != NULL)
    delete[] cof06;
if(cof16 != NULL)
    delete[] cof16;
if(cof07 != NULL)
    delete[] cof07;
if(cof17 != NULL)
    delete[] cof17;
if(cof08 != NULL)
    delete[] cof08;
if(cof18 != NULL)
    delete[] cof18;

}
return 0;
}