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
#include <math.h>
#include <header/config.hpp>
#include <header/constants.hpp>
#include <header/globals.hpp>
#include <header/particleDistribution.hpp>
#include <header/boxSize.hpp>
#include <header/greenTable.hpp>
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
double *cof01,*cof11,*cof02,*cof12,*cof03,
*cof13,*cof04,*cof14,*cof05,*cof15,*cof06,*cof16,*cof07,
*cof17,*cof08,*cof18;
const int nGreen = 10000;

shearratei = configuration.getDynincrshrate(); //! shear-rate -> Arquivo de configuracao
per = (4.0 / 3.0) * configuration.getBrownianpecletnum(); //! rotational Peclet number -> Arquivo de configuracao

phi = configuration.getVolumefracpart();
numParticles = configuration.getNumpart();
numRealizations = configuration.getNumreal();
totalRealParticle = numParticles * numRealizations;

if(configuration.getStatanalysis()){
    npast = 2;
}else{
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
particleDistribution(configuration.getMonopolidisp(), numRealizations, numParticles, diam, betaVec);

// Calculating the size of the simulation box based on the 
// number of particles and on the volume fraction defined
// by the user in the simconfig.dat

boxSize(numParticles, configuration.getVolumefracpart(), configuration.getBoxaspectratio(), configuration.getInitialspheraggr());

double qsi = 1.0 * (pow(M_PI,0.5)  / (pow(l * l * h,(1.0/3.0))));


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
 
// periodicStructure();
}

return 0;
}