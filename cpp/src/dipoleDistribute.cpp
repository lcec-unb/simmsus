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


/***************************************************************
* Subroutine responsible for creating an initial distribution 
* of the dipole moments of all particles in all realizations
****************************************************************/
#include <header/globals.hpp>
#include <header/randomic.hpp>
#include <math.h>
#include <iostream>
#include <iomanip>


void distributeDipole(bool orderedDipoles, double percentNonMagPart, bool mixMagNonMagPart){

int e, i, j; 
double modip;
int total = 3 * numParticles * numRealizations;
double *nr = new double[total]{};
int *nr0Index = new int[total]{};
int *nr1Index = new int[total]{};
int *nr2Index = new int[total]{};
randomic(-1.0,1.0,(total),nr);
radomicAccess(0,total - 1,(total),nr0Index);
radomicAccess(0,total - 1,(total),nr1Index);
radomicAccess(0,total - 1,(total),nr2Index);

// If dipoles are distributed in an ordered way

if(orderedDipoles) {

    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            DI0[j * numParticles + i] = 0.0;
            DI1[j * numParticles + i] = 1.0;
            DI2[j * numParticles + i] = 0.0;
        }
    }
}
else{
// For a random dipole distribution
e = percentNonMagPart * numParticles;

// If we are mixing magnetic particles with non magnetic ones...
if(mixMagNonMagPart){
    for(j = 0; j < numRealizations; j++){
        for(i = e + 1; i < numParticles; i++){ 
            DI0[j * numParticles + i] = 0.0;
            DI1[j * numParticles + i] = 0.0;
            DI2[j * numParticles + i] = 0.0;  
        }
    }

    for(j = 0; j < numRealizations; j++){
        for(i = 0; i = e; i++){ 
            DI0[j * e + i] = nr[nr0Index[j * numParticles + i]]; 
            DI1[j * e + i] = nr[nr1Index[j * numParticles + i]];
            DI2[j * e + i] = nr[nr2Index[j * numParticles + i]];            
        }
    }
 

// Normalizing the vectors
    for(j = 0; j < numRealizations; j++){
        for(i = e + 1; i < numParticles; i++){ 
            DI0[j * numParticles + i] = 0.0;
            DI1[j * numParticles + i] = 0.0;
            DI2[j * numParticles + i] = 0.0;
        }
    }

    for(j = 0; j < numRealizations; j++){
        for(i = 0; i = e ; i++){  
            modip = pow(pow(DI0[j * e + i],2.0) + (pow(DI1[j * e + i],2.0)) + (pow(DI2[j * e + i],2.0)),0.5);
            DI0[j * e + i] /= modip;
            DI1[j * e + i] /= modip;
            DI2[j * e + i] /= modip;
        }
    }
}
else{
// If all the particles are magnetic particles, {...
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            DI0[j * numParticles + i] = nr[nr0Index[j * numParticles + i]]; 
            DI1[j * numParticles + i] = nr[nr1Index[j * numParticles + i]];
            DI2[j * numParticles + i] = nr[nr2Index[j * numParticles + i]];            
        }
    }
//    std::cout << "Chegou aqui!" << "\n";
   // Normalizing the vectors
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            modip = pow(pow(DI0[j * numParticles + i],2.0)+ pow(DI1[j * numParticles + i],2.0)+ pow(DI2[j * numParticles + i],2.0),0.5);
            DI0[j * numParticles + i] /= modip;
            DI1[j * numParticles + i] /= modip;
            DI2[j * numParticles + i] /= modip;            
        }
    }
    }
}

}
