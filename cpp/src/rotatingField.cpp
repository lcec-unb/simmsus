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

// !*************************************************!
// ! 		     SIMMSUS			  ! 
// !SUBROUTINE: rotating_field			  !					         
// !Last update: 16/07/2023			  !
// !*************************************************!

// !*************************************************!
// ! Subroutine resposible for computing magnetic    !
// ! torques due to an external rotating field	  !
// !*************************************************!
#include <iostream>
#include <numeric>
#include <math.h>
#include <header/globals.hpp>

void rotatingField(double alfa, double wt, bool mistura, bool gravidade, bool shear, double brownianmpecletnum){

int auxiliary = 1, i, j;
double percentual;

// ! These subroutine may consider a magnetic fluidized 
// ! bed mixing magnetic and non-magnetic particles, 
// ! the logical variable "mistura" activates this po-
// ! ssibility. If mistura = TRUE then we have to set 
// ! zero dipoles for a certain percentage "percentual" 
// ! of the particles.

if(mistura){
    auxiliary = (percentual* numParticles) + 1;
}

for(j = 0; j < numRealizations; j++){
    for(i = auxiliary; i < numParticles; i++){ 

    // ! Rotating field on the upper side

    if(gravidade) {
        TORQUES20[j * numParticles + i] = (3.0 * alfa / 4.0 * brownianmpecletnum) * (DI1[j * numParticles + i] * cos(wt) - (DI2[j * numParticles + i]*sin(wt)));
        TORQUES21[j * numParticles + i] = (3.0 * alfa / 4.0 * brownianmpecletnum) * (-(DI0[j * numParticles + i] * cos(wt)));
        TORQUES22[j * numParticles + i] = (3.0 * alfa / 4.0 * brownianmpecletnum) * ((DI0[j * numParticles + i] * sin(wt)));
    }
    else{
        if(shear) {
            TORQUES20[j * numParticles + i]=(3.0 * alfa / 4.0 * brownianmpecletnum) * (DI1[j * numParticles + i]* cos(wt) - (DI2[j * numParticles + i]*sin(wt)));
            TORQUES21[j * numParticles + i]=(3.0 * alfa / 4.0 * brownianmpecletnum) * (-(DI0[j * numParticles + i]* cos(wt)));
            TORQUES22[j * numParticles + i]=(3.0 * alfa / 4.0 * brownianmpecletnum) * ((DI0[j * numParticles + i]* sin(wt)));
        }else{
            TORQUES20[j * numParticles + i]=(3.0 * alfa / 4.0) * (DI1[j * numParticles + i] * cos(wt) - (DI2[j * numParticles + i] * sin(wt)));
            TORQUES21[j * numParticles + i]=(3.0 * alfa / 4.0) * (-(DI0[j * numParticles + i] * cos(wt)));
            TORQUES22[j * numParticles + i]=(3.0 * alfa / 4.0) * ((DI0[j * numParticles + i] * sin(wt)));
        }
    }
}
} 

}