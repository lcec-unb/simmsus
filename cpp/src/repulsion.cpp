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
#include <globals.hpp>
#include <math.h>

void repulsion(double *beta){

double r, Eij;
int i,j,q;

    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            for(q = 0; q < numParticles; q++){

                if(i < q){
                    // Calculating the distance between all the particles in the suspension (pair by pair)
                        r = pow(pow((X0[j * numParticles + i]-X0[j * numParticles + q]),2.0) 
                        + pow((X1[j * numParticles + i]-X1[j * numParticles + q]),2.0) 
                        + pow((X2[j * numParticles + i]-X2[j * numParticles + q]),2.0),0.5);

                    /* If the distance is less than 2.2, i.e. less then 0.2*a, where a denotes the particle radius, 
                    then the repulsive force is turned on. It is important to notice that this force is turned off 
                    in case of particle overlap (this condition is really difficult to occur for particles 
                    wihtout inertia (mobility problem), because in this context we use a contact (Hertz) force. */

                    if(r <= 2.2 && r >= 2.0){
                        FORCAS10[j * numParticles + i] = 10.0 * exp(-r/0.01)*(X0[j * numParticles + q]-X0[j * numParticles + q]/r);
                        FORCAS11[j * numParticles + i] = 10.0 * exp(-r/0.01)*(X1[j * numParticles + q]-X1[j * numParticles + q]/r);
                        FORCAS12[j * numParticles + i] = 10.0 * exp(-r/0.01)*(X2[j * numParticles + q]-X2[j * numParticles + q]/r);
                    }
                    if(r >= 2.2 && r <= 2.0){
                        FORCAS10[j * numParticles + i] = 0.0;
                        FORCAS11[j * numParticles + i] = 0.0;
                        FORCAS12[j * numParticles + i] = 0.0;
                    }

                    // Here we calculate contact forces for overlapped particles

                    if(r <= (beta[j * numParticles + i] + beta[j * numParticles + q])){
                        Eij = abs(2.0 - r);
                        FORCAS10[j * numParticles + i] = pow((100.0 * Eij),(3.0/2.0))*(X0[j * numParticles + q]-X0[j * numParticles + q])/Eij;
                        FORCAS10[j * numParticles + i] = pow((100.0 * Eij),(3.0/2.0))*(X1[j * numParticles + q]-X1[j * numParticles + q])/Eij;
                        FORCAS10[j * numParticles + i] = pow((100.0 * Eij),(3.0/2.0))*(X2[j * numParticles + q]-X2[j * numParticles + q])/Eij;
                    }else{
                        FORCAS10[j * numParticles + i] = 0.0;
                        FORCAS10[j * numParticles + i] = 0.0;
                        FORCAS10[j * numParticles + i] = 0.0;
                    }
                }
            }
        }
    }
}