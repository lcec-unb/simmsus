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

/*************************************************
 Subroutine responsible for creating the initial 
 particle distribution for all simultaneous 	  
 numerical experiments.			  
						  
 Several possible initial configurations are im- 
 plemented here. These include:		  
						  
 1 - An initial random distribution		  
 2 - An ordered NxNxN array of particles 	  
 3 - An initial spherical aggregate		  
*************************************************/
#include <iostream>
#include <header/globals.hpp>
#include <header/randomic.hpp>
#include <header/initialCondition.hpp>
#include <math.h>
using namespace std;

void initialCondition(bool initialSphericalAggregate){

//Calculating the aggregate's radius   
double ragreg = pow((numParticles / phi),(1.0/3.0));

//Inserting the aggregate in a giant lattice to avoid periodicity condition 
l = 100.0 * ragreg;
h = l;
double r, nr1, nr2, nr3, modrand;
int total = 3 * numParticles * numRealizations;
double *nr = new double[total]{};
double *nra = new double[total]{};
double *nrs = new double[3]{};
int *nr0Index = new int[total]{};
int *nr1Index = new int[total]{};
int *nr2Index = new int[total]{};
randomic(-1.0,1.0,(total),nr);
radomicAccess(0,total - 1,(total),nr0Index);
radomicAccess(0,total - 1,(total),nr1Index);
radomicAccess(0,total - 1,(total),nr2Index);
randomic(-1.0,1.0,(total),nr);
randomic(0.0,1.0,(total),nra);

// If we are starting a new simulation, { we will check if the initial condition is a spherical aggregate
if(initialSphericalAggregate) {

    
    // Defining the center of the aggregate

    double xcentro = l / 2.0;
    double ycentro = l / 2.0;
    double zcentro = h - 2.0 * ragreg;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = l/2.0 - ragreg;
    xmax = xmin + 2.0 * ragreg;

    for(int j = 0; j < numRealizations; j++){
        for(int i = 0; i < numParticles; i++){
            // Range of X position of the particles in the aggregate
            // verificar se precisa somar 1 ao i e j -> Fortran C++
            X0[j * numParticles + i] = xcentro + (xmax - xmin) * 0.5 * nr[nr0Index[j * numParticles + i]];

            // Using the circle equation to distribute the particles in the y direction (sphere projection on a xy plane)
            double yc = pow(pow(ragreg,2.0) - pow(X0[j * numParticles + i] - xcentro,2.0),0.5);
            ymin = ycentro - yc;
            ymax = ycentro + yc;

            X1[j * numParticles + i] = ycentro +(ymax - ymin) * 0.5 * nr[nr1Index[j * numParticles + i]];

            // Using the circle equation to define the z positions considering now the sphere projection on planes zy and zx
            double zc = pow(pow(ragreg,2.0) - pow(X0[j * numParticles + i] - xcentro,2.0) - pow(X1[j * numParticles + i] - ycentro,2.0),0.5);    
            zmin = zcentro - zc;
            zmax = zcentro + zc;

            X2[j * numParticles + i] = zcentro +(zmax - zmin) * 0.5 * nr[nr2Index[j * numParticles + i]];
        }
    }
}else{

    // If you want to generate a random initial distribution, { you must call the random number generation routine...

    for(int j = 0; j < numRealizations; j++){
        for(int i = 0; i < numParticles; i++){
            X0[j * numParticles + i] = l * nra[nr0Index[j * numParticles + i]]; 
            X1[j * numParticles + i] = l * nra[nr1Index[j * numParticles + i]]; 
            X2[j * numParticles + i] = h * nra[nr2Index[j * numParticles + i]]; 
        }
    }

    // Now that we already created an initial random configuration, we must check for particle overlaps

    // We call the random number generation subroutine
    randomic(-1.0,1.0,3, nrs);
    for(int k = 0; k < numRealizations; k++){
            for(int i = 0; i < numParticles; i++){
                for(int j = 0; j < numParticles; j++){
                    if(i != j){

                        // Calculate the distance between all the pair of particles in the suspension
                        r = pow(pow(X0[k* numParticles + i] - X0[k* numParticles + j],2.0) + pow(X1[k* numParticles + i]-X1[k* numParticles + j],2.0) + pow(X2[k* numParticles + i] - X2[k* numParticles + j],2.0),0.5);

                        // If at any point the distance is smaller { 0.01*a, { we give a Brownian kick in the overlaped particles
                        if(r <= 2.01){

                        // Creating a small random vector used to "kick" one of the overlaped particles

                        nr1 = nrs[0];
                        nr2 = nrs[1];
                        nr3 = nrs[2];

                        modrand = pow(pow(nr1,2.0) + pow(nr2,2.0) + pow(nr3,2.0),0.5);

                        // Normalizing this vector

                        nr1 /= modrand;
                        nr2 /= modrand;
                        nr3 /= modrand;

                        nr1 *= 0.25;
                        nr2 *= 0.25;
                        nr3 *= 0.25;

                        // Changing the position of one of the overlaped particles with a Brownian "kick"

                        X0[k* numRealizations + i] += nr1;
                        X1[k* numRealizations + i] += nr2;
                        X2[k* numRealizations + i] += nr3;

                        }
                    // Calculate the new distance between the "problematic" particles
                    r = pow(pow(X0[k* numParticles + i] - X0[k* numParticles + j],2.0) + pow(X1[k* numParticles + i]-X1[k* numParticles + j],2.0) + pow(X2[k* numParticles + i] - X2[k* numParticles + j],2.0),0.5);


                    // If at any point the particles are still overlaped we give new 
                    // Brownian kicks in the particles that present this problem until the situation is solved

                //     if(r <= 2.01){
                //     //para checar o que fazer aqui dentro que nao seja 

                //     // go to 103
                //     }

                // // Imposing conditions to avoid particles outside the box

                // if(X0[k* numRealizations + i] < 0) {
                //     X0[k* numRealizations + i] = abs(X0[k* numRealizations + i]);
                //     // go to 103
                // }
                // if(X1[k* numRealizations + i] < 0) {
                //     X1[k* numRealizations + i] = abs(X1[k* numRealizations + i]);
                //     // go to 103
                // }
                // if(X2[k* numRealizations + i] < 0) {
                //     X2[k* numRealizations + i] = abs(X2[k* numRealizations + i]);
                //     // go to 103
                // }

                // if(X0[k* numRealizations + i] > l) {
                //     X0[k* numRealizations + i] -= l;
                //     // go to 103
                // }
                // if(X1[k* numRealizations + i] > l) {
                //     X1[k* numRealizations + i] -= l;
                //     // go to 103
                // }
                // if(X2[k* numRealizations + i] > l) {
                //     X2[k* numRealizations + i] -= l;
                //     // go to 103
                // }
                }
                }
                // Lets define the initial velocity of the particles as the Stokes velocity (mobility problem)
               
                for(int j = 0; j < numRealizations; j++){
                    for(int i = 0; i < numParticles; i++){
                        U0[j * numParticles + i] = 0.0;
                        U1[j * numParticles + i] = 0.0;
                        U2[j * numParticles + i] = -1.0;
                    }
                } 
}
    }
}
}
