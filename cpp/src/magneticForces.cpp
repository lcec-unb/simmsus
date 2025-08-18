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
#include <numeric>
#include <math.h>
#include <header/randomic.hpp>
#include <header/globals.hpp>
#include <header/magneticForces.hpp>

using namespace std;

void magneticForces()
{
    double alpha2 = configuration->getAlpha2();
    double percentual = configuration->getPercentnonmagpart();
    bool mistura = configuration->getMixmagnonmagpart();
    bool gravidade = configuration->getSedimentation();
    bool shear = configuration->getTurnonshrate();
    double brownianmpecletnum = configuration->getBrownianpecletnum();
    int i, j, q;
    int auxiliary;
    double r, rij[3], modrij, lambda;
    double termo0, termo1, termo2, termo3, sumAux;
    auxiliary = 1;
    double *aux0, *aux1, *aux2;
    aux0 = new double[numRealizations * numParticles]{};
    aux1 = new double[numRealizations * numParticles]{};
    aux2 = new double[numRealizations * numParticles]{};

    if (mistura)
    {
        auxiliary = (percentual * numParticles) + 1;
    }

    if (gravidade)
    {
        lambda = alpha2 * 6.0 / brownianmpecletnum;
    }
    else
    {
        if (shear)
        {
            lambda = alpha2 * 6.0 / brownianmpecletnum;
        }
        else
        {
            lambda = alpha2 * 6.0;
        }
    }

    for (j = 0; j < numRealizations; j++)
    {
        for (i = auxiliary; i < numParticles; i++)
        {
            for (q = auxiliary; q < numParticles; q++)
            {
                // vetorizar
                if (i != q)
                {
                    //! Checking the distance between the particles
                    r = pow(pow((X0[j * numParticles + i] - X0[j * numParticles + q]), 2.0) + pow((X1[j * numParticles + i] - X1[j * numParticles + q]), 2.0) + pow((X2[j * numParticles + i] - X2[j * numParticles + q]), 2.0), 0.5);

                    // If particles are to close we turn off magnetic interactions to avoid overlap
                    if (r <= 2.2)
                    {
                        aux0[j * numParticles + q] = 0.0;
                        aux1[j * numParticles + q] = 0.0;
                        aux2[j * numParticles + q] = 0.0;
                    }
                    else
                    {
                        //! Calculating the vector R_{ij} which connects a particle i to a particle j
                        rij[0] = X0[j * numParticles + i] - X0[j * numParticles + q];
                        rij[1] = X1[j * numParticles + i] - X1[j * numParticles + q];
                        rij[2] = X2[j * numParticles + i] - X2[j * numParticles + q];
                        //! Normalizing vector R_{ij}
                        modrij = pow(pow(rij[0], 2.0) + pow(rij[1], 2.0) + pow(rij[2], 2.0), 0.5);
                        rij[0] = rij[0] / modrij;
                        rij[1] = rij[1] / modrij;
                        rij[2] = rij[2] / modrij;

                        termo0 = (DI0[j * numParticles + i] * DI0[j * numParticles + q]);
                +(DI1[j * numParticles + i]*DI1[j * numParticles + q])+(DI2[j * numParticles + i]*DI2[j * numParticles + q]))*rij[0];
                termo1 = ((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) + (DI2[j * numParticles + i] * rij[2])) * DI0[j * numParticles + q];
                termo2 = ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2])) * DI0[j * numParticles + i];
                termo3 = (((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) + (DI2[j * numParticles + i] * rij[2])) * ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2]))) * rij[0];

                aux0[j * numParticles + q] = (lambda / (1.0 * pow(r, 4.0))) * (termo0 + termo1 + termo2 - 5.0 * termo3);

                termo0 = ((DI0[j * numParticles + i] * DI0[j * numParticles + q]) + (DI1[j * numParticles + i] * DI1[j * numParticles + q]) + (DI2[j * numParticles + i] * DI2[j * numParticles + q])) * rij[1];
                termo1 = ((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) + (DI2[j * numParticles + i] * rij[2])) * DI1[j * numParticles + q];
                termo2 = ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2])) * DI1[j * numParticles + i];
                termo3 = (((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) + (DI2[j * numParticles + i] * rij[2])) * ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2]))) * rij[1];

                aux1[j * numParticles + q] = (lambda / (1.0 * pow(r, 4.0))) * (termo0 + termo1 + termo2 - 5.0 * termo3);

                termo0 = ((DI0[j * numParticles + i] * DI0[j * numParticles + q]) + (DI1[j * numParticles + i] * DI1[j * numParticles + q]) + (DI2[j * numParticles + i] * DI2[j * numParticles + q])) * rij[2];
                termo1 = ((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) + (DI2[j * numParticles + i] * rij[2])) * DI2[j * numParticles + q];
                termo2 = ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2])) * DI2[j * numParticles + i];
                termo3 = (((DI0[j * numParticles + i] * rij[0]) + (DI1[j * numParticles + i] * rij[1]) +
                           (DI2[j * numParticles + i] * rij[2])) *
                          ((DI0[j * numParticles + q] * rij[0]) + (DI1[j * numParticles + q] * rij[1]) + (DI2[j * numParticles + q] * rij[2]))) *
                         rij[2];

                aux2[j * numParticles + q] = (lambda / (1.0 * pow(r, 4.0))) * (termo0 + termo1 + termo2 - 5.0 * termo3);
                    }
                }
            }
            sumAux = reduce(aux0[0], aux0[numRealizations * numParticles], aux0[0]);
            FORCAS40[j * numParticles + i] = sumAux;
            sumAux = reduce(aux1[0], aux1[numRealizations * numParticles], aux1[0]);
            FORCAS41[j * numParticles + i] = sumAux;
            sumAux = reduce(aux2[0], aux2[numRealizations * numParticles], aux2[0]);
            FORCAS42[j * numParticles + i] = sumAux;
            aux0[j * numParticles + i] = 0.0;
            aux1[j * numParticles + i] = 0.0;
            aux2[j * numParticles + i] = 0.0;
        }
    }
}