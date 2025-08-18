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
#include <header/magneticTorque.hpp>


void magneticTorque()
{

    bool mistura = configuration->getMixmagnonmagpart();
    bool gravidade = configuration->getSedimentation();
    bool shear = configuration->getTurnonshrate();
    double brownianmpecletnum = configuration->getBrownianpecletnum();

    int auxiliary = 1;
    double r, rij[3], modrij, lambda, alpha2;
    double termo1, termo2, termo3, percentual;
    double aux1[numRealizations][numParticles];
    double aux2[numRealizations][numParticles];
    double aux3[numRealizations][numParticles];

    if (gravidade)
    {
        lambda = alpha2 * 24.0;
    }
    else
    {
        if (shear)
        {
            lambda = alpha2 * 3.0 / brownianmpecletnum;
        }
        else
        {
            lambda = alpha2 * 9.0 / 2.0;
        }
    }

    // ! These subroutine may consider a magnetic fluidized
    // ! bed mixing magnetic and non-magnetic particles,
    // ! the logical variable "mistura" activates this po-
    // ! ssibility. If mistura = TRUE { we have to set
    // ! zero dipoles for a certain percentage "percentual"
    // ! of the particles.

    if (mistura)
    {
        auxiliary = (percentual * numParticles) + 1;
    }

    for (int j = 0; j < numRealizations; j++)
    {
        for (int i = auxiliary; i < numParticles; i++)
        {
            for (int q = auxiliary; q < numParticles; q++)
            {
                if (i != q)
                {

                    // ! Calculating the distance between all the pairs of particles
                    r = pow(pow(X0[j * numParticles + i] - X0[j * numParticles + q], 2.0) +
                                pow((X1[j * numParticles + i] - X1[j * numParticles + q]), 2.0) +
                                pow((X2[j * numParticles + i] - X2[j * numParticles + q]), 2.0),
                            0.5);
                    if (r <= 2.0)
                    {
                        aux1[j][q] = 0.0;
                        aux2[j][q] = 0.0;
                        aux3[j][q] = 0.0;
                    }
                }
                else
                {
                    // ! Calculating the vector that connects a particle i to a particle j
                    rij[0] = X0[j * numParticles + i] - X0[j * numParticles + q];
                    rij[1] = X1[j * numParticles + i] - X1[j * numParticles + q];
                    rij[2] = X2[j * numParticles + i] - X2[j * numParticles + q];

                    // ! Normalizing this vector
                    modrij = pow(pow(rij[0], 2.0) + pow(rij[1], 2.0) + pow(rij[2], 2.0), 0.5);
                    rij[0] /= modrij;
                    rij[1] /= modrij;
                    rij[2] /= modrij;

                    termo2 = ((DI0[j * numParticles + i] * rij[0]) +
                              (DI1[j * numParticles + q] * rij[1]) +
                              (DI2[j * numParticles + q] * rij[2]));

                    termo1 = ((DI1[j * numParticles + i] * DI2[j * numParticles + q]) -
                              (DI2[j * numParticles + i] * DI1[j * numParticles + q]));

                    termo3 = ((DI1[j * numParticles + i] * rij[2]) -
                              (DI2[j * numParticles + i] * rij[1]));

                    aux1[j][q] = (lambda / (1.0 * pow(r, 3.0))) * (((-1.0 / 3.0) * termo1) + (termo2 * termo3));

                    termo1 = (DI2[j * numParticles + i] * DI0[j * numParticles + q]) - (DI0[j * numParticles + i] * DI2[j * numParticles + q]);
                    termo3 = ((DI2[j * numParticles + i] * rij[0]) - (DI0[j * numParticles + i] * rij[2]));

                    aux2[j][q] = (lambda / (1.0 * pow(r, 3.0))) * (((-1.0 / 3.0) * termo1) + (termo2 * termo3));

                    termo1 = ((DI0[j * numParticles + i] * DI1[j * numParticles + q]) - (DI1[j * numParticles + i] * DI0[j * numParticles + q]));
                    termo3 = ((DI0[j * numParticles + i] * rij[1]) - (DI1[j * numParticles + i] * rij[0]));

                    aux3[j][q] = (lambda / (1.0 * pow(r, 3.0))) * (((-1.0 / 3.0) * termo1) + (termo2 * termo3));
                }
            }
        }
        TORQUES10[j * numParticles + i] = sum(aux1(j, :));
        TORQUES11[j * numParticles + i] = sum(aux2(j, :));
        TORQUES12[j * numParticles + i] = sum(aux3(j, :));

        aux1 = 0.0;
        aux2 = 0.0;
        aux3 = 0.0;
        r = 0.0;
        modrij = 0.0;
    }
}
}