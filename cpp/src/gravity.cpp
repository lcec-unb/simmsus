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
#include <randomic.hpp>
#include <cmath>
#include <globals.hpp>

using namespace std;

void gravity(double *beta)
{
    bool gravidade = configuration->getSedimentation();
    int i, j;

    if (gravidade)
    {

        for (j = 0; j < numRealizations; j++)
        {
            for (i = 0; i < numParticles; i++)
            {

                FORCAS30[j * numParticles + i] = 0.0;
                FORCAS31[j * numParticles + i] = 0.0;
                FORCAS32[j * numParticles + i] = pow(-beta[j * numParticles + i], 3.0);
            }
        }
    }
    else
    {
        for (j = 0; j < numRealizations; j++)
        {
            for (i = 0; i < numParticles; i++)
            {
                FORCAS30[j * numParticles + i] = 0.0;
                FORCAS31[j * numParticles + i] = 0.0;
                FORCAS32[j * numParticles + i] = 0.0;
            }
        }
    }
}