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
#include <math.h>
#include <globals.hpp>
#include <header/brownian.hpp>

void brownian(double *beta)
{
    double dt = configuration->getNumtimestep();
    double Pe = configuration->getBrownianpecletnum();
    bool gravidade = configuration->getSedimentation();
    bool shear = configuration->getTurnonshrate();
    bool torque = configuration->getSte();
    int total = 3 * numParticles * numRealizations;

    double *nr = new double[total];
    randomic(-1.0, 1.0, total, nr);
    double *diarand = new double[total];
    double nr1, nr2, nr3;
    int i, j;

    if (gravidade || shear)
    {
        for (j = 0; j < numRealizations; j++)
        {
            for (i = 0; i < numParticles; i++)
            {
                nr1 = nr[(i * 2 + (i - 2) + (numParticles * 3 * (j - 1)))];
                nr2 = nr[(i * 2 + (i - 1) + (numParticles * 3 * (j - 1)))];
                nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j - 1)))];

                double modrand = pow(pow(nr1, 2.0) + pow(nr2, 2.0) + pow(nr3, 2.0), 0.5);
                nr1 /= modrand;
                nr2 /= modrand;
                nr3 /= modrand;
                double resPow = (pow(beta[j * numParticles + i], (-2.0))) * (pow((6.0 / (Pe * dt)), 0.5));
                FORCAS60[j * numParticles + i] = resPow * nr1;
                FORCAS61[j * numParticles + i] = resPow * nr2;
                FORCAS62[j * numParticles + i] = resPow * nr3;
            }
        }
    }
    else
    {
        for (j = 0; j < numRealizations; j++)
        {
            for (i = 0; i < numParticles; i++)
            {
                nr1 = nr[(i * 2 + (i - 2) + (numParticles * 3 * (j - 1)))];
                nr2 = nr[(i * 2 + (i - 1) + (numParticles * 3 * (j - 1)))];
                nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j - 1)))];

                double modrand = pow(pow(nr1, 2.0) + pow(nr2, 2.0) + pow(nr3, 2.0), 0.5);
                nr1 /= modrand;
                nr2 /= modrand;
                nr3 /= modrand;
                double resPow = (pow(beta[j * numParticles + i], (-2.0))) * (pow((6.0 / dt), 0.5));
                FORCAS60[j * numParticles + i] = resPow * nr1;
                FORCAS61[j * numParticles + i] = resPow * nr2;
                FORCAS62[j * numParticles + i] = resPow * nr3;
            }
        }
    }

    if (torque)
    {
        randomic(-1.0, 1.0, (3 * numParticles * numRealizations), diarand);

        if (gravidade || shear)
        {
            for (j = 0; j < numRealizations; j++)
            {
                for (i = 0; i < numParticles; i++)
                {
                    nr1 = nr[(i * 2 + (i - 2) + (numParticles * 3 * (j - 1)))];
                    nr2 = nr[(i * 2 + (i - 1) + (numParticles * 3 * (j - 1)))];
                    nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j - 1)))];
                    double modrand = pow(pow(nr1, 2.0) + pow(nr2, 2.0) + pow(nr3, 2.0), 0.5);
                    nr1 /= modrand;
                    nr2 /= modrand;
                    nr3 /= modrand;
                    double resPow = pow(beta[j * numParticles + i], (-2.0)) * pow(9.0 / (2.0 * Pe * dt), 0.5);
                    TORQUES30[j * numParticles + i] = resPow * nr1;
                    TORQUES31[j * numParticles + i] = resPow * nr2;
                    TORQUES32[j * numParticles + i] = resPow * nr3;
                }
            }
        }
        else
        {
            for (j = 0; j < numRealizations; j++)
            {
                for (i = 0; i < numParticles; i++)
                {
                    nr1 = nr[(i * 2 + (i - 2) + (numParticles * 3 * (j - 1)))];
                    nr2 = nr[(i * 2 + (i - 1) + (numParticles * 3 * (j - 1)))];
                    nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j - 1)))];
                    double modrand = pow(pow(nr1, 2.0) + pow(nr2, 2.0) + pow(nr3, 2.0), 0.5);
                    nr1 /= modrand;
                    nr2 /= modrand;
                    nr3 /= modrand;
                    double resPow = pow(beta[j * numParticles + i], (-2.0)) * pow(9.0 / (2.0 * dt), 0.5);
                    TORQUES30[j * numParticles + i] = resPow * nr1;
                    TORQUES31[j * numParticles + i] = resPow * nr2;
                    TORQUES32[j * numParticles + i] = resPow * nr3;
                }
            }
        }
    }