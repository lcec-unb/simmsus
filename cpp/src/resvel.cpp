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

#include <header/resvel.hpp>

double resvel(double a, double b, double c, double d)
{
    double k1 = b * (-a + d) / c;
    double k2 = b * ((-a - 0.5 * k1) + d) / c;
    double k3 = b * ((-a - 0.5 * k2) + d) / c;
    double k4 = b * ((-a - k3) + d) / c;

    return a += (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}