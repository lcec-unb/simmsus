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
// ! Subroutine resposible for calculating the parti-!
// ! cles velocity using a 4th order Runge-Kutta.    !
// !*************************************************!

// ! IT IS IMPORTANT TO NOTICE THAT THIS SUBROUTINE 
// ! DOES NOT APPLY FOR A MOBILITY PROBLEM. IT ONLY 
// ! MAKES SENSE WHEN WE CONSIDER PARTICLE INERTIA 
// ! (RESISTANCE FORMULATION)
// double a ! velocity component
// double b ! time step
// double c ! Stokes number (null for zero inertia)
// double d ! sum of forces in a given direction
// double k1,k2,k3,k4 ! internal variables

#include <header/resvel.hpp>

double resvel(double a, double b, double c, double d){
    double k1 = b * (-a + d) / c;
    double k2 = b * ((-a - 0.5 * k1)+d) / c;
    double k3 = b * ((-a - 0.5 * k2)+d) / c;
    double k4 = b * ((-a - k3) + d) / c;

    return a += (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}