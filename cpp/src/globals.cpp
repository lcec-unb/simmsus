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
#include <header/globals.hpp>

int nnr;
int npast;
int numRealizations;
int numParticles;
int totalRealParticle;
int k;
double l, h, qsi;
int numBoxes, numRecBoxes;
double shearratei;
double per, phi;
bool periodicity;
double *diam;
double *betaVec;
double *X0, *X1, *X2;
double *W0, *W1, *W2;
double *U0, *U1, *U2;
double *XI0, *XI1, *XI2;
double *ILF0, *ILF1, *ILF2;
double *ILR0, *ILR1, *ILR2;
double *DI0, *DI1, *DI2;
double *FORCAS10, *FORCAS11, *FORCAS12;
double *FORCAS20, *FORCAS21, *FORCAS22;
double *FORCAS30, *FORCAS31, *FORCAS32;
double *FORCAS40, *FORCAS41, *FORCAS42;
double *FORCAS50, *FORCAS51, *FORCAS52;
double *FORCAS60, *FORCAS61, *FORCAS62;
double *TORQUES10, *TORQUES11, *TORQUES12;
double *TORQUES20, *TORQUES21, *TORQUES22;
double *TORQUES30, *TORQUES31, *TORQUES32;
double *FT0, *FT1, *FT2;
double *TT0, *TT1, *TT2;
double *AUX0, *AUX1, *AUX2, *AUX3;
double *HAUX10, *HAUX11, *HAUX12;
double *HAUX20, *HAUX21, *HAUX22;
double *C0, *C1, *C2;
double *H10, *H11, *H12;
double *H20, *H21, *H22;
double *AUXT0, *AUXT1, *AUXT2;
double *T10, *T11, *T12;
double *T20, *T21, *T22;
double *F10, *F11, *F12;
double *F20, *F21, *F22;
double *AUXF0, *AUXF1, *AUXF2;
double *AUXT0, *AUXT1, *AUXT2;
double *cof01, *cof11, *cof02;
double *cof12, *cof03, *cof13;
double *cof04, *cof14, *cof05; 
double *cof15, *cof06, *cof16;
double *cof07, *cof17, *cof08;
double *cof18;