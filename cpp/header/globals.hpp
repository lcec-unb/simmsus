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
#include <math.h>
#include <header/config.hpp>

#ifndef SRC_HEADERS_GLOBALS_HPP_
#define SRC_HEADERS_GLOBALS_HPP_

extern Configuration *configuration;
extern int nnr;
extern int npast;
extern int numRealizations;
extern int numParticles;
extern int num;
extern double qsi;
extern int k;
extern double l;
extern double h;
extern int totalRealParticle;
extern int numBoxes;
extern int numRecBoxes;
extern double shearratei; //= shearrate; //! shear-rate -> Arquivo de configuracao
extern double per; // = (4.0 / 3.0) * Pe; //! rotational Peclet number -> Arquivo de configuracao
extern double phi;
extern bool periodicity;
extern double *diam;
extern double *betaVec;
extern double *X0, *X1, *X2;
extern double *U0, *U1, *U2;
extern double *W0, *W1, *W2;
extern double *XI0, *XI1, *XI2;
extern double *ILF0, *ILF1, *ILF2;
extern double *ILR0, *ILR1, *ILR2;
extern double *DI0, *DI1, *DI2;
extern double *FORCAS10, *FORCAS11, *FORCAS12;
extern double *FORCAS20, *FORCAS21, *FORCAS22;
extern double *FORCAS30, *FORCAS31, *FORCAS32;
extern double *FORCAS40, *FORCAS41, *FORCAS42;
extern double *FORCAS50, *FORCAS51, *FORCAS52;
extern double *FORCAS60, *FORCAS61, *FORCAS62;
extern double *TORQUES10, *TORQUES11, *TORQUES12;
extern double *TORQUES20, *TORQUES21, *TORQUES22;
extern double *TORQUES30, *TORQUES31, *TORQUES32;
extern double *FT0, *FT1, *FT2;
extern double *TT0, *TT1, *TT2;
extern double *AUX0, *AUX1, *AUX2, *AUX3;
extern double *HAUX10, *HAUX11, *HAUX12;
extern double *HAUX20, *HAUX21, *HAUX22;
extern double *C0, *C1, *C2;
extern double *H10, *H11, *H12;
extern double *H20, *H21, *H22;
extern double *AUXT0, *AUXT1, *AUXT2;
extern double *T10, *T11, *T12;
extern double *T20, *T21, *T22;
extern double *F10, *F11, *F12;
extern double *F20, *F21, *F22;
extern double *AUXF0, *AUXF1, *AUXF2;
extern double *cof01, *cof11, *cof02;
extern double *cof12, *cof03, *cof13;
extern double *cof04, *cof14, *cof05; 
extern double *cof15, *cof06, *cof16;
extern double *cof07, *cof17, *cof08;
extern double *cof18;
#endif /* SRC_HEADERS_GLOBALS_HPP_ */