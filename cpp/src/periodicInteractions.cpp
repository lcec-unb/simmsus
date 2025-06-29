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

// *************************************************!
//  		     SIMMSUS			  ! 
// SUBROUTINE: periodic_interactions		  !					         
// Last update: 16/07/2023			  !
// *************************************************!

// *************************************************!
//  Suroutine responsible for computing properly all!
//  periodic interactions. 			  !
// 						  !
//  This subroutine is called whenever the program  !
//  identifies the presence of any periodic interac-!
//  tion. These periodic interactions can be of the !
//  following kinds:				  !
// 						  !
//  1 - Long range hydrodynamic interactions	  !
//  2 - Long range periodic magnetic torques        !
//  3 - Long range periodic magnetic forces  	  !
// 						  !
//  This is probably the most expensive and complex !
//  subroutine in this code.			  !
// *************************************************!
#include <iostream>
#include <numeric>
#include <globals.hpp>
#include <math.h>
 #include <cstring>
using namespace std;
//Transformar o objeto de configuração em global

void periodicInteractions(double dt, double shearrate, double alpha2, double lambda, bool shear, bool ligaih, bool tmagper, bool fmagper, bool gravidade, double brownianmpecletnum){
// We start by computing first the sums in the real space 
// Here, nb denotes the number of physical boxes. This
// number can be different from nbr, which denotes the
// nubmer of reciprocal boxes. In simconfig.dat we 
// recommend the user to set nb=125 and nbr=27. These
// values ensure a good precision when computing the 
// average sedimentation velocity of a suspension of
// spheres in CreeM_PIng-flow

// This loop works like this: we fix a realization and
// a given particle and the make a sweep computing the
// long-range interactions between this particles and
// all the other particles in the surrounding boxes.
// We {, extend this procedure to all real particles
// in all the simultaneous numerical experiments.

int q,i,s,j,d;
double rij[3], konda[3], knormal[3], rn[3], auxf[3][3], auxt[3][3];
double modrij, kr, kr2, modk;
double diferenca_interpol1;
double mobilidade_self[3][3];
double mobs0, mobs1, mobs2;
double coeficiente1,coeficiente2, coeficiente3;
double coeficiente4, coeficiente5, coeficiente6;
double coeficiente7, coeficiente8;
double termo1, termo2, termo3, termo4, termo5;
double hidrodinamica_aux1[3][3];
double mobilidade1[3][3], mobilidade2[3][3];
double forcareciproca[3][3], torquereciproco[3][3];
double contribuicao_reciproco[3][3], contribuicao_self[3][3], contribuicao_fisico[3][3]; //numpart * numrealizations
double aux_periodico;
double *forcareal0 = new double[numBoxes];
double *forcareal1 = new double[numBoxes];
double *forcareal2 = new double[numBoxes];
double *torquereal0 = new double[numBoxes];
double *torquereal1 = new double[numBoxes];
double *torquereal2 = new double[numBoxes];
double *FREC0 = new double[numRecBoxes];
double *FREC1 = new double[numRecBoxes];
double *FREC2 = new double[numRecBoxes];
double *TREC0 = new double[numRecBoxes];
double *TREC1 = new double[numRecBoxes];
double *TREC2 = new double[numRecBoxes];



for (q == 0; q < numRealizations; q++){
    for (i == 0; i < numParticles; i++){
        for (s == 0; s < numBoxes; s++){
            for (j == 0; j < numParticles; j++){    

                int aux_periodico = abs(ILF0[s]) + abs(ILF1[s]) + abs(ILF2[s]);

                // Checking the distance between a particle "i"  and other particle "j" in the real boxes (physical and images)

                // Calculating the "r" vector (distance)

                if (!shear)
                {
                    rij[0] = X0[q * numParticles + i] - (X0[q * numParticles + j] + (ILF0[s] * l));
                    rij[1] = X1[q * numParticles + i] - (X1[q * numParticles + j] + (ILF1[s] * l));
                    rij[2] = X2[q * numParticles + i] - (X2[q * numParticles + j] + (ILF2[s] * h));
                }
                else
                {
                    rij[0] = X0[q * numParticles + i] - (X0[q * numParticles + j] + (ILF0[s] * l));
                    rij[1] = X1[q * numParticles + i] - (X1[q * numParticles + j] + (ILF1[s] * (k * ILF2[s] * shearrate * dt))); // arrumar shearrate
                    rij[2] = X2[q * numParticles + i] - (X2[q * numParticles + j] + (ILF2[s] * h));
                }

                // Normalizing the "r" vector -> r/|r|

                modrij = (pow(rij[0], 2.0) + pow(rij[1], 2.0) + pow(pow(rij[2], 2.0), 0.5));

                rn[0] = rij[0] / modrij;
                rn[1] = rij[1] / modrij;
                rn[2] = rij[2] / modrij;

                d = 1.0 + (9999.0 * (modrij - 2)) / (pow(3.0, 0.5) * l * pow(numBoxes, (1.0 / 3.0)));
                diferenca_interpol1 = (modrij - cof01[d]);

                // ************ SUMS IN THE PHYSICAL SPACE IN THE REAL BOX (WHERE THE REAL PARTICLES ARE) AND "IMAGE" BOXES***********!
                // ******************************* HYDRODYNAMIC INTERACTIONS *********************************************************!

                // Building the self-mobility matrix

                if (aux_periodico == 0.0)
                {
                    if (ligaih)
                    {
                        if (i == j)
                        {
                            mobs0 = 1.0 - (6.0 * (pow(M_PI, -0.5)) * qsi) + ((40.0 / 3.0) * (pow(M_PI, -0.5)) * (pow(qsi, 3.0)));
                            mobs1 = 1.0 - (6.0 * (pow(M_PI, -0.5)) * qsi) + ((40.0 / 3.0) * (pow(M_PI, -0.5)) * (pow(qsi, 3.0)));
                            mobs2 = 1.0 - (6.0 * (pow(M_PI, -0.5)) * qsi) + ((40.0 / 3.0) * (pow(M_PI, -0.5)) * (pow(qsi, 3.0)));
                        }
                    }

                    // Building the called mobility matrix 1 (lattice sum in the physical space)

                    if (ligaih)
                    {
                        if (modrij > (2.0))
                        {

                            // Interpolating the pre-calculated Green functions

                            if (diferenca_interpol1 > 0.0)
                            {
                                coeficiente1 = cof11[d] + ((modrij - cof01[d]) / (cof01[d + 1] - cof01[d])) * (cof11[d + 1] - cof11[d]);
                                coeficiente2 = cof12[d] + ((modrij - cof02[d]) / (cof02[d + 1] - cof02[d])) * (cof12[d + 1] - cof12[d]);
                            }
                            else
                            {
                                coeficiente1 = cof11[d - 1] + ((modrij - cof01[d - 1]) / (cof01[d] - cof01[d - 1])) * (cof11[d] - cof11[d - 1]);
                                coeficiente2 = cof12[d - 1] + ((modrij - cof02[d - 1]) / (cof02[d] - cof02[d - 1])) * (cof12[d] - cof12[d - 1]);
                            }

                            mobilidade1[0][0] = coeficiente1 + coeficiente2 * (rn[0] * rn[0]);
                            mobilidade1[0][1] = coeficiente2 * (rn[0] * rn[1]);
                            mobilidade1[0][2] = coeficiente2 * (rn[0] * rn[2]);

                            mobilidade1[1][0] = coeficiente2 * (rn[1] * rn[0]);
                            mobilidade1[1][1] = coeficiente1 + coeficiente2 * (rn[1] * rn[1]);
                            mobilidade1[1][2] = coeficiente2 * (rn[1] * rn[2]);

                            mobilidade1[2][0] = coeficiente2 * (rn[2] * rn[0]);
                            mobilidade1[2][1] = coeficiente2 * (rn[2] * rn[1]);
                            mobilidade1[2][2] = coeficiente1 + coeficiente2 * (rn[2] * rn[2]);

                            HAUX10[j] = mobilidade1[0][0] * FT0[q * numParticles + j] + mobilidade1[0][1] * FT1[q * numParticles + j] + mobilidade1[0][2] * FT2[q * numParticles + j];
                            HAUX11[j] = mobilidade1[1][0] * FT0[q * numParticles + j] + mobilidade1[1][1] * FT1[q * numParticles + j] + mobilidade1[1][2] * FT2[q * numParticles + j];
                            HAUX12[j] = mobilidade1[2][0] * FT0[q * numParticles + j] + mobilidade1[2][1] * FT1[q * numParticles + j] + mobilidade1[2][2] * FT2[q * numParticles + j];
                        } // ending the "if" of the condition (modrij>2.0)
                    } // ending the "if" of the condition (ligaih)

                    // ****************** MAGNETIC TORQUES **************************************!

                    if (tmagper)
                    {
                        if (modrij > 2.0)
                        {

                            //! Interpolating the pre-calculated Green functions

                            if (diferenca_interpol1 > 0.0)
                            {
                                coeficiente4 = cof14[d] + ((modrij - cof04[d]) / (cof04[d + 1] - cof04[d])) * (cof14[d + 1] - cof14[d]);
                                coeficiente5 = cof15[d] + ((modrij - cof05[d]) / (cof05[d + 1] - cof05[d])) * (cof15[d + 1] - cof15[d]);
                            }
                            else
                            {
                                coeficiente4 = cof14[d - 1] + ((modrij - cof04[d - 1]) / (cof04[d] - cof04[d - 1])) * (cof14[d] - cof14[d - 1]);
                                coeficiente5 = cof15[d - 1] + ((modrij - cof05[d - 1]) / (cof05[d] - cof05[d - 1])) * (cof15[d] - cof15[d - 1]);
                            }

                            // Computing periodic torques due to magnetic interactions in the real space
                            // !eps=10.0

                            // !if(aux_periodico==0.0){
                            // !termo5= 4.0*M_PI/((2.0+eps)*(l**3.0)) ! Surface term (only present in the physical lattice)
                            // !else
                            // !termo5=0.0
                            // !}

                            termo5 = 0.0;

                            termo2 = ((DI0[q * numParticles + j] * rij[0]) + (DI1[q * numParticles + j] * rij[1]) + (DI2[q * numParticles + j] * rij[2]));

                            if (gravidade)
                            {
                                lambda = alpha2 * 6.0 * M_PI / brownianmpecletnum;
                            }
                            else
                            {
                                lambda = alpha2 * 6.0 * M_PI;
                            }

                            termo1 = (DI2[q * numParticles + i] * DI1[q * numParticles + j]) - (DI1[q * numParticles + i] * DI2[q * numParticles + j]);
                            termo3 = (DI1[q * numParticles + i] * rij[2]) - (DI2[q * numParticles + i] * rij[1]);
                            AUXT0[j] = lambda * ((termo1 * coeficiente4) + (termo2 * termo3 * coeficiente5) + termo5 * termo1);

                            termo1 = (DI0[q * numParticles + i] * DI2[q * numParticles + j]) - (DI2[q * numParticles + i] * DI0[q * numParticles + j]);
                            termo3 = (DI2[q * numParticles + i] * rij[0]) - (DI0[q * numParticles + i] * rij[2]);
                            AUXT1[j] = lambda * ((termo1 * coeficiente4) + (termo2 * termo3 * coeficiente5) + termo5 * termo1);

                            termo1 = (DI1[q * numParticles + i] * DI0[q * numParticles + j]) - (DI0[q * numParticles + i] * DI1[q * numParticles + j]);
                            termo3 = (DI0[q * numParticles + i] * rij[1]) - (DI1[q * numParticles + i] * rij[0]);
                            AUXT2[j] = lambda * ((termo1 * coeficiente4) + (termo2 * termo3 * coeficiente5) + termo5 * termo1);

                        } // ending the "if" of the condition (modrij>2.0)
                    } // ending the "if" of the condition (tmagper)

                    // *********************************** MAGNETIC FORCES *****************************************!
                    if (fmagper)
                    {
                        if (modrij > (2.0))
                        {

                            // Interpolating the pre-calculated Green functions

                            if (diferenca_interpol1 > 0.0)
                            {
                                coeficiente6 = cof16[d] + ((modrij - cof06[d]) / (cof06[d + 1] - cof06[d])) * (cof16[d + 1] - cof16[d]);
                            }
                            else
                            {
                                coeficiente6 = cof16[d - 1] + ((modrij - cof06[d - 1]) / (cof06[d] - cof06[d - 1])) * (cof16[d] - cof16[d - 1]);
                            }

                            //! Computing periodic forces due to magnetic interactions in the real space

                            if (gravidade)
                            {
                                lambda = alpha2 * 8.0 * M_PI / brownianmpecletnum;
                            }
                            else
                            {
                                lambda = alpha2 * 8.0 * M_PI;
                            }

                            termo1 = ((DI0[q * numParticles + i] * DI0[q * numParticles + j]) +
                                      (DI1[q * numParticles + i] * DI1[q * numParticles + j]) +
                                      (DI2[q * numParticles + i] * DI2[q * numParticles + j])) *
                                     rij[0];

                            termo3 = ((DI0[q * numParticles + j] * rij[0]) +
                                      (DI1[q * numParticles + j] * rij[1]) +
                                      (DI2[q * numParticles + j] * rij[2])) *
                                     DI0[q * numParticles + i];

                            termo2 = ((DI0[q * numParticles + i] * rij[0]) +
                                      (DI1[q * numParticles + i] * rij[1]) +
                                      (DI2[q * numParticles + i] * rij[2])) *
                                     DI0[q * numParticles + j];

                            termo4 = (((DI0[numParticles + i]) * rij[0]) + (DI1[numParticles + i] * rij[1]) +
                                      (DI2[numParticles + i] * rij[2])) *
                                     ((DI0[q * numParticles + j] * rij[0]) +
                                      (DI1[q * numParticles + j] * rij[1]) + (DI2[q * numParticles + j] * rij[2])) *
                                     rij[0];

                            AUXF0[j] = lambda * (((termo1 + termo2 + termo3) * coeficiente5) - termo4 * coeficiente6);

                            termo1 = ((DI0[q * numParticles + i] * DI0[q * numParticles + j]) +
                                      (DI1[q * numParticles + i] * DI1[q * numParticles + j]) +
                                      (DI2[q * numParticles + i] * DI2[q * numParticles + j])) *
                                     rij[1];

                            termo3 = ((DI0[q * numParticles + j] * rij[0]) +
                                      (DI1[q * numParticles + j] * rij[1]) +
                                      (DI2[q * numParticles + i] * rij[2])) *
                                     DI1[q * numParticles + i];

                            termo2 = ((DI0[q * numParticles + i] * rij[0]) +
                                      (DI1[q * numParticles + i] * rij[1]) +
                                      (DI2[q * numParticles + i] * rij[2])) *
                                     DI1[q * numParticles + j];

                            termo4 = (((DI0[numParticles + i] * rij[0]) + (DI1[numParticles + i] * rij[1]) +
                                       (DI2[numParticles + i] * rij[2])) *
                                      ((DI0[q * numParticles + j] * rij[0]) +
                                       (DI1[q * numParticles + i] * rij[1]) +
                                       (DI2[q * numParticles + j] * rij[2]))) *
                                     rij[1];

                            AUXF1[j] = lambda * (((termo1 + termo2 + termo3) * coeficiente5) - termo4 * coeficiente6);

                            termo1 = ((DI0[q * numParticles + i] * DI0[q * numParticles + j]) +
                                      (DI1[q * numParticles + i] * DI1[q * numParticles + j]) +
                                      (DI2[q * numParticles + i] * DI2[q * numParticles + i])) *
                                     rij[2];

                            termo3 = ((DI0[q * numParticles + j] * rij[0]) + (DI1[q * numParticles + j] * rij[1]) + (DI2[q * numParticles + i] * rij[2])) * DI2[q * numParticles + i];

                            termo2 = ((DI0[q * numParticles + i] * rij[0]) + (DI1[q * numParticles + i] * rij[1]) + (DI2[q * numParticles + i] * rij[2])) * DI2[q * numParticles + j];

                            termo4 = (((DI0[numParticles + i] * rij[0]) + (DI1[numParticles + i] * rij[1]) +
                                       (DI2[numParticles + i] * rij[2])) *
                                      ((DI0[q * numParticles + j] * rij[0]) +
                                       (DI1[q * numParticles + j] * rij[1]) + (DI2[q * numParticles + j] * rij[2]))) *
                                     rij[2];

                            AUXF2[j] = lambda * (((termo1 + termo2 + termo3) * coeficiente5) - termo4 * coeficiente6);
                        } //! ending the "if" of the condition  (modrij>2.0)
                    } //! ending the "if" of the condition (fmagper)
                }
                // ************************** SUMING EVERYTHING *********************************!

                // Computing the sum of the previous vectors and matrices

                if (ligaih)
                {
                    H10[s] = reduce(HAUX10[0], HAUX10[numParticles], HAUX10[0]);
                    H11[s] = reduce(HAUX11[0], HAUX11[numParticles], HAUX11[0]);
                    H12[s] = reduce(HAUX12[0], HAUX12[numParticles], HAUX12[0]);
                }

                if (fmagper)
                {
                    forcareal0[s] = reduce(AUX0[0], AUX0[numRecBoxes], AUX0[0]);
                    forcareal1[s] = reduce(AUX1[0], AUX1[numRecBoxes], AUX1[0]);
                    forcareal2[s] = reduce(AUX2[0], AUX2[numRecBoxes], AUX2[0]);
                }

                if (tmagper)
                {
                    torquereal0[s] = reduce(AUXT0[0], AUXT0[numRecBoxes], AUXT0[0]);
                    torquereal1[s] = reduce(AUXT1[0], AUXT1[numRecBoxes], AUXT1[0]);
                    torquereal2[s] = reduce(AUXT2[0], AUXT2[numRecBoxes], AUXT2[0]);
                    memset(AUXT0, 0, sizeof(AUXT0));
                    memset(AUXT1, 0, sizeof(AUXT1));
                    memset(AUXT2, 0, sizeof(AUXT2));
                }
            }

            // Computing the contributions of the real space interactions in the velocities of each particle

            if (ligaih)
            {
                U0[q * numParticles + i] = mobilidade_self[0][0] * FT0[q * numParticles + i] + reduce(H10[0], H10[numRecBoxes], H10[0]);
                U1[q * numParticles + i] = mobilidade_self[1][1] * FT1[q * numParticles + i] + reduce(H11[0], H11[numRecBoxes], H11[0]);
                U2[q * numParticles + i] = mobilidade_self[2][2] * FT2[q * numParticles + i] + reduce(H12[0], H12[numRecBoxes], H12[0]);
            }

            // Computing the contributions of the real space interactions in the magnetic forces and torques acting on each particle

            if (fmagper)
            {
                FORCAS40[q * numParticles + i] = reduce(forcareal0[0], forcareal0[numRecBoxes], forcareal0[0]);
                FORCAS41[q * numParticles + i] = reduce(forcareal1[0], forcareal1[numRecBoxes], forcareal1[0]);
                FORCAS42[q * numParticles + i] = reduce(forcareal2[0], forcareal2[numRecBoxes], forcareal2[0]);
            }

            if (tmagper)
            {
                TORQUES10[q * numParticles + i] = reduce(torquereal0[0], torquereal0[numRecBoxes], torquereal0[0]);
                TORQUES11[q * numParticles + i] = reduce(torquereal1[0], torquereal1[numRecBoxes], torquereal1[0]);
                TORQUES12[q * numParticles + i] = reduce(torquereal2[0], torquereal2[numRecBoxes], torquereal2[0]);
            }

            contribuicao_self[q][i] = mobilidade_self[2][2] * FT2[q * numParticles + i];
            contribuicao_fisico[q][i] = reduce(H12[0], H12[numBoxes], H12[0]);
        }

        // ************************************* SUMS IN THE RECIPROCAL SPACE ***********************************!
        for (q == 0; q < numRealizations; q++)
        {
            for (i == 0; i < numParticles; i++)
            {
                for (s == 0; s < numBoxes; s++)
                {
                    for (j == 0; j < numParticles; j++)
                    {

                        aux_periodico = abs(ILR0[s]) + abs(ILR1[s]) + abs(ILR2[s]);

                        if (aux_periodico != 0.0)
                        {

                            rij[0] = X0[q * numParticles + i] - X0[q * numParticles + j];
                            rij[1] = X1[q * numParticles + i] - X1[q * numParticles + j];
                            rij[2] = X2[q * numParticles + i] - X2[q * numParticles + j];

                            // Calculating the normalized "r" vector = r/|r|

                            modrij = pow(pow(rij[0], 2.0) + pow(rij[1], 2.0) + pow(rij[2], 2.0), 0.5);

                            rn[0] = rij[0] / modrij;
                            rn[1] = rij[1] / modrij;
                            rn[2] = rij[2] / modrij;

                            if (modrij != 0)
                            {

                                // Computing the wave number vector based on the lattice index

                                konda[0] = ILR0[s] * 2.0 * M_PI / l;
                                konda[1] = ILR1[s] * 2.0 * M_PI / l;
                                konda[2] = ILR2[s] * 2.0 * M_PI / h;

                                // Calculating the wave number module

                                modk = pow(pow(konda[0], 2.0) + pow(konda[1], 2.0) + pow(konda[2], 2.0), 0.5);

                                // Normalized wave number vector

                                knormal[0] = konda[0] / modk;
                                knormal[1] = konda[1] / modk;
                                knormal[2] = konda[2] / modk;

                                // ********************** HYDRODYNAMIC INTERACTIONS *********************************************************!

                                if (ligaih)
                                {

                                    // Interpolating the pre-calculated Green function

                                    // call interpola_reciproco(nbr,cof3,cof3,coeficiente3,modk,l)

                                    // Calculating the mobility matrix for the reciprocal space contribution

                                    mobilidade2[0][0] = coeficiente3 * (1.0 - knormal[0] * knormal[0]);
                                    mobilidade2[0][1] = coeficiente3 * (-knormal[0] * knormal[1]);
                                    mobilidade2[0][2] = coeficiente3 * (-knormal[0] * knormal[2]);
                                    mobilidade2[1][0] = coeficiente3 * (-knormal[1] * knormal[0]);
                                    mobilidade2[1][1] = coeficiente3 * (1.0 - knormal[1] * knormal[1]);
                                    mobilidade2[1][2] = coeficiente3 * (-knormal[1] * knormal[2]);
                                    mobilidade2[2][0] = coeficiente3 * (-knormal[2] * knormal[0]);
                                    mobilidade2[2][1] = coeficiente3 * (-knormal[2] * knormal[1]);
                                    mobilidade2[2][2] = coeficiente3 * (1.0 - knormal[2] * knormal[2]);
                                    kr = (knormal[0] * rn[0]) + (knormal[1] * rn[1]) + (knormal[2] * rn[2]);

                                    kr2 = cos((konda[0] * rij[0]) + (konda[1] * rij[1]) + (konda[2] * rij[2]));

                                    HAUX20[j] = kr * (mobilidade2[0][0] * FT0[q * numParticles + j] +
                                                      mobilidade2[0][1] * FT1[q * numParticles + j] +
                                                      mobilidade2[0][2] * FT2[q * numParticles + j]);
                                    HAUX21[j] = kr * (mobilidade2[1][0] * FT0[q * numParticles + j] +
                                                      mobilidade2[1][1] * FT1[q * numParticles + j] +
                                                      mobilidade2[1][2] * FT2[q * numParticles + j]);
                                    HAUX22[j] = kr * (mobilidade2[2][0] * FT0[q * numParticles + j] +
                                                      mobilidade2[2][1] * FT1[q * numParticles + j] +
                                                      mobilidade2[2][2] * FT2[q * numParticles + j]);
                                } //! closing the "if" of the condition (ligaih)

                                //!**************************** MAGNETIC TORQUES ***********************************************!

                                if (tmagper)
                                {

                                    kr2 = 1.0 * ((konda[0] * rij[0]) + (konda[1] * rij[1]) + (konda[2] * rij[2]));

                                    // Interpolating the pre-calculated Green function

                                    // call interpola_reciproco(nbr,cof3,cof7,coeficiente7,modk,l);

                                    // Computing the reciprocal space contribution
                                    termo2 = ((DI0[q * numParticles + j] * konda[0]) + (DI1[q * numParticles + j] * konda[1]) + (DI2[q * numParticles + j] * konda[2]));

                                    // !if(brownianmpecletnum==0.0){
                                    // !lambda=alpha2*8.0
                                    // !else
                                    // !lambda=alpha2*8.0/Per
                                    // !}

                                    if (gravidade)
                                    {
                                        lambda = alpha2 * 6.0 * M_PI / brownianmpecletnum;
                                    }
                                    else
                                    {
                                        lambda = alpha2 * 6.0 * M_PI;
                                    }

                                    termo1 = (DI1[q * numParticles + i] * konda[2]) - (DI2[q * numParticles + i] * konda[1]);
                                    AUXT0[j] = lambda * (coeficiente7 * (termo1 * termo2) * cos(kr2)); // daria problema eliminando parenteses?
                                    termo1 = (DI2[q * numParticles + i] * konda[0]) - (DI0[q * numParticles + i] * konda[2]);
                                    AUXT1[j] = lambda * (coeficiente7 * (termo1 * termo2) * cos(kr2));
                                    termo1 = (DI0[q * numParticles + i] * konda[1]) - (DI1[q * numParticles + i] * konda[0]);
                                    AUXT2[j] = lambda * (coeficiente7 * (termo1 * termo2) * cos(kr2));
                                } //! closing the "if" of the condition (tmagper)

                                // ********************************* MAGNETIC FORCES **********************************************!
                                if (fmagper)
                                {

                                    // Interpolating the pre-calculated Green function

                                    // call interpola_reciproco(nbr,cof3,cof8,coeficiente8,modk,l)

                                    // Computing the reciprocal space contribution

                                    kr = (konda[0] * rij[0]) + (konda[1] * rij[1]) + (konda[2] * rij[2]);
                                    termo1 = ((DI0[q * numParticles + j] * konda[0]) + (DI1[q * numParticles + j] * konda[1]) + (DI2[q * numParticles + j] * konda[2]));
                                    termo2 = ((DI0[q * numParticles + j] * konda[0]) + (DI1[q * numParticles + j] * konda[1]) + (DI2[q * numParticles + j] * konda[2]));

                                    if (brownianmpecletnum == 0.0)
                                    {
                                        lambda = alpha2 * 8.0;
                                    }
                                    else
                                    {
                                        lambda = alpha2 * 8.0 / brownianmpecletnum;
                                    }

                                    double auxResult = lambda * (coeficiente8 * termo1 * termo2 * sin(2.0 * M_PI * kr / l));
                                    AUXF0[j] = auxResult * konda[0];
                                    AUXF1[j] = auxResult * konda[1];
                                    AUXF2[j] = auxResult * konda[2];
                                } // closing the "if" of the condition (fmagper)

                            } // closing the "if" of the condition (modrij!=0)
                        } // closing the "if" of the condition (auxperiodico!=0)
                    }

                    // Making the final sum for the physical and reciprocal contribution for all possible interactions

                    if (ligaih)
                    {
                        H20[s] = reduce(HAUX20[0], HAUX20[numParticles], HAUX20[0]);
                        H21[s] = reduce(HAUX21[0], HAUX21[numParticles], HAUX21[0]);
                        H22[s] = reduce(HAUX22[0], HAUX22[numParticles], HAUX22[0]);
                    }

                    if (fmagper)
                    {
                        forcareciproca[s][0] = reduce(AUXF0[0], AUXF0[numRecBoxes], AUXF0[0]);
                        forcareciproca[s][1] = reduce(AUXF1[0], AUXF1[numRecBoxes], AUXF1[0]);
                        forcareciproca[s][2] = reduce(AUXF2[0], AUXF2[numRecBoxes], AUXF2[0]);
                    }

                    if (tmagper)
                    {
                        torquereciproco[s][0] = reduce(AUXT0[0], AUXT0[numRecBoxes], AUXT0[0]);
                        torquereciproco[s][1] = reduce(AUXT1[0], AUXT1[numRecBoxes], AUXT1[0]);
                        torquereciproco[s][2] = reduce(AUXT2[0], AUXT2[numRecBoxes], AUXT2[0]);
                    }
                }

                if (ligaih)
                {
                    // Computing the velocities of the particles, now with the reciprocal space sum contribution
                    U0[q * numParticles + i] += reduce(H20[0], H20[numRecBoxes], H20[0]);
                    U1[q * numParticles + i] += reduce(H21[0], H21[numRecBoxes], H21[0]);
                    U2[q * numParticles + i] += reduce(H22[0], H22[numRecBoxes], H22[0]);
                    contribuicao_reciproco[q][i] = reduce(H22[0], H22[numRecBoxes], H22[0]);
                }

                if (fmagper)
                {
                    // Computing the magnetic forces acting on the particles, now with the reciprocal space sum contribution
                    FORCAS40[q * numParticles + i] += reduce(FREC0[0], FREC0[numRecBoxes], FREC0[0]);
                    FORCAS41[q * numParticles + i] += reduce(FREC1[0], FREC1[numRecBoxes], FREC1[0]);
                    FORCAS42[q * numParticles + i] += reduce(FREC2[0], FREC2[numRecBoxes], FREC2[0]);
                }

                if (tmagper)
                {
                    // Computing the magnetic torques acting on the particles, now with the reciprocal space sum contribution
                    // reduce(aux0[j * numParticles], aux0[j * numParticles + (numParticles -1)], 0);
                    TORQUES10[q * numParticles + i] += reduce(TREC0[0], TREC0[numRecBoxes], TREC0[0]);
                    TORQUES11[q * numParticles + i] += reduce(TREC1[0], TREC1[numRecBoxes], TREC1[0]);
                    TORQUES12[q * numParticles + i] += reduce(TREC2[0], TREC2[numRecBoxes], TREC2[0]);
                }
            }
        }
    }
}