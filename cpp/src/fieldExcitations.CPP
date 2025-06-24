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
#include <header/fieldExcitations.hpp>
#include <math.h>
#include <iostream>
#include <iomanip>
!*************************************************!
! Subroutine resonsible for pre-calculating the	  !
! magnetic field excitation and shear-rate ramp in!
! case the user wants to build a bifurcation dia- ! 
! gram of the magnetization response		  !
!						  !
! The following possibilities are considered for  !
! the magnetic field H(t)			  !	
!						  !	
! 1 - H(t) as a solution of the Duffing oscilator !
! 2 - H(t) = cos(w1*t) + cos(w2*t) 		  !
! 3 - H(t) = sin(w*t) [oscillatory field]	  !
! 4 - H(t) = sin(i*w*t) [i varying periodically]  !
! 5 - H(t) = sin(w*t)ex + cos(w*t)ey [rotating]   !
!					          !
! Field excitation number 4 implies that the field!
! frequency will periodically change as the simula!
! tion evolves. This is useful to build bifurca-  !
! tion diagrams. 				  !
! 						  !
! This subroutine also produces an oscillatory    !
! shear ramp.				   	  !
!*************************************************!
void fieldExcitations(){   
    
    if(oscilacampo) then    
    
    ! Condições iniciais para o campo H e sua derivada H'=y
    
    y(1)=freqcampo
    campo(1)=0.0
    tempototal(1)=0.0
    
    ! Calculando o intervalo de tempo de atuação de cada frequência (para construção dos diagramas de bifurcação)
    
    intervalo=npast*dt/nfreq
    multiplofreq=0
    
    do k=2,npast
    ! MONTANDO O CAMPO MAGNÉTICO APLICADO CONSIDERANDO A MUDANÇA DA FREQUÊNCIA A CADA PERÍODO PRÉ-DETERMINADO PELO USUÁRIO    
    
    tempototal(k) = tempototal(k-1) + dt 
    campo(k) = sin(freqcampo*tempototal(k)) 
    y(k) = freqcampo*cos(freqcampo*tempototal(k))    
    end do
else{

campo(k)=1;
}

}

}