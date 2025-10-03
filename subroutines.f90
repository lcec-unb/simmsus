!*************************************************!
! 		     SIMMSUS			  ! 
!MODULE:subroutines				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

module subroutines
use variables
contains

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: allocatevariables			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine allocatevariables
  use iso_fortran_env, only: error_unit
  implicit none
  integer :: istat

  if (bifshear) then
    if (allocated(gpvetor)) deallocate(gpvetor, stat=istat)
    allocate(gpvetor(npast), stat=istat); if (istat/=0) goto 900
  end if

  if (allocated(trap)) deallocate(trap, stat=istat)
  allocate(trap(npast), stat=istat); if (istat/=0) goto 900

  if (allocated(X)) deallocate(X, stat=istat)
  allocate(X(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(U)) deallocate(U, stat=istat)
  allocate(U(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(W)) deallocate(W, stat=istat)
  allocate(W(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(FORCAS)) deallocate(FORCAS, stat=istat)
  allocate(FORCAS(6,rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(TORQUES)) deallocate(TORQUES, stat=istat)
  allocate(TORQUES(3,rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(FT)) deallocate(FT, stat=istat)
  allocate(FT(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(Tt)) deallocate(Tt, stat=istat)
  allocate(Tt(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(Di)) deallocate(Di, stat=istat)
  allocate(Di(rea,N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(aux1)) deallocate(aux1, stat=istat)
  allocate(aux1(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(aux2)) deallocate(aux2, stat=istat)
  allocate(aux2(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(aux3)) deallocate(aux3, stat=istat)
  allocate(aux3(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(aux4)) deallocate(aux4, stat=istat)
  allocate(aux4(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(nr)) deallocate(nr, stat=istat)
  allocate(nr(nnr), stat=istat); if (istat/=0) goto 900

  if (allocated(hidrodinamica_aux1)) deallocate(hidrodinamica_aux1, stat=istat)
  allocate(hidrodinamica_aux1(N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(hidrodinamica_aux2)) deallocate(hidrodinamica_aux2, stat=istat)
  allocate(hidrodinamica_aux2(N,3), stat=istat); if (istat/=0) goto 900

  if (allocated(contribuicao_self)) deallocate(contribuicao_self, stat=istat)
  allocate(contribuicao_self(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(contribuicao_fisico)) deallocate(contribuicao_fisico, stat=istat)
  allocate(contribuicao_fisico(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(contribuicao_reciproco)) deallocate(contribuicao_reciproco, stat=istat)
  allocate(contribuicao_reciproco(rea,N), stat=istat); if (istat/=0) goto 900

  if (tmagper) then
    if (allocated(auxt)) deallocate(auxt, stat=istat)
    allocate(auxt(N,3), stat=istat); if (istat/=0) goto 900
  end if

  if (fmagper) then
    if (allocated(auxf)) deallocate(auxf, stat=istat)
    allocate(auxf(N,3), stat=istat); if (istat/=0) goto 900
  end if

  if (allocated(XI)) deallocate(XI, stat=istat)
  allocate(XI(nb,rea,N,3), stat=istat); if (istat/=0) goto 900

  if (leito) then
    if (allocated(usistema)) deallocate(usistema, stat=istat)
    allocate(usistema(rea,3), stat=istat); if (istat/=0) goto 900
  end if

  if (grafmag) then
    if (allocated(magtempo)) deallocate(magtempo, stat=istat)
    allocate(magtempo(3,npast), stat=istat); if (istat/=0) goto 900

    if (allocated(flutmag)) deallocate(flutmag, stat=istat)
    allocate(flutmag(N,rea), stat=istat); if (istat/=0) goto 900
  end if

  if (allocated(tempototal)) deallocate(tempototal, stat=istat)
  allocate(tempototal(npast), stat=istat); if (istat/=0) goto 900

  if (agregado_inicial) then
    if (allocated(centro_massa)) deallocate(centro_massa, stat=istat)
    allocate(centro_massa(rea,3), stat=istat); if (istat/=0) goto 900
  end if

  if (allocated(DIAM)) deallocate(DIAM, stat=istat)
  allocate(DIAM(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(beta)) deallocate(beta, stat=istat)
  allocate(beta(rea,N), stat=istat); if (istat/=0) goto 900

  if (allocated(diarand)) deallocate(diarand, stat=istat)
  allocate(diarand(rea*N), stat=istat); if (istat/=0) goto 900

  if (allocated(campo)) deallocate(campo, stat=istat)
  allocate(campo(npast), stat=istat); if (istat/=0) goto 900

  if (allocated(y)) deallocate(y, stat=istat)
  allocate(y(npast), stat=istat); if (istat/=0) goto 900

  return

900 continue
  write(error_unit,*) 'Allocation has failed in subroutine allocatevariables(). Check your dimensions/entry.'
  stop 2
end subroutine allocatevariables

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: particle_distribution		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine particle_distribution
implicit none
integer i, j

if(polidispersidade)then

call randomica(1.5,2.5,diarand,(N*rea),8)

do j=1,rea
do i=1,N
DIAM(j,i)=diarand((j-1)*N + i)
end do
end do

else

do j=1,rea
do i=1,N
DIAM(j,i)=2.0
end do
end do

end if

do j=1,rea
do i=1,n
beta(j,i) = DIAM(j,i)/DIAM(j,1)
end do
end do

end subroutine particle_distribution

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: box_size				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine box_size
implicit none

if(agregado_inicial) then

ragreg=(N/phi)**(1.0/3.0)
l=100.0*ragreg
h=l

else
 l=((N/(razao*phi))*(4.0*3.1416)/(3.0))**(1.0/3.0)
 razao2=razao
 h=razao2*l
end if

end subroutine box_size

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: gera_arquivos			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating all the fi- !
! les used in the simulation 			  !
!*************************************************!
subroutine gera_arquivos(a,b,c)
implicit none

logical a,b
integer c,i

! Subtitle of the units used for file storage

! From 1 to rea -> posicao.plt
! From rea+1 to 2*rea -> velocidade.plt
! From 2*rea+1 to 3*rea -> dipolo.plt
! From 3*rea+1 to 4*rea -> dip_parcial.plt
! 5*rea -> mag_tempo.plt
! 100*rea-> particula_teste.plt
! 200*rea-> distribuicao_agregados.plt
! 300*rea -> duffing_field.plt
! 1093*j para j=1,rea -> phi_local.plt
! 2012*j para j=1,rea -> energia.plt
! From 400*rea +1 to 400*rea + nfreq


! Creating a file that will be used to plot the trajectory of an arbitrary test particle

open(100*c,file='test_particle.plt')
write(100*c,*) 'Variables="X","Y","Z","Dx","Dy","Dz","t"'


if(bifurcation.or.bifshear) then
do i=1,nfreq
write(rea_char, '(I3)') i
open (400*c+i,file='magnetization'//rea_char//'.plt')
write(400*c+i,*) 'Variables="H","dH/dt","Mx","My","Mz","dMx","dMy","dMz","t"'
end do
end if

if(duffing) then
open(300*c,file='duffing_field.plt')
write(300*c,*)'Variables="t","H(t)","dH/dt"'
end if

if(fator)then
open(6*rea,file='structure_factor.plt')
end if

! Creating the files for each realization

if(a)then
do i=1,c
write(rea_char, '(I3)') i
if(.not.ovito)then
open (i,file='position'//rea_char//'.plt')
else
open (i,file='position'//rea_char//'.xyz')
end if
if(.not.ovito)then
if(gravadipolo) then
write(i,*) 'Variables="X","Y","Z","Dx","Dy","Dz","D"'
else
write(i,*) 'Variables="X","Y","Z","D"'
end if
end if
end do
end if

if(b)then
do i=1,c
write(rea_char, '(I3)') i
open (c+i,file='velocity'//rea_char//'.plt')
write(c+i,*) 'Variables="U","V","W"'
end do
end if

if(grafmag) then
  open(5*c,file='magnetization.plt')
  write(5*c,*)'Variables="H","dH/dt","Mx","My","Mz","dMx","dMy","dMz","t"'
end if

if(printphi) then
do i=1,rea
write(rea_char, '(I3)') i
open (1093*i,file='local_phi'//rea_char//'.plt')
write(1093*i,*) 'Variables="X","Y","Z","Phi"'
end do
end if

! Right now all the necessary files for storaging the position, 
! velocity, magnetization and other important physical variables 
! have been created!

write(*,*) 'All output files have been created'
write(*,*) ''

end subroutine gera_arquivos

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: tabelagreen			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!************************************************!
! Subroutine resposible for pre-calculating the  !
! Green functions used in the calculation of pe- !
! riodic interactions, due to hydrodynamic for-  !
! ces and magnetic forces and torques.           ! 
! This procedure is done to decrease the compu-  !
! tational cost, since theses functions are pre- !
! calculated only once and are simply called any-!
! time a force or torque on a particle is calcu- !
! lated through the code in future loops.        !
!************************************************!

subroutine tabelagreen(a,e,f,g,y)
  implicit none
  ! arguments
  real     a        ! qsi (parâmetro de convergência)
  real     e, y     ! e = l, y = h
  integer  f, g     ! f = nb, g = nbr

  ! local variables
  integer i
  real RR1, RR2, RR3, RR4, RR5, RR7
  real SIG, SIG2, SIG3, SIG4, SIG5, SIG7
  real b(2,10000), c(2,10000), d(2,10000)
  real tb(2,10000), tc(2,10000), td(2,10000), tf(2,10000)
  real te(2,10000)

  ! Making a table with all the coefficients
  do i = 1, 10000

    b(1,i) = 2.0 + (i-1)*(((3.0**0.5)*e*(f**(1.0/3.0)))/9999.0)
    c(1,i) = b(1,i)

    if (tmagper) then
      tb(1,i) = b(1,i)
      tc(1,i) = b(1,i)
    end if
    if (fmagper) then
      td(1,i) = b(1,i)
    end if

    ! Calculating the powers of all the possible distances between the particles
    RR1 = b(1,i)
    RR2 = b(1,i)**2.0
    RR3 = b(1,i)**3.0
    RR4 = b(1,i)**4.0
    RR5 = b(1,i)**5.0
    RR7 = b(1,i)**7.0

    SIG  = qsi
    SIG2 = qsi**2.0
    SIG3 = qsi**3.0
    SIG4 = qsi**4.0
    SIG5 = qsi**5.0
    SIG7 = qsi**7.0

    ! Computing the functions used to calculate hydrodynamic interactions (in real space)
    if (ligaih) then
      b(2,i) = (0.75/RR1 + 0.5/RR3)*ERFC(SIG*RR1) + (4.0*SIG7*RR4 + 3.0*SIG3*RR2 - &
               & 20.0*SIG5*RR2 - 4.50*SIG + 14.0*SIG3 + SIG/RR2)*EXP(-SIG2*RR2)/SQRT(pi)
      c(2,i) = (0.75/RR1 - 1.5/RR3)*ERFC(SIG*RR1) + (-4.0*SIG7*RR4 - 3.0*SIG3*RR2 + &
               & 16.0*SIG5*RR2 + 1.50*SIG - 2.0*SIG3 - 3.0*SIG/RR2)*EXP(-SIG2*RR2)/SQRT(pi)
    end if

    ! Computing the functions used to calculate periodic magnetic torques (in real space)
    if (tmagper) then
      tb(2,i) = (ERFC(SIG*RR1) + ((2.0*SIG*RR1)/SQRT(pi))*EXP(-SIG2*RR2))/RR3
      tc(2,i) = (3.0*ERFC(SIG*RR1) + ((2.0*SIG*RR1)/SQRT(pi))*(3.0 + 2.0*SIG2*RR2)*EXP(-SIG2*RR2))/RR5
    end if

    ! Computing the functions used to calculate periodic magnetic forces (in real space)
    if (fmagper) then
      td(2,i) = (15.0*ERFC(SIG*RR1) + ((2.0*SIG*RR1)/SQRT(pi))*(15.0 + (10.0*SIG2*RR2 + 4.0*SIG4*RR4)) &
                 *EXP(-SIG2*RR2))/RR7
    end if

    ! Same procedure done now for the reciprocal space
    ! Calculating all possible wave numbers
    if (g**(1.0/3.0) .eq. 5.0) then
      d(1,i) = (2.0*pi/e) + (i-1)*(((4.0*(3.0**0.5)*pi/e) - (2.0*pi/e))/9999.0)
    end if
    if (g**(1.0/3.0) .eq. 3.0) then
      d(1,i) = (2.0*pi/e) + (i-1)*(((2.0*(3.0**0.5)*pi/e) - (2.0*pi/e))/9999.0)
    end if

    ! Determining the 2th and 4th powers of these possible wave numbers
    RR2 = d(1,i)**2.0
    RR4 = d(1,i)**4.0

    ! Computing the functions used to calculate hydrodynamic interactions (in reciprocal space)
    if (ligaih) then
      d(2,i) = (1.0/(e*e*y))*(1.0 - RR2/3.0)*(1.0 + 0.25*RR2/SIG2 + 0.125*RR4/SIG4) * &
               6.0*pi*EXP(-0.25*RR2/SIG2)/RR2
    end if

    ! Computing the functions used to calculate periodic magnetic torques (in reciprocal space)
    if (tmagper) then
      te(1,i) = d(1,i)
      te(2,i) = -(1.0/(e*e*y))*(4.0*pi*EXP(-((pi/e)**2.0)*RR2/SIG2))
    end if

    ! Computing the functions used to calculate periodic magnetic forces (in reciprocal space)
    if (fmagper) then
      tf(1,i) = d(1,i)
      tf(2,i) = ((8.0*(pi**2.0))/((e**4.0)*RR2))*EXP(-((pi**2.0)*RR2)/(SIG2*(e**2.0)))/RR2
    end if

    cof1(1,i) = b(1,i)
    cof1(2,i) = b(2,i)
    cof2(1,i) = c(1,i)
    cof2(2,i) = c(2,i)
    cof3(1,i) = d(1,i)
    cof3(2,i) = d(2,i)

    if (tmagper) then
      cof4(1,i) = tb(1,i)
      cof4(2,i) = tb(2,i)
      cof5(1,i) = tc(1,i)
      cof5(2,i) = tc(2,i)
      cof7(1,i) = te(1,i)
      cof7(2,i) = te(2,i)
    end if

    if (fmagper) then
      cof6(1,i) = td(1,i)
      cof6(2,i) = td(2,i)
      cof8(1,i) = tf(1,i)
      cof8(2,i) = tf(2,i)
    end if

  end do

end subroutine tabelagreen


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: estrutura_periodica		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Suroutine responsible for creating a periodic   !
! strucutre (lattice) containing copies of a phy- !
! sical cell with the interacting particles repli-!
! cated through the real and reciprocal spaces    !
! using the Ewald summation technique (1921) des- !
! cribed in details in the works of Beenakker     !
! (1986), Abade (2002), Gontijo (2013)            !
!*************************************************!

subroutine estrutura_periodica(k)
implicit none

integer, intent(in) :: k
integer auxper(5), auxper2(5), auxper3(5)
integer aux, auxk, a, b, c, i, s

aux = nint(nb**(1.0/3.0))
auxk = nint(nbr**(1.0/3.0))

! Creating the lattices indeces in the real space

auxper=0

if(aux == 5)then
auxper(1)=1
auxper(2)=2
auxper(3)=3
auxper(4)=4
auxper(5)=5

auxper2(1)=0
auxper2(2)=5
auxper2(3)=10
auxper2(4)=15
auxper2(5)=20

auxper3(1)=0
auxper3(2)=25
auxper3(3)=50
auxper3(4)=75
auxper3(5)=100
end if

if(aux == 3)then
auxper(1)=1
auxper(2)=2
auxper(3)=3

auxper2(1)=0
auxper2(2)=3
auxper2(3)=6

auxper3(1)=0
auxper3(2)=9
auxper3(3)=18
end if

do a=1,aux
do b=1,aux
do c=1,aux

! Number of physical boxes

s=auxper3(a)+ auxper2(b) +auxper(c)

if(aux == 5) then
ILF(s,1)=a-3
ILF(s,2)=b-3
ILF(s,3)=c-3
end if

if(aux == 3) then
ILF(s,1)=a-2
ILF(s,2)=b-2
ILF(s,3)=c-2
end if

end do
end do
end do

! Creating the initial configuration of all the physical lattices

do a=1,nb
do b=1,rea
do i=1,N

XI(a,b,i,1)= X(b,i,1) + ILF(a,1)*l
XI(a,b,i,2)= X(b,i,2) + ILF(a,2)*l
XI(a,b,i,3)= X(b,i,3) + ILF(a,3)*h

end do
end do
end do

if(k.eq.0)then
open(872,file='typical_initial_lattice.plt')
write(872,*)'Variables= "X1","X2","X3"'
do a=1,nb
do b=1,N
write(872,*) XI(a,1,b,1),XI(a,1,b,2),XI(a,1,b,3)
end do
end do
end if

! Creating the lattice's indeces in the reciprocal space

if(auxk == 5)then
auxper(1)=1
auxper(2)=2
auxper(3)=3
auxper(4)=4
auxper(5)=5

auxper2(1)=0
auxper2(2)=5
auxper2(3)=10
auxper2(4)=15
auxper2(5)=20

auxper3(1)=0
auxper3(2)=25
auxper3(3)=50
auxper3(4)=75
auxper3(5)=100
end if

if(auxk == 3)then
auxper(1)=1
auxper(2)=2
auxper(3)=3

auxper2(1)=0
auxper2(2)=3
auxper2(3)=6

auxper3(1)=0
auxper3(2)=9
auxper3(3)=18
end if

do a=1,auxk
do b=1,auxk
do c=1,auxk

! Number of reciprocal lattices

s=auxper3(a)+ auxper2(b) +auxper(c)

if(auxk == 5) then
ILR(s,1)=a-3
ILR(s,2)=b-3
ILR(s,3)=c-3
end if

if(auxk == 3) then
ILR(s,1)=a-2
ILR(s,2)=b-2
ILR(s,3)=c-2
end if

end do
end do
end do

end subroutine estrutura_periodica

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: dipole_distribution		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating an initial  !
! distribution of the dipole moments of all parti-!
! cles in all realizations			  !
!*************************************************!

subroutine dipole_distribution(a,b,c)
  implicit none
  ! arguments
  integer b, c             ! b = number of realizations, c = number of particles
  real    a(b,c,3)         ! Dipoles

  ! variáveis locais
  real    d, e
  integer f
  integer i, j

  ! If dipoles are distributed in an ordered way
  if (dipolo_ordenado) then

    do j = 1, rea
      do i = 1, N
        a(j,i,1) = 0.0
        a(j,i,2) = 1.0
        a(j,i,3) = 0.0
      end do
    end do

  else
    ! For a random dipole distribution
    e = percentual*c
    f = e

    call randomica(-1.0, 1.0, nr, (3*N*rea), 3)

    ! If we are mixing magnetic particles with non magnetic ones...
    if (mistura) then
      !$OMP PARALLEL DO
      do j = 1, b
        do i = 1, f
          a(j,i,1) = 0.0
          a(j,i,2) = 0.0
          a(j,i,3) = 0.0
        end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO
      do j = 1, b
        do i = f+1, c
          a(j,i,1) = nr((i*2+(i-2)+(c*3*(j-1))))
          a(j,i,2) = nr((i*2+(i-1)+(c*3*(j-1))))
          a(j,i,3) = nr((i*2+(i)+(c*3*(j-1))))
        end do
      end do
      !$OMP END PARALLEL DO

      ! Normalizing the vectors
      !$OMP PARALLEL DO
      do j = 1, b
        do i = 1, f
          a(j,i,1) = 0.0
          a(j,i,2) = 0.0
          a(j,i,3) = 0.0
        end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO
      do j = 1, b
        do i = f+1, N
          d = ((a(j,i,1)**2.0) + (a(j,i,2)**2.0) + (a(j,i,3)**2.0))**0.5
          a(j,i,1) = a(j,i,1)/d
          a(j,i,2) = a(j,i,2)/d
          a(j,i,3) = a(j,i,3)/d
        end do
      end do
      !$OMP END PARALLEL DO

      Di = a

    else
      ! If all the particles are magnetic particles, then...
      !$OMP PARALLEL DO
      do j = 1, b
        do i = 1, c
          a(j,i,1) = nr((i*2+(i-2)+(c*3*(j-1))))
          a(j,i,2) = nr((i*2+(i-1)+(c*3*(j-1))))
          a(j,i,3) = nr((i*2+(i)+(c*3*(j-1))))
        end do
      end do
      !$OMP END PARALLEL DO

      ! Normalizing the vectors
      !$OMP PARALLEL DO
      do j = 1, rea
        do i = 1, N
          d = ((a(j,i,1)**2.0) + (a(j,i,2)**2.0) + (a(j,i,3)**2.0))**0.5
          a(j,i,1) = a(j,i,1)/d
          a(j,i,2) = a(j,i,2)/d
          a(j,i,3) = a(j,i,3)/d
        end do
      end do
      !$OMP END PARALLEL DO

      Di = a
    end if
  end if

end subroutine dipole_distribution


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: initial_condition			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating the initial !
! particle distribution for all simultaneous 	  !
! numerical experiments.			  !
!						  !
! Several possible initial configurations are im- !
! plemented here. These include:		  !
!						  !
! 1 - An initial random distribution		  !
! 2 - An ordered NxNxN array of particles 	  !
! 3 - An initial spherical aggregate		  !
!*************************************************!

subroutine initial_condition
  implicit none
  
  ! local variables
  integer :: auxiliar1, loop, loop2
  integer :: i, j, k
  integer :: a, b, c
  integer, parameter :: MAX_KICKS = 200
  integer :: tries
  real    :: r, reale, reale2
  real    :: nr1, nr2, nr3, modrand, scalefac
  real    :: eps

  eps = max(1.0e-12, tiny(1.0))

  ! If you have to continue a previous simulation, then get into this "if clause"
  if (continua) then

    open(100*rea, file='particula_teste.plt', STATUS='OLD')

509 FORMAT(F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)

    do i = 1, rea
      write(rea_char, '(I3)') i
      open (rea+i, file='velocidade'//rea_char//'.plt', STATUS='OLD')
    end do

    do j = 1, rea
      read (rea+j,'(A)') linha1
      do k = 1, (iter/n2)
        read(rea+j,'(A)') linha2
        do i = 1, N
          read(rea+j,509) U(j,i,1), U(j,i,2), U(j,i,3)
        end do
      end do
    end do

    do i = 1, rea
      write(rea_char, '(I3)') i
      open (i, file='posicao'//rea_char//'.plt', STATUS='OLD')
    end do

    do j = 1, rea
      read (j,'(A)') linha1
      do k = 1, (iter/n2)
        read(j,'(A)') linha2
        do i = 1, N
          read(j,666) X(j,i,1), X(j,i,2), X(j,i,3), Di(j,i,1), Di(j,i,2), Di(j,i,3), DIAM(j,i)
        end do
      end do
    end do

  else

    if (agregado_inicial) then

      open(666, file='vel_tempo.plt')
      write(666,*) 'Variables="t","Ux","Uy","Uz"'

      ragreg = (N/phi)**(1.0/3.0)
      l = 100.0*ragreg
      h = l

      call randomica(-1.0, 1.0, nr, (3*N*rea), 5)

      xcentro = l/2.0
      ycentro = l/2.0
      zcentro = h - 2.0*ragreg

      do j = 1, rea
        do i = 1, N
          xmin = l/2.0 - ragreg
          xmax = xmin + 2.0*ragreg

          X(j,i,1) = xcentro + (xmax-xmin)*0.5*nr((i*2+(i-2)+(N*3*(j-1))))

          ymin = ycentro - (ragreg**2.0 - (X(j,i,1)-xcentro)**2.0)**0.5
          ymax = ycentro + (ragreg**2.0 - (X(j,i,1)-xcentro)**2.0)**0.5

          X(j,i,2) = ycentro + (ymax-ymin)*0.5*nr((i*2+(i-1)+(N*3*(j-1))))

          zmin = zcentro - (ragreg**2.0 - (X(j,i,1)-xcentro)**2.0 - (X(j,i,2)-ycentro)**2.0)**0.5
          zmax = zcentro + (ragreg**2.0 - (X(j,i,1)-xcentro)**2.0 - (X(j,i,2)-ycentro)**2.0)**0.5

          X(j,i,3) = zcentro + (zmax-zmin)*0.5*nr((i*2+(i)+(N*3*(j-1))))
        end do
      end do

      write(*,*) 'Aggregate radius =', ragreg
      write(*,*) 'Box length =', l
      write(*,*) 'Box height =', h

      call randomica(-1.0, 1.0, nr, (3*N*rea), 6)
      
      do k = 1, rea	
      centro_massa(k,1) = sum(X(k,:,1))/N
      centro_massa(k,2) = sum(X(k,:,2))/N
      centro_massa(k,3) = sum(X(k,:,3))/N
      end do
      
    else

      if (ordenado) then
        !$OMP PARALLEL DO
        do j = 1, rea
          loop  = 0
          loop2 = 0
          do i = 1, N
            reale    = (i-1.0)/(N**(1.0/3.0))
            inteiro  = (i-1)/(N**(1.0/3.0))
            reale2   = (loop+1.0)/(N**(1.0/3.0))
            inteiro2 = (loop+1)/(N**(1.0/3.0))

            if (reale  .ne. 0.0 .and. reale  .eq. inteiro ) loop  = loop  + 1
            if (reale2 .ne. 1.0 .and. reale2 .eq. inteiro2) loop2 = loop2 + 1

            auxiliar1 = i/(N**(2.0/3.0))

            a = i - loop*(N**(1.0/3.0))
            b = 1 + loop - auxiliar1*(N**(1.0/3.0))
            c = i/(N**(2.0/3.0))

            X(j,i,1) = l/(2.0*(N**(1.0/3.0))) + a*(l - (l/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)
            X(j,i,2) = l/(2.0*(N**(1.0/3.0))) + b*(l - (l/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)
            X(j,i,3) = h/(2.0*(N**(1.0/3.0))) + c*(h - (h/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)
          end do
        end do
        !$OMP END PARALLEL DO

      else
        call randomica(0.0, 1.0, nr, (3*N*rea), 1)

        !$OMP PARALLEL DO
        do j = 1, rea
          do i = 1, N
            X(j,i,1) = (l) * nr((i*2+(i-2)+(N*3*(j-1))))
            X(j,i,2) = (l) * nr((i*2+(i-1)+(N*3*(j-1))))
            X(j,i,3) = (h) * nr((i*2+(i)+(N*3*(j-1))))
          end do
        end do
        !$OMP END PARALLEL DO

        ! Checking overlaps
        !$OMP PARALLEL DO
        do k = 1, rea
          do i = 1, N
            do j = 1, N
              if (i .ne. j) then
                tries = 0
                do
                  r = ((X(k,i,1)-X(k,j,1))**2.0 + (X(k,i,2)-X(k,j,2))**2.0 + (X(k,i,3)-X(k,j,3))**2.0)**0.5
                  if (r .gt. 2.01) exit

                  ! new brownian "kick" 
                  call randomica(-1.0, 1.0, nr, 3, 2)
                  nr1 = nr(1); nr2 = nr(2); nr3 = nr(3)
                  modrand = (nr1**2.0 + nr2**2.0 + nr3**2.0)**0.5

                  if (modrand .le. eps) then
                    ! small kick
                    nr1 = 0.25; nr2 = 0.0; nr3 = 0.0
                  else
                    scalefac = 0.25/modrand
                    nr1 = nr1*scalefac; nr2 = nr2*scalefac; nr3 = nr3*scalefac
                  end if

                  X(k,i,1) = X(k,i,1) + nr1
                  X(k,i,2) = X(k,i,2) + nr2
                  X(k,i,3) = X(k,i,3) + nr3

                  ! keep inside the box
                  if (X(k,i,1) .lt. 0.0) X(k,i,1) = abs(X(k,i,1))
                  if (X(k,i,2) .lt. 0.0) X(k,i,2) = abs(X(k,i,2))
                  if (X(k,i,3) .lt. 0.0) X(k,i,3) = abs(X(k,i,3))

                  if (X(k,i,1) .gt. l) X(k,i,1) = X(k,i,1) - l
                  if (X(k,i,2) .gt. l) X(k,i,2) = X(k,i,2) - l
                  if (X(k,i,3) .gt. h) X(k,i,3) = X(k,i,3) - h

                  tries = tries + 1
                  if (tries >= MAX_KICKS) exit   ! avoid infinity loop
                end do
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        !End of overlap checking
      end if

      !$OMP PARALLEL DO
      do j = 1, rea
        do i = 1, N
          U(j,i,1) = 0.0
          U(j,i,2) = 0.0
          U(j,i,3) = -1.0
        end do
      end do
      !$OMP END PARALLEL DO
    end if

  end if

end subroutine initial_condition



!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: field_excitations			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

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

subroutine field_excitations
  implicit none
  integer  k

510 FORMAT(F30.4,F30.4,F30.4)

  if (duffing) then
    ! Initial conditions for field H and its derivative H'=y
    y(1)         = 0.0
    campo(1)     = 1.0
    tempototal(1)= 0.0

    write(300*rea,510) tempototal(1), campo(1), y(1)

    ! 4th order Runge–Kutta to solve H(t) and H'(t)
    do k = 2, npast
      k1 = C4*cos(freqcampo*tempototal(k-1)) - (C1*y(k-1) + C2*campo(k-1) + C3*(campo(k-1)**3.0))
      k2 = C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt))-(C1*(y(k-1)+ 0.5*dt*k1)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
      k3 = C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt))-(C1*(y(k-1)+ 0.5*dt*k2)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
      k4 = C4*cos(freqcampo*(tempototal(k-1)+dt))     -(C1*(y(k-1)+ dt*k3)   +C2*campo(k-1)+C3*(campo(k-1)**3.0))

      g1 = y(k-1)
      g2 = y(k-1) + dt*0.5*g1
      g3 = y(k-1) + dt*0.5*g2
      g4 = y(k-1) + dt*g3

      y(k)      = y(k-1)      + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
      campo(k)  = campo(k-1)  + (dt/6.0)*(g1 + 2.0*g2 + 2.0*g3 + g4)
      tempototal(k) = tempototal(k-1) + dt 

      write(300*rea,510) tempototal(k), campo(k), y(k)
    end do

    print *,'******************************************************************************'
    print *,'*                                                                            *'
    print *,'*             DUFFING HARMONIC EXCITATION SUCCESSFULLY PRE-CALCULATED        *'
    print *,'*                                                                            *'
    print *,'******************************************************************************'
    print *,''

  end if

  ! Double frequency excitation
  if (beating) then
    y(1)         = 0.0
    campo(1)     = 2.0
    tempototal(1)= 0.0

    do k = 2, npast
      tempototal(k) = tempototal(k-1) + dt 
      campo(k) = cos(freqcampo*tempototal(k)) + cos(freqbeat*tempototal(k))
      y(k)     = -freqcampo*sin(freqcampo*tempototal(k)) - freqbeat*sin(freqbeat*tempototal(k))
    end do

    print *,'******************************************************************************'
    print *,'*                                                                            *'
    print *,'*           BEATING PATTERN FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED     *'
    print *,'*                                                                            *'
    print *,'******************************************************************************'
    print *,''

  end if

  ! Oscillating field (not duffing nor beating (double-frequency))
  if (oscilacampo) then
    if (.not. duffing) then
      if (.not. beating) then
        y(1)         = freqcampo
        campo(1)     = 0.0
        tempototal(1)= 0.0

        intervalo    = npast*dt/nfreq
        multiplofreq = 0

        do k = 2, npast
          if (bifurcation) then
            tempototal(k) = tempototal(k-1) + dt 
            contfreqinteiro1 = tempototal(k-1)/intervalo
            contfreqinteiro2 = tempototal(k)/intervalo
            if (contfreqinteiro1 .ne. contfreqinteiro2) multiplofreq = multiplofreq + 1
            campo(k) = sin((freqcampo + (bifmax - freqcampo)*multiplofreq)*tempototal(k))
            y(k)     = (freqcampo + (bifmax - freqcampo)*multiplofreq)* &
                        cos((freqcampo + (bifmax - freqcampo)*multiplofreq)*tempototal(k))
          else
            tempototal(k) = tempototal(k-1) + dt 
            campo(k) = sin(freqcampo*tempototal(k)) 
            y(k)     = freqcampo*cos(freqcampo*tempototal(k))
          end if
        end do
        multiplofreq = 0

        print *,'******************************************************************************'
        print *,'*                                                                            *'
        print *,'*                HARMONIC FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED       *'
        print *,'*                                                                            *'
        print *,'******************************************************************************'
        print *,''

      end if
    end if
  end if

  ! Shear ramp (shear rate bifurcation)
  do k = 2, npast
    if (bifshear) then
      tempototal(k) = tempototal(k-1) + dt 
      contfreqinteiro1 = tempototal(k-1)/intervalo
      contfreqinteiro2 = tempototal(k)/intervalo
      if (contfreqinteiro1 .ne. contfreqinteiro2) multiplofreq = multiplofreq + 1
      gpvetor(k) = shearrate + multiplofreq*((bifmax - shearrate)/nfreq) 
    end if
  end do
  multiplofreq = 0

end subroutine field_excitations
 

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: campo_phi				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating a tecplot   !	
! format file with a 3D field of the local volume !
! fraction of particles inside the simulation box !
! domain. This subroutine divides the simulation  !
! box space into 10x10x10 cells, counts the number!
! of particles in each cell and divide by the num-!
! ber of particles within the suspension space.   !
!*************************************************!

subroutine campo_phi(a,s)
  implicit none
  ! argumentos
  integer :: a, s

  ! variáveis locais
  real    :: dx, dy, dz
  integer :: nc, ncubo
  integer :: i, j, k, o, p, auxiliar
  integer, allocatable :: auxper(:), auxper2(:), auxper3(:)
  real,    allocatable :: philocal(:,:)
  real    :: reali, realj, realk, realn

  nc     = 10
  ncubo  = nc**3        ! inteiro**inteiro → inteiro (evita conversão implícita)
  realn  = nc

  allocate(philocal(a,ncubo))
  allocate(auxper(nc))
  allocate(auxper2(nc))
  allocate(auxper3(nc))

  dx = l/realn
  dy = l/realn
  dz = h/realn

  philocal = 0.0

  do i = 1, nc
    auxper(i)  = i
    auxper2(i) = (i-1)*nc
    auxper3(i) = (i-1)*(nc**2)
  end do

  do i = 1, nc
    do j = 1, nc
      do k = 1, nc
        do p = 1, rea
          do o = 1, N
            if (X(p,o,1) .ge. (i-1)*dx .and. X(p,o,1) .le. i*dx) then
              if (X(p,o,2) .ge. (j-1)*dy .and. X(p,o,2) .le. j*dy) then
                if (X(p,o,3) .ge. (k-1)*dz .and. X(p,o,3) .le. k*dz) then
                  auxiliar = auxper3(k) + auxper2(j) + auxper(i)
                  philocal(p,auxiliar) = philocal(p,auxiliar) + 1
                end if
              end if
            end if
          end do
        end do
      end do
    end do
  end do

  do p = 1, rea
    ! (Se quiser, troque "10" por nc no texto:)
    write(1093*p,*) 'zone F=POINT,I=10,J=10,K=10'
    do i = 1, nc
      do j = 1, nc
        do k = 1, nc
          auxiliar = auxper3(k) + auxper2(j) + auxper(i)

          reali = i; realj = j; realk = k
          write(1093*p,*) (reali-1.0)*dx*(1.0 + 1.2/realn), &
                           (realj-1.0)*dy*(1.0 + 1.2/realn), &
                           (realk-1.0)*dz*(1.0 + 1.2/realn),  &
                           philocal(p,auxiliar)
        end do
      end do
    end do
  end do

  deallocate(philocal)
  deallocate(auxper, auxper2, auxper3)
end subroutine campo_phi


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: brownian				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for implementing Brownian!
! forces and torques				  !
!*************************************************!

subroutine brownian(k)
  implicit none
  ! argumentos
  integer, intent(in) :: k

  ! variáveis locais
  integer :: i, j
  real    :: nr1, nr2, nr3, modrand

  call randomica(-1.0, 1.0, nr, (3*N*rea), 3 + k)

  !$OMP PARALLEL DO
  do j = 1, rea
    do i = 1, N

      nr1 = nr((i*2 + (i-2) + (N*3*(j-1))))
      nr2 = nr((i*2 + (i-1) + (N*3*(j-1))))
      nr3 = nr((i*2 + (i)   + (N*3*(j-1))))

      modrand = sqrt(nr1**2.0 + nr2**2.0 + nr3**2.0)
      nr1 = nr1 / modrand
      nr2 = nr2 / modrand
      nr3 = nr3 / modrand

      if (gravidade) then
        FORCAS(6,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr1
        FORCAS(6,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr2
        FORCAS(6,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr3
      else if (shear) then
        FORCAS(6,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr1
        FORCAS(6,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr2
        FORCAS(6,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(6.0/(Pe*dt)) * nr3
      else
        FORCAS(6,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(6.0/dt) * nr1
        FORCAS(6,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(6.0/dt) * nr2
        FORCAS(6,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(6.0/dt) * nr3
      end if

    end do
  end do
  !$OMP END PARALLEL DO

  if (torque) then
    call randomica(-1.0, 1.0, nr, (3*N*rea), npast + k)

    !$OMP PARALLEL DO
    do j = 1, rea
      do i = 1, N
        nr1 = nr((i*2 + (i-2) + (N*3*(j-1))))
        nr2 = nr((i*2 + (i-1) + (N*3*(j-1))))
        nr3 = nr((i*2 + (i)   + (N*3*(j-1))))

        modrand = sqrt(nr1**2.0 + nr2**2.0 + nr3**2.0)
        nr1 = nr1 / modrand
        nr2 = nr2 / modrand
        nr3 = nr3 / modrand

        if (gravidade) then
          TORQUES(3,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr1
          TORQUES(3,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr2
          TORQUES(3,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr3
        else if (shear) then
          TORQUES(3,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr1
          TORQUES(3,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr2
          TORQUES(3,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*Pe*dt)) * nr3
        else
          TORQUES(3,j,i,1) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*dt)) * nr1
          TORQUES(3,j,i,2) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*dt)) * nr2
          TORQUES(3,j,i,3) = (beta(j,i)**(-2.0)) * sqrt(9.0/(2.0*dt)) * nr3
        end if

      end do
    end do
    !$OMP END PARALLEL DO
  end if
end subroutine brownian

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: gravity				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for implementing gravita-!
! tional forces					  !
!*************************************************!
subroutine gravity
  implicit none
  integer :: i, j

  !$OMP PARALLEL DO
  do j = 1, rea
    do i = 1, N
      if (gravidade) then
        ! Gravitational forces
        FORCAS(3,j,i,1) = 0.0
        FORCAS(3,j,i,2) = 0.0
        FORCAS(3,j,i,3) = -beta(j,i)**3.0
      else
        ! Gravitational forces
        FORCAS(3,j,i,1) = 0.0
        FORCAS(3,j,i,2) = 0.0
        FORCAS(3,j,i,3) = 0.0
      end if
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine gravity


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: repulsion				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing repulsive  !
! forces (lubrication) and contact forces for  	  !
! overlapped particles using the Hertz model	  !
!*************************************************!
subroutine repulsion
  implicit none
  integer :: i, j, q
  real    :: r, Eij

  !$OMP PARALLEL DO
  do j = 1, rea
    do i = 1, N
      do q = 1, N
        if (i .ne. q) then
          ! pairwise distance
          r = sqrt( (X(j,i,1)-X(j,q,1))**2.0 + (X(j,i,2)-X(j,q,2))**2.0 + (X(j,i,3)-X(j,q,3))**2.0 )

          ! short-range repulsion
          
          if (r .le. 2.2) then
            if (r .ge. 2.0) then
              FORCAS(1,j,i,1) = 10.0*exp(-r/0.01) * (X(j,i,1) - X(j,q,1)/r)
              FORCAS(1,j,i,2) = 10.0*exp(-r/0.01) * (X(j,i,2) - X(j,q,2)/r)
              FORCAS(1,j,i,3) = 10.0*exp(-r/0.01) * (X(j,i,3) - X(j,q,3)/r)
            end if
          end if
          if (r .ge. 2.2) then
            if (r .le. 2.0) then
              FORCAS(1,j,i,1) = 0.0
              FORCAS(1,j,i,2) = 0.0
              FORCAS(1,j,i,3) = 0.0
            end if
          end if

          ! contact (Hertz) for overlapped particles
          if (r .le. (beta(j,i)+beta(j,q))) then
            Eij = abs(2.0 - r)
            if (abs(Eij) .gt. 1.0e-05) then
              FORCAS(2,j,i,1) = (100.0*Eij**(3.0/2.0)) * (X(j,i,1)-X(j,q,1))/Eij
              FORCAS(2,j,i,2) = (100.0*Eij**(3.0/2.0)) * (X(j,i,2)-X(j,q,2))/Eij
              FORCAS(2,j,i,3) = (100.0*Eij**(3.0/2.0)) * (X(j,i,3)-X(j,q,3))/Eij
            end if
          else
            FORCAS(2,j,i,1) = 0.0
            FORCAS(2,j,i,2) = 0.0
            FORCAS(2,j,i,3) = 0.0
          end if

        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine repulsion


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: forca_magnetica			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing long-range !
! non-periodic dipolar interactions between the	  !
! particles					  !
!*************************************************!

subroutine forca_magnetica
  implicit none
  integer :: auxiliary
  integer :: i, j, q
  real    :: r, termo1, termo2, termo3, termo4

  if (mistura) then
    auxiliary = (percentual*N) + 1
  else
    auxiliary = 1
  end if

  if (gravidade) then
    lambda = alpha2*6.0/Pe
  else if (shear) then
    lambda = alpha2*6.0/Pe
  else
    lambda = alpha2*6.0
  end if

  do j = 1, rea
    do i = auxiliary, N
      do q = auxiliary, N
        if (i .ne. q) then
          ! distância entre partículas
          r = sqrt( (X(j,i,1)-X(j,q,1))**2.0 + (X(j,i,2)-X(j,q,2))**2.0 + (X(j,i,3)-X(j,q,3))**2.0 )

          ! desliga interação magnética quando próximas
          if (r .le. 2.2) then
            aux1(j,q) = 0.0
            aux2(j,q) = 0.0
            aux3(j,q) = 0.0
          else
            ! vetor R_ij
            rij(1) = X(j,i,1) - X(j,q,1)
            rij(2) = X(j,i,2) - X(j,q,2)
            rij(3) = X(j,i,3) - X(j,q,3)
            ! normalização
            modrij = sqrt(rij(1)**2.0 + rij(2)**2.0 + rij(3)**2.0)
            rij(1) = rij(1)/modrij
            rij(2) = rij(2)/modrij
            rij(3) = rij(3)/modrij

            termo1 = (Di(j,i,1)*Di(j,q,1) + Di(j,i,2)*Di(j,q,2) + Di(j,i,3)*Di(j,q,3)) * rij(1)
            termo2 = (Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * Di(j,q,1)
            termo3 = (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3)) * Di(j,i,1)
            termo4 = ((Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * &
                      (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3))) * rij(1)
            aux1(j,q) = (lambda/(r**4.0)) * (termo1 + termo2 + termo3 - 5.0*termo4)

            termo1 = (Di(j,i,1)*Di(j,q,1) + Di(j,i,2)*Di(j,q,2) + Di(j,i,3)*Di(j,q,3)) * rij(2)
            termo2 = (Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * Di(j,q,2)
            termo3 = (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3)) * Di(j,i,2)
            termo4 = ((Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * &
                      (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3))) * rij(2)
            aux2(j,q) = (lambda/(r**4.0)) * (termo1 + termo2 + termo3 - 5.0*termo4)

            termo1 = (Di(j,i,1)*Di(j,q,1) + Di(j,i,2)*Di(j,q,2) + Di(j,i,3)*Di(j,q,3)) * rij(3)
            termo2 = (Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * Di(j,q,3)
            termo3 = (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3)) * Di(j,i,3)
            termo4 = ((Di(j,i,1)*rij(1) + Di(j,i,2)*rij(2) + Di(j,i,3)*rij(3)) * &
                      (Di(j,q,1)*rij(1) + Di(j,q,2)*rij(2) + Di(j,q,3)*rij(3))) * rij(3)
            aux3(j,q) = (lambda/(r**4.0)) * (termo1 + termo2 + termo3 - 5.0*termo4)
          end if
        end if
      end do

      FORCAS(4,j,i,1) = sum(aux1(j,:))
      FORCAS(4,j,i,2) = sum(aux2(j,:))
      FORCAS(4,j,i,3) = sum(aux3(j,:))

      aux1 = 0.0
      aux2 = 0.0
      aux3 = 0.0
    end do
  end do

end subroutine forca_magnetica

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: campo_externo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing magnetic	  !
! forces between an external magnetic field and	  !
! each particle in suspension.  		  !
!*************************************************!

subroutine campo_externo
  implicit none
  integer :: auxiliary
  integer :: i, j

 !$OMP PARALLEL DO
do j=1,rea
do i=auxiliary,N
! posicao_campo = 1 -> Applied field in the lower wall
if(posicao_campo.eq.1)then
if(Pe.eq.0.0) then
if(X(j,i,3).ge.(2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*(2.0**3.0))
end if
else
if(X(j,i,3).ge.(2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*(2.0**3.0))
end if
end if
end if

! posicao_campo = 2 -> Applied field on the upper wall

if(posicao_campo.eq.2) then
if(Pe.eq.0.0) then
if(X(j,i,3).le.(h-2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((h-X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((h-2.0)**3.0))
end if
else
if(X(j,i,3).le.(h-2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((h-X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((h-2.0)**3.0))
end if
end if
end if

! posicao_campo = 3 -> Applied field on the right side

if(posicao_campo.eq.3) then
if(Pe.eq.0.0) then
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(1.0*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(1.0*((l-2.0)**3.0))
end if
else
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(Pe*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(Pe*((l-2.0)**3.0))
end if
end if
end if

! posicao_campo = 4 -> Applied field on the left side

if(posicao_campo.eq.4) then
if(Pe.eq.0.0)then
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(1.0*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(1.0*((l-2.0)**3.0))
end if
else
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(Pe*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(Pe*((l-2.0)**3.0))
end if
end if
end if

end do
end do
 !$OMP END PARALLEL DO

end subroutine campo_externo



!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: periodic_interactions		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Suroutine responsible for computing properly all!
! periodic interactions. 			  !
!						  !
! This subroutine is called whenever the program  !
! identifies the presence of any periodic interac-!
! tion. These periodic interactions can be of the !
! following kinds:				  !
!						  !
! 1 - Long range hydrodynamic interactions	  !
! 2 - Long range periodic magnetic torques        !
! 3 - Long range periodic magnetic forces  	  !
!						  !
! This is probably the most expensive and complex !
! subroutine in this code.			  !
!*************************************************!

subroutine periodic_interactions(k)
  implicit none
  integer, intent(in) :: k

  !------------------- parâmetros de corte (agora independentes) -------------------
  ! cutoff: ativa/desativa os cortes. Se .false., soma completa.
  logical :: cutoff
  ! Y1 controla o corte no real: rcut = Y1/qsi
  ! Y2 controla o corte no recíproco: kcut = 2*qsi*Y2
  real :: Y1, Y2
  real :: rcut, rcut2, kcut
  !---------------------------------------------------------------------------------

  ! variáveis locais
  real    :: mobilidade_self(3,3), mobilidade1(3,3), mobilidade2(3,3)
  real    :: coeficiente1, coeficiente2, coeficiente3, coeficiente4
  real    :: coeficiente5, coeficiente6, coeficiente7, coeficiente8
  real    :: rn(3)
  integer :: aux_periodico
  integer :: q, i, s, j, a, d
  integer, parameter :: dmin = 2, dmax = 9999
  real    :: eps
  real    :: Tx, Ty, Tz, T2, diag_caixa, r2, modrij
  real    :: kr, kr2, modk
  real    :: termo1, termo2, termo3, termo4, termo5
  real    :: konda(3), knormal(3)

  ! valores numéricos robustos
  eps = max(1.0e-12, tiny(1.0))

  ! ----------------- escolha do corte -----------------
  cutoff = .true.   ! <<< ajuste aqui: .true. usa corte; .false. faz soma completa
  Y1     = 5.0      ! <<< varrer real: 4,6,8,10 ...
  Y2     = 50.0      ! <<< varrer recíproco: 3,5,7,9 ...

  if (cutoff) then
     rcut  = Y1 / qsi
     kcut  = 2.0 * qsi * Y2
  else
     rcut  = huge(1.0)
     kcut  = huge(1.0)
  end if
  rcut2 = rcut * rcut
  ! ----------------------------------------------------

  diag_caixa = sqrt(l*l + l*l + h*h)

  ! We start by computing first the sums in the real space 
  ! Here, nb denotes the number of physical boxes. This
  ! number can be different from nbr, which denotes the
  ! nubmer of reciprocal boxes. In simconfig.dat we 
  ! recommend the user to set nb=125 and nbr=27. These
  ! values ensure a good precision when computing the 
  ! average sedimentation velocity of a suspension of
  ! spheres in Creeping-flow

  ! This loop works like this: we fix a realization and
  ! a given particle and the make a sweep computing the
  ! long-range interactions between this particles and
  ! all the other particles in the surrounding boxes.
  ! We then, extend this procedure to all real particles
  ! in all the simultaneous numerical experiments.

  do q=1,rea
    do i=1,N

      ! Building the self-mobility matrix (independe de j,s)
      if (ligaih) then
        mobilidade_self = 0.0
        do a=1,3
          mobilidade_self(a,a)= 1.0 - (6.0*(pi**(-0.5))*qsi) + ((40.0/3.0)*(pi**(-0.5))*(qsi**3.0))
        end do
      end if

      do s=1,nb
        do j=1,N

          aux_periodico = abs(ILF(s,1)) + abs(ILF(s,2)) + abs(ILF(s,3))

          ! Checking the distance between a particle "i" and other particle "j" in the real boxes (physical and images)  
          ! Calculating the "r" vector (distance)

          if (.not. shear) then
            Tx = ILF(s,1)*l
            Ty = ILF(s,2)*l
            Tz = ILF(s,3)*h
          else
            Tx = ILF(s,1)*l
            Ty = (ILF(s,2)*l) + (k*ILF(s,3)*shearrate*dt)
            Tz = ILF(s,3)*h
          end if

          ! culling rápido por caixa imagem (se muito longe, não alcança rcut na pior hipótese)
          T2 = Tx*Tx + Ty*Ty + Tz*Tz
          if (sqrt(T2) - diag_caixa > rcut) cycle

          ! vetor r_ij
          rij(1) = X(q,i,1) - (X(q,j,1) + Tx)
          rij(2) = X(q,i,2) - (X(q,j,2) + Ty)
          rij(3) = X(q,i,3) - (X(q,j,3) + Tz)

          ! Normalizing the "r" vector -> r/|r|
          r2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
          if (r2 <= 4.0) cycle         ! ignora esfera interna (r <= 2a) — igual ao seu >2.0
          if (r2 >  rcut2) cycle       ! corte no real

          modrij = sqrt(r2)

          if (modrij .ge. eps) then
            rn(1) = rij(1)/modrij
            rn(2) = rij(2)/modrij
            rn(3) = rij(3)/modrij
          else
            rn = 0.0
          end if

          d = 1 + int(9999.0*(modrij-2.0)/((3.0**0.5)*l*(nb**(1.0/3.0))))
          if (d < dmin) d = dmin
          if (d > dmax) d = dmax
          diferenca_interpol1 = (modrij - cof1(1,d))

          !************ SUMS IN THE PHYSICAL SPACE IN THE REAL BOX (WHERE THE REAL PARTICLES ARE) AND "IMAGE" BOXES***********!

          !******************************* HYDRODYNAMIC INTERACTIONS *********************************************************!

          ! Building the self-mobility matrix
          if (aux_periodico.eq.0 .and. i.eq.j) then
            ! já calculada fora
          end if

          ! Building the called mobility matrix 1 (lattice sum in the physical space)
          if (ligaih) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente1 = cof1(2,d) + ((modrij-cof1(1,d))/(cof1(1,d+1)-cof1(1,d))) * (cof1(2,d+1)-cof1(2,d))
                coeficiente2 = cof2(2,d) + ((modrij-cof2(1,d))/(cof2(1,d+1)-cof2(1,d))) * (cof2(2,d+1)-cof2(2,d))
              else
                coeficiente1 = cof1(2,d-1) + ((modrij-cof1(1,d-1))/(cof1(1,d)-cof1(1,d-1))) * (cof1(2,d)-cof1(2,d-1))
                coeficiente2 = cof2(2,d-1) + ((modrij-cof2(1,d-1))/(cof2(1,d)-cof2(1,d-1))) * (cof2(2,d)-cof2(2,d-1))
              end if

              mobilidade1(1,1)=coeficiente1 + coeficiente2*(rn(1)*rn(1))
              mobilidade1(1,2)=coeficiente2*(rn(1)*rn(2))
              mobilidade1(1,3)=coeficiente2*(rn(1)*rn(3))

              mobilidade1(2,1)=coeficiente2*(rn(2)*rn(1))
              mobilidade1(2,2)=coeficiente1 + coeficiente2*(rn(2)*rn(2))
              mobilidade1(2,3)=coeficiente2*(rn(2)*rn(3))

              mobilidade1(3,1)=coeficiente2*(rn(3)*rn(1))
              mobilidade1(3,2)=coeficiente2*(rn(3)*rn(2))
              mobilidade1(3,3)=coeficiente1 + coeficiente2*(rn(3)*rn(3))

              hidrodinamica_aux1(j,1)=mobilidade1(1,1)*FT(q,j,1)+mobilidade1(1,2)*FT(q,j,2)+mobilidade1(1,3)*FT(q,j,3)
              hidrodinamica_aux1(j,2)=mobilidade1(2,1)*FT(q,j,1)+mobilidade1(2,2)*FT(q,j,2)+mobilidade1(2,3)*FT(q,j,3)
              hidrodinamica_aux1(j,3)=mobilidade1(3,1)*FT(q,j,1)+mobilidade1(3,2)*FT(q,j,2)+mobilidade1(3,3)*FT(q,j,3)
            end if
          end if

          !****************** MAGNETIC TORQUES **************************************!
          if (tmagper) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente4=cof4(2,d) + ((modrij-cof4(1,d))/(cof4(1,d+1)-cof4(1,d)))*(cof4(2,d+1)-cof4(2,d))
                coeficiente5=cof5(2,d) + ((modrij-cof5(1,d))/(cof5(1,d+1)-cof5(1,d)))*(cof5(2,d+1)-cof5(2,d))
              else
                coeficiente4=cof4(2,d-1) + ((modrij-cof4(1,d-1))/(cof4(1,d)-cof4(1,d-1)))*(cof4(2,d)-cof4(2,d-1))
                coeficiente5=cof5(2,d-1) + ((modrij-cof5(1,d-1))/(cof5(1,d)-cof5(1,d-1)))*(cof5(2,d)-cof5(2,d-1))
              end if

              ! Computing periodic torques due to magnetic interactions in the real space
              termo5 = 0.0
              termo2 = (Di(q,j,1)*rij(1)) + (Di(q,j,2)*rij(2)) + (Di(q,j,3)*rij(3))

              if (gravidade) then
                lambda=alpha2*6.0*pi/Pe
              else
                lambda=alpha2*6.0*pi
              end if

              termo1 = (Di(q,i,3)*Di(q,j,2))-(Di(q,i,2)*Di(q,j,3))
              termo3 = (Di(q,i,2)*rij(3))-(Di(q,i,3)*rij(2))
              auxt(j,1) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)

              termo1 = (Di(q,i,1)*Di(q,j,3))-(Di(q,i,3)*Di(q,j,1))
              termo3 = (Di(q,i,3)*rij(1))-(Di(q,i,1)*rij(3))
              auxt(j,2) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)

              termo1 = (Di(q,i,2)*Di(q,j,1))-(Di(q,i,1)*Di(q,j,2))
              termo3 = (Di(q,i,1)*rij(2))-(Di(q,i,2)*rij(1))
              auxt(j,3) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)
            end if
          end if

          !*********************************** MAGNETIC FORCES *****************************************!
          if (fmagper) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente6=cof6(2,d) + ((modrij-cof6(1,d))/(cof6(1,d+1)-cof6(1,d)))*(cof6(2,d+1)-cof6(2,d))
              else
                coeficiente6=cof6(2,d-1) + ((modrij-cof6(1,d-1))/(cof6(1,d)-cof6(1,d-1)))*(cof6(2,d)-cof6(2,d-1))
              end if

              if (gravidade) then
                lambda=alpha2*8.0*pi/Pe
              else
                lambda=alpha2*8.0*pi
              end if

              ! (di.dj)rij
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(1)
              ! (dj.rij)di
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,1)
              ! (di.rij)dj
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,1)
              ! ((di.rij)(dj.rij))rij
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(1)

              auxf(j,1)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

              ! y
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(2)
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,2)
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,2)
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(2)

              auxf(j,2)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

              ! z
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(3)
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,3)
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,3)
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(3)

              auxf(j,3)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)
            end if
          end if

        end do

        !************************** SUMING EVERYTHING *********************************!

        ! Computing the sum of the previous vectors and matrices
        if (ligaih) then
          hidro1(s,1)=sum(hidrodinamica_aux1(:,1))
          hidro1(s,2)=sum(hidrodinamica_aux1(:,2))
          hidro1(s,3)=sum(hidrodinamica_aux1(:,3))
        end if

        if (fmagper) then
          forcareal(s,1)=sum(auxf(:,1))
          forcareal(s,2)=sum(auxf(:,2))
          forcareal(s,3)=sum(auxf(:,3))
        end if

        if (tmagper) then
          torquereal(s,1)=sum(auxt(:,1))
          torquereal(s,2)=sum(auxt(:,2))
          torquereal(s,3)=sum(auxt(:,3))
          auxt = 0.0
        end if

      end do

      ! Computing the contributions of the real space interactions in the velocities of each particle
      if (ligaih) then
        U(q,i,1)=mobilidade_self(1,1)*FT(q,i,1) + sum(hidro1(:,1))
        U(q,i,2)=mobilidade_self(2,2)*FT(q,i,2) + sum(hidro1(:,2))
        U(q,i,3)=mobilidade_self(3,3)*FT(q,i,3) + sum(hidro1(:,3))
      end if

      hidro1=0.0

      ! Computing the contributions of the real space interactions in the magnetic forces and torques acting on each particle
      if (fmagper) then
        FORCAS(4,q,i,1)=sum(forcareal(:,1))
        FORCAS(4,q,i,2)=sum(forcareal(:,2))
        FORCAS(4,q,i,3)=sum(forcareal(:,3))
      end if

      if (tmagper) then
        TORQUES(1,q,i,1)=sum(torquereal(:,1))
        TORQUES(1,q,i,2)=sum(torquereal(:,2))
        TORQUES(1,q,i,3)=sum(torquereal(:,3))
      end if

      contribuicao_self(q,i)   = mobilidade_self(3,3)*FT(q,i,3)
      contribuicao_fisico(q,i) = sum(hidro1(:,3))

    end do
  end do

  !************************************* SUMS IN THE RECIPROCAL SPACE ***********************************!
  do q=1,rea
    do i=1,N
      do s=1,nbr
        do j=1,N

          aux_periodico = abs(ILR(s,1)) + abs(ILR(s,2)) + abs(ILR(s,3))

          if (aux_periodico.ne.0) then

            rij(1)= X(q,i,1)-X(q,j,1)
            rij(2)= X(q,i,2)-X(q,j,2)
            rij(3)= X(q,i,3)-X(q,j,3)

            ! Calculating the normalized "r" vector = r/|r|
            modrij = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))

            if (modrij .ne. 0.0) then
              rn(1)=rij(1)/modrij
              rn(2)=rij(2)/modrij
              rn(3)=rij(3)/modrij
            else
              rn = 0.0
            end if

            if (modrij.ne.0.0) then

              ! Computing the wave number vector based on the lattice index 
              konda(1)=ILR(s,1)*2.0*pi/l
              konda(2)=ILR(s,2)*2.0*pi/l
              konda(3)=ILR(s,3)*2.0*pi/h

              ! Calculating the wave number module
              modk = sqrt(konda(1)*konda(1) + konda(2)*konda(2) + konda(3)*konda(3))

              ! -------------------- CORTE NO RECÍPROCO (independente) --------------------
              if (modk > kcut) cycle
              ! ---------------------------------------------------------------------------

              ! Normalized wave number vector
              knormal(1)=konda(1)/modk
              knormal(2)=konda(2)/modk
              knormal(3)=konda(3)/modk

              !********************** HYDRODYNAMIC INTERACTIONS *********************************************************!
              if (ligaih) then

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof3,coeficiente3,modk,l)

                ! Calculating the mobility matrix for the reciprocal space contribution
                mobilidade2(1,1)=coeficiente3*(1.0 - knormal(1)*knormal(1))
                mobilidade2(1,2)=coeficiente3*(-knormal(1)*knormal(2))
                mobilidade2(1,3)=coeficiente3*(-knormal(1)*knormal(3))
                mobilidade2(2,1)=coeficiente3*(-knormal(2)*knormal(1))
                mobilidade2(2,2)=coeficiente3*(1.0 - knormal(2)*knormal(2))
                mobilidade2(2,3)=coeficiente3*(- knormal(2)*knormal(3))
                mobilidade2(3,1)=coeficiente3*(-knormal(3)*knormal(1))
                mobilidade2(3,2)=coeficiente3*(-knormal(3)*knormal(2))
                mobilidade2(3,3)=coeficiente3*(1.0 - knormal(3)*knormal(3))

                kr  = (knormal(1)*rn(1)) + (knormal(2)*rn(2)) + (knormal(3)*rn(3))
                kr2 = cos((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

                hidrodinamica_aux2(j,1)=kr*(mobilidade2(1,1)*FT(q,j,1)+mobilidade2(1,2)*FT(q,j,2)+mobilidade2(1,3)*FT(q,j,3))
                hidrodinamica_aux2(j,2)=kr*(mobilidade2(2,1)*FT(q,j,1)+mobilidade2(2,2)*FT(q,j,2)+mobilidade2(2,3)*FT(q,j,3))
                hidrodinamica_aux2(j,3)=kr*(mobilidade2(3,1)*FT(q,j,1)+mobilidade2(3,2)*FT(q,j,2)+mobilidade2(3,3)*FT(q,j,3))
              end if

              !**************************** MAGNETIC TORQUES ***********************************************!
              if (tmagper) then

                kr2 = 1.0*((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof7,coeficiente7,modk,l)

                ! Computing the reciprocal space contribution
                termo2 = (Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3))

                if (gravidade) then
                  lambda=alpha2*6.0*pi/Pe
                else
                  lambda=alpha2*6.0*pi
                end if

                termo1=(Di(q,i,2)*konda(3))-(Di(q,i,3)*konda(2))
                auxt(j,1)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
                termo1=(Di(q,i,3)*konda(1))-(Di(q,i,1)*konda(3))
                auxt(j,2)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
                termo1=(Di(q,i,1)*konda(2))-(Di(q,i,2)*konda(1))
                auxt(j,3)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
              end if

              !********************************* MAGNETIC FORCES **********************************************!
              if (fmagper) then

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof8,coeficiente8,modk,l)

                ! Computing the reciprocal space contribution
                kr     = (konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3))
                termo1 = (Di(q,i,1)*konda(1))+(Di(q,i,2)*konda(2))+(Di(q,i,3)*konda(3))
                termo2 = (Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3))

                if (Pe.eq.0.0) then
                  lambda=alpha2*8.0
                else
                  lambda=alpha2*8.0/Pe
                end if

                auxf(j,1)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(1)
                auxf(j,2)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(2)
                auxf(j,3)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(3)
              end if

            end if
          end if

        end do

        ! Making the final sum for the physical and reciprocal contribution for all possible interactions
        if (ligaih) then
          hidro2(s,1)=sum(hidrodinamica_aux2(:,1))
          hidro2(s,2)=sum(hidrodinamica_aux2(:,2))
          hidro2(s,3)=sum(hidrodinamica_aux2(:,3))
        end if

        if (fmagper) then
          forcareciproca(s,1)=sum(auxf(:,1))
          forcareciproca(s,2)=sum(auxf(:,2))
          forcareciproca(s,3)=sum(auxf(:,3))
        end if

        if (tmagper) then
          torquereciproco(s,1)=sum(auxt(:,1))
          torquereciproco(s,2)=sum(auxt(:,2))
          torquereciproco(s,3)=sum(auxt(:,3))
        end if

      end do

      if (ligaih) then
        ! Computing the velocities of the particles, now with the reciprocal space sum contribution
        U(q,i,1)= U(q,i,1) + sum(hidro2(:,1))
        U(q,i,2)= U(q,i,2) + sum(hidro2(:,2))
        U(q,i,3)= U(q,i,3) + sum(hidro2(:,3))
        contribuicao_reciproco(q,i)=sum(hidro2(:,3))
      end if

      if (fmagper) then
        ! Computing the magnetic forces acting on the particles, now with the reciprocal space sum contribution
        FORCAS(4,q,i,1)=FORCAS(4,q,i,1)+sum(forcareciproca(:,1))
        FORCAS(4,q,i,2)=FORCAS(4,q,i,2)+sum(forcareciproca(:,2))
        FORCAS(4,q,i,3)=FORCAS(4,q,i,3)+sum(forcareciproca(:,3))
      end if

      if (tmagper) then
        ! Computing the magnetic torques acting on the particles, now with the reciprocal space sum contribution
        TORQUES(1,q,i,1)=TORQUES(1,q,i,1)+sum(torquereciproco(:,1))
        TORQUES(1,q,i,2)=TORQUES(1,q,i,2)+sum(torquereciproco(:,2))
        TORQUES(1,q,i,3)=TORQUES(1,q,i,3)+sum(torquereciproco(:,3))
      end if

    end do
  end do
end subroutine periodic_interactions





subroutine periodic_interactions_old(k)

  implicit none
  integer, intent(in) :: k

  !------------------- parâmetros de corte (coerentes real/recíproco) -------------------
  ! Y controla ambos: rcut = Y/qsi e kcut = 2*qsi*Y
  real :: Y
  real :: rcut, rcut2, kcut
  !--------------------------------------------------------------------------------------

  ! variáveis locais
  real    :: mobilidade_self(3,3), mobilidade1(3,3), mobilidade2(3,3)
  real    :: coeficiente1, coeficiente2, coeficiente3, coeficiente4
  real    :: coeficiente5, coeficiente6, coeficiente7, coeficiente8
  real    :: rn(3)
  integer :: aux_periodico
  integer :: q, i, s, j, a, d
  integer, parameter :: dmin = 2, dmax = 9999
  real    :: eps
  real    :: Tx, Ty, Tz, T2, diag_caixa, r2, modrij
  real    :: kr, kr2, modk
  real    :: termo1, termo2, termo3, termo4, termo5
  real    :: konda(3), knormal(3)

  ! valores numéricos robustos
  eps = max(1.0e-12, tiny(1.0))

  ! ----------------- escolha do corte -----------------
  ! ajuste Y para varrer a convergência: 4–8 é um bom começo
  Y     = 5.0
  rcut  = Y / qsi
  rcut2 = rcut * rcut
  kcut  = 2.0 * qsi * Y
  ! ----------------------------------------------------

  diag_caixa = sqrt(l*l + l*l + h*h)

  ! We start by computing first the sums in the real space 
  ! Here, nb denotes the number of physical boxes. This
  ! number can be different from nbr, which denotes the
  ! nubmer of reciprocal boxes. In simconfig.dat we 
  ! recommend the user to set nb=125 and nbr=27. These
  ! values ensure a good precision when computing the 
  ! average sedimentation velocity of a suspension of
  ! spheres in Creeping-flow

  ! This loop works like this: we fix a realization and
  ! a given particle and the make a sweep computing the
  ! long-range interactions between this particles and
  ! all the other particles in the surrounding boxes.
  ! We then, extend this procedure to all real particles
  ! in all the simultaneous numerical experiments.

  do q=1,rea
    do i=1,N

      ! Building the self-mobility matrix (independe de j,s)
      if (ligaih) then
        mobilidade_self = 0.0
        do a=1,3
          mobilidade_self(a,a)= 1.0 - (6.0*(pi**(-0.5))*qsi) + ((40.0/3.0)*(pi**(-0.5))*(qsi**3.0))
        end do
      end if

      do s=1,nb
        do j=1,N

          aux_periodico = abs(ILF(s,1)) + abs(ILF(s,2)) + abs(ILF(s,3))

          ! Checking the distance between a particle "i" and other particle "j" in the real boxes (physical and images)  
          ! Calculating the "r" vector (distance)

          if (.not. shear) then
            Tx = ILF(s,1)*l
            Ty = ILF(s,2)*l
            Tz = ILF(s,3)*h
          else
            Tx = ILF(s,1)*l
            Ty = (ILF(s,2)*l) + (k*ILF(s,3)*shearrate*dt)
            Tz = ILF(s,3)*h
          end if

          ! culling rápido por caixa imagem (se muito longe, não alcança rcut na pior hipótese)
          T2 = Tx*Tx + Ty*Ty + Tz*Tz
          if (sqrt(T2) - diag_caixa > rcut) cycle

          ! vetor r_ij
          rij(1) = X(q,i,1) - (X(q,j,1) + Tx)
          rij(2) = X(q,i,2) - (X(q,j,2) + Ty)
          rij(3) = X(q,i,3) - (X(q,j,3) + Tz)

          ! Normalizing the "r" vector -> r/|r|
          r2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
          if (r2 <= 4.0) cycle         ! ignora esfera interna (r <= 2a) — igual ao seu >2.0
          if (r2 >  rcut2) cycle       ! corte no real

          modrij = sqrt(r2)

          if (modrij .ge. eps) then
            rn(1) = rij(1)/modrij
            rn(2) = rij(2)/modrij
            rn(3) = rij(3)/modrij
          else
            rn = 0.0
          end if

          d = 1 + int(9999.0*(modrij-2.0)/((3.0**0.5)*l*(nb**(1.0/3.0))))
          if (d < dmin) d = dmin
          if (d > dmax) d = dmax
          diferenca_interpol1 = (modrij - cof1(1,d))

          !************ SUMS IN THE PHYSICAL SPACE IN THE REAL BOX (WHERE THE REAL PARTICLES ARE) AND "IMAGE" BOXES***********!

          !******************************* HYDRODYNAMIC INTERACTIONS *********************************************************!

          ! Building the self-mobility matrix
          if (aux_periodico.eq.0 .and. i.eq.j) then
            ! já calculada fora
          end if

          ! Building the called mobility matrix 1 (lattice sum in the physical space)
          if (ligaih) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente1 = cof1(2,d) + ((modrij-cof1(1,d))/(cof1(1,d+1)-cof1(1,d))) * (cof1(2,d+1)-cof1(2,d))
                coeficiente2 = cof2(2,d) + ((modrij-cof2(1,d))/(cof2(1,d+1)-cof2(1,d))) * (cof2(2,d+1)-cof2(2,d))
              else
                coeficiente1 = cof1(2,d-1) + ((modrij-cof1(1,d-1))/(cof1(1,d)-cof1(1,d-1))) * (cof1(2,d)-cof1(2,d-1))
                coeficiente2 = cof2(2,d-1) + ((modrij-cof2(1,d-1))/(cof2(1,d)-cof2(1,d-1))) * (cof2(2,d)-cof2(2,d-1))
              end if

              mobilidade1(1,1)=coeficiente1 + coeficiente2*(rn(1)*rn(1))
              mobilidade1(1,2)=coeficiente2*(rn(1)*rn(2))
              mobilidade1(1,3)=coeficiente2*(rn(1)*rn(3))

              mobilidade1(2,1)=coeficiente2*(rn(2)*rn(1))
              mobilidade1(2,2)=coeficiente1 + coeficiente2*(rn(2)*rn(2))
              mobilidade1(2,3)=coeficiente2*(rn(2)*rn(3))

              mobilidade1(3,1)=coeficiente2*(rn(3)*rn(1))
              mobilidade1(3,2)=coeficiente2*(rn(3)*rn(2))
              mobilidade1(3,3)=coeficiente1 + coeficiente2*(rn(3)*rn(3))

              hidrodinamica_aux1(j,1)=mobilidade1(1,1)*FT(q,j,1)+mobilidade1(1,2)*FT(q,j,2)+mobilidade1(1,3)*FT(q,j,3)
              hidrodinamica_aux1(j,2)=mobilidade1(2,1)*FT(q,j,1)+mobilidade1(2,2)*FT(q,j,2)+mobilidade1(2,3)*FT(q,j,3)
              hidrodinamica_aux1(j,3)=mobilidade1(3,1)*FT(q,j,1)+mobilidade1(3,2)*FT(q,j,2)+mobilidade1(3,3)*FT(q,j,3)
            end if
          end if

          !****************** MAGNETIC TORQUES **************************************!
          if (tmagper) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente4=cof4(2,d) + ((modrij-cof4(1,d))/(cof4(1,d+1)-cof4(1,d)))*(cof4(2,d+1)-cof4(2,d))
                coeficiente5=cof5(2,d) + ((modrij-cof5(1,d))/(cof5(1,d+1)-cof5(1,d)))*(cof5(2,d+1)-cof5(2,d))
              else
                coeficiente4=cof4(2,d-1) + ((modrij-cof4(1,d-1))/(cof4(1,d)-cof4(1,d-1)))*(cof4(2,d)-cof4(2,d-1))
                coeficiente5=cof5(2,d-1) + ((modrij-cof5(1,d-1))/(cof5(1,d)-cof5(1,d-1)))*(cof5(2,d)-cof5(2,d-1))
              end if

              ! Computing periodic torques due to magnetic interactions in the real space
              termo5 = 0.0
              termo2 = (Di(q,j,1)*rij(1)) + (Di(q,j,2)*rij(2)) + (Di(q,j,3)*rij(3))

              if (gravidade) then
                lambda=alpha2*6.0*pi/Pe
              else
                lambda=alpha2*6.0*pi
              end if

              termo1 = (Di(q,i,3)*Di(q,j,2))-(Di(q,i,2)*Di(q,j,3))
              termo3 = (Di(q,i,2)*rij(3))-(Di(q,i,3)*rij(2))
              auxt(j,1) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)

              termo1 = (Di(q,i,1)*Di(q,j,3))-(Di(q,i,3)*Di(q,j,1))
              termo3 = (Di(q,i,3)*rij(1))-(Di(q,i,1)*rij(3))
              auxt(j,2) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)

              termo1 = (Di(q,i,2)*Di(q,j,1))-(Di(q,i,1)*Di(q,j,2))
              termo3 = (Di(q,i,1)*rij(2))-(Di(q,i,2)*rij(1))
              auxt(j,3) = lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)
            end if
          end if

          !*********************************** MAGNETIC FORCES *****************************************!
          if (fmagper) then
            if (modrij.gt.2.0) then

              ! Interpolating the pre-calculated Green functions
              if (diferenca_interpol1.gt.0.0) then
                coeficiente6=cof6(2,d) + ((modrij-cof6(1,d))/(cof6(1,d+1)-cof6(1,d)))*(cof6(2,d+1)-cof6(2,d))
              else
                coeficiente6=cof6(2,d-1) + ((modrij-cof6(1,d-1))/(cof6(1,d)-cof6(1,d-1)))*(cof6(2,d)-cof6(2,d-1))
              end if

              if (gravidade) then
                lambda=alpha2*8.0*pi/Pe
              else
                lambda=alpha2*8.0*pi
              end if

              ! (di.dj)rij
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(1)
              ! (dj.rij)di
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,1)
              ! (di.rij)dj
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,1)
              ! ((di.rij)(dj.rij))rij
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(1)

              auxf(j,1)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

              ! y
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(2)
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,2)
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,2)
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(2)

              auxf(j,2)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

              ! z
              termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(3)
              termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,3)
              termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,3)
              termo4=(((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3))) * ((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(3)

              auxf(j,3)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)
            end if
          end if

        end do

        !************************** SUMING EVERYTHING *********************************!

        ! Computing the sum of the previous vectors and matrices
        if (ligaih) then
          hidro1(s,1)=sum(hidrodinamica_aux1(:,1))
          hidro1(s,2)=sum(hidrodinamica_aux1(:,2))
          hidro1(s,3)=sum(hidrodinamica_aux1(:,3))
        end if

        if (fmagper) then
          forcareal(s,1)=sum(auxf(:,1))
          forcareal(s,2)=sum(auxf(:,2))
          forcareal(s,3)=sum(auxf(:,3))
        end if

        if (tmagper) then
          torquereal(s,1)=sum(auxt(:,1))
          torquereal(s,2)=sum(auxt(:,2))
          torquereal(s,3)=sum(auxt(:,3))
          auxt = 0.0
        end if

      end do

      ! Computing the contributions of the real space interactions in the velocities of each particle
      if (ligaih) then
        U(q,i,1)=mobilidade_self(1,1)*FT(q,i,1) + sum(hidro1(:,1))
        U(q,i,2)=mobilidade_self(2,2)*FT(q,i,2) + sum(hidro1(:,2))
        U(q,i,3)=mobilidade_self(3,3)*FT(q,i,3) + sum(hidro1(:,3))
      end if

      hidro1=0.0

      ! Computing the contributions of the real space interactions in the magnetic forces and torques acting on each particle
      if (fmagper) then
        FORCAS(4,q,i,1)=sum(forcareal(:,1))
        FORCAS(4,q,i,2)=sum(forcareal(:,2))
        FORCAS(4,q,i,3)=sum(forcareal(:,3))
      end if

      if (tmagper) then
        TORQUES(1,q,i,1)=sum(torquereal(:,1))
        TORQUES(1,q,i,2)=sum(torquereal(:,2))
        TORQUES(1,q,i,3)=sum(torquereal(:,3))
      end if

      contribuicao_self(q,i)   = mobilidade_self(3,3)*FT(q,i,3)
      contribuicao_fisico(q,i) = sum(hidro1(:,3))

    end do
  end do

  !************************************* SUMS IN THE RECIPROCAL SPACE ***********************************!
  do q=1,rea
    do i=1,N
      do s=1,nbr
        do j=1,N

          aux_periodico = abs(ILR(s,1)) + abs(ILR(s,2)) + abs(ILR(s,3))

          if (aux_periodico.ne.0) then

            rij(1)= X(q,i,1)-X(q,j,1)
            rij(2)= X(q,i,2)-X(q,j,2)
            rij(3)= X(q,i,3)-X(q,j,3)

            ! Calculating the normalized "r" vector = r/|r|
            modrij = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))

            if (modrij .ne. 0.0) then
              rn(1)=rij(1)/modrij
              rn(2)=rij(2)/modrij
              rn(3)=rij(3)/modrij
            else
              rn = 0.0
            end if

            if (modrij.ne.0.0) then

              ! Computing the wave number vector based on the lattice index 
              konda(1)=ILR(s,1)*2.0*pi/l
              konda(2)=ILR(s,2)*2.0*pi/l
              konda(3)=ILR(s,3)*2.0*pi/h

              ! Calculating the wave number module
              modk = sqrt(konda(1)*konda(1) + konda(2)*konda(2) + konda(3)*konda(3))

              ! -------------------- CORTE NO RECÍPROCO --------------------
              if (modk > kcut) cycle
              ! ------------------------------------------------------------

              ! Normalized wave number vector
              knormal(1)=konda(1)/modk
              knormal(2)=konda(2)/modk
              knormal(3)=konda(3)/modk

              !********************** HYDRODYNAMIC INTERACTIONS *********************************************************!
              if (ligaih) then

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof3,coeficiente3,modk,l)

                ! Calculating the mobility matrix for the reciprocal space contribution
                mobilidade2(1,1)=coeficiente3*(1.0 - knormal(1)*knormal(1))
                mobilidade2(1,2)=coeficiente3*(-knormal(1)*knormal(2))
                mobilidade2(1,3)=coeficiente3*(-knormal(1)*knormal(3))
                mobilidade2(2,1)=coeficiente3*(-knormal(2)*knormal(1))
                mobilidade2(2,2)=coeficiente3*(1.0 - knormal(2)*knormal(2))
                mobilidade2(2,3)=coeficiente3*(- knormal(2)*knormal(3))
                mobilidade2(3,1)=coeficiente3*(-knormal(3)*knormal(1))
                mobilidade2(3,2)=coeficiente3*(-knormal(3)*knormal(2))
                mobilidade2(3,3)=coeficiente3*(1.0 - knormal(3)*knormal(3))

                kr  = (knormal(1)*rn(1)) + (knormal(2)*rn(2)) + (knormal(3)*rn(3))
                kr2 = cos((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

                hidrodinamica_aux2(j,1)=kr*(mobilidade2(1,1)*FT(q,j,1)+mobilidade2(1,2)*FT(q,j,2)+mobilidade2(1,3)*FT(q,j,3))
                hidrodinamica_aux2(j,2)=kr*(mobilidade2(2,1)*FT(q,j,1)+mobilidade2(2,2)*FT(q,j,2)+mobilidade2(2,3)*FT(q,j,3))
                hidrodinamica_aux2(j,3)=kr*(mobilidade2(3,1)*FT(q,j,1)+mobilidade2(3,2)*FT(q,j,2)+mobilidade2(3,3)*FT(q,j,3))
              end if

              !**************************** MAGNETIC TORQUES ***********************************************!
              if (tmagper) then

                kr2 = 1.0*((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof7,coeficiente7,modk,l)

                ! Computing the reciprocal space contribution
                termo2 = (Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3))

                if (gravidade) then
                  lambda=alpha2*6.0*pi/Pe
                else
                  lambda=alpha2*6.0*pi
                end if

                termo1=(Di(q,i,2)*konda(3))-(Di(q,i,3)*konda(2))
                auxt(j,1)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
                termo1=(Di(q,i,3)*konda(1))-(Di(q,i,1)*konda(3))
                auxt(j,2)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
                termo1=(Di(q,i,1)*konda(2))-(Di(q,i,2)*konda(1))
                auxt(j,3)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))
              end if

              !********************************* MAGNETIC FORCES **********************************************!
              if (fmagper) then

                ! Interpolating the pre-calculated Green function
                call interpola_reciproco(nbr,cof3,cof8,coeficiente8,modk,l)

                ! Computing the reciprocal space contribution
                kr     = (konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3))
                termo1 = (Di(q,i,1)*konda(1))+(Di(q,i,2)*konda(2))+(Di(q,i,3)*konda(3))
                termo2 = (Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3))

                if (Pe.eq.0.0) then
                  lambda=alpha2*8.0
                else
                  lambda=alpha2*8.0/Pe
                end if

                auxf(j,1)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(1)
                auxf(j,2)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(2)
                auxf(j,3)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(3)
              end if

            end if
          end if

        end do

        ! Making the final sum for the physical and reciprocal contribution for all possible interactions
        if (ligaih) then
          hidro2(s,1)=sum(hidrodinamica_aux2(:,1))
          hidro2(s,2)=sum(hidrodinamica_aux2(:,2))
          hidro2(s,3)=sum(hidrodinamica_aux2(:,3))
        end if

        if (fmagper) then
          forcareciproca(s,1)=sum(auxf(:,1))
          forcareciproca(s,2)=sum(auxf(:,2))
          forcareciproca(s,3)=sum(auxf(:,3))
        end if

        if (tmagper) then
          torquereciproco(s,1)=sum(auxt(:,1))
          torquereciproco(s,2)=sum(auxt(:,2))
          torquereciproco(s,3)=sum(auxt(:,3))
        end if

      end do

      if (ligaih) then
        ! Computing the velocities of the particles, now with the reciprocal space sum contribution
        U(q,i,1)= U(q,i,1) + sum(hidro2(:,1))
        U(q,i,2)= U(q,i,2) + sum(hidro2(:,2))
        U(q,i,3)= U(q,i,3) + sum(hidro2(:,3))
        contribuicao_reciproco(q,i)=sum(hidro2(:,3))
      end if

      if (fmagper) then
        ! Computing the magnetic forces acting on the particles, now with the reciprocal space sum contribution
        FORCAS(4,q,i,1)=FORCAS(4,q,i,1)+sum(forcareciproca(:,1))
        FORCAS(4,q,i,2)=FORCAS(4,q,i,2)+sum(forcareciproca(:,2))
        FORCAS(4,q,i,3)=FORCAS(4,q,i,3)+sum(forcareciproca(:,3))
      end if

      if (tmagper) then
        ! Computing the magnetic torques acting on the particles, now with the reciprocal space sum contribution
        TORQUES(1,q,i,1)=TORQUES(1,q,i,1)+sum(torquereciproco(:,1))
        TORQUES(1,q,i,2)=TORQUES(1,q,i,2)+sum(torquereciproco(:,2))
        TORQUES(1,q,i,3)=TORQUES(1,q,i,3)+sum(torquereciproco(:,3))
      end if

    end do
  end do

end subroutine periodic_interactions_old

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: interpola_reciproco		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for interpolating the va- !
! lues of the pre-calculated the Green functions  !
! used in the calculation of periodic interactions!
! due to hydrodynamic forces, magnetic forces and !
! torques (only in the reciprocal space)          !
!*************************************************!

subroutine interpola_reciproco(a,b,c,d,e,f)
  implicit none
  ! a = number of reciprocal lattices (e.g., 27 or 125)
  ! b = cof3   (grade de k de referência)   [2 x 10000]
  ! c = cofx   (tabela a interpolar)        [2 x 10000]
  ! d = coeficientex (saída: valor interpolado)
  ! e = |k| (modk)
  ! f = l (comprimento da caixa)

  integer, intent(in)  :: a
  real   , intent(in)  :: b(2,10000), c(2,10000), e, f
  real   , intent(out) :: d

  integer :: nside, g
  real    :: kmin, kmax, t, diferenca_interpol2
  real    :: two_pi_over_l, root3

  two_pi_over_l = 2.0*pi / f
  root3         = sqrt(3.0)

  ! nside = a^(1/3) arredondado: esperamos 3 (27) ou 5 (125)
  nside = nint(a**(1.0/3.0))

  ! Definição do intervalo tabulado de |k| para a malha recíproca usada
  select case (nside)
  case (5)      ! 125 = 5^3
    kmin = two_pi_over_l                                ! 1ª casca: 1*2π/L
    kmax = (4.0*two_pi_over_l)*root3                    ! casca máxima ~ 4√3 * 2π/L
  case (3)      ! 27 = 3^3
    kmin = two_pi_over_l
    kmax = (2.0*two_pi_over_l)*root3
  case default
    ! fallback conservador (evita divisão por zero e mantém funcional)
    kmin = two_pi_over_l
    kmax = (real(nside-1)*two_pi_over_l)*root3
    if (kmax <= kmin) kmax = kmin + two_pi_over_l
  end select

  ! Parametrização no intervalo [0,1] com clamp para robustez
  if (e <= kmin) then
    t = 0.0
  else if (e >= kmax) then
    t = 1.0
  else
    t = (e - kmin) / (kmax - kmin)
  end if

  ! Escolha de índice garantindo que usaremos (g-1,g) ou (g,g+1) válidos:
  ! Aqui geramos g em [2,9999] para sempre termos g-1 >= 1 e g+1 <= 10000.
  g = 2 + int( t * 9997.0 )
  if (g < 2)     g = 2
  if (g > 9999)  g = 9999

  diferenca_interpol2 = e - b(1,g)

  if (diferenca_interpol2 >= 0.0) then
    ! interpola entre g e g+1
    d = c(2,g) + ((e - c(1,g)) / (c(1,g+1) - c(1,g))) * (c(2,g+1) - c(2,g))
  else
    ! interpola entre g-1 e g
    d = c(2,g-1) + ((e - c(1,g-1)) / (c(1,g) - c(1,g-1))) * (c(2,g) - c(2,g-1))
  end if

end subroutine interpola_reciproco


subroutine interpola_reciproco_old(a,b,c,d,e,f)
  implicit none

! a = number of reciprocal lattices
! b = cof3
! c = cofx
! d = coeficientex
! e = norm of k  (modk)
! f = l (length of the box)

  integer :: a
  real    :: b(2,10000), c(2,10000), d
  integer :: g
  real    :: modk, l, diferenca_interpol2, e, f

  if(a**(1.0/3.0).eq.5.0)then
    g = 1 + int( 9999.0*((e-(2*pi/f))/(((4*pi/f)*(3.0**0.5))-(2*pi/f))) )
  end if
  if(nbr**(1.0/3.0).eq.3.0)then
    g = 1 + int( 9999.0*((e-(2*pi/f))/(((2*pi/f)*(3.0**0.5))-(2*pi/f))) )
  end if
  diferenca_interpol2 = (e - b(1,g))

  if(g-1.eq.0) then
    g = g + 1
  end if

  if(diferenca_interpol2.gt.0.0)then
    d = c(2,g) + ((e-c(1,g))/(c(1,g+1)-c(1,g)))*(c(2,g+1)-c(2,g))
  else
    d = c(2,g-1) + ((e-c(1,g-1))/(c(1,g)-c(1,g-1)))*(c(2,g)-c(2,g-1))
  end if

end subroutine interpola_reciproco_old


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resvel				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the parti-!
! cles velocity using a 4th order Runge-Kutta.    !
!*************************************************!

! IT IS IMPORTANT TO NOTICE THAT THIS SUBROUTINE 
! DOES NOT APPLY FOR A MOBILITY PROBLEM. IT ONLY 
! MAKES SENSE WHEN WE CONSIDER PARTICLE INERTIA 
! (RESISTANCE FORMULATION)

subroutine resvel(a,b,c,d)
  implicit none
  real :: a                      ! velocity component
  real :: b                      ! time step
  real :: c                      ! Stokes number (null for zero inertia)
  real :: d                      ! sum of forces in a given direction
  real :: k1, k2, k3, k4         ! internal variables

  k1 = b*(-a + d)/c
  k2 = b*((-a - 0.5*k1) + d)/c
  k3 = b*((-a - 0.5*k2) + d)/c
  k4 = b*((-a - k3) + d)/c

  a = a + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

end subroutine resvel

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: respos				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the parti-!
! cles position.				  !
!*************************************************!

subroutine respos(a,b,c)
  implicit none
  real :: a                      ! position
  real :: b                      ! time-step
  real :: c                      ! velocity component

  a = a + b*c
end subroutine respos

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: writting_files			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for writting output files !
! containing the positions, velocities and orien- !
! tations of all particles in all realiations     !
!*************************************************!

subroutine writting_files(k,k_real)
  implicit none
  integer :: k, teste2
  real    :: k_real, teste1
  integer :: i, j

509 FORMAT(F30.4,F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)
if(continua)then
teste1=k/n3
teste2=k/n2

if(teste1.eq.teste2) then

if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt
write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

else

if(k.eq.1) then
if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt
write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

teste1=k/n3
teste2=k/n2

if(k.ne.1) then

if(teste1.eq.teste2) then
if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,*) 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt

write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'  

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if
end if
end if
end if

end subroutine writting_files

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: torque_magnetico			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to dipolar interactions when this   !
! computation is performed in a non-periodic way  !
!*************************************************!
subroutine torque_magnetico
  implicit none
  integer :: auxiliary
  integer :: i, j, q
  real    :: r, modrij
  real    :: termo1, termo2, termo3

if(gravidade)then
lambda=alpha2*24.0
else
if(shear) then
lambda=alpha2*3.0/Pe
else
lambda=alpha2*9.0/2.0
end if
end if

! These subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

if(mistura)then
 auxiliary = int(percentual*N) + 1
else
 auxiliary=1
end if

do j=1,rea
do i=auxiliary,N
do q=auxiliary,N
if(i.ne.q) then
! Calculating the distance between all the pairs of particles
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5
if(r.le.2.0) then
!if(r.le.2.0) then
aux1(j,q)=0.0
aux2(j,q)=0.0
aux3(j,q)=0.0
else

! Calculating the vector that connects a particle i to a particle j
rij(1)=X(j,i,1)-X(j,q,1)
rij(2)=X(j,i,2)-X(j,q,2)
rij(3)=X(j,i,3)-X(j,q,3)
! Normalizing this vector
modrij=((rij(1)**2.0)+(rij(2)**2.0)+(rij(3)**2.0))**0.5
rij(1)=rij(1)/modrij
rij(2)=rij(2)/modrij
rij(3)=rij(3)/modrij

termo2=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))

termo1=((Di(j,i,2)*Di(j,q,3))-(Di(j,i,3)*Di(j,q,2)))
termo3=((Di(j,i,2)*rij(3))-(Di(j,i,3)*rij(2)))

aux1(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

termo1=((Di(j,i,3)*Di(j,q,1))-(Di(j,i,1)*Di(j,q,3)))
termo3=((Di(j,i,3)*rij(1))-(Di(j,i,1)*rij(3)))

aux2(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

termo1=((Di(j,i,1)*Di(j,q,2))-(Di(j,i,2)*Di(j,q,1)))
termo3=((Di(j,i,1)*rij(2))-(Di(j,i,2)*rij(1)))

aux3(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

end if
end if
end do
TORQUES(1,j,i,1)=sum(aux1(j,:))
TORQUES(1,j,i,2)=sum(aux2(j,:))
TORQUES(1,j,i,3)=sum(aux3(j,:))

aux1=0.0
aux2=0.0
aux3=0.0
r=0.0
modrij=0.0

end do
end do

end subroutine torque_magnetico


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: torque_externo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to an external field		  !
!*************************************************!

subroutine torque_externo(alfa)
  implicit none
  integer :: auxiliary
  real    :: alfa
  integer :: i, j

! This subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

  if(mistura)then
    auxiliary = int(percentual*N) + 1
  else
    auxiliary = 1
  end if

 !$OMP PARALLEL DO
  do j=1,rea
    do i=auxiliary,N

! posicao_campo = 1 -> Applied field on the lower wall

      if(posicao_campo.eq.1)then
        if(gravidade) then
          TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/(4.0*Pe))
          TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,1)/(4.0*Pe))
          TORQUES(2,j,i,3)=0.0
        else
          if(shear) then
            TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/(4.0*Pe))
            TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,1)/(4.0*Pe))
            TORQUES(2,j,i,3)=0.0
          else
            TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/4.0)
            TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,1)/4.0)
            TORQUES(2,j,i,3)=0.0
          end if
        end if
      end if

! posicao_campo = 2 ->  Applied field on the upper wall

      if(posicao_campo.eq.2)then
        if(gravidade) then
          TORQUES(2,j,i,1)=( 3.0*alfa*Di(j,i,2)/(4.0*Pe))
          TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/(4.0*Pe))
          TORQUES(2,j,i,3)=0.0
        else
          if(shear) then
            TORQUES(2,j,i,1)=( 3.0*alfa*Di(j,i,2)/(4.0*Pe))
            TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/(4.0*Pe))
            TORQUES(2,j,i,3)=0.0
          else
            TORQUES(2,j,i,1)=( 3.0*alfa*Di(j,i,2)/4.0)
            TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/4.0)
            TORQUES(2,j,i,3)=0.0
          end if
        end if
      end if

! posicao_campo = 3 ->  Applied field on the right side

      if(posicao_campo.eq.3) then
        if(gravidade) then
          TORQUES(2,j,i,1)=0.0
          TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/(4.0*Pe))
          TORQUES(2,j,i,3)=( 3.0*alfa*Di(j,i,2)/(4.0*Pe))
        else
          if(shear) then
            TORQUES(2,j,i,1)=0.0
            TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/(4.0*Pe))
            TORQUES(2,j,i,3)=( 3.0*alfa*Di(j,i,2)/(4.0*Pe))
          else
            TORQUES(2,j,i,1)=0.0
            TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/4.0)
            TORQUES(2,j,i,3)=( 3.0*alfa*Di(j,i,2)/4.0)
          end if
        end if
      end if

! posicao_campo = 4 ->  Applied field on the left side

      if(posicao_campo.eq.4) then
        if(gravidade) then
          TORQUES(2,j,i,1)=0.0
          TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,3)/(4.0*Pe))
          TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/(4.0*Pe))
        else
          if(shear) then
            TORQUES(2,j,i,1)=0.0
            TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,3)/(4.0*Pe))
            TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/(4.0*Pe))
          else
            TORQUES(2,j,i,1)=0.0
            TORQUES(2,j,i,2)=( 3.0*alfa*Di(j,i,3)/4.0)
            TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/4.0)
          end if
        end if
      end if

    end do
  end do
 !$OMP END PARALLEL DO
end subroutine torque_externo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: rotating_field			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to an external rotating field	  !
!*************************************************!

subroutine rotating_field(alfa, wt)
  implicit none
  integer :: auxiliary
  real    :: alfa, wt
  integer :: i, j

! This subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

  if(mistura)then
    auxiliary = int(percentual*N) + 1
  else
    auxiliary = 1
  end if

 !$OMP PARALLEL DO
  do j=1,rea
    do i=auxiliary,N

! Rotating field on the upper side

      if(gravidade) then
        TORQUES(2,j,i,1)=(3.0*alfa/(4.0*Pe))*( Di(j,i,2)*cos(wt) - Di(j,i,3)*sin(wt) )
        TORQUES(2,j,i,2)=(3.0*alfa/(4.0*Pe))*( - Di(j,i,1)*cos(wt) )
        TORQUES(2,j,i,3)=(3.0*alfa/(4.0*Pe))*(   Di(j,i,1)*sin(wt) )
      else
        if(shear) then
          TORQUES(2,j,i,1)=(3.0*alfa/(4.0*Pe))*( Di(j,i,2)*cos(wt) - Di(j,i,3)*sin(wt) )
          TORQUES(2,j,i,2)=(3.0*alfa/(4.0*Pe))*( - Di(j,i,1)*cos(wt) )
          TORQUES(2,j,i,3)=(3.0*alfa/(4.0*Pe))*(   Di(j,i,1)*sin(wt) )
        else
          TORQUES(2,j,i,1)=(3.0*alfa/4.0)*( Di(j,i,2)*cos(wt) - Di(j,i,3)*sin(wt) )
          TORQUES(2,j,i,2)=(3.0*alfa/4.0)*( - Di(j,i,1)*cos(wt) )
          TORQUES(2,j,i,3)=(3.0*alfa/4.0)*(   Di(j,i,1)*sin(wt) )
        end if
      end if
    end do
  end do
 !$OMP END PARALLEL DO

end subroutine rotating_field


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resvel_sem_inercia			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the velo- !
! city of the particles in a non inertial context !
!*************************************************!

subroutine resvel_sem_inercia(a,b,c)
  implicit none
  real :: a   ! velocity component
  real :: b   ! sum of forces in a given direction
  real :: c   ! 0 or 1, 0 = no gravity, 1 = gravity

  a = -c + b
end subroutine resvel_sem_inercia

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resomega				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the angu- !
! lar velocity of the particles 		  !
!*************************************************!

subroutine resomega(a,b,c,d)
  implicit none
  real :: a  ! rotational velocitty in a given direction
  real :: b  ! time step
  real :: c  ! Rotational Stokes number
  real :: d  ! Sum of torques in a given direction
  real :: k1, k2, k3, k4  ! Internal variables

  k1 = b*(-a + d)/c
  k2 = b*((-a - 0.5*k1) + d)/c
  k3 = b*((-a - 0.5*k2) + d)/c
  k4 = b*((-a - k3) + d)/c

  a  = a + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
end subroutine resomega


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: evoldip				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for evolving the dipole	  !
! orientation for each particle			  !
!*************************************************!

subroutine evoldip(a,b,c,d,e,f)
  implicit none
  real :: a  ! Dipole moment orientation to be calculated
  real :: b  ! Particle dipole moment in the j direction
  real :: c  ! Particle dipole moment in the k direction
  real :: d  ! Angular velocity in the j direction
  real :: e  ! Angular velocity in the k direction
  real :: f  ! Time-step

! d(di)/dt = omega_i x di for particle i

  a = a + (d*c - e*b)*f
end subroutine evoldip


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: randomica				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating all random!
! numbers used in all simulations		  !
!*************************************************!

subroutine randomica(a,b,c,n,d)
implicit none
  ! a,b = random number range
  real,    intent(in)  :: a, b
  integer, intent(in)  :: n, d
  ! c = generated random sequence (size = n)
  real,    intent(out) :: c(*)

  integer :: m, i, e
  integer :: f(8)
  integer, allocatable :: seedrand(:)
  
  call random_seed(size = m)
  if (m > 0) then
    allocate(seedrand(m))
  end if

  call date_and_time(values=f)
  call system_clock(count=e)

  if (m > 0) then
    do i = 1, m
      seedrand(i) = 47*d + f(8)*i*d*12000 + e*3*d 
    end do
    call random_seed(put = seedrand)
    deallocate(seedrand)
  end if

  call random_number(c(1:n))
  c(1:n) = a + (b - a) * c(1:n)

end subroutine randomica


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: ordem_crescente			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for allocating numbers in !
! ascending order				  !
!*************************************************!

subroutine ordem_crescente(k, s, n)
  implicit none
  integer, intent(in)    :: n
  real,    intent(inout) :: k(n), s(n)
  integer :: j, q
  real    :: temp1, temp2, temp3, temp4, aux

  do j = 1, n-1
    do q = j+1, n
      aux = k(q) - k(j)
      if (aux .lt. 0.0) then
        temp1 = k(q);  temp2 = k(j)
        temp3 = s(q);  temp4 = s(j)
        k(j) = temp1;  k(q) = temp2
        s(j) = temp3;  s(q) = temp4
      end if
    end do
  end do
end subroutine ordem_crescente


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: media_dipolo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating an ensem- !
! ble average of the particles orientation to cal-!
! culate the suspension's magnetization		  !
!*************************************************!

subroutine media_dipolo(U,N,rea,media,direcao)
  implicit none
  integer :: N, rea
  integer :: i, j, k
  integer :: direcao
  real    :: U(rea,N,3)
  real    :: Di(rea,N,3)
  real    :: media

  Di = U

  media = sum(Di(:,:,direcao)) / (N*rea)
end subroutine media_dipolo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: acopla_particulas			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for sticking particles to !
! each other as soon as they form an aggregate.   !
! This subroutine prevents aggregate breakups.    !
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine acopla_particulas(X,N,D,rea,Di)
  implicit none
  integer :: N, D, rea
  integer :: i, j, k, s, iter, iter2
  integer :: LOOP, AGRE, PARTAGRE
  integer :: AGREGPART(N/2)
  integer :: PARTICULATESTE
  integer :: PARTPORAGRE
  integer :: agregtotal(N)
  integer :: parteste(N)
  integer :: agreteste((N/2),N)
  real    :: X(rea,N,3), dt
  integer :: G(N)
  integer :: A(N/2,N)
  real    :: aux1, aux2, aux3, aux4, aux5
  real    :: Di(rea,N,3)

  G = 0

  if (D .eq. 1) then
    agregtotal = 0
  end if

  LOOP = 0
  AGRE = 0
  PARTAGRE = 0
  PARTICULATESTE = 0
  AGREGPART = 0
  parteste = 0
  agreteste = 0

  i = 1
  j = 2

10 continue
  if (parteste(i) .ne. 0) then
    LOOP = 1
  end if

  if (parteste(j) .ne. 0) then
    LOOP = 1
  end if

  if (i .ne. j) then
    r = (((X(D,i,1)-X(D,j,1))**2.0) + ((X(D,i,2)-X(D,j,2))**2.0) + ((X(D,i,3)-X(D,j,3))**2.0))**0.5

    if (r .le. 2.06) then
      if (LOOP .eq. 0) then
        AGRE      = AGRE + 1
        PARTAGRE  = PARTAGRE + 1
        parteste(i) = i
        parteste(j) = j

        agreteste(AGRE,i) = i
        agreteste(AGRE,j) = j

        AGREGPART(AGRE)   = PARTAGRE
      else
        do k = 1, N/2
          if (agreteste(k,i) .eq. i) then
            AGRE = k
          end if
        end do
        parteste(i) = i
        parteste(j) = j
        PARTAGRE = AGREGPART(AGRE) + 1
        AGREGPART(AGRE) = PARTAGRE
        agreteste(AGRE,j) = j
      end if

    end if
  end if

  j = j + 1

  if (j .eq. N) then
    i = i + 1
    j = i + 1
    LOOP = 0
    PARTAGRE = 1
  end if

  if (i .ne. N-1) then
    LOOP = 0
    PARTAGRE = 1
    go to 10
  end if

  do k = 1, AGRE
    do i = 2, N
      A(k,i) = agreteste(k,i)
    end do
  end do

! Coupling a variable of particles that belong to the same aggregate
! This quantity can be force, torque or dipole orientation of the 
! particles. This guarantee that the particles inside the same 
! aggregate structure will move and rotate together.

  do k = 1, AGRE

    iter = 1

    do i = 2, N
      if (A(k,i) .ne. 0) then
        if (iter .eq. 1) then
          aux3 = Di(D,i,1)
          aux4 = Di(D,i,2)
          aux5 = Di(D,i,3)
          iter = iter + 1
        else
          Di(D,i,1) = aux3
          Di(D,i,2) = aux4
          Di(D,i,3) = aux5
        end if
      end if

    end do
  end do

end subroutine acopla_particulas


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: analise_de_agregados		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for analyzing the number  !
! of aggregates and couting the particle distri-  !
! bution inside aggregates. It works together with!
! subroutine fdp_agregados			  ! 
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine analise_de_agregados(X,N,D,rea,s)
  implicit none
  integer :: N, D, rea, s
  integer :: i, j, k
  integer :: LOOP, AGRE, PARTAGRE
  integer :: AGREGPART(N/2)
  integer :: PARTICULATESTE
  integer :: PARTPORAGRE
  integer :: agregtotal(N)
  integer :: parteste(N)
  integer :: agreteste((N/2),N)
  real    :: X(rea,N,3), dt
  integer :: G(N)
  real    :: aux1, aux2
  real    :: r

  G = 0

  if (D.eq.1) then
    agregtotal = 0
  end if

  LOOP = 0
  AGRE = 0
  PARTAGRE = 0
  PARTICULATESTE = 0
  AGREGPART = 0
  parteste = 0
  agreteste = 0

  i = 1
  j = 2

10 if (parteste(i).ne.0) then
     LOOP = 1
   end if

  if (parteste(j).ne.0) then
    LOOP = 1
  end if

  if (i.ne.j) then
    r = (((X(D,i,1)-X(D,j,1))**2.0) + ((X(D,i,2)-X(D,j,2))**2.0) + ((X(D,i,3)-X(D,j,3))**2.0))**0.5
    if (r.le.2.06) then
      if (LOOP.eq.0.0) then
        AGRE = AGRE + 1
        PARTAGRE = PARTAGRE + 1
        parteste(i) = i
        parteste(j) = j

        agreteste(AGRE,i) = i
        agreteste(AGRE,j) = j

        AGREGPART(AGRE) = PARTAGRE
      else
        do k = 1, N/2
          if (agreteste(k,i).eq.i) then
            AGRE = k
          end if
        end do
        parteste(i) = i
        parteste(j) = j
        PARTAGRE = AGREGPART(AGRE) + 1
        AGREGPART(AGRE) = PARTAGRE
        agreteste(AGRE,j) = j
      end if
    end if
  end if

  j = j + 1

  if (j.eq.N) then
    i = i + 1
    j = i + 1
    LOOP = 0
    PARTAGRE = 1
  end if

  if (i.ne.N-1) then
    LOOP = 0
    PARTAGRE = 1
    go to 10
  end if

! Counting the number of aggregates with N particles

  do i = 2, N
    do k = 1, AGRE
      if (AGREGPART(k).eq.i) then
        G(i) = G(i) + 1
      end if
    end do
  end do

  do i = 2, N
    agregtotal(i) = agregtotal(i) + G(i)
  end do

  if (D.eq.rea) then
    call fdp_agregados(agregtotal, N, s)
  end if
end subroutine analise_de_agregados

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: fdp_agregados			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Acessory subroutine used in subroutine:	  !
! analise_de_agregados				  ! 
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine fdp_agregados(A,N,k)
  implicit none
  integer :: N, i, k
  integer :: A(N)
  integer :: numero_agregados
  real    :: aux2, aux1
  real    :: f(N)

  numero_agregados = sum(A)
  aux1 = numero_agregados

  do i = 2, N
    f(i) = A(i) ! / aux1
  end do

  write(200*rea,*) f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10), k*dt
end subroutine fdp_agregados


!************************************************
end module subroutines
