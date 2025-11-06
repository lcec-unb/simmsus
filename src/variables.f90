module variables
implicit none
! Declaring all global variables used in the code


! INTEGERS (NON ALLOCATABLE)
integer N, razao,iter,posicao_campo
integer nb,nbr,npast,rea, nnr
integer inteiro, inteiro2
integer n2,auxiliar_continua,DeAllocateStatus
integer nfreq, contfreqinteiro1, contfreqinteiro2, multiplofreq

! REAL (NON ALLOCATABLE)
real, parameter :: pi = acos(-1.0)
real razao2, diferenca_interpol1, diferenca_interpol2
real lambda, kreal, aux_real, bifmax
real derivada1, derivada2, derivada3, intervalo
real modulodipolo,kr2,dt,percentual,qsi, modk, kr
real rij(3), freqmax,C1, C2, C3, C4, k_real
real l,h,phi,nr1, nr2, nr3, Pe, Per, ragreg
real r, alpha2,alpha,xcentro, ycentro, zcentro
real xmin, xmax, ymin, ymax, zmin, zmax
real modrij, modrand,auxcont, freq, shearratei
real termo1,termo2,termo3,termo4,termo5,eps
real tempo,dist,difusao
real reale, reale2, Str, freqcampo, freqbeat
real n3, shearrate, St
real UMEDIA(3),SIGMA(3), konda(3), knormal(3)
real k1,k2,k3,k4,g1,g2,g3,g4
real ILF(125,3), ILR(27,3)
real hidro1(125,3), hidro2(27,3)
real cof1(2,10000), cof2(2,10000)
real cof3(2,10000), cof4(2,10000)
real cof5(2,10000), cof6(2,10000)
real cof7(2,10000), cof8(2,10000)
real torquereal(125,3),torquereciproco(27,3)
real forcareal(125,3),forcareciproca(27,3)

! REAL (ALLOCATABLES)
real, allocatable :: gpvetor(:)
real, allocatable :: campo(:)
real, allocatable :: y(:)
real, allocatable :: trap(:)
real, allocatable :: Di(:,:,:)
real, allocatable :: centro_massa(:,:)
real, allocatable :: aux1(:,:),aux2(:,:),aux3(:,:), aux4(:,:)
real, allocatable :: auxt(:,:),auxf(:,:)
real, allocatable :: contribuicao_self(:,:)
real, allocatable :: contribuicao_fisico(:,:)
real, allocatable :: contribuicao_reciproco(:,:)
real, allocatable :: contribuicoes(:,:)
real, allocatable :: FORCAS(:,:,:,:)
real, allocatable :: TORQUES(:,:,:,:)
real, allocatable :: FT(:,:,:)
real, allocatable :: Tt(:,:,:)
real, allocatable :: X(:,:,:)
real, allocatable :: W(:,:,:)
real, allocatable :: XI(:,:,:,:)
real, allocatable :: U(:,:,:)
real, allocatable :: V(:,:,:)
real, allocatable :: var(:,:)
real, allocatable :: errovar(:,:)
real, allocatable :: flutmag(:,:)
real, allocatable :: auxiliarcor(:,:)
real, allocatable :: nr(:)
real, allocatable :: tempototal(:)
real, allocatable :: aux_erro_vel(:,:,:)
real, allocatable :: aux_erro_var(:,:,:)
real, allocatable :: variancia(:,:,:)
real, allocatable :: flut(:,:,:,:)
real, allocatable :: auxcor(:,:,:)
real, allocatable :: errocor(:,:)
real, allocatable :: funcaor(:,:)
real, allocatable :: dif(:,:)
real, allocatable :: velmedia(:,:,:)
real, allocatable :: vmedia(:,:)
real, allocatable :: errovmedia(:,:)
real, allocatable :: fmedia(:,:)
real, allocatable :: errofmedia(:,:)
real, allocatable :: flut_tempo(:,:)
real, allocatable :: usistema(:,:)
real, allocatable :: DIAM(:,:)
real, allocatable :: auxflut(:,:,:)
real, allocatable :: diarand(:)
real, allocatable :: hidrodinamica_aux1(:,:)
real, allocatable :: hidrodinamica_aux2(:,:)
real, allocatable :: magtempo(:,:)
real, allocatable :: beta(:,:)
real, allocatable :: difusaoaux(:,:)

! LOGICAL
logical beating,posicao,fator,velocidade,estatistica, dipolo_ordenado
logical estatica,ordenado,torque,externo,gravidade,leito,mistura
logical gravadipolo,grafmag,browniano,tmagper,fmagper,ligaih
logical periodicidade,shear,oscillatory,agregado_inicial,continua
logical polidispersidade,printphi,oscilacampo,inertia,rotating
logical duffing,bifurcation, bifshear, ovito,magpart

! STRINGS
 character(3) rea_char
 character(23) linha1
 character(20) linha2
 character(25) linha3
 
end module variables
