program simmsus
use variables

implicit none
  real    :: ti, tf, tpros


print *,'*****************************************'
print *,'*                                       *'
print *,'*               SIMMSUS                 *'
print *,'*  SIMULATION OF MAGNETIC SUSPENSIONS   *'
print *,'*_______________________________________*'
print *,'*                                       *'
print *,'*   PROF. RAFAEL GABLER GONTIJO, PhD    *'
print *,'*_______________________________________*'
print *,'*                                       *'
print *,'*      IN DEVELOPMENT SINCE 2009        *'
print *,'*_______________________________________*'
print *,'*                                       *'
print *,'*       LAST UPDATE: 12/09/2025         *'
print *,'*                                       *'
print *,'*****************************************'
print *,''
print *,'*****************************************'
print *,'*                                       *'
print *,'*        Numerical simulation of        *'
print *,'* magnetic suspensions of hard spheres  *'  
print *,'*_______________________________________*'
print *,'*                                       *'
print *,'*    Langevin and Stokesian Dynamics    *'
print *,'*                                       *'
print *,'*****************************************'
print *,''

! Collecting simulation data
 call input
 
! Start to count the simulation time
 call cpu_time(ti)
 
! Calling the main subroutine responsible for executing the simulation
 call main
 
! Calling the statistical analysis subroutine
if(estatistica) then
 call statistics
end if

! Stops counting the simulation time
 call cpu_time(tf)
 tpros=tf-ti
 
print *, 'TOTAL SIMULATION TIME:',tpros,'SECONDS'
write(*,*) ''
end program simmsus
