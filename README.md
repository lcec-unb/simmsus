# What is simmsus?

SIMMSUS is a research code written in FORTRAN that simulates the motion of a system of interacting particles. These particles can be simulated in different scenarios and may interact through different physical mechanisms. The code was initially developed for studying the physics of magnetic spherical particles in suspensions in order to better understand the properties of magnetic fluids (ferrofluids).

## History of the project

SIMMSUS is a research code written in FORTRAN that simulates the motion of a system of interacting particles. These particles can be simulated in different scenarios and may interact through different physical mechanisms. The code was initially developed for studying the physics of magnetic spherical particles in suspensions in order to better understand the properties of magnetic fluids (ferrofluids).

Since the beggining of the development, the code was designed to consider the effects of Brownian motion, long range dipole-dipole and hydrodynamic interactions. The first physics simulated through SIMMSUS was the study of the alignment of the particles in the direction of an applied steady-state magnetic field to see wether the code was capable of capturing the behavior of the equilibrium magnetization predicted by theoretical asymptotic models available on the literature (Ivanov and Kusnetsova, 2001).

The comparison between the predictions of SIMMSUS and the asymptotic models served as a first validation of the numerical code. After this initial validation, we have implemented several additional features on the code and have used SIMMSUS to explore different physical situations regarding the behavior of magnetic fluids. Many of the studies conducted through SIMMSUS have been published in prestigious academic Journals, such as Physics of Fluids, Journal of Magnetism and Magnetic Materials, Powder Technology and Mechanics Research Communications.

In its present version SIMMSUS is capable of simulating Brownian and non-Brownian suspensions and the user can turn on/off specific physical mechanisms in order to build a customized scenario. For example the user can mix magnetic particles with non-magnetic particles, turn on/off the gravitational field, apply an oscillatory or steady-state magnetic field over the particles, apply a simple or oscillatory shear over the particles, simulate different initial configurations (ordered, random or spherical distributions), simulate mono or polydisperse particles (with the same radius or with different radius) and the list goes on.

Regarding the application of time-dependent magnetic fields, we have conducted a rigorous study on the dynamical susceptibility response of ferrofluids using SIMMSUS and have validated the dynamical solution of the rotational motion of the dipole moments of the particles provided by the code by comparing the numerical solution with an asymptotic theoretical model (Berkov et al, 2009). The code seems to predict with excellent precision the dynamical behavior of the internal structure of ferrofluids.

## General structure of the source-code

The code contains 8 files, from which 7 are necessary to produce the executable file (program) after the compilation and 1 is a configuration file.
These are the following files:

1. simmsus.f90
2. input.f90
3. main.f90
4. subroutines.f90
5. statistics.f90
6. variables.f90
7. makefile

In order to compile the source-code you must open a Linux terminal and run the make command. For the first time you should run the make command twice. In the first time the compiler will produce the modules from the files variables.f90 and subroutines.f90 and the second run will use these recently created modules to produce the executable file simmsus.ex.

Since the origins of its earlier development, SIMMSUS has been tested and compiled exclusively on Linux using the Intel Fortran Compiler and up to this version of the manual we have no idea if it compiles or runs using other compilers or operational systems.

In order to run a simulation, the user will need only 2 files: the executable file (simmsus.ex) and the input configuration file simconfig.dat. The user can drag these two files to a new folder to run a new set of numerical simulations defined on the configuration file simconfig.dat. 

The configuration file is a text file with several questions that the user should answer in order to direct
the path that SIMMSUS should cross in order to produce a specific set of simulations for the intended physics. We will talk about this file later. Figure (1) shows how the source-code files are related to each other.

The file simmsus.f90 basically calls other files from the source code structure in a straightfoward order. The first call is related to the reading of the
simulation data. The file input.f90 reads the configuration file simconfig.dat and storage in logical, real and integer variables all the information necessary to run a set of simultaneous configurations. SIMMSUS may run a single simulation or perform several simultaneous numerical experiments varying the initial configuration of the particles and the set of random numbers used to emulate Brownian forces and torques. Therefore, the use may configure in the simconfig.dat file all the informations regarding all the simultaneous numerical experiments. 

We will talk about the simconfig.dat file later on this
README file. It is important to notice that the input.f90 file uses the variables module, written in the variables.f90 file. We will also talk about this module later on this manual.

After reading the configuration file, simmsus.f90 calls the main.f90. This file calls the necessary subroutines in a specific order in order to run the simulations. Here, the dynamical variables are allocated in the computer memory, the initial condition is generated and the main loop occurs. Here the physics is processed and several output files are created along the simulation
process. The main.f90 file uses the variables and subroutine modules. 

The last one contains all the subroutines used on the code. After the processing stage, simmsus.f90 checks if the user has chosen to perform a statistical analysis of the calculated results and if so calls the statistics.f90 file, which also uses the variables module.

The makefile contains the compilation instructions. Here the name of the source code files, the compilation command (and consequently the compiler) is specified and any optional compilation flags are also defined. We will talk about the general logic and structure of each of these files in the next subsections of this README file.

## simmsus.f90

