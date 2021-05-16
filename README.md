# Cahn-Hilliard-Equation-Solver-FDM-

!! filename = CH.f90
!! Vahid Attari
!! Created: 30 Feb. 2016
!! Modified: ....
!! Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
!!
!! Acknowledgements:  Based on Cahn-Hilliard 1965 paper
!!
!! Purpose:
!!   - Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases
!!     to self consistantly model the Spinodal Composition Phenomenon
!!   - Dynamic coupling is not provided ...
!!   
!! General Algorithm function:
!!
!!   1. Retrieve parameter data from file "parameters.dat"
!!   2. Assess thermodynamics of the associated system 
!!   3. Reads initial phase distribution from "phase.dat" file
!!   4. Calculate Phase Evolution with time integration
!!      -  Nucleate Phases
!!      -  Resolve boundary conditions 
!!         -- Periodic boundaries in all directions
!!      -  Solve differential equations via 9-stencil finite difference
!!      -  Update phase information and concentration data
!!
!! Compilation instructions: >> make
!!    - Manual: >>  ifort -o a.out CH.f90
!!
!! Execution: >> ./a.out 
!!                                     
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!====================================================================================

!   This code simulates the early stage of Spinodal Decomposition...

!   Cahn- Hilliard solver 
!   with periodic boundary conditions.
!   Nonlinear term: f(u) = u - u**3
