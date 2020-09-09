@def title = "CFD in Practice"
@def mintoclevel=2
@def maxtoclevel=3

# CFD in Practice 

\tableofcontents <!-- you can use \toc as well -->

## Introduction

* [Introduction and Review](/chap1/)
* [Mesh Generation](/chap2/)
* [Discontinuous Galerkin Method](/chap4/)
* [Incompressible Solver](/Incompressible/)

Lagrangian approach:
* vortex method
* smoothed particle hydrodynamics (SPH)
* Lattice Boltzmann method (LBM)

Eulerian approach:
* finite difference
* finite volume
* finite element

## Experiences

* Half indices are nice conceptually, but an implementation obviously should use integer indices. A standard convention is needed in each code. Different programming languages have different indices flexibility and conventions.
* Sticking to dense 3D arrays for your first fluid solver, and prototyping a sparse grid simulation in 2D before doing it in 3D.
* That last point bears repeating, and emphasizing, redundantly and multiple times as necessary. Always prototype a 3D solver in 2D first. Even once you get to a working 3D solver, keep your 2D prototype around and up-to-date so that you can easily jump back to it when working out new features or old bugs. In fact, if it's at all meaningful (which for incompressible flow is not always the case) prototype in 1D before 2D.

## Programming Languages

The languages involved in this note are:
* Julia
* C++
* MATLAB
* Python
* Fortran

The original class demo are given in Python; I wrote my projects originally in MATLAB and Fortran; later they were reimplemented in Julia; C++ is used for MAC grid from the book *Fluid Simulation for Computer Graphics*.

## References

1. CFD notes for AE423, AE523 and AE623.
2. *Fluid Simulation for Computer Graphics*
