#pragma once

#include <stdio.h>

//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();  
//  update positions and velocities using Velocity Verlet algorithm 
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
//  sets potentialEnergy to total potential energy from particle coordinates
double VelocityVerletAndPotential(double dt, int iter, FILE *fp,double* potentialEnergy);  
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
void computeAccelerations();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();