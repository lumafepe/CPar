#pragma once

#include <stdio.h>
#include <immintrin.h>



#define MAXPART 5004 // Maximum array size.
/* Class for SIMD operations */


/**
 * Represents properties of an atom used in a molecular dynamics simulation.
 */
typedef struct AtomProperties {
    double VolFac;   // Volume scaling factor.
    double PressFac; // Pressure scaling factor.
    double TempFac;  // Temperature scaling factor.
    double timefac;  // Time scaling factor.
} AtomProperties;

/**
 * Set of helper functions used on main().
 */

/**
 * Retrieve simulation title.
 *
 * This function retrieves via user input the simulation title and prepends it to
 * the filenames to be outputted as a result of the simulation.
 *
 * @param prefix Output variable that will contain the simulation title.
 * @param tfn Trajectory output filename.
 * @param ofn Generic output filename.
 * @param afn Average values output filename.
 */
void getTitle(char *prefix, char *tfn, char *ofn, char *afn);

/**
 * Retrieve preferred gas type.
 *
 * This function shows a menu with the possible values for the gas parameter
 * and retrieves via user input the chosen one.
 *
 * @note If no value is provided by the user, it defaults to Argon.
 *
 * @param atomType Output variable that will contain the gas type.
 */
void getGasType(char *atomType);

/**
 * Get atom properties from the chosen gas.
 *
 * This function gets the atom properties corresponding the chosen atom type.
 *
 * @param VolFac Output variable, volume scaling factor.
 * @param TempFac Output variable, pressure scaling factor.
 * @param PressFac Output variable, temperature scaling factor.
 * @param timefac Output variable, time scaling factor.
 *
 * @param atype Gas chosen.
 */
void getAtomProperties(double *VolFac, double *TempFac, double *PressFac, double *timefac, char* atype);

/**
 * This function is used to get the simulation parameters from the user.
 *
 * @param rho Pointer to the number density in moles/m^3.
 * @param Vol Pointer to the volume of the gas.
 * @param VolFac Volume conversion factor.
 * @param TempFac Temperature conversion factor.
 */
void getParameters(double *rho, double *Vol, double VolFac, double TempFac);

/**
 * @brief This function is used to print the results of the simulation.
 *
 * @param tfn Pointer to the trajectory file name.
 * @param ofn Pointer to the output file name.
 * @param afn Pointer to the averages file name.
 * @param Tavg Average temperature in Kelvin.
 * @param Pavg Average pressure in Pascal.
 * @param gc Gas constant in J * mol^-1 K^-1.
 * @param VolNat Total volume in m^3.
 * @param Z The compressibility (unit-less).
 */
void printResults(char *tfn, char *ofn, char *afn, double Tavg, double Pavg, double gc, double VolNat, double Z);

/**
 * Initialize particle positions and velocities.
 *
 * This function initializes the positions of particles in a cubic lattice with spacing 'pos' and
 * also calls the 'initializeVelocities' function to set the initial velocities.
 *
 * @return None, but updates the positions and velocities of particles.
 */
void initialize();

/**
 * Calculate the mean squared velocity of particles in the system.
 *
 * This function computes the average squared velocity of particles by summing the squared velocities
 * of all particles and dividing by the total number of particles (N).
 *
 * @return The mean squared velocity of particles.
 */
double MeanSquaredVelocity();

/**
 * Calculate the kinetic energy of the system.
 *
 * This function computes the kinetic energy of the system by summing the kinetic energy of each particle.
 *
 * @return The total kinetic energy of the system.
 */
double Kinetic();

/**
 * Set accelerations to zero for all particles.
 *
 * This function initializes the accelerations for all particles to zero. It is typically used
 * at the beginning of a simulation to ensure a clean start for the computation of accelerations.
 *
 * @return None, but updates the accelerations in the 'a' array.
 */
void setAccelerationToZero();

/**
 * Calculate the Lennard-Jones force between particles at a given distance.
 *
 * This function calculates the Lennard-Jones force between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the force.
 *
 * @param InvrSqd The 1 / distance² between particles. 
 *
 * @param InvrSqd3 The 1 / distance⁶ between particles.
 *
 * @return The Lennard-Jones force.
 */
double lennardJonesForceVector(double rSqd,double InvrSqd3);

/**
 * Calculate the Lennard-Jones potential energy at a given distance.
 *
 * This function calculates the Lennard-Jones potential energy between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the potential energy.
 *
 * @param InvrSqd3 The 1\ distance⁶ between particles.
 *
 * @return The Lennard-Jones potential energy.
 */
double potentialEnergy(double InvrSqd3);


/**
 * Perform non-SIMD operations on arrays for molecular dynamics simulation.
 *
 * This function calculates forces and potential energy contributions between particles without utilizing
 * SIMD (Single Instruction, Multiple Data) operations. It processes particles one by one until the processing
 * index is aligned with a multiple of four.
 *
 * @param i - Index for the central particle.
 * @param j - Starting index for processing elements.
 * @param updates - 2D array representing updates to particles' accelerations.
 *
 * @note This function assumes that the arrays `r`, `v`, and `a` are globally defined.
 *
 * @return double - The total potential energy calculated during the simulation.
 */
inline double loopNonSIMD(int i, int j, double updates[][MAXPART]);



/**
 * Calculates accelerations and the Lennard-Jones potential energy using AVX instructions.
 *
 * This function computes the accelerations of particles and the Lennard-Jones potential energy
 * for a system of particles. It uses AVX (Advanced Vector Extensions) instructions for
 * vectorized calculations to improve performance. The function updates accelerations based on
 * the forces and calculates the potential energy. It handles both vectorized and non-vectorized
 * calculations.
 *
 * @return The total Lennard-Jones potential energy for the system.
 */
double computeAccelerationsAndPotentialVector();

/**
 * Calculates the accelerations of particles based on the Lennard-Jones potential.
 *
 * This function computes the accelerations of each particle in a system of particles
 * using the Lennard-Jones potential. It calculates the forces between particles and
 * then computes the accelerations based on the forces and particle masses. The function
 * handles vectorized operations for efficiency.
 */
void computeAccelerations();

/**
 * Perform a Velocity Verlet time integration step and calculate the potential energy
 * and pressure from elastic collisions with walls.
 *
 * This function performs a Velocity Verlet time integration step, updating particle positions
 * and velocities, and computes the potential energy and pressure contributions due to elastic
 * collisions with walls.
 *
 * @param dt              Time step size for the integration.
 * @param potentialEnergy Pointer to a variable to store the calculated potential energy.
 *
 * @return The pressure contribution from elastic collisions with walls.
 */
double VelocityVerletAndPotential(double dt, double* potentialEnergy);

/**
 * Initializes particle velocities for the simulation.
 *
 * This function sets the velocities of N particles in a simulation according to
 * specific procedures, ensuring the center of mass velocity is zero and that the
 * average velocity is scaled to match the desired initial temperature (Tinit).
 *
 * @note The function assumes that necessary variables like N, m, v, and Tinit
 * are declared and initialized before calling this function.
 *
 * @param None
 *
 * @return None
 */
void initializeVelocities();

/**
 * Generate a random number from a Gaussian distribution (mean=0, variance=1).
 *
 * This function implements the Box-Muller transform to generate random numbers
 * from a standard Gaussian distribution. It keeps track of a previously generated
 * value for efficiency and generates a new one every other call.
 *
 * @note The function uses the 'rand()' function for random number generation,
 * so make sure to seed the random number generator before using this function.
 *
 * @return A random number from a Gaussian distribution with a mean of 0 and a
 *         variance of 1.
 */
double gaussdist();
