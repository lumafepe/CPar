#pragma once

#include <stdio.h>
#include <immintrin.h>

/* Macros namespace for vectorization using AVX compiler intrinsics. */

#define vector __m256d // 256 bits available, we're using 4 packet doubles.

/* Set, load and store operations */
#define load _mm256_loadu_pd
#define set _mm256_set1_pd
#define store _mm256_storeu_pd

/* Arithmetic operations */
#define add _mm256_add_pd // a + b
#define sub _mm256_sub_pd // a - b
#define mul _mm256_mul_pd // a * m
#define div _mm256_div_pd // a / b

/* Custom operations */
#define mul3(a, b, c) _mm256_mul_pd(a, _mm256_mul_pd(b, c)) // Multiplies 3 vectors.
#define mul4(a, b, c, d) _mm256_mul_pd(_mm256_mul_pd(a, b), _mm256_mul_pd(c, d)) // Multiplies 4 vectors.


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
void getTitle(char *prefix, char *tfn, char *ofn, char *afn);
void getGasType(char *atomType);
void getAtomProperties(double *VolFac, double *TempFac, double *PressFac, double *timefac, char* atype);
void getParameters(double *rho, double *Vol, double VolFac, double TempFac);

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
 * @param rSqd The squared distance between particles.
 *
 * @return The Lennard-Jones force.
 */
vector lennardJonesForceVector(double rSqd);

/**
 * Calculate the Lennard-Jones force between particles at a given distance.
 *
 * This function calculates the Lennard-Jones force between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the force.
 *
 * @param rSqd The squared distance between particles.
 *
 * @note Uses AVX compiler intrinsics to perform the calculation.
 *
 * @return The Lennard-Jones force.
 */
vector lennardJonesForceVector(vector rSqd);

/**
 * Calculate the Lennard-Jones potential energy at a given distance.
 *
 * This function calculates the Lennard-Jones potential energy between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the potential energy.
 *
 * @param rSqd The squared distance between particles.
 *
 * @return The Lennard-Jones potential energy.
 */
double potentialEnergy(double rSqd);

/**
 * Calculate the Lennard-Jones potential energy at a given distance.
 *
 * This function calculates the Lennard-Jones potential energy between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the potential energy.
 *
 * @param rSqd The squared distance between particles.
 *
 * @note Uses AVX compiler intrinsics to perform the calculation.
 *
 * @return The Lennard-Jones potential energy.
 */
vector potentialEnergyVector(vector rSqd);

/**
 * Sums the elements of a vector and returns the result.
 *
 * This function takes a vector and calculates the sum of its elements.
 * It first copies the elements of the input vector to a temporary array
 * and then iterates through the array to compute the sum.
 *
 * @param vect The input vector to be summed.
 *
 * @return The sum of the elements in the vector.
 */
double sumVector(vector vect);

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
