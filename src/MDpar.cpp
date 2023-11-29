#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "MDpar.h"

/* Vector class implementation. */

Vector::Vector() {
    value = _mm256_setzero_pd();
}

Vector::Vector(double initialValue) {
    value = _mm256_set1_pd(initialValue);
}

Vector::Vector(__m256d initialValue) {
    value = initialValue;
}

Vector::Vector(double a, double b, double c, double d) {
    value = _mm256_set_pd(d,c,b,a);
}

Vector::Vector(double* array) {
    value = _mm256_load_pd(array);
}

void Vector::store(double* array) const {
    _mm256_store_pd(array,value);
}

__m256d Vector::getValue() const {
    return value;
}

Vector Vector::operator+(const Vector& other) const {
    return Vector(_mm256_add_pd(value,other.getValue()));
}

Vector Vector::operator-(const Vector& other) const {
    return Vector(_mm256_sub_pd(value,other.getValue()));
}

Vector Vector::operator*(const Vector& other) const {
    return Vector(_mm256_mul_pd(value,other.getValue()));
}

Vector Vector::operator/(const Vector& other) const {
    return Vector(_mm256_div_pd(value,other.getValue()));
}

double Vector::sum() const {
    return value[0] + value[1] + value[2] + value[3];
}

/* Simulation global parameters. */

int N = 5000; // Number of particles.

// Lennard-Jones parameters in natural units!
double sigma = 1.,
       sigma_over_6 = sigma;

Vector V1(1.0),
       V24(24.0),
       V48(48.0),
       Vsigma6(sigma_over_6);

double epsilon = 1.,
       epsilon_times_4 = 4 * epsilon;

double m = 1.,
       kB = 1.,
       NA = 6.022140857e23,
       kBSI = 1.38064852e-23;  // m^2 * kg / (s^2 * K)

double L, // Size of the box, which will be specified in natural units.
       Tinit; // Initial Temperature in Natural Units.


double r[3][MAXPART] __attribute__((aligned(32))), // Position array.
       v[3][MAXPART] __attribute__((aligned(32))), // Velocity array.
       a[3][MAXPART] __attribute__((aligned(32))); // Acceleration array.

const char* fileHeaders[] = {
    "  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n",
    "  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n",
    " --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n"
}; // Output file table headers.


/**
 * Molecular Dynamics simulation of a gas system.
 *
 * This function represents the main program for a Molecular Dynamics simulation of a gas system.
 * It initializes simulation parameters, performs the simulation, and calculates various properties
 * such as temperature, pressure, kinetic energy, potential energy, gas constant, and compressibility.
 * The simulation results are written to output files, and averages of temperature and pressure
 * are calculated and reported.
 *
 * @return 0 on successful execution.
 */
int main()
{
    double Pavg = 0.,
           Tavg = 0.;

    double dt, Temp, Press;
    double KE, PE, mvs, gc, Z;

    // Get and set title for calculation while also creating the respective output files.
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];
    getTitle(prefix, tfn, ofn, afn);

    FILE *tfp = fopen(tfn,"w"), // The MD trajectory, coordinates of every particle at each time step.
         *ofp = fopen(ofn,"w"), // Output of other quantities (T, P, gc, etc.) at every time step.
         *afp = fopen(afn,"w"); // Average T, P, gc, etc. from the simulation.

    fprintf(tfp,"%i\n", N); // Write down the number of particles to the trajectory file.

    // Get and set gas atom for the calculation (defaults to Argon).
    char atype[10];
    getGasType(atype);

    // Get and set atom properties.
    double VolFac, TempFac, PressFac, timefac;
    getAtomProperties(&VolFac, &TempFac, &PressFac, &timefac, atype);

    // Get and set simulation parameters, like temperature and density.
    double Vol, rho;
    getParameters(&rho, &Vol, VolFac, TempFac);

    int NumTime; // Amount of time-steps for the simulation.
    if (strcmp(atype, "He") == 0) {
        dt = 0.2e-14 / timefac; // dt in natural units of time (2 f.s.)
        NumTime = 50000;
    }
    else {
        dt = 0.5e-14 / timefac; // dt in natural units of time (5 f.s)
        NumTime = 200;
    }

    // Put all the atoms in simple crystal lattice and give them random velocities
    // that corresponds to the initial temperature we have specified.
    initialize();

    // Based on their positions, calculate the initial intermolecular forces.
    // The accelerations of each particle will be defined from the forces and their
    // mass, and this will allow us to update their positions via Newton's law.
    computeAccelerations();

    fprintf(ofp, "%s", fileHeaders[0]);
    printf("\n\nProgress: [");

    // Start the simulation.
    int i, tenp = floor(NumTime / 10);
    for (i = 0; i < NumTime + 1; i++) {

        for (int u = 1; u <= 10; u++) {
            if (i == u * tenp) {
                if (u == 10)
                    printf(" 100 ]\n\n");
                else
                    printf(" %d |", 10 * u);
                break;
            }
        }
        fflush(stdout);

        Press = VelocityVerletAndPotential(dt, &PE); // Update positions and velocities according to Newton's Law.
        Press *= PressFac; // Compute the pressure as the sum of momentum changes from wall collisions / time-step,

        mvs = MeanSquaredVelocity(); // Calculate velocity squared.
        KE = Kinetic(); // Calculate Kinetic energy.
        Temp = m * mvs / (3 * kB) * TempFac; // Calculate temperature from kinetic theory.
        gc = NA * Press * (Vol * VolFac) / (N * Temp); // Calculate gas constant.
        Z  = Press * (Vol * VolFac) / (N * kBSI * Temp); // Calculate compressibility.

        Tavg += Temp;
        Pavg += Press;

        fprintf(
            ofp,
            "  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",
            (i * dt * timefac), Temp, Press, KE, PE, KE + PE
        ); // Write values to output file.
    }

    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here.
    Pavg /= NumTime;
    Tavg /= NumTime;

    Z = Pavg * (Vol * VolFac) / (N * kBSI * Tavg);
    gc = NA * Pavg * (Vol * VolFac) / (N * Tavg);

    fprintf(afp,"%s", fileHeaders[1]);
    fprintf(afp,"%s", fileHeaders[2]);
    fprintf(
        afp,
        "  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e         %i\n",
        i * dt * timefac, Tavg, Pavg, gc, Z, Vol * VolFac, N
    );

    printResults(tfn, ofn, afn, Tavg, Pavg, gc, Vol * VolFac, Z);
    
    fclose(tfp); fclose(ofp); fclose(afp);
    return 0;
}

/**
 * Initialize particle positions and velocities.
 *
 * This function initializes the positions of particles in a cubic lattice with spacing 'pos' and
 * also calls the 'initializeVelocities' function to set the initial velocities.
 *
 * @return None, but updates the positions and velocities of particles.
 */
void initialize() {

    // Number of atoms in each direction.
    int n = int(ceil(pow(N, 1.0/3)));
    
    // Spacing between atoms along a given direction.
    double pos = L / n;
    
    // Index for number of particles assigned positions.
    int p = 0;

    // Initialize positions.
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            for (int k=0; k < n; k++) {
                if (p<N) {
                    r[0][p] = (i + 0.5) * pos;
                    r[1][p] = (j + 0.5) * pos;
                    r[2][p] = (k + 0.5) * pos;
                }
                p++;
            }
        }
    }
    
    // Call function to initialize velocities.
    initializeVelocities();    
}   


/**
 * Calculate the mean squared velocity of particles in the system.
 *
 * This function computes the average squared velocity of particles by summing the squared velocities
 * of all particles and dividing by the total number of particles (N).
 *
 * @return The mean squared velocity of particles.
 */
double MeanSquaredVelocity() { 

    double s = 0.;
    for (int j = 0; j < 3;j++)
        for (int i = 0; i < N; i++)   
            s += v[j][i] * v[j][i];

    return s / N;
}

/**
 * Calculate the kinetic energy of the system.
 *
 * This function computes the kinetic energy of the system by summing the kinetic energy of each particle.
 *
 * @return The total kinetic energy of the system.
 */
double Kinetic() { //Write Function here!  
    
    double v2, kin = 0.;

    for (int i = 0; i < N; i++) {
        v2 = 0.;
        for (int j=0; j<3; j++) v2 += v[j][i]*v[j][i];
        kin += 0.5 * m * v2;
    }
    
    return kin;
}


/**
 * Set accelerations to zero for all particles.
 *
 * This function initializes the accelerations for all particles to zero. It is typically used
 * at the beginning of a simulation to ensure a clean start for the computation of accelerations.
 *
 * @return None, but updates the accelerations in the 'a' array.
 */
void setAccelerationToZero() {
    Vector zero = Vector();

    for (int k = 0; k < 3; k++) {
        int i = 0;
        for (; i < N - (N % 4); i += 4) zero.store(&a[k][i]); /* vector acceleration array. */
        for (; i < N; i++) a[k][i] = 0; /* Double acceleration array. */
    }
}

/**
 * Calculate the Lennard-Jones force between particles at a given distance.
 *
 * This function calculates the Lennard-Jones force between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the force.
 *
 * @param InvrSqd The 1 \ distance² between particles.
 * 
 * @param InvrSqd3 The 1 \ distance⁶ between particles.
 *
 * @return The Lennard-Jones force.
 */
double lennardJonesForce(double InvrSqd,double InvrSqd3) {
    double InvrSqd4 = InvrSqd3 * InvrSqd,
           InvrSqd7 = InvrSqd4 * InvrSqd3;

    return 24 * (InvrSqd7 + InvrSqd7 - InvrSqd4);
}

/* Vectorized version of the above function. */
Vector lennardJonesForceVector(Vector InvrSqdV,Vector InvrSqdV3) {

    Vector InvrSqdV4 = InvrSqdV * InvrSqdV3; // InvrSqd ^ 4
    return InvrSqdV4 * ((V48*InvrSqdV3)-V24); // 24 * (2 * InvrSqd7 - InvrSqd4)
}

/**
 * Calculate the Lennard-Jones potential energy at a given distance.
 *
 * This function calculates the Lennard-Jones potential energy between particles based on the squared distance
 * between them. It uses the 12-6 Lennard-Jones potential to compute the potential energy.
 *
 * @param InvrSqd3 The 1/ distance⁶ between particles.
 *
 * @return The Lennard-Jones potential energy.
 */
double potentialEnergy(double InvrSqd3) {
    double term2 = sigma_over_6 * InvrSqd3;
    double term1 = term2 * term2;

    return term1 - term2;
}

/* Vectorized version of the above function. */
Vector potentialEnergyVector(Vector InvrSqdV3) {

    Vector term2 = Vsigma6 * InvrSqdV3; // sigma_over_6 / (rSqd ^ 3)

    Vector term1 = term2 * term2; // term2 * term2

    return term1 - term2; // term1 - term2
}


/**
 * Perform non-SIMD operations on vectors for molecular dynamics simulation.
 *
 * @param i - Index for particles.
 * @param j - Starting index for processing elements.
 * @param updates - 2D array representing updates to particles.
 *
 * @return double - The total potential energy calculated during the simulation.
 */
inline double loopNonSIMD(int i, int j, double updates[][MAXPART]) {

    double ri[3] = {
            r[0][i],
            r[1][i],
            r[2][i]
    };
    double ai[3] = {0,0,0};
    double rSqd, f, rij[3], pot = 0;

    for (; j % 4 != 0; j++) {

        for (int k = 0; k < 3; k++)
            rij[k] = ri[k] - r[k][j]; // Calculate relative positions and squared distances

        rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // Calculate inverse squared distance and its cube
        double InvrSqd = 1 / rSqd;
        double InvrSqd3 = InvrSqd * InvrSqd * InvrSqd;

        f = lennardJonesForce(InvrSqd,InvrSqd3);

        // Update accelerations and apply updates to the array
        for (int k = 0; k < 3; k++){
            rij[k] = rij[k] * f;
            ai[k] = ai[k] + rij[k];
            updates[k][j] = updates[k][j] - rij[k];
        }
        pot = pot + potentialEnergy(InvrSqd3); // Accumulate potential energy
    }

    for (int k = 0; k < 3; k++)
        updates[k][i] += ai[k]; // Accumulate acceleration contributions into the updates array

    return pot;
}

/**
 * Perform SIMD operations on vectors for molecular dynamics simulation.
 *
 * @param i - Index for particles.
 * @param j - Starting index for processing elements.
 * @param updates - 2D array representing updates to particles.
 *
 * @return Vector - The total potential energy calculated during the simulation.
 */
inline Vector loopSIMD(int i, int j, double updates[][MAXPART]) {

    Vector ri[3] = {
            Vector(r[0][i]),
            Vector(r[1][i]),
            Vector(r[2][i])
    };

    Vector ai[3] = {
            Vector(),
            Vector(),
            Vector()
    };

    Vector rSqd, f, rij[3], pot = Vector();

    for (; j < N; j += 4) {

        for (int k = 0; k < 3; k++)
            rij[k] = ri[k] - Vector(&r[k][j]);

        rSqd = rij[0] * rij[0] + rij[1] * rij[1] +rij[2] * rij[2];

        // Calculate inverse squared distance and its cube
        Vector InvrSqd = V1 / rSqd;
        Vector InvrSqd3 = InvrSqd * InvrSqd * InvrSqd;

        // Calculate forces using Lennard-Jones potential
        f = lennardJonesForceVector(InvrSqd,InvrSqd3);

        // Update accelerations and apply updates to the array
        for (int k = 0; k < 3; k++) {
            rij[k] = rij[k] * f;
            ai[k] = ai[k] + rij[k];
            (Vector(&updates[k][j]) - rij[k]).store(&updates[k][j]);
        }

        // Accumulate potential energy contributions
        pot = pot + potentialEnergyVector(InvrSqd3);
    }

    for (int k = 0; k < 3; k++)
        updates[k][i] += ai[k].sum(); // Accumulate acceleration contributions into the updates array

    return pot;
}

// Sum vectors in parallel.
#pragma omp declare reduction(+ : Vector : omp_out = omp_out + omp_in) initializer(omp_priv = Vector())

/**
 * Calculates accelerations and the Lennard-Jones potential energy using AVX instructions.
 *
 * This function computes the accelerations of particles and the Lennard-Jones potential energy
 * for a system of particles. It uses AVX (Advanced vector Extensions) instructions for
 * vectorized calculations to improve performance. The function updates accelerations based on
 * the forces and calculates the potential energy. It handles both vectorized and non-vectorized
 * calculations.
 *
 * @return The total Lennard-Jones potential energy for the system.
 */
double computeAccelerationsAndPotentialVector() {

    Vector totalPotV;
    double totalPot = 0.;

    setAccelerationToZero();
    #pragma omp parallel reduction(+:a)
    {
        double updates[3][MAXPART] __attribute__((aligned(32)))={0};
        #pragma omp for schedule(dynamic) reduction(+:totalPot) reduction(+:totalPotV)
        for (int i = 0; i < N - 1; i++){
            totalPot = totalPot + loopNonSIMD(i,i + 1,updates);
            totalPotV = totalPotV + loopSIMD(i,((i + 4)/4)*4,updates);
        }
        #pragma omp simd aligned(updates:32)
        for (int i=0;i<N;i++) {
            a[0][i] += updates[0][i];
            a[1][i] += updates[1][i];
            a[2][i] += updates[2][i];
        }
    }
    return 2 * epsilon_times_4 * (totalPotV.sum()+ totalPot);
}

/**
 * Calculates the accelerations of particles based on the Lennard-Jones potential.
 *
 * This function computes the accelerations of each particle in a system of particles
 * using the Lennard-Jones potential. It calculates the forces between particles and
 * then computes the accelerations based on the forces and particle masses. The function
 * handles vectorized operations for efficiency.
 */
void computeAccelerations() {

    double ai[3], ri[3], rij[3],
           f, rSqd;

    setAccelerationToZero();

    for (int i = 0; i < N - 1; i++)
    {
        // Store the position of the particle i and set the acceleration to zero.
        for (int k = 0; k < 3; k++) {
            ri[k] = r[k][i];
            ai[k] = 0.;
        }

        // Handling non vectorizable values.
        for (int j = i + 1; j < N; j++) {
            for (int k = 0; k < 3; k++) rij[k] = ri[k] - r[k][j];
        
            rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
            double InvrSqd = 1/rSqd;
            double InvrSqd3 = 1/(rSqd * rSqd * rSqd);
            f = lennardJonesForce(InvrSqd,InvrSqd3);

            for (int k = 0; k < 3; k++){
                rij[k] *= f;
                ai[k] += rij[k];
                a[k][j] -= rij[k];
            }
        }

        // Sum + Reduce of all accelerations.
        for (int k = 0; k < 3; k++) 
            a[k][i] += ai[k];
    }
}


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
double VelocityVerletAndPotential(double dt, double *potentialEnergy) {

    double psum = 0.;

    /*
    Compute accelerations from forces at current position 
    This call was removed (commented) for pedagogical reasons computeAccelerations();
    Update positions and velocity with current velocity and acceleration.
    */

    for (int j = 0; j < 3; j++) {
        for (int i=0; i < N; i++) {
            r[j][i] += v[j][i] * dt + 0.5 * a[j][i] * dt * dt;
            v[j][i] += 0.5 * a[j][i] * dt;
        }

    }

    // Update accelerations from updated positions.
    *potentialEnergy = computeAccelerationsAndPotentialVector();

    // Update velocity with updated acceleration.
    for (int j = 0; j < 3; j++) {

        for (int i = 0; i < N; i++){
            v[j][i] += 0.5 * a[j][i] * dt;
        }
        
        // Elastic walls.
        for (int i = 0; i < N; i++) {
            if (r[j][i] < 0.|| r[j][i] >= L) {
                v[j][i] *= -1.; // Elastic walls.
                psum += 2 * m * fabs(v[j][i]) / dt;  // Contribution to pressure from "left" walls.
            }
        }
    }
    
    return psum / (6 * L * L);
}


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
void initializeVelocities() {
    
    int i, j;
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) v[j][i] = gaussdist(); // Pull a number from a Gaussian Distribution.
    }
    
    /*
    Compute center of mass velocity according to the following formula:
    Vcm = sum_i^N m*v_i / sum_i^N M
    */
    double vCM[3] = {0, 0, 0};
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) vCM[j] += m * v[j][i];
    }
    
    for (i = 0; i < 3; i++) vCM[i] /= N*m;

    /*
    Subtract out the center-of-mass velocity from the
    velocity of each particle... effectively set the
    center of mass velocity to zero so that the system does
    not drift in space!
    */
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) v[j][i] -= vCM[j];
    }
    
    /*
    Now we want to scale the average velocity of the system
    by a factor which is consistent with our initial temperature, Tinit
    */
    double vSqdSum = 0., lambda;
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) vSqdSum += v[j][i]*v[j][i];
    }
    
    lambda = sqrt(3 * (N-1) * Tinit / vSqdSum);
    
    for (i = 0; i < N; i++) {
        for (j=0; j<3; j++) v[j][i] *= lambda;
    }
}


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
double gaussdist() {

    static bool available = false;
    static double gset;

    double fac, rsq, v1, v2;

    if (!available) {
        
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        return v2 * fac;
    } else {
        
        available = false;
        return gset;
        
    }
}

/* Helper functions used to retrieve user input and show results. ================================================= */

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
void getTitle(char *prefix, char *tfn, char *ofn, char *afn) {

    printf("Welcome!\n");
    printf("Calculation title: ");

    scanf("%s", prefix);

    /* Trajectory file */
    strcpy(tfn, prefix);
    strcat(tfn, "_traj.xyz");

    /* Output file */
    strcpy(ofn, prefix);
    strcat(ofn, "_output.txt");

    /* Average file */
    strcpy(afn, prefix);
    strcat(afn, "_average.txt");

    printf("> Title is set to '%s'\n\n", prefix);
}

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
void getGasType(char *atomType) {

    printf("Which noble gas would you like to simulate ?\n");
    printf("+ Helium - He\n");
    printf("+ Neon - Ne\n");
    printf("+ Argon - Ar\n");
    printf("+ Krypton - Kr\n");
    printf("+ Xeon - Xe\n");

    printf("Atom type (defaults to Ar): ");
    scanf("%s", atomType);

    printf("Using the gas '%s'\n\n", atomType);
}

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
void getAtomProperties(double *VolFac, double *TempFac, double *PressFac, double *timefac, char* atype) {

    AtomProperties atomPropertiesMap[] = {
            {1.8399744000000005e-29, 8152287.336171632, 10.864459551225972, 1.7572698825166272e-12}, // He
            {2.0570823999999997e-29, 27223022.27659913, 40.560648991243625, 2.1192341945685407e-12}, // Ne
            {3.7949992920124995e-29, 51695201.06691862, 142.0950000000000, 2.09618e-12},             // Ar
            {4.5882712000000004e-29, 59935428.40275003, 199.1817584391428, 8.051563913585078e-13},   // Kr
            {5.4872e-29, 70527773.72794868, 280.30305642163006, 9.018957925790732e-13}               // Xe
    };

    int atomIndex = 2; // Default type is Argon.
    long unsigned i, len = sizeof(atomPropertiesMap) / sizeof(atomPropertiesMap[0]);
    for (i = 0; i < len; i++) {
        if (strcmp(atype, "He") == 0) { atomIndex = 0; break; }
        else if (strcmp(atype, "Ne") == 0) { atomIndex = 1; break; }
        else if (strcmp(atype, "Ar") == 0) { atomIndex = 2; break; }
        else if (strcmp(atype, "Kr") == 0) { atomIndex = 3; break; }
        else if (strcmp(atype, "Xe") == 0) { atomIndex = 4; break; }
    }

    *VolFac = atomPropertiesMap[atomIndex].VolFac;
    *TempFac = atomPropertiesMap[atomIndex].TempFac;
    *PressFac = atomPropertiesMap[atomIndex].PressFac;
    *timefac = atomPropertiesMap[atomIndex].timefac;
}

/**
 * This function is used to get the simulation parameters from the user.
 *
 * @param rho Pointer to the number density in moles/m^3.
 * @param Vol Pointer to the volume of the gas.
 * @param VolFac Volume conversion factor.
 * @param TempFac Temperature conversion factor.
 */
void getParameters(double *rho, double *Vol, double VolFac, double TempFac) {

    printf("You will now enter the simulation parameters:\n");
    printf("Initial temperature (Kº) of the gas: ");
    scanf("%lf", &Tinit);

    if (Tinit < 0.) {
        printf("Absolute temperature must be a positive number.\n");
        exit(EXIT_FAILURE);
    }

    printf("\nEnter the number density in moles/m^3: ");
    scanf("%lf", rho);

    Tinit /= TempFac; // Convert temperature into natural units.
    *Vol = N / (*rho * NA);
    *Vol /= VolFac; // Convert volume into natural units.

    // Limiting N to MAXPART for practical reasons.
    if (N >= MAXPART) {
        printf("\nThe maximum number of particles is %d, please adjust your input file accordingly.\n", MAXPART);
        exit(EXIT_FAILURE);
    }

    // Check to see if the volume makes sense - is it too small?
    // Remember VDW radius of the particles is 1 natural unit of length
    // and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    // will be initialized with an inter-particle separation equal to 2xVDW radius.
    if (*Vol < N) {
        printf("\nYour density is very high, simulations with density greater than 1 particle may diverge.\n");
        exit(EXIT_FAILURE);
    }

    L = pow(*Vol,(1./3)); // Length of the box in natural units.
}

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
void printResults(char *tfn, char *ofn, char *afn, double Tavg, double Pavg, double gc, double VolNat, double Z) {

    printf("> To animate your simulation, open the file '%s' with VMD after the simulation completes.\n", tfn);
    printf("> To analyze instantaneous data about your molecule, open the file '%s' with your favorite text editor.\n\n", ofn);

    printf("> The following thermodynamic averages will be computed and written to the file '%s':\n", afn);
    printf("    > AVERAGE TEMPERATURE (K):                 %15.5f\n", Tavg);
    printf("    > AVERAGE PRESSURE  (Pa):                  %15.5f\n", Pavg);
    printf("    > PV/nT (J * mol^-1 K^-1):                 %15.5f\n", gc);
    printf("    > PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n", 100 * fabs(gc - 8.3144598) / 8.3144598);
    printf("    > THE COMPRESSIBILITY (unit-less):          %15.5f \n", Z);
    printf("    > TOTAL VOLUME (m^3):                      %10.5e \n", VolNat);
    printf("    > NUMBER OF PARTICLES (unit-less):          %d \n", N);
}