#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#include "MDcuda.h"

/* Simulation global parameters. */

int N = 5000; // Number of particles.

// Lennard-Jones parameters in natural units!
double sigma = 1.,
       sigma_over_6 = sigma;


double epsilon = 1.,
       epsilon_times_4 = 4 * epsilon;

double m = 1.,
       kB = 1.,
       NA = 6.022140857e23,
       kBSI = 1.38064852e-23;  // m^2 * kg / (s^2 * K)

double L, // Size of the box, which will be specified in natural units.
       Tinit; // Initial Temperature in Natural Units.


double r[3][MAXPART], // Position array.
       v[3][MAXPART], // Velocity array.
       a[3][MAXPART]; // Acceleration array.
double *d_r,*d_v,*d_a,*d_result,*d_reduction;

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
    cudaMalloc(&d_a,sizeof(double)*MAXPART*3);
    cudaMalloc(&d_r,sizeof(double)*MAXPART*3);
    cudaMalloc(&d_v,sizeof(double)*MAXPART*3);
    cudaMalloc(&d_result,sizeof(double));

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

    cudaMemcpy(d_v, v, 3 * MAXPART * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_a, a, 3 * MAXPART * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, r, 3 * MAXPART * sizeof(double), cudaMemcpyHostToDevice);

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
    
    cudaFree(d_r);
    cudaFree(d_a);
    cudaFree(d_v);
    cudaFree(d_result);

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


__global__ void calculateSquareVelocity(double* v, int N, double* result) {
    extern __shared__ double tempResult[];
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    tempResult[threadIdx.x]=0;

    if (i < N) {
        for (int j = 0; j < 3; j++) {
            tempResult[threadIdx.x] += v[j * MAXPART + i] * v[j * MAXPART + i];
        }
    }

    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            tempResult[threadIdx.x] += tempResult[threadIdx.x + stride];
        }
        __syncthreads();
    }

    // Store the sum in d_pot
    if (threadIdx.x == 0) {
        result[blockIdx.x] = tempResult[0];
    }
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

    dim3 threadsPerBlock(256);
    dim3 blocksPerGrid((N + threadsPerBlock.x - 1) / threadsPerBlock.x);

    cudaMalloc(&d_reduction,sizeof(double)*blocksPerGrid.x);
    double result[blocksPerGrid.x];

    calculateSquareVelocity<<<blocksPerGrid, threadsPerBlock,threadsPerBlock.x*sizeof(double)>>>(d_v,N,d_reduction);
    
    cudaMemcpy(&result, d_reduction, sizeof(double)*blocksPerGrid.x, cudaMemcpyDeviceToHost);

    cudaFree(d_reduction);
    for (int i=1;i<blocksPerGrid.x;i++)
        result[0]+=result[i];
    return result[0] / N;
}

__global__ void calculateKinetic(double* v, int N, double m, double* result) {
    extern __shared__ double tempResult[];
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    tempResult[threadIdx.x]=0;

    if (i < N) 
        for (int j = 0; j < 3; j++)
            tempResult[threadIdx.x] += v[j * MAXPART + i] * v[j * MAXPART + i];


    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            tempResult[threadIdx.x] += tempResult[threadIdx.x + stride];
        }
        __syncthreads();
    }

    // Store the sum in d_pot
    if (threadIdx.x == 0) {
        result[blockIdx.x] = tempResult[0];
    }
    
}


/**
 * Calculate the kinetic energy of the system.
 *
 * This function computes the kinetic energy of the system by summing the kinetic energy of each particle.
 *
 * @return The total kinetic energy of the system.
 */
double Kinetic() { //Write Function here!
    dim3 threadsPerBlock(256);
    dim3 blocksPerGrid((N + threadsPerBlock.x - 1) / threadsPerBlock.x);

    cudaMalloc(&d_reduction,sizeof(double)*blocksPerGrid.x);
    double result[blocksPerGrid.x];

    calculateKinetic<<<blocksPerGrid, threadsPerBlock,threadsPerBlock.x*sizeof(double)>>>(d_v, N, m, d_reduction);

    cudaMemcpy(result, d_reduction, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_reduction);
    for (int i=1;i<blocksPerGrid.x;i++)
        result[0]+=result[i];

    return result[0];
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
__device__ double cudaLennardJonesForce(double InvdSqd, double InvdSqd3) {
    double InvdSqd4 = InvdSqd3 * InvdSqd,
           InvdSqd7 = InvdSqd4 * InvdSqd3;
    return 24 * (InvdSqd7 + InvdSqd7 - InvdSqd4);
}

// CUDA device function
__device__ double cudaPotentialEnergy(double sigma_over_6,double InvdSqd3) {
    double term2 = sigma_over_6 * InvdSqd3;
    double term1 = term2 * term2;
    return term1 - term2;
}


__global__ void cudaComputeAccelerationsAndPotentialVector(double* d_r, double* d_a, int N, double epsilon_times_4, double sigma_over_6, double* d_pot) {
    extern __shared__ double pot[];

    pot[threadIdx.x]=0;

    __syncthreads();

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N) {
        double ri[3] = {
            d_r[0 * MAXPART + i],
            d_r[1 * MAXPART + i],
            d_r[2 * MAXPART + i]
        };
        double rij[3],ai[3]={0};

        for (int j = 0; j < N; j++) {
            if (j==i) continue;
            for (int k = 0; k < 3; k++)
                rij[k] = ri[k] - d_r[k * MAXPART + j];

            double rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            double InvdSqd = 1.0 / rSqd;
            double InvdSqd3 = InvdSqd * InvdSqd * InvdSqd;

            double f = cudaLennardJonesForce(InvdSqd, InvdSqd3);

            for (int k = 0; k < 3; k++) 
                ai[k]+=rij[k] * f;
            if (j>i)
                pot[threadIdx.x] += cudaPotentialEnergy(sigma_over_6, InvdSqd3); // Accumulate potential energy
        }

        // Update accelerations in global memory
        for (int k = 0; k < 3; k++) {
            d_a[k * MAXPART + i] = ai[k];
        }
    }
    //reduction
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            pot[threadIdx.x] += pot[threadIdx.x + stride];
        }
        __syncthreads();
    }

    // Store the sum in d_pot
    if (threadIdx.x == 0) {
        d_pot[blockIdx.x] = pot[0];
    }

}


double computeAccelerationsAndPotentialVector() {
    cudaMemset(d_a, 0, sizeof(double) * 3 * MAXPART);
    

    dim3 threadsPerBlock(64);
    dim3 blocksPerGrid((N + threadsPerBlock.x) / threadsPerBlock.x);

    cudaMalloc(&d_reduction,sizeof(double)*blocksPerGrid.x);
    double totalPot[blocksPerGrid.x];

    cudaComputeAccelerationsAndPotentialVector<<<blocksPerGrid, threadsPerBlock,threadsPerBlock.x*sizeof(double)>>>(d_r, d_a, N, epsilon_times_4, sigma_over_6, d_reduction);

    cudaMemcpy(&totalPot, d_reduction, sizeof(double)*blocksPerGrid.x, cudaMemcpyDeviceToHost);

    cudaFree(d_reduction);
    for (int i=1;i<blocksPerGrid.x;i++)
        totalPot[0]+=totalPot[i];

    return 2 * epsilon_times_4 * totalPot[0];
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
    for(int i=0;i<3;i++)
        for (int j=0;j<N;j++)
            a[i][j]=0;

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

__global__ void updateMotion(double* r, double* v, double* a, int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N) {
        for (int j = 0; j < 3; j++) {
            r[j * MAXPART + i] += v[j * MAXPART + i] * dt + 0.5 * a[j * MAXPART + i] * dt * dt;
            v[j * MAXPART + i] += 0.5 * a[j * MAXPART + i] * dt;
        }
    }
}
__global__ void updateVelocities(double* r, double* v, double* a, double L, double dt, double m, int N, double* psum) {
    extern __shared__ double tempResult[];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    tempResult[threadIdx.x]=0;
    if (i < N) {
        for (int j = 0; j < 3; j++) {
            v[j * MAXPART + i] += 0.5 * a[j * MAXPART + i] * dt;

            // Elastic walls.
            if (r[j * MAXPART + i] < 0. || r[j * MAXPART + i] >= L) {
                v[j * MAXPART + i] *= -1.0; // Elastic walls.
                tempResult[threadIdx.x] += 2 * m * fabs(v[j * MAXPART + i]) / dt;  // Contribution to pressure from "left" walls.
            }
        }
    }

    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            tempResult[threadIdx.x] += tempResult[threadIdx.x + stride];
        }
        __syncthreads();
    }

    // Store the sum in d_pot
    if (threadIdx.x == 0) {
        psum[blockIdx.x] = tempResult[0];
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
    /*
    Compute accelerations from forces at current position 
    This call was removed (commented) for pedagogical reasons computeAccelerations();
    Update positions and velocity with current velocity and acceleration.
    */

    dim3 threadsPerBlock(256);
    dim3 blocksPerGrid((N + threadsPerBlock.x - 1) / threadsPerBlock.x);

    updateMotion<<<blocksPerGrid, threadsPerBlock>>>(d_r, d_v, d_a, N, dt);


    // Update accelerations from updated positions.
    *potentialEnergy = computeAccelerationsAndPotentialVector();


    cudaMalloc(&d_reduction,sizeof(double)*blocksPerGrid.x);
    double result[blocksPerGrid.x];

    updateVelocities<<<blocksPerGrid, threadsPerBlock,threadsPerBlock.x*sizeof(double)>>>(d_r, d_v, d_a,L,dt,m, N, d_reduction);

    cudaMemcpy(result, d_reduction, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_reduction);
    for (int i=1;i<blocksPerGrid.x;i++)
        result[0]+=result[i];
 
    return result[0] / (6 * L * L);
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