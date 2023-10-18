/*
 MD.c - a simple molecular dynamics program for simulating real gas properties of Lennard-Jones particles.
 
 Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 Electronic Contact:  foleyj10@wpunj.edu
 Mail Contact:   Prof. Jonathan Foley
 Department of Chemistry, William Paterson University
 300 Pompton Road
 Wayne NJ 07470
 
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "MD.h"


// Number of particles
int N;

// Lennard-Jones parameters in natural units!
double sigma = 1.;
double sigma_over_6 = sigma;

double epsilon = 1.;
double epsilon_times_4 = 4 * epsilon;

double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

// Size of box, which will be specified in natural units
double L;

// Initial Temperature in Natural Units
double Tinit;  //2;
// Vectors!
//
const int MAXPART=5001;
//  Position
double r[3][MAXPART] __attribute__((aligned(32)));
//  Velocity
double v[3][MAXPART] __attribute__((aligned(32)));
//  Acceleration
double a[3][MAXPART] __attribute__((aligned(32)));

// atom type
char atype[10];


int main()
{
    
    //  variable delcarations
    int i;
    double dt, Vol, Temp, Press, Pavg, Tavg, rho;
    double VolFac, TempFac, PressFac, timefac;
    double KE, PE, mvs, gc, Z;
    char trash[10000], prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *infp, *tfp, *ofp, *afp;
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  WELCOME TO WILLY P CHEM MD!\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
    scanf("%s",prefix);
    strcpy(tfn,prefix);
    strcat(tfn,"_traj.xyz");
    strcpy(ofn,prefix);
    strcat(ofn,"_output.txt");
    strcpy(afn,prefix);
    strcat(afn,"_average.txt");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  TITLE ENTERED AS '%s'\n",prefix);
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    /*     Table of values for Argon relating natural units to SI units:
     *     These are derived from Lennard-Jones parameters from the article
     *     "Liquid argon: Monte carlo and molecular dynamics calculations"
     *     J.A. Barker , R.A. Fisher & R.O. Watts
     *     Mol. Phys., Vol. 21, 657-673 (1971)
     *
     *     mass:     6.633e-26 kg          = one natural unit of mass for argon, by definition
     *     energy:   1.96183e-21 J      = one natural unit of energy for argon, directly from L-J parameters
     *     length:   3.3605e-10  m         = one natural unit of length for argon, directly from L-J parameters
     *     volume:   3.79499-29 m^3        = one natural unit of volume for argon, by length^3
     *     time:     1.951e-12 s           = one natural unit of time for argon, by length*sqrt(mass/energy)
     ***************************************************************************************/
    
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Edit these factors to be computed in terms of basic properties in natural units of
    //  the gas being simulated
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
    printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    scanf("%s",atype);
    
    if (strcmp(atype,"He")==0) {
        
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
        
    }
    else if (strcmp(atype,"Ne")==0) {
        
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
        
    }
    else if (strcmp(atype,"Ar")==0) {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        //strcpy(atype,"Ar");
        
    }
    else if (strcmp(atype,"Kr")==0) {
        
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
        
    }
    else if (strcmp(atype,"Xe")==0) {
        
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
        
    }
    else {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype,"Ar");
        
    }
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n                     YOU ARE SIMULATING %s GAS! \n",atype);
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
    scanf("%lf",&Tinit);
    // Make sure temperature is a positive number!
    if (Tinit<0.) {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }
    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;
    
    
    printf("\n\n  ENTER THE NUMBER DENSITY IN moles/m^3\n");
    printf("  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3\n");
    printf("  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3\n");
    
    scanf("%lf",&rho);
    
    N = 10*216;
    Vol = N/(rho*NA);
    
    Vol /= VolFac;
    
    //  Limiting N to MAXPART for practical reasons
    if (N>=MAXPART) {
        
        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
        
    }
    //  Check to see if the volume makes sense - is it too small?
    //  Remember VDW radius of the particles is 1 natural unit of length
    //  and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    //  will be initialized with an interparticle separation equal to 2xVDW radius
    if (Vol<N) {
        
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n",N,Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }
    // Vol = L*L*L;
    // Length of the box in natural units:
    L = pow(Vol,(1./3));
    
    //  Files that we can write different quantities to
    tfp = fopen(tfn,"w");     //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn,"w");     //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn,"w");    //  Average T, P, gc, etc from the simulation
    
    int NumTime;
    if (strcmp(atype,"He")==0) {
        
        // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
        dt = 0.2e-14/timefac;
        //  We will run the simulation for NumTime timesteps.
        //  The total time will be NumTime*dt in natural units
        //  And NumTime*dt multiplied by the appropriate conversion factor for time in seconds
        NumTime=50000;
    }
    else {
        dt = 0.5e-14/timefac;
        NumTime=200;
        
    }
    
    //  Put all the atoms in simple crystal lattice and give them random velocities
    //  that corresponds to the initial temperature we have specified
    initialize();
    
    //  Based on their positions, calculate the ininial intermolecular forces
    //  The accellerations of each particle will be defined from the forces and their
    //  mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();
    
    
    // Print number of particles to the trajectory file
    fprintf(tfp,"%i\n",N);
    
    //  We want to calculate the average Temperature and Pressure for the simulation
    //  The variables need to be set to zero initially
    Pavg = 0;
    Tavg = 0;
    
    
    int tenp = floor(NumTime/10);
    fprintf(ofp,"  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
    printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
    for (i=0; i<NumTime+1; i++) {
        
        //  This just prints updates on progress of the calculation for the users convenience
        if (i==tenp) printf(" 10 |");
        else if (i==2*tenp) printf(" 20 |");
        else if (i==3*tenp) printf(" 30 |");
        else if (i==4*tenp) printf(" 40 |");
        else if (i==5*tenp) printf(" 50 |");
        else if (i==6*tenp) printf(" 60 |");
        else if (i==7*tenp) printf(" 70 |");
        else if (i==8*tenp) printf(" 80 |");
        else if (i==9*tenp) printf(" 90 |");
        else if (i==10*tenp) printf(" 100 ]\n");
        fflush(stdout);
        
        
        // This updates the positions and velocities using Newton's Laws
        // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
        // which is a Kinetic Theory of gasses concept of Pressure
        Press = VelocityVerletAndPotential(dt, &PE);
        Press *= PressFac;
        
        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  Now we would like to calculate something about the system:
        //  Instantaneous mean velocity squared, Temperature, Pressure
        //  Potential, and Kinetic Energy
        //  We would also like to use the IGL to try to see if we can extract the gas constant
        mvs = MeanSquaredVelocity();
        KE = Kinetic();
        
        // Temperature from Kinetic Theory
        Temp = m*mvs/(3*kB) * TempFac;
        
        // Instantaneous gas constant and compressibility - not well defined because
        // pressure may be zero in some instances because there will be zero wall collisions,
        // pressure may be very high in some instances because there will be a number of collisions
        gc = NA*Press*(Vol*VolFac)/(N*Temp);
        Z  = Press*(Vol*VolFac)/(N*kBSI*Temp);
        
        Tavg += Temp;
        Pavg += Press;
        
        fprintf(ofp,"  %.11e  %.11e  %.11e %.11e  %.11e  %.11e \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
        
        
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);
    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
    fprintf(afp,"  %.11e  %.11e       %.11e     %.11e       %.11e        %.11e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
    printf("\n  AVERAGE TEMPERATURE (K):                 %.11e\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %.11e\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %.11e\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %.11e\n",100*fabs(gc-8.3144598)/8.3144598);
    printf("\n  THE COMPRESSIBILITY (unitless):          %.11e \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %.11e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);
    
    fclose(tfp);
    fclose(ofp);
    fclose(afp);
    
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
    vector zero = set(0.);

    for (int k = 0; k < 3; k++) {
        int i = 0;
        for (; i < N - (N % 4); i += 4) store(&a[k][i], zero); /* vector acceleration array. */
        for (; i < N; i++) a[k][i] = 0; /* Double acceleration array. */
    }
}

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
double lennardJonesForce(double rSqd) {
    double InvrSqd = 1. / rSqd,
           InvrSqd4 = InvrSqd * InvrSqd * InvrSqd * InvrSqd,
           InvrSqd7 = InvrSqd4 * InvrSqd * InvrSqd * InvrSqd;

    return 24 * (2 * InvrSqd7 - InvrSqd4);
}

/* Vectorized version of the above function. */
vector lennardJonesForceVector(vector rSqd) {

    vector InvrSqd = div(set(1.), rSqd); // 1 / rSqd
    vector InvrSqd4 = mul4(InvrSqd, InvrSqd, InvrSqd, InvrSqd); // InvrSqd ^ 4
    vector InvrSqd7 = mul4(InvrSqd4, InvrSqd, InvrSqd, InvrSqd); // InvrSqd ^ 7

    return mul(
        set(24.),
        sub(
            add(
                InvrSqd7,
                InvrSqd7
            ),
            InvrSqd4
        )
    ); // 24 * (2 * InvrSqd7 - InvrSqd4)
}

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
double potentialEnergy(double rSqd) {
    double term2 = sigma_over_6 / (rSqd * rSqd * rSqd);
    double term1 = term2 * term2;

    return term1 - term2;
}

/* Vectorized version of the above function. */
vector potentialEnergyVector(vector rSqd) {

    vector term2 = div(
        set(sigma_over_6),
        mul3(rSqd, rSqd, rSqd)
    ); // sigma_over_6 / (rSqd ^ 3)

    vector term1 = mul(term2, term2); // term2 * term2

    return sub(term1, term2); // term1 - term2
}

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
double sumVector(vector vect) {
    double sum = 0.0,
           vectArray[4];

    store(vectArray, vect); // vectArray[0] = vect[0] ; ... ; vectArray[3] = vect[3]

    for (int i = 0; i < 4; i++) sum += vectArray[i];
    return sum;
}

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

    double ai[3], ri[3], rij[3],
           f, rSqd, pot = 0.0; // Non vectorized variables.

    vector rijV[3], rijVsqd[3], aiV[3],
           potential = set(0.); // Vectorized variables.

    setAccelerationToZero();

    for (int i = 0; i < N - 1; i++)
    {
        // Store the position of the particle i and set the acceleration to zero.
        for (int k = 0; k < 3; k++){
            ri[k] = r[k][i];
            ai[k] = 0.;
            aiV[k] = set(0.);
        }

        int j = i + 1;
        int lastValueVectorizable = N - ((N - j) % 4); // We need to divide the loop into vectorizable values and non-vectorizable ones.

        for (; j < lastValueVectorizable; j += 4) {

            // Difference in each coordinate between particle i and j.
            for (int k = 0; k < 3; k++){
                rijV[k] = sub(set(ri[k]), load(&r[k][j])); // ri[k] - r[k][j]
                rijVsqd[k] = mul(rijV[k],rijV[k]); // rijV[k] * rijV[k]
            }

            // Squared of the distance between particle i and j.
            vector rSqdV = add(
                add(rijVsqd[0], rijVsqd[1]),
                rijVsqd[2]
            ); // rijVsqd[0] + rijVsqd[1] + rijVsqd[2]

            // Forces applied to particle i and j.
            vector fv = lennardJonesForceVector(rSqdV);

            for (int k = 0; k < 3; k++) {
                rijV[k] = mul(rijV[k], fv); // rijV[k] * fv
                aiV[k] = add(aiV[k], rijV[k]); // aiV[k] + rijV[k]

                store(&a[k][j], sub(load(&a[k][j]), rijV[k])); // a[k][j] = a[k][j] - rijV[k]
            }

            potential = add(potential, potentialEnergyVector(rSqdV)); // Update potential energy.
        }

        // Now we handle the values that didn't make into the vectorization.
        // The same logic applied above, works here.
        for (; j < N; j++) {

            for (int k = 0; k < 3; k++) rij[k] = ri[k] - r[k][j];
        
            rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            f = lennardJonesForce(rSqd);

            for (int k = 0; k < 3; k++){
                rij[k] *= f;
                ai[k] += rij[k];
                a[k][j] -= rij[k];
            }

            pot += potentialEnergy(rSqd);
        }
        for (int k = 0; k < 3; k++) {
            a[k][i] += sumVector(aiV[k]);
            a[k][i] += ai[k];
        }
    }

    return 2 * epsilon_times_4 * (sumVector(potential) + pot);
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

    vector rijV[3], rijVsqd[3], aiV[3];

    setAccelerationToZero();

    for (int i = 0; i < N - 1; i++)
    {
        // Store the position of the particle i and set the acceleration to zero.
        for (int k = 0; k < 3; k++) {
            ri[k] = r[k][i];
            ai[k] = 0.;
            aiV[k] = set(0.);
        }

        int j = i + 1;
        int lastValueVectorizable = N - ( (N - j) % 4);

        // Handling vectorizable values.
        for (; j < lastValueVectorizable; j += 4) {

            // Difference in each coordinate between particle i and j.
            for (int k = 0; k < 3; k++){
                rijV[k] = sub(set(ri[k]), load(&r[k][j])); // ri[k] - r[k][j]
                rijVsqd[k] = mul(rijV[k], rijV[k]); // rijV[k] * rijV[k]
            }

            // Squared of the distance between particle i and j.
            vector rSqdV = add(
                add(rijVsqd[0], rijVsqd[1]),
                rijVsqd[2]
            ); // rijVsqd[0] + rijVsqd[1] + rijVsqd[2]

            // Forces applied to particle i and j.
            vector fv = lennardJonesForceVector(rSqdV);

            for (int k = 0; k < 3; k++) {
                rijV[k] = mul(rijV[k],fv); // rijV[k] * fv
                aiV[k] = add(aiV[k],rijV[k]); // aiV[k] * rijV[k]
                store(&a[k][j],sub(load(&a[k][j]),rijV[k])); // a[k][j] = a[k][j] - rijV[k]
            }
        }

        // Handling non vectorizable values.
        for (; j < N; j++) {
            for (int k = 0; k < 3; k++) rij[k] = ri[k] - r[k][j];
        
            rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
            f = lennardJonesForce(rSqd);

            for (int k = 0; k < 3; k++){
                rij[k] *= f;
                ai[k] += rij[k];
                a[k][j] -= rij[k];
            }
        }

        // Sum + Reduce of all accelerations.
        for (int k = 0; k < 3; k++) {
            a[k][i] += sumVector(aiV[k]);
            a[k][i] += ai[k];
        }
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
    
    int i, j;
    double psum = 0.;

    /*
    Compute accelerations from forces at current position 
    This call was removed (commented) for pedagogical reasons computeAccelerations();
    Update positions and velocity with current velocity and acceleration.
    */
    vector dtV = set(dt);
    vector dtV05 = set(dt*0.5);
    for (j = 0; j < 3; j++){
        int i = 0;
        for (; i < N-(N%4); i+=4) {
            vector vji = load(&v[j][i]);
            vector aji = load(&a[j][i]);
            vector rji = load(&r[j][i]);
            store(&r[j][i],add(rji,mul(dtV,add(vji,mul(dtV05,aji)))));
            store(&v[j][i],add(vji,mul(dtV05,aji)));
        }
        for (;i<N;i++){
            r[j][i] += v[j][i] * dt + 0.5 * a[j][i] * dt * dt;
            v[j][i] += 0.5 * a[j][i] * dt;
        }
    }
    //  Update accellerations from updated positions
    *potentialEnergy=computeAccelerationsAndPotentialVector();
    //  Update velocity with updated acceleration
    for (j=0; j<3; j++){
        int i=0;
        for (; i<N-(N%4); i+=4){
            vector vji = load(&v[j][i]);
            vector aji = load(&a[j][i]);
            store(&v[j][i],add(vji,mul(dtV05,aji)));
        }
        for (; i<N; i++){
            v[j][i] += 0.5*a[j][i]*dt;
        }
        
        // Elastic walls.
        for (i = 0; i < N; i++) {
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
