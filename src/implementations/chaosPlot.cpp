/*
 * chaosPlot.cpp
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "..\headers\chaosPlot.h"
#include "..\headers\planetsAndProbe.h"
#include "..\headers\utils.h"

// Main function
// Outputs lyapunov exponents to a csv for plotting a 2d Chaos Map of the Earth Sun Moon JW Telescope system
void chaosPlot(double dt, int totalSteps, int resolution[2], int NUMPOINTS, double RADIUS)
{
	int probeNum = resolution[0] * resolution[1];
	// x and z positions for each probe
	double *positions = (double *)malloc(2 * probeNum * sizeof(double));
	// Get circle perturbations for probes
	double **perturbations = (double **)malloc(NUMPOINTS * sizeof(double *)); for (int i = 0; i < NUMPOINTS; ++i) { perturbations[i] = (double *)malloc(2 * sizeof(double)); }
	returnPerturbationsCircle(perturbations, NUMPOINTS);
	// Get positions for all probes
	partitionPositions(positions, resolution);
	// Calculate lyapunovs: Either multiProbeChaos or Orbital Separation
	double *lyapunovs = multiProbeChaos(dt, totalSteps, positions, probeNum, perturbations, NUMPOINTS, RADIUS);

	// TODO: I NEED TO BE RENORMALIZING, THERE ARE MANY OVERFLOWS
	// double *lyapunovs = orbitalSeparationMultiProbe(dt, totalSteps, positions, probeNum, perturbations, NUMPOINTS, RADIUS);
	// Write to file
	multiLyapunovToFile(lyapunovs, probeNum);
	// Free stuff
	// Free lyapunovs array
	free(lyapunovs); free(positions); for (int i = 0; i < NUMPOINTS; ++i) { free(perturbations[i]); } free(perturbations);
}

// Function implementations utilized by main function chaosPlot

// Returns an array of x, z probe positions distributed across L1 to L2 for all probes
void partitionPositions(double *positions, int resolution[2])
{
	// L1 is -1.5 mil km, L2 is 1.5 mil km. L1 to L2 to 3 mil km
	double xIncr = 3000000000.0 / (double)(resolution[0] - 1); double zIncr = 3000000000.0 / (double)(resolution[1] - 1);
	for (int i = 0; i < resolution[0] * resolution[1]; ++i)
	{
		positions[2 * i] = -1500000000.0 + xIncr * (i % resolution[0]);
		positions[2 * i + 1] = -1500000000.0 + zIncr * (i / resolution[0]);
	}
}

// Submain functions

// Extracts and returns array of min and max lyapunov exponents for a number of probes at various positions. Renormalizes at each timestep
// Orbital Separation method for 2d chaos plot

// FIX THIS FUNCTION: BROKEN. IMPLEMENT PROBES FOR GRAM SCHMIDT
double *orbitalSeparationMultiProbe(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int NUMPOINTS, double RADIUS)
{
	// Initialize Stuff
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe *probes = (Probe *)malloc(probeNum * sizeof(Probe));
	// Store data: Min and Max Lyapunov exponents for each probe
	double *lyapunovs = (double *)malloc(2 * probeNum * sizeof(double));
	// Sum of log of relative separation every timestep for each probes perturbations: ORBITAL SEPARATION METHOD FOR CALCULATING LYAPUNOVS
	// Max lyapunov will be extracted from this at the end
	double **LfSums = (double **)malloc(probeNum * sizeof(double *));
	// Paralellogram: THESE NEED TO BE PROBES
	double **orthoNormal = (double **)malloc(probeNum * sizeof(double *));
	// Min lyapunov will be extracted from this at the end
	double *areaSums = (double *)malloc(probeNum * sizeof(double));
	double initialArea = RADIUS * RADIUS;
	// Initialize perturbations for every probe
	for (int i = 0; i < probeNum; i++)
	{
		// Allocate mem to LfSums
		LfSums[i] = (double *)malloc(NUMPOINTS * sizeof(double));
		// Allocate mem for orthonormal vectors: [x1, z1, x2, z2]
		orthoNormal[i] = (double *)malloc(4 * sizeof(double));
		// Initialize Probe and Perturbed Probes for each. Set LF initial val to 0 as it is a running sum
		// Initialize orthoNormal vectors, areaSums[i] = 0.0
		probes[i].x = probePositions[2 * i]; probes[i].y = 0.0; probes[i].z = probePositions[2 * i + 1]; probes[i].vx = 0.0; probes[i].vy = 0.0; probes[i].vz = 30173.943; probes[i].perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
		orthoNormal[i][0] = RADIUS; orthoNormal[i][1] = 0; orthoNormal[i][2] = 0; orthoNormal[i][3] = RADIUS; areaSums[i] = 0;
		for (int j = 0; j < NUMPOINTS; ++j) { probes[i].perturbed[j].x = probes[i].x + perturbations[j][0]; probes[i].perturbed[j].y = 0.0; probes[i].perturbed[j].z = probes[i].z + perturbations[j][1]; copyProbeToPerturbed(&probes[i].perturbed[j], probes[i]); LfSums[i][j] = 0.0; }
	}
	// Loop
	int c = 0;
	while (c < totalSteps)
	{
		multiChaosRK4(&e, &s, probes, probeNum, dt, NUMPOINTS);

		// Calculate Lf for timestep: and add to LfSums for all probes. Lfk = ln(dk / e) between perturbation and probe. d is distance at time k, e is initial distance
		// Move perturbation back to initial distance but preserve direction: yk* = yk + e * (ck / dk) where ck = (yk* - yk) (renormalize)
		// Calculate difference in area between evolved parallelogram and initial
		// DEBUGGING: maybe not broken?
//		std::setprecision(5);
		for (int i = 0; i < probeNum; ++i)
		{
			// DEBUGGING
//			std::cout << "Probe " << i << " perturbation: dk / e, log(dk / e)" << std::endl;
			for (int j = 0; j < NUMPOINTS; ++j)
			{
				double ck[2] = {probes[i].perturbed[j].x - probes[i].x, probes[i].perturbed[j].z - probes[i].z};
				double dist = sqrt((ck[0] * ck[0]) + (ck[1] * ck[1]));

				// DEBUGGING
//				std::cout << (double)(dist / RADIUS) << " " << log(dist / RADIUS) << std::endl;
				LfSums[i][j] += log(dist / RADIUS);
				probes[i].perturbed[j].x = probes[i].x + RADIUS * (ck[0] / dist); probes[i].perturbed[j].z = probes[i].z + RADIUS * (ck[1] / dist);

				// Gram Schmidt:
				//gramSchmidtRenormalize(orthoNormal, areaSums, probeNum, initialArea);
				//gramSchmidtRenormalize(Probe *probes, double *areaSums, int probeNum, double initialArea, int NUMPOINTS)
			}
//			std::cout << std::endl;
		}
		c++;
	}

	// Extract lyapunovs:
	extractOrbitalMultiLyapunov(probes, probeNum, lyapunovs, LfSums, areaSums, dt * totalSteps, NUMPOINTS);

	// Free stuff
	for (int i = 0; i < probeNum; ++i) { free(probes[i].perturbed); free(LfSums[i]); free(orthoNormal[i]); } free(probes); free(LfSums); free(orthoNormal); free(areaSums);

	return lyapunovs;

}

// Extracts and returns array of min and max lyapunov exponents for a number of probes at various positions
// NO RENORMALIZATION
double *multiProbeChaos(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int NUMPOINTS, double RADIUS)
{
	// Initialize Stuff
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe *probes = (Probe *)malloc(probeNum * sizeof(Probe));

	// Store data: Min and Max Lyapunov exponents for each probe
	double *lyapunovs = (double *)malloc(2 * probeNum * sizeof(double));
	double initialArea = RADIUS * RADIUS;

	// Initialize perturbations for every probe
	for (int i = 0; i < probeNum; i++)
	{
		// Initialize Probe and Perturbed Probes for each
		probes[i].x = probePositions[2 * i]; probes[i].y = 0.0; probes[i].z = probePositions[2 * i + 1]; probes[i].vx = 0.0; probes[i].vy = 0.0; probes[i].vz = 30173.943; probes[i].perturbed = (Probe *)malloc((2 + NUMPOINTS) * sizeof(Probe));
		for (int j = 0; j < NUMPOINTS; ++j) { probes[i].perturbed[j].x = probes[i].x + perturbations[j][0]; probes[i].perturbed[j].y = 0.0; probes[i].perturbed[j].z = probes[i].z + perturbations[j][1]; copyProbeToPerturbed(&probes[i].perturbed[j], probes[i]); }
		// Final two perturbed probes are the orthonormal vectors for gram schmidt
		// V1: <1, 0, 0>, V2: <0, 0, 1>
		copyProbe(&probes[i].perturbed[NUMPOINTS], probes[i]); probes[i].perturbed[NUMPOINTS].x += RADIUS; copyProbe(&probes[i].perturbed[NUMPOINTS + 1], probes[i]); probes[i].perturbed[NUMPOINTS + 1].z += RADIUS;

	}
	// Count
	int c = 0;
	// Loop
	while (c < totalSteps)
	{
		// Either Euler or RK4
		// multiChaosEuler(&e, &s, probes, probeNum, dt, NUMPOINTS);
		// RK4: PROBABLY BROKEN
		multiChaosRK4(&e, &s, probes, probeNum, dt, NUMPOINTS);
		c++;
	}
	// Extract lyapunovs
	extractMultiLyapunov(probes, probeNum, initialArea, lyapunovs, dt * totalSteps, NUMPOINTS, RADIUS);

	// Free stuff
	for (int i = 0; i < probeNum; ++i) { free(probes[i].perturbed); } free(probes);

	return lyapunovs;
}

// Function Implementations utilized by submain functions

// Classic Runge Kutta for multiple probes with perturbations
void multiChaosRK4(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int NUMPOINTS)
{
	// Initialize Planet arrays
	Planet k2Planets[2], k3Planets[2], k4Planets[2];
	// Intialize Probe arrays
	Probe *k2Probes = (Probe *)malloc(probeNum * sizeof(Probe)); Probe *k3Probes = (Probe *)malloc(probeNum * sizeof(Probe)); Probe *k4Probes = (Probe *)malloc(probeNum * sizeof(Probe));
	// Copy data into planet arrays for all k
	copyPlanet(&k2Planets[0], *e); copyPlanet(&k2Planets[1], *s); copyPlanet(&k3Planets[0], *e); copyPlanet(&k3Planets[1], *s); copyPlanet(&k4Planets[0], *e); copyPlanet(&k3Planets[1], *s);
	// Initialize kAccels
	double **k1Accels = (double **)malloc((1 + probeNum) * sizeof(double *)); double **k2Accels = (double **)malloc((1 + probeNum) * sizeof(double *)); double **k3Accels = (double **)malloc((1 + probeNum) * sizeof(double *)); double **k4Accels = (double **)malloc((1 + probeNum) * sizeof(double *));
	// Earth kAccels: x and z. Only need two doubles
	k1Accels[0] = (double *)malloc(2 * sizeof(double)); k2Accels[0] = (double *)malloc(2 * sizeof(double)); k3Accels[0] = (double *)malloc(2 * sizeof(double)); k4Accels[0] = (double *)malloc(2 * sizeof(double));
	// Copy data into probeCopies and their perturbations for all k
	for (int i = 0; i < probeNum; ++i)
	{
		// Allocate memory for Probe kAccels. Indexes 0 and 1 are for probe x and z, the rest are perturbations x and z.
		k1Accels[i + 1] = (double *)malloc((2 + 2 * (2 + NUMPOINTS)) * sizeof(double)); k2Accels[i + 1] = (double *)malloc((2 + 2 * (2 + NUMPOINTS)) * sizeof(double)); k3Accels[i + 1] = (double *)malloc((2 + 2 * (2 + NUMPOINTS)) * sizeof(double)); k4Accels[i + 1] = (double *)malloc((2 + 2 * (2 + NUMPOINTS)) * sizeof(double));
		// Copy for Probe copies and initialize perturbed arrays for each probe
		copyProbe(&k2Probes[i], probes[i]); k2Probes[i].perturbed = (Probe *)malloc((2 + NUMPOINTS) * sizeof(Probe)); copyProbe(&k3Probes[i], probes[i]); k3Probes[i].perturbed = (Probe *)malloc((2 + NUMPOINTS) * sizeof(Probe)); copyProbe(&k4Probes[i], probes[i]); k4Probes[i].perturbed = (Probe *)malloc((2 + NUMPOINTS) * sizeof(Probe));
		// Copy for Perturbation copies
		for (int j = 0; j < NUMPOINTS + 2; ++j) { copyProbe(&k2Probes[i].perturbed[j], probes[i].perturbed[j]); copyProbe(&k3Probes[i].perturbed[j], probes[i].perturbed[j]); copyProbe(&k4Probes[i].perturbed[j], probes[i].perturbed[j]); }
	}

	// k1 step
	multiChaosCalcAccel(*e, *s, probes, probeNum, k1Accels, NUMPOINTS);
	// k2 Step
	rk4StepMultiChaos(&k2Planets[0], k2Probes, *e, probes, probeNum, k1Accels, dt / 2, NUMPOINTS);
	multiChaosCalcAccel(k2Planets[0], k2Planets[1], k2Probes, probeNum, k2Accels, NUMPOINTS);
	// k3 Step
	rk4StepMultiChaos(&k3Planets[0], k3Probes, k2Planets[0], k2Probes, probeNum, k2Accels, dt / 2, NUMPOINTS);
	multiChaosCalcAccel(k3Planets[0], k3Planets[1], k3Probes, probeNum, k3Accels, NUMPOINTS);
	// k4 Step
	rk4StepMultiChaos(&k4Planets[0], k4Probes, k3Planets[0], k3Probes, probeNum, k3Accels, dt, NUMPOINTS);
	multiChaosCalcAccel(k4Planets[0], k4Planets[1], k4Probes, probeNum, k4Accels, NUMPOINTS);
	// rk4 Step
	// Move objects
	e->vx += (k1Accels[0][0] + 2 * k2Accels[0][0] + 2 * k3Accels[0][0] + k4Accels[0][0]) * dt / 6.0; e->vz += (k1Accels[0][1] + 2 * k2Accels[0][1] + 2 * k3Accels[0][1] + k4Accels[0][1]) * dt / 6.0; e->x += (e->vx + 2 * k2Planets[0].vx + 2 * k3Planets[0].vx + k4Planets[0].vx) * dt / 6.0;  e->z += (e->vz + 2 * k2Planets[0].vz + 2 * k3Planets[0].vz + k4Planets[0].vz) * dt / 6.0;
	// Probes and perturbations
	for (int i = 0; i < probeNum; ++i)
	{
		probes[i].vx += (k1Accels[i + 1][0] + 2 * k2Accels[i + 1][0] + 2 * k3Accels[i + 1][0] + k4Accels[i + 1][0]) * dt / 6.0; probes[i].vz += (k1Accels[i + 1][1] + 2 * k2Accels[i + 1][1] + 2 * k3Accels[i + 1][1] + k4Accels[i + 1][1]) * dt / 6.0; probes[i].x += (probes[i].vx + 2 * k2Probes[i].vx + 2 * k3Probes[i].vx + k4Probes[i].vx) * dt / 6.0; probes[i].z += (probes[i].vz + 2 * k2Probes[i].vz + 2 * k3Probes[i].vz + k4Probes[i].vz) * dt / 6.0;
		// Perturbations
		for (int j = 0; j < NUMPOINTS + 2; ++j)
		{
			probes[i].perturbed[j].vx += (k1Accels[i + 1][2 * j + 2] + 2 * k2Accels[i + 1][2 * j + 2] + 2 * k3Accels[i + 1][2 * j + 2] + k4Accels[i + 1][2 * j + 2]) * dt / 6.0; probes[i].perturbed[j].vz += (k1Accels[i + 1][2 * j + 3] + 2 * k2Accels[i + 1][2 * j + 3] + 2 * k3Accels[i + 1][2 * j + 3] + k4Accels[i + 1][2 * j + 3]) * dt / 6.0; probes[i].perturbed[j].x += (probes[i].perturbed[j].vx + 2 * k2Probes[i].perturbed[j].vx + 2 * k3Probes[i].perturbed[j].vx + k4Probes[i].perturbed[j].vx) * dt / 6.0; probes[i].perturbed[j].z += (probes[i].perturbed[j].vz + 2 * k2Probes[i].perturbed[j].vz + 2 * k3Probes[i].perturbed[j].vz + k4Probes[i].perturbed[j].vz) * dt / 6.0;
		}
	}

	// Free stuff
	// Probes and Perturbations and Accels
	free(k1Accels[0]); free(k2Accels[0]); free(k3Accels[0]); free(k4Accels[0]); for (int i = 0; i < probeNum; ++i) { free(k2Probes[i].perturbed); free(k3Probes[i].perturbed); free(k4Probes[i].perturbed); free(k1Accels[i + 1]); free(k2Accels[i + 1]); free(k3Accels[i + 1]); free(k4Accels[i + 1]); } free(k1Accels); free(k2Accels); free(k3Accels); free(k4Accels); free(k2Probes); free(k3Probes); free(k4Probes);
}

// Euler method for multiple probes and perturbations
void multiChaosEuler(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int NUMPOINTS)
{
	double **accels = (double **)malloc((1 + probeNum) * sizeof(double *)); accels[0] = (double *)malloc(2 * sizeof(double)); for (int i = 1; i < 1 + probeNum; ++i) { accels[i] = (double *)malloc((2 + 2 * (NUMPOINTS + 2)) * sizeof(double)); }
	multiChaosCalcAccel(*e, *s, probes, probeNum, accels, NUMPOINTS);

	// DEBUGGING
//	if (row == 4 && loopCount == 17)
//	{
//		std::setprecision(5);
//		std::cout << "MULTI CHAOS EULER" << std::endl;
//		std::cout << "Loop " << loopCount << std::endl;
//		std::cout << "Printing position of probe" << std::endl;
//		std::cout << "x: " << probes[2].x << " y: " << probes[2].y << "z: " << probes[2].z << std::endl << std::endl;
//		std::cout << "Printing velocities of probe" << std::endl;
//		std::cout << "vx: " << probes[2].vx << " vy: " << probes[2].vy << " vz: " << probes[2].vz << std::endl << std::endl;
//		std::cout << "Printing accelerations on probe" << std::endl;
//		std::cout << "ax: " << accels[3][0] << " ay: " << 0.0 << "az: " << accels[3][1] << std::endl << std::endl;
//		std::cout << "Printing accelerations on earth" << std::endl;
//		std::cout << "ax: " << accels[0][0] << " ay: " << 0.0 << "az: " << accels[0][1] << std::endl << std::endl;
//	}
	e->vx += accels[0][0] * dt; e->vz += accels[0][1] * dt;
	e->x += e->vx * dt; e->z += e->vz * dt;
	for (int i = 0; i < probeNum; ++i)
	{
		probes[i].vx += accels[i + 1][0] * dt; probes[i].vz += accels[i + 1][1] * dt;
		probes[i].x += probes[i].vx * dt; probes[i].z += probes[i].vz * dt;
		for (int j = 0; j < NUMPOINTS + 2; ++j)
		{
			probes[i].perturbed[j].vx += accels[i + 1][2 * j + 2] * dt; probes[i].perturbed[j].vz += accels[i + 1][2 * j + 3] * dt;
			probes[i].perturbed[j].x += probes[i].perturbed[j].vx * dt; probes[i].perturbed[j].z += probes[i].perturbed[j].vz * dt;
		}
	}
}


// Calculates accelerations. Handles multiple probes with perturbations
void multiChaosCalcAccel(Planet e, Planet s, Probe *probes, int probeNum, double **accels, int NUMPOINTS)
{
	double G = 6.6743e-11;

	// Sun on Earth: x and z
	double r21[2] = {e.x - s.x, e.z - s.z}; double magr21 = sqrt((r21[0] * r21[0]) + (r21[1] * r21[1])); double r21Unit[2] = {r21[0] / magr21, r21[1] / magr21}; double gravity = (-1 * G * s.m) / (magr21 * magr21);
	accels[0][0] = gravity * r21Unit[0]; accels[0][1] = gravity * r21Unit[1];

	// Calculate for probes and perturbations
	for (int i = 0; i < probeNum; ++i)
	{
		// Sun on James Webb telescope
		double r31[2] = {probes[i].x - s.x, probes[i].z - s.z}; double magr31 = sqrt((r31[0] * r31[0]) + (r31[1] * r31[1])); double r31Unit[2] = {r31[0] / magr31, r31[1] / magr31}; double gravity31 = (-1 * G * s.m) / (magr31 * magr31);
		// Earth on JW
		double r32[2] = {probes[i].x - e.x, probes[i].z - e.z};
		double magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1]));
		double r32Unit[2] = {r32[0] / magr32, r32[1] / magr32};
		double gravity32 = (-1 * G * e.m) / (magr32 * magr32);
		// Net acceleration
		accels[i + 1][0] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]); accels[i + 1][1] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);

		// DEBUGGING
//		if (row == 4 && i == 2)
//		{
//			std::cout << "Multi chaos calc force" << std::endl;
//			std::cout << "Sun on earth gravity" << std::endl;
//			std::cout << gravity << "Vectors x: " << accels[0][0] << " y: " << 0.0 << " z: " << accels[0][1] << std::endl;
//			std::cout << "Sun on probe gravity" << std::endl;
//			std::cout << gravity31 << "Vectors x: " << gravity31 * r31Unit[0] << " y: " << 0.0 << " z: " << gravity31 * r31Unit[1] << std::endl;
//			std::cout << "Earth on probe gravity" << std::endl;
//			std::cout << gravity32 << "Vectors x: " << gravity32 * r32Unit[0] << " y: " << 0.0 << " z: " << gravity32 * r32Unit[1] << std::endl;
//			std::cout << "Earth probe vectors" << std::endl;
//			std::cout << "Vector r32: x:" << r32[0] << " y: " << 0.0 << " z: " << r32[1] << std::endl;
//			std::cout << "Vector Magnitude: " << magr32 << std::endl;
//			std::cout << "Unit vector x: " << r32Unit[0] << " y: " << 0.0 << " z: " << r32Unit[1] << std::endl;
//			std::cout << "Net gravity on probe" << std::endl;
//			std::cout << "Vectors x: " << accels[i + 1][0] << " y: " << 0.0 << " z: " << accels[i + 1][1] << std::endl << std::endl;
//		}

		// Calc for all perturbations of probe as well as two gram schmidt probes
		for (int j = 0; j < NUMPOINTS + 2; ++j)
		{
			// Sun on perturbation
			r31[0] = probes[i].perturbed[j].x - s.x; r31[1] = probes[i].perturbed[j].z - s.z; magr31 = sqrt((r31[0] * r31[0]) + (r31[1] * r31[1])); r31Unit[0] = r31[0] / magr31; r31Unit[1] = r31[1] / magr31; gravity31 = (-1 * G * s.m) / (magr31 * magr31);
			// Earth on Perturbation
			r32[0] = probes[i].perturbed[j].x - e.x; r32[1] = probes[i].perturbed[j].z - e.z; magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1])); r32Unit[0] = r32[0] / magr32; r32Unit[1] = r32[1] / magr32; gravity32 = (-1 * G * e.m) / (magr32 * magr32);
			// Net
			accels[i + 1][2 * j + 2] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]); accels[i + 1][2 * j + 3] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);
		}
	}
}

// Moves planets and probes one step forward. Helper function to multiChaosRK4
void rk4StepMultiChaos(Planet *eCopy, Probe *probeCopies, Planet e, Probe *probes, int probeNum, double **accels, double dt, int NUMPOINTS)
{
	// First index is earth accels
	eCopy->vx += accels[0][0] * dt; eCopy->vz += accels[0][1] * dt; eCopy->x += e.vx * dt; eCopy->z += e.vz * dt;
	// Step for probes and perturbations
	for (int i = 0; i < probeNum; ++i)
	{
		// Probe step
		probeCopies[i].vx += accels[i + 1][0] * dt; probeCopies[i].vz += accels[i + 1][1] * dt; probeCopies[i].x += probes[i].vx * dt; probeCopies[i].z += probes[i].vz * dt;
		// Perturbations step
		for (int j = 0; j < NUMPOINTS + 2; ++j) { probeCopies[i].perturbed[j].vx += accels[i + 1][2 * j + 2] * dt; probeCopies[i].perturbed[j].vz += accels[i + 1][2 * j + 3] * dt; probeCopies[i].perturbed[j].x += probes[i].perturbed[j].vx * dt; probeCopies[i].perturbed[j].z += probes[i].perturbed[j].vz * dt; }
	}
}

// Extracts lyapunov exponents from multiple probes and deposits them into an array
void extractMultiLyapunov(Probe *probes, int probeNum, double initialArea, double *lyapunovs, double t, int NUMPOINTS, double RADIUS)
{
	// DEBUGGING
	std::cout << "MULTI CHAOS LYAPUNOV" << std::endl;
	std::cout << "Printing final position of probe of interest" << std::endl;
	std::cout << "X pos: " << probes[22].x << " Z pos: " << probes[22].z << std::endl;
	std::cout << std::endl;

	// Loop over every probe
	for (int i = 0; i < probeNum; ++i)
	{
		// Maximum lyapunov calculation
		// Calculate distance from every perturbation to center probe
		double v[2] = {probes[i].perturbed[0].x - probes[i].x, probes[i].perturbed[0].z - probes[i].z};
		double dist = sqrt((v[0] * v[0]) + (v[1] * v[1]));

		// DEBUGGING
//		if (i == 22)
//		{
//			std::cout << "Probe perturbation 0 is at pos: x: " << probes[i].perturbed[0].x << " z: " << probes[i].perturbed[0].z << std::endl;
//		}
		// Search for max and min distance and index
		double max = dist;
		int ind = 0;
		for (int j = 1; j < NUMPOINTS; ++j)
		{
			v[0] = probes[i].perturbed[j].x - probes[i].x; v[1] = probes[i].perturbed[j].z - probes[i].z;
			dist = sqrt((v[0] * v[0]) + (v[1] * v[1]));
			if (dist > max) { max = dist; ind = j; }
			// DEBUGGING
			if (i == 22)
			{
				std::cout << "Probe perturbation " << j << " is at pos: x: " << probes[i].perturbed[j].x << " z: " << probes[i].perturbed[j].z << std::endl;
				std::cout << "Distance to probe is: " << dist << std::endl;
			}
		}

		// DEBUGGING
		if (i == 22)
		{
			std::cout << "Printing gram schmidt perturbations" << std::endl;
			std::cout << "V1: pos: x: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS].x << " z: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS].z << std::endl;
			std::cout << "V2: pos: x: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS + 1].x << " z: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS + 1].z << std::endl;
			std::cout << "Printing max distance perturbation vector: x: " << std::setprecision(20) << probes[i].perturbed[ind].x << " z: " << std::setprecision(20) << probes[i].perturbed[ind].z << std::endl;
			double v1[2] = {probes[i].perturbed[NUMPOINTS].x - probes[i].x, probes[i].perturbed[NUMPOINTS].z - probes[i].z};
			double v2[2] = {probes[i].perturbed[NUMPOINTS + 1].x - probes[i].x, probes[i].perturbed[NUMPOINTS + 1].z - probes[i].z};
			double maxv[2] = {probes[i].perturbed[ind].x - probes[i].x, probes[i].perturbed[ind].z - probes[i].z};
			double magv1 = sqrt((v1[0] * v1[0]) + (v1[1] * v1[1]));
			double magv2 = sqrt((v2[0] * v2[0]) + (v2[1] * v2[1]));
			double magMaxv = sqrt((maxv[0] * maxv[0]) + (maxv[1] * maxv[1]));
			std::cout << "Max dist found in perturbation is: " << std::setprecision(10) << max << std::endl;
			std::cout << "Magnitude of max vector is: " << std::setprecision(10) << magMaxv << std::endl;
			// Make them unit vectors for calculating angle to avoid overflow
			double u1[2] = {v1[0] / magv1, v1[1] / magv1};
			double u2[2] = {v2[0] / magv2, v2[1] / magv2};
			double umax[2] = {maxv[0] / magMaxv, maxv[1] / magMaxv};
			double u1Dotumax = (u1[0] * umax[0]) + (u1[1] * umax[1]);
			double u2Dotumax = (u2[0] * umax[0]) + (u2[1] * umax[1]);
			double a1 = acos(u1Dotumax); double a2 = acos(u2Dotumax);
			std::cout << "Printing angle between v1 and max vector" << std::endl;
			std::cout << std::setprecision(10) << a1 << std::endl;
			std::cout << "Printing angle between v2 and max vector" << std::endl;
			std::cout << std::setprecision(10) << a2 << std::endl;
		}
		// Calculate and store max lyapunov
		lyapunovs[2 * i + 1] = (1.0 / t) * log(max / RADIUS);
		// Calculate and store min lyapunov
		double logAreaDif = gramSchmidt(probes, initialArea, i, NUMPOINTS, ind);
		lyapunovs[2 * i] = (logAreaDif / t) - lyapunovs[2 * i + 1];
	}

	// DEBUGGING
	std::setprecision(5);
	for (int i = 0; i < probeNum; ++i)
	{
		std::cout << "Probe " << i << std::endl;
		std::cout << "Final pos: x: " << probes[i].x << " z: " << probes[i].z << std::endl;
		std::cout << "Lyapunovs: Lyapunov min: " << lyapunovs[2 * i] << " Lyapunov max: " << lyapunovs[2 * i + 1] << std::endl;
	}
}
// Gram Schmidt: Needed for calculation of minimum lyapunov. Returns difference in area between evolved parallelogram and original rectange. Calculates for single probe
double gramSchmidt(Probe *probes, double initialArea, int i, int NUMPOINTS, int maxInd)
{
	// Find angle between vectors: Final two perturbation probes
	double v1[2] = {probes[i].perturbed[NUMPOINTS].x - probes[i].x, probes[i].perturbed[NUMPOINTS].z - probes[i].z};
	double v2[2] = {probes[i].perturbed[NUMPOINTS + 1].x - probes[i].x, probes[i].perturbed[NUMPOINTS + 1].z - probes[i].z};
	double magv1 = sqrt((v1[0] * v1[0]) + (v1[1] * v1[1]));
	double magv2 = sqrt((v2[0] * v2[0]) + (v2[1] * v2[1]));
	// Make them unit vectors for calculating angle to avoid overflow
	double u1[2] = {v1[0] / magv1, v1[1] / magv1};
	double u2[2] = {v2[0] / magv2, v2[1] / magv2};
	double u1Dotu2 = (u1[0] * u2[0]) + (u1[1] * u2[1]);
	double a = acos(u1Dotu2);
	if (i == 22)
	{
		std::cout << std::endl << "GRAM SCHMIDT" << std::endl;
		std::cout << "V1: x: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS].x << " z: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS].z << std::endl;
		std::cout << "V2: x: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS + 1].x << " z: " << std::setprecision(20) << probes[i].perturbed[NUMPOINTS + 1].z << std::endl;
		std::cout << "V1mag: " << std::setprecision(20) << magv1 << std::endl;
		std::cout << "V2mag: " << std::setprecision(20) << magv2 << std::endl;
		std::cout << "v1Dotv2: " << std::setprecision(20) << u1Dotu2 << std::endl;
		std::cout << "Angle between vectors is: " << a << std::endl;
		double maxv[2] = {probes[i].perturbed[maxInd].x - probes[i].x, probes[i].perturbed[maxInd].z - probes[i].z};
		double magMaxv = sqrt((maxv[0] * maxv[0]) + (maxv[1] * maxv[1]));
		double umax[2] = {maxv[0] / magMaxv, maxv[1] / magMaxv};


		gramSchmidtToFile(u1, u2, umax);
	}
	// Return area
	return log(magv1 * magv2 * sin(a) / initialArea);
}
// Gram Schmidt: Needed for calculation of minimum lyapunov. Calculates difference in area and renormalizes to RADIUS orthogonal vectors that preserve direction
void gramSchmidtRenormalize(Probe *probes, double *areaSums, int probeNum, double initialArea, int NUMPOINTS)
{
	for (int i = 0; i < probeNum; ++i)
	{
		// Find angle between vectors
		double magv1 = sqrt((probes[i].perturbed[NUMPOINTS].x * probes[i].perturbed[NUMPOINTS].x) + (probes[i].perturbed[NUMPOINTS].z * probes[i].perturbed[NUMPOINTS].z));
		double magv2 = sqrt((probes[i].perturbed[NUMPOINTS + 1].x * probes[i].perturbed[NUMPOINTS + 1].x) + (probes[i].perturbed[NUMPOINTS + 1].z * probes[i].perturbed[NUMPOINTS + 1].z));
		double v1Dotv2 = (probes[i].perturbed[NUMPOINTS].x * probes[i].perturbed[NUMPOINTS + 1].x) + (probes[i].perturbed[NUMPOINTS].z * probes[i].perturbed[NUMPOINTS + 1].z);
		double a = acos(v1Dotv2 / (magv1 * magv2));
		// Calculate area
		areaSums[i] += log(magv1 * magv2 * sin(a) / initialArea);
		// Gram schmidt: u1 = v1 / magv1: y2 = v2 - (v2 dot u1)*u1: u2 = y2 / magy2
		double u1[2] = {probes[i].perturbed[NUMPOINTS].x / magv1, probes[i].perturbed[NUMPOINTS].z / magv1};
		double v2Dotu1 = (probes[i].perturbed[NUMPOINTS + 1].x * u1[0]) + (probes[i].perturbed[NUMPOINTS + 1].z * u1[1]);
		double y2[2] = {probes[i].perturbed[NUMPOINTS + 1].x - (v2Dotu1 * u1[0]), probes[i].perturbed[NUMPOINTS + 1].z - (v2Dotu1 * u1[1])};
		double magy2 = sqrt((y2[0] * y2[0]) + (y2[1] * y2[1]));
		probes[i].perturbed[NUMPOINTS].x = probes[i].x + u1[0]; probes[i].perturbed[NUMPOINTS].z = probes[i].z + u1[1]; probes[i].perturbed[NUMPOINTS + 1].x = (y2[0] / magy2) + probes[i].x; probes[i].perturbed[NUMPOINTS + 1].z = (y2[1] / magy2) + probes[i].z;
	}
}

// Extracts lyapunov exponents for multiple probes using orbit separation method
void extractOrbitalMultiLyapunov(Probe *probes, int probeNum, double *lyapunovs, double **LfSums, double *areaSums, double t, int NUMPOINTS)
{
	for (int i = 0; i < probeNum; ++i)
	{
		// Find min and max Ls
		double max = LfSums[i][0];
		for (int j = 1; j < NUMPOINTS; ++j)
		{
			if (LfSums[i][j] > max) { max = LfSums[i][j]; }
		}
		// Calculate max lyapunov and min lyapunov: From eq: 4.40, 4.41 in Chaos Dynamic Systems
		max /= t;
		// Calculate min lyapunov
		lyapunovs[2 * i] = (areaSums[i] / t) - max; lyapunovs[2 * i + 1] = max;
	}

	// Print results
	std::setprecision(5);
	for (int i = 0; i < probeNum; ++i)
	{
		std::cout << "Probe " << i << " Lyapunov min: " << lyapunovs[2 * i] << " Lyapunov max: " << lyapunovs[2 * i + 1] << std::endl;
	}
}

// Writes to csv file lyapunov exponent mins, maxes, and maxes - mins.
void multiLyapunovToFile(double *data, int dataLen)
{
	std::ofstream file;
	file.open("lyapunovChaosPlot.csv");

	file << "min,max,max-min" << std::endl;
	for (int i = 0; i < dataLen; ++i)
	{
		file << std::to_string(data[2 * i]) << "," << std::to_string(data[2 * i + 1]) << "," << std::to_string(data[2 * i + 1] - data[2 * i]) << std::endl;
	}
	file.close();
}

// Outputs normalized initial and final gram schmidt vectors for specific probe to csv for plotting. Called by gram schmidt function.
void gramSchmidtToFile(double u1f[2], double u2f[2], double umax[2])
{
	// DEBUGGING
	std::cout << "Unit vector v1f: x: " << std::setprecision(10) << u1f[0] << " z: " << std::setprecision(10) << u1f[1] << std::endl;
	std::cout << "Unit vector v2f: x: " << std::setprecision(10) << u2f[0] << " z: " << std::setprecision(10) << u2f[1] << std::endl;
	std::ofstream file1, file2;						// File1 is initial, File2 is final
	file1.open("gramSchmidtInitial.csv"); file2.open("gramSchmidtFinal.csv");
	// Probe 23 of 25 is at L1 lagrange point if chaos plot is 5x5 resolution
	file1 << "x,z" << std::endl; file2 << "x,z" << std::endl;
	// Initial x,z : V1i x,z: V2i x,z
	file1 << std::to_string(1.0) << "," << std::to_string(0.0) << std::endl << std::to_string(0.0) << "," << std::to_string(1.0) << std::endl << std::to_string(umax[0]) << "," << std::to_string(umax[1]) << std::endl;
	// Probe final x,z: V1f x,z: v2f x,z
	file2 << std::to_string(u1f[0]) << "," << std::to_string(u1f[1]) << std::endl << std::to_string(u2f[0]) << "," << std::to_string(u2f[1]) << std::endl << std::to_string(umax[0]) << "," << std::to_string(umax[1]) << std::endl;
	file1.close(); file2.close();
}


