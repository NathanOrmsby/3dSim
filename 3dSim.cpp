//============================================================================
// Name        : 3dSim.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "sphereGenerator.h"


// DEFINED GLOBALS
#define bufSize 9
#define NUMPOINTS 360
#define ITERATIONS 1				// Plotting data per ITERATIONS timesteps
#define RADIUS 500000.0   			// Radius of generated circle or sphere

typedef struct
{
	double x, y, z, vx, vy, vz;
	double m;
} Planet;

typedef struct Probe
{
	double x, y, z, vx, vy, vz;
	Probe *perturbed;
} Probe;

//-------------------------------------------------------------------------------------------
// Declarations
void calcForce(Planet e, Planet s, double **forces);
void chaosCalcForce(Planet e, Planet s, Probe jw, double *accels);
void rk4(Planet *e, Planet *s, Planet *jw, double dt);
void chaosRK4(Planet *e, Planet *s, Probe *jw, double dt);
void multiChaosRK4(Planet *e, Planet *s, Probe *probes, int probeNum, double dt);
void euler(Planet *e, Planet *s, Planet *jw, double dt);
void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt, int loopCount);
void multiChaosEuler(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int loopCount, int row);
void rk4StepSingleChaos(Planet *eCopy, Probe *jwCopy, Planet e, Probe jw, double *accels, double dt);
void rk4StepMultiChaos(Planet *eCopy, Probe *probeCopies, Planet e, Probe *probes, int probeNum, double **accels, double dt);
void copyPlanet(Planet *o, Planet i);
bool orbitsToFile(std::string fileName, double **data, int dataLen);
void perturbationsToFile(double **data, int dataLen, bool offset);
void multiLyapunovToFile(double **data, int dataLen);
void initializePerturbations(Probe *p, bool sphere);
void copyProbeToPerturbed(Probe *o, Probe i);
void copyProbe(Probe *o, Probe i);
void singleProbeChaos(double dt, int totalSteps, bool offset, bool sphere);
double **multiProbeChaos(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int row);
void chaosPlot(double dt, int totalSteps, int resolution[2]);
void extractSingleLyapunov(double **data, int dataLen, double lyapunov[2], bool sphere, double t);
void extractMultiLyapunov(Probe *probes, int probeNum, double **lyapunovs, double t);
void partitionPositions(double *positions, int resolution[2], int row);\
void returnPerturbations(double **perturberations);
void multiLyapunovToArr(double **data, int dataLen, double **lyapunovs, int probeNum);

int main() {

	// Generate unit sphere or circle if not already done
	bool offset = false;
	bool sphere = false;

	if (sphere) { writeSphereToFile(NUMPOINTS, RADIUS); }
	else { writeCircleToFile(NUMPOINTS, RADIUS); }

	// Timestep
	double dt = 500;
	int totalSteps = 18;

	// 2D Chaos Plot
	int resolution[2] = {5, 5};
	// I THINK THIS IS BROKEN
	chaosPlot(dt, totalSteps, resolution);
	// Debugging
	singleProbeChaos(dt, totalSteps, offset, sphere);
	return 0;

	// Initialization
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe jw; jw.x = 0; jw.y = 0; jw.z = 1500000000; jw.vx = 0; jw.vy = 0; jw.vz = 28000; jw.perturbed = nullptr;


	// Debugging
	initializePerturbations(&jw, false);
	return 0;

	// Loop
	int c = 0;
	int stop = 50;

	// fileName
	std::string fname = "data.csv";


	// Storage buffer
	double **data = (double **)malloc(stop * sizeof(double *));
	for (int i = 0; i < stop; ++i)
	{
		data[i] = (double *)malloc(bufSize * sizeof(double));
	}

	// DEBUGGING
//	std::cout << "AT START" << std::endl;
//	std::cout << "Object positions" << std::endl;
//	std::cout << "Earth" << std::endl;
//	std::cout << e.x << " " << e.y << " " << e.z << std::endl;
//	std::cout << "Sun" << std::endl;
//	std::cout << s.x << " " << s.y << " " << s.z << std::endl;
//	std::cout << "James Webb" << std::endl;
//	std::cout << jw.x << " " << jw.y << " " << jw.z << std::endl;

	while (c < stop)
	{
		// Deposit data
		// std::cout << c << std::endl;
		data[c][0] = e.x;
		data[c][1] = e.y;
		data[c][2] = e.z;
		data[c][3] = s.x;
		data[c][4] = s.y;
		data[c][5] = s.z;
		data[c][6] = jw.x;
		data[c][7] = jw.y;
		data[c][8] = jw.z;
		// rk4(&e, &s, &jw, dt);
//		euler(&e, &s, &jw, dt);
		c++;
	}

	// Write data to file
	orbitsToFile(fname, data, stop);

	// DEBUGGING
//	std::cout << "Object positions" << std::endl;
//	std::cout << "Earth" << std::endl;
//	std::cout << e.x << " " << e.y << " " << e.z << std::endl;
//	std::cout << "Sun" << std::endl;
//	std::cout << s.x << " " << s.y << " " << s.z << std::endl;
//	std::cout << "James Webb" << std::endl;
//	std::cout << jw.x << " " << jw.y << " " << jw.z << std::endl;
}

//------------------------------------------------------------------------------------------------------
// Main functions

// Produce 2d chaos plot over the distance L1 to L2 with resolution X x X probes
void chaosPlot(double dt, int totalSteps, int resolution[2])
{
	// Data for storing lyapunov exponents
	double **data = (double **)malloc(10 * resolution[0] * sizeof(double *)); for (int i = 0; i < 10 * resolution[0]; ++i) { data[i] = (double *)malloc(2 * sizeof(double)); }
	int dataLen = 0;
	// Do the plot one row at a time
	double *positions = (double *)malloc(2 * resolution[0] * sizeof(double));				// x and z positions for each probe
	// Get circle perturbations for probes
	double **perturbations = (double **)malloc(NUMPOINTS * sizeof(double *)); for (int i = 0; i < NUMPOINTS; ++i) { perturbations[i] = (double *)malloc(2 * sizeof(double)); }
	returnPerturbations(perturbations);
	// std::cout << "Made it here 1" << std::endl;

	for (int row = 0; row < resolution[1]; ++row)
	{
		// Get positions for probes in row
		// BROKEN
		partitionPositions(positions, resolution, row);

		// DEBUGGING
		std::cout << "Printing positions of probes" << std::endl;
		for (int i = 0; i < resolution[0]; ++i)
		{
			std::cout << positions[2 * i] << " , " << positions[2 * i + 1] << std::endl;
		}
		std::cout << std::endl << std::endl;

		// Do simulation for probes and get lyapunovs
		// BROKEN: MULTI PROBE CHAOS IS BROKEN
		double **lyapunovs = multiProbeChaos(dt, totalSteps, positions, resolution[0], perturbations, row);
		// Dump lyapunovs into data arr, increment data length
		multiLyapunovToArr(data, dataLen, lyapunovs, resolution[0]);
		dataLen += resolution[0];

		// Write to file at certain length, Clear data array
		if (dataLen == 10 * resolution[0]) { multiLyapunovToFile(data, dataLen); dataLen = 0; }

		// Free lyapunovs array
		for (int i = 0; i < resolution[0]; ++i) { free(lyapunovs[i]); } free(lyapunovs);
	}

	// Write residual data to file if there is any
	if (dataLen > 0) { multiLyapunovToFile(data, dataLen); }

	// Free stuff
	for (int i = 0; i < 10 * resolution[0]; ++i) { free(data[i]); } free(data); free(positions); for (int i = 0; i < NUMPOINTS; ++i) { free(perturbations[i]); } free(perturbations);
}

// Extracts and returns array of min and max lyapunov exponents for a number of probes at various positions
double **multiProbeChaos(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int row)
{
	// Initialize Stuff
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe *probes = (Probe *)malloc(probeNum * sizeof(Probe));

	// Store data: Min and Max Lyapunov exponents for each probe
	double **lyapunovs = (double **)malloc(probeNum * sizeof(double *));

	// Initialize perturbations for every probe
	// MAYBE BROKEN
	for (int i = 0; i < probeNum; i++)
	{
		lyapunovs[i] = (double *)malloc(2 * sizeof(double));				// Allocate space for lyapunov exponents
		// Initialize Probe and Perturbed Probes for each
		probes[i].x = probePositions[2 * i]; probes[i].y = 0.0; probes[i].z = probePositions[2 * i + 1]; probes[i].vx = 0.0; probes[i].vy = 0.0; probes[i].vz = 30173.943; probes[i].perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
		for (int j = 0; j < NUMPOINTS; ++j) { probes[i].perturbed[j].x = probes[i].x + perturbations[j][0]; probes[i].perturbed[j].y = 0.0; probes[i].perturbed[j].z = probes[i].z + perturbations[j][1]; copyProbeToPerturbed(&probes[i].perturbed[j], probes[i]); }
	}
	// Count
	int c = 0;
	// Loop
	while (c < totalSteps)
	{
		// Try euler
		multiChaosEuler(&e, &s, probes, probeNum, dt, c, row);
		// RK4: PROBABLY BROKEN
		//multiChaosRK4(&e, &s, probes, probeNum, dt);
		c++;
	}

	// Extract lyapunovs
	extractMultiLyapunov(probes, probeNum, lyapunovs, dt * totalSteps);

	return lyapunovs;

}
// Generates a plot and extracts lyapunov exponents for a single probe placed in L2
void singleProbeChaos(double dt, int totalSteps, bool offset, bool sphere)
{
	// Initialize
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe jw; jw.x = 0; jw.y = 0; jw.z = 1500000000; jw.vx = 0; jw.vy = 0; jw.vz = 30173.943; jw.perturbed = nullptr;

	initializePerturbations(&jw, sphere);

	// Store data: Store all data of every ITERATIONS step
	int n = (totalSteps + ITERATIONS - 1) / ITERATIONS;
	double **data = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i)
	{
		data[i] = (double *)malloc((3 + 3 * NUMPOINTS) * sizeof(double));
	}

	// Count, Data Count
	int c = 0;
	int dc = 0;

	// Loop
	while (c < totalSteps)
	{
		// Store data of JW and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			data[dc][0] = jw.x; data[dc][1] = jw.y; data[dc][2] = jw.z;
			for (int i = 0; i < NUMPOINTS; ++i) { data[dc][3 * i + 3] = jw.perturbed[i].x; data[dc][3 * i + 4] = jw.perturbed[i].y; data[dc][3 * i + 5] = jw.perturbed[i].z; }
			dc++;
		}

		// EULER
		chaosEuler(&e, &s, &jw, dt, c);

		// RK4
		//chaosRK4(&e, &s, &jw, dt);
		c++;
	}

	// Extract lyapunov exponents
	double lyapunov[2];
	extractSingleLyapunov(data, n, lyapunov, sphere, dt * totalSteps);

	// Write data to file
	perturbationsToFile(data, n, offset);

	// Free stuff
	for (int i = 0; i < n; ++i)
	{
		free(data[i]);
	}
	free(data);
}
//-------------------------------------------------------------------------------------------------------

void extractSingleLyapunov(double **data, int dataLen, double lyapunov[2], bool sphere, double t)
{
	// Calculate final vector distance to every point on ellipsoid
	double *dist = (double *)malloc(NUMPOINTS * sizeof(double));
	if (sphere)
	{
		for (int i = 0; i < NUMPOINTS; ++i)
		{
			double v[3] = {data[dataLen - 1][3 * i + 3] - data[dataLen - 1][0], data[dataLen - 1][3 * i + 4] - data[dataLen - 1][1], data[dataLen - 1][3 * i + 5] - data[dataLen - 1][2]};
			dist[i] = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
		}
	}
	else
	{
		for (int i = 0; i < NUMPOINTS; ++i)
		{
			// Only need x and z values
			double v[2] = {data[dataLen - 1][3 * i + 3] - data[dataLen - 1][0], data[dataLen - 1][3 * i + 5] - data[dataLen - 1][2]};
			dist[i] = sqrt((v[0] * v[0]) + (v[1] * v[1]));
		}
	}

	// Find max and min
	double max = dist[0];
	double min = dist[0];
	for (int i = 1; i < NUMPOINTS; ++i)
	{
		if (dist[i] > max)
		{
			max = dist[i];
		}
		if (dist[i] < min)
		{
			min = dist[i];
		}
	}

	std::cout << "Printing single probe initial position" << std::endl;
	std::cout << "x: " << data[0][0] << " z: " << data[0][2] << std::endl;

	std::cout << "Printing single probe final position" << std::endl;
	std::cout << "x: " << data[dataLen - 1][0] << " z: " << data[dataLen - 1][2] << std::endl;

	std::cout << "Printing single probe chaos max and min final vector distance" << std::endl;
	std::cout << "Max: " << max << " Min: " << min << std::endl;

	// Calculate Lyapunov max and min
	lyapunov[0] = (1.0 / t) * log(min / RADIUS);
	lyapunov[1] = (1.0 / t) * log(max / RADIUS);

	std::cout << "Printing single probe chaos max and min lyapunov" << std::endl;
	std::cout << "Max lyapunov: " << lyapunov[1] << std::endl;
	std::cout << "Min lyapunov: " << lyapunov[0] << std::endl;


}

void extractMultiLyapunov(Probe *probes, int probeNum, double **lyapunovs, double t)
{
	// Loop over every probe
	for (int i = 0; i < probeNum; ++i)
	{
		// Calculate distance from every perturbation to center probe
		double v[2] = {probes[i].perturbed[0].x - probes[i].x, probes[i].perturbed[0].z - probes[i].z};
		double dist = sqrt((v[0] * v[0]) + (v[1] * v[1]));
		// Search for max and min
		double max = dist; double min = dist;
		for (int j = 1; j < NUMPOINTS; ++j)
		{
			v[0] = probes[i].perturbed[j].x - probes[i].x; v[1] = probes[i].perturbed[j].z - probes[i].z;
			dist = sqrt((v[0] * v[0]) + (v[1] * v[1]));
			if (dist > max) { max = dist; } if (dist < min) { min = dist; }
		}
		// Calculate lyapunov exponents: [min, max]
		lyapunovs[i][0] = (1.0 / t) * log(min / RADIUS); lyapunovs[i][1] = (1.0 / t) * log(max / RADIUS);

		// DEBUGGING
		std::cout << "Printing final positions and lyapunovs" << std::endl;
		std::cout << "X pos: " << probes[i].x << " Z pos: " << probes[i].z << " Min lyapunov: " << lyapunovs[i][0] << " Max lyapunov: " << lyapunovs[i][1] << std::endl;
		std::cout << std::endl << std::endl;
	}
}

void calcForce(Planet e, Planet s, Planet jw, double *forces)
{
	double G = 6.6743e-11;

	// Earth Sun
	double r21[3] = {e.x - s.x, e.y - s.y, e.z - s.z};
	double dx = e.x - s.x, dy = e.y - s.y, dz = e.z - s.z;
	double magr21 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	double r21Unit[3] = {r21[0] / magr21, r21[1] / magr21, r21[2] / magr21};
	double gravity = (-1 * G * s.m * e.m) / (magr21 * magr21);

	forces[0] = gravity * r21Unit[0];
	forces[1] = gravity * r21Unit[1];
	forces[2] = gravity * r21Unit[2];
	forces[3] = 0;
	forces[4] = 0;
	forces[5] = 0;

	// James Webb
	double r31[3] = {jw.x - s.x, jw.y - s.y, jw.z - s.z};
	dx = jw.x - s.x; dy = jw.y - s.y; dz = jw.z - s.z;
	double magr31 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	double r31Unit[3] = {r31[0] / magr31, r31[1] / magr31, r31[2] / magr31};
	double gravity31 = (-1 * G * s.m * jw.m) / (magr31 * magr31);

	double r32[3] = {jw.x - e.x, jw.y - e.y, jw.z - e.z};
	double magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1]) + (r32[2] * r32[2]));
	double r32Unit[3] = {r32[0] / magr32, r32[1] / magr32, r32[2] / magr32};
	double gravity32 = (-1 * G * e.m * jw.m) / (magr32 * magr32);

	forces[6] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]);
	forces[7] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);
	forces[8] = (gravity31 * r31Unit[2]) + (gravity32 * r32Unit[2]);
}

// Calculates accelerations for earth, sun, james webb, and all its perturbations
void chaosCalcForce(Planet e, Planet s, Probe jw, double *accels)
{
	double G = 6.6743e-11;

	// Sun on Earth
	double r21[3] = {e.x - s.x, e.y - s.y, e.z - s.z};
	double magr21 = sqrt((r21[0] * r21[0]) + (r21[1] * r21[1]) + (r21[2] * r21[2]));
	double r21Unit[3] = {r21[0] / magr21, r21[1] / magr21, r21[2] / magr21};
	double gravity = (-1 * G * s.m) / (magr21 * magr21);
	accels[0] = gravity * r21Unit[0]; accels[1] = gravity * r21Unit[1]; accels[2] = gravity * r21Unit[2];
	// Sun on James Webb telescope
	double r31[3] = {jw.x - s.x, jw.y - s.y, jw.z - s.z};
	double magr31 = sqrt((r31[0] * r31[0]) + (r31[1] * r31[1]) + (r31[2] * r31[2]));
	double r31Unit[3] = {r31[0] / magr31, r31[1] / magr31, r31[2] / magr31};
	double gravity31 = (-1 * G * s.m) / (magr31 * magr31);
	// Earth on JW
	double r32[3] = {jw.x - e.x, jw.y - e.y, jw.z - e.z};
	double magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1]) + (r32[2] * r32[2]));
	double r32Unit[3] = {r32[0] / magr32, r32[1] / magr32, r32[2] / magr32};
	double gravity32 = (-1 * G * e.m) / (magr32 * magr32);

	accels[3] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]);
	accels[4] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);
	accels[5] = (gravity31 * r31Unit[2]) + (gravity32 * r32Unit[2]);

	// Perturbation Loop
	for (int i = 0; i < NUMPOINTS; ++i)
	{
		// Sun on Perturbation
		double r31[3] = {jw.perturbed[i].x - s.x, jw.perturbed[i].y - s.y, jw.perturbed[i].z - s.z};
		double magr31 = sqrt((r31[0] * r31[0]) + (r31[1] * r31[1]) + (r31[2] * r31[2]));
		double r31Unit[3] = {r31[0] / magr31, r31[1] / magr31, r31[2] / magr31};
		double gravity31 = (-1 * G * s.m) / (magr31 * magr31);
		// Earth on Perturbation
		double r32[3] = {jw.perturbed[i].x - e.x, jw.perturbed[i].y - e.y, jw.perturbed[i].z - e.z};
		double magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1]) + (r32[2] * r32[2]));
		double r32Unit[3] = {r32[0] / magr32, r32[1] / magr32, r32[2] / magr32};
		double gravity32 = (-1 * G * e.m) / (magr32 * magr32);

		accels[3 * i + 6] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]);
		accels[3 * i + 7] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);
		accels[3 * i + 8] = (gravity31 * r31Unit[2]) + (gravity32 * r32Unit[2]);
	}
}

void multiChaosCalcAccel(Planet e, Planet s, Probe *probes, int probeNum, double **accels)
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
		double r32[2] = {probes[i].x - e.x, probes[i].z - e.z}; double magr32 = sqrt((r32[0] * r32[0]) + (r32[1] * r32[1])); double r32Unit[2] = {r32[0] / magr32, r32[1] / magr32}; double gravity32 = (-1 * G * e.m) / (magr32 * magr32);
		// Net acceleration
		accels[i + 1][0] = (gravity31 * r31Unit[0]) + (gravity32 * r32Unit[0]); accels[i + 1][1] = (gravity31 * r31Unit[1]) + (gravity32 * r32Unit[1]);

		// Calc for all perturbations of probe
		for (int j = 0; j < NUMPOINTS; ++j)
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

void rk4(Planet *e, Planet *s, Planet *jw, double dt)
{

	// k1 step
	// Get net force at t0
	// Calculate the net force at the starting position
	double k1Force[9];
	calcForce(*e, *s, *jw, k1Force);

	// K1 velocity
	double k1_vel[6];
	k1_vel[0] = e->vx; k1_vel[1] = e->vy; k1_vel[2] = e->vz;
	k1_vel[3] = s->vx; k1_vel[4] = s->vy; k1_vel[5] = s->vz;
	k1_vel[6] = jw->vx; k1_vel[7] = jw->vy; k1_vel[8] = jw->vz;

	// k2 step
	// Make new mass_list for k2 state.
	Planet k2Planets[3];

	// Copy properties of mass_list
	copyPlanet(&k2Planets[0], *e);
	copyPlanet(&k2Planets[1], *s);
	copyPlanet(&k2Planets[2], *jw);

	// Push k2 state forward dt / 2 using k1 as net force

	// Update velocities and positions: THIS USES k1 / 2
	// Update positions
	for (int i = 0; i < 3; ++i)
	{
		k2Planets[i].vx += (k1Force[3 * i] / k2Planets[i].m) * 0.5 * dt; k2Planets[i].vy += (k1Force[3 * i + 1] / k2Planets[i].m) * 0.5 * dt; k2Planets[i].vz += (k1Force[3 * i + 2] / k2Planets[i].m) * 0.5 * dt;
		k2Planets[i].x += k1_vel[3 * i] * 0.5 * dt; k2Planets[i].y += k1_vel[3 * i + 1] * 0.5 * dt; k2Planets[i].z += k1_vel[3 * i + 2] * 0.5 * dt;
	}

	// Calculate k2 force
	double k2Force[9];
	calcForce(k2Planets[0], k2Planets[1], k2Planets[2], k2Force);

	// K2 velocity
	double k2_vel[9];
	k2_vel[0] = k2Planets[0].vx; k2_vel[1] = k2Planets[0].vy; k2_vel[2] = k2Planets[0].vz;
	k2_vel[3] = k2Planets[1].vx; k2_vel[4] = k2Planets[1].vy; k2_vel[5] = k2Planets[1].vz;
	k2_vel[6] = k2Planets[2].vx; k2_vel[7] = k2Planets[2].vy; k2_vel[8] = k2Planets[2].vz;

	// k3 step
	// Make new mass_list for k2 state.
	Planet k3Planets[3];

	// Copy properties of mass_list
	copyPlanet(&k3Planets[0], *e);
	copyPlanet(&k3Planets[1], *s);
	copyPlanet(&k3Planets[2], *jw);

	// Push k2 state forward dt / 2 using k1 as net force

	// Update velocities and positions: THIS USES k1 / 2
	for (int i = 0; i < 3; ++i)
	{
		k3Planets[i].vx += (k2Force[3 * i] / k3Planets[i].m) * 0.5 * dt; k3Planets[i].vy += (k2Force[3 * i + 1] / k3Planets[i].m) * 0.5 * dt; k3Planets[i].vz += (k2Force[3 * i + 2] / k3Planets[i].m) * 0.5 * dt;
		k3Planets[i].x += k2_vel[3 * i] * 0.5 * dt; k3Planets[i].y += k2_vel[3 * i + 1] * 0.5 * dt; k3Planets[i].z += k2_vel[3 * i + 2] * 0.5 * dt;
	}

	// Calculate k3 force
	double k3Force[9];
	calcForce(k3Planets[0], k3Planets[1], k3Planets[2], k3Force);

	// k3 vel
	double k3_vel[9];
	k3_vel[0] = k3Planets[0].vx; k3_vel[1] = k3Planets[0].vy; k3_vel[2] = k3Planets[0].vz;
	k3_vel[3] = k3Planets[1].vx; k3_vel[4] = k3Planets[1].vy; k3_vel[5] = k3Planets[1].vz;
	k3_vel[6] = k3Planets[2].vx; k3_vel[7] = k3Planets[2].vy; k3_vel[8] = k3Planets[2].vz;

	// k4 step
	// Make new mass_list for k2 state.
	Planet k4Planets[3];

	// Copy properties of mass_list
	copyPlanet(&k4Planets[0], *e);
	copyPlanet(&k4Planets[1], *s);
	copyPlanet(&k4Planets[2], *jw);

	// Push k4 state forward dt using k3 as net force

	// Update velocities: THIS USES k3
	for (int i = 0; i < 3; ++i)
	{
		k4Planets[i].vx += (k3Force[3 * i] / k4Planets[i].m) * dt; k4Planets[i].vy += (k3Force[3 * i + 1] / k4Planets[i].m) * dt; k4Planets[i].vz += (k3Force[3 * i + 2] / k4Planets[i].m) * dt;
		k4Planets[i].x += k3_vel[3 * i] * dt; k4Planets[i].y += k3_vel[3 * i + 1] * dt; k4Planets[i].z += k3_vel[3 * i + 2] * dt;
	}

	// Calculate k4 force
	double k4Force[9];
	calcForce(k4Planets[0], k4Planets[1], k4Planets[2], k4Force);

	// k4 vel
	double k4_vel[6];
	k4_vel[0] = k4Planets[0].vx; k4_vel[1] = k4Planets[0].vy; k4_vel[2] = k4Planets[0].vz;
	k4_vel[3] = k4Planets[1].vx; k4_vel[4] = k4Planets[1].vy; k4_vel[5] = k4Planets[1].vz;
	k4_vel[6] = k4Planets[2].vx; k4_vel[7] = k4Planets[2].vy; k4_vel[8] = k4Planets[2].vz;

	// Calculate the weighted rk4 net force (k1_force + 2*k2_force + 2*k3_force + k4_force)
	double rk4Force[9];
	double rk4_vel[9];
	for (int i = 0; i < 9; ++i)
	{
		rk4Force[i] = k1Force[i] + 2 * k2Force[i] + 2 * k3Force[i] + k4Force[i];
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

	// Calculate the weighted rk4 velocity (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

	// Push the initial state forward dt using the weighted average rk4 force as net force

	// Update velocities: THIS USES k3
	e->vx += (rk4Force[0] / e->m / 6.0) * dt; e->vy += (rk4Force[1] / e->m / 6.0) * dt; e->vz += (rk4Force[2] / e->m / 6.0) * dt;
	s->vx += (rk4Force[3] / s->m / 6.0) * dt; s->vy += (rk4Force[4] / s->m / 6.0) * dt; s->vz += (rk4Force[5] / s->m / 6.0) * dt;
	jw->vx += (rk4Force[6] / jw->m / 6.0) * dt; jw->vy += (rk4Force[7] / jw->m / 6.0) * dt; jw->vz += (rk4Force[8] / jw->m / 6.0) * dt;


	// Update positions: Using k3 velocity. Initial position velocity
	e->x += (rk4_vel[0] / 6.0) * dt; e->y += (rk4_vel[1] / 6.0) * dt; e->z += (rk4_vel[2] / 6.0) * dt;
	s->x += (rk4_vel[3] / 6.0) * dt; s->y += (rk4_vel[4] / 6.0) * dt; s->z += (rk4_vel[5] / 6.0) * dt;
	jw->x += (rk4_vel[6] / 6.0) * dt; jw->y += (rk4_vel[7] / 6.0) * dt; jw->z += (rk4_vel[8] / 6.0) * dt;
}

void euler(Planet *e, Planet *s, Planet *jw, double dt)
{
	double forces[9];

	calcForce(*e, *s, *jw, forces);

	e->vx += forces[0] / e->m * dt; e->vy += forces[1] / e->m * dt; e->vz += forces[2] / e->m * dt;
	s->vx += forces[3] / s->m * dt; s->vy += forces[4] / s->m * dt; s->vz += forces[5] / s->m * dt;
	jw->vx += forces[6] / jw->m * dt; jw->vy += forces[7] / jw->m * dt; jw->vz += forces[8] / jw->m * dt;

	e->x += e->vx * dt; e->y += e->vy * dt; e->z += e->vz * dt;
	s->x += s->vx * dt; s->y += s->vy * dt; s->z += s->vz * dt;
	jw->x += jw->vx * dt; jw->y += jw->vy * dt; jw->z += jw->vz * dt;
}

// Chaos euler only applies acceleration, doesn't calculate it
void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt, int loopCount)
{
	// Assign memory to acceleration array for Euler
	double *accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double));
	chaosCalcForce(*e, *s, *jw, accels);

	std::setprecision(5);
	std::cout << "SINGLE CHAOS EULER" << std::endl;
	std::cout << "Loop " << loopCount << std::endl;
	std::cout << "Printing position of probe" << std::endl;
	std::cout << "x: " << jw->x << " y: " << jw->y << "z: " << jw->z << std::endl << std::endl;
	std::cout << "Printing velocities of probe" << std::endl;
	std::cout << "vx: " << jw->vx << " vy: " << jw->vy << " vz: " << jw->vz << std::endl << std::endl;
	std::cout << "Printing accelerations on probe" << std::endl;
	std::cout << "ax: " << accels[3] << " ay: " << accels[4] << "az: " << accels[5] << std::endl << std::endl;
	std::cout << "Printing accelerations on earth" << std::endl;
	std::cout << "ax: " << accels[0] << " ay: " << accels[1] << "az: " << accels[2] << std::endl << std::endl;


	e->vx += accels[0] * dt; e->vy += accels[1] * dt; e->vz += accels[2] * dt;
	e->x += e->vx * dt; e->y += e->vy * dt; e->z += e->vz * dt;
	jw->vx += accels[3] * dt; jw->vy += accels[4] * dt; jw->vz += accels[5] * dt;
	jw->x += jw->vx * dt; jw->y += jw->vy * dt; jw->z += jw->vz * dt;

	for (int i = 0; i < NUMPOINTS; ++i)
	{
		jw->perturbed[i].vx += accels[3 * i + 6] * dt; jw->perturbed[i].vy += accels[3 * i + 7] * dt; jw->perturbed[i].vz += accels[3 * i + 8] * dt;
		jw->perturbed[i].x += jw->perturbed[i].vx * dt; jw->perturbed[i].y += jw->perturbed[i].vy * dt; jw->perturbed[i].z += jw->perturbed[i].vz * dt;
	}

	free(accels);
}

void multiChaosEuler(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int loopCount, int row)
{
	double **accels = (double **)malloc((1 + probeNum) * sizeof(double *)); accels[0] = (double *)malloc(2 * sizeof(double)); for (int i = 1; i < 1 + probeNum; ++i) { accels[i] = (double *)malloc((2 + 2 * NUMPOINTS) * sizeof(double)); }
	multiChaosCalcAccel(*e, *s, probes, probeNum, accels);

	// DEBUGGING
	if (row == 4)
	{
		std::setprecision(5);
		std::cout << "MULTI CHAOS EULER" << std::endl;
		std::cout << "Loop " << loopCount << std::endl;
		std::cout << "Printing position of probe" << std::endl;
		std::cout << "x: " << probes[2].x << " y: " << probes[2].y << "z: " << probes[2].z << std::endl << std::endl;
		std::cout << "Printing velocities of probe" << std::endl;
		std::cout << "vx: " << probes[2].vx << " vy: " << probes[2].vy << " vz: " << probes[2].vz << std::endl << std::endl;
		std::cout << "Printing accelerations on probe" << std::endl;
		std::cout << "ax: " << accels[3][0] << " ay: " << 0.0 << "az: " << accels[3][1] << std::endl << std::endl;
		std::cout << "Printing accelerations on earth" << std::endl;
		std::cout << "ax: " << accels[0][0] << " ay: " << 0.0 << "az: " << accels[0][1] << std::endl << std::endl;
	}
	e->vx += accels[0][0] * dt; e->vz += accels[0][1] * dt;
	e->x += e->vx * dt; e->z += e->vz * dt;
	for (int i = 1; i < probeNum; ++i)
	{
		probes[i].vx += accels[i][0] * dt; probes[i].vz += accels[i][1] * dt;
		probes[i].x += probes[i].vx * dt; probes[i].z += probes[i].vz * dt;
		for (int j = 0; j < NUMPOINTS; ++j)
		{
			probes[i].perturbed[j].vx += accels[i + 1][2 * j + 2] * dt; probes[i].perturbed[j].vz += accels[i + 1][2 * j + 3] * dt;
			probes[i].perturbed[j].x += probes[i].perturbed[j].vx * dt; probes[i].perturbed[j].z += probes[i].perturbed[j].vz * dt;
		}
	}
}

// Moves objects forward using kAccels and kVels. Step for rk4
void rk4StepSingleChaos(Planet *eCopy, Probe *jwCopy, Planet e, Probe jw, double *accels, double dt)
{
	// DEBUGGING
//	std::cout << "Before" << std::endl;
//	std::cout << "JW pos: " << jwCopy->x << " " << jwCopy->y << " " << jwCopy->z << " vel: " << jwCopy->vx << " " << jwCopy->vy << " " << jwCopy->vz << std::endl;
	eCopy->vx += accels[0] * dt; eCopy->vy += accels[1] * dt; eCopy->vz += accels[2] * dt;
	jwCopy->vx += accels[3] * dt; jwCopy->vy += accels[4] * dt; jwCopy->vz += accels[5] * dt;
	eCopy->x += e.vx * dt; eCopy->y += e.vy * dt; eCopy->z += e.vz * dt;
	jwCopy->x += jw.vx * dt; jwCopy->y += jw.vy * dt; jwCopy->z += jw.vz * dt;
	// DEBUGGING
//	std::cout << "After" << std::endl;
//	std::cout << "JW pos: " << jwCopy->x << " " << jwCopy->y << " " << jwCopy->z << " vel: " << jwCopy->vx << " " << jwCopy->vy << " " << jwCopy->vz << std::endl;
//	std::cout << std::endl;

	for (int i = 0; i < NUMPOINTS; ++i)
	{
		jwCopy->perturbed[i].vx += accels[3 * i + 6] * dt; jwCopy->perturbed[i].vy += accels[3 * i + 7] * dt; jwCopy->perturbed[i].vz += accels[3 * i + 8] * dt;
		jwCopy->perturbed[i].x += jw.perturbed[i].vx * dt; jwCopy->perturbed[i].y += jw.perturbed[i].vy * dt; jwCopy->perturbed[i].z += jw.perturbed[i].vz * dt;
	}
}

void rk4StepMultiChaos(Planet *eCopy, Probe *probeCopies, Planet e, Probe *probes, int probeNum, double **accels, double dt)
{
	// First index is earth accels
	eCopy->vx += accels[0][0] * dt; eCopy->vz += accels[0][1] * dt; eCopy->x += e.vx * dt; eCopy->z += e.vz * dt;
	// Step for probes and perturbations
	for (int i = 0; i < probeNum; ++i)
	{
		// Probe step
		probeCopies[i].vx += accels[i + 1][0] * dt; probeCopies[i].vz += accels[i + 1][1] * dt; probeCopies[i].x += probes[i].vx * dt; probeCopies[i].z += probes[i].vz * dt;
		// Perturbations step
		for (int j = 0; j < NUMPOINTS; ++j) { probeCopies[i].perturbed[j].vx += accels[i + 1][2 * j + 2] * dt; probeCopies[i].perturbed[j].vz += accels[i + 1][2 * j + 3] * dt; probeCopies[i].perturbed[j].x += probes[i].perturbed[j].vx * dt; probeCopies[i].perturbed[j].z += probes[i].perturbed[j].vz * dt; }
	}
}

// RK4 for a single probe with perturbations
void chaosRK4(Planet *e, Planet *s, Probe *jw, double dt)
{
	// Initialize Object arrays
	Planet k2Planets[2], k3Planets[2], k4Planets[2];
	Probe k2Probe, k3Probe, k4Probe;
	// Copy data into planet arrays for all k
	copyPlanet(&k2Planets[0], *e); copyPlanet(&k2Planets[1], *s); copyPlanet(&k3Planets[0], *e); copyPlanet(&k3Planets[1], *s); copyPlanet(&k4Planets[0], *e); copyPlanet(&k3Planets[1], *s);
	// Copy data into probes for all k
	copyProbe(&k2Probe, *jw); copyProbe(&k3Probe, *jw); copyProbe(&k4Probe, *jw);
	// Assign memory for perturbed probes for all k
	k2Probe.perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe)); k3Probe.perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe)); k4Probe.perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
	// Initialize probe perturbations for all k copies
	for (int i = 0; i < NUMPOINTS; ++i)
	{
		copyProbe(&k2Probe.perturbed[i], jw->perturbed[i]); copyProbe(&k3Probe.perturbed[i], jw->perturbed[i]); copyProbe(&k4Probe.perturbed[i], jw->perturbed[i]);
	}

	// Initialize KAccels
	double *k1Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k2Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k3Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k4Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double));

	// k1 step
	// Put accelerations in array
	chaosCalcForce(*e, *s, *jw, k1Accels);

	// k2 Step
	// Move k2 object copies dt / 2 using k1 Accelerations and k1 Velocities

	// DEBUGGING
//	std::cout << "Debugging k2 Step" << std::endl;
	rk4StepSingleChaos(&k2Planets[0], &k2Probe, *e, *jw, k1Accels, dt / 2.0);
	// Put data in arrays
	chaosCalcForce(k2Planets[0], k2Planets[1], k2Probe, k2Accels);

	// k3 Step
	// Move k3 object copies dt / 2 using k2 Accelerations and k2 Velocities
	rk4StepSingleChaos(&k3Planets[0], &k3Probe, k2Planets[0], k2Probe, k2Accels, dt / 2.0);
	// Put data in arrays
	chaosCalcForce(k3Planets[0], k3Planets[1], k3Probe, k3Accels);

	// k4 Step
	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
	rk4StepSingleChaos(&k4Planets[0], &k4Probe, k3Planets[0], k3Probe, k3Accels, dt);
	// Put data in arrays
	chaosCalcForce(k4Planets[0], k4Planets[1], k4Probe, k4Accels);

	// rk4 Step
	// Move objects

	// DEBUGGING
//	std::cout << "Before: " << std::endl;
//	std::cout << "Earth pos: " << e->x << " " << e->y << " " << e->z << " vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
//	std::cout << "JW pos: " << jw->x << " " << jw->y << " " << jw->z << " vel: " << jw->vx << " " << jw->vy << " " << jw->vz << std::endl;
//	std::cout << "JW Copies" << std::endl;
//	std::cout << "k2Probe pos: " << k2Probe.x << " " << k2Probe.y << " " << k2Probe.z << " vel: " << k2Probe.vx << " " << k2Probe.vy << " " << k2Probe.vz << std::endl;
//	std::cout << "k3Probe pos: " << k3Probe.x << " " << k3Probe.y << " " << k3Probe.z << " vel: " << k3Probe.vx << " " << k3Probe.vy << " " << k3Probe.vz << std::endl;
//	std::cout << "k4Probe pos: " << k4Probe.x << " " << k4Probe.y << " " << k4Probe.z << " vel: " << k4Probe.vx << " " << k4Probe.vy << " " << k4Probe.vz << std::endl;
//	std::cout << std::endl << std::endl;
	e->vx += (k1Accels[0] + 2 * k2Accels[0] + 2 * k3Accels[0] + k4Accels[0]) * dt / 6.0; e->vy += (k1Accels[1] + 2 * k2Accels[1] + 2 * k3Accels[1] + k4Accels[1]) * dt / 6.0; e->vz += (k1Accels[2] + 2 * k2Accels[2] + 2 * k3Accels[2] + k4Accels[2]) * dt / 6.0;
	jw->vx += (k1Accels[3] + 2 * k2Accels[3] + 2 * k3Accels[3] + k4Accels[3]) * dt / 6.0; jw->vy += (k1Accels[4] + 2 * k2Accels[4] + 2 * k3Accels[4] + k4Accels[4]) * dt / 6.0; jw->vz += (k1Accels[5] + 2 * k2Accels[5] + 2 * k3Accels[5] + k4Accels[5]) * dt / 6.0;
	e->x += (e->vx + 2 * k2Planets[0].vx + 2 * k3Planets[0].vx + k4Planets[0].vx) * dt / 6.0; e->y += (e->vy + 2 * k2Planets[0].vy + 2 * k3Planets[0].vy + k4Planets[0].vy) * dt / 6.0; e->z += (e->vz + 2 * k2Planets[0].vz + 2 * k3Planets[0].vz + k4Planets[0].vz) * dt / 6.0;
	jw->x += (jw->vx + 2 * k2Probe.vx + 2 * k3Probe.vx + k4Probe.vx) * dt / 6.0; jw->y += (jw->vy + 2 * k2Probe.vy + 2 * k3Probe.vy + k4Probe.vy) * dt / 6.0; jw->z += (jw->vz + 2 * k2Probe.vz + 2 * k3Probe.vz + k4Probe.vz) * dt / 6.0;
//	std::cout << "After: " << std::endl;
//	std::cout << "Earth pos: " << e->x << " " << e->y << " " << e->z << " vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
//	std::cout << "JW pos: " << jw->x << " " << jw->y << " " << jw->z << " vel: " << jw->vx << " " << jw->vy << " " << jw->vz << std::endl;
	// Move perturbations
	for (int i = 0; i < NUMPOINTS; ++i)
	{
		jw->perturbed[i].vx += (k1Accels[3 * i + 6] + 2 * k2Accels[3 * i + 6] + 2 * k3Accels[3 * i + 6] + k4Accels[3 * i + 6]) * dt / 6.0; jw->perturbed[i].vy += (k1Accels[3 * i + 7] + 2 * k2Accels[3 * i + 7] + 2 * k3Accels[3 * i + 7] + k4Accels[3 * i + 7]) * dt / 6.0; jw->perturbed[i].vz += (k1Accels[3 * i + 8] + 2 * k2Accels[3 * i + 8] + 2 * k3Accels[3 * i + 8] + k4Accels[3 * i + 8]) * dt / 6.0;
		jw->perturbed[i].x += (jw->perturbed[i].vx + 2 * k2Probe.perturbed[i].vx + 2 * k3Probe.perturbed[i].vx + k4Probe.perturbed[i].vx) * dt / 6.0; jw->perturbed[i].y += (jw->perturbed[i].vy + 2 * k2Probe.perturbed[i].vy + 2 * k3Probe.perturbed[i].vy + k4Probe.perturbed[i].vy) * dt / 6.0; jw->perturbed[i].z += (jw->perturbed[i].vz + 2 * k2Probe.perturbed[i].vz + 2 * k3Probe.perturbed[i].vz + k4Probe.perturbed[i].vz) * dt / 6.0;
	}

	// Free stuff
	free(k1Accels); free(k2Accels); free(k3Accels); free(k4Accels);
}

void multiChaosRK4(Planet *e, Planet *s, Probe *probes, int probeNum, double dt)
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
		k1Accels[i + 1] = (double *)malloc((2 + 2 * NUMPOINTS) * sizeof(double)); k2Accels[i + 1] = (double *)malloc((2 + 2 * NUMPOINTS) * sizeof(double)); k3Accels[i + 1] = (double *)malloc((2 + 2 * NUMPOINTS) * sizeof(double)); k4Accels[i + 1] = (double *)malloc((2 + 2 * NUMPOINTS) * sizeof(double));
		// Copy for Probe copies and initialize perturbed arrays for each probe
		copyProbe(&k2Probes[i], probes[i]); k2Probes[i].perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe)); copyProbe(&k3Probes[i], probes[i]); k3Probes[i].perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe)); copyProbe(&k4Probes[i], probes[i]); k4Probes[i].perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
		// Copy for Perturbation copies
		for (int j = 0; j < NUMPOINTS; ++j) { copyProbe(&k2Probes[i].perturbed[j], probes[i].perturbed[j]); copyProbe(&k3Probes[i].perturbed[j], probes[i].perturbed[j]); copyProbe(&k4Probes[i].perturbed[j], probes[i].perturbed[j]); }
	}

	// k1 step
	multiChaosCalcAccel(*e, *s, probes, probeNum, k1Accels);
	// k2 Step
	rk4StepMultiChaos(&k2Planets[0], k2Probes, *e, probes, probeNum, k1Accels, dt / 2);
	multiChaosCalcAccel(k2Planets[0], k2Planets[1], k2Probes, probeNum, k2Accels);
	// k3 Step
	rk4StepMultiChaos(&k3Planets[0], k3Probes, k2Planets[0], k2Probes, probeNum, k2Accels, dt / 2);
	multiChaosCalcAccel(k3Planets[0], k3Planets[1], k3Probes, probeNum, k3Accels);
	// k4 Step
	rk4StepMultiChaos(&k4Planets[0], k4Probes, k3Planets[0], k3Probes, probeNum, k3Accels, dt);
	multiChaosCalcAccel(k4Planets[0], k4Planets[1], k4Probes, probeNum, k4Accels);
	// rk4 Step
	// Move objects
	e->vx += (k1Accels[0][0] + 2 * k2Accels[0][0] + 2 * k3Accels[0][0] + k4Accels[0][0]) * dt / 6.0; e->vz += (k1Accels[0][1] + 2 * k2Accels[0][1] + 2 * k3Accels[0][1] + k4Accels[0][1]) * dt / 6.0; e->x += (e->vx + 2 * k2Planets[0].vx + 2 * k3Planets[0].vx + k4Planets[0].vx) * dt / 6.0;  e->z += (e->vz + 2 * k2Planets[0].vz + 2 * k3Planets[0].vz + k4Planets[0].vz) * dt / 6.0;
	// Probes and perturbations
	for (int i = 0; i < probeNum; ++i)
	{
		probes[i].vx += (k1Accels[i + 1][0] + 2 * k2Accels[i + 1][0] + 2 * k3Accels[i + 1][0] + k4Accels[i + 1][0]) * dt / 6.0; probes[i].vz += (k1Accels[i + 1][1] + 2 * k2Accels[i + 1][1] + 2 * k3Accels[i + 1][1] + k4Accels[i + 1][1]) * dt / 6.0; probes[i].x += (probes[i].vx + 2 * k2Probes[i].vx + 2 * k3Probes[i].vx + k4Probes[i].vx) * dt / 6.0; probes[i].z += (probes[i].vz + 2 * k2Probes[i].vz + 2 * k3Probes[i].vz + k4Probes[i].vz) * dt / 6.0;
		// Perturbations
		for (int j = 0; j < NUMPOINTS; ++j)
		{
			probes[i].perturbed[j].vx += (k1Accels[i + 1][2 * j + 2] + 2 * k2Accels[i + 1][2 * j + 2] + 2 * k3Accels[i + 1][2 * j + 2] + k4Accels[i + 1][2 * j + 2]) * dt / 6.0; probes[i].perturbed[j].vz += (k1Accels[i + 1][2 * j + 3] + 2 * k2Accels[i + 1][2 * j + 3] + 2 * k3Accels[i + 1][2 * j + 3] + k4Accels[i + 1][2 * j + 3]) * dt / 6.0; probes[i].perturbed[j].x += (probes[i].perturbed[j].vx + 2 * k2Probes[i].perturbed[j].vx + 2 * k3Probes[i].perturbed[j].vx + k4Probes[i].perturbed[j].vx) * dt / 6.0; probes[i].perturbed[j].z += (probes[i].perturbed[j].vz + 2 * k2Probes[i].perturbed[j].vz + 2 * k3Probes[i].perturbed[j].vz + k4Probes[i].perturbed[j].vz) * dt / 6.0;
		}
	}

	// Free stuff
	// Probes and Perturbations and Accels
	free(k1Accels[0]); free(k2Accels[0]); free(k3Accels[0]); free(k4Accels[0]); for (int i = 0; i < probeNum; ++i) { free(k2Probes[i].perturbed); free(k3Probes[i].perturbed); free(k4Probes[i].perturbed); free(k1Accels[i + 1]); free(k2Accels[i + 1]); free(k3Accels[i + 1]); free(k4Accels[i + 1]); } free(k1Accels); free(k2Accels); free(k3Accels); free(k4Accels); free(k2Probes); free(k3Probes); free(k4Probes);
}


//-------------------------------------------------------------------------------------------------------------------------------
// Utility functions

// Returns an array of x, z probe positions distributed across L1 to L2
void partitionPositions(double *positions, int resolution[2], int row)
{
	// L1 to L2 is 3 mil km
	double xIncr = 3000000000.0 / (double)(resolution[0] - 1); double zIncr = 3000000000.0 / (double)(resolution[1] - 1);
	double zPos = -1500000000.0 + zIncr * row;
	for (int i = 0; i < resolution[0]; ++i)
	{
		positions[2 * i] = -1500000000.0 + xIncr * i;
		positions[2 * i + 1] = zPos;
	}
}

// Transfers from lyapunov arr to data array to write to file later
void multiLyapunovToArr(double **data, int dataLen, double **lyapunovs, int probeNum)
{
	for (int i = dataLen; i < probeNum + dataLen; ++i) { data[i][0] = lyapunovs[i - dataLen][0]; data[i][1] = lyapunovs[i - dataLen][1]; }
}

bool orbitsToFile(std::string fileName, double **data, int dataLen)
{
	std::ofstream file;
	file.open(fileName);

	file << "ex" << "," << "ey" << "," << "ez" << "," << "sx" << "," << "sy" << "," << "sz" << "," << "jwx" << "," << "jwy" << "," << "jwz" << "," << std::endl;
	for (int i = 0; i < dataLen; ++i)
	{
		file << std::to_string(data[i][0]) << "," << std::to_string(data[i][1]) << "," << std::to_string(data[i][2]) << ","
			 << std::to_string(data[i][3]) << "," << std::to_string(data[i][4]) << "," << std::to_string(data[i][5]) << ","
		     << std::to_string(data[i][6]) << "," << std::to_string(data[i][7]) << "," << std::to_string(data[i][8]) << std::endl;
	}
	file.close();
	return true;
}

void multiLyapunovToFile(double **data, int dataLen)
{
	std::ofstream file;
	file.open("lyapunovChaosPlot.csv", std::ios_base::app);

	for (int i = 0; i < dataLen; ++i)
	{
		file << std::to_string(data[i][0]) << "," << std::to_string(data[i][1]) << "," << std::to_string(data[i][1] - data[i][0]) << std::endl;
	}
	file.close();
}
// Writes data points for evolution of single James Webb with perturbations to file.
void perturbationsToFile(double **data, int dataLen, bool offset)
{
	std::string fileName1 = "singleChaosData\\probe.csv";
	std::string fileName2 = "singleChaosData\\perturbations.csv";
	std::ofstream file1;
	std::ofstream file2;
	file1.open(fileName1);
	file2.open(fileName2);

	std::cout << "Data length: " << dataLen << std::endl;

	file1 << "x" << "," << "y" << "," << "z" << std::endl;
	file2 << "x" << "," << "y" << "," << "z" << std::endl;

	if (offset)
	{
		for (int i = 0; i < dataLen; ++i)
		{
			file1 << std::to_string(data[i][0] - data[0][0]) << "," << std::to_string(data[i][1] - data[0][1]) << "," << std::to_string(data[i][2] - data[0][2]) << std::endl;

			// DEBUGGING
			for (int j = 3; j < 3 * NUMPOINTS + 3; j += 3)
			{
				file2 << std::to_string(data[i][j] - data[0][0]) << "," << std::to_string(data[i][j + 1] - data[0][1]) << "," << std::to_string(data[i][j + 2] - data[0][2]) << std::endl;
			}
		}
	}
	else
	{
		for (int i = 0; i < dataLen; ++i)
		{
			file1 << std::to_string(data[i][0]) << "," << std::to_string(data[i][1]) << "," << std::to_string(data[i][2]) << std::endl;
			for (int j = 3; j < 3 * NUMPOINTS + 3; j += 3)
			{
				file2 << std::to_string(data[i][j]) << "," << std::to_string(data[i][j + 1]) << "," << std::to_string(data[i][j + 2]) << std::endl;
			}

		}

		// DEBUGGING
//		file1 << std::to_string(data[0][0]) << "," << std::to_string(data[0][1]) << "," << std::to_string(data[0][2] - data[0][2]) << std::endl;
//		for (int j = 3; j < 3 * NUMPOINTS + 3; j += 3)
//		{
//			file2 << std::to_string(data[0][j]) << "," << std::to_string(data[0][j + 1]) << "," << std::to_string(data[0][j + 2] - data[0][2]) << std::endl;
//		}
	}
	file1.close();
	file2.close();
}

// Returns an array structure of the perturbations read from file. Faster to read from array for repeated use. Struct: data[NUMPOINTS][x, z]
void returnPerturbations(double **perturbations)
{
	// Read from file
	FILE *f = fopen("unitCircle.csv", "r");
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);

	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);

	// Read from buffer into struct
	int pointCount = 0;
	char arr[15];
	int arrCount = 0;
	int commaCount = 0;
	int i = 0;

	while (pointCount < NUMPOINTS)
	{
		if (buffer[i] == '\n')
		{
			arr[arrCount] = '\0';
			perturbations[pointCount][1] = std::atof(arr);
			pointCount++;
			arrCount = 0;
			commaCount = 0;

		}
		else if (buffer[i] == ',')
		{
			arr[arrCount] = '\0';
			if (commaCount == 0)
			{
				perturbations[pointCount][0] = std::atof(arr);
			}
			arrCount = 0;
			commaCount++;
		}
		else
		{
			arr[arrCount] = buffer[i];
			arrCount++;
		}
		++i;
	}

	free(buffer);
}

void freePerturbations(double **perturbations)
{
	for (int i = 0; i < NUMPOINTS; ++i) { free(perturbations[i]); }
	free(perturbations);
}

// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid
void initializePerturbations(Probe *p, bool sphere)
{
	// Initialize all perturbation probes
	p->perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
	for (int i = 0; i < NUMPOINTS; ++i) { copyProbeToPerturbed(&p->perturbed[i], *p); }

	// Read from file
	FILE *f;
	if (sphere) { f = fopen("unitSphere.csv", "r"); }
	else { f = fopen("unitCircle.csv", "r"); }
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);

	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);

	// Read from buffer into perturbed probes
	int pointCount = 0;
	char arr[15];
	int arrCount = 0;
	int commaCount = 0;
	int i = 0;

	if (sphere)
	{
		while (pointCount < NUMPOINTS)
		{
			if (buffer[i] == '\n')
			{
				arr[arrCount] = '\0';
				p->perturbed[pointCount].z = std::atof(arr) + p->z;
				pointCount++;
				arrCount = 0;
				commaCount = 0;

			}
			else if (buffer[i] == ',')
			{
				arr[arrCount] = '\0';
				if (commaCount == 0)
				{
					p->perturbed[pointCount].x = std::atof(arr) + p->x;
				}
				else
				{
					p->perturbed[pointCount].y = std::atof(arr) + p->y;
				}
				arrCount = 0;
				commaCount++;
			}
			else
			{
				arr[arrCount] = buffer[i];
				arrCount++;
			}
			++i;
		}
	}
	else
	{
		while (pointCount < NUMPOINTS)
		{
			if (buffer[i] == '\n')
			{
				arr[arrCount] = '\0';
				p->perturbed[pointCount].z = std::atof(arr) + p->z;
				pointCount++;
				arrCount = 0;
				commaCount = 0;

			}
			else if (buffer[i] == ',')
			{
				arr[arrCount] = '\0';
				if (commaCount == 0)
				{
					p->perturbed[pointCount].x = std::atof(arr) + p->x;
				}
				else
				{
					p->perturbed[pointCount].y = 0.0;
				}
				arrCount = 0;
				commaCount++;
			}
			else
			{
				arr[arrCount] = buffer[i];
				arrCount++;
			}
			++i;
		}
	}

	// DEBUGGING
//	std::cout << "Printing JW telescope location" << std::endl;
//	std::cout << p->x << " " << p->y << " " << p->z << std::endl;
//
//	std::cout << "Printing perturbations" << std::endl;
//	for (int i = 0; i < NUMPOINTS; ++i)
//	{
//		std::cout << p->perturbed[i].x << " " << p->perturbed[i].y << " " << p->perturbed[i].z << std::endl;
//	}
//	std::cout << std::endl;
//	std::cout << "Velocities: " << p->perturbed[0].vx << " " << p->perturbed[0].vy << " " << p->perturbed[0].vz << std::endl;
//	std::cout << "Last" << std::endl;
//	std::cout << p->perturbed[NUMPOINTS - 1].x << " " << p->perturbed[NUMPOINTS - 1].y << " " << p->perturbed[NUMPOINTS - 1].z << std::endl;
//	std::cout << "Velocities: " << p->perturbed[NUMPOINTS - 1].vx << " " << p->perturbed[NUMPOINTS - 1].vy << " " << p->perturbed[NUMPOINTS - 1].vz << std::endl;

	free(buffer);
}

// Copies identical data from
void copyPlanet(Planet *o, Planet i)
{
	o->m = i.m; o->x = i.x; o->y = i.y; o->z = i.z; o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}
// Copies velocity data from a probe to its perturbed probe. Perturbed probe does not have further perturbations
void copyProbeToPerturbed(Probe *o, Probe i)
{
	o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}

// Copies data from one probe to another, doesn't worry about perturbed pointer
void copyProbe(Probe *o, Probe i)
{
	o->x = i.x; o->y = i.y; o->z = i.z; o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}
