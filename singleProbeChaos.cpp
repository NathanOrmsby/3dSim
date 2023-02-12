/*
 * singleProbeChaos.cpp
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "..\headers\singleProbeChaos.h"
#include "..\headers\planetsAndProbe.h"
#include "..\headers\utils.h"

// Main function singleProbeChaos: Plots an evolution of a single probe and its perturbations (spherical or circular xz) over a specified time
void singleProbeChaos(double dt, int totalSteps, bool offset, bool sphere, int NUMPOINTS, double RADIUS, int ITERATIONS)
{
	// Initialize Planets and Probe
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Probe jw; jw.x = 0; jw.y = 0; jw.z = 1.5e9; jw.vx = 0; jw.vy = 0; jw.vz = 30173.943; jw.perturbed = nullptr;
	initializePerturbations(&jw, sphere, NUMPOINTS);
	// Store data: Store all data of every ITERATIONS step
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
	double **data = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i) { data[i] = (double *)malloc((3 + 3 * NUMPOINTS) * sizeof(double)); }
	// Loop
	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{
		// Store data of JW and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			data[dc][0] = jw.x; data[dc][1] = jw.y; data[dc][2] = jw.z;
			for (int i = 0; i < NUMPOINTS; ++i)
			{
				data[dc][3 * i + 3] = jw.perturbed[i].x; data[dc][3 * i + 4] = jw.perturbed[i].y; data[dc][3 * i + 5] = jw.perturbed[i].z;
			}
			dc++;
		}
		// Numerical Integrator: Euler or RK4
		chaosEuler(&e, &s, &jw, dt, NUMPOINTS);
		// chaosRK4(&e, &s, &jw, dt, NUMPOINTS);
		c++;
	}
	// Store latest info
	data[n - 1][0] = jw.x; data[n - 1][1] = jw.y; data[n - 1][2] = jw.z;
	for (int i = 0; i < NUMPOINTS; ++i) { data[n - 1][3 * i + 3] = jw.perturbed[i].x; data[n - 1][3 * i + 4] = jw.perturbed[i].y; data[n - 1][3 * i + 5] = jw.perturbed[i].z; }
	// Extract lyapunov exponents
	double lyapunov[2];
	extractSingleLyapunov(data, n, lyapunov, sphere, dt * totalSteps, NUMPOINTS, RADIUS);
	// Write data to file
	perturbationsToFile(data, n, offset, NUMPOINTS);
	// Free the data arr
	for (int i = 0; i < n; ++i) { free(data[i]); }
	free(data);
}
// Implementations for functions used by main function singleProbeChaos

// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid
void initializePerturbations(Probe *p, bool sphere, int NUMPOINTS)
{
	// Initialize all perturbation probes
	p->perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
	for (int i = 0; i < NUMPOINTS; ++i) { copyProbeToPerturbed(&p->perturbed[i], *p); }

	// Read from file
	FILE *f;
	// Desktop path
	if (sphere) { f = fopen("D:\\3dsim\\unitSphere.csv", "r"); }
	else { f = fopen("D:\\3dsim\\unitCircle.csv", "r"); }
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);
	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);
	// Read from buffer into perturbed probes
	int pointCount = 0;
	char arr[30];
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
	free(buffer);
}

// RK4 for a single probe with perturbations
// Maybe debugged?
void chaosRK4(Planet *e, Planet *s, Probe *jw, double dt, int NUMPOINTS)
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
	for (int i = 0; i < NUMPOINTS; ++i) { copyProbe(&k2Probe.perturbed[i], jw->perturbed[i]); copyProbe(&k3Probe.perturbed[i], jw->perturbed[i]); copyProbe(&k4Probe.perturbed[i], jw->perturbed[i]); }
	// Initialize KAccels
	double *k1Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k2Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k3Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double)); double *k4Accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double));
	// k1 step
	// Put accelerations in array
	chaosCalcForce(*e, *s, *jw, k1Accels, NUMPOINTS);

	// k2 Step
	// Move k2 object copies dt / 2 using k1 Accelerations and k1 Velocities
	rk4StepSingleChaos(&k2Planets[0], &k2Probe, *e, *jw, k1Accels, dt / 2.0, NUMPOINTS);
	// Put data in arrays
	chaosCalcForce(k2Planets[0], k2Planets[1], k2Probe, k2Accels, NUMPOINTS);

	// k3 Step
	// Move k3 object copies dt / 2 using k2 Accelerations and k2 Velocities
	rk4StepSingleChaos(&k3Planets[0], &k3Probe, k2Planets[0], k2Probe, k2Accels, dt / 2.0, NUMPOINTS);
	// Put data in arrays
	chaosCalcForce(k3Planets[0], k3Planets[1], k3Probe, k3Accels, NUMPOINTS);

	// k4 Step
	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
	rk4StepSingleChaos(&k4Planets[0], &k4Probe, k3Planets[0], k3Probe, k3Accels, dt, NUMPOINTS);
	// Put data in arrays
	chaosCalcForce(k4Planets[0], k4Planets[1], k4Probe, k4Accels, NUMPOINTS);

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

	// BUG: Division isnt commutative
	e->vx += (k1Accels[0] + 2 * k2Accels[0] + 2 * k3Accels[0] + k4Accels[0]) * dt / 6.0; e->vy += (k1Accels[1] + 2 * k2Accels[1] + 2 * k3Accels[1] + k4Accels[1]) * dt / 6.0; e->vz += (k1Accels[2] + 2 * k2Accels[2] + 2 * k3Accels[2] + k4Accels[2]) * dt / 6.0;
	e->x += (e->vx + 2 * k2Planets[0].vx + 2 * k3Planets[0].vx + k4Planets[0].vx) * dt / 6.0; e->y += (e->vy + 2 * k2Planets[0].vy + 2 * k3Planets[0].vy + k4Planets[0].vy) * dt / 6.0; e->z += (e->vz + 2 * k2Planets[0].vz + 2 * k3Planets[0].vz + k4Planets[0].vz) * dt / 6.0;
	jw->vx += (k1Accels[3] + 2 * k2Accels[3] + 2 * k3Accels[3] + k4Accels[3]) * dt / 6.0; jw->vy += (k1Accels[4] + 2 * k2Accels[4] + 2 * k3Accels[4] + k4Accels[4]) * dt / 6.0; jw->vz += (k1Accels[5] + 2 * k2Accels[5] + 2 * k3Accels[5] + k4Accels[5]) * dt / 6.0;
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

// Euler method for single probe with perturbations
void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt, int NUMPOINTS)
{
	// Assign memory to acceleration array for Euler
	double *accels = (double *)malloc((6 + 3 * NUMPOINTS) * sizeof(double));
	chaosCalcForce(*e, *s, *jw, accels, NUMPOINTS);

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

// Calculates accelerations for earth, sun, single james webb, and all its perturbations
void chaosCalcForce(Planet e, Planet s, Probe jw, double *accels, int NUMPOINTS)
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

	// DEBUGGING
//	std::cout << "Single chaos calc force" << std::endl;
//	std::cout << "Sun on earth gravity" << std::endl;
//	std::cout << gravity << "Vectors x: " << accels[0] << " y: " << accels[1] << " z: " << accels[2] << std::endl;
//	std::cout << "Sun on probe gravity" << std::endl;
//	std::cout << gravity31 << "Vectors x: " << gravity31 * r31Unit[0] << " y: " << gravity31 * r31Unit[1] << " z: " << gravity31 * r31Unit[2] << std::endl;
//	std::cout << "Earth on probe gravity" << std::endl;
//	std::cout << gravity32 << "Vectors x: " << gravity32 * r32Unit[0] << " y: " << gravity32 * r32Unit[1] << " z: " << gravity32 * r32Unit[2] << std::endl;
//	std::cout << "Earth probe vectors" << std::endl;
//	std::cout << "Vector r32: x:" << r32[0] << " y: " << r32[1] << " z: " << r32[2] << std::endl;
//	std::cout << "Vector Magnitude: " << magr32 << std::endl;
//	std::cout << "Unit vector x: " << r32Unit[0] << " y: " << r32Unit[1] << " z: " << r32Unit[2] << std::endl;
//	std::cout << "Net gravity on probe" << std::endl;
//	std::cout << "Vectors x: " << accels[3] << " y: " << accels[4] << " z: " << accels[5] << std::endl;

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

// Moves objects forward using kAccels and kVels. Step for rk4
// PRETTY SURE THIS IS FINE
void rk4StepSingleChaos(Planet *eCopy, Probe *jwCopy, Planet e, Probe jw, double *accels, double dt, int NUMPOINTS)
{
	// DEBUGGING
//	std::cout << "Before" << std::endl;
//	std::cout << "JW pos: " << jwCopy->x << " " << jwCopy->y << " " << jwCopy->z << " vel: " << jwCopy->vx << " " << jwCopy->vy << " " << jwCopy->vz << std::endl;
	eCopy->vx += accels[0] * dt; eCopy->vy += accels[1] * dt; eCopy->vz += accels[2] * dt;
	eCopy->x += e.vx * dt; eCopy->y += e.vy * dt; eCopy->z += e.vz * dt;
	jwCopy->vx += accels[3] * dt; jwCopy->vy += accels[4] * dt; jwCopy->vz += accels[5] * dt;
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

// Extracts two lyapunov exponents from a single probe with perturbations
void extractSingleLyapunov(double **data, int dataLen, double lyapunov[2], bool sphere, double t, int NUMPOINTS, double RADIUS)
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
			std::cout << "Probe perturbation " << i << " is at pos: x: " << data[dataLen -1][3 * i + 3] << " z: " << data[dataLen -1][3 * i + 5] << std::endl;
			std::cout << "Distance to probe is: " << std::setprecision(5) << dist[i] / RADIUS << std::endl;
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

// Writes data points for evolution of single James Webb with perturbations to file.
void perturbationsToFile(double **data, int dataLen, bool offset, int NUMPOINTS)
{
	std::string fileName1 = "..\\singleChaosData\\probe.csv";
	std::string fileName2 = "..\\singleChaosData\\perturbations.csv";
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





