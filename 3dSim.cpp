//============================================================================
// Name        : 3dSim.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "sphereGenerator.h"

#define bufSize 9
#define NUMPOINTS 360
// Plotting data per ITERATIONS timesteps
#define ITERATIONS 1

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

void calcForce(Planet e, Planet s, double **forces);
void rk4(Planet *e, Planet *s, Planet *jw, double dt);
void euler(Planet *e, Planet *s, Planet *jw, double dt);
void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt);
void copyPlanet(Planet a, Planet *b);
bool orbitsToFile(std::string fileName, double **data, int dataLen);
void perturbationsToFile(double **data, int dataLen, bool offset);
void initializePerturbations(Probe *p, bool sphere);
void copyProbe(Probe *o, Probe i);
void singleProbeChaos(double dt, int totalSteps, bool offset, bool sphere);

int main() {

	// Generate unit sphere or circle if not already done
	writeSphereToFile(NUMPOINTS, 5000.0);
	//writeCircleToFile(NUMPOINTS, 500);

	// Timestep
	double dt = 1;
	int totalSteps = 3;

	// Debugging
	bool offset = true;
	bool sphere = true;
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

	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{

		// Store data of JW and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			data[dc][0] = jw.x;
			data[dc][1] = jw.y;
			data[dc][2] = jw.z;
			for (int i = 0; i < NUMPOINTS; ++i)
			{
				data[dc][3 * i + 3] = jw.perturbed[i].x;
				data[dc][3 * i + 4] = jw.perturbed[i].y;
				data[dc][3 * i + 5] = jw.perturbed[i].z;
			}
			dc++;
		}

		// Do euler and loop
		chaosEuler(&e, &s, &jw, dt);
		c++;
	}

	// Write data to file
	perturbationsToFile(data, n, offset);

	// Free
	for (int i = 0; i < n; ++i)
	{
		free(data[i]);
	}
	free(data);

}
//-------------------------------------------------------------------------------------------------------

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
	for (int i = 0; i < (NUMPOINTS / 3); ++i)
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

void rk4(Planet *e, Planet *s, Planet *jw, double dt)
{

	// k1 step
	// Get net force at t0
	// Calculate the net force at the starting position
	double k1Force[9];
	calcForce(*e, *s, *jw, k1Force);


//	std::cout << "Earth Pos: " << e->x << " " << e->y << " " << e->z << " " << "vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
//	std::cout << "Sun Pos: " << s->x << s->y << s->z << "vel: " << s->vx << s->vy << s->vz << std::endl;
//
//	std::cout << "K1 forces" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k1Force[i] << std::endl;
//	}

	// K1 velocity
	double k1_vel[6];
	k1_vel[0] = e->vx; k1_vel[1] = e->vy; k1_vel[2] = e->vz;
	k1_vel[3] = s->vx; k1_vel[4] = s->vy; k1_vel[5] = s->vz;
	k1_vel[6] = jw->vx; k1_vel[7] = jw->vy; k1_vel[8] = jw->vz;

//	std::cout << "K1 velocities" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k1_vel[i] << std::endl;
//	}

	// k2 step
	// Make new mass_list for k2 state.
	Planet k2Planets[3];

	// Copy properties of mass_list
	copyPlanet(*e, &k2Planets[0]);
	copyPlanet(*s, &k2Planets[1]);
	copyPlanet(*jw, &k2Planets[2]);

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

//	std::cout << "K2 forces" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k2Force[i] << std::endl;
//	}

	// K2 velocity
	double k2_vel[9];
	k2_vel[0] = k2Planets[0].vx; k2_vel[1] = k2Planets[0].vy; k2_vel[2] = k2Planets[0].vz;
	k2_vel[3] = k2Planets[1].vx; k2_vel[4] = k2Planets[1].vy; k2_vel[5] = k2Planets[1].vz;
	k2_vel[6] = k2Planets[2].vx; k2_vel[7] = k2Planets[2].vy; k2_vel[8] = k2Planets[2].vz;

//	std::cout << "K2 velocities" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k2_vel[i] << std::endl;
//	}

	// k3 step
	// Make new mass_list for k2 state.
	Planet k3Planets[3];

	// Copy properties of mass_list
	copyPlanet(*e, &k3Planets[0]);
	copyPlanet(*s, &k3Planets[1]);
	copyPlanet(*jw, &k3Planets[2]);

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

//	std::cout << "K3 forces" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k3Force[i] << std::endl;
//	}

	// k3 vel
	double k3_vel[9];
	k3_vel[0] = k3Planets[0].vx; k3_vel[1] = k3Planets[0].vy; k3_vel[2] = k3Planets[0].vz;
	k3_vel[3] = k3Planets[1].vx; k3_vel[4] = k3Planets[1].vy; k3_vel[5] = k3Planets[1].vz;
	k3_vel[6] = k3Planets[2].vx; k3_vel[7] = k3Planets[2].vy; k3_vel[8] = k3Planets[2].vz;

//	std::cout << "K3 velocities" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k3_vel[i] << std::endl;
//	}


	// k4 step
	// Make new mass_list for k2 state.
	Planet k4Planets[3];

	// Copy properties of mass_list
	copyPlanet(*e, &k4Planets[0]);
	copyPlanet(*s, &k4Planets[1]);
	copyPlanet(*jw, &k4Planets[2]);

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

//	std::cout << "K4 forces" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k4Force[i] << std::endl;
//	}

	// k4 vel
	double k4_vel[6];
	k4_vel[0] = k4Planets[0].vx; k4_vel[1] = k4Planets[0].vy; k4_vel[2] = k4Planets[0].vz;
	k4_vel[3] = k4Planets[1].vx; k4_vel[4] = k4Planets[1].vy; k4_vel[5] = k4Planets[1].vz;
	k4_vel[6] = k4Planets[2].vx; k4_vel[7] = k4Planets[2].vy; k4_vel[8] = k4Planets[2].vz;


//	std::cout << "K4 velocities" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << k4_vel[i] << std::endl;
//	}


	// Calculate the weighted rk4 net force (k1_force + 2*k2_force + 2*k3_force + k4_force)
	double rk4Force[9];
	double rk4_vel[9];
	for (int i = 0; i < 9; ++i)
	{
		rk4Force[i] = k1Force[i] + 2 * k2Force[i] + 2 * k3Force[i] + k4Force[i];
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

//	std::cout << "rK4 forces" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << rk4Force[i] << std::endl;
//	}
//
//	std::cout << "rK4 velocities" << std::endl;
//	for (int i = 0; i < 6; i++)
//	{
//		std::cout << rk4_vel[i] << std::endl;
//	}

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

//	std::cout << "Earth Pos: " << e->x << " " << e->y << " " << e->z << " " << "vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
//	std::cout << "Sun Pos: " << s->x << " " << s->y << " " << s->z << " " << "vel: " << s->vx << " " << s->vy << " " << s->vz << std::endl;
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

void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt)
{
	double *accels = (double *)malloc(6 + (3 * NUMPOINTS) * sizeof(double));
	chaosCalcForce(*e, *s, *jw, accels);

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

void copyPlanet(Planet a, Planet *b)
{
	b->m = a.m; b->x = a.x; b->y = a.y; b->z = a.z; b->vx = a.vx; b->vy = a.vy; b->vz = a.vz;
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


// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid
void initializePerturbations(Probe *p, bool sphere)
{
	// Initialize all perturbation probes
	p->perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
	for (int i = 0; i < NUMPOINTS; ++i) { copyProbe(&p->perturbed[i], *p); }

	// Read from file
	FILE *f;
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
	std::cout << "Printing JW telescope location" << std::endl;
	std::cout << p->x << " " << p->y << " " << p->z << std::endl;

	std::cout << "Printing perturbations" << std::endl;
	for (int i = 0; i < NUMPOINTS; ++i)
	{
		std::cout << p->perturbed[i].x << " " << p->perturbed[i].y << " " << p->perturbed[i].z << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Velocities: " << p->perturbed[0].vx << " " << p->perturbed[0].vy << " " << p->perturbed[0].vz << std::endl;
	std::cout << "Last" << std::endl;
	std::cout << p->perturbed[NUMPOINTS - 1].x << " " << p->perturbed[NUMPOINTS - 1].y << " " << p->perturbed[NUMPOINTS - 1].z << std::endl;
	std::cout << "Velocities: " << p->perturbed[NUMPOINTS - 1].vx << " " << p->perturbed[NUMPOINTS - 1].vy << " " << p->perturbed[NUMPOINTS - 1].vz << std::endl;
	free(buffer);
}

// Copies identical data from a probe to its perturbed probe. Perturbed probe does not have further perturbations
void copyProbe(Probe *o, Probe i)
{
	o->vx = i.vx;
	o->vy = i.vy;
	o->vz = i.vz;
	o->perturbed = nullptr;
}
