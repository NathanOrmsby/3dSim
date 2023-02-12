/*
 * orbitsToFile.cpp
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */


#include <iostream>
#include <math.h>
#include <fstream>
#include "..\headers\planetsAndProbe.h"
#include "..\headers\orbitsToFile.h"
#include "..\headers\utils.h"

// Main function
void orbitsToFile(double dt, int totalSteps)
{
	// Initialization: Maybe make a from file feature?
	Planet e; e.x = 0; e.y = 0; e.z = 0; e.vx = 0; e.vy = 0; e.vz = 29787.7; e.m = 5.972e24;
	Planet s; s.x = 149597870700; s.y = 0; s.z = 0; s.vx = 0; s.vy = 0; s.vz = 0; s.m = 1.989e30;
	Planet jw; jw.x = 0; jw.y = 0; jw.z = 1500000000; jw.vx = 0; jw.vy = 0; jw.vz = 28000;
	// Loop
	int c = 0;
	int stop = 50;
	// fileName
	std::string fname = "data.csv";
	// Storage buffer
	int bufSize = 9;
	double **data = (double **)malloc(stop * sizeof(double *));
	for (int i = 0; i < stop; ++i)
	{
		data[i] = (double *)malloc(bufSize * sizeof(double));
	}

	while (c < stop)
	{
		// Deposit data
		data[c][0] = e.x;
		data[c][1] = e.y;
		data[c][2] = e.z;
		data[c][3] = s.x;
		data[c][4] = s.y;
		data[c][5] = s.z;
		data[c][6] = jw.x;
		data[c][7] = jw.y;
		data[c][8] = jw.z;
		// Numerical Integration: Euler or rk4.
		rk4(&e, &s, &jw, dt);
//		euler(&e, &s, &jw, dt);
		c++;
	}
	orbitsToFile(fname, data, stop);
}


// Implementations of functions utilized by orbitsToFile
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

// Euler for only planets
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

// Writes planetary orbits to file
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

