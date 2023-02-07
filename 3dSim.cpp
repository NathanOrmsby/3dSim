//============================================================================
// Name        : 3dSim.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <math.h>

typedef struct
{
	double x, y, z, vx, vy, vz;
	double m;
} Planet;

//-------------------------------------------------------------------------------------------

void calcForce(Planet e, Planet s, double *forces);
void rk4(Planet *e, Planet *s, Planet *jw, double dt);
void copyPlanet(Planet a, Planet *b);

int main() {

	// Timestep
	double dt = 1 / 1000;
	dt = 50;

	Planet e;
	e.x = 0; e.y = 0; e.z = 0;
	e.vx = 0; e.vy = 0; e.vz = 29.7877e3;
	e.m = 5.972e24;
	Planet s;
	s.x = 149597870700; s.x = 147120163000; s.y = 0; s.z = 0;
	s.vx = 0; s.vy = 0; s.vz = 0;
	s.m = 1.989e30;
	Planet jw;
	jw.x = 1850000; jw.y = 0; jw.z = 0;
	jw.vx = 0; jw.vy = 0; jw.vz = 30059.3;
	jw.m = 6200;

	// Loop
	int c = 0;
	while (c < 100)
	{
		rk4(&e, &s, &jw, dt);
		c++;
	}
	std::cout << "Earth Pos: " << e.x << " " << e.y << " " << e.z << " " << "vel: " << e.vx << " " << e.vy << " " << e.vz << std::endl;
	std::cout << "Sun Pos: " << s.x << " " << s.y << " " << s.z << " " << "vel: " << s.vx << " " << s.vy << " " << s.vz << std::endl;
	std::cout << "Sun Pos: " << s.x << " " << s.y << " " << s.z << " " << "vel: " << s.vx << " " << s.vy << " " << s.vz << std::endl;
}

//-------------------------------------------------------------------------------------------------------

void calcForce(Planet e, Planet s, Planet jw, double *forces)
{
	double G = 6.6743e-11;

	// Earth Sun
	double r12 = sqrt((s.x - e.x) * (s.x - e.x) + (s.y - e.y) * (s.y - e.y) + (s.z - e.z) * (s.z - e.z));
	double r21 = sqrt((e.x - s.x) * (e.x - s.x) + (e.y - s.y) * (e.y - s.y) + (e.z - s.z) * (e.z - s.z));
	double n = ((G * e.m * s.m) / (r12 * r12 * r12));
	forces[0] = n * (s.x - e.x);
	forces[1] = n * (s.y - e.y);
	forces[2] = n * (s.z - e.z);
	n = ((G * e.m * s.m) / (r21 * r21 * r21));
	forces[3] = n * (e.x - s.x);
	forces[4] = n * (e.y - s.y);
	forces[5] = n * (e.z - s.z);

	// James Webb
	double r13 = sqrt((jw.x - s.x) * (jw.x - s.x) + (jw.y - s.y) * (jw.y - s.y) + (jw.z - s.z) * (jw.z - s.z));
	double r23 = sqrt((jw.x - e.x) * (jw.x - e.x) + (jw.y - e.y) * (jw.y - e.y) + (jw.z - e.z) * (jw.z - e.z));
	double n13 = ((G * s.m * jw.m) / (r13 * r13 * r13));
	double n23 = ((G * e.m * jw.m) / (r23 * r23 * r23));
	forces[6] = (n13 * (jw.x - s.x)) + (n23 * (jw.x - e.x));
	forces[7] = (n13 * (jw.y - s.y)) + (n23 * (jw.y - e.y));
	forces[8] = (n13 * (jw.z - s.z)) + (n23 * (jw.z - e.z));
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

void copyPlanet(Planet a, Planet *b)
{
	b->m = a.m; b->x = a.x; b->y = a.y; b->z = a.z; b->vx = a.vx; b->vy = a.vy; b->vz = a.vz;
}
