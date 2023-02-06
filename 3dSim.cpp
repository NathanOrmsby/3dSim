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
void rk4(Planet *e, Planet *s, double dt);

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

	// Loop
	int c = 0;
	while (c < 100)
	{
		rk4(&e, &s, dt);
		c++;
	}
}

void calcForce(Planet e, Planet s, double *forces)
{
	double G = 6.6743e-11;

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
}

void rk4(Planet *e, Planet *s, double dt)
{

	// k1 step
	// Get net force at t0
	// Calculate the net force at the starting position
	double k1Force[6];
	calcForce(*e, *s, k1Force);


	std::cout << "Earth Pos: " << e->x << " " << e->y << " " << e->z << " " << "vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
	std::cout << "Sun Pos: " << s->x << s->y << s->z << "vel: " << s->vx << s->vy << s->vz << std::endl;

	std::cout << "K1 forces" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k1Force[i] << std::endl;
	}

	// K1 velocity
	double k1_vel[6];
	k1_vel[0] = e->vx; k1_vel[1] = e->vy; k1_vel[2] = e->vz;
	k1_vel[3] = s->vx; k1_vel[4] = s->vy; k1_vel[5] = s->vz;

	std::cout << "K1 velocities" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k1_vel[i] << std::endl;
	}

	// k2 step
	// Make new mass_list for k2 state.
	Planet k2Planets[2];

	// Copy properties of mass_list

	k2Planets[0].m = e->m; k2Planets[0].x = e->x; k2Planets[0].y = e->y; k2Planets[0].z = e->z; k2Planets[0].vx = e->vx; k2Planets[0].vy = e->vy; k2Planets[0].vz = e->vz;
	k2Planets[1].m = s->m; k2Planets[1].x = s->x; k2Planets[1].y = s->y; k2Planets[1].z = s->z; k2Planets[1].vx = s->vx; k2Planets[1].vy = s->vy; k2Planets[1].vz = s->vz;

	// Push k2 state forward dt / 2 using k1 as net force

	// Update velocities: THIS USES k1 / 2
	k2Planets[0].vx += (k1Force[0] / k2Planets[0].m) * 0.5 * dt; k2Planets[0].vy += (k1Force[1] / k2Planets[0].m) * 0.5 * dt; k2Planets[0].vz += (k1Force[2] / k2Planets[0].m) * 0.5 * dt;
	k2Planets[1].vx += (k1Force[3] / k2Planets[1].m) * 0.5 * dt; k2Planets[1].vy += (k1Force[4] / k2Planets[1].m) * 0.5 * dt; k2Planets[1].vz += (k1Force[5] / k2Planets[1].m) * 0.5 * dt;

	// Update positions: Using k1 velocity. Initial position velocity
	k2Planets[0].x += k1_vel[0] * 0.5 * dt; k2Planets[0].y += k1_vel[1] * 0.5 * dt; k2Planets[0].z += k1_vel[2] * 0.5 * dt;
	k2Planets[1].x += k1_vel[3] * 0.5 * dt; k2Planets[1].y += k1_vel[4] * 0.5 * dt; k2Planets[1].z += k1_vel[5] * 0.5 * dt;

	// Calculate k2 force
	double k2Force[6];
	calcForce(k2Planets[0], k2Planets[1], k2Force);

	std::cout << "K2 forces" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k2Force[i] << std::endl;
	}

	// K2 velocity
	double k2_vel[6];
	k2_vel[0] = k2Planets[0].vx; k2_vel[1] = k2Planets[0].vy; k2_vel[2] = k2Planets[0].vz;
	k2_vel[3] = k2Planets[1].vx; k2_vel[4] = k2Planets[1].vy; k2_vel[5] = k2Planets[1].vz;

	std::cout << "K2 velocities" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k2_vel[i] << std::endl;
	}

	// k3 step
	// Make new mass_list for k2 state.
	Planet k3Planets[2];

	// Copy properties of mass_list
	k3Planets[0].m = e->m; k3Planets[0].x = e->x; k3Planets[0].y = e->y; k3Planets[0].z = e->z; k3Planets[0].vx = e->vx; k3Planets[0].vy = e->vy; k3Planets[0].vz = e->vz;
	k3Planets[1].m = s->m; k3Planets[1].x = s->x; k3Planets[1].y = s->y; k3Planets[1].z = s->z; k3Planets[1].vx = s->vx; k3Planets[1].vy = s->vy; k3Planets[1].vz = s->vz;

	// Push k2 state forward dt / 2 using k1 as net force

	// Update velocities: THIS USES k1 / 2
	k3Planets[0].vx += (k2Force[0] / k3Planets[0].m) * 0.5 * dt; k3Planets[0].vy += (k2Force[1] / k3Planets[0].m) * 0.5 * dt; k3Planets[0].vz += (k2Force[2] / k3Planets[0].m) * 0.5 * dt;
	k3Planets[1].vx += (k2Force[3] / k3Planets[1].m) * 0.5 * dt; k3Planets[1].vy += (k2Force[4] / k3Planets[1].m) * 0.5 * dt; k3Planets[1].vz += (k2Force[5] / k3Planets[1].m) * 0.5 * dt;

	// Update positions: Using k1 velocity. Initial position velocity
	k3Planets[0].x += k2_vel[0] * 0.5 * dt; k3Planets[0].y += k2_vel[1] * 0.5 * dt; k3Planets[0].z += k2_vel[2] * 0.5 * dt;
	k3Planets[1].x += k2_vel[3] * 0.5 * dt; k3Planets[1].y += k2_vel[4] * 0.5 * dt; k3Planets[1].z += k2_vel[5] * 0.5 * dt;

	// Calculate k3 force
	double k3Force[6];
	calcForce(k3Planets[0], k3Planets[1], k3Force);

	std::cout << "K3 forces" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k3Force[i] << std::endl;
	}

	// k3 vel
	double k3_vel[6];
	k3_vel[0] = k3Planets[0].vx; k3_vel[1] = k3Planets[0].vy; k3_vel[2] = k3Planets[0].vz;
	k3_vel[3] = k3Planets[1].vx; k3_vel[4] = k3Planets[1].vy; k3_vel[5] = k3Planets[1].vz;

	std::cout << "K3 velocities" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k3_vel[i] << std::endl;
	}


	// k4 step
	// Make new mass_list for k2 state.
	Planet k4Planets[2];

	// Copy properties of mass_list
	k4Planets[0].m = e->m; k4Planets[0].x = e->x; k4Planets[0].y = e->y; k4Planets[0].z = e->z; k4Planets[0].vx = e->vx; k4Planets[0].vy = e->vy; k4Planets[0].vz = e->vz;
	k4Planets[1].m = s->m; k4Planets[1].x = s->x; k4Planets[1].y = s->y; k4Planets[1].z = s->z; k4Planets[1].vx = s->vx; k4Planets[1].vy = s->vy; k4Planets[1].vz = s->vz;

	// Push k4 state forward dt using k3 as net force

	// Update velocities: THIS USES k3
	k4Planets[0].vx += (k3Force[0] / k4Planets[0].m) * dt; k4Planets[0].vy += (k3Force[1] / k4Planets[0].m) * dt; k4Planets[0].vz += (k3Force[2] / k4Planets[0].m) * dt;
	k4Planets[1].vx += (k3Force[3] / k4Planets[1].m) * dt; k4Planets[1].vy += (k3Force[4] / k4Planets[1].m) * dt; k4Planets[1].vz += (k3Force[5] / k4Planets[1].m) * dt;

	// Update positions: Using k3 velocity. Initial position velocity
	k4Planets[0].x += k3_vel[0]  * dt; k4Planets[0].y += k3_vel[1]  * dt; k4Planets[0].z += k3_vel[2]  * dt;
	k4Planets[1].x += k3_vel[3]  * dt; k4Planets[1].y += k3_vel[4]  * dt; k4Planets[1].z += k3_vel[5]  * dt;

	// Calculate k4 force
	double k4Force[6];
	calcForce(k4Planets[0], k4Planets[1], k4Force);

	std::cout << "K4 forces" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k4Force[i] << std::endl;
	}

	// k4 vel
	double k4_vel[6];
	k4_vel[0] = k4Planets[0].vx; k4_vel[1] = k4Planets[0].vy; k4_vel[2] = k4Planets[0].vz;
	k4_vel[3] = k4Planets[1].vx; k4_vel[4] = k4Planets[1].vy; k4_vel[5] = k4Planets[1].vz;

	std::cout << "K4 velocities" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << k4_vel[i] << std::endl;
	}


	// Calculate the weighted rk4 net force (k1_force + 2*k2_force + 2*k3_force + k4_force)
	double rk4Force[6];
	double rk4_vel[6];
	for (int i = 0; i < 6; i++)
	{
		rk4Force[i] = k1Force[i] + 2 * k2Force[i] + 2 * k3Force[i] + k4Force[i];
	}
	for (int i = 0; i < 6; i++)
	{
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

	std::cout << "rK4 forces" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << rk4Force[i] << std::endl;
	}

	std::cout << "rK4 velocities" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::cout << rk4_vel[i] << std::endl;
	}

	// Calculate the weighted rk4 velocity (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

	// Push the initial state forward dt using the weighted average rk4 force as net force

	// Update velocities: THIS USES k3
	e->vx += (rk4Force[0] / e->m / 6.0) * dt; e->vy += (rk4Force[1] / e->m / 6.0) * dt; e->vz += (rk4Force[2] / e->m / 6.0) * dt;
	s->vx += (rk4Force[3] / s->m / 6.0) * dt; s->vy += (rk4Force[4] / s->m / 6.0) * dt; s->vz += (rk4Force[5] / s->m / 6.0) * dt;

	// Update positions: Using k3 velocity. Initial position velocity
	e->x += (rk4_vel[0] / 6.0) * dt; e->y += (rk4_vel[1] / 6.0) * dt; e->z += (rk4_vel[2] / 6.0) * dt;
	s->x += (rk4_vel[3] / 6.0) * dt; s->y += (rk4_vel[4] / 6.0) * dt; s->z += (rk4_vel[5] / 6.0) * dt;

	std::cout << "Earth Pos: " << e->x << " " << e->y << " " << e->z << " " << "vel: " << e->vx << " " << e->vy << " " << e->vz << std::endl;
	std::cout << "Sun Pos: " << s->x << " " << s->y << " " << s->z << " " << "vel: " << s->vx << " " << s->vy << " " << s->vz << std::endl;
}
