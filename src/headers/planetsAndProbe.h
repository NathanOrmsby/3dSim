/*
 * planetsAndProbe.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef PLANETSANDPROBE_H_
#define PLANETSANDPROBE_H_

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



#endif /* PLANETSANDPROBE_H_ */
