/*
 * orbitsToFile.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef HEADERS_ORBITSTOFILE_H_
#define HEADERS_ORBITSTOFILE_H_

#include "planetsAndProbe.h"

// Main function: Outputs orbital path data to csv
void orbitsToFile(double dt, int totalSteps);

// Functions utilized by orbitsToFile

// Calculates gravitational forces between only planets
void calcForce(Planet e, Planet s, Planet jw, double *forces);
// rk4 for only planets
void rk4(Planet *e, Planet *s, Planet *jw, double dt);
// Euler for only planets
void euler(Planet *e, Planet *s, Planet *jw, double dt);
// Initializes the perturbations about a probe
void initializePerturbations(Probe *p, bool sphere);
// Writes planetary orbits to file
bool orbitsToFile(std::string fileName, double **data, int dataLen);


#endif /* HEADERS_ORBITSTOFILE_H_ */
