/*
 * singleProbeChaos.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEPROBECHAOS_H_
#define HEADERS_SINGLEPROBECHAOS_H_

#include "planetsAndProbe.h"

// Main function: Extracts lyapunov exponents for a single probe with perturbations. Also outputs data to csv for plottong of the evolution of the spherical or circular perturbations through time.
void singleProbeChaos(double dt, int totalSteps, bool offset, bool sphere, int NUMPOINTS, double RADIUS, int ITERATIONS);
// Functions utilized by main function singleProbeChaos

// Initializes the perturbations about a probe reading directly from file
void initializePerturbations(Probe *p, bool sphere, int NUMPOINTS);
// RK4 for a single probe with perturbations
// Maybe debugged?
void chaosRK4(Planet *e, Planet *s, Probe *jw, double dt, int NUMPOINTS);
// Euler for a single probe with perturbations
void chaosEuler(Planet *e, Planet *s, Probe *jw, double dt, int NUMPOINTS);
// Calculates accelerations for earth, sun, single james webb, and all its perturbations
void chaosCalcForce(Planet e, Planet s, Probe jw, double *accels, int NUMPOINTS);
// Moves objects forward using kAccels and kVels. Step for rk4
// PRETTY SURE THIS IS FINE
void rk4StepSingleChaos(Planet *eCopy, Probe *jwCopy, Planet e, Probe jw, double *accels, double dt, int NUMPOINTS);
// Extracts two lyapunov exponents from a single probe with perturbations
void extractSingleLyapunov(double **data, int dataLen, double lyapunov[2], bool sphere, double t, int NUMPOINTS, double RADIUS);
// Writes data points for evolution of single James Webb with perturbations to file.
void perturbationsToFile(double **data, int dataLen, bool offset, int NUMPOINTS);
#endif /* HEADERS_SINGLEPROBECHAOS_H_ */
