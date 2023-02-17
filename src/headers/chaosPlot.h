/*
 * chaosPlot.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef HEADERS_CHAOSPLOT_H_
#define HEADERS_CHAOSPLOT_H_

#include "planetsAndProbe.h"
// Main function
// Outputs lyapunov exponents to a csv for plotting a 2d Chaos Map of the Earth Sun Moon JW Telescope system
void chaosPlot(double dt, int totalSteps, int resolution[2], int NUMPOINTS, double RADIUS);

// Returns an array of x, z probe positions distributed across L1 to L2 for all probes
void partitionPositions(double *positions, int resolution[2]);

// Submain functions utilized by main function

// Extracts and returns array of min and max lyapunov exponents for a number of probes at various positions. Renormalizes at each timestep
// Orbital Separation method for 2d chaos plot
double *orbitalSeparationMultiProbe(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int NUMPOINTS, double RADIUS);

// Extracts and returns array of min and max lyapunov exponents for a number of probes at various positions
// NO RENORMALIZATION
double *multiProbeChaos(double dt, int totalSteps, double *probePositions, int probeNum, double **perturbations, int NUMPOINTS, double RADIUS);


// Functions utilized by submain functions

// Classic Runge Kutta for multiple probes with perturbations
void multiChaosRK4(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int NUMPOINTS);
// Calculates accelerations. Handles multiple probes with perturbations
void multiChaosCalcAccel(Planet e, Planet s, Probe *probes, int probeNum, double **accels, int NUMPOINTS);
// Moves planets and probes one step forward. Helper function to multiChaosRK4
void rk4StepMultiChaos(Planet *eCopy, Probe *probeCopies, Planet e, Probe *probes, int probeNum, double **accels, double dt, int NUMPOINTS);
// Euler method for multiple probes and perturbations
void multiChaosEuler(Planet *e, Planet *s, Probe *probes, int probeNum, double dt, int NUMPOINTS);
// Not actually Gram Schmidt: Needed for calculation of minimum lyapunov. Returns difference in area between evolved parallelogram and original rectange. Calculates for single probe
double gramSchmidt(Probe *probes, double initialArea, int i, int NUMPOINTS, int maxInd);
// Gram Schmidt: Needed for calculation of minimum lyapunov. Calculates difference in area and renormalizes to RADIUS orthogonal vectors that preserve direction
void gramSchmidtRenormalize(Probe *probes, double *areaSums, int probeNum, double initialArea, int NUMPOINTS);
// Extracts lyapunov exponents from multiple probes and deposits them into an array: DOES NOT CALCULATE MIN USING GRAM SCHMIDT RN: FIX THIS
void extractMultiLyapunov(Probe *probes, int probeNum, double initialArea, double *lyapunovs, double t, int NUMPOINTS, double RADIUS);
// Extracts lyapunov exponents for multiple probes using orbit separation method
void extractOrbitalMultiLyapunov(Probe *probes, int probeNum, double *lyapunovs, double **LfSums, double *areaSums, double t, int NUMPOINTS);
// Writes to csv file lyapunov exponent mins, maxes, and maxes - mins.
void multiLyapunovToFile(double *data, int dataLen);
// Outputs normalized initial and final gram schmidt vectors as well as the vector in the direction of maximum exponential divergence for specific probe to csv for plotting. Called by gram schmidt function.
void gramSchmidtToFile(double u1f[2], double u2f[2], double umax[2]);

#endif /* HEADERS_CHAOSPLOT_H_ */
