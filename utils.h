/*
 * utils.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef HEADERS_UTILS_H_
#define HEADERS_UTILS_H_

// Utility functions used by multiple main functions

// Returns array struct of perturbation vectors about a probe for a generated circle
void returnPerturbationsCircle(double **perturbations, int NUMPOINTS);
// Copies data from Planet i to Planet o
void copyPlanet(Planet *o, Planet i);
// Copies only velocity data from a probe to its perturbed probe. Perturbed probe does not have further perturbations
void copyProbeToPerturbed(Probe *o, Probe i);
// Copies data from one probe to another, doesn't worry about perturbed pointer
void copyProbe(Probe *o, Probe i);





#endif /* HEADERS_UTILS_H_ */
