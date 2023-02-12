//============================================================================
// Name        : 3dSim.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

// Custom Headers
#include "headers\sphereGenerator.h"
#include "headers\planetsAndProbe.h"
#include "headers\utils.h"
#include "headers\chaosPlot.h"
#include "headers\orbitsToFile.h"
#include "headers\singleProbeChaos.h"


//-------------------------------------------------------------------------------------------


int main(void) {

	// Generate unit sphere or circle if not already done
	bool offset = false;
	bool sphere = false;
	int NUMPOINTS = 20;
	double RADIUS = 50.0;

	if (sphere) { writeSphereToFile(NUMPOINTS, RADIUS); }
	else { writeCircleToFile(NUMPOINTS, RADIUS); }

	// Timestep
	double dt = 500;
	int totalSteps = 29000;

	// 2D Chaos Plot
	int resolution[2] = {5, 5};
	chaosPlot(dt, totalSteps, resolution, NUMPOINTS, RADIUS);
	// Debugging Single rk4 vs euler

	// Single Probe
	int ITERATIONS = 1;
	//singleProbeChaos(dt, totalSteps, offset, sphere, NUMPOINTS, RADIUS, ITERATIONS);
	return 0;
}
