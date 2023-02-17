/*
 * utils.cpp
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#include <iostream>
#include <math.h>
#include "..\headers\planetsAndProbe.h"
#include "..\headers\utils.h"

// Utility functions utilized by multiple main function implementations

// Returns an array structure of circular perturbations in xz plane read from file. Faster to read from array for use over multiple probes. Struct: data[NUMPOINTS][x, z]
void returnPerturbationsCircle(double **perturbations, int NUMPOINTS)
{
	// Read from file
	FILE *f = fopen("unitCircle.csv", "r");
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);

	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);

	// Read from buffer into struct
	int pointCount = 0;
	char arr[30];
	int arrCount = 0;
	int commaCount = 0;
	int i = 0;

	while (pointCount < NUMPOINTS)
	{
		if (buffer[i] == '\n')
		{
			arr[arrCount] = '\0';
			perturbations[pointCount][1] = std::atof(arr);
			pointCount++;
			arrCount = 0;
			commaCount = 0;

		}
		else if (buffer[i] == ',')
		{
			arr[arrCount] = '\0';
			if (commaCount == 0)
			{
				perturbations[pointCount][0] = std::atof(arr);
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

	free(buffer);
}

// Copies identical data from
void copyPlanet(Planet *o, Planet i)
{
	o->m = i.m; o->x = i.x; o->y = i.y; o->z = i.z; o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}
// Copies velocity data from a probe to its perturbed probe. Perturbed probe does not have further perturbations
void copyProbeToPerturbed(Probe *o, Probe i)
{
	o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}

// Copies data from one probe to another, doesn't worry about perturbed pointer
void copyProbe(Probe *o, Probe i)
{
	o->x = i.x; o->y = i.y; o->z = i.z; o->vx = i.vx; o->vy = i.vy; o->vz = i.vz;
}




