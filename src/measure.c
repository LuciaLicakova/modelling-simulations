#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
#ifdef FAST
extern int nnsize;
#endif

// Measure potential energy, and pressure and perhaps kinetic energy
void measure(Par *par, double *atoms, Measured *val)
{
  int i, j, d;
  double dist[D];	// Distance between particles i and j.
  double f[D];		// Force between particles
  double r2;		// Distance square between particles
  double *pos = atoms;
#ifdef VEL
  double *vel = atoms + D * par->n;
#endif

  double epot = 0;
  double virial = 0.0;

  // Potential energy
  for (i = 0; i < par->n; i++) {
    for (j = i + 1; j < par->n; j++) {
      // Calculate distance squared between particles i and j
      // Fixed this (1)  r2 = dist2( ... );	// The function dist2 is in "common.c"
      // it takes L (the size of the simulation cell), p1 and p2 (particle positions)
      // If the value describes how the simulation is set up, it belongs in par. If it describes what the individual particles are doing at this exact moment (positions, velocities), it belongs in data arrays like pos or vel.
      // pos points to the start of the array (particle 0)
      // each particle has D coordinates, so to reach particle i, we need to shift by D*i positions
      r2 = dist2(par->L, pos + D*i, pos + D*j, dist);
      if (r2 < CUT * CUT) {
	epot += pair_interaction(r2);
	virial += sqrt(r2) * force_magnitude(r2);
      }
    }
  }

  epot /= par->n;	// We want potential energy per particle
  val->epot = epot;	// Return information in the struct 'val'
  double etot = epot;	// Total energy (ifdef VEL: add ekin)

  // Calculate pressure from static variables
  // pressure = (N * T + virial / Dimensionality) / volume
  double pressure = (par->n * par->T + virial / D) / par->vol;
  // Use the struct val to return information
  val->pressure = pressure;
  
#ifdef VEL
  // Kinetic energy. 
  double ekin = 0.0;
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      ekin += par->mass * vel[D * i + d] * vel[D * i + d] / 2.0;

  ekin /= par->n;	// We want kinetic energy per particle
  val->ekin = ekin;
  etot += ekin;
#endif

  val->etot = etot;
  // include energy squared
  val->etot2 = etot * etot;

  return;
}
