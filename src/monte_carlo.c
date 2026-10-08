#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
extern int nnsize;

int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double ediff;
  double newpos[D], oldpos[D], dist_old[D], dist_new[D];

#ifdef FAST
  check_neighbor_list(par, pos);
#endif
  // Loop over the particles
  for (i = 0; i < par->n; i++) {
    // Suggest a new position: r'_i = r_i + b * delta
    // where delta is a random vector
    for (d = 0; d < D; d++) {
      oldpos[d] = pos[D * i + d];
      newpos[d] = oldpos[d] + par->b * dran_sign();
      // ensure new position is within the simulation box
      newpos[d] = inrange(newpos[d], par->L[d]);
    }

    // Calculate energy difference Delta_U = U_newpos - U_oldpos
    // Computing the total potential energy would take O(N^2) operations
    // the change in energy is calculated by considering only the interactions
    // between the moved particle i and all other particles j
    // Delta_U = sum_{j != i} [u(r'_i - r_j) - u(r_i - r_j)]
    // Delta_U = energy after the move - energy before the move
    ediff = 0.0;
    for (j = 0; j < par->n; j++) {
      if (j == i) continue;
      // When a particle goes off the right edge, it reappears on the left edge
      // if L=10, particle A is at x = 0.1 and particle B is at x = 9.9, the distance is -9.8
      // but since the box wraps around, B is just across the wall from A, only 0.2 units apart
      // hence we need to wrap the distance component into the range [-L/2, L/2] to measure the distance the short way round

      double r2_old = 0.0, r2_new = 0.0;
      for (d = 0; d < D; d++) {
        // interaction energy between particle i and all other particles j before i moves
        double d_old = oldpos[d] - pos[D * j + d];
        // shift the coordinate difference so it is always wrapped into the range [-L/2, L/2]
        d_old -= par->L[d] * round(d_old / par->L[d]); // Minimum image convention
        r2_old += d_old * d_old;
        // interaction energy between particle i and all other particles j at its proposed new position
        double d_new = newpos[d] - pos[D * j + d];
        d_new -= par->L[d] * round(d_new / par->L[d]);
        r2_new += d_new * d_new;
      }

      // the Lennard-Jones interaction: u(r) = 4 * eps * ((sigma/r)^12 - (sigma/r)^6))
      // assume eps = sigma = 1
      if (r2_old < CUT * CUT) {
        double r6_inv = 1.0 / (r2_old * r2_old * r2_old);
        double u_old = 4.0 * (r6_inv * r6_inv - r6_inv);
        ediff -= u_old;
      }
      if (r2_new < CUT * CUT) {
        double r6_inv = 1.0 / (r2_new * r2_new * r2_new);
        double u_new = 4.0 * (r6_inv * r6_inv - r6_inv);
        ediff += u_new;
      }
    }

    // alpha = min(1, e^(-Delta_U/T))
    double alpha = 1.0;
    if (ediff > 0.0) {
      alpha = exp(-ediff / par->T);
    }
    // Metropolis acceptance criterion
    // Generate random number xi in [0, 1) and accept if xi < alpha
    if (dran() < alpha) {
      for (d = 0; d < D; d++) {
        // accept the move: update the position of particle i
        pos[D * i + d] = newpos[d];
      }
      naccept++;
    }
  }
						 
  return naccept;
}