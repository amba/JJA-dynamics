#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <error.h>
#include <errno.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#define N_x  200
#define N_y  200
#define SCALAR_TYPE float
#define SINE sinf

#define EPS  0.01
static SCALAR_TYPE phases[N_x*N_y]; // phi(i,j) = phases(i + N_y * j)

// static fields are all initialized to zero
// static SCALAR_TYPE currents_x[(N_x-1)*N_y];
// static SCALAR_TYPE *currents_y[N_x*(N_y-1)];


static int run_steps(int num_steps) {
  for (int n = 0; n < num_steps; ++n) {
    for (int i = 0; i < N_x; ++i) {
      for (int j = 0; j < N_y; ++j) {
        SCALAR_TYPE div_I_super = 0;
        SCALAR_TYPE phase = phases[i + N_y * j];
        if (i < N_x-1) {
          div_I_super += SINE(phases[i+1 + N_y*j] - phase);
        }
        if (i > 0) {
          div_I_super += SINE(phases[i-1 + N_y*j] - phase);
        }
        if (j < N_y - 1) {
          div_I_super += SINE(phases[i + N_y*(j+1)] - phase);
        }
        if (j > 0) {
          div_I_super += SINE(phases[i + N_y*(j-1)] - phase);
        }
        phases[i + N_y *j] -= EPS * div_I_super;
      }
    }
  }
}

int
main (int argc, char **argv)
{
  printf("Running ground state annealer\n");
  printf("N_x = %d, N_y = %d\n", N_x, N_y);
  run_steps(100000);
  return 0;
}
