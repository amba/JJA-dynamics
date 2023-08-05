#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#define N_x  100
#define N_y  100
#define SCALAR_TYPE float
#define SINE sinf
#define COS cosf
#define PI 3.141592

static const SCALAR_TYPE EPS  = 0.4;
static SCALAR_TYPE frustration = 0.25* 2 * PI;


static SCALAR_TYPE phases[N_x*N_y]; // phi(i,j) = phases(i + N_y * j)

// ^
// |
//
// j (row)
//    i (column)   - >


// static fields are all initialized to zero
// static SCALAR_TYPE currents_x[(N_x-1)*N_y];
// static SCALAR_TYPE *currents_y[N_x*(N_y-1)];

static void random_init() {
    for (int i = 0; i < N_x; ++i) {
      for (int j = 0; j < N_y; ++j) {
        // generate random number in [-pi,pi)
        phases[i + N_y * j] = 2*PI * ((float) rand() / (float) RAND_MAX) - PI;
      }
    }
}

static void run_step(SCALAR_TYPE temp) {
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      SCALAR_TYPE div_I_super = 0;
      SCALAR_TYPE phase = phases[i + N_y * j];
      if (i < N_x-1) {
        div_I_super += SINE(phases[i+1 + N_y*j] - phase - frustration * j);
      }
      if (i > 0) {
        div_I_super += SINE(phases[i-1 + N_y*j] - phase + frustration *j);
      }
      if (j < N_y - 1) {
        div_I_super += SINE(phases[i + N_y*(j+1)] - phase);
      }
      if (j > 0) {
        div_I_super += SINE(phases[i + N_y*(j-1)] - phase);
      }
      phases[i + N_y *j] += EPS * div_I_super;
      if (temp > 0) {
        // float local_temp = temp *(1 + temp_scale_factor * (abs(i - N_x/2) + abs(j - N_y/2)) / ((float) N_x/2+ (float) N_y/2));
        //printf("i=%d, j=%d, local_temp = %g\n",i,j, local_temp);
        phases[i + N_y * j] += 2*temp*((float) rand() / (float) RAND_MAX) - temp;
      }
    }
  }
}

static SCALAR_TYPE free_energy() {
  SCALAR_TYPE f = 0;
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      SCALAR_TYPE phase = phases[i + N_y * j];
      if (i < N_x-1) {
        f += 1-COS(phases[i+1 + N_y*j] - phase - frustration * j);
      }
      if (j < N_y - 1) {
        f += 1-COS(phases[i + N_y*(j+1)] - phase);
      }
    }
  }
  return f;
}

static void save_to_file(char *basename) {
  char output_filename[300];


  //
  // save phases to basename.dat
  //
  
  sprintf(output_filename, "%s.%s", basename, ".dat");
  FILE *file = fopen(output_filename, "w");
  fprintf(file, "#\t\tcolumn(i)\t\trow(j)\t\t\t\tphase\n");
  for (int j = 0; j < N_y; ++j) {
    for (int i = 0; i < N_x; ++i) {
      fprintf(file, "%d\t\t%d\t\t%.5g\n",i,j,phases[i+N_y*j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  
  //
  // save I_x to basename_Ix.dat
  //

  sprintf(output_filename, "%s_Ix.dat", basename);
  file = fopen(output_filename, "w");
  fprintf(file, "#\t\tcolumn(i)\t\trow(j)\t\t\tI_x\n");
  for (int j = 0; j < N_y; ++j) {
    for (int i = 0; i < N_x-1; ++i) {
      SCALAR_TYPE I_x =  SINE(phases[i+1 + N_y*j] - phases[i+N_y*j] - frustration * j);
      fprintf(file, "%d\t\t%d\t\t%.5g\n",i,j,I_x);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  
  
  //
  // save I_y to basename_Ix.dat
  //
  sprintf(output_filename, "%s_Iy.dat", basename);
  file = fopen(output_filename, "w");
  fprintf(file, "#\t\tcolumn(i)\t\trow(j)\t\t\tI_y\n");
  for (int j = 0; j < N_y-1; ++j) {
    for (int i = 0; i < N_x; ++i) {
      SCALAR_TYPE I_y =  SINE(phases[i + N_y*(j+1)] - phases[i + N_y*j]);
      fprintf(file, "%d\t\t%d\t\t%.5g\n",i,j,I_y);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  
}

int
main (int argc, char **argv)
{
  // parse command line options
  int num_steps = 10;
  float T_start = 1;

  int c;
  while ((c = getopt(argc, argv, "f:n:t:")) != -1)
    switch (c) {
    case 'f':
      frustration = 2 * PI * strtof(optarg, NULL);
      break;
    case 'n':
      num_steps = strtol(optarg, NULL, 10);
      break;
    case 't':
      T_start = strtof(optarg, NULL);
      break;
    case '?':
      if (optopt == 'f' || optopt == 'n' || optopt == 't')
        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf (stderr,
                 "Unknown option character `\\x%x'.\n",
                 optopt);
      return 1;
    default:
      abort ();
    }
    
  random_init();
  printf("Running ground state annealer\n");
  printf("N_x = %d, N_y = %d\n", N_x, N_y);
  for (int i = 0; i < num_steps; ++i) {
    float temp = T_start * ( 1 - (float) i / num_steps);
    if (i % 100 == 0)
      printf("i = %10d, temp = %10g, f = %.13g\n", i, temp, free_energy());
    run_step(temp);
  }  
  save_to_file("output");
  return 0;
}
