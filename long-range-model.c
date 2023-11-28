#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/stat.h> // for mkdir
#include <time.h>

static int N_x = 100;
static int N_y = 100;
#define SINE sinf
#define COS cosf
#define PI 3.141592

static const float EPS  = 0.4;
static float frustration = 0;
float xchi = 1.0;


static float *phases;
static float phi0_y = 0;



static uint32_t random_state = 1;

// ^
// |
//
// j (row)
//    i (column)   - >


// static fields are all initialized to zero
// static float currents_x[(N_x-1)*N_y];
// static float *currents_y[N_x*(N_y-1)];

static inline uint32_t lcg_parkmiller() {
  // faster than the glibc rand() function
  random_state = (uint64_t)random_state * 48271 % 0x7fffffff;
  return random_state;
}


static inline float random_phase(float temp) {
  // return random float in interval [-temp, temp]
  return temp * (2.0f * (float) lcg_parkmiller() / 0x7fffffff - 1);
}

static void random_init() {
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      // generate random number in [-pi,pi)
      phases[i + N_y * j] = 2*PI * random_phase(0.5); //((float) rand() / (float) RAND_MAX) - PI;
    }
  }
}

static void normalize_phases() {
  // normalize all phases to range [-pi,pi)
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      phases[i + N_y * j] = remainderf(phases[i + N_y * j], 2*PI);
    }
  }
}

static inline void run_step(float temp) {
  for (int i1 = 0; i1 < N_x; ++i1) {
    for (int j1 = 0; j1 < N_y; ++j1) {
      float phase1 = phases[i1 + N_x * j1];
      float div_I_super = 0;
      for (int i2 = 0; i2 < N_x; ++i2) {
	for (int j2 = 0; j2 < N_y; ++j2) {
	  float phase2 = phases[i2 + N_x * j2];
	  float r_squared = (i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2);
	  div_I_super += SINE(phase2 - phase1) * expf(-r_squared / xchi*xchi);
	}
      }
      phases[i1 + N_x *j1] += EPS * div_I_super;
      phases[i1 + N_x * j1] += random_phase(temp);
	
    }
  }
}



static inline float free_energy() {
  float f = 0;
  for (int i1 = 0; i1 < N_x; ++i1) {
    for (int j1 = 0; j1 < N_y; ++j1) {
      float phase1 = phases[i1 + N_x * j1];
      for (int i2 = 0; i2 < N_x; ++i2) {
	for (int j2 = 0; j2 < N_y; ++j2) {
	  float phase2 = phases[i2 + N_x * j2];
	  float r_squared = (i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2);
	  f += (1-COS(phase2 - phase1)) * expf(-r_squared / xchi*xchi);
	}
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
  
  sprintf(output_filename, "%s.dat", basename);
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
      float I_x =  SINE(phases[i+1 + N_y*j] - phases[i+N_y*j] - frustration * j);
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
      float I_y =  SINE(phases[i + N_y*(j+1)] - phases[i + N_y*j]);
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
  int num_steps = 10000;
  float T_start = 1;
  int c;
  while ((c = getopt(argc, argv, "n:t:N:")) != -1)
    switch (c) {
      /* case 'f': */
      /*   frustration = 2 * PI * strtof(optarg, NULL); */
      /*   break; */
    case 'n':
      num_steps = strtol(optarg, NULL, 10);
      break;
    case 't':
      T_start = strtof(optarg, NULL);
      break;
    case 'N':
      N_x = strtol(optarg, NULL, 10);
      N_y = N_x;
      break;
      /* case 'p': */
      /*   phi0_y = PI * strtof(optarg, NULL); */
      /*   break; */
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
  char output_dir[200];
  char output_file_voltage[300];
  FILE *file_IV;
  /* create output directory */
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  snprintf(output_dir,
	   sizeof(output_dir),
	   "%d-%02d-%02d_%02d-%02d-%02d_frustration_sweep_Nx=%d_Ny=%d_annealing-steps=%d",
	   tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec,
	   N_x,
	   N_y,
	   num_steps);
  /* mkdir(output_dir, 0777); */
  /* printf("output dir: %s\n", output_dir); */
  /* snprintf(output_file_voltage, sizeof(output_file_voltage), "%s/IV.dat", output_dir); */
  char output_filename[100];
  snprintf(output_filename, sizeof(output_filename), "output_Nx=%d_Ny=%d_nsteps=%d_Tstart=%g", N_x, N_y, num_steps, T_start);
  phases = (float *) calloc(N_x * N_y, sizeof(float));
  random_init();

  for (int i = 0; i < num_steps; ++i) {
    float temp = T_start * ( 1 - (float) i / num_steps);
    //if (i % 100 == 0)
    printf("i = %10d, temp = %10g, f = %.13g\n", i, temp, free_energy());
    run_step(temp);
  }
  save_to_file(output_filename);
}
