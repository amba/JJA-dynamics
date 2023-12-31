#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/stat.h> // for mkdir
#include <time.h>

static int N_x = 200;
static int N_y = 200;
#define SINE sinf
#define COS cosf
#define PI 3.141592

static const float EPS  = 0.4;
static float frustration = 0;


static float *phases;
static float phi0_y = 0.2*PI;
// [N_x*N_y + 2]; // phi(i,j) = phases(i + N_y * j)
// φ_L = phases[Nx*Ny] (left electrode)
// φ_R = phases[Nx*Ny + 1] (right electrode)



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
    phases[N_x*N_y] = 2*PI*random_phase(0.5);
    phases[N_x*N_y+1] = 2*PI*random_phase(0.5);
}

static void normalize_phases() {
  // normalize all phases to range [-pi,pi)
    for (int i = 0; i < N_x; ++i) {
      for (int j = 0; j < N_y; ++j) {
        phases[i + N_y * j] = remainderf(phases[i + N_y * j], 2*PI);
      }
    }
    phases[N_x*N_y] = remainderf(phases[N_x*N_y], 2*PI);
    phases[N_x*N_y+1] = remainderf(phases[N_x*N_y+1], 2*PI);
}

static inline void run_step(float temp, float I_bias) {
  float phi_L = phases[N_x*N_y];
  float phi_R = phases[N_x*N_y+1];
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      float div_I_super = 0;
      float phase = phases[i + N_y * j];

      // x-direction
      if (i == 0)  {
        div_I_super += SINE(phases[i+1 + N_y*j] - phase - frustration * j);
        div_I_super += SINE(phi_L - phase + frustration *j);
      }
      else if (i == N_x-1) {
        div_I_super += SINE(phases[i-1 + N_y*j] - phase + frustration *j);
        div_I_super += SINE(phi_R - phase - frustration * j);
      }
      else {
        div_I_super += SINE(phases[i+1 + N_y*j] - phase - frustration * j);
        div_I_super += SINE(phases[i-1 + N_y*j] - phase + frustration *j);
      }

      // y-direction
      if (j < N_y - 1) {
        div_I_super += SINE(phases[i + N_y*(j+1)] - phase - phi0_y);
      }
      if (j > 0) {
        div_I_super -= SINE(phase - phases[i + N_y*(j-1)] - phi0_y);
      }
      phases[i + N_y *j] += EPS * div_I_super;
      phases[i + N_y * j] += random_phase(temp);
    }
  }
  // update electrodes
  float IL = 0;
  float IR = 0;
  float div_IL = 0;
  float div_IR = 0;
  // use Newton method to find the solutions of I_L = I_bias and I_R = I_bias
  for (int j = 0; j < N_y; ++j) {
    IL += SINE(phases[N_y*j] - phi_L - frustration * j);
    IR += SINE(phi_R - phases[N_y*j + N_x-1] - frustration * j);
    div_IL -= COS(phases[N_y*j] - phi_L -frustration*j);
    div_IR += COS(phi_R - phases[N_y*j+N_x  -1] - frustration * j);
  }
  /* printf("IL = %.3g, IR = %.3g\n", IL, IR); */
  /* printf("φL = %.5g, φR = %.5g\n", phi_L / (2*PI), phi_R / (2*PI)); */
  phases[N_x*N_y] += (I_bias - IL) / div_IL;
  phases[N_x*N_y+1] += (I_bias - IR) / div_IR;
}

static inline float free_energy() {
  float f = 0;
  float phi_L = phases[N_x*N_y];
  float phi_R = phases[N_x*N_y+1];
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      float phase = phases[i + N_y * j];

      if (i == 0) 
        f += 1 -COS(phase - phi_L - frustration * j);

      if (i < N_x-1) 
        f += 1-COS(phases[i+1 + N_y*j] - phase - frustration * j);
      else
        f += 1 - COS(phi_R - phase - frustration * j);
      
      if (j < N_y - 1) 
        f += 1-COS(phases[i + N_y*(j+1)] - phase - phi0_y);
      
    }
  }
  return f;
}

static void save_to_file(char *basename) {
  char output_filename[300];


  //
  // save phases to basename.dat
  //
  printf("saving data to %s\n", basename);
  sprintf(output_filename, "%s.%s", basename, "dat");
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
  float phi_L = phases[N_x*N_y];
  float phi_R = phases[N_x*N_y+1];

  
  sprintf(output_filename, "%s_Ix.dat", basename);
  file = fopen(output_filename, "w");
  fprintf(file, "#\t\tcolumn(i)\t\trow(j)\t\t\tI_x\n");
  for (int j = 0; j < N_y; ++j) {
    float I_x = SINE(phases[N_y*j] - phi_L - frustration * j);
    fprintf(file, "%d\t\t%d\t\t%.5g\n",0,j,I_x);
    for (int i = 0; i < N_x-1; ++i) {
      I_x =  SINE(phases[i+1 + N_y*j] - phases[i+N_y*j] - frustration * j);
      fprintf(file, "%d\t\t%d\t\t%.5g\n",i+1,j,I_x);
    }
    I_x = SINE(phi_R - phases[N_y*j  +N_x-1] - frustration * j);
    fprintf(file, "%d\t\t%d\t\t%.5g\n",N_x+1,j,I_x);
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
      float I_y =  SINE(phases[i + N_y*(j+1)] - phases[i + N_y*j] - phi0_y);
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
  int num_steps = 20000;
  //  float T_start = 1;
  int c;
  while ((c = getopt(argc, argv, "f:n:t:N:")) != -1)
    switch (c) {
    case 'f':
      frustration = 2 * PI * strtof(optarg, NULL);
      break;
    case 'n':
      num_steps = strtol(optarg, NULL, 10);
      break;
    /* case 't': */
    /*   T_start = strtof(optarg, NULL); */
    /*   break; */
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
  char output_dir[100];
  char output_file[200];
  //  FILE *file_IV;
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
  mkdir(output_dir, 0777);
  printf("output dir: %s\n", output_dir);
  snprintf(output_file, sizeof(output_file), "%s/phases.dat", output_dir); 
  /* file_IV = fopen(output_file_voltage, "w"); */
  /* fprintf(file_IV, "#\t\tf\t\tj_bias\t\tV\t\tΔφ_end\n"); */
  
  phases = (float *) calloc(N_x * N_y + 2, sizeof(float));
  float temp = 0;
  float I_bias = 0.1 * N_y;
  printf("N_x = %d, N_y = %d\n", N_x, N_y);
  random_init();
  for (int i = 0; i < num_steps ; ++i) {
    run_step(temp, I_bias);
    if (i % 100 == 0) {
      printf("i = %d, F = %g\n", i, free_energy());
      printf("phi_L = %g\n", phases[N_x*N_y]);
      printf("phi_R = %g\n", phases[N_x*N_y+1]);
    }
    
      
  }
  save_to_file(output_file);
}
