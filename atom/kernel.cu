#include "kernel.h"
#include <stdio.h>
#define TX 32
#define TY 32
#define LEN 5.f
#define TIME_STEP 0.005f
#define FINAL_TIME 10.f
#define PI           3.14159265358979323846  /* pi */
#define SCALE_RGB (255.f / LEN) * 30 

// scale coordinates onto [-LEN, LEN]
__device__
float scale(int i, int w) { return 2 * LEN*(((1.f*i)/w) - 0.5f); }

// function for right-hand side of y-equation
__device__
float f(float x, float y, float param, float sys) {
  if (sys == 1) return x - 2 * param*y; // negative stiffness
  if (sys == 2) return -x + param*(1 - x*x)*y; //van der Pol
  else return -x - 2 * param*y;
}

// explicit Euler solver
__device__
float2 euler(float x, float y, float dt, float tFinal,
             float param, float sys) {
  float dx = 0.f, dy = 0.f;
  for (float t = 0; t < tFinal; t += dt) {
    dx = dt*y;
    dy = dt*f(x, y, param, sys);
    x += dx;
    y += dy;
  }
  return make_float2(x, y);
}

__device__
unsigned char clip(float x){ return x > 255 ? 255 : (x < 0 ? 0 : x); }


// Source of equations: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Map%3A_Physical_Chemistry_for_the_Biosciences_(Chang)/11%3A_Quantum_Mechanics_and_Atomic_Structure/11.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom

__device__
float density(float ra, float th, float ph, int Z, int nlm) {

//  const float a0 = 5.29177210903E-11;
  const float a0 = 1.5;
  const float p = (Z * ra) / a0;

  float result = 0.0;

  if (nlm == 100) {
    result = 1 / sqrt(PI);
    result *= pow (Z / a0, 3.0/2.0) * exp(-1 * p);
  } else if (nlm == 210) { 
    result = 1 / sqrt(PI * 32);
    result *= pow (Z / a0, 3.0/2.0);
    result *= p * exp(-1 * p / 2);
    result *= cos(th);
  } else if (nlm == 211) { 
    result = 1 / sqrt(PI * 64);
    result *= pow (Z / a0, 3.0/2.0);
    result *= p * exp(-1 * p / 2);
    result *= sin(th) * exp(ph);
  } else if (nlm == 322) { 
    result = 1 / sqrt(PI * 162);
    result *= pow (Z / a0, 3.0/2.0);
    result *= pow(p,2) * exp(-1 * p / 3);
    result *= pow(sin(th),2) * exp(ph * 2);
  } else if (nlm == 321) { 
    result = 1 / sqrt(PI * 81);
    result *= pow (Z / a0, 3.0/2.0);
    result *= pow(p,2) * exp(-1 * p / 3.0);
    result *= sin(th) * cos(th) * exp(ph);
  } else if (nlm == 320) { 
    result = 1 / sqrt(PI * 81);
    result *= pow (Z / a0, 3.0/2.0);
    result *= pow(p,2) * exp(-1.0 * p / 3.0);
    result *= (3 * pow(cos(th),2)) - 1;
  } else if (nlm == 310) { 
    result = (1.0 / 81.0) * sqrt(2.0 / PI);
    result *= pow (Z / a0, 3.0/2.0);
    result *= (6 * ra - pow(p,2.0)) * exp(-1 * p / 3.0);
    result *= cos(th);
  }
  return result; 
}

// kernel function to compute decay and shading
__global__
void stabImageKernel(uchar4 *d_out, int w, int h, float p, int s, float z0, int nlm) {
  const int c = blockIdx.x*blockDim.x + threadIdx.x;
  const int r = blockIdx.y*blockDim.y + threadIdx.y;
  if ((c >= w) || (r >= h)) return; // Check if within image bounds
  const int i = c + r*w; // 1D indexing
  
  const float x0 = scale(c, w);
  const float y0 = scale(r, h);
  const float dist_0 = 1 / sqrt(x0*x0 + y0*y0);
  const float dist_1 = sin(x0*x0 + y0*y0);
  const float dist_2 = cos(x0*x0);
  const float2 pos = euler(x0, y0, TIME_STEP, FINAL_TIME, p, s);
  const float dist_f = sqrt(pos.x*pos.x + pos.y*pos.y);

  const float ra = sqrt(x0*x0 + y0*y0 + z0*z0);
  const float th = atan2(y0 , x0);
  const float ph = atan2(sqrt(x0*x0 + y0*y0) , z0);
  
  float dens = 0;

  if (nlm >= 300) 
    dens = density(ra, th, ph, 9, nlm) + density(ra, th, ph, 4, 210);
  else 
    dens = density(ra, th, ph, 5, nlm);

  // printf("density: %.2f", dens);

  // assign colors based on distance from origin
  const float dist_r = dist_f / dist_0;
  // d_out[i].x = clip(abs(ph) * 255); // red ~ growth
  d_out[i].x = (dens >= 0) ? clip(dens * SCALE_RGB) : 0; // red ~ growth
  d_out[i].y = ((c == w / 2) || (r == h / 2)) ? 255 : 0; // axes
  d_out[i].z = (dens < 0) ? clip(dens * -1 * SCALE_RGB) : 0; // blue ~ 1/growth
  d_out[i].w = 255;
}

void kernelLauncher(uchar4 *d_out, int w, int h, float p, int s, float z0, int nlm) {
  const dim3 blockSize(TX, TY);
  const dim3 gridSize = dim3((w + TX - 1)/TX, (h + TY - 1)/TY);
  stabImageKernel<<<gridSize, blockSize >>>(d_out, w, h, p, s, z0, nlm);
}
