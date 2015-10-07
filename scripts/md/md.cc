#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<iostream>

//parameters

#define L 30 //lattice length
#define dx 10 //lattice spacing
#define Lx (L/dx) //renormalised length
#define N (Lx*Lx*Lx) //total number of particles
#define NSTEP 100000 //total number of time steps

using namespace std;

//constants

double const vmax = 2.7; //T_LJ = 2.4
double const dt = 0.01;
double const rcut = 3.0;

//variables

double x[N];
double y[N];
double z[N];
double vx[N];
double vy[N];
double vz[N];
double fx[N];
double fy[N];
double fz[N];
double fx0[N];
double fy0[N];
double fz0[N];
double xi[N];
double yi[N];
double zi[N];
double T; //reduced temperature
double D; //diffusion constant

//functions

void init();
void move();
int overlap(int i);
double get_distance2(int i, int j);
void force();

FILE *file;

double get_distance2(int i, int j) //squared distance
{
 double xmin, ymin, zmin, diff, r2;
      diff = x[i] - x[j];
      xmin = (diff - L*round(diff/L));
      diff = y[i] - y[j];
      ymin = (diff - L*round(diff/L));
      diff = z[i] - z[j];
      zmin = (diff - L*round(diff/L));

 r2 = xmin*xmin + ymin*ymin + zmin*zmin;

 return r2;
}

void init() //initialisation
{

 /*Position initialisation*/
 
 int l;
 
 for(int i = 0; i < Lx; i++){
    for(int j = 0; j < Lx; j++){
      for(int k = 0; k < Lx; k++){
	l = (i*Lx*Lx + j*Lx + k);
	x[l] = i*dx;
	y[l] = j*dx;
	z[l] = k*dx;
      }
    }  
  }dran

 /*Velocity initialisation*/
 
 double vxc, vyc, vzc = 0.;

 for(int i = 0; i < N; i++){
    vx[i] = vmax*(2.*drand48() - 1.);
    vy[i] = vmax*(2.*drand48() - 1.);
    vz[i] = vmax*(2.*drand48() - 1.);
    
    vxc += vx[i];
    vyc += vy[i];
    vzc += vz[i];
  }

 vxc /= N;
 vyc /= N;
 vzc /= N;

 for(int i = 0; i < N; i++){
    vx[i] -= vxc;
    vy[i] -= vyc;
    vz[i] -= vzc;
  }
  printf("%f \t %f \t %f \n", vxc/vmax, vyc/vmax, vzc/vmax);
  
  for(int i = 0; i < N; i++){ //initial positions
    xi[i] = x[i];
    yi[i] = y[i];
    zi[i] = z[i];
   //printf("%d \t %f \t %f \t %f \n", i, x[i], y[i], z[i]);
  }
}

int overlap(int i) //overlap function
{
  double r2;
  
  for(int j = 0; j < N; j++){
    r2 = get_distance2(i,j);
  
    if(i != j && r2 == 0.){
      return 1; //overlap
    }
  }
  
  return 0;
}

void force() //total force
{ 
  double r2, r6, r12, xmin, ymin, zmin, diff, f;
  
  for(int i = 0; i < N; i++){
    fx[i] = 0;
    fy[i] = 0;
    fz[i] = 0;
  }
  
  for(int i = 0; i < N; i++){
    for(int j = i+1; j < N; j++){
      diff = x[i] - x[j];
      xmin = (diff - L*round(diff/L));
      diff = y[i] - y[j];
      ymin = (diff - L*round(diff/L));
      diff = z[i] - z[j];
      zmin = (diff - L*round(diff/L));
      
      r2 = xmin*xmin + ymin*ymin + zmin*zmin;
      
      if(r2 < rcut*rcut){
      	r6 = r2*r2*r2;
      	r12 = r6*r6;
      	f = (48./(r12*r2) - 24./(r6*r2)); //LJ-potential
      	fx[i] += f*xmin;
      	fy[i] += f*ymin;
      	fz[i] += f*zmin;
      	fx[j] -= f*xmin;
      	fy[j] -= f*ymin;
      	fz[j] -= f*zmin;
      }
    }
  }
}

void move() //one time step move dt
{
  //double fx1, fy1, fz1;
  
  force();
  
  for(int i = 0; i < N; i++){
    x[i] += dt*vx[i] + 0.5*dt*dt*fx[i]; //updated positions
    y[i] += dt*vy[i] + 0.5*dt*dt*fy[i];
    z[i] += dt*vz[i] + 0.5*dt*dt*fz[i];
    
    fx0[i] = fx[i];
    fy0[i] = fy[i];
    fz0[i] = fz[i];
  }
  
  force();
  
  for(int i = 0; i < N; i++){
    vx[i] += 0.5*dt*(fx0[i] + fx[i]); //updated velocities
    vy[i] += 0.5*dt*(fy0[i] + fy[i]);
    vz[i] += 0.5*dt*(fz0[i] + fz[i]);
  }
}

int main(int argc, char *argv[])
{
  
  double r2, v2, r20, rms;

  int seed = time(0);

  srand48(seed);

  init();
 
  file = fopen("temperature_50.dat","w"); 
 
  for(int k = 0; k < NSTEP; k++){
     move();
    if(k % 500 == 0){
      v2 = 0.;
      for(int i = 0; i < N; i++){
        v2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      }
      
      T = 1./3.*v2/N;
     
      fprintf(file, "%f \t %f \n", k*dt, T);
    }
  }

  fclose(file);

  r2, r20 = 0.;

  for(int i = 0; i < N; i++){
    r2 += x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
    
    r20 += xi[i]*xi[i] + yi[i]*yi[i] + zi[i]*zi[i];
  }

  rms = (r2 - r20)/N;

  D = rms/(6.*NSTEP*dt);

  printf("Temperature: T = %f \t r2 = %f \t rms = %f \t Diffusion: D =  %f \n", T, r20/N, sqrt(rms), D);

  return 0;

}


