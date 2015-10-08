#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<iostream>

//parameters

// #define L 50 //lattice length
// #define dx 10 //lattice spacing
// #define Lx (L/dx) //renormalised length
// #define N (Lx*Lx*Lx) //total number of particles
// #define NSTEP 500000 //total number of time steps

#define L 50 //lattice length
#define dx 10 //lattice spacing
#define Lx (L/dx) //renormalised length
#define N (Lx*Lx*Lx) //total number of particles
#define NSTEP 1000000 //total number of time steps

using namespace std;

//constants

double const vmax     = 2.7; //T_LJ = 2.4
double const dt       = 0.001;
double const rcut     = 3.0;
double const rcut2    = rcut*rcut;
double const halfdtdt = 0.5*dt*dt;
double const halfdt   = 0.5*dt;


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
double fx1[N];
double fy1[N];
double fz1[N];
double xi[N];
double yi[N];
double zi[N];
double T; //reduced temperature
double D; //diffusion constant

// file objects

FILE *file_x, *file_y, *file_z, *file_vx, *file_vy, *file_vz, *file_fx, *file_fy, *file_fz;


//functions

void init();
void move();
void force();

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
  }

 /*Velocity initialisation*/
 
 double vxcm, vycm, vzcm = 0.;

 for(int i = 0; i < N; i++){

    vx[i] = vmax*(2.*drand48()-1.);
    vy[i] = vmax*(2.*drand48()-1.);
    vz[i] = vmax*(2.*drand48()-1.);

    vxcm += vx[i]; // summing velocities (for subtracting center of mass V)
    vycm += vy[i];
    vzcm += vz[i];

  }

 /*Subtract center of mass velocities*/

 vxcm /= N; // averaging velocities (for subtracting center of mass V)
 vycm /= N;
 vzcm /= N;

 for(int i = 0; i < N; i++){
    vx[i] -= vxcm;
    vy[i] -= vycm;
    vz[i] -= vzcm;
  }

  printf("\t %f \t %f \t %f \n", vxcm/vmax, vycm/vmax, vzcm/vmax);
 
/* */ 
  for(int i = 0; i < N; i++){ //initial positions
    xi[i] = x[i];
    yi[i] = y[i];
    zi[i] = z[i];
   //printf("%d \t %f \t %f \t %f \n", i, x[i], y[i], z[i]);
  }
}


void force() //total force
{ 
  double r2, xmin, ymin, zmin, diff, f;
  
  for(int i = 0; i<N; i++){
    fx[i] = 0;
    fy[i] = 0;
    fz[i] = 0;
  }
  
  for(int i = 0; i<N; i++){
    for(int j = i+1; j<N; j++){
      // @test:0034:replacing "diff" with the raw x[i]-x[i] calcs
      // diff = x[i] - x[j];
      // xmin = (diff - L*round(diff/L));
      // diff = y[i] - y[j];
      // ymin = (diff - L*round(diff/L));
      // diff = z[i] - z[j];
      // zmin = (diff - L*round(diff/L));
      // @test:0034: ^
      xmin = ((x[i]-x[j]) - L*round((x[i]-x[j])/L));
      ymin = ((y[i]-y[j]) - L*round((y[i]-y[j])/L));
      zmin = ((z[i]-z[j]) - L*round((z[i]-z[j])/L));

      r2 = xmin*xmin + ymin*ymin + zmin*zmin;
      
      if(r2 < rcut2){
        // @test:@0005: trying to use pow() instead of the below three lines of code (commented)
        // r6 = r2*r2*r2;
        // r12 = r6*r6;
        // f = (48./(r12*r2) - 24./(r6*r2)); //LJ-potential
        // @test:@0005: testing below 3 lines instead of above  
        f = (48./(pow(r2,7)) - 24./(pow(r2,4))); //LJ-potential

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
  //force();        // fx, fy and fz stores forces for f(t)  // @test:0021:moving force out of the loop

  // MOVE_R 
  for(int i = 0; i < N; i++){
    x[i] += dt*vx[i] + halfdtdt*fx[i]; //updated positions
    y[i] += dt*vy[i] + halfdtdt*fy[i];
    z[i] += dt*vz[i] + halfdtdt*fz[i];
    
    fx1[i] = fx[i]; // fx1 stores forces at f(t)
    fy1[i] = fy[i]; // fy1 stores forces at f(t)
    fz1[i] = fz[i]; // fz1 stores forces at f(t)
  }

  //fprintf(file, "%f \t %f \n", k*dt, T);  
  
  force();          // fx, fy and fz stores forces for f(t+1)
  
  // MOVE_V
  for(int i = 0; i < N; i++){
    vx[i] += halfdt*(fx1[i] + fx[i]); //updated velocities
    vy[i] += halfdt*(fy1[i] + fy[i]);
    vz[i] += halfdt*(fz1[i] + fz[i]);
  }
}

int main(int argc, char *argv[])
{
  
  double r2, v2, r20, rms, time_spent;


  // Code timing:start
  clock_t begin, end;
  begin = clock();

  // @test:@0011:using 1 as seed for testing for runtime    @todo: make the srand48 truely random later by using time(0) as the seed
  // srand48(time(0));  // seed for random initial particle velocities

  // @test:@0011: ^
  srand48(1);  // seed for random initial particle velocities


  printf("Initiating particle's coordinates and velocities...\n");
  
  init();            // initialize 

  // @todo: uncomment to save data to file... 
  // file = fopen("temperature_50.dat","w"); 

  printf("Simulating particle's dynamics...\n");


  // void init_file(char title){

  //   // e.g. title: Coordinates X

  //   file_x = fopen("x.dat","w"); 
  //   fprintf(file_x, "%s\n",title);
  //   for (int i = 0; i < N; i++){
  //     fprintf(file_x, "Molecule %i\t", i); 
  //   }
  //   fprintf(file_x, "\n");

  // }

  file_x  = fopen("../../data/x.dat","w"); 
  file_y  = fopen("../../data/y.dat","w"); 
  file_z  = fopen("../../data/z.dat","w"); 
  file_vx = fopen("../../data/vx.dat","w"); 
  file_vy = fopen("../../data/vy.dat","w"); 
  file_vz = fopen("../../data/vz.dat","w"); 
  file_fx = fopen("../../data/fx.dat","w"); 
  file_fy = fopen("../../data/fy.dat","w"); 
  file_fz = fopen("../../data/fz.dat","w"); 

  fprintf(file_x,  "X(t) Coordinates\n");
  fprintf(file_y,  "Y(t) Coordinates\n");
  fprintf(file_z,  "Z(t) Coordinates\n");
  fprintf(file_vx, "X(t) Velocities\n");
  fprintf(file_vy, "Y(t) Velocities\n");
  fprintf(file_vz, "Z(t) Velocities\n");
  fprintf(file_fx, "X(t) Forces\n");
  fprintf(file_fy, "Y(t) Forces\n");
  fprintf(file_fz, "Z(t) Forces\n");

  for (int i = 0; i < N; i++){
    fprintf(file_x, "Molecule %i\t", i); 
    fprintf(file_y, "Molecule %i\t", i);
    fprintf(file_z, "Molecule %i\t", i);
    fprintf(file_vx, "Molecule %i\t", i);
    fprintf(file_vy, "Molecule %i\t", i);
    fprintf(file_vz, "Molecule %i\t", i);
    fprintf(file_fx, "Molecule %i\t", i);
    fprintf(file_fy, "Molecule %i\t", i);
    fprintf(file_fz, "Molecule %i\t", i);
  }

  fprintf(file_x, "\n");
  fprintf(file_y, "\n");
  fprintf(file_z, "\n");
  fprintf(file_vx, "\n");
  fprintf(file_vy, "\n");
  fprintf(file_vz, "\n");
  fprintf(file_fx, "\n");
  fprintf(file_fy, "\n");
  fprintf(file_fz, "\n");

  //file_T = fopen("temperature_50.dat","w"); 

  force();

  for(int k = 0; k < NSTEP; k++){

    if(k % 500 == 0){
      printf("\t %f percent\n",((float)k/NSTEP)*100.0);  
    }
    
    move();

    for (int i = 0; i < N; i++){
      fprintf(file_x, "%f\t",  x[i]); 
      fprintf(file_y, "%f\t",  y[i]);
      fprintf(file_z, "%f\t",  z[i]);
      fprintf(file_vx, "%f\t", vx[i]);
      fprintf(file_vy, "%f\t", vy[i]);
      fprintf(file_vz, "%f\t", vz[i]);
      fprintf(file_fx, "%f\t", fx[i]);
      fprintf(file_fy, "%f\t", fy[i]);
      fprintf(file_fz, "%f\t", fz[i]);
    }
    fprintf(file_x, "\n");
    fprintf(file_y, "\n");
    fprintf(file_z, "\n");
    fprintf(file_vx, "\n");
    fprintf(file_vy, "\n");
    fprintf(file_vz, "\n");
    fprintf(file_fx, "\n");
    fprintf(file_fy, "\n");
    fprintf(file_fz, "\n");

    // @todo: uncomment up to "fprintf" to get temperature
    // if(k % 500 == 0){
    //   v2 = 0.;
    //   for(int i = 0; i < N; i++){
    //     v2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    //   }
      
    //   T = 1./3.*v2/N;
     
    //   fprintf(file, "%f \t %f \n", k*dt, T);
    // }
  }

  // @todo: uncomment to save data to file...
  // fclose(file);
  fclose(file_x);
  fclose(file_y);
  fclose(file_z);
  fclose(file_vx);
  fclose(file_vy);
  fclose(file_vz);
  fclose(file_fx);
  fclose(file_fy);
  fclose(file_fz);

  /* Calculate diffusion coefficient   @todo: make this it's own function */

  r2, r20 = 0.;
  for(int i = 0; i < N; i++){
    r2  += x[i]*x[i]   + y[i]*y[i]   + z[i]*z[i];  
    r20 += xi[i]*xi[i] + yi[i]*yi[i] + zi[i]*zi[i];
  }
  rms = (r2 - r20)/N;
  D = rms/(6.*NSTEP*dt);

  printf("Temperature: T = %f \t r2 = %f \t rms = %f \t Diffusion: D =  %f \n", T, r20/N, sqrt(rms), D);

  // Code timing:end
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("Runtime: %f seconds \n",time_spent);

  return 0;

}