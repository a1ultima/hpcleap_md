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
#define NSTEP 1000 //total number of time steps

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
double fx0[N];
double fy0[N];
double fz0[N];
double xi[N];
double yi[N];
double zi[N];
double T; //reduced temperature
double D; //diffusion constant

// file objects

//FILE *file_x, *file_y, *file_z, *file_vx, *file_vy, *file_vz, *file_fx, *file_fy, *file_fz, *file_T;

FILE *file[9];
const char *dirs[]    ={"../../data/x.dat" , "../../data/y.dat","../../data/z.dat" ,"../../data/vx.dat","../../data/vy.dat","../../data/vz.dat","../../data/fx.dat","../../data/fy.dat","../../data/fz.dat" };
const char *headers[] ={"#X(t) Coordinates\n#Time (t)\t","#Y(t) Coordinates\n#Time (t)\t","#Z(t) Coordinates\n#Time (t)\t","#X(t) Velocities\n#Time (t)\t","#Y(t) Velocities\n#Time (t)\t","#Z(t) Velocities\n#Time (t)\t","#X(t) Forces\n#Time (t)\t","#Y(t) Forces\n#Time (t)\t","#Z(t) Forces\n#Time (t)\t"};

//functions

void init();
void move();
void force();
void temperature();
void diffusion_coefficient(double *, double *, double *);


// main program

int main(int argc, char *argv[])
{
  
  double v2, rms, D, time_spent, r20;  // double r2, v2, r20, rms, time_spent;

  // Code timing:start
  clock_t begin, end;
  begin = clock();

  // @test:@0011:using 1 as seed for testing for runtime    @todo: make the srand48 truely random later by using time(0) as the seed
  // srand48(time(0));  // seed for random initial particle velocities

  // @test:@0011: ^
  srand48(1);  // seed for random initial particle velocities


  printf("Initiating particle's coordinates and velocities...\n");
  
  init();            // initialize 

  printf("Simulating particle's dynamics...\n");

  // open files and initiate headers  @file
  for(int fi=0; fi<9; fi++ )
  {
    file[fi] = fopen(dirs[fi],"w"); 
    fprintf(file[fi],"%s",headers[fi]);

    for (int i = 0; i < N; i++){
      fprintf(file[fi],  "Molecule %i\t", i); 
    }
    fprintf(file[fi],  "\n");
  }
    

  // Initial force calculation
  force();

  // Increment time (t + dt)
  for(int k = 0; k < NSTEP; k++){

    if(k % 500 == 0){
      printf("\t %f percent\n",((float)k/NSTEP)*100.0);  
    }
    
    // Move particle positions & velocities
    move();
    
    // Write data to file[fi]  @file
    double *data[] = {x,y,z,vx,vy,vz,fx,fy,fz};

    for(int fi=0; fi<9; fi++ )
    {
      fprintf(file[fi], "%f\t", k*dt); // time 
      for (int i = 0; i<N; i++){ // for each molecule
        fprintf(file[fi],  "%f\t", data[fi][i]); 
      }
      fprintf(file[fi],  "\n");
    }
    
    // fprintf(file[1], "%f\t", x[i]); 
    // fprintf(file[2], "%f\t", y[i]);
    // fprintf(file[3], "%f\t", z[i]);
    // fprintf(file[4], "%f\t", vx[i]);
    // fprintf(file[5], "%f\t", vy[i]);
    // fprintf(file[6], "%f\t", vz[i]);
    // fprintf(file[7], "%f\t", fx[i]);
    // fprintf(file[8], "%f\t", fy[i]);
    // fprintf(file[9], "%f\t", fz[i]);

  //   fprintf(file_x, "%f\t", k*dt);
  //   fprintf(file_y, "%f\t", k*dt);
  //   fprintf(file_z, "%f\t", k*dt);
  //   fprintf(file_vx, "%f\t", k*dt);
  //   fprintf(file_vy, "%f\t", k*dt);
  //   fprintf(file_vz, "%f\t", k*dt);
  //   fprintf(file_fx, "%f\t", k*dt);
  //   fprintf(file_fy, "%f\t", k*dt);
  //   fprintf(file_fz, "%f\t", k*dt);

  //   for (int i = 0; i < N; i++){
  //     fprintf(file_x,  "%f\t", x[i]); 
  //     fprintf(file_y,  "%f\t", y[i]);
  //     fprintf(file_z,  "%f\t", z[i]);
  //     fprintf(file_vx, "%f\t", vx[i]);
  //     fprintf(file_vy, "%f\t", vy[i]);
  //     fprintf(file_vz, "%f\t", vz[i]);
  //     fprintf(file_fx, "%f\t", fx[i]);
  //     fprintf(file_fy, "%f\t", fy[i]);
  //     fprintf(file_fz, "%f\t", fz[i]);
  //   }
  //   fprintf(file_x, "\n");
  //   fprintf(file_y, "\n");
  //   fprintf(file_z, "\n");
  //   fprintf(file_vx, "\n");
  //   fprintf(file_vy, "\n");
  //   fprintf(file_vz, "\n");
  //   fprintf(file_fx, "\n");
  //   fprintf(file_fy, "\n");
  //   fprintf(file_fz, "\n");

  //   // Calculate Temperature 
  //   temperature();

  //   fprintf(file_T, "%f \t %f \n", k*dt, T);
  }

  // fclose(file_x);
  // fclose(file_y);
  // fclose(file_z);
  // fclose(file_vx);
  // fclose(file_vy);
  // fclose(file_vz);
  // fclose(file_fx);
  // fclose(file_fy);
  // fclose(file_fz);
  // fclose(file_T);
  
  // close files  @file
  for(int fi=0; fi<9; fi++ ){
    fclose(file[fi]);
  }

  diffusion_coefficient( &D, &rms, &r20 );

  printf("Temperature: T = %f \t r2 = %f \t rms = %f \t Diffusion: D =  %f \n", T, r20/N, sqrt(rms), D);

  // Code timing:end
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("Runtime: %f seconds \n",time_spent);

  return 0;
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

      xmin = ((x[i]-x[j]) - L*round((x[i]-x[j])/L));
      ymin = ((y[i]-y[j]) - L*round((y[i]-y[j])/L));
      zmin = ((z[i]-z[j]) - L*round((z[i]-z[j])/L));

      r2 = xmin*xmin + ymin*ymin + zmin*zmin;
      
      if(r2 < rcut2){

        f = (48./(pow(r2,7)) - 24./(pow(r2,4))); //LJ-potential

        fx[i] += f*xmin; // @ask:mugdha @mugdha:ask: if the functions cannot affect the outer-environment then how can these fx[i]s be affected
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
  // Molecules' positions (x,y,z), then velocities (vx,vy,vz) are moved according to forces

  // MOVE_R: change molecules' positions according to forces 
  for(int i = 0; i < N; i++){
    x[i] += dt*vx[i] + halfdtdt*fx[i]; //updated positions
    y[i] += dt*vy[i] + halfdtdt*fy[i];
    z[i] += dt*vz[i] + halfdtdt*fz[i];
    
    fx0[i] = fx[i]; // fx0 stores forces at f(t)
    fy0[i] = fy[i]; // fy0 stores forces at f(t)
    fz0[i] = fz[i]; // fz0 stores forces at f(t)
  }
  
  force();          // fx, fy and fz stores forces for f(t+1)
  
  // MOVE_V: change molecules' velocities according to forces
  for(int i = 0; i < N; i++){
    vx[i] += halfdt*(fx0[i] + fx[i]); //updated velocities
    vy[i] += halfdt*(fy0[i] + fy[i]);
    vz[i] += halfdt*(fz0[i] + fz[i]);
  }
}

void temperature()
{
  // Temperatures calculated over all molecules
  int v2 = 0.;
  for(int i = 0; i < N; i++){
    v2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  T = 1./3.*v2/N;
}

void diffusion_coefficient(double *D, double *rms, double *r20)
{
  // Diffusion coefficients calculated over all molecules
  double r2 = 0;
  *r20 = 0;

  for(int i = 0; i < N; i++){
    r2  += x[i]*x[i]   + y[i]*y[i]   + z[i]*z[i];  
    *r20 += xi[i]*xi[i] + yi[i]*yi[i] + zi[i]*zi[i];
  }
  *rms = (r2 - *r20)/N;  
  *D   = *rms/(6.*NSTEP*dt); 

}