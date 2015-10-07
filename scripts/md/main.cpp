#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>

using namespace std;

double randGen();

int main()
{

    // random seed
    srand((unsigned)time(NULL)*(unsigned)time(NULL));

    // Initialize parameters
    int n       = 1;
    int NSTEP   = 10;
    int vmax    = 2.7;
    int rcut    = 3;
    int L       = 3;
    int dt      = 0.001;

    int max_n = int(pow(n,1./3.));

    // x,y,z coordinates of particles in a 3D lattice
    double x[n][NSTEP];
    double y[n][NSTEP];
    double z[n][NSTEP];

    // filling coordinate lattice with 0s
    for(int i=0;i<n;i++){
        for(int j=0;j<NSTEP;j++){
            x[i][j]=0;
            y[i][j]=0;
            z[i][j]=0;
            //cout << x[i][j] << " ";
        }
        //cout << endl;
     }

    // Velocity array
    double vx[n][NSTEP];
    double vy[n][NSTEP];
    double vz[n][NSTEP];

    // fill velocity arrays with
    for(int i=0;i<n;i++){
        vx[i][0] = vmax*(2*randGen()-1);
        cout << vx[i][0] << endl;
        vy[i][0] = vmax*(2*randGen()-1);
        cout << vx[i][0] << endl;
        vz[i][0] = vmax*(2*randGen()-1);
        cout << vx[i][0] << endl;
    }

    // subtract centre of mass velocity
//    vxcm = vx[0:,0].mean()
//    vycm = vy[0:,0].mean()
//    vzcm = vz[0:,0].mean()

    int avg_vx=0;
    int avg_vy=0;
    int avg_vz=0;

    for(int i=0;i<=n;i++){

        avg_vx += vx[i][0];
        avg_vx += vy[i][0];
        avg_vx += vz[i][0];

        // cout << "new iteration" << endl;
        // cout << avg_vx << endl;
        // cout << avg_vy << endl;
        // cout << avg_vz << endl;
    }


    return 0;
}

// Print an array
//    for(int i=0;i<2;i++){
//        for(int j=0;j<3;j++){
//            cout << x[i][j] << " ";
//        }
//        cout << endl;
//    }

double randGen(){
    // random number generator [0,1]
    double r = (double)rand() / RAND_MAX;
    return r;
}


