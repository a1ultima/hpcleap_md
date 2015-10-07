#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>


using namespace std;

double randGen();

int main()
{

    double moo=pow(9,1./3.);

    cout << moo << endl;

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


