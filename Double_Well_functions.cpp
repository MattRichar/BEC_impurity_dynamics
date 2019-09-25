#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <cmath>
#include <fstream>
#include <complex>
#include "Schrodinger.hpp"
#include <valarray>

using namespace std;

//If I want to change parameter values, I need to do it in both the RK4 function and the main program.
//This is the function for the RK4 algorithm if it was real. I want to keep this since it compiles with the real RK4 step function created. 
valarray<double> f(valarray<double> X,double t )
{
   valarray<double> F( X.size() );

   int N=(X.size()-2)/2;
   int i,j;
   double Dn=0.0;
   double v=1.0,tilt=7.0,J=0.5;

   for (i=0;i<2*N+2;i++){
     if (i==0){
       Dn=N/2-i;
       F[i]=((v*t + tilt)*Dn + tilt/2)*X[i]-J/2*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn-1)),0.5)*X[i+1]-J*X[i+N+1];
     }
     
     if (i<N+1 && i>0){
       Dn=N/2-i;
       F[i]=((v*t + tilt)*Dn + tilt/2)*X[i]-J/2*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn-1)),0.5)*X[i+1]-J/2*pow((N + 2*(Dn+1))*(N - 2*(Dn+1) + 2),0.5)*X[i-1]-J*X[i+N+1];
     }
     
     if (i>=N+1 && i<2*N+1){
       Dn=N/2-i+N+1;
       F[i]=((-v*t + tilt)*Dn - tilt/2)*X[i]-J/2*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn - 1)),0.5)*X[i+1]-J/2*pow((N + 2*(Dn + 1))*(N - 2*(Dn + 1) + 2),0.5)*X[i-1]-J*X[i-(N+1)];
     }
     
     if (i==2*N+1){
       Dn=N/2-i+N+1;
       F[i]=((-v*t + tilt)*Dn - tilt/2)*X[i]-J/2*pow((N + 2*(Dn + 1))*(N - 2*(Dn + 1) + 2),0.5)*X[i-1]-J*X[i-(N+1)];
     }

   }
   
   return F;
}

//This is the right hand side of the RK4 algorithm to the double-well Hamiltonian i.e. -iH(t)
valarray<complex<double>> g(valarray<complex<double>> X, double t ,double v, double tilt, double J)
{
  valarray<complex<double>> F( X.size() );

   int N=(X.size()-2)/2;
   int i,j;
   double Dn=0.0;
   complex<double> I=0.0+1.0i,n=N+0.0*i;

   for (i=0;i<2*N+2;i++){
     if (i==0){
       Dn=N/2-i;
       F[i]=-I*( ((v*t + tilt)*Dn*pow(J*n,-1) + 0.5*tilt*pow(J,-1))*(X[i])-0.5*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn-1)),0.5)*X[i+1]*pow(n,-1) -X[i+N+1] );
     }
     
     if (i<N+1 && i>0){
       Dn=N/2-i;
       F[i]=-I*( ((v*t + tilt)*Dn*pow(J*n,-1) + 0.5*tilt*pow(J,-1))*(X[i])-0.5*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn-1)),0.5)*X[i+1]*pow(n,-1) - 0.5*pow((N + 2*(Dn+1))*(N - 2*(Dn+1) + 2),0.5)*X[i-1]*pow(n,-1) -X[i+N+1] );
     }
     
     if (i>=N+1 && i<2*N+1){
       Dn=N/2-i+N+1;
       F[i]=-I*( ((-v*t + tilt)*Dn*pow(J*n,-1) - 0.5*tilt*pow(J,-1))*(X[i])-0.5*pow((N + 2*(Dn - 1) + 2)*(N - 2*(Dn - 1)),0.5)*X[i+1]*pow(n,-1) - 0.5*pow((N + 2*(Dn + 1))*(N - 2*(Dn + 1) + 2),0.5)*X[i-1]*pow(n,-1) -X[i-(N+1)] );
     }
     
     if (i==2*N+1){
       Dn=N/2-i+N+1;
       F[i]=-I*( ((-v*t + tilt)*Dn*pow(J*n,-1) - 0.5*tilt*pow(J,-1))*(X[i])-0.5*pow((N + 2*(Dn + 1))*(N - 2*(Dn + 1) + 2),0.5)*X[i-1]*pow(n,-1) -X[i-(N+1)] );
     }

   }
   
   return F;
}

// void Gen_Ham(int N, double t, double v, double tilt, double J,double Hamt[N][N])
// {
//   int i,j;
//   double Dn;
  
//   for (i=0; i<2*N+2; i++){
//    for  (j=0; j<2*N+2; j++){
//        //First "Quadrant"
//        if (i<N+1 && j<N+1){
//  	 Dn=N/2-i;
//  	 if (i==j) {
//  	   Hamt[i][i]=(v*t+tilt)*Dn + tilt/2;
//  	 }
	 
//  	 if (i==j-1){
//  	   Hamt[i][j]=-J/4*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5);
//  	 }

//  	 if (i==j+1){
//  	   Hamt[i][j]=-J/4*pow((N + 2*Dn+2)*(N - 2*Dn),0.5);
//  	 }
//        }

//        // Off diagonal blocks i.e. Second and third "Quadrants"

//        if ((i>=N+1 && j<N+1)||(i<N+1 && j>=N+1)){
//  	 if (i==j+N+1) {
//  	   Hamt[i][j]=-J/2;
//  	 }
//  	 if (j==i+N+1) {
//  	   Hamt[i][j]=-J/2;
//  	 }
//        }

//        //Fourth "Quadrant"
//        if (i>=N+1 && j>=N+1){
//  	 Dn=N/2-i+N+1;
//  	 if (i==j) {
//  	   Hamt[i][i]=(-v*t+tilt)*Dn - tilt/2;
//  	 }
	 
//  	 if (i==j-1){
//  	   Hamt[i][j]=-J/4*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5);
//  	 }

//  	 if (i==j+1){
//  	   Hamt[i][j]=-J/4*pow((N + 2*Dn + 2)*(N - 2*Dn),0.5);
//  	 }
//        }
     
//    }
//   }//End of construction.


// }

