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

//Anytime you use a class construction from Eigen, it needs to be declared here.

void RK4_step( double h, double &t, valarray<double> &X, valarray<double> (*f)(valarray<double>X, double t))
{
  int N = X.size();
  valarray<double> k1(N),k2(N),k3(N),k4(N);

  k1=h*(f(X,t));

  k2=h*(f(X+0.5*k1,t+0.5*h));

  k3=h*(f(X+0.5*k2,t+0.5*h));

  k4=h*(f(X+k3,t+h));

  X+=(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);

}

//This RK4 solves the Shrodinger equation for the double-well system. I introduce the parameters of v, tilt, and J. If I want to generalize this to some function, I need to adopt the notation for an introduced function such as in RK4_step. 
void RK4_Im_Shrod_step(double h, double &t, valarray<complex<double>> &X, double v, double tilt, double J)
{

  int N = X.size();
  valarray<complex<double>> k1(N),k2(N),k3(N),k4(N);
  complex<double> I=0.0+1.0i, hI=h+0.0i, a=0.5+0.0i, b=1.0/6.0+0.0i, c=2.0+0.0i;  //For some reason, valarray gets angry at me unless I put the coefficients in the RK4 in a complex form for the complex arguments instead of explicity writing it out as a complex number. 

  k1=hI*(g(X,t,v,tilt,J));

  k2=hI*(g(X+a*k1,t+0.5*h,v,tilt,J));

  k3=hI*(g(X+a*k2,t+0.5*h,v,tilt,J));

  k4=hI*(g(X+k3,t+0.5*h,v,tilt,J));

  X+=b*(k1+c*k2+c*k3+k4);

}

