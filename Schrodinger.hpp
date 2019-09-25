#include <stdio.h> // (*@ \label{code:Q2 header} @*)
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <fstream>
#include <complex>
#include <valarray>

#ifndef Schrod_H_  /* Include guard */
#define Schrod_H_

using namespace std;

// vector<double> RK4_step(int N, double h, double t,vector<double> X,  double *f(vector<double> X, double t));
valarray<double> f(valarray<double> X, double t);
void RK4_step( double h, double &t, valarray<double> &X, valarray<double> (*f)(valarray<double>X, double t));
void RK4_Im_Shrod_step(double h, double &t, valarray<complex<double>> &X, double v, double tilt, double J);
valarray<complex<double>> g(valarray<complex<double>> X, double t, double v, double tilt, double J);
//void Gen_Ham(int N, double t, double v, double tilt, double J, double Hamt[N][N]);

#endif // Schrod_H_
