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
#include <iomanip>


//Anytime you use a class construction from Eigen, it needs to be declared here. 

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

using namespace std;
using namespace Eigen;

int main(int argc, char ** argv)
{
  // Call to the schedular for parameters
 if(argc<3){ cout << "not enough params" << endl; return 0; }  
 int N=atoi( argv[1] );
 int m=atoi( argv[2] );

 //Setting the parameters for the Hamiltonian and dynamics simulation.
 int i,j,Nt;
 double Dn,v,t,Dt,lambda1,W,tmax,tmin,density,Wcrit,tbefore,DW,k,q,l,uncertainty;
 double tilt=0.000,J=1.0,Wmax=4.0,Wmin,error=pow(10.0,-4);
 Wcrit=2*J;//critical value for tilt=0, J=2
 Wmin=0.0;
 //Proper notation for Eigen
 MatrixXd H(2*N+2,2*N+2),Ht(2*N+2,2*N+2),HW(2*N+2,2*N+2);
 VectorXcd psi(2*N+2), E0t(2*N+2), E1t(2*N+2), E2t(2*N+2), E3t(2*N+2), E4t(2*N+2), E5t(2*N+2), E6t(2*N+2), E7t(2*N+2), E8t(2*N+2), E9t(2*N+2), psiW(2*N+2), psibefore(2*N+2), psibefore1(2*N+2), psiHHt(2*N+2),psiHt(2*N+2), E0W(2*N+2), E0Wbefore(2*N+2);

 // psi=es.eigenvectors().col(0);
 // lambda1=es.eigenvalues()[0];

  //!!!!!!!!!!!!!!I am scaling each of the parameters for the impurity by a factor of N in order to make every time in the Hamiltonian of O(1)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  //Here I am initializing the system and finding the initial vector.
  for (i=0; i<2*N+2; i++){
   for  (j=0; j<2*N+2; j++){
       //First "Quadrant"
       if (i<N+1 && j<N+1){
 	 Dn=N/2-i;
 	 if (i==j) {
 	   H(i,i)=(Wmin+tilt)*(Dn)*pow(J*N,-1) + 0.5*tilt*pow(J,-1);
 	 }
	 
 	 if (i==j-1){
 	   H(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	 }

 	 if (i==j+1){
 	   H(i,j)=-0.5*pow((N + 2*Dn+2)*(N - 2*Dn),0.5)*pow(N,-1);
 	 }
       }

       // Off diagonal blocks i.e. Second and third "Quadrants"

       if ((i>=N+1 && j<N+1)||(i<N+1 && j>=N+1)){
 	 if (i==j+N+1) {
 	   H(i,j)=-1;
 	 }
 	 if (j==i+N+1) {
 	   H(i,j)=-1;
 	 }
       }

       //Fourth "Quadrant"
       if (i>=N+1 && j>=N+1){
 	 Dn=N/2-i+N+1;
 	 if (i==j) {
 	   H(i,i)=(-Wmin+tilt)*(Dn)*pow(J*N,-1) - 0.5*tilt*pow(J,-1);
 	 }
	 
 	 if (i==j-1){
 	   H(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	 }

 	 if (i==j+1){
 	   H(i,j)=-0.5*pow((N + 2*Dn + 2)*(N - 2*Dn),0.5)*pow(N,-1);
 	 }
       }
     
   }
  }

  for (i=0;i<=2*N+1;i++){
    for (j=0;j<=2*N+1;j++){
      cout << H(i,j) << "  " ;

    }
    cout << endl;
  }
    
  

 SelfAdjointEigenSolver<MatrixXd> es(H);
 psi=es.eigenvectors().col(0);

 cout << psi << endl;
 Dn=0.0;
 //t=tmin;

 vector<complex<double>> psiv(psi.data(), psi.data() + psi.rows()*psi.cols());
 complex<double>* ptr=&psiv[0]; //This is how I convert from vectorXd to vector
 Map<VectorXcd> psitest(ptr,2*N+2);

 valarray<complex<double>> psival(psiv.size()),psi0(psiv.size()),psiWval(psiv.size()); //This converts from vector to valarray
 copy(begin(psiv), end(psiv), begin(psival));

 psi0=psival;


 //This part of the code involves solving for the eigenvectors and eigenvalues of the matrix for specific values of W. This part only has a computation time based on the size of the system.
 ofstream eigs;
 ofstream wavefunctions;
 ofstream difs;
 //ofstream vectors;
 eigs.open ("eigenvalues.txt");
 wavefunctions.open ("wavefunctions.txt");
 difs.open ("difs.txt");
 //vectors.open ("eigvecs.txt");

 // Setting the precision for 
 eigs.precision(15);
 difs.precision(10);
 eigs.precision(15);
 difs.precision(10);
 wavefunctions.precision(10);
 //vectors.precision(10);
       

 W=0;
 DW=0.01;
 while(W<Wmax){
   E0Wbefore=E0W;
   
   if(W==Wmin) E0Wbefore=psi;
   
   for (i=0; i<2*N+2; i++){
     for  (j=0; j<2*N+2; j++){
       //First "Quadrant"
       if (i<N+1 && j<N+1){
 	 Dn=N/2-i;
 	 if (i==j) {
 	   HW(i,i)=(W+tilt)*(Dn)*pow(J*N,-1) + 0.5*tilt*pow(J,-1);
 	 }
	 
 	 if (i==j-1){
 	   HW(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	 }

 	 if (i==j+1){
 	   HW(i,j)=-0.5*pow((N + 2*Dn+2)*(N - 2*Dn),0.5)*pow(N,-1);
 	 }
       }

       // Off diagonal blocks i.e. Second and third "Quadrants"

       if ((i>=N+1 && j<N+1)||(i<N+1 && j>=N+1)){
 	 if (i==j+N+1) {
 	   HW(i,j)=-1;
 	 }
 	 if (j==i+N+1) {
 	   HW(i,j)=-1;//pow(N,0.5);
 	 }
       }

       //Fourth "Quadrant"
       if (i>=N+1 && j>=N+1){
 	 Dn=N/2-i+N+1;
 	 if (i==j) {
 	   HW(i,i)=(-W+tilt)*(Dn)*pow(J*N,-1) - 0.5*tilt*pow(J,-1);
 	 }
	 
 	 if (i==j-1){
 	   HW(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	 }

 	 if (i==j+1){
 	   HW(i,j)=-0.5*pow((N + 2*Dn + 2)*(N - 2*Dn),0.5)*pow(N,-1);
 	 }
       }
     
     }
   }

   SelfAdjointEigenSolver<MatrixXd> esW(HW);
   E0W=esW.eigenvectors().col(0);

   wavefunctions << abs(E0W(N+1)*E0W(N+1))+abs(E0W(N+2)*E0W(N+2))*pow(2,-0.5) << "  " << W << endl;
  
   eigs << esW.eigenvalues()[0]/(N*N) << "  "  << esW.eigenvalues()[1]/(N*N) << "  " << esW.eigenvalues()[2]/(N*N) << "  " << esW.eigenvalues()[3]/(N*N) << "  " << esW.eigenvalues()[4]/(N*N) << "  " << esW.eigenvalues()[5]/(N*N) << "  " <<  abs(esW.eigenvalues()[1]- esW.eigenvalues()[0]) << "  " <<  abs((E0W.dot(E0Wbefore))*(E0Wbefore.dot(E0W))) << "  " << 2.0*(1.0 - abs((E0W.dot(E0Wbefore))*(E0Wbefore.dot(E0W))))/(DW*DW) << "  " << W << endl;

   difs <<  abs(esW.eigenvalues()[1]- esW.eigenvalues()[0]) << "  " <<  abs(esW.eigenvalues()[2]- esW.eigenvalues()[0]) << "  " <<  abs(esW.eigenvalues()[3]- esW.eigenvalues()[0]) << "  " <<  abs(esW.eigenvalues()[4]- esW.eigenvalues()[0]) <<  "  " << abs(esW.eigenvalues()[5]- esW.eigenvalues()[0]) << "  " <<  abs(esW.eigenvalues()[6]- esW.eigenvalues()[0]) << "  " << W << endl;

   W+=DW;
 }

 ofstream probabilities;
 ofstream probabilitiesv;
 //ofstream accuracy;
 ofstream defects;
 //ofstream scaling;
 ofstream fidelity;
 ofstream fisher;

 // Defining the value of v*=vmin from the system size for the BEC-Impurity double well so that it scales nicely on a logarithmic scale.
 // I am multiplying by nine tenths afterwards so that I start just before the minimum value to show the finite size scaling. 
 double vmin = pow( double(N), -5.0/3.0 )*0.9;
 v = 0.00005 * pow( vmin / 0.00005, double(m)/20 );


 //Here I am creating a bunch of txt files for each specific parameter value of v (which is the quench speed).
 char probsTXT[100]; sprintf( probsTXT, "probs_N%i_v%.7f.txt", N, v);//sprintf fills in the appropriate parameter values.
 char probVTXT[100]; sprintf( probVTXT, "probV_N%i_v%.7f.txt", N, v);
 char defecTXT[100]; sprintf( defecTXT, "defec_N%i_v%.7f.txt", N, v);
 char fidelTXT[100]; sprintf( fidelTXT, "fidel_N%i_v%.7f.txt", N, v);
 char fishTXT[100];  sprintf( fishTXT, "fish_N%i_v%.7f.txt",N,v);

 probabilities.open  (probsTXT);
 probabilitiesv.open (probVTXT);
 defects.open        (defecTXT);
 fidelity.open       (fidelTXT);
 fisher.open         (fishTXT);
 
 //probabilitiesv.open ("probabilitiesv.txt");
 //accuracy.open ("accuracy.txt");
 //defects.open ("defects.txt");
 //scaling.open ("scaling.txt");
 //fidelity.open ("fidelity.txt");

 probabilitiesv.precision(15);
 fidelity.precision(15);
 //scaling.precision(10);
 defects.precision(15);
 probabilities.precision(15);
 fisher.precision(15);

 complex<double> prob1,prob2,prob3, prob4,prob5,prob6,prob7,prob8, prob9,prob10,probW1,probW2,probW3,prob1before;
 
 //W=v*t
 //Here I am applying a linear quench on W
 Dt=0.0001;//This needs must be less than 0.01, otherwise the wavefunction breaks down.
 DW=0.01;
   
   tmin=Wmin/v;
   tmax=Wmax/v;
   t=tmin;
   W=v*t;
   k=(10.0*Wmin);
   q=0;
   l=1;
   
   while(t<=tmax){

     psibefore=psi;

     if((l-1)*0.01 - v*Dt<v*t && v*t<=(l-1)*0.01 + v*Dt){
       psibefore1=psi;
       cout << psibefore[1] << "  " << v*t << "  " << v <<  endl;
       
     }


     RK4_Im_Shrod_step(Dt, t, psival, v, tilt, J);//Solving the Schrodinger equation.
     //RK4_Im_Shrod_step(Dt, W, psiWval, v, tilt, J);

     if(t==tmin) {psival=psi0;//The initial value has to be set since the RK4 works one step ahead.

     }

     for(i=0;i<2*N+2;i++){
       //psiW[i]=psiWval[i];
       psi[i]=psival[i]; //This loop converts from a valarray to vectorXcd.
     }

     //Dynamic fidelity loop. Here I am dividing the susceptibility by DW since that is v*v*Dt*Dt
     if((l)*0.01 - v*Dt<v*t && v*t<(l)*0.01 + v*Dt){
       
       fidelity << abs((psi.dot(psibefore))*(psibefore.dot(psi))) << "  " << 2.0*(1.0 -  abs((psi.dot(psibefore))*(psibefore.dot(psi))))/(DW*DW) << "  " << abs((psi.dot(psibefore1))*(psibefore1.dot(psi)))<< "  " << 2.0*(1.0 - abs((psi.dot(psibefore1))*(psibefore1.dot(psi))))/(Dt*Dt) << "  " << v*t << "  " << t << "  " << v << endl;
       //cout << "help1 "<< psi[1] << "  " << t << endl;
       l++;

     }
     

     //if(v==0.01 && t<=1.01) cout << psi << "  " << t << endl;
      

     //cout << psibefore[1] << "  " << psi[1] << "  " <<  abs((psi.dot(psibefore))*(psibefore.dot(psi))) << endl;

     if(k*0.1-error <= v*t && v*t<= k*0.1 + error){//This is allows me to pick out specific values I want to record for W
       for (i=0; i<2*N+2; i++){//Constructing the time dependent Hamiltonian
 	 for  (j=0; j<2*N+2; j++){
 	   //First "Quadrant"
 	   if (i<N+1 && j<N+1){
 	     Dn=N/2-i;
 	     if (i==j) {
 	       Ht(i,i)=(v*t+tilt)*(Dn)*pow(J*N,-1) + 0.5*tilt*pow(J,-1);
 	     }
	 
 	     if (i==j-1){
 	       Ht(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	     }

 	     if (i==j+1){
 	       Ht(i,j)=-0.5*pow((N + 2*Dn+2)*(N - 2*Dn),0.5)*pow(N,-1);
 	     }
 	   }

 	   // Off diagonal blocks i.e. Second and third "Quadrants"

 	   if ((i>=N+1 && j<N+1)||(i<N+1 && j>=N+1)){
 	     if (i==j+N+1) {
 	       Ht(i,j)=-1;//pow(N,0.5);
 	     }
 	     if (j==i+N+1) {
 	       Ht(i,j)=-1;//pow(N,0.5);
 	     }
 	   }

 	   //Fourth "Quadrant"
 	   if (i>=N+1 && j>=N+1){
 	     Dn=N/2-i+N+1;
 	     if (i==j) {
 	       Ht(i,i)=(-v*t+tilt)*(Dn)*pow(J*N,-1) - 0.5*tilt*pow(J,-1);
 	     }
	 
 	     if (i==j-1){
 	       Ht(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
 	     }

 	     if (i==j+1){
 	       Ht(i,j)=-0.5*pow((N + 2*Dn + 2)*(N - 2*Dn),0.5)*pow(N,-1);
 	     }
 	   }
     
 	 }
       } //End of time dependant Hamiltonian construction.

       SelfAdjointEigenSolver<MatrixXd> est(Ht);

       E0t=est.eigenvectors().col(0);
       E1t=est.eigenvectors().col(1);
       E2t=est.eigenvectors().col(2);
       E3t=est.eigenvectors().col(3);
       E4t=est.eigenvectors().col(4);
       E5t=est.eigenvectors().col(5);
       E6t=est.eigenvectors().col(6);
       E7t=est.eigenvectors().col(7);
       E8t=est.eigenvectors().col(8);
       E9t=est.eigenvectors().col(9);
       psiHHt = Ht*Ht*psi;
       psiHt = Ht*psi;

       uncertainty = abs( psi.dot(psiHHt)-(psi.dot(psiHt))*(psi.dot(psiHt)) );
			    

       prob1 = (psi.dot(E0t))*(E0t.dot(psi));
       prob1before = (psibefore.dot(E0t))*(E0t.dot(psibefore));
       //probW1 = (psiW.dot(E0t))*(E0t.dot(psiW));
       density = (1 - abs(prob1) );
       prob2 = (psi.dot(E1t))*(E1t.dot(psi));
       prob3 = (psi.dot(E2t))*(E2t.dot(psi));
       prob4 = (psi.dot(E3t))*(E3t.dot(psi));
       prob5 = (psi.dot(E4t))*(E4t.dot(psi));
       prob6 = (psi.dot(E5t))*(E5t.dot(psi));
       prob7 = (psi.dot(E6t))*(E6t.dot(psi));
       prob8 = (psi.dot(E7t))*(E7t.dot(psi));
       prob9 = (psi.dot(E8t))*(E8t.dot(psi));
       prob10 = (psi.dot(E9t))*(E9t.dot(psi));
				   
       
       //cout<<E0t << "  " << psi << "  " << t << "  " << v << endl;
       
       //if(v*t<=Wmax) {
 	 //}
	 
       //accuracy << abs(prob1) << "  " << t << endl;
       //if(v==0.01) cout << abs(prob1) << "  " << t << endl;

       //cout << abs(psi.dot(psi)) << "  " << t << endl;

       //if(Wmin <= v*t && v*t <= 2.0){
       defects << density << "  " << v*t << "  " << 1/v << "  " << log(density) << "  " << log(1/v) << endl;
       fisher << uncertainty << "  " << v*t << "  " << v << "  " << Dt << "  " << N << endl;
       //scaling << log(density) << "  " << log (1/v) << "  " << v*t << endl;

       probabilities << abs(prob1) << "  "  << abs(prob2) << "  " << abs(prob3) << "  "<< abs(prob4) << "  "<< abs(prob5) << "  "<< abs(prob6) << "  " << abs(prob7) << "  " << abs(prob8) << "  " << abs(prob9) << "  " << abs(prob10) << "  " << v*t << "  " << v  << endl; 
       //}

       //if(t<=200){

       //}
       //if(v*t>Wmax)break; //This only works since v<0.1;

       //cout << v*t << "  " << k << "  " <<abs(v*t-k*0.1) << "  " << t << "  " << v << endl;
       //cout << psibefore[1] << "  " << tbefore << "  " << psi[1] << "  " << t << "  " << abs(prob1-prob1before) << "  " << abs(prob1-(psi.dot(E0t))*(E0t.dot(psi))) << endl;
       k++;

     }

     //picking off values for t.
     if(q*0.1-error <= t && t <= q*0.1+error && t<=20+error){

       for (i=0; i<2*N+2; i++){//Constructing the time dependent Hamiltonian
     	 for  (j=0; j<2*N+2; j++){
     	   //First "Quadrant"
     	   if (i<N+1 && j<N+1){
     	     Dn=N/2-i;
     	     if (i==j) {
     	       Ht(i,i)=(v*t+tilt)*(Dn)*pow(J*N,-1) + 0.5*tilt*pow(J,-1);
     	     }
      
     	     if (i==j-1){
     	       Ht(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
     	     }

     	     if (i==j+1){
     	       Ht(i,j)=-0.5*pow((N + 2*Dn+2)*(N - 2*Dn),0.5)*pow(N,-1);
     	     }
     	   }

     	   // Off diagonal blocks i.e. Second and third "Quadrants"

     	   if ((i>=N+1 && j<N+1)||(i<N+1 && j>=N+1)){
     	     if (i==j+N+1) {
     	       Ht(i,j)=-1;//pow(N,0.5);
     	     }
     	     if (j==i+N+1) {
     	       Ht(i,j)=-1;//pow(N,0.5);
     	     }
     	   }

     	   //Fourth "Quadrant"
     	   if (i>=N+1 && j>=N+1){
     	     Dn=N/2-i+N+1;
     	     if (i==j) {
     	       Ht(i,i)=(-v*t+tilt)*(Dn)*pow(J*N,-1) - 0.5*tilt*pow(J,-1);
     	     }
      
     	     if (i==j-1){
     	       Ht(i,j)=-0.5*pow((N + 2*Dn)*(N - 2*Dn + 2),0.5)*pow(N,-1);
     	     }

     	     if (i==j+1){
     	       Ht(i,j)=-0.5*pow((N + 2*Dn + 2)*(N - 2*Dn),0.5)*pow(N,-1);
     	     }
     	   }
     
     	 }
       } //End of time dependant Hamiltonian construction.

       SelfAdjointEigenSolver<MatrixXd> est(Ht);

       E0t=est.eigenvectors().col(0);

       prob1 = (psi.dot(E0t))*(E0t.dot(psi));
     
       probabilitiesv << abs(prob1) << "  " << t << endl;
       q++;
     }

     tbefore=t;
     t+=Dt; //This must come after the RK4 call.
     

   }//End of while

 probabilitiesv.close();
 //accuracy.close();
 defects.close();
 //scaling.close();
 fidelity.close();
 //difs.close();
 //wavefunctions.close();
 //eigs.close();
 // vectors.close();
 fisher.close();
}

