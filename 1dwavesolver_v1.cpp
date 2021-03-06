#include <stdio.h>
#include <complex.h> 
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <complex> 
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <string>

const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
const double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
const int nx = 3000;
const double dx = 1e-8;
const int nt = 10000;
const double dt = 1e-17;
const double sigma = dx*1e-1;
const double lamb0 = 500e-6;
double shift = -nx/3*dx;
double omega = 2*M_PI*c/lamb0;
const double I0 = 1e9;
double A0 = sqrt(2*I0/(eps0*c));


using namespace std;

complex<double> f(double x, complex<double> k, double s){
    return A0 * exp(-2*(x-s)*(x-s)/sigma) * cos(k*x-s);
    //return A0 * exp(-2*(x-s)*(x-s)/sigma);
    //return A0 * cos(k*x);
}
complex<double> g(double x, complex<double> k, double omega, double s){
    return -A0 * c * (4 * (x-s)/sigma * exp(-2*(x-s)*(x-s)/sigma) * cos(k*x-s) + k * exp(-2*(x-s)*(x-s)/sigma) * sin(k*x-s));
    //return 4 * A0 * c * (x-s)/sigma * exp(-2*(x-s)*(x-s)/sigma);
    //return -A0 * omega  * sin(k*x);
    //return A0 * omega * exp(-2*x*x/sigma) * cos(k*x); 
}

void initialconds(double *x, complex<double> *kx, complex<double> *u, complex<double> *ut, double omega){
    for(int i=0; i<nx; i++){
        x[i] = (i-nx/2)*dx;
        kx[i] = 2*M_PI/lamb0;
        u[i] = f(x[i], kx[i],shift);
        ut[i] = g(x[i], kx[i], omega,shift);
    }
}


int main(){
    double *x = (double *)malloc(nx*sizeof(double));
    complex<double> *kx = (complex<double> *)malloc(nx*sizeof(complex<double>));
    complex<double> *u = (complex<double> *)malloc(nx*sizeof(complex<double>));
    complex<double> *ut = (complex<double> *)malloc(nx*sizeof(complex<double>));
    double alpha = c * dt / dx;

    complex<double> **E = (complex<double> **)malloc(nx*sizeof(complex<double> *));
    for(int i=0; i<nx; i++){
        E[i] = (complex<double> *)malloc((nt+1)*sizeof(complex<double>));
    }

    initialconds(x, kx, u, ut, omega);

    FILE *fp;

    fp = fopen("initialdistribution.csv", "w+");
    fprintf(fp, "x,Re[k],Im[k],Re[u],Im[u]\n");
    for(int i=0; i<nx; i++){
        fprintf(fp, "%f,%f,%f,%f,%f\n",x[i],real(kx[i]),imag(kx[i]),real(u[i]),imag(u[i]));
    }
    fclose(fp);

    for(int i=0; i<nx; i++){
        complex<double> fi = f(x[i],kx[i],shift);
        complex<double> gi = g(x[i],kx[i],omega,shift);
        complex<double> fim1;
        complex<double> fip1;

        if (i == 0)
        {
            fim1 = 0.0;
            fip1 = f(x[i+1],kx[i+1],shift);
        }
        else if (i == nx - 1)
        {
            fim1 = f(x[i-1],kx[i-1],shift);
            fip1 = 0.0;
        }
        else 
        {
            fim1 = f(x[i-1],kx[i-1],shift);
            fip1 = f(x[i+1],kx[i+1],shift);
        }
        E[i][1] = dt * gi + (1 - alpha * alpha) * fi + 0.5 * alpha * alpha * (fim1 + fip1);
        E[i][0] = E[i][1] - 2 * dt * gi;
    }
    
    
    for(int j=2; j<nt; j++)
    {
        for(int i=0; i<nx; i++)
        {
            if (i == 0)
            {
                E[i][j] = -E[i][j-2] + 2 * (1 - alpha * alpha) * E[i][j-1] + alpha * alpha * E[i+1][j-1];
            }
            else if (i == nx - 1)
            {
                E[i][j] = -E[i][j-2] + 2 * (1 - alpha * alpha) * E[i][j-1] + alpha * alpha * E[i-1][j-1];
            }
            else
            {
                E[i][j] = -E[i][j-2] + 2 * (1 - alpha * alpha) * E[i][j-1] + alpha * alpha * (E[i+1][j-1] + E[i-1][j-1]);
            }
            // std:cout << E[i][j] << "\n";     
        }
    }


    
    for(int j=0; j<nt; j=j+500)
    {
        char buffer[32];
        snprintf(buffer, sizeof(char) * 32, "data/Ep_%i.csv", j);
        fp = fopen(buffer, "w+");
        fprintf(fp, "x,Re[E],Im[E]\n");
        for(int i=0; i<nx; i++)
        { 
            fprintf(fp, "%f,%f,%f\n",x[i],real(E[i][j]),imag(E[i][j]));
        }    
        fprintf(fp, "\n");
        fclose(fp);
    }
    
    

    return 0;
}


