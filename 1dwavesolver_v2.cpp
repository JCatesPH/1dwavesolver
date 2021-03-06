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
const int nx = 5000;
const double dx = 1e-6;
const int nt = 10000;
const double dt = 1e-15;
const double sigma = dx*1e-1;
const double lamb0 = 5e-7;
double k0 = 2*M_PI/lamb0;
const int nk = 5;
double dk = k0*1e-1;
double shift = nx/10*dx;
double omega = 2*M_PI*c/lamb0;
const double I0 = 1e3;
double A0 = sqrt(2*I0/(eps0*c));



using namespace std;

complex<double> f(double x, complex<double> *k, double s){
    complex<double> Im (0.0,1.0);
    complex<double> sum = 0.0;
    for (int i=0; i<nk; i++){
        sum = sum + A0 * cos(k[i] * x);
    }
    return sum;
}
complex<double> g(double x, complex<double> *k, double omega, double s){
    complex<double> Im (0.0,1.0);
    complex<double> sum = 0.0;
    for (int i=0; i<nk; i++){
        sum = sum - A0 * k[i] * c * sin(k[i] * x);
    }
    return sum;
}

void initialconds(double *x, complex<double> *kx, complex<double> *u, complex<double> *ut, double omega, double s){
    for (int i=0; i<nk; i++){
        kx[i] = k0 - (nx / 2 - i) * dk;
    }
    for(int i=0; i<nx; i++){
        x[i] = (i-nx/4)*dx;
        u[i] = f(x[i], kx, shift);
        ut[i] = g(x[i], kx, omega, shift);
    }
    
}


int main(){
    double *x = (double *)malloc(nx*sizeof(double));
    complex<double> *kx = (complex<double> *)malloc(nk*sizeof(complex<double>));
    complex<double> *u = (complex<double> *)malloc(nx*sizeof(complex<double>));
    complex<double> *ut = (complex<double> *)malloc(nx*sizeof(complex<double>));
    double alpha = c * dt / dx;

    complex<double> *Et2 = (complex<double> *)malloc(nx*sizeof(complex<double>));
    complex<double> *Et1 = (complex<double> *)malloc(nx*sizeof(complex<double>));
    complex<double> *Et0 = (complex<double> *)malloc(nx*sizeof(complex<double>));

    cout << "k0 = " << k0 << "\n";
    cout << "nx = " << nx << "\n";
    cout << "dx = " << dx << "\n";
    cout << "nt = " << nt << "\n";
    cout << "dt = " << dt << "\n";

    initialconds(x, kx, u, ut, omega, shift);

    FILE *fp;

    fp = fopen("initialdistribution.csv", "w+");
    fprintf(fp, "x,Re[E],Im[E]\n");
    for(int i=0; i<nx; i++){
        fprintf(fp, "%f,%f,%f,\n",x[i],real(u[i]),imag(u[i]));
    }
    fclose(fp);

    for(int i=0; i<nx; i++){
        complex<double> fi = f(x[i],kx,shift);
        complex<double> gi = g(x[i],kx,omega,shift);
        complex<double> fim1;
        complex<double> fip1;

        if (i == 0)
        {
            fim1 = 0.0;
            fip1 = f(x[i+1],kx,shift);
        }
        else if (i == nx - 1)
        {
            fim1 = f(x[i-1],kx,shift);
            fip1 = 0.0;
        }
        else 
        {
            fim1 = f(x[i-1],kx,shift);
            fip1 = f(x[i+1],kx,shift);
        }
        Et1[i] = dt * gi + (1 - alpha * alpha) * fi + 0.5 * alpha * alpha * (fim1 + fip1);
        Et2[i] = Et1[i] - 2 * dt * gi;
    }
    
    
    for(int j=0; j<nt; j++)
    {
        for(int i=0; i<nx; i++)
        {
            if (i == 0)
            {
                Et0[i] = -Et2[i] + 2 * (1 - alpha * alpha) * Et1[i] + alpha * alpha * Et1[i+1];
            }
            else if (i == nx - 1)
            {
                Et0[i] = -Et2[i] + 2 * (1 - alpha * alpha) * Et1[i] + alpha * alpha * Et1[i-1];
            }
            else
            {
                Et0[i] = -Et2[i] + 2 * (1 - alpha * alpha) * Et1[i] + alpha * alpha * (Et1[i+1] + Et1[i-1]);
            }

            // std:cout << E[i][j] << "\n";     
        }
        for(int i=0; i<nx; i++)
        {
            Et2[i] = Et1[i];
            Et1[i] = Et0[i];
        }

        if (j % 200 == 0){
            char buffer[32];
            snprintf(buffer, sizeof(char) * 32, "data/Ep_%i.csv", j);
            fp = fopen(buffer, "w+");
            fprintf(fp, "x,Re[E],Im[E]\n");
            for(int i=0; i<nx; i++)
            { 
                fprintf(fp, "%f,%f,%f\n",x[i],real(Et0[i]),imag(Et0[i]));
            }    
            fprintf(fp, "\n");
            fclose(fp);
        }
    }
    

    return 0;
}


