#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"
/*
#define MATHLIB_STANDALONE 0 */ /*It is essential to have this before the call 
                               to the Rmath's header file because this decides
                               the definitions to be set. */

#include<Rmath.h>

#include "overlap.h"

/*
void EigValDec(int size, double *W, double **A, double (*determinant));
int vecMin(double *x, int p, double (*min));
int vecMax(double *x, int p, double (*max));
void Anull(double **X, int ax, int bx);
void anull(double *X, int p);
void XAXt(double **X, int p, double **A, double **Res);
void cpy(double **a, int nrows, int ncols, double **b);
void cpy2(double **a, int nrows, int ncols, double ***b, int k);
*/


/* genSigma :
   generates covariance matrix based on (p + 1) observations (unstable covariance matrix) */

void genSigma(int p, double **VC){

	int i,j,k,n;
	double **x, *mu;

	n = p + 1;

	MAKE_MATRIX(x, n, p);
	MAKE_VECTOR(mu, p);

	anull(mu, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			x[i][j] = rnorm(0.0, 1.0);
			mu[j] = mu[j] + x[i][j];
		}
        }

	for (j=0;j<p;j++){
		mu[j] = mu[j] / n;
	}	
	
	Anull(VC, p, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			for (k=0; k<p; k++){
				VC[j][k] = VC[j][k] + (x[i][j] - mu[j]) * (x[i][k] - mu[k]);
			}
		}
        }

	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			VC[j][k] = VC[j][k] / (n - 1);
		}
	}

	FREE_MATRIX(x);
	FREE_VECTOR(mu);

}


/* genSigmaEcc :
   generates covariance matrix  with prespecified eccentricity */

void genSigmaEcc(int p, int K, double emax, double ***S){

	int i, k;

	double dtmt, minL, maxL, e;
	double *Eig;
	double **VC, **L, **R;

	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(VC, p, p);
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(R, p, p);

	for (k=0; k<K; k++){
		genSigma(p, VC);
		cpy2(VC, p, p, S, k);		

		EigValDec(p, Eig, VC, &dtmt);

		i = vecMin(Eig, p, &minL);
		i = vecMax(Eig, p, &maxL);

		e = pow(1 - minL / maxL, 0.5);


		if (e > emax){

			Anull(L, p, p);

			for (i=0; i<p; i++){
				Eig[i] = maxL * (1 - emax * emax * (maxL - Eig[i]) / (maxL - minL));
				L[i][i] = Eig[i];
			}

			XAXt(VC, p, L, R);
			cpy2(R, p, p, S, k);

		}

	}

	FREE_MATRIX(VC);
	FREE_MATRIX(L);
	FREE_MATRIX(R);

}



/* genSphSigma :
   generates spherical covariance matrix */

void genSphSigma(int p, int K, double ***S){

	int i, k;
	double r;
	double **L;

	MAKE_MATRIX(L, p, p);		

	Anull(L, p, p);
	
	for (k=0; k<K; k++){

		r = runif(0.0, 1.0);
		for (i=0; i<p; i++){			
			L[i][i] = r;
		}

		cpy2(L, p, p, S, k);
	
	}

	FREE_MATRIX(L);

}



/* genSphSigma :
   generates matrix of means */

void genMu(int p, int K, double **Mu, double Ubound){
		
	int i, k;
	
	if (Ubound <= 0) Ubound = 1.0;
	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
		
			Mu[k][i] = runif(0.0, Ubound);

		}
	}

}



/* genPi :
   generates mixing proportions */

void genPi(int K, double PiLow, double *Pi){

	int flag, k;
	double s;

	flag = 0;


	if ((PiLow >= 1) | (PiLow <= 0)){
/*		printf("Warning: PiLow is out of range... generated equal mixing proportions...\n"); */
		for (k=0; k<K; k++){
			Pi[k] = 1.0 / K;
		}
	} else {
		s = 0.0;
		for (k=0; k<K; k++){
			Pi[k] = rgamma(1.0, 1.0);
			s += Pi[k];
		}
		for (k=0; k<K; k++){
			Pi[k] = PiLow + Pi[k] / s * (1 - K * PiLow);
			if (Pi[k] < PiLow){
				flag = 1;
				break;
			}
		}
		if (flag == 1){
/*			printf("Warning: PiLow is too high... generated equal mixing proportions...\n"); */
			for (k=0; k<K; k++){
				Pi[k] = 1.0 / K;
			}
		}
		
		
	}

}
