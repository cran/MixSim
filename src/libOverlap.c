#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#define Inf 1e+140

#include "overlap.h"


void ComputePars(int p, int K, double *Pi, double **Mu, double ***S, double ***li, double ***di, double **const1){

	int i, j, k;

	double dtmt;

	double *m1, *m2, *Eig, *detS;
	double **Ga, **L, **Ga2, **L2, **Si;
	double ***Sinv, ***Sh;


	MAKE_VECTOR(detS, K);
	MAKE_3ARRAY(Sinv, K, p, p);
	MAKE_3ARRAY(Sh, K, p, p);
	
	MAKE_VECTOR(m1, p);
	MAKE_VECTOR(m2, p);
	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(Ga, p, p);	
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(Ga2, p, p);	
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(Si, p, p);

	for (k=0; k<K; k++){

		cpy1(S, k, p, p, Ga);
		EigValDec(p, Eig, Ga, &dtmt);

		detS[k] = dtmt;

		Anull(L, p, p);
		for (i=0; i<p; i++){
			L[i][i] = pow(Eig[i], 0.5);
		}
		XAXt2(Ga, p, L, Sh, k);

		for (i=0; i<p; i++){
			L[i][i] = 1 / Eig[i];
		}
		XAXt2(Ga, p, L, Sinv, k);

	}



	for (i=0; i<(K-1); i++){
		for (j=i+1; j<K; j++){

			cpy1(Sh, i, p, p, Ga);
			cpy1(Sinv, j, p, p, L);
			cpy1(Sinv, i, p, p, L2);
			for (k=0; k<p; k++){
				m1[k] = Mu[i][k] - Mu[j][k];
				m2[k] = -m1[k];
			}


			XAXt(Ga, p, L, Si);

			EigValDec(p, Eig, Si, &dtmt);
			for (k=0; k<p; k++){
				li[i][j][k] = Eig[k];
			}

			multiply(L2, p, p, Ga, p, p, Ga2);
			matxvec(Ga2, p, p, m1, p, Eig);
			tA(Si, p, p, Ga2);
			matxvec(Ga2, p, p, Eig, p, m1);
			for (k=0; k<p; k++){
				di[i][j][k] = m1[k];
			}



			cpy1(Sh, j, p, p, Ga2);


			XAXt(Ga2, p, L2, Si);

			EigValDec(p, Eig, Si, &dtmt);
			for (k=0; k<p; k++){
				li[j][i][k] = Eig[k];
			}

			multiply(L, p, p, Ga2, p, p, Ga);
			matxvec(Ga, p, p, m2, p, Eig);
			tA(Si, p, p, Ga);
			matxvec(Ga, p, p, Eig, p, m2);
			for (k=0; k<p; k++){
				di[j][i][k] = m2[k];
			}


			const1[i][j] = log((Pi[j]*Pi[j]) / (Pi[i]*Pi[i]) * detS[i]/detS[j]);
			const1[j][i] = -const1[i][j];


		}
	}


	FREE_VECTOR(detS);
	FREE_3ARRAY(Sinv);
	FREE_3ARRAY(Sh);
	
	FREE_VECTOR(m1);
	FREE_VECTOR(m2);
	FREE_VECTOR(Eig);
	FREE_MATRIX(Ga);
	FREE_MATRIX(L);
	FREE_MATRIX(Ga2);
	FREE_MATRIX(L2);
	FREE_MATRIX(Si);

}



void GetOmegaMap(double c, int p, int K, double ***li, double ***di, double **const1, int *fix, double *pars, int lim, double asympt, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax){


	int i, j, k;
	double Cnst1, t, s, TotalOmega, OmegaOverlap;
	double eps, acc, sigma;
	
	double *Li, *Di, *ncp, *coef, *ldprod, *const2;
	double *trace;
	int ifault;
	int *df;

	MAKE_VECTOR(Li, p);
	MAKE_VECTOR(Di, p);
	MAKE_VECTOR(coef, p);
	MAKE_VECTOR(ldprod, p);
	MAKE_VECTOR(ncp, p);
	MAKE_VECTOR(const2, p);
	MAKE_VECTOR(df, p);

	MAKE_VECTOR(trace, 7);

	sigma = 0.0;

	eps = pars[0];
	acc = pars[1];

	TotalOmega = 0.0;
	(*MaxOmega) = 0.0;

   	for (k=0; k<p; k++){
		df[k] = 1;
	}

	i = 0;
	j = 1;


	if (asympt == 0){

		while (i < (K-1)){

			if (fix[i] == 1){
				
				for (k=0; k<p; k++){
					Di[k] = di[i][j][k];
				}

				if (fix[j] == 1){
					for (k=0; k<p; k++){
						Li[k] = li[i][j][k];
					}
					Cnst1 = const1[i][j];
				} else {
			      		for (k=0; k<p; k++){
						Li[k] = li[i][j][k] / c;
					}	
					Cnst1 = const1[i][j] - p * log(c);
				}

			} else {

				for (k=0; k<p; k++){
					Di[k] = di[i][j][k] / pow(c, 0.5);
				}

				if (fix[j] == 1){
					for (k=0; k<p; k++){
						Li[k] = c * li[i][j][k];
					}
					Cnst1 = const1[i][j] + p * log(c);
				} else {
			      		for (k=0; k<p; k++){
						Li[k] = li[i][j][k];
					}	
					Cnst1 = const1[i][j];
				}

			}


			s = 0;
			for (k=0; k<p; k++){
				coef[k] = Li[k] - 1.0;
				ldprod[k] = Li[k] * Di[k];
				const2[k] = ldprod[k] * Di[k] / coef[k];
				s = s + const2[k];
				ncp[k] = pow(ldprod[k] / coef[k], 2);
			}
			t = s + Cnst1;

			OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);



			if (fix[j] == 1){
				
				for (k=0; k<p; k++){
					Di[k] = di[j][i][k];
				}

				if (fix[i] == 1){
					for (k=0; k<p; k++){
						Li[k] = li[j][i][k];
					}
					Cnst1 = const1[j][i];
				} else {
			      		for (k=0; k<p; k++){
						Li[k] = li[j][i][k] / c;
					}	
					Cnst1 = const1[j][i] - p * log(c);
				}

			} else {

				for (k=0; k<p; k++){
					Di[k] = di[j][i][k] / pow(c, 0.5);
				}

				if (fix[i] == 1){
					for (k=0; k<p; k++){
						Li[k] = c * li[j][i][k];
					}
					Cnst1 = const1[j][i] + p * log(c);
				} else {
			      		for (k=0; k<p; k++){
						Li[k] = li[j][i][k];
					}	
					Cnst1 = const1[j][i];
				}

			}

			s = 0;
			for (k=0; k<p; k++){
				coef[k] = Li[k] - 1.0;
				ldprod[k] = Li[k] * Di[k];
				const2[k] = ldprod[k] * Di[k] / coef[k];
				s = s + const2[k];
				ncp[k] = pow(ldprod[k] / coef[k], 2);
			}
			t = s + Cnst1;

			OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
	


			OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
			TotalOmega = TotalOmega + OmegaOverlap;

			if (OmegaOverlap > (*MaxOmega)){
				(*MaxOmega) = OmegaOverlap;
				rcMax[0] = i;
				rcMax[1] = j;
			}

			
			if (j < (K - 1)){
				j = j + 1;
			} else {
				i = i + 1;
				j = i + 1;
			}

		}

	}






	if (asympt == 1){

		while (i < (K-1)){

			if (fix[i] == 1){
				
				if (fix[j] == 1){
					
					s = 0;
					for (k=0; k<p; k++){
						Di[k] = di[i][j][k];
						Li[k] = li[i][j][k];
						coef[k] = Li[k] - 1.0;
						ldprod[k] = Li[k] * Di[k];
						const2[k] = ldprod[k] * Di[k] / coef[k];
						s = s + const2[k];
						ncp[k] = pow(ldprod[k] / coef[k], 2);
					}
					t = s + const1[i][j];

					OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

				} else {

					OmegaMap[i][j] = 0.0;

				}

			} else {

				if (fix[j] == 1){

					OmegaMap[i][j] = 0.0;

				} else {
			      		for (k=0; k<p; k++){
						coef[k] = li[i][j][k] - 1.0;
						ncp[k] = 0.0;
					}
					t = const1[i][j];

					OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);					

				}
			}




			if (fix[j] == 1){
				
				if (fix[i] == 1){
					
					s = 0;
					for (k=0; k<p; k++){
						Di[k] = di[j][i][k];
						Li[k] = li[j][i][k];
						coef[k] = Li[k] - 1.0;
						ldprod[k] = Li[k] * Di[k];
						const2[k] = ldprod[k] * Di[k] / coef[k];
						s = s + const2[k];
						ncp[k] = pow(ldprod[k] / coef[k], 2);
					}
					t = s + const1[j][i];
				
					OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

				} else {

					OmegaMap[j][i] = 0.0;

				}

			} else {

				if (fix[i] == 1){

					OmegaMap[j][i] = 0.0;

				} else {
			      		for (k=0; k<p; k++){
						coef[k] = li[j][i][k] - 1.0;
						ncp[k] = 0.0;
					}
					t = const1[j][i];
					
					OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

				}
			}


			OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
			TotalOmega = TotalOmega + OmegaOverlap;

			if (OmegaOverlap > (*MaxOmega)){
				(*MaxOmega) = OmegaOverlap;
				rcMax[0] = i;
				rcMax[1] = j;
			}

			
			if (j < (K - 1)){
				j = j + 1;
			} else {
				i = i + 1;
				j = i + 1;
			}

		}

	}

	(*BarOmega) = TotalOmega / (K * (K - 1) / 2.0);

	for (k=0; k<K; k++){
		OmegaMap[k][k] = 1.0;
	}


	FREE_VECTOR(Li);
	FREE_VECTOR(Di);
	FREE_VECTOR(coef);
	FREE_VECTOR(ldprod);
	FREE_VECTOR(ncp);
	FREE_VECTOR(const2);
	FREE_VECTOR(df);

	FREE_VECTOR(trace);
	
}


void ExactOverlap(int p, int K, double *Pi, double **Mu, double ***S, double *pars, int lim, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax){

	double c, Balpha, Malpha;
	int asympt;

	int  *fix;
	double **const1;
	double ***li, ***di;

	MAKE_VECTOR(fix, K);

	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);

	ComputePars(p, K, Pi, Mu, S, li, di, const1);

	anulli(fix, K);

	c = 1.0;
	asympt = 0;
	GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);

	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;

	FREE_VECTOR(fix);

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);

}


/* FIND MULTIPLIER C ON THE INTERVAL (lower, upper) */

void FindC(double lower, double upper, double Omega, int method, int p, int K, double ***li, double ***di, double **const1, int *fix, double *pars, int lim, double (*c), double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax){

	double diff, eps, acc;
	int sch, asympt, stopIter;

	eps = pars[0];
	acc = pars[1];

	diff = Inf;
	stopIter = 1000;

	sch = 0;

	while (fabs(diff) > eps){

		(*c) = (lower + upper) / 2.0;

		asympt = 0;
		GetOmegaMap((*c), p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &(*BarOmega), &(*MaxOmega), rcMax);
	
		if (method == 0){

			if ((*BarOmega) < Omega){ /* clusters are too far */
				lower = (*c);
			} else {
				upper = (*c);
			}
			
			diff = (*BarOmega) - Omega;

		} else {

			if ((*MaxOmega) < Omega){ /* clusters are too far */
				lower = (*c);
			} else {
				upper = (*c);
			}
			
			diff = (*MaxOmega) - Omega;			

		}

		sch = sch + 1;

		if (sch == stopIter){
			(*c) = 0.0;
/*			printf("Error: required overlap was not reached in %i iterations...\n", stopIter); */
			break;
		}

	}

}



void OmegaClust(double Omega, int method, int p, int K, double PiLow, double Ubound, double emax, double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail)){

	int asympt, sch;
	double c, diff, lower, upper, eps, acc, Balpha, Malpha;

	int *fix;
	double **const1;
	double ***li, ***di;


	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);

	MAKE_VECTOR(fix, K);

	anulli(fix, K);

	eps = pars[0];
	acc = pars[1];

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);
		

	sch = 0;

	do{

		(*fail) = 0;

		/* generate parameters */
/*		printf("Simulating dataset...\n"); */

		genPi(K, PiLow, Pi);

		genMu(p, K, Mu, Ubound);
		if (sph == 0){
			genSigmaEcc(p, K, emax, S);
		} else {
			genSphSigma(p, K, S);
		}

		/* prepare parameters */
		
		ComputePars(p, K, Pi, Mu, S, li, di, const1);

		/* check if desired overlap is reachable */

		asympt = 1;
		c = 0.0;
		GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);

		if (method == 0){
			diff = Balpha - Omega;
		} else {
			diff = Malpha - Omega;
		}

		if (diff < -eps){ /* not reachable */
/*			printf("Warning: the desired overlap cannot be reached...\n"); */
			(*fail) = 1;
		} else {

			lower = 0.0;
			upper = pow(2.0, 2);

			do{
				FindC(lower, upper, Omega, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, rcMax);
				
				lower = upper;
				upper = upper * upper;
				if (upper > 1000000){
/*					printf("Warning: the desired overlap cannot be reached...\n");
					printf("Simulating another dataset...\n"); */
					(*fail) = 1;			
					break;
				} /* (upper > 1000000) prevents nonstopping loops */

			} while(c == 0);

		}


		if ((*fail) == 0){

			/* correct covariances by multiplier C */

			cxS(p, K, S, c);

			break;

		}

		sch = sch + 1;

		if (sch == resN){
			printf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
			(*fail) = 1;
			break;
		}

	} while ((*fail) == 1);

	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);

	FREE_VECTOR(fix);


}





void OmegaBarOmegaMax(int p, int K, double PiLow, double Ubound, double emax, double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail)){


	int i, j, k, asympt, sch, rowN, colN, method;
	double c, diff, lower = 0, upper = 0, eps, acc, Balpha, Malpha;

	int *fix, *fix2;
	double **const1, **const12, **OmegaMap2;
	double ***li, ***di, ***li2, ***di2;

	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);
	MAKE_3ARRAY(li2, 2, 2, p);
	MAKE_3ARRAY(di2, 2, 2, p);
	MAKE_MATRIX(const12, 2, 2);
	MAKE_MATRIX(OmegaMap2, 2, 2);

	MAKE_VECTOR(fix, K);
	MAKE_VECTOR(fix2, 2);

	anulli(fix, K);
	anulli(fix2, 2);

	eps = pars[0];
	acc = pars[1];

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);

	(*fail) = 1;

	if ((Malpha < Balpha) | (Malpha > Balpha * K * (K - 1) / 2.0)){ /* wrong parameters*/

		printf("Error: incorrect values of average and maximum overlaps...\n");

	} else {


		sch = 0;
	
		do{

			/* generate parameters */
/*			printf("Simulating dataset...\n");	*/

			genPi(K, PiLow, Pi);
			genMu(p, K, Mu, Ubound);
			if (sph == 0){
				genSigmaEcc(p, K, emax, S);
			} else {
				genSphSigma(p, K, S);
			}

			/* prepare parameters */
		
			ComputePars(p, K, Pi, Mu, S, li, di, const1);

			/* check if maximum overlap is reachable */

			asympt = 1;
			c = 0.0;
			GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);

			diff = Malpha - (*MaxOmega);

			if (diff >= -eps){ /* reachable */

				lower = 0.0;
				upper = pow(2.0, 10);

				do{

					/* find C for two currently largest clusters */
					
					rowN = rcMax[0];
					colN = rcMax[1];

					for (i=0; i<2; i++){
						for (j=0; j<2; j++){
							for (k=0; k<p; k++){
								li2[i][j][k] = li[rcMax[i]][rcMax[j]][k];
								di2[i][j][k] = di[rcMax[i]][rcMax[j]][k];
							}
							const12[i][j] = const1[rcMax[i]][rcMax[j]];
						}
					}				

					Malpha = (*MaxOmega);
					method = 1;
					FindC(lower, upper, Malpha, method, p, 2, li2, di2, const12, fix2, pars, lim, &c, OmegaMap2, &Balpha, &Malpha, rcMax);

					if (c == 0){ /* abnormal termination */
/*						printf("Warning: the desired overlap cannot be reached...\n"); */
						(*fail) = 1;
						break;
					}

					asympt = 0;
					GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);
					upper = c;

			       		diff = Balpha - (*BarOmega);
					if (diff < -eps){ /* BarOmega is not reachable */
/*						printf("Warning: the desired overlap cannot be reached...\n"); */
						(*fail) = 1;
						break;
					}

					diff = Malpha - (*MaxOmega);
					if (diff < eps){ /* MaxOmega has been reached */
						(*fail) = 0;
						break;
					}

				} while ((*fail) != 0);

			}	



			if ((*fail) == 0){ /* OmegaMax is reached and OmegaBar is reachable */

				/* correct covariances by multiplier C */

				cxS(p, K, S, c);

				ComputePars(p, K, Pi, Mu, S, li, di, const1);

				fix[rcMax[0]] = 1;
				fix[rcMax[1]] = 1;
				upper = 1;

				Balpha = (*BarOmega);
				method = 0;
				FindC(lower, upper, Balpha, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, rcMax);

				/* correct covariances by multiplier C */

				for (k=0; k<K; k++){
					for (i=0; i<p; i++){
						for (j=0; j<p; j++){
							if (fix[k] == 0){
								S[k][i][j] = c * S[k][i][j];
							}
						}
					}
				}

				break;

			}

			
			sch = sch + 1;

			if (sch == resN){
				printf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
				(*fail) = 1;
				break;
			}

		} while ((*fail) == 1);

	}

	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);
	FREE_3ARRAY(li2);
	FREE_3ARRAY(di2);
	FREE_MATRIX(const12);
	FREE_MATRIX(OmegaMap2);

	FREE_VECTOR(fix);
	FREE_VECTOR(fix2);	


}
