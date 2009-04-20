#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#define Inf 1e+140


void fprintOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax){

	int i, j;
	
	FILE *fout;

	fout = fopen("OUTPUT/overMap.dat", "w");

	for (i=0; i<K; i++){		
		for (j=0; j<K; j++){
			fprintf(fout, "%lf ", OmegaMap[i][j]);
		}
		fprintf(fout, "\n");
	}

	fclose(fout);

	fout = fopen("OUTPUT/overBarMax.dat", "w");

	fprintf(fout, "%lf %lf %i %i\n", BarOmega, MaxOmega, rcMax[0], rcMax[1]);

	fclose(fout);

}



void printOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax){

	int i, j;

	printf("Map of misclassification probabilities:\n");
	for (i=0; i<K; i++){		
		for (j=0; j<K; j++){
			printf("[%i][%i]:%lf ", i, j, OmegaMap[i][j]);
		}
		printf("\n");
	}

	printf("Average Overlap: %lf\n", BarOmega);
	printf("Maximum Overlap: %lf", MaxOmega);
	printf(" (components: %i and %i)\n", rcMax[0], rcMax[1]);

}


void fprintParameters(int p, int K, double *Pi, double **Mu, double ***S){

	int i, j, k;

	FILE *fout;

	fout = fopen("OUTPUT/Pi.dat", "w");
	for (k=0; k<K; k++){		
		fprintf(fout, "%lf\n",Pi[k]);
	}
	fclose(fout);

	fout = fopen("OUTPUT/Mu.dat", "w");
	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			fprintf(fout, "%lf ", Mu[k][i]);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);

	fout = fopen("OUTPUT/LTSigma.dat", "w");
	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			for (j=0; j<(i+1); j++){
				fprintf(fout, "%lf ", S[k][i][j]);
			}
		}
		fprintf(fout, "\n");
	}
	fclose(fout);

}


void printParameters(int p, int K, double *Pi, double **Mu, double ***S){

	int i, j, k;

	printf("\nMixture parameters:\n");
	printf("Pi:\n");
	for (k=0; k<K; k++){		
		printf("%lf ",Pi[k]);
	}
	printf("\n");
	printf("Mu:\n");
	for (k=0; k<K; k++){
		printf("[%i]: ",k);
		for (i=0; i<p; i++){
			printf("%lf ", Mu[k][i]);
		}
		printf("\n");
	}
	printf("Sigma:\n");
	for (k=0; k<K; k++){
		printf("[%i]:\n",k);
		for (i=0; i<p; i++){
			printf("  ");
			for (j=0; j<p; j++){
				printf("%lf ", S[k][i][j]);
			}
			printf("\n");
		}
	}

}
