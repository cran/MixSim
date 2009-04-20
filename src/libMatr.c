#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#include "overlap.h"


/* Multiplies matrices a and b and puts the result in c which should be
 pre-allocated.   exit() will be called if a and b are incompatible
*/

void multiply(double **a, int arows, int acols,
	      double **b, int brows, int bcols, double **c)
{
  int i, j, k;
  
  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[i][j] = 0;
      for (k=0; k<acols; k++)
	c[i][j] += a[i][k] * b[k][j];
    }

}

void multiply2(double **a, int arows, int acols,
	       double **b, int brows, int bcols, double ***c, int m)
{
  int i, j, k;
  
  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[m][i][j] = 0;
      for (k=0; k<acols; k++)
	c[m][i][j] += a[i][k] * b[k][j];
    }

}

/* Multiplies matrix a and vector x and puts the result in y which should be
 pre-allocated.   exit() will be called if a and x are incompatible
*/
void matxvec(double **a, int arows, int acols,
		double *x, int xrows, double *y)
{
  int i, k;
  
  for (i=0; i<arows; i++){
    y[i] = 0;
    for (k=0; k<acols; k++){
      y[i] += a[i][k] * x[k];
    }
  }

}


/*copy matrix A to matrix B*/
void cpy(double **a, int nrows, int ncols, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j];
    }
  }

}


/*copy matrix A[k] to matrix B */
void cpy1(double ***a, int k, int nrows, int ncols, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[k][i][j];
    }
  }

}

/*copy matrix A to matrix B[k] */
void cpy2(double **a, int nrows, int ncols, double ***b, int k)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[k][i][j]=a[i][j];
    }
  }

}


int vecMin(double *x, int p, double (*min)){

	int i, minN;

	(*min) = x[0];
	minN = 0;

	for (i=0;i<p;i++){
		if (x[i] < (*min)){
			(*min) = x[i];
			minN = i;
		}
	}

	return minN;

}


int vecMax(double *x, int p, double (*max)){

	int i, maxN;

	(*max) = x[0];
	maxN = 0;

	for (i=0;i<p;i++){
		if (x[i] > (*max)){
			(*max) = x[i];
			maxN = i;
		}
	}

	return maxN;

}




int vec11vecSQ(double *y, int p, double **Res){
	
	int i,j;

	for (i=0;i<p;i++){
		for (j=0;j<p;j++){
			Res[i][j] = y[i]*y[j];
		}
	}
	
	return 0;
}


double vecNNvec(int p, double *y, double *x){
	
	int i;
	double Res;

	Res = 0;
	for (i=0;i<p;i++){
		Res = Res + y[i]*x[i];
	}

	return Res;
}


int mat_(int a, int b,double **Res, double **Y){
	
	int i,j;

	for (i=0;i<a;i++){
		for (j=0;j<b;j++){
			Res[i][j] = Res[i][j] - Y[i][j];
		}
	}
	
	return 0;
}



int vecsum(int a, int b,double **OO, double *Res){
	
	int i,j;

	for (i=0;i<a;i++){
		Res[i] = 0;
		for (j=0;j<b;j++){
			Res[i] = Res[i] + OO[i][j];
		}
	}
	
	return 0;
}


int MatrixProd(double **OO, int p, int m, double **Res){
     
     int i,j,k;

     for (i=0; i<p; i++){
         for (j=0; j<p; j++){
             Res[i][j]=0;
             for (k=0; k<m; k++){
                 Res[i][j]=Res[i][j]+OO[i][k]*OO[j][k];
             }
         }     
     }
     
     return 0;
}


int Kronecker(double **A, int a1, int a2, double **B, int b1, int b2, double **Res){

  int inda1, inda2, indb1, indb2, indRes1, indRes2;
  int i;
  int n;


  n = a1 * b1 * a2 * b2;

  indRes1 = 0;
  indRes2 = -1;

  inda1 = 0;
  inda2 = 0;
  indb1 = 0;
  indb2 = -1;

  for (i=0; i<n; i++){

    indb2++;
    indRes2++;

    if (indb2 == b2){

      indb2 = 0;
      inda2++;

      if (inda2 == a2){
	
	inda2 = 0;
	indb1++;
	indRes1++;
	indRes2 = 0;

	if (indb1 == b1){

	  indb1 = 0;
	  inda1++;

	}

      }
    }

    Res[indRes1][indRes2] = A[inda1][inda2] * B[indb1][indb2];

  }

  return 0;

}


int Gmat(int p, int m, double **Res){
     
     int a,b,i,i1,i2,n,ind;

	 n = 0;

     for (a=0; a<p; a++){
         for (b=0; b<p; b++){
         	
         	if (a < b){
         		i1 = b;
         		i2 = a;
         	} else {
         		i1 = a;
         		i2 = b;
         	}
         	
         	ind = m - (p - i2) * (p - i2 + 1) / 2 + i1 - i2;
         	
         	for (i=0; i<m; i++){
         	
         		if (i != ind ){
         			Res[n][i] = 0;	
         		} else {
         			Res[n][i] = 1;
         		}
         		
         	}
         	
         	n++;

         }     
     }
     
     return 0;
}


void tA(double **A, int a, int b, double **Res){

	int i,j;

   	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			Res[i][j] = A[j][i];
		}
	}
	
}


int ZXY(double **Z, int az, int bz, double **X, int ax, int bx, double **Y, int ay, int by, double **Res){

	double **Res1;

	MAKE_MATRIX(Res1, az, bx);	

	multiply(Z, az, bz, X, ax, bx, Res1);
	multiply(Res1, az, bx, Y, ay, by, Res);

	FREE_MATRIX(Res1);
 
        return 0;
    
}



void XAXt(double **X, int p, double **A, double **Res){

	double **Res1, **Res2;

	MAKE_MATRIX(Res1, p, p);	
	MAKE_MATRIX(Res2, p, p);

	tA(X, p, p, Res2);

	multiply(X, p, p, A, p, p, Res1);
	multiply(Res1, p, p, Res2, p, p, Res);

	FREE_MATRIX(Res1);
 	FREE_MATRIX(Res2);

}


/* Computes X %*% A %*% t(X) and writes results into Res[k,,] */

void XAXt2(double **X, int p, double **A, double ***Res, int k){

	double **Res1, **Res2;

	MAKE_MATRIX(Res1, p, p);	
	MAKE_MATRIX(Res2, p, p);

	tA(X, p, p, Res2);

	multiply(X, p, p, A, p, p, Res1);
	multiply2(Res1, p, p, Res2, p, p, Res, k);

	FREE_MATRIX(Res1);
 	FREE_MATRIX(Res2);

}


void Anull(double **X, int ax, int bx){
     
     int i, j;

     for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		X[i][j] = 0.0;
	 }
     }
}


void anull(double *x, int p){
     
     int i;

     for (i=0; i<p; i++){
	     x[i] = 0.0;
     }
}


void Anulli(int **X, int ax, int bx){
     
     int i, j;

     for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		X[i][j] = 0;
	 }
     }
}


void anulli(int *x, int p){
     
     int i;

     for (i=0; i<p; i++){
	     x[i] = 0;
     }
}


int asvector(double **X, int ax, int bx, double *ResVec){
     
    int i,j,k;

    k = 0;

    for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		 	ResVec[k] = X[i][j];
		 	k++;
         }
     }

    return 0;
    
}


void cxS(int p, int K, double ***S, double c){

	int i, j, k;

	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			for (j=0; j<p; j++){
				S[k][i][j] = c * S[k][i][j];
			}
		}
	}

}
