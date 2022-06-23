/*
	Bayesian Functional GWAS --- MCMC (bfGWAS:MCMC)
    Copyright (C) 2016  Jingjing Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"
#include <gsl/gsl_rng.h>

#include "bvsrm.h"

using namespace std;
using namespace Eigen;

extern "C" void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA, float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC);
extern "C" void spotrf_(char *UPLO, int *N, float *A, int *LDA, int *INFO);
extern "C" void spotrs_(char *UPLO, int *N, int *NRHS, float *A, int *LDA, float *B, int *LDB, int *INFO);
extern "C" void ssyev_(char* JOBZ, char* UPLO, int *N, float *A, int *LDA, float *W, float *WORK, int *LWORK, int *INFO);
extern "C" void ssyevr_(char* JOBZ, char *RANGE, char* UPLO, int *N, float *A, int *LDA, float *VL, float *VU, int *IL, int *IU, float *ABSTOL, int *M, float *W, float *Z, int *LDZ, int *ISUPPZ, float *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);
extern "C" void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
extern "C" void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
extern "C" void dsyev_(char* JOBZ, char* UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);

extern "C" void dsysv_( char* uplo, int* n, int* nrhs, double* a, int* lda,
                      int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info );

extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                      double* b, int* ldb, int* info );

extern "C" void dsyevr_(char* JOBZ, char *RANGE, char* UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);


//cholesky decomposition, A is distroyed
void lapack_float_cholesky_decomp (gsl_matrix_float *A)
{
	int N=A->size1, LDA=A->size1, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_decomp."<<endl; return;}
	
	spotrf_(&UPLO, &N, A->data, &LDA, &INFO);
	if (INFO!=0) {cout<<"Cholesky decomposition unsuccessful in lapack_cholesky_decomp."<<endl; return;}	
	
	return;
}

//cholesky decomposition, A is distroyed
void lapack_cholesky_decomp (gsl_matrix *A)
{
	int N=A->size1, LDA=A->size1, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_decomp."<<endl; return;}
	
	dpotrf_(&UPLO, &N, A->data, &LDA, &INFO);
	if (INFO!=0) {cout<<"Cholesky decomposition unsuccessful in lapack_cholesky_decomp."<<endl; return;}	
	
	return;
}


//cholesky solve, A is decomposed,
void lapack_cholesky_solve (gsl_matrix *A, const gsl_vector *b, gsl_vector *x)
{
    int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
    char UPLO='L';
    
    if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_solve."<<endl; return;}
    
    gsl_vector_memcpy (x, b);
    dpotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
    if (INFO!=0) {cout<<"Cholesky solve unsuccessful in lapack_cholesky_solve."<<endl; return;}
    
    return;
}


//cholesky solve, A is decomposed, 
void lapack_float_cholesky_solve (gsl_matrix_float *A, const gsl_vector_float *b, gsl_vector_float *x)
{
	int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_solve."<<endl; return;}
	
	gsl_vector_float_memcpy (x, b);
	spotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
	if (INFO!=0) {cout<<"Cholesky solve unsuccessful in lapack_cholesky_solve."<<endl; return;}	
	
	return;
}




void lapack_sgemm (char *TransA, char *TransB, float alpha, const gsl_matrix_float *A, const gsl_matrix_float *B, float beta, gsl_matrix_float *C)
{
	int M, N, K1, K2, LDA=A->size1, LDB=B->size1, LDC=C->size2;
	
	if (*TransA=='N' || *TransA=='n') {M=A->size1; K1=A->size2;}
	else if (*TransA=='T' || *TransA=='t') {M=A->size2; K1=A->size1;}
	else {cout<<"need 'N' or 'T' in lapack_sgemm"<<endl; return;}
	
	if (*TransB=='N' || *TransB=='n') {N=B->size2; K2=B->size1;}
	else if (*TransB=='T' || *TransB=='t')  {N=B->size1; K2=B->size2;}
	else {cout<<"need 'N' or 'T' in lapack_sgemm"<<endl;  return;}
	
	if (K1!=K2) {cout<<"A and B not compatible in lapack_sgemm"<<endl; return;}
	if (C->size1!=(size_t)M || C->size2!=(size_t)N) {cout<<"C not compatible in lapack_sgemm"<<endl; return;}
	
	gsl_matrix_float *A_t=gsl_matrix_float_alloc (A->size2, A->size1);
	gsl_matrix_float_transpose_memcpy (A_t, A);
	gsl_matrix_float *B_t=gsl_matrix_float_alloc (B->size2, B->size1);
	gsl_matrix_float_transpose_memcpy (B_t, B);
	gsl_matrix_float *C_t=gsl_matrix_float_alloc (C->size2, C->size1);
	gsl_matrix_float_transpose_memcpy (C_t, C);
	
	sgemm_(TransA, TransB, &M, &N, &K1, &alpha, A_t->data, &LDA, B_t->data, &LDB, &beta, C_t->data, &LDC);
	gsl_matrix_float_transpose_memcpy (C, C_t);
	
	gsl_matrix_float_free (A_t);
	gsl_matrix_float_free (B_t);
	gsl_matrix_float_free (C_t);
	return;
}



void lapack_dgemm (char *TransA, char *TransB, double alpha, const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
{
	int M, N, K1, K2, LDA=A->size1, LDB=B->size1, LDC=C->size2;
	
	if (*TransA=='N' || *TransA=='n') {M=A->size1; K1=A->size2;}
	else if (*TransA=='T' || *TransA=='t') {M=A->size2; K1=A->size1;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl; return;}
	
	if (*TransB=='N' || *TransB=='n') {N=B->size2; K2=B->size1;}
	else if (*TransB=='T' || *TransB=='t')  {N=B->size1; K2=B->size2;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl;  return;}
	
	if (K1!=K2) {cout<<"A and B not compatible in lapack_dgemm"<<endl; return;}
	if (C->size1!=(size_t)M || C->size2!=(size_t)N) {cout<<"C not compatible in lapack_dgemm"<<endl; return;}
	
	gsl_matrix *A_t=gsl_matrix_alloc (A->size2, A->size1);
	gsl_matrix_transpose_memcpy (A_t, A);
	gsl_matrix *B_t=gsl_matrix_alloc (B->size2, B->size1);
	gsl_matrix_transpose_memcpy (B_t, B);
	gsl_matrix *C_t=gsl_matrix_alloc (C->size2, C->size1);
	gsl_matrix_transpose_memcpy (C_t, C);

	dgemm_(TransA, TransB, &M, &N, &K1, &alpha, A_t->data, &LDA, B_t->data, &LDB, &beta, C_t->data, &LDC);

	gsl_matrix_transpose_memcpy (C, C_t);
	
	gsl_matrix_free (A_t);
	gsl_matrix_free (B_t);
	gsl_matrix_free (C_t);
	return;
}



//eigen value decomposition, matrix A is destroyed, float seems to have problem with large matrices (in mac)
void lapack_float_eigen_symmv (gsl_matrix_float *A, gsl_vector_float *eval, gsl_matrix_float *evec, const size_t flag_largematrix)
{
	if (flag_largematrix==1) {
		int N=A->size1, LDA=A->size1, INFO, LWORK=-1;
		char JOBZ='V', UPLO='L';
				
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		//	float temp[1];
		//	ssyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, temp, &LWORK, &INFO);
		//	if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}
		//	LWORK=(int)temp[0];
		
		LWORK=3*N;
		float *WORK=new float [LWORK];	
		ssyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, WORK, &LWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}
		
		gsl_matrix_float_view A_sub=gsl_matrix_float_submatrix(A, 0, 0, N, N);
		gsl_matrix_float_memcpy (evec, &A_sub.matrix);
		gsl_matrix_float_transpose (evec);
		
		delete [] WORK;
	} else {	
		int N=A->size1, LDA=A->size1, LDZ=A->size1, INFO, LWORK=-1, LIWORK=-1;
		char JOBZ='V', UPLO='L', RANGE='A';
		float ABSTOL=1.0E-7;
		
		//VL, VU, IL, IU are not referenced; M equals N if RANGE='A'
		float VL=0.0, VU=0.0;
		int IL=0, IU=0, M;
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_float_eigen_symmv."<<endl; return;}
		
		int *ISUPPZ=new int [2*N];
				
		float WORK_temp[1];
		int IWORK_temp[1];
		ssyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK_temp, &LWORK, IWORK_temp, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_float_eigen_symmv."<<endl; return;}	
		LWORK=(int)WORK_temp[0]; LIWORK=(int)IWORK_temp[0];	
		 
		//LWORK=26*N;
		//LIWORK=10*N;
		float *WORK=new float [LWORK];
		int *IWORK=new int [LIWORK];
		
		ssyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_float_eigen_symmv."<<endl; return;}
		
		gsl_matrix_float_transpose (evec);
		
		delete [] ISUPPZ;
		delete [] WORK;
		delete [] IWORK;
	}
	
	
	return;
}



//eigen value decomposition, matrix A is destroyed
void lapack_eigen_symmv (gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec, const size_t flag_largematrix)
{
	if (flag_largematrix==1) {
		int N=A->size1, LDA=A->size1, INFO, LWORK=-1;
		char JOBZ='V', UPLO='L';		
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		//	double temp[1];
		//	dsyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, temp, &LWORK, &INFO);
		//	if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}		
		//	LWORK=(int)temp[0];
		
		LWORK=3*N;
		double *WORK=new double [LWORK];	
		dsyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, WORK, &LWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}
		
		gsl_matrix_view A_sub=gsl_matrix_submatrix(A, 0, 0, N, N);
		gsl_matrix_memcpy (evec, &A_sub.matrix);
		gsl_matrix_transpose (evec);
		
		delete [] WORK;
	} else {	
		int N=A->size1, LDA=A->size1, LDZ=A->size1, INFO, LWORK=-1, LIWORK=-1;
		char JOBZ='V', UPLO='L', RANGE='A';
		double ABSTOL=1.0E-7;
		
		//VL, VU, IL, IU are not referenced; M equals N if RANGE='A'
		double VL=0.0, VU=0.0;
		int IL=0, IU=0, M;
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		int *ISUPPZ=new int [2*N];
		
		double WORK_temp[1];
		int IWORK_temp[1];

		dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK_temp, &LWORK, IWORK_temp, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}	
		LWORK=(int)WORK_temp[0]; LIWORK=(int)IWORK_temp[0];	

		//LWORK=26*N;
		//LIWORK=10*N;
		double *WORK=new double [LWORK];
		int *IWORK=new int [LIWORK];
		
		dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}

		gsl_matrix_transpose (evec);
		
		delete [] ISUPPZ;
		delete [] WORK;
		delete [] IWORK;
	}
	
	return;
}

//DO NOT set eigen values to be positive
double EigenDecomp (gsl_matrix *G, gsl_matrix *U, gsl_vector *eval, const size_t flag_largematrix)
{
//#ifdef WITH_LAPACK
//	lapack_eigen_symmv (G, eval, U, flag_largematrix);
//#else
	gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc (G->size1);
	gsl_eigen_symmv (G, eval, U, w);
	gsl_eigen_symmv_free (w);	
//#endif
	/*
	for (size_t i=0; i<eval->size; ++i) {
		if (gsl_vector_get (eval, i)<1e-10) {
//			cout<<gsl_vector_get (eval, i)<<endl;
			gsl_vector_set (eval, i, 0);			
		}
	}
	*/
	//calculate track_G=mean(diag(G))	
	double d=0.0;
	for (size_t i=0; i<eval->size; ++i) {
		d+=gsl_vector_get(eval, i);
	}
	d/=(double)eval->size;
	
	return d;
}


//DO NOT set eigen values to be positive
double EigenDecomp (gsl_matrix_float *G, gsl_matrix_float *U, gsl_vector_float *eval, const size_t flag_largematrix)
{
#ifdef WITH_LAPACK
	lapack_float_eigen_symmv (G, eval, U, flag_largematrix);
#else
	//gsl doesn't provide float precision eigen decomposition; plus, float precision eigen decomposition in lapack may not work on OS 10.4
	//first change to double precision
	gsl_matrix *G_double=gsl_matrix_alloc (G->size1, G->size2);
	gsl_matrix *U_double=gsl_matrix_alloc (U->size1, U->size2);
	gsl_vector *eval_double=gsl_vector_alloc (eval->size);
	for (size_t i=0; i<G->size1; i++) {
		for (size_t j=0; j<G->size2; j++) {
			gsl_matrix_set(G_double, i, j, gsl_matrix_float_get(G, i, j));
		}
	}	
	gsl_eigen_symmv_workspace *w_space=gsl_eigen_symmv_alloc (G->size1);
	gsl_eigen_symmv (G_double, eval_double, U_double, w_space);
	gsl_eigen_symmv_free (w_space);	
	
	//change back to float precision
	for (size_t i=0; i<G->size1; i++) {
		for (size_t j=0; j<G->size2; j++) {
			gsl_matrix_float_set(K, i, j, gsl_matrix_get(G_double, i, j));
		}
	}
	for (size_t i=0; i<U->size1; i++) {
		for (size_t j=0; j<U->size2; j++) {
			gsl_matrix_float_set(U, i, j, gsl_matrix_get(U_double, i, j));
		}
	}
	for (size_t i=0; i<eval->size; i++) {
		gsl_vector_float_set(eval, i, gsl_vector_get(eval_double, i));
	}	
	
	//delete double precision matrices
	gsl_matrix_free (G_double);
	gsl_matrix_free (U_double);
	gsl_vector_free (eval_double);
#endif
	/*
	for (size_t i=0; i<eval->size; ++i) {
		if (gsl_vector_float_get (eval, i)<1e-10) {
			gsl_vector_float_set (eval, i, 0);
		}
	}
	*/
	//calculate track_G=mean(diag(G))	
	double d=0.0;
	for (size_t i=0; i<eval->size; ++i) {
		d+=gsl_vector_float_get(eval, i);
	}
	d/=(double)eval->size;
	
	return d;
}

//JY add topdm function
void topdm(gsl_matrix *Omega){
    
    double evali;    
    size_t size_o = Omega->size1;
    if(size_o != Omega->size2) {cout << "Omega size error"; return;}
    
    else{
    gsl_matrix *Omega_temp = gsl_matrix_alloc(size_o, size_o);
    gsl_matrix_memcpy(Omega_temp, Omega);
    gsl_matrix *D = gsl_matrix_alloc(size_o, size_o);
    gsl_matrix *U = gsl_matrix_alloc(size_o, size_o);
    gsl_matrix_set_all(D, 0.0);
    gsl_matrix_set_all(U, 0.0);
        
    gsl_vector *eval_temp = gsl_vector_alloc(size_o);
    gsl_matrix *evec_temp = gsl_matrix_alloc(size_o, size_o);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(size_o);
    gsl_eigen_symmv(Omega_temp, eval_temp, evec_temp, w);
        
        double eval_sum = 0.0;
        size_t posin = 0;
        for (size_t i=0; i<size_o; i++) {
            evali = gsl_vector_get(eval_temp, i);
            if (evali > 0) {
                eval_sum+=gsl_vector_get(eval_temp, i);
                posin++;
            }
            else continue;
        }
        eval_sum /= (double)posin;
        
    size_t neig = 0;
        
        double EPS = eval_sum * 0.0000001;
        for(size_t i = 0; i < size_o; i++){
            evali = gsl_vector_get(eval_temp, i);
            if (evali < EPS)
                {
                    gsl_matrix_set(D, i, i, EPS); neig++;
                   // cout << "eigen value = EPS = " << EPS << endl;
                }
            else {gsl_matrix_set(D, i, i, evali); }
            
    }
    //cout << "first of D = "<<gsl_matrix_get(D, 0, 0)<< "last of D = " << gsl_matrix_get(D, (size_o-1), (size_o-1)) << endl;
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec_temp, D, 0.0, U);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, evec_temp, 0.0, Omega);
    
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(Omega_temp);
    gsl_matrix_free(D);
    gsl_matrix_free(U);
    gsl_matrix_free(evec_temp);
    gsl_vector_free(eval_temp);
    
  //if(neig > 0 ) cout << "topdm success with number of negative eigen = " << neig << endl;
    }
    return;
}
//end

void CholeskyInverse(gsl_matrix *XtX){
    //topdm(XtX);
    gsl_linalg_cholesky_decomp(XtX);
    gsl_linalg_cholesky_invert(XtX);
    return;
}



//Eigensolve start
void EigenInverse(gsl_matrix *XtX){
    size_t s_size = XtX -> size1;
    MatrixXd XtX_eig(s_size, s_size);
    MatrixXd XtX_eiginv(s_size, s_size);
    MatrixXd Identity_Y(s_size, s_size);
    Identity_Y.setIdentity();
 
    for(size_t i=0; i < s_size; ++i){
        for (size_t j=i; j<s_size; ++j) {
            XtX_eig(i, j) = gsl_matrix_get(XtX, i, j);
            if (j != i)
                XtX_eig(j, i) = gsl_matrix_get(XtX, j, i);
        }
    }
    
    XtX_eiginv = XtX_eig.llt().solve(Identity_Y);
    //Standard Cholesky decomposition (LL^T) of a matrix
    
    for(size_t i=0; i < s_size; ++i){
        for (size_t j=i; j<s_size; ++j) {
            gsl_matrix_set(XtX, i, j, XtX_eiginv(i, j));
            if (j != i)
                gsl_matrix_set(XtX, j, i, XtX_eiginv(j, i));
        }
    }
    return;
}

void EigenSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta)
{
    size_t s_size = Xty -> size;
    
    MatrixXd XtX_eig(s_size, s_size);
    VectorXd Xty_eig(s_size);
    VectorXd beta_eig(s_size);
    
    for(size_t i=0; i < s_size; ++i){
        Xty_eig(i) = gsl_vector_get(Xty, i);
        XtX_eig(i, i) = gsl_matrix_get(XtX, i, i);
        for (size_t j=(i+1); j<s_size; ++j) {
            XtX_eig(i, j) = gsl_matrix_get(XtX, i, j);
            XtX_eig(j, i) = gsl_matrix_get(XtX, j, i);
        }
    }
    
    beta_eig = XtX_eig.jacobiSvd(ComputeThinU | ComputeThinV).solve(Xty_eig);
    //beta_eig = XtX_eig.fullPivLu().solve(Xty_eig);
    //beta_eig = XtX_eig.fullPivHouseholderQr().solve(Xty_eig);

    
	for(size_t i=0; i < s_size; ++i){
        gsl_vector_set(beta, i, beta_eig(i));
    }
    
	return;
}


void EigenSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta, const double &lambda)
{
    size_t s_size = Xty -> size;
    
     MatrixXd XtX_eig(s_size, s_size);
     VectorXd Xty_eig(s_size);
     VectorXd beta_eig(s_size);
    
    
   // gsl_rng * r;
   // const gsl_rng_type * T;
    //gsl_rng_env_setup();
    //T = gsl_rng_default;
    //r = gsl_rng_alloc (T);
    //cout << "lambda" << lambda << endl;
    //gsl_rng_set(r, (unsigned)lambda);
    
    
     for(size_t i=0; i < s_size; ++i){
         Xty_eig(i) = gsl_vector_get(Xty, i);
         //cout << "add " << lambda * gsl_rng_uniform(r) << endl;
         
         XtX_eig(i, i) = gsl_matrix_get(XtX, i, i) + lambda;
         //XtX_eig(i, i) = gsl_matrix_get(XtX, i, i) + lambda * gsl_rng_uniform(r);
         
         for (size_t j=(i+1); j<s_size; ++j) {
                 XtX_eig(i, j) = gsl_matrix_get(XtX, i, j);
                 XtX_eig(j, i) = gsl_matrix_get(XtX, j, i);
             }
     }
    //Eigen::JacobiSVD.setThreshold(0.00001);
    beta_eig = XtX_eig.jacobiSvd(ComputeThinU | ComputeThinV).solve(Xty_eig);
   // beta_eig = XtX_eig.fullPivLu().solve(Xty_eig);
    //beta_eig = XtX_eig.fullPivHouseholderQr().solve(Xty_eig);
    
	for(size_t i=0; i < s_size; ++i){
        gsl_vector_set(beta, i, beta_eig(i));
    }
    
    //gsl_rng_free(r);
}

void GSLSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta_hat, const double lambda)
{
    int s_size = Xty -> size;
    gsl_matrix *XtXtemp = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(XtXtemp, XtX);
    gsl_vector_view XtXtemp_diag = gsl_matrix_diagonal(XtXtemp);
    gsl_vector_add_constant(&XtXtemp_diag.vector, lambda);
    
    gsl_permutation * p = gsl_permutation_alloc(s_size);
    //int s = s_size;
    gsl_linalg_LU_decomp(XtXtemp, p, &s_size);
    gsl_linalg_LU_solve(XtXtemp, p, Xty, beta_hat);
    gsl_permutation_free(p);
    gsl_matrix_free(XtXtemp);
    
}

//JY write eigen solve

double CalcLogdet(gsl_matrix *Omega)
{
    double logdet_O=0.0;
    gsl_matrix *OmegaTemp = gsl_matrix_alloc(Omega->size1, Omega->size2);
    gsl_matrix_memcpy(OmegaTemp, Omega);
    
    topdm(OmegaTemp);
    gsl_linalg_cholesky_decomp(OmegaTemp);
    //cout << "cholesky decompose Omega success" << endl;
    
    for (size_t i=0; i< OmegaTemp->size1; ++i) {
        logdet_O+=log(gsl_matrix_get (OmegaTemp, i, i));
    }
    logdet_O*=2.0;
    
    gsl_matrix_free(OmegaTemp);
    return logdet_O;
}


double LapackCholSolve(gsl_matrix *Omega, const gsl_vector *Xty, gsl_vector *OiXty)
{
    double logdet_O=0.0;
    topdm(Omega);
    
    //cholesky decomposition, A is distroyed
    lapack_cholesky_decomp(Omega);
    for (size_t i=0; i<Omega->size1; ++i) {
        logdet_O+=log(gsl_matrix_get (Omega, i, i));
    }
    logdet_O*=2.0;
    
    return logdet_O;
    
    lapack_cholesky_solve(Omega, Xty, OiXty);
    
    return logdet_O;
}

double LapackLogDet(const gsl_matrix *Omega)
{
    double logdet_O=0.0;
    gsl_matrix *OmegaTemp = gsl_matrix_alloc(Omega->size1, Omega->size2);
    gsl_matrix_memcpy(OmegaTemp, Omega);
    
    topdm(OmegaTemp); //cholesky decomposition, A is distroyed
    lapack_cholesky_decomp(OmegaTemp);
    
    for (size_t i=0; i< OmegaTemp->size1; ++i) {
        logdet_O+=log(gsl_matrix_get (OmegaTemp, i, i));
    }
    logdet_O*=2.0;
    
    gsl_matrix_free(OmegaTemp);
    return logdet_O;

}


double CholeskySolve(gsl_matrix *Omega, const gsl_vector *Xty, gsl_vector *OiXty)
{
    /*size_t s_size = Xty->size;
    double lambda = 0.0;
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(Omega, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.001;
    //cout << "lambda = " << lambda << endl;
    EigenSolve(Omega, Xty, OiXty, lambda); */
    //cout << "lambda = " << lambda << endl;
    //EigenSolve(Omega, Xty, OiXty);
    
	double logdet_O=0.0;
    topdm(Omega);
    gsl_linalg_cholesky_decomp(Omega);
    //cout << "cholesky decompose Omega success" << endl;
	
	for (size_t i=0; i<Omega->size1; ++i) {
		logdet_O+=log(gsl_matrix_get (Omega, i, i));
	}
	logdet_O*=2.0;
    
    //gsl_linalg_cholesky_solve(Omega, Xty, OiXty);
    
	return logdet_O;
}

void Ginv(gsl_matrix *XtX_gtemp){
    
    size_t n_row = XtX_gtemp->size1;
    gsl_vector *svd_work = gsl_vector_alloc(n_row);
    gsl_vector *S = gsl_vector_alloc(n_row);
    gsl_matrix *U = gsl_matrix_alloc(n_row, n_row);
    gsl_matrix *V = gsl_matrix_alloc(n_row, n_row);
    gsl_linalg_SV_decomp(XtX_gtemp, V, S, svd_work);
    gsl_matrix_memcpy(U, XtX_gtemp);
    //cout << "SVD decomposition success"<< endl;
    //cout << "last value of U = " << gsl_matrix_get(XtX_gtemp, s_size-1, s_size-1) << endl;
    //cout << "first Singular value = " << gsl_vector_get(S, 0) << endl;
    //cout << "last Singular value = " << gsl_vector_get(S, (s_size-1)) << endl;
    //cout << "last value of V = " << gsl_matrix_get(V, s_size-1, s_size-1) << endl;
    
    gsl_matrix * SI = gsl_matrix_calloc (n_row, n_row);
    gsl_matrix_set_all(SI, 0.0);
	for (size_t i = 0; i < n_row; i++) {
        double tol = 0.000001 * gsl_vector_get(S, 0) * 0.5 * sqrt(n_row + n_row + 1);
        if (gsl_vector_get (S, i) > tol)
            gsl_matrix_set (SI, i, i, (1.0 / gsl_vector_get (S, i)));
	}
	//THE PSEUDOINVERSE//
	gsl_matrix * SIpUT = gsl_matrix_alloc (n_row, n_row);
    gsl_matrix_set_all(SIpUT, 0.0);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, SI, U, 0.0, SIpUT);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, V, SIpUT, 0.0, XtX_gtemp);
    //cout << "persudo inverse success" << endl;
    
	gsl_matrix_free(SI);
	gsl_matrix_free(SIpUT);
    gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
    gsl_vector_free(svd_work);
    
    return;
}
/*
 void lapack_cholesky_decomp (gsl_matrix *A)
 {
 int N=A->size1, LDA=A->size1, INFO;
 char UPLO='L';
 
 if (N!=(int)A->size2) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_decomp."<<endl; return;}
 
 dpotrf_(&UPLO, &N, A->data, &LDA, &INFO);
 if (INFO!=0) {cout<<"Cholesky decomposition unsuccessful in lapack_cholesky_decomp."<<endl; return;}
 
 return;
 }
 
 
 //cholesky solve, A is decomposed,
 void lapack_cholesky_solve (gsl_matrix *A, const gsl_vector *b, gsl_vector *x)
 {
 int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
 char UPLO='L';
 
 if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_solve."<<endl; return;}
 
 gsl_vector_memcpy (x, b);
 dpotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
 if (INFO!=0) {cout<<"Cholesky solve unsuccessful in lapack_cholesky_solve."<<endl; return;}
 
 return;
 }
 */

/*int LapackSolve(gsl_matrix *A, gsl_vector *b, gsl_vector *x){
    
    int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
    int ipiv[N];
    
    if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in LapackSolve."<<endl; exit(-1);}
    
    gsl_vector_memcpy (x, b);
    dgesv_( &N, &NRHS, A->data, &LDA, ipiv, x->data, &LDB, &INFO);
    if (INFO!=0) {
        //PrintVector(x);
        cout<<"INFO ="<< INFO << ": Lapack solve unsuccessful in LapackSolve."<<endl;
    }
    
    return INFO;
}*/

int LapackSolve(const gsl_matrix *XtX, const gsl_vector *b, gsl_vector *x){
    
    gsl_matrix *A = gsl_matrix_alloc(XtX->size1, XtX->size2);
    gsl_matrix_memcpy(A, XtX);
    
    int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
    int ipiv[N];
    
    if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in LapackSolve."<<endl; exit(-1);}
    
    gsl_vector_memcpy (x, b);
    dgesv_( &N, &NRHS, A->data, &LDA, ipiv, x->data, &LDB, &INFO);
    
    if (INFO!=0) {
     //   PrintVector(x);
        cout<<"INFO ="<< INFO << ": Lapack solve unsuccessful in LapackSolve."<<endl;
    }
    
    gsl_matrix_free(A);
    return INFO;
}


double CholeskySolve(gsl_matrix_float *Omega, gsl_vector_float *Xty, gsl_vector_float *OiXty)
{
	double logdet_O=0.0;
	
//#ifdef WITH_LAPACK
//	lapack_float_cholesky_decomp(Omega);
//	for (size_t i=0; i<Omega->size1; ++i) {
//		logdet_O+=log(gsl_matrix_float_get (Omega, i, i));
//	}
//	logdet_O*=2.0;
//	lapack_float_cholesky_solve(Omega, Xty, OiXty);
//#else
    
	gsl_matrix *Omega_double=gsl_matrix_alloc (Omega->size1, Omega->size2);
	double d;
	for (size_t i=0; i<Omega->size1; ++i) {
		for (size_t j=0; j<Omega->size2; ++j) {
			d=(double)gsl_matrix_float_get (Omega, i, j);
			gsl_matrix_set (Omega_double, i, j, d);
		}
	}
	
	int status = gsl_linalg_cholesky_decomp(Omega_double);
	if(status == GSL_EDOM) {
		cout << "## non-positive definite matrix" << endl;
		//		exit(0); 
	}	
	
	for (size_t i=0; i<Omega->size1; ++i) {
		for (size_t j=0; j<Omega->size2; ++j) {
			d=gsl_matrix_get (Omega_double, i, j);
			if (j==i) {logdet_O+=log(d);}
			gsl_matrix_float_set (Omega, i, j, (float)d);
		}
	}
	logdet_O*=2.0;	
	
	gsl_vector_float_memcpy (OiXty, Xty);
	gsl_blas_strsv(CblasLower, CblasNoTrans, CblasNonUnit, Omega, OiXty); 
	gsl_blas_strsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, OiXty); 	
	//	gsl_linalg_cholesky_solve(XtX, Xty, iXty);
	
	gsl_matrix_free (Omega_double);
//#endif
	
	return logdet_O;
}	


//LU decomposition
void LUDecomp (gsl_matrix *LU, gsl_permutation *p, int *signum)
{
	gsl_linalg_LU_decomp (LU, p, signum);
	return;
}

void LUDecomp (gsl_matrix_float *LU, gsl_permutation *p, int *signum)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	gsl_linalg_LU_decomp (LU_double, p, signum);
	
	//copy float matrix to double
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_float_set (LU, i, j, gsl_matrix_get(LU_double, i, j));
		}
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	return;
}


//LU invert
void LUInvert (const gsl_matrix *LU, const gsl_permutation *p, gsl_matrix *inverse)
{
	gsl_linalg_LU_invert (LU, p, inverse);
	return;
}

void LUInvert (const gsl_matrix_float *LU, const gsl_permutation *p, gsl_matrix_float *inverse)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	gsl_matrix *inverse_double=gsl_matrix_alloc (inverse->size1, inverse->size2);
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	gsl_linalg_LU_invert (LU_double, p, inverse_double);
	
	//copy float matrix to double
	for (size_t i=0; i<inverse->size1; i++) {
		for (size_t j=0; j<inverse->size2; j++) {
			gsl_matrix_float_set (inverse, i, j, gsl_matrix_get(inverse_double, i, j));
		}
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	gsl_matrix_free (inverse_double);
	return;
}

//LU lndet
double LULndet (gsl_matrix *LU)
{
	double d;
	d=gsl_linalg_LU_lndet (LU);
	return d;
}

double LULndet (gsl_matrix_float *LU)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	double d;
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	d=gsl_linalg_LU_lndet (LU_double);
	
	//copy float matrix to double
	/*
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_float_set (LU, i, j, gsl_matrix_get(LU_double, i, j));
		}
	}
	*/
	//free matrix
	gsl_matrix_free (LU_double);
	return d;
}


//LU solve
void LUSolve (const gsl_matrix *LU, const gsl_permutation *p, const gsl_vector *b, gsl_vector *x)
{
	gsl_linalg_LU_solve (LU, p, b, x);
	return;
}

void LUSolve (const gsl_matrix_float *LU, const gsl_permutation *p, const gsl_vector_float *b, gsl_vector_float *x)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	gsl_vector *b_double=gsl_vector_alloc (b->size);
	gsl_vector *x_double=gsl_vector_alloc (x->size);	
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	for (size_t i=0; i<b->size; i++) {
		gsl_vector_set (b_double, i, gsl_vector_float_get(b, i));
	}
	
	for (size_t i=0; i<x->size; i++) {
		gsl_vector_set (x_double, i, gsl_vector_float_get(x, i));
	}
	
	//LU decomposition
	gsl_linalg_LU_solve (LU_double, p, b_double, x_double);
	
	//copy float matrix to double
	for (size_t i=0; i<x->size; i++) {
		gsl_vector_float_set (x, i, gsl_vector_get(x_double, i));
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	gsl_vector_free (b_double);
	gsl_vector_free (x_double);
	return;
}


