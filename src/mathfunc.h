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

#ifndef __MATHFUNC_H__                
#define __MATHFUNC_H__

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef short int int16;

using namespace std;

double VectorVar (const gsl_vector *v);
void CenterMatrix (gsl_matrix *G);
void CenterMatrix (gsl_matrix *G, gsl_vector *w);
void ScaleMatrix (gsl_matrix *G);
double CenterVector (gsl_vector *y);
void CalcUtX (const gsl_matrix *U, gsl_matrix *UtX);
void CalcUtX (const gsl_matrix *U, const gsl_matrix *X, gsl_matrix *UtX);
void CalcUtX (const gsl_matrix *U, const gsl_vector *x, gsl_vector *Utx);
double CalcHWE (const size_t n_hom1, const size_t n_hom2, const size_t n_ab);
void Kronecker(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H);
void KroneckerSym(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H);

double VectorSum (const vector<double> &v);
double VectorMean (const vector<double> &v);

void PrintVector(const gsl_vector * x);
void PrintVector(const vector <double> &x);
void PrintVector(const vector <size_t> &x);
void PrintVector(const double *x);
void PrintVector(const vector <double> &x, const size_t s);
void PrintVector(const uchar *x, const size_t length);
void PrintVector(const gsl_vector * x, const size_t s);
void PrintMatrix(const gsl_matrix * X, const size_t nrow, const size_t ncol);

void expVector(vector<double> &expvec, vector<double> &logvec);
void CalcXVbeta(gsl_matrix *X, const gsl_vector * sigma_vec);

bool comp_sigma (pair<double, size_t> a, pair<double, size_t> b);
bool comp_lr (pair<size_t, double> a, pair<size_t, double> b);
bool comp_res (pair<size_t, double> a, pair<size_t, double> b);

void NormRes(gsl_vector * z_res);

#endif
