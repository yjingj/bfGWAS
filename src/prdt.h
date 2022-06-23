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

#ifndef __PRDT_H__                
#define __PRDT_H__


#include <vector>
#include <map>
#include <string.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#include "param.h"


using namespace std;
typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;

using namespace std;

class PRDT {
	
public:
	// IO related parameters
	size_t a_mode;
	size_t d_pace;
	
	string file_bfile;
	string file_geno;
	string file_out;
	
	vector<vector<bool> > indicator_pheno;
	vector<bool> indicator_cvt;
	vector<bool> indicator_idv;
	vector<SNPINFO> snpInfo;
	map<string, double> mapRS2est;
	
	size_t n_ph;
	size_t np_obs, np_miss;
	size_t ns_total;
	size_t ns_test;
	
	double time_eigen;
	
	// Main functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	void WriteFiles (gsl_vector *y_prdt);
	void WriteFiles (gsl_matrix *Y_full);
	void AddBV (gsl_matrix *G, const gsl_vector *u_hat, gsl_vector *y_prdt);
	void AnalyzeBimbam (gsl_vector *y_prdt);
	void AnalyzePlink (gsl_vector *y_prdt);
	void MvnormPrdt (const gsl_matrix *Y_hat, const gsl_matrix *H, gsl_matrix *Y_full);
};


#endif







