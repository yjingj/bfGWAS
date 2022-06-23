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

#ifndef __LM_H__                
#define __LM_H__

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
 
#include "param.h"
#include "io.h"
#include "ReadVCF.h"


using namespace std;
typedef short int int16;
typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;

class LM {
	
public:
	// IO related parameters
	int a_mode;				//analysis mode, 50+1/2/3/4 for Frequentist tests
	size_t d_pace;		//display pace
	
	string file_bfile;
	string file_geno;
	string file_out;
	string file_vcf;
		
	// Summary statistics
	size_t ni_total, ni_test;	//number of individuals
	size_t ns_total, ns_test;	//number of snps
	size_t n_cvt;
	double time_opt;		//time spent
	double pheno_var, pheno_mean;
	
	vector<bool> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
	vector<bool> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	
	vector<SNPINFO> snpInfo;		//record SNP information
	
	// Not included in PARAM
	vector<SUMSTAT> sumStat;		//Output SNPSummary Data
	
	// Main functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	void AnalyzePlink (const gsl_matrix *W, const gsl_vector *y);
	void AnalyzeGeno (const gsl_matrix *W, const gsl_vector *y, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID);
	void AnalyzeVCF (const gsl_matrix *W, const gsl_vector *y, string &GTfield, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID); // NEED to be rewritten
	void WriteFiles ();
};

void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, double &trace_G, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag);

void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, const vector<double> &SNPmean, vector<double> &Gvec, vector<double>&XtX_diagvec, const vector<SNPPOS> &snp_pos, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag, const vector<pair<int, double> > &UcharTable);

void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, const vector<double> &SNPmean, vector<double> &Gvec, vector<double> &XtX_diagvec, vector<double> &Z_scores, vector<double> &SE_beta, vector<double> &pval_lrt, const vector<SNPPOS> &snp_pos, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag, const vector<pair<int, double> > &UcharTable);

#endif
