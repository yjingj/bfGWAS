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

#ifndef __PARAM_H__                
#define __PARAM_H__

#include <vector>
#include <map>
#include <set>
#include <iostream>     // std::cout, std::endl
#include <iomanip>
#include <limits>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#include "VcfFileReader.h"
#include "StringBasics.h"
#include "StringHash.h"
#include "MemoryAllocators.h"


using namespace std;
typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;


class SNPINFO {
public:
	string chr;
	string rs_number;
	double cM;
	long int base_position;
	string a_minor;
	string a_major;
	int n_miss;
	double missingness;
	double maf;
    
    vector<bool> indicator_func;
    vector<double> weight;
    double weight_i;
    
    void printMarker();
};

struct genMarker
{
    string rs;
    string chr;
    long int bp;
    string Ref;
    string Alt;
    
    void iniRecord(VcfRecord& record);
    void printMarker();
    
};

//JY
struct SNPPOS{
    size_t pos;
    string rs;
    string chr;
    long int bp;
    string a_minor;
	string a_major;
	double maf;
    vector<bool> indicator_func;
    //vector<double> weight;
    //double weight_i;
    
    void printMarker();
};

//JY
void printSNPInfo(vector<SNPPOS> &snp_pos, int numMarker);
void CalcWeight(const vector<bool> &indicator_func, vector<double> &weight, const double weight_i);

//results for lmm
class SUMSTAT {
public:
	double beta;				//REML estimator for beta
	double se;				//SE for beta  
	double lambda_remle;		//REML estimator for lambda
	double lambda_mle;		//MLE estimator for lambda
	double p_wald;			//p value from a Wald test
	double p_lrt;				//p value from a likelihood ratio test
	double p_score;			//p value from a score test

	double score_u; // for meta-analysis
	double score_v;
};

//results for mvlmm
class MPHSUMSTAT {
public:
	vector<double> v_beta;	//REML estimator for beta
	double p_wald;			//p value from a Wald test
	double p_lrt;				//p value from a likelihood ratio test
	double p_score;			//p value from a score test
	vector<double> v_Vg;	//estimator for Vg, right half
	vector<double> v_Ve;	//estimator for Ve, right half
	vector<double> v_Vbeta;	//estimator for Vbeta, right half
};


//hyper-parameters for bslmm
class HYPBSLMM {
public:
    double h;
    double rho;
    double logp;
    
    vector<double> theta;
    vector<double> log_theta; // log(theta) for each function type
    vector<double> subvar; // variance for each function type
    vector<double> rho_vec;
    vector<size_t> m_gamma; // # of selected SNP of each function type
    
    double rv;
    double sigma_b2; //
	double pve;
	double pge;
	size_t n_gamma;
    
};


class PARAM {
public:
    //multiple functionrelated parameters
    size_t n_type;
    vector<size_t> mFunc; // # of variants of each variant type
    double e; //hyper parameter in the prior gamma distribution
    double vscale;
    map<string, int> mapFunc2Code;
    int iniType;
    bool FIXHYP;
    bool saveSNP;
    bool saveGeno;
    bool saveLD;

    string iniSNPfile;
    string hypfile;
    double rv;
    double pheno_mean, pheno_var;
    
	// IO related parameters
    size_t UnCompBufferSize;
    vector <size_t> CompBuffSizeVec;
    bool Compress_Flag;
    StringIntHash sampleID2vcfInd;
    map<string, size_t>  PhenoID2Ind;
    vector<size_t> SampleVcfPos;
    vector<string> VcfSampleID; // size=total sample #
    vector<string> VcfSampleID_test; // sample id for ni_test in order of X[i][j]
    vector<string> InputSampleID; //size = ni_total
    vector<pair<int, double> > UcharTable;
    vector<double> SNPmean;
    
	bool mode_silence;
	int a_mode;				//analysis mode, 1/2/3/4 for Frequentist tests
	int k_mode;				//kinship read mode: 1: n by n matrix, 2: id/id/k_value; 
	size_t d_pace;		
	
    string file_func_code; //coded all unique variant function types
    string file_sample;
    string file_vcf;
    string GTfield;
    
	string file_bfile;
	string file_geno;
	string file_pheno;
	string file_anno;		
	string file_kin;
	string file_out;
	string file_snps;		//file containing analyzed snps or genes

	
	// QC related parameters	
	double miss_level;
	double maf_level;	
	double hwe_level;
	
	// BVSRM MCMC related parameters
    size_t win;
    size_t nadd_accept, ndel_accept, nswitch_accept, nother_accept;
    size_t nadd, ndel, nswitch, nother;
    
	double h_min, h_max, h_scale;			//priors for h
	double rho_min, rho_max, rho_scale;		//priors for rho
	double logp_min, logp_max, logp_scale;		//priors for log(pi)
	size_t s_min, s_max;			//minimum and maximum number of gammas
	size_t w_step;					//number of warm up/burn in iterations
	size_t s_step;					//number of sampling iterations
	size_t n_accept;				//number of acceptance
	size_t n_mh;					//number of MH steps within each iteration
    size_t region_pip;              //number of MCMC with SNPs>0
	long int randseed;
	double trace_G;
	
	HYPBSLMM cHyp_initial;
		
	// Summary statistics
	bool error;
	size_t ni_total, ni_test, ni_cvt;	//number of individuals
	size_t np_obs, np_miss;		//number of observed and missing phenotypes
    
	size_t ns_total, ns_test;	//number of snps
	size_t n_cvt;
    
	size_t ng_total, ng_test;	//number of genes
	size_t ni_control, ni_case;	//number of controls and number of cases
	double time_total;		//record total time
	double time_G;			//time spent on reading files the second time and calculate K
	double time_opt;		//time spent on optimization iterations/or mcmc
	double time_Omega;		//time spent on calculating Omega
	double time_hyp;		//time spent on sampling hyper-parameters, in PMM
	double time_Proposal;  //time spend on constructing the proposal distribution (i.e. the initial lmm or lm analysis)

	// Data
    
	vector<double>  pheno;	//a vector record all phenotypes, NA replaced with -9

	vector<bool>  indicator_pheno;			//a matrix record when a phenotype is missing for an individual; 0 missing, 1 available

	vector<bool> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis

	vector<bool> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	

    map<int, string> mapCode2Func; // map unique code to a unique function type
	map<string, int> mapID2num;		//map small ID number to number, from 0 to n-1
	map<string, string> mapRS2chr;		//map rs# to chromosome location
	map<string, long int> mapRS2bp;		//map rs# to base position
	map<string, double> mapRS2cM;		//map rs# to cM
	map<string, double> mapRS2est;			//map rs# to parameters
	
	vector<SNPINFO> snpInfo;		//record SNP information
	set<string> setSnps;			//a set of snps for analysis
	
	//constructor
	PARAM();
	
	//functions
	void ReadFiles ();		
	void CheckParam (); 
	void CheckData ();	
	void PrintSummary ();
	void ReadGenotypes (uchar **X, gsl_matrix *K, const bool calc_K);
	void WriteGenotypes(uchar **X);
    
	//void CheckCvt ();
	void ProcessPheno();
    void ReorderPheno(gsl_vector *y);
    
	void CopyPheno (gsl_vector *y);
	void CalcKin (gsl_matrix *matrix_kin);

    
};


#endif

