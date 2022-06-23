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

#ifndef __BVSRM_H__                
#define __BVSRM_H__

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector_double.h>

#include "param.h"

using namespace std;




class BVSRM {

public:
    //multiple function related parameters
    size_t n_type;
    vector<size_t> mFunc; // # of variants of each variant type
    double e, e_shape, e_rate; //hyper parameter in the prior gamma distribution
    double vscale;
    map<string, int> mapFunc2Code;
    int iniType;
    bool FIXHYP;
    bool saveSNP;
    bool saveLD;
    string iniSNPfile;
    string hypfile;
    vector< pair<size_t, size_t> > SNPorder_vec; //<pos, rank>
    vector< pair<size_t, size_t> > SNPrank_vec; //<pos, order>
    double GV, rv, tau, logrv;
    vector<double> SNPsd, XtX_diagvec;
    vector<double> SNPmean;
    
    
    double h, h_global;
    vector <double> rho_vec;
    vector <double> Gvec;
    vector <double> theta; // global hyper parameter
    vector <double> log_theta;
    vector <double> log_qtheta;
    vector <double> theta_total;
    vector <double> subvar, inv_subvar, log_subvar, sumbeta2; // global hyper parameter

    
	// IO related parameters
    size_t UnCompBufferSize;
    vector <size_t> CompBuffSizeVec;
    bool Compress_Flag;
	int a_mode;
	size_t d_pace;
    vector<pair<int, double> > UcharTable;
    vector<string> InputSampleID; //size = ni_total
    vector<string> VcfSampleID_test; // sample id for ni_test in order of X[i][j]
	
	string file_bfile;
	string file_geno;
    string file_vcf;
	string file_out;
	
	// LMM related parameters
	double l_min;
	double l_max;
	size_t n_region;
	double pve_null;
	double pheno_mean;
    double pheno_var;
	
	// BSLMM MCMC related parameters
    
    //JY added win, Wvar, ns_neib;
    size_t win, ns_neib;
    size_t nadd_accept, ndel_accept, nswitch_accept, nother_accept;
    size_t nadd, ndel, nswitch, nother;
    int Switch_Flag;
    
    
	double h_min, h_max, h_scale;			//priors for h
	double rho_min, rho_max, rho_scale;		//priors for rho
	double logp_min, logp_max, logp_scale;		//priors for log(pi)
	size_t s_min, s_max;			//minimum and maximum number of gammas
	size_t w_step;					//number of warm up/burn in iterations
	size_t s_step;					//number of sampling iterations
	size_t r_pace;					//record pace
	size_t w_pace;					//write pace
	size_t n_accept;				//number of acceptance
	size_t n_mh;					//number of MH steps within each iteration
    size_t region_pip;              //number of MCMC with SNPs>0
	double geo_mean;				//mean of the geometric distribution
	long int randseed;
	double trace_G;	
	
	HYPBSLMM cHyp_initial;

	// Summary statistics
	size_t ni_total, ns_total;	//number of total individuals and snps
	size_t ni_test, ns_test;	//number of individuals and snps used for analysis
	size_t n_cvt;				//number of covariates
	double time_UtZ;
	double time_Omega;		//time spent on optimization iterations
	double time_Proposal;        //time spent on constructing the proposal distribution for gamma (i.e. lmm or lm analysis)
	vector<bool> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis

	vector<bool> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	
	vector<SNPINFO> snpInfo;		//record SNP information
	
	// Not included in PARAM
	gsl_rng *gsl_r;
	gsl_ran_discrete_t *gsl_t; //JY added dynamic gsl_s
    
	map<size_t, size_t> mapRank2pos;
	map<size_t, size_t> mapOrder2pos; // JY: map order index to snp position
    map<size_t, size_t> mapPos2Order; // JY: map position to snp order
    map<size_t, size_t> mapPos2Rank; // JY: map position to snp rank
    map<size_t, size_t> mapRank2Order; // JY: map rank index to snp order
    map<size_t, size_t> mapOrder2Rank; // JY: map order index to snp rank

	
	// ************* Main Functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);	

    
    //************* MCMC related functions
    void MCMC (uchar **X, const gsl_vector *y, bool original_method);

    void WriteMCMC(const vector<string> &snps_mcmc);

    void WriteHyptemp(gsl_vector *LnPost, vector<double> &em_gamma);

    void WriteParam(vector<pair<double, double> > &beta_g, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr, const vector<double> &Z_scores, const vector<double> &SE_beta, const vector<double> pval_lrt);

    // Initialize the model
    void CalcPgamma (double *p_gamma, size_t p_gamma_top);

    void setHyp(double theta_temp, double subvar_temp);

    void CreateSnpPosVec(vector<SNPPOS> &snp_pos);

    void InitialMCMC ( uchar **UtX, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos);

    void set_mgamma(class HYPBSLMM &cHyp, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos);

    void WriteIniSNP (const vector<size_t> &rank, const vector<SNPPOS> &snp_pos);
    
    void WriteIniRank (const vector<string> &iniRank);

    // MH-related functions
    double ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, uchar **X, const gsl_vector *z, const gsl_matrix *Xgamma_old, const gsl_matrix *XtX_old, const gsl_vector *Xtz_old, const double &ztz, int &flag_gamma, gsl_matrix *Xgamma_new, gsl_matrix *XtX_new, gsl_vector *Xtz_new);

    double CalcLikegamma(const class HYPBSLMM &cHyp);
            
    void SetXgamma (gsl_matrix *Xgamma, uchar **X, vector<size_t> &rank);

    void SetXgamma ( uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, const vector<size_t> &rank_new, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new);

    void SetXgammaDel(const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const vector<size_t> &rank_old, size_t col_id, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new);

    void SetXgammaAdd (uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, size_t ranki, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new);

    void CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty);
    
    double CalcLikelihood (const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag);   

    double CalcPosterior (const double yty, class HYPBSLMM &cHyp);

    double CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag, double &loglike);

    void CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX_gamma, const gsl_vector *Xtz_gamma, gsl_vector *z_res, const size_t &s_size, const double &ztz);
    
    // Others
    gsl_ran_discrete_t * MakeProposal(const size_t &o, double *p_BF, uchar **X, const gsl_vector *z_res, const map<size_t, int> &mapRank2in);

    double CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec, size_t posj);

    void getSubVec(gsl_vector *sigma_subvec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos);
    
    bool ColinearTest(uchar ** X, const gsl_matrix * Xtemp, const gsl_matrix * XtX_temp, size_t r_add, size_t s_size);
    
    void WriteGenotypeFile(uchar **X, const vector<SNPPOS> &snp_pos);

    void WriteFGWAS_InputFile(const vector<SNPPOS> &snp_pos, const vector<double> &Z_scores, const vector<double> &SE_beta);

    void CalcLD(uchar **X, gsl_matrix * LD);

    double CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2);

    void SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z);

    void CalcCC_PVEnZ (gsl_vector *z_hat, class HYPBSLMM &cHyp);

    void CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp);

};



#endif


