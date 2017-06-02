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
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include <limits>

#include "lapack.h"
#include "param.h"
#include "bvsrm.h"
#include "lm.h"
#include "mathfunc.h"

using namespace std;

// define to_string function : convert to string
template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

void BVSRM::CopyFromParam (PARAM &cPar) 
{
    VcfSampleID_test = cPar.VcfSampleID_test;

    rv = cPar.rv;
    n_type = cPar.n_type;
    mFunc = cPar.mFunc;
    e = cPar.e;
    vscale = cPar.vscale;
    FIXHYP = cPar.FIXHYP;
    saveSNP = cPar.saveSNP;
    saveLD = cPar.saveLD;
    iniType = cPar.iniType;
    iniSNPfile = cPar.iniSNPfile;
    hypfile = cPar.hypfile;
    SNPmean = cPar.SNPmean;

    
    UnCompBufferSize = cPar.UnCompBufferSize;
    CompBuffSizeVec = cPar.CompBuffSizeVec;
    Compress_Flag = cPar.Compress_Flag;
    win = cPar.win;
    nadd_accept = cPar.nadd_accept;
    ndel_accept= cPar.ndel_accept;
    nswitch_accept= cPar.nswitch_accept;
    nother_accept= cPar.nother_accept;
    nadd= cPar.nadd;
    ndel= cPar.ndel;
    nswitch= cPar.nswitch;
    nother= cPar.nother;
    
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
    file_vcf = cPar.file_vcf;
	file_out=cPar.file_out;
	
	l_min=cPar.h_min;	
	l_max=cPar.h_max;  
	pheno_mean=cPar.pheno_mean;
	
	time_UtZ=0.0;
	time_Omega=0.0;
	n_accept=0;
    region_pip = 0;
    Switch_Flag = 0;
	
	h_min=cPar.h_min;	
	h_max=cPar.h_max;  
	h_scale=cPar.h_scale;
	rho_min=cPar.rho_min;	
	rho_max=cPar.rho_max;  
	rho_scale=cPar.rho_scale;
	logp_min=cPar.logp_min;	
	logp_max=cPar.logp_max;  
	logp_scale=cPar.logp_scale;
	
	s_min=cPar.s_min;
	s_max=cPar.s_max;
	w_step=cPar.w_step;
	s_step=cPar.s_step;
	n_mh=cPar.n_mh;
	randseed=cPar.randseed;
	trace_G=cPar.trace_G;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	
	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;
	
	return;
}


void BVSRM::CopyToParam (PARAM &cPar) 
{
	cPar.time_Omega=time_Omega;
	cPar.time_Proposal=time_Proposal;
	cPar.cHyp_initial=cHyp_initial;
	cPar.n_accept=n_accept;
	cPar.pheno_mean=pheno_mean;
    cPar.pheno_var = pheno_var;
	cPar.randseed=randseed;
    
    cPar.nadd_accept=nadd_accept;
    cPar.ndel_accept=ndel_accept;
    cPar.nswitch_accept=nswitch_accept;
    cPar.nother_accept=nother_accept;
    cPar.nadd=nadd;
    cPar.ndel=ndel;
    cPar.nswitch=nswitch;
    cPar.nother=nother;
    cPar.region_pip = region_pip;
	
	return;
}

bool comp_vec (size_t a, size_t b)
{
	return (a < b);
}

bool comp_pi (pair<string, double> a, pair<string, double> b)
{
    return (a.second > b.second);
}

//JY add to write initial significant SNP id out
void BVSRM::WriteIniRank (const vector<string> &iniRank)
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".inirank.txt";
    
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	for (size_t i=0; i<iniRank.size(); ++i) {

			outfile<<iniRank[i]<<endl;

	}
	
	outfile.clear();
	outfile.close();
	return;
}

void BVSRM::WriteIniSNP (const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".iniSNP";
    //cout << "write iniSNP at " << file_str << endl;
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for (size_t i=0; i<rank.size(); ++i) {
        outfile<< snp_pos[SNPrank_vec[rank[i]].second].rs <<endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
}

void BVSRM::WriteMCMC(const vector<string> &snps_mcmc){
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".mcmc";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for (size_t i=0; i<snps_mcmc.size(); ++i) {
        outfile<< snps_mcmc[i] <<endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
}


void BVSRM::WriteParam(vector<pair<double, double> > &beta_g, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr, const vector<double> &Z_scores, const vector<double> &SE_beta, const vector<double> pval_lrt)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".paramtemp";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    //outfile<<"markerID"<<"\t"<<"chr"<<"\t" <<"bp"<<"\t" <<"REF"<<"\t" <<"ALT"<<"\t" << "maf" << "\t" << "Func_code"<< "\t" <<"beta"<<"\t"<<"gamma" << "\t" <<"Zscore" << "\t" << "SE_beta" << "\t" << "LRT" << "\t" << "pval_lrt"  << "\t" << "rank" << endl;
    
    size_t pos;
    vector< pair<string , double> > pi_vec;
    string rs;
    double pi_temp;
    
    for (size_t i=0; i<ns_test; ++i) {
        
        // save the data along the order of all variants, snp_pos is sorted by order
        rs = snp_pos[i].rs;
        outfile<<rs<<"\t"<<snp_pos[i].chr<<"\t"<<snp_pos[i].bp<<"\t"<< snp_pos[i].a_major<<"\t"<<snp_pos[i].a_minor<<"\t" ;
        outfile << scientific << setprecision(6)  << snp_pos[i].maf << "\t";
        
        for (size_t j=0; j < n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                outfile << j << "\t";
                break;
            }
            else if(j == (n_type - 1)) outfile << "NA" << "\t";
        }
        
        pos = snp_pos[i].pos;
        //beta_g is saved by position
        if (beta_g[pos].second!=0) {
            pi_temp = beta_g[pos].second/(double)s_step;
            outfile << beta_g[pos].first/beta_g[pos].second<< "\t" << pi_temp <<"\t";
        }
        else {
            pi_temp = 0.0;
            outfile << 0.0 << "\t" << 0.0 << "\t";
        }
        pi_vec.push_back(make_pair(rs, pi_temp));
        
        if ( SNPorder_vec[i].first  != pos) {
            cerr << "ERROR: SNPorder_vec[i].first not equal to pos... \n";
            exit(-1);
        }
        outfile << scientific << setprecision(6) << Z_scores[pos] << "\t" << SE_beta[pos] << "\t" << pos_loglr[SNPorder_vec[i].second].second << "\t"<< pval_lrt[pos] << "\t" ;
        outfile << mapOrder2Rank[i];
        outfile << endl;
    }
    outfile.clear();	
    outfile.close();
}

void BVSRM::WriteFGWAS_InputFile(const vector<SNPPOS> &snp_pos, const vector<double> &Z_scores, const vector<double> &SE_beta)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".ifgwas";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    //header of the space-delimited text file
    outfile<<"SNPID"<<" "<<"CHR"<<" " <<"POS"<<" " << "F" << " " << "Z"<< " " << "SE" << " " << "N" << " " << "nonsynonymous" << endl;
    
    size_t pos;
    
    for (size_t i=0; i<ns_test; ++i) {
        // save the data along the order of all variants, snp_pos is sorted by order
        pos = snp_pos[i].pos;

        outfile<< snp_pos[i].rs <<" "<< snp_pos[i].chr<<" " <<snp_pos[i].bp << " ";
        
        outfile << scientific << setprecision(6)  << snp_pos[i].maf << " " << Z_scores[pos] << " " << SE_beta[pos] << " ";

        outfile << ni_test << " ";

        if(snp_pos[i].indicator_func[0]) outfile << 1 ;
        else outfile << 0;
        
        outfile << endl;
    }
    
    outfile.clear();    
    outfile.close();
}


void BVSRM::WriteGenotypeFile(uchar **X, const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".geno";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    //write header with VcfSampleID_test 
    outfile<<"ID"<<"\t"<<"CHROM"<<"\t" <<"POS"<<"\t" << "REF"<< "\t" << "ALT"  << "\t";
    for (size_t i=0; i<ni_test; i++) {
        if (i ==(ni_test-1)) {
            outfile << VcfSampleID_test[i] << endl;
        }
        else outfile << VcfSampleID_test[i] << "\t";
    } 

    size_t pos;
    string rs;
    double geno_j;
    uchar c;
    
    for (size_t i=0; i<ns_test; ++i) {
        
        pos = snp_pos[i].pos;
        // save the data along the order of all variants, snp_pos is sorted by order
        rs = snp_pos[i].rs;
        outfile<< rs <<"\t"<< snp_pos[i].chr<<"\t" <<snp_pos[i].bp << "\t";
        outfile << snp_pos[i].a_major << "\t" <<snp_pos[i].a_minor << "\t";
        
        /*for (size_t j=0; j < n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                outfile << j << "\t";
                break;
            }
            else if(j == (n_type - 1)) outfile << "NA" << "\t";
        }*/       
        //outfile << scientific << setprecision(6)  << snp_pos[i].maf << "\t";
        
        for (size_t j=0; j < ni_test; j++) {
            c = X[pos][j];
            geno_j = UcharTable[(int)c].second;

            if(geno_j < 0.0 || geno_j > 2.0){
                cout << "ERROR: genotype = " << geno_j << endl;
                exit(-1);
            }else{
                if (geno_j == 0.0) geno_j = 0;
                else if (geno_j == 2.0) geno_j = 2;
                else if (geno_j == 1.0) geno_j = 1;

                if (j == (ni_test-1))
                    outfile << fixed << setprecision(2)  << geno_j << endl;
                else
                    outfile << fixed << setprecision(2) << geno_j << "\t";
            }
        }
    }
    
    outfile.clear();
    outfile.close();
}

void BVSRM::CalcPgamma (double *p_gamma, size_t p_gamma_top)
{
    double p, q;

    if(p_gamma_top < 50){
        p_gamma_top=100;
    } else if(p_gamma_top > 300){
        p_gamma_top=300;
    }

  if((ns_test-p_gamma_top) < 0){
        p = 1.0 / double(ns_test);
        for (size_t i=0; i<ns_test; ++i) {
            p_gamma[i] = p;
        }
  } else{
    p = 0.9 / double(p_gamma_top);
    q = 0.1 / ((double)(ns_test-p_gamma_top));
    for (size_t i=0; i<ns_test; ++i) {
        if(i < p_gamma_top) p_gamma[i] = p;
            else p_gamma[i] = q;
    }
  }
    return;
}


//currently used
void BVSRM::SetXgamma (gsl_matrix *Xgamma, uchar **X, vector<size_t> &rank)
{
	size_t pos;
	for (size_t i=0; i<rank.size(); ++i) {
		pos=SNPrank_vec[rank[i]].first;
		gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
        getGTgslVec(X, &Xgamma_col.vector, pos, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
    }	
	return;
}

//currently used
void BVSRM::SetXgamma (uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, const vector<size_t> &rank_new, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    double d;
    // cout << "X_add set start" << endl;
    
    //rank_old and rank_new are sorted already inside PorposeGamma
    //calculate vectors rank_remove and rank_add
    //  size_t v_size=max(rank_old.size(), rank_new.size());
    //make sure that v_size is larger than repeat
    size_t v_size=20;
    vector<size_t> rank_remove(v_size), rank_add(v_size), rank_union(s_max+v_size);
    vector<size_t>::iterator it;
    
    it=set_difference (rank_old.begin(), rank_old.end(), rank_new.begin(), rank_new.end(), rank_remove.begin());
    rank_remove.resize(it-rank_remove.begin());
    stable_sort(rank_remove.begin(), rank_remove.end(), comp_vec);
   // cout << "rank_remove: "; PrintVector(rank_remove);
    
    it=set_difference (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_add.begin());
    rank_add.resize(it-rank_add.begin());
    stable_sort(rank_add.begin(), rank_add.end(), comp_vec);
  //  cout << "rank_add: "; PrintVector(rank_add);
    
    it=set_union (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_union.begin());
    rank_union.resize(it-rank_union.begin());
    stable_sort(rank_union.begin(), rank_union.end(), comp_vec);
   // cout << "rank_union: "; PrintVector(rank_union);
    
    //map rank_remove and rank_add
    map<size_t, int> mapRank2in_remove, mapRank2in_add;
    for (size_t i=0; i<rank_remove.size(); i++) {
        mapRank2in_remove[rank_remove[i]]=1;
    }
    for (size_t i=0; i<rank_add.size(); i++) {
        mapRank2in_add[rank_add[i]]=1;
    }
    
    //obtain the subset of matrix/vector
    gsl_matrix_const_view Xold_sub=gsl_matrix_const_submatrix(X_old, 0, 0, X_old->size1, rank_old.size());
    gsl_matrix_const_view XtXold_sub=gsl_matrix_const_submatrix(XtX_old, 0, 0, rank_old.size(), rank_old.size());
    gsl_vector_const_view Xtyold_sub=gsl_vector_const_subvector(Xty_old, 0, rank_old.size());
    
    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(X_new, 0, 0, X_new->size1, rank_new.size());
    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
    gsl_vector_view Xtynew_sub=gsl_vector_subvector(Xty_new, 0, rank_new.size());
    
    if (rank_remove.size()==0 && rank_add.size()==0) {
        gsl_matrix_memcpy(&Xnew_sub.matrix, &Xold_sub.matrix);
        gsl_matrix_memcpy(&XtXnew_sub.matrix, &XtXold_sub.matrix);
        gsl_vector_memcpy(&Xtynew_sub.vector, &Xtyold_sub.vector);
        //cout << "rank_old = rank_new; " << "Xgamma_new set success" << endl;
    } else {
        size_t i_old, j_old, i_new, j_new, i_add, j_add, i_flag, j_flag;
        if (rank_add.size()==0) {
            //only delete a snp
            i_old=0; i_new=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}
                
                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                
                d=gsl_vector_get (Xty_old, i_old);
                gsl_vector_set (Xty_new, i_new, d);
                
                j_old=i_old; j_new=i_new;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}
                    
                    d=gsl_matrix_get(XtX_old, i_old, j_old);
                    
                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}
                    
                    j_old++; j_new++;
                }
                i_old++; i_new++;
            }
            //cout << "X_add = NULL; " << "Xgamma_new set success" << endl;
        } else {
            //rank_add has length > 0
            gsl_matrix *X_add=gsl_matrix_alloc(X_old->size1, rank_add.size() );
            gsl_matrix *XtX_aa=gsl_matrix_alloc(X_add->size2, X_add->size2);
            gsl_matrix *XtX_ao=gsl_matrix_alloc(X_add->size2, rank_old.size());
            gsl_vector *Xty_add=gsl_vector_alloc(X_add->size2);
            
            //get X_add
            SetXgamma (X_add, X, rank_add);
            //cout << "X_add set success" << endl;
            
            //get t(X_add)X_add and t(X_add)X_temp
            clock_t time_start=clock();
            
            //somehow the lapack_dgemm does not work here
            //#ifdef WITH_LAPACK
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_add, 0.0, XtX_aa);
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_old, 0.0, XtX_ao);
            
            //#else
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, X_add, 0.0, XtX_aa);
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, &Xold_sub.matrix, 0.0, XtX_ao);
            //#endif
            gsl_blas_dgemv(CblasTrans, 1.0, X_add, y, 0.0, Xty_add);
            
            time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            
            //save to X_new, XtX_new and Xty_new
            i_old=0; i_new=0; i_add=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}
                if (mapRank2in_add.count(rank_new[i_new])!=0) {i_flag=1; //within x_add
                } else {i_flag=0; //within x_common
                }
                
                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                if (i_flag==1) {
                    gsl_vector_view Xcopy_col=gsl_matrix_column(X_add, i_add);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                } else {
                    gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                }
                //  cout << "Xgamma_new set success" << endl;
                
                if (i_flag==1) {
                    d=gsl_vector_get (Xty_add, i_add);
                } else {
                    d=gsl_vector_get (Xty_old, i_old);
                }
                gsl_vector_set (Xty_new, i_new, d);
                // cout << "Xty_new set success" << endl;
                
                j_old=i_old; j_new=i_new; j_add=i_add;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}
                    if (mapRank2in_add.count(rank_new[j_new])!=0) {j_flag=1;} else {j_flag=0;}
                    
                    if (i_flag==1 && j_flag==1) {
                        d=gsl_matrix_get(XtX_aa, i_add, j_add);          
                    } else if (i_flag==1 && j_flag==0) {
                        d=gsl_matrix_get(XtX_ao, i_add, j_old);
                    } else if (i_flag==0 && j_flag==1) {
                        d=gsl_matrix_get(XtX_ao, j_add, i_old);
                    } else {
                        d=gsl_matrix_get(XtX_old, i_old, j_old);
                    }
                    
                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}
                    
                    j_new++;
                    if (j_flag==1) {j_add++;} else {j_old++;}
                }
                //cout << "XtX_new success" << endl;
                i_new++; if (i_flag==1) {i_add++;} else {i_old++;}
            }
            // cout << "X_gamma set success" << endl;
            
            gsl_matrix_free(X_add);
            gsl_matrix_free(XtX_aa);
            gsl_matrix_free(XtX_ao);
            gsl_vector_free(Xty_add);
        }
        
    }
    
    rank_remove.clear();
    rank_add.clear();
    rank_union.clear();
    mapRank2in_remove.clear();
    mapRank2in_add.clear();
    
    return;
}


//currently used in propose gamma
void BVSRM::SetXgammaDel (const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const vector<size_t> &rank_old, size_t col_id, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    size_t s_size = rank_old.size();
    size_t s2;
    
    if (col_id==0) {
        s2 = s_size-1;
        gsl_matrix_const_view X2old = gsl_matrix_const_submatrix(X_old, 0, 1, ni_test, s2);
        gsl_matrix_const_view XtX22_sub = gsl_matrix_const_submatrix(XtX_old, 1, 1, s2, s2);
        gsl_vector_const_view Xty2_sub = gsl_vector_const_subvector(Xty_old, 1, s2);
        
        gsl_matrix_view Xnew2_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, s2);
        gsl_matrix_view XtXnew22_sub=gsl_matrix_submatrix(XtX_new, 0, 0, s2, s2);
        gsl_vector_view Xtynew2_sub=gsl_vector_subvector(Xty_new, 0, s2);
        
        gsl_matrix_memcpy(&Xnew2_sub.matrix, &X2old.matrix);
        gsl_matrix_memcpy(&XtXnew22_sub.matrix, &XtX22_sub.matrix);
        gsl_vector_memcpy(&Xtynew2_sub.vector, &Xty2_sub.vector);
    }
    else if(col_id == (s_size-1)){

        gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, col_id);
        gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, col_id, col_id);
        gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, col_id);
        
        gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, col_id);
        gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, col_id, col_id);
        gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, col_id);
        
        gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
        gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
        gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);
    }
    else{
        s2 = s_size - col_id - 1;
        gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, col_id);
        gsl_matrix_const_view X2old = gsl_matrix_const_submatrix(X_old, 0, col_id+1, ni_test, s2);
        
        gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, col_id, col_id);
        gsl_matrix_const_view XtX12_sub = gsl_matrix_const_submatrix(XtX_old, 0, col_id+1, col_id, s2);
        gsl_matrix_const_view XtX21_sub = gsl_matrix_const_submatrix(XtX_old, col_id+1, 0, s2, col_id);
        gsl_matrix_const_view XtX22_sub = gsl_matrix_const_submatrix(XtX_old, col_id+1, col_id+1, s2, s2);
        
        gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, col_id);
        gsl_vector_const_view Xty2_sub = gsl_vector_const_subvector(Xty_old, col_id+1, s2);
        
        gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, col_id);
        gsl_matrix_view Xnew2_sub=gsl_matrix_submatrix(X_new, 0, col_id, ni_test, s2);
        
        gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, col_id, col_id);
        gsl_matrix_view XtXnew12_sub=gsl_matrix_submatrix(XtX_new, 0, col_id, col_id, s2);
        gsl_matrix_view XtXnew21_sub=gsl_matrix_submatrix(XtX_new, col_id, 0, s2, col_id);
        gsl_matrix_view XtXnew22_sub=gsl_matrix_submatrix(XtX_new, col_id, col_id, s2, s2);
        
        gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, col_id);
        gsl_vector_view Xtynew2_sub=gsl_vector_subvector(Xty_new, col_id, s2);
        
        gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
        gsl_matrix_memcpy(&Xnew2_sub.matrix, &X2old.matrix);
        
        gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
        gsl_matrix_memcpy(&XtXnew12_sub.matrix, &XtX12_sub.matrix);
        gsl_matrix_memcpy(&XtXnew21_sub.matrix, &XtX21_sub.matrix);
        gsl_matrix_memcpy(&XtXnew22_sub.matrix, &XtX22_sub.matrix);
        
        gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);
        gsl_vector_memcpy(&Xtynew2_sub.vector, &Xty2_sub.vector);
    }

}

void BVSRM::SetXgammaAdd (uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, size_t ranki, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    double xty;
    size_t s_size = rank_old.size();
    size_t pos = SNPrank_vec[ranki].first;
    
    if (s_size==0) {
        cerr << "setXgammaAdd rank_old has size 0\n";
        exit(-1);
    }
    //copy rank_old
    gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, s_size);
    gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, s_size);
    
    gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, s_size);
    gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, s_size, s_size);
    gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, s_size);
    
    gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
    gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
    gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);
    
    //create ranki
    gsl_vector_view xvec = gsl_matrix_column(X_new, s_size);
    getGTgslVec(X, &xvec.vector, pos, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);

    gsl_vector_view Xtx_col = gsl_matrix_subcolumn(XtX_new, s_size, 0, s_size);
    gsl_vector_view Xtx_row = gsl_matrix_subrow(XtX_new, s_size, 0, s_size);
    gsl_blas_dgemv(CblasTrans, 1.0, &X1old.matrix, &xvec.vector, 0.0, &Xtx_col.vector);
    gsl_vector_memcpy(&Xtx_row.vector, &Xtx_col.vector);
    gsl_matrix_set(XtX_new, s_size, s_size, XtX_diagvec[pos]);
    
    gsl_blas_ddot(&xvec.vector, y, &xty);
    gsl_vector_set(Xty_new, s_size, xty);
    
}


//end of currently used

double BVSRM::CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2) 
{
	double pve, var_y;	
	
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *Xty=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *OiXty=gsl_vector_alloc (UtXgamma->size2);

	gsl_matrix_set_identity (Omega);
	gsl_matrix_scale (Omega, 1.0/sigma_a2); 

#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", 1.0, UtXgamma, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, UtXgamma, UtXgamma, 1.0, Omega);	
#endif
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma, Uty, 0.0, Xty);

	CholeskySolve(Omega, Xty, OiXty);
	
	gsl_blas_ddot (Xty, OiXty, &pve);
	gsl_blas_ddot (Uty, Uty, &var_y);
	
	pve/=var_y;
	
	gsl_matrix_free (Omega);
	gsl_vector_free (Xty);
	gsl_vector_free (OiXty);

	return pve;
}

bool comp_snp(const SNPPOS& lhs, const SNPPOS& rhs){
    return (lhs.chr.compare(rhs.chr) < 0) || ((lhs.chr.compare(rhs.chr) == 0) && (lhs.bp < rhs.bp));
}


void BVSRM::setHyp(double theta_temp, double subvar_temp){
        
    // Default initial values   
    cout << "rv from command line = " << rv << endl;
    theta.assign(n_type, 1.0e-6);
    subvar.assign(n_type, 10.0);

    //cout << "load fixed hyper parameter values from : " << hypfile << endl;
    string line;
    char *pch, *nch;
    size_t group_idx=0;
    
    if (hypfile.empty()) {
        cout << "Did not specify hypefile, use default values ...\n";
    }
    else{
        ifstream infile(hypfile.c_str(), ifstream::in);
        if(!infile) {
            cout << "Error opening file " << hypfile << endl; 
            exit(-1);
        }
        // cout << "load hyp from hypfile... " << hypfile << endl;
        while (!safeGetline(infile, line).eof()) {
            if ((line[0] == 'h') || (line[0] == '#') || (line[0] == 'p')) {
                continue;
            }
            else{
                pch = (char *)line.c_str();
                nch = strchr(pch, '\t');
                theta[group_idx] = strtod(pch, NULL);
                if(nch == NULL){
                 cerr << "Need input initial hyper parameter value for sigma2 \n";
                    } else{pch = nch+1;}
                if(group_idx < n_type)  {
                    subvar[group_idx] = strtod(pch, NULL);
                }
                group_idx++;
            }
        }
        infile.clear();
        infile.close();
    }

    log_theta.clear();
    log_qtheta.clear();
    for(size_t i=0; i < n_type; i++){
        log_theta.push_back(log(theta[i]));
        log_qtheta.push_back(log(1.0 - theta[i]));
    }
    
    cout << "Initial causal probability per category = "; PrintVector(theta);
    cout << "Initial effect-size variance per category = "; PrintVector(subvar);
    //cout << "log_qtheta: "; PrintVector(log_qtheta); 

}


//InitialMCMC currently used
void BVSRM::InitialMCMC (uchar **X, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos)

{
    //double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
    double q_genome=gsl_cdf_chisq_Qinv(5e-8, 1);
    //cout << "significant chisquare value : " << q_genome << endl;
    cHyp.n_gamma=0;
    for (size_t i=0; i<pos_loglr.size(); ++i) {
        if (2.0*pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
    }
    //cout << "number of snps before adjust = " << cHyp.n_gamma << endl;
    if (cHyp.n_gamma<30) {cHyp.n_gamma=30;}
    if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
    if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}
    
    
    if (!iniSNPfile.empty() && iniType == 0) {
        
        ifstream infile(iniSNPfile.c_str(), ifstream::in);
        if(!infile) {cout << "Error opening file " << iniSNPfile << endl; exit(-1);}
        string lineid;
        rank.clear();
        size_t orderj, rankj;
        
        cout << "Start loading initial snp IDs from " << iniSNPfile << "\n";
        
        while (!safeGetline(infile, lineid).eof()) {
            
            orderj = 0;
            for (size_t i=0; i < snp_pos.size(); i++) {
                if (snp_pos[i].rs.compare(lineid) == 0) {
                    orderj=i;
                    rankj = SNPorder_vec[orderj].second;
                    rank.push_back(rankj);
                    //cout << lineid << " with rank = " << rankj;
                    //snp_pos[orderj].printMarker();
                    break;
                }
            }
        }
        infile.close();
        infile.clear();
        
        if (rank.size() == 0) {
            for (size_t i=0; i<cHyp.n_gamma; ++i) {
                rank.push_back(i);
            }
        } //take rank 0 if tracked no SNPs from the iniSNPfile
        
        cHyp.n_gamma = rank.size();
    }
    else if(iniType == 0) {iniType = 1;}
    else if(iniType == 1) {
        cout << "Start with top variants.\n";
        rank.clear();
        for (size_t i=0; i<cHyp.n_gamma; ++i) {
            rank.push_back(i);
        }
    } // Start with most significant variants from SVT
    else if(iniType == 2) {
        cout << "Start with top SVT SNPs/ColinearTest.\n";
        size_t posr, j=0, i=0;
        double xtx;
        
        rank.clear();
        rank.push_back(0);
        posr = SNPrank_vec[0].first;
        xtx = XtX_diagvec[posr];
        // cout << "rank added: " << 0 << ", ";
        
        gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
        gsl_vector * xvec = gsl_vector_alloc(ni_test);
        getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); //get geno column
        
        gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
        gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);
        
        do{
       // for (size_t i=1; i < cHyp.n_gamma; ++i){
            gsl_matrix_set_col(Xr, i, xvec);
            gsl_matrix_set(XtXr, i, i, xtx);
            
            if (i>0) {
                gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
                gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
                gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);
                
                gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, i, 0, i);
                gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, i, 0, i);
                gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
                gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);
                
            }
            
            if ((rank.size() < cHyp.n_gamma)) {
            
                do{
                    j++; // Consider rank j
                    //cout << "consider rank j" << j << endl;
                } while( (j < ns_test) && ColinearTest(X, Xr, XtXr, j, rank.size()));
                
                rank.push_back(j);
                posr = SNPrank_vec[j].first;
                xtx = XtX_diagvec[posr];
                getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); // get geno column
            }
            i++;
            
        }while(i < (cHyp.n_gamma));
        
        //PrintMatrix(XtXr, cHyp.n_gamma, cHyp.n_gamma);
        
        gsl_matrix_free(Xr);
        gsl_matrix_free(XtXr);
        gsl_vector_free(xvec);
        gsl_vector_free(Xtxvec);

    } // Start with most significant variants from SVT
    else if(iniType == 3){
        cout << "Start with Step-wise selected variants. \n";
    vector<pair<size_t, double> > rank_loglr;
    size_t posr, radd;

    size_t topMarkers=500;
    if(ns_test<500){topMarkers = ns_test;}
        
    for (size_t i=1; i<topMarkers; ++i) {
        rank_loglr.push_back(make_pair(i, pos_loglr[i].second));
    }
    cout << endl;
    
    rank.clear();
    rank.push_back(0);
    posr = SNPrank_vec[0].first;
    //cout << "rank added: " << 0 << " with LRT "<< pos_loglr[0].second << "," ;

    gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
    gsl_vector * xvec = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); //get geno column
    
    gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
    gsl_vector * Xtyr = gsl_vector_alloc(cHyp.n_gamma);
    
    gsl_vector * yres = gsl_vector_alloc(ni_test);
    gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);
    
    double xty, yty;
    gsl_blas_ddot(Uty, Uty, &yty);
    
    for (size_t i=1; i < s_max; ++i){
        gsl_matrix_set_col(Xr, (i-1), xvec);
        gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
        
        gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
        gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);
        
        gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, (i-1), 0, i);
        gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, (i-1), 0, i);
        gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
        gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);
        
        gsl_blas_ddot(xvec, Uty, &xty);
        gsl_vector_set(Xtyr, (i-1), xty);
        
        // calculate conditional yres
        CalcRes(Xr, Uty, XtXr, Xtyr, yres, i, yty);
        for (size_t j=0; j<rank_loglr.size(); ++j) {
            posr = SNPrank_vec[rank_loglr[j].first].first;
            getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); // get geno column
            rank_loglr[j].second = CalcLR(yres, xvec, posr);
        }
        stable_sort (rank_loglr.begin(), rank_loglr.end(), comp_lr); //sort the initial rank.
        
        if (rank_loglr[0].second > q_genome) {
            radd = rank_loglr[0].first;
           // cout << "XtXr : "; PrintMatrix(XtXr, rank.size(), rank.size());
            if (ColinearTest(X, Xr, XtXr, radd, rank.size())) {
                continue;
            }
            else {
                posr = SNPrank_vec[radd].first;
                getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
                rank.push_back(radd);
                rank_loglr.erase(rank_loglr.begin());
                //cout << "rank added: " << radd << " with LRT "<< rank_loglr[0].second << "," ;
            }
        }
        else break;
    }
        cHyp.n_gamma = rank.size();
        //cout <<"Initial XtX: \n"; PrintMatrix(XtXr, cHyp.n_gamma, cHyp.n_gamma);
        gsl_matrix_free(Xr);
        gsl_matrix_free(XtXr);
        gsl_vector_free(Xtyr);
        gsl_vector_free(xvec);
        gsl_vector_free(yres);
        gsl_vector_free(Xtxvec);
    }
    //cout << "number of snps = " << cHyp.n_gamma << endl;
    //stable_sort (rank.begin(), rank.end(), comp_vec); //sort the initial rank.
    cout << "Starting model has variants with ranks: \n"; PrintVector(rank);
    
    cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
    cHyp.h=pve_null;
    
    if (cHyp.logp==0) {cHyp.logp=-0.000001;}
    if (cHyp.h==0) {cHyp.h=0.1;}
    
    gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
    SetXgamma (UtXgamma, X, rank);
   // cout<<"initial value of h = "<< h <<endl;
    
    double sigma_a2;
    if (trace_G!=0) {
        sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp));
    } else {
        sigma_a2=cHyp.h*1.0/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
    }
    if (sigma_a2==0) {sigma_a2=0.025;}
   // cout << "initial sigma_a2 = " << sigma_a2 << endl;
    
    //cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
    gsl_matrix_free (UtXgamma);
    
    if (cHyp.rho>1.0) {cHyp.rho=1.0;}
    if (cHyp.h<h_min) {cHyp.h=h_min;}
    if (cHyp.h>h_max) {cHyp.h=h_max;}
    //if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
    //if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
    if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
    if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}
    
    //cout << "start setHyp... \n";
    setHyp(((double)cHyp.n_gamma/(double)ns_test), sigma_a2);
    cHyp.theta = theta;
    cHyp.log_theta = log_theta;
    cHyp.subvar = subvar; // initial subvar vector
    
    // cout<<"initial value of h = "<< h <<endl;
    // cout<<"initial value of rho = "<<cHyp.rho<<endl;
    // cout<<"initial value of theta_vec = "; PrintVector(theta);
    // cout << "initial value of sub-variance_vec = "; PrintVector(subvar);
    cout<<"Initially selected number of variants in the model = "<<cHyp.n_gamma<<endl;
    
    return;
}

// for the new model
double BVSRM::CalcPosterior (const double yty, class HYPBSLMM &cHyp)
{
    double logpost=0.0;
    
    //for quantitative traits, calculate pve and pge
    //pve and pge for case/control data are calculted in CalcCC_PVEnZ
    if (a_mode==11) {
        cHyp.pve=0.0;
     //   cHyp.pge=1.0;
    }
    //calculate likelihood
   // if (a_mode==11) {logpost-=0.5*(double)ni_test*log(yty);}
   // else {logpost-=0.5*yty;}
    
    // calculate gamma likelihood
    //logpost += CalcLikegamma(cHyp);
    
    return logpost;
}

//Used for EM-Block
void BVSRM::getSubVec(gsl_vector *sigma_subvec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    size_t order_i;
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = SNPrank_vec[rank[i]].second;
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                gsl_vector_set(sigma_subvec, i, subvar[j]);
                continue;
            }
        }
    }
}

//set sigma_subvec and mgamma vectoer and trace vector Gvec

//used in EM-block
void BVSRM::set_mgamma(class HYPBSLMM &cHyp, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    size_t order_i;
    
    cHyp.m_gamma.assign(n_type, 0);
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = SNPrank_vec[rank[i]].second;
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                cHyp.m_gamma[j]++;
                continue;
            }
        }
    }
}


//used in EM_BLock
double BVSRM::CalcLikegamma(const class HYPBSLMM &cHyp)
{
    double loglikegamma = 0.0;
    
    for (size_t i=0; i < n_type; i++) {
        loglikegamma +=  log_theta[i] * ((double)cHyp.m_gamma[i]) + ((double)(mFunc[i] - cHyp.m_gamma[i])) * log_qtheta[i];
    }
    return loglikegamma;
}

// for the new model
double BVSRM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag, double &loglike)
{
    //conditioning on hyper parameters: subvar, log_theta
    double logpost=0.0;
    double d;
    double logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;
    
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);
    
    //calculate Omega
    gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
   // cout << "Omega : "; PrintMatrix(Omega, 5, 5);
    CalcXVbeta(Omega, &sigma_sub.vector);
    gsl_vector_view Omega_diag = gsl_matrix_diagonal(Omega);
    gsl_vector_add_constant(&Omega_diag.vector, 1.0);
    
    if(LapackSolve(Omega, &Xty_sub.vector, beta_hat)!=0)
       EigenSolve(Omega, &Xty_sub.vector, beta_hat);
    logdet_O=LapackLogDet(Omega);
    
    //cout << "beta_hat from solve : "; PrintVector(beta_hat);
    gsl_vector_mul(beta_hat, &sigma_sub.vector);
    //cout << "beta_hat: "; PrintVector(beta_hat);
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
    
    double bxy;
    gsl_blas_ddot (&Xty_sub.vector, beta_hat, &bxy);
    double R2 = bxy / yty;
    
    
   /* double lambda = 0.0;
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(Omega, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.01;

    int k=0;*/
     
    if (R2 > 1.0 || R2 < -0.0) {
        
        //cout << "R2 in CalcPosterior = " << R2 << endl;
        Error_Flag=1;
        /*
        WriteMatrix(&Xgamma_sub.matrix, "_X");
        WriteMatrix(&XtX_sub.matrix, "_XtX");
        WriteMatrix(Omega, "_Omega");
        WriteVector(&sigma_sub.vector, "_sigma");
        WriteVector(&Xty_sub.vector, "_Xty");
        WriteVector(beta_hat, "_beta");
        // Error_Flag = 1;
        //cerr << "Error in calcPosterior: P_yy = " << P_yy << endl;
        //cout << "est beta_hat: "; PrintVector(beta_hat);
        // exit(-1);
        
        while (R2 > 1.1 || R2 < -0.1) {
            
            cout << "add lambda : negative R2 in calcposterior..." << R2 << "; k = " << k<< endl;
            
         gsl_vector_add_constant(&Omega_diag.vector, lambda);
         if(LapackSolve(Omega, &Xty_sub.vector, beta_hat)!=0)
             EigenSolve(Omega, &Xty_sub.vector, beta_hat);
         gsl_blas_ddot (&Xty_sub.vector, beta_hat, &bxy);
            R2 = bxy / yty;
            k++;
            cout << "now R2 = " << R2 << "; k = " << k << endl;
        
            if (k > 9) {
                cout << "reached k = " << k << endl;
                break;
            }
        }
        
        if (R2 > 1.0 || R2 < -0.0) {
            bxy=0.0;
            cHyp.pve=0.0;
            
            gsl_vector_set_zero(&beta_sub.vector);
            Error_Flag=1;
        }
        else{
            Error_Flag=0;
            gsl_vector_memcpy(&beta_sub.vector, beta_hat);
            gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);
            gsl_blas_ddot (Xb, Xb, &d);
            if (a_mode==11) {
                cHyp.pve=d/(double)ni_test;
                //cHyp.pve/=cHyp.pve+1.0/tau;
                // cHyp.pge=1.0;
            }
        }*/
    }

    else{
        Error_Flag=0;
        gsl_vector_memcpy(&beta_sub.vector, beta_hat);
        gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);
        gsl_blas_ddot (Xb, Xb, &d);
        if (a_mode==11) {
            cHyp.pve=d/(double)ni_test;
            //cHyp.pve/=cHyp.pve+1.0/tau;
            // cHyp.pge=1.0;
        }
    }
    
    logpost = tau * bxy;
    loglike = -0.5 * ((double)cHyp.n_gamma * logrv + (double)cHyp.m_gamma[0] * log_subvar[0] + (double)cHyp.m_gamma[1] * log_subvar[1] - logpost);
    
    logpost = -0.5 * (logdet_O - logpost);

    gsl_matrix_free (Omega);
    gsl_vector_free (beta_hat);
    
    return logpost;
}

// for the new model
//calculate likelihood P(Y | gamma, subvar, theta)
double BVSRM::CalcLikelihood (const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag)
{
    //double sigma_a2=cHyp.h/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test * trace_G);
    
    double loglike=0.0;
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;
    
  if (s_size == 0) {
        //calculate likelihood if ngamma=0
        if (a_mode==11) {loglike-=0.5*(double)ni_test*log(yty);}
        else {loglike-=0.5*yty;}
  }
    
  else{
   // gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
      gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
      gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
      gsl_vector *beta_hat=gsl_vector_alloc (s_size);
      gsl_vector *Xty_temp=gsl_vector_alloc (s_size);
      
      //calculate Omega
      gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
      CalcXVbeta(Omega, &sigma_sub.vector);
      gsl_matrix_set_identity (M_temp);
      gsl_matrix_add (Omega, M_temp);
      
      //calculate beta_hat
      gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
      logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);	//solve Omega * beta_hat = Xty for beta_hat
      // Omega was inverted here
      // logdet_0 = det(Omega)
      gsl_vector_mul(beta_hat, &sigma_sub.vector);
      gsl_blas_ddot (Xty_temp, beta_hat, &d);
      P_yy-=d;
      if (P_yy <= 0) {
          Error_Flag = 1;
          cout << "Error in calclikelihood: P_yy = " << P_yy << endl;
          cout << "h = "<<setprecision(6) << cHyp.h << "; rho: " << cHyp.rho_vec[0] << ", " << cHyp.rho_vec[1];
          cout << "theta: "<<setprecision(6) << exp(cHyp.log_theta[0]) << ", " << exp(cHyp.log_theta[1]);
          cout << "; subvar: "<<setprecision(6) << cHyp.subvar[0] << ", " << cHyp.subvar[1];
          cout << "beta_hat: "; PrintVector(beta_hat);
          
          cout << "set beta_hat to 0\n";
          gsl_vector_set_zero(beta_hat);
          P_yy = yty;
         // exit(-1);
      }
      else {Error_Flag = 0;}
    
    loglike = -0.5 * logdet_O;
    if (a_mode==11) {loglike -= 0.5 * (double)ni_test * log(P_yy);}
    else {loglike -= 0.5*P_yy;}
    
    gsl_matrix_free (Omega);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (Xty_temp);
  }
    
    return loglike;
}


//currently used
double BVSRM::ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, uchar **X, const gsl_vector *z, const gsl_matrix *Xgamma_old, const gsl_matrix *XtX_old, const gsl_vector *Xtz_old, const double &ztz, int &flag_gamma, gsl_matrix *Xgamma_new, gsl_matrix *XtX_new, gsl_vector *Xtz_new)
{
    map<size_t, int> mapRank2in;
    double unif, logp = 0.0;
    size_t r_add, r_remove, col_id, r;
    
    rank_new.clear();
    if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"size wrong"<<endl;}
    if (cHyp_old.n_gamma!=0) {
        for (size_t i=0; i<rank_old.size(); ++i) {
            r=rank_old[i];
            rank_new.push_back(r);
            mapRank2in[r]=1;
        }
    }
    cHyp_new.n_gamma=cHyp_old.n_gamma;
    
    //for (size_t i=0; i<repeat; ++i) {
        
        unif=gsl_rng_uniform(gsl_r);
        
        if (unif < 0.33 && cHyp_new.n_gamma<s_max) {flag_gamma=1;}
        else if (unif>=0.33 && unif < 0.67 && cHyp_new.n_gamma>s_min) {flag_gamma=2;}
        else if (unif>=0.67 && cHyp_new.n_gamma>0 && cHyp_new.n_gamma<ns_test) {flag_gamma=3;}
        else {flag_gamma=0;}
        
        if(flag_gamma==1)  {//add a snp;
            // cout << "add a snp" << endl;
            
          /*  if (rank_old.size() > 0) {
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0) \
                    || (ColinearTest(X, Xgamma_old, XtX_old, r_add, rank_old.size()))
                    );
            }
            else {*/
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0));
           // }
            
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            
            mapRank2in[r_add]=1;
            rank_new.push_back(r_add);
            cHyp_new.n_gamma++;
            logp += -log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);
            
            if (rank_old.size()>0) {
                SetXgammaAdd(X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, r_add, Xgamma_new, XtX_new, Xtz_new);
            }
            else{
                SetXgamma (Xgamma_new, X, rank_new);
                CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
            }
           // cout << "XtX from set calcXtX: \n"; PrintMatrix(XtX_new, rank_new.size(), rank_new.size());
           // cout << "succesfully added a snp" << endl;

        }
        else if (flag_gamma==2) {//delete a snp;
            //  cout << "delete a snp" << endl;
            
            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];
            
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            prob_total+=p_gamma[r_remove];
            
            mapRank2in.erase(r_remove);
            rank_new.erase(rank_new.begin()+col_id);
            logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
            cHyp_new.n_gamma--;
            
            if (rank_new.size() > 0) {
                SetXgammaDel(Xgamma_old, XtX_old, Xtz_old, rank_old, col_id, Xgamma_new, XtX_new, Xtz_new);
            }
            // cout << "succesfully deleted a snp" << endl;
        }
        else if (flag_gamma==3) {//switch a snp;
            // cout << "switch a snp" << endl;
            long int o_add, o_remove;
            long int o_rj, o_aj;
            size_t j_add, j_remove, o;
            
            gsl_ran_discrete_t *gsl_s, *gsl_a; //JY added dynamic gsl_s
            double *p_BFr = new double[ns_neib];
            double *p_BFa = new double[ns_neib];
            
            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];//careful with the proposal
            if(mapRank2in.count(r_remove) == 0) {cout << "wrong proposal of r_remove;" << endl; exit(-1);}
            o_remove = SNPrank_vec[r_remove].second;
            rank_new.erase(rank_new.begin()+col_id);
            size_t s_size = rank_new.size();
            mapRank2in.erase(r_remove);
            
            //cout << "s_size = "<< s_size << "; s_size+1 = " << s_size+1 << endl;
            gsl_matrix *Xgamma_temp=gsl_matrix_alloc (ni_test, s_size+1);
            gsl_matrix *XtX_gamma=gsl_matrix_alloc (s_size+1, s_size+1);
            gsl_vector *Xtz_gamma=gsl_vector_alloc (s_size+1);
            gsl_vector *z_res = gsl_vector_alloc(ni_test);
            
            //cout << "Switch step rank_old:"; PrintVector(rank_old);
            //cout <<"XtX_old: "; PrintMatrix(XtX_old, rank_old.size(), rank_old.size());
            //cout << "temp rank_new:"; PrintVector(rank_new);
            if (s_size > 0) {
                //SetXgamma (Xgamma_temp, X, rank_new);
                //CalcXtX (Xgamma_temp, z, s_size, XtX_gamma, Xtz_gamma);
                //cout << "XtX from set calcXtX: \n"; PrintMatrix(XtX_gamma, s_size, s_size);
                
                SetXgammaDel(Xgamma_old, XtX_old, Xtz_old, rank_old, col_id, Xgamma_temp, XtX_gamma, Xtz_gamma);
                //cout << "XtX from set XgammaDel success: \n";
                //PrintMatrix(XtX_gamma, s_size, s_size);
                
                CalcRes(Xgamma_temp, z, XtX_gamma, Xtz_gamma, z_res, s_size, ztz);
                gsl_s = MakeProposal(o_remove, p_BFr, X, z_res, mapRank2in);
            }
            else {
                gsl_s = MakeProposal(o_remove, p_BFr, X, z, mapRank2in);
            }
            
            j_add = gsl_ran_discrete(gsl_r, gsl_s);
            o_add = (o_remove - win) + j_add;
            if((o_add < 0) || (o_add >= (long int)ns_test) || (o_add == (long int)o_remove))
                cout << "ERROR proposing switch snp"; //new snp != removed snp
            r_add = SNPorder_vec[(size_t)o_add].second;
            
            //cout << "XtX from set Xgamma: \n"; PrintMatrix(XtX_gamma, s_size, s_size);
            if (s_size>0) {
                
                /*if (ColinearTest(X, Xgamma_temp, XtX_gamma, r_add, s_size)) {
                    flag_gamma=-1;
                    //cout << "Failed colinear test in switch" << endl;
                }
                else{*/
                    gsl_a = MakeProposal(o_add, p_BFa, X, z_res, mapRank2in);
                    
                    double prob_total_remove=1.0;
                    double prob_total_add=1.0;
                    
                    for (size_t ii=0; ii<rank_new.size(); ++ii) {
                        r = rank_new[ii];
                        o = SNPrank_vec[r].second;
                        o_rj = ((long int)o - o_remove) + win;
                        o_aj = ((long int)o - o_add) + win;
                        if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
                        if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
                    }
                    
                    j_remove = o_remove - o_add + win;
                    logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                    logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)
                    
                    SetXgammaAdd(X, Xgamma_temp, XtX_gamma, Xtz_gamma, z, rank_new, r_add, Xgamma_new, XtX_new, Xtz_new);
                   // cout << "XtX from setXgammaAdd success: \n";
                    //PrintMatrix(XtX_new, s_size+1, s_size+1);
                    
                    mapRank2in[r_add]=1;
                    rank_new.push_back(r_add);
                    
                    gsl_ran_discrete_free(gsl_a);
                //}
            }
            else{
                //construct gsl_s, JY
                //cout << "o_add = " << o_add <<  "; r_add = "<<r_add << endl;
                gsl_a = MakeProposal(o_add, p_BFa, X, z, mapRank2in);
                
                double prob_total_remove=1.0;
                double prob_total_add=1.0;
                
                for (size_t ii=0; ii<rank_new.size(); ++ii) {
                    r = rank_new[ii];
                    o = SNPrank_vec[r].second;
                    o_rj = ((long int)o - o_remove) + win;
                    o_aj = ((long int)o - o_add) + win;
                    if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
                    if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
                }
                
                j_remove = o_remove - o_add + win;
                logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)
                
                mapRank2in[r_add]=1;
                rank_new.push_back(r_add);
                SetXgamma (Xgamma_new, X, rank_new);
                CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
                
                gsl_ran_discrete_free(gsl_a);
            }
            
            gsl_matrix_free(Xgamma_temp);
            gsl_matrix_free(XtX_gamma);
            gsl_vector_free(Xtz_gamma);
            gsl_vector_free(z_res);
            
            gsl_ran_discrete_free(gsl_s);
            
            delete[] p_BFr;
            delete[] p_BFa;
            // cout << "successfully switched a snp" << endl;
        }
        
        else {logp+=0.0;}//do not change
    //}
    mapRank2in.clear();
    return logp;
}


void BVSRM::WriteHyptemp(gsl_vector *LnPost, vector<double> &em_gamma){
    
    double em_logpost = 0.0, logpost_max =  gsl_vector_max(LnPost);
    //cout << "logpost_max = " << logpost_max << endl;
    for (size_t i=0; i < s_step; i++) {
        em_logpost += exp(gsl_vector_get(LnPost, i) - logpost_max);
    }
    em_logpost /= double(s_step);
    em_logpost = log(em_logpost) + logpost_max;
    
    //save E(file_out, lnpost, GV, rv, n[i], Gvec[i], m[i], sigma2[i])
    string file_hyp;
    file_hyp = "./output/" + file_out;
    file_hyp += ".hyptemp";

    ofstream outfile_hyp;

    // write *.hyptemp
    outfile_hyp.open (file_hyp.c_str(), ofstream::out);
    if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}

    outfile_hyp << file_out << "\t";
    outfile_hyp << scientific << setprecision(6) << em_logpost << "\t" << (GV / (double)s_step) << "\t" << rv ;
    for(size_t i=0; i < n_type; i++){
        outfile_hyp << "\t" << mFunc[i] ;
        outfile_hyp << scientific << setprecision(6) << "\t" << Gvec[i] ;
        outfile_hyp << "\t" << (em_gamma[i] / (double)s_step) ;
        if(em_gamma[i] > 0)
            {outfile_hyp << "\t" << sumbeta2[i] / em_gamma[i] ;}
        else {outfile_hyp << "\t" << sumbeta2[i];}
    }
    outfile_hyp << endl;

    outfile_hyp.clear();
    outfile_hyp.close();
    
}



//Current MCMC function
void BVSRM::MCMC (uchar **X, const gsl_vector *y, bool original_method) {
    
    if (original_method) {
        cout << "Run MCMC...\n";
    }
    // cout << "# of unique function types = " << n_type << endl;
    
    //new model related
    gsl_vector *sigma_subvec_old = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_old);
    gsl_vector *sigma_subvec_new = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_new);
    gsl_vector *LnPost = gsl_vector_alloc(s_step); //save logPost...

    vector<double> em_gamma(n_type, 0.0); 
        //save sum of m_q, beta2_q;
    GV = 0.0; 
    sumbeta2.assign(n_type, 0.0);
    
    gsl_vector *Xb_new=gsl_vector_alloc (ni_test);
    gsl_vector *Xb_old=gsl_vector_alloc (ni_test);
    gsl_vector *z_hat=gsl_vector_alloc (ni_test);
    gsl_vector *z=gsl_vector_alloc (ni_test);
    
    gsl_matrix *Xgamma_old=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_old=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_old=gsl_vector_alloc (s_max);
    gsl_vector *beta_old=gsl_vector_alloc (s_max);
    
    gsl_matrix *Xgamma_new=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_new=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_new=gsl_vector_alloc (s_max);
    gsl_vector *beta_new=gsl_vector_alloc (s_max);
    
    double ztz=0.0;
    gsl_vector_memcpy (z, y);
    double mean_z = CenterVector (z); // center phenotype in case

    gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    pheno_var = ztz / ((double)(ni_test-1)) ;

    //cout << "ztz = " << ztz << "; phenotype variance = " << ztz / pheno_var <<endl;
    gsl_vector_scale(z, 1.0 / sqrt(pheno_var)); // standardize phenotype z
    gsl_blas_ddot(z, z, &ztz); // calculate ztz for one more time after standardization
    
    //Initialize variables for MH
    double logPost_new, logPost_old, loglike_new, loglike_old, loglikegamma;
    double logMHratio;
    vector<size_t> rank_new, rank_old;
    class HYPBSLMM cHyp_old, cHyp_new;
    bool Error_Flag=0;
    
    if (a_mode==13) {
        pheno_mean=0.0;
    }
    vector<pair<double, double> > beta_g; //save beta estimates
    for (size_t i=0; i<ns_test; i++) {
        beta_g.push_back(make_pair(0.0, 0.0));
    }
        
    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);
    
    // Jingjing add a vector of "snpPos" structs snp_pos
    vector<SNPPOS> snp_pos;
    CreateSnpPosVec(snp_pos); //ordered by position here
    //cout << "1 / SNP_sd  = "; PrintVector(SNPsd, 10);
    
    vector<pair<size_t, double> > pos_loglr;
    vector<double> pval_lrt;
    vector<double> Z_scores;
    vector<double> SE_beta;

    cout << "Calculating Z_scores, standard errors of effect-sizes, LRT statistics, pvals ... \n";
    MatrixCalcLmLR (X, z, pos_loglr, ns_test, ni_test, SNPsd, SNPmean, Gvec, XtX_diagvec, Z_scores, SE_beta, pval_lrt, snp_pos, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); //calculate trace_G or Gvec, Z_scores, SE_beta
 //calculate trace_G or Gvec
    trace_G = VectorSum(Gvec) / double(ns_test);
    cout << "Trace of Genotype Matrix = " << trace_G << endl;

    for(size_t i=0; i < n_type; i++){
        Gvec[i] /= mFunc[i];
    } 
    //cout << "Avg trace of genotypes per category : " ; PrintVector(Gvec); 

    //calculat LD matrix and save in the output
    if(saveLD){
        cout << "start calculating LD matrix ... \n";
        gsl_matrix *LD = gsl_matrix_alloc(ns_test, ns_test);
        CalcLD(X, LD);
        cout << "start saving LD matrix ... \n";
        WriteMatrix(LD, ".LD");
        gsl_matrix_free(LD);
    }
    // end of calculate LD matrix
    
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // order snp_pos by chr/bp
    stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio
    
    //Jingjing's edit, create maps between rank and order
    size_t pos;
    for (size_t i=0; i<ns_test; ++i) {
        mapRank2pos[i]=pos_loglr[i].first;
        mapPos2Rank[pos_loglr[i].first] = i;
        
        mapOrder2pos[i] = snp_pos[i].pos;
        mapPos2Order[snp_pos[i].pos] = i;
    }
    
    for (size_t i=0; i<ns_test; ++i) {
        pos = mapRank2pos[i];
        mapRank2Order[i] = mapPos2Order[pos];
        
        pos = mapOrder2pos[i];
        mapOrder2Rank[i] = mapPos2Rank[pos];
    }
    
    SNPorder_vec.clear();
    SNPrank_vec.clear();
    for (size_t i=0; i<ns_test; i++) {
        SNPorder_vec.push_back(make_pair(snp_pos[i].pos, mapOrder2Rank[i]));
        SNPrank_vec.push_back(make_pair(pos_loglr[i].first, mapRank2Order[i]));
    }
    
    
    //end of Jingjing's edit
    
    //Calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
    gsl_rng_env_setup();
    const gsl_rng_type * gslType;
    gslType = gsl_rng_default;
    if (randseed<0)
    {
        time_t rawtime;
        time (&rawtime);
        tm * ptm = gmtime (&rawtime);
        
        randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
    }
    gsl_r = gsl_rng_alloc(gslType);
    gsl_rng_set(gsl_r, randseed);
    double *p_gamma = new double[ns_test];

    size_t p_gamma_top=0;
    for(size_t i=0; i < pval_lrt.size(); i++){
        if (pval_lrt[i] < 5e-8) p_gamma_top++;
    }
    cout << "Number of variants with p-value < 5e-8 : " << p_gamma_top << endl;

    // CalcPgamma (p_gamma); // calculate discrete distribution for gamma
    CalcPgamma (p_gamma, p_gamma_top);

    gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma); // set up proposal function for gamma
    
    //Initial parameters
    cout << "Start initializing MCMC ... \n";
    InitialMCMC (X, z, rank_old, cHyp_old, pos_loglr, snp_pos); // Initialize rank and cHyp
   
    //cout<< "Residual Variance Proportion = " << rv << endl;
    rv *= ztz / ((double)(ni_test-1)); tau = 1.0 / rv; // assume residual varaince = 80% phenotype variance
    logrv = log(2.0 * M_PI * rv);
    //cout<< "Fix Residual Variance = " << rv << endl;
    //cout << "tau = " << tau << "; log(2pi*rv) = " <<logrv << endl;
    
    inv_subvar.assign(n_type, 0.0), log_subvar.assign(n_type, 0.0);
    for(size_t i=0; i < n_type; i++){
        inv_subvar[i] = (1.0 / subvar[i]); 
        log_subvar[i] = (log(subvar[i])); 
    }
    //cout << "inv_subvar = "; PrintVector(inv_subvar);
    //cout << "log_subvar = "; PrintVector(log_subvar);
    
    if (cHyp_old.n_gamma > 0) {
        SetXgamma (Xgamma_old, X, rank_old);
        CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
    }
    
    //cout << "Set m_gamma... \n";
    set_mgamma(cHyp_old, rank_old, snp_pos);
    cout << "Initial number of selected variants per category : ";
    PrintVector(cHyp_old.m_gamma); 
    //cout << "Set sigma_subvec... \n";
    getSubVec(sigma_subvec_old, rank_old, snp_pos);
    //PrintVector(sigma_subvec_old, rank_old.size());
    
    cHyp_initial=cHyp_old;
    gsl_vector_memcpy(sigma_subvec_new, sigma_subvec_old);
    
    //Calculate first loglikelihood
    //cout << "first calculating logpost ... \n";
    if (cHyp_old.n_gamma==0) {
        loglikegamma = CalcLikegamma(cHyp_old);
        logPost_old = CalcPosterior (ztz, cHyp_old) + loglikegamma;
        loglike_old = loglikegamma;
    }
    else {
        loglikegamma = CalcLikegamma(cHyp_old);
        logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag, loglike_old) + loglikegamma;
        loglike_old += loglikegamma;
    }
    if (!Error_Flag) {
       // cout <<  "logPost_old = " << logPost_old << endl;
    }
    else {
        cerr << "Failed at initialMCMC...\n";
        exit(-1);
    }
    //cout <<  "Initial logPost_old = " << logPost_old << endl;
    
    //calculate centered z_hat, and pve
    if (a_mode==13) {
        if (cHyp_old.n_gamma==0) {
            CalcCC_PVEnZ (z_hat, cHyp_old);
        }
        else {
            CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
        }
    }
    
    //Start MCMC
    size_t k_save_sample=0;
    w_pace=1000;
    int accept; // accept_theta; naccept_theta=0,
    size_t total_step=w_step+s_step;
    size_t repeat=1;
    int flag_gamma=0;
    double accept_percent, betai; // accept_theta_percent;

    
    cHyp_new = cHyp_old;
    rank_new = rank_old;

    vector <string> snps_mcmc; // save locations of included snps per iteration
    string snps_mcmc_temp;
    size_t order_i;

    for (size_t t=0; t<total_step; ++t) {
        
       if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));
                cout << endl;
            }
        //		if (t>10) {break;}
        
        if (a_mode==13) {
            SampleZ (y, z_hat, z); //sample z
            mean_z=CenterVector (z);
            gsl_blas_ddot(z,z,&ztz);
            
            //First proposal need to be revised
            if (cHyp_old.n_gamma==0) {
                loglikegamma = CalcLikegamma(cHyp_old);
                logPost_old = CalcPosterior (ztz, cHyp_old) + loglikegamma;
                loglike_old = loglikegamma;
            } else {
                gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
                gsl_vector_view Xtz_sub=gsl_vector_subvector(Xtz_old, 0, rank_old.size());
                gsl_blas_dgemv (CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0, &Xtz_sub.vector);
                loglikegamma = CalcLikegamma(cHyp_old);
                logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag, loglike_old) + loglikegamma;
                loglike_old += loglikegamma;
            }
        }

        //////// Set repeat number
        //if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
        //else {repeat=1;}
        //cout << "n_mh = " << n_mh << endl;
        
        for (size_t i=0; i<n_mh; ++i) {

            //cout << "\n \n propose gamam...\n";
            //cout << "old rank: "; PrintVector(rank_old);
            //repeat = 1;
            logMHratio = ProposeGamma (rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat, X, z, Xgamma_old, XtX_old, Xtz_old, ztz, flag_gamma, Xgamma_new, XtX_new, Xtz_new); //JY
           // rank_new.clear(); cHyp_new.n_gamma=0;
            //cout << "propose new rank: "; PrintVector(rank_new);
            //cout << "flag_gamma = " << flag_gamma << endl;
            //cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
            //cout << "propose gamma logMHratio = "<<logMHratio << "; MHratio = " << exp(logMHratio) << endl;
            
        if (flag_gamma > 0) {

            if(flag_gamma==1) nadd++;
            else if(flag_gamma==2) ndel++;
            else  nswitch++;            
            
            if (rank_new.size() > 0) {
                set_mgamma(cHyp_new, rank_new, snp_pos);
                getSubVec(sigma_subvec_new, rank_new, snp_pos);
                loglikegamma = CalcLikegamma(cHyp_new);
                logPost_new = CalcPosterior (Xgamma_new, XtX_new, Xtz_new, ztz, Xb_new, beta_new, cHyp_new, sigma_subvec_new, Error_Flag, loglike_new) + loglikegamma;
                loglike_new += loglikegamma;
            }
            else {
                cHyp_new.m_gamma.assign(n_type, 0);
                loglikegamma = CalcLikegamma(cHyp_new);
                logPost_new = CalcPosterior (ztz, cHyp_new) + loglikegamma;
                loglike_new = loglikegamma;
            }
           // cout << "new m_gamma: " << cHyp_new.m_gamma[0] << ", "<< cHyp_new.m_gamma[1]<< endl;

             // cout << "Calcposterior success." << endl;
            if (!Error_Flag) {
                logMHratio += logPost_new-logPost_old;
                //cout <<"logPost_old = " << logPost_old<< "; logPost_new = "<< logPost_new<< "\n logMHratio = " << logMHratio<< "; MHratio = " << exp(logMHratio) << endl;
                if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio)
                    { accept=1; if (flag_gamma < 4) n_accept++;}
                else {accept=0;}
            }
            else{
                accept=0;
            }
        }
        else{
            nother++;
            accept = 0;
        }
            
            //cout << "accept = " << accept << endl;
            
            if (accept==1) {
                    if(flag_gamma==1) nadd_accept++;
                    else if(flag_gamma==2) ndel_accept++;
                    else if(flag_gamma==3) nswitch_accept++;
                    else nother_accept++;
                
                    logPost_old=logPost_new;
                    loglike_old = loglike_new;
                    cHyp_old.n_gamma = cHyp_new.n_gamma;
               // cout << "cHyp_old.m_gamma = "; PrintVector(cHyp_old.m_gamma);
                    cHyp_old.m_gamma = cHyp_new.m_gamma;
                //cout << "cHyp_new.m_gamma = "; PrintVector(cHyp_new.m_gamma);
                //cout << "cHyp_old.m_gamma = "; PrintVector(cHyp_old.m_gamma);
                    cHyp_old.pve = cHyp_new.pve;
                    //cHyp_old.rv = cHyp_new.rv;
                    //cHyp_old.pge = cHyp_new.pge;
                    gsl_vector_memcpy (Xb_old, Xb_new);
                rank_old.clear();
                for (size_t i=0; i<rank_new.size(); i++) {
                    rank_old.push_back(rank_new[i]);
                }
                if (rank_old.size() != rank_new.size()) {
                    cerr << "Error: rank_old size != rank_new size\n";
                    exit(-1);
                }
                //cout << "Accept proposal: "; PrintVector(rank_old);
                
                if(rank_old.size()>0){
                    gsl_vector_view sigma_oldsub=gsl_vector_subvector(sigma_subvec_old, 0, rank_old.size());
                    gsl_vector_view sigma_newsub=gsl_vector_subvector(sigma_subvec_new, 0, rank_old.size());
                    gsl_vector_memcpy(&sigma_oldsub.vector, &sigma_newsub.vector);
                
                    gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXold_sub=gsl_matrix_submatrix(XtX_old, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtzold_sub=gsl_vector_subvector(Xtz_old, 0, rank_new.size());
                    gsl_vector_view betaold_sub=gsl_vector_subvector(beta_old, 0, rank_new.size());
                    
                    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(Xgamma_new, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtznew_sub=gsl_vector_subvector(Xtz_new, 0, rank_new.size());
                    gsl_vector_view betanew_sub=gsl_vector_subvector(beta_new, 0, rank_new.size());
                    
                    gsl_matrix_memcpy(&Xold_sub.matrix, &Xnew_sub.matrix);
                    gsl_matrix_memcpy(&XtXold_sub.matrix, &XtXnew_sub.matrix);
                    gsl_vector_memcpy(&Xtzold_sub.vector, &Xtznew_sub.vector);
                    gsl_vector_memcpy(&betaold_sub.vector, &betanew_sub.vector);
                    
                    gsl_vector_memcpy(Xb_old, Xb_new);
                }
                //else{
                  //  gsl_vector_set_zero(Xb_old); // set Xb = 0
                //}
            } else {
                cHyp_new.n_gamma = cHyp_old.n_gamma;
                cHyp_new.m_gamma = cHyp_old.m_gamma;
                rank_new.clear();
            }
          //  cout << "copy data from new propose -> old " << endl;
        } //end of n_mh
        
        //calculate z_hat, and pve
        if (a_mode==13) {
            if (cHyp_old.n_gamma==0) {
                CalcCC_PVEnZ (z_hat, cHyp_old);
            }
            else {
                CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
            }
            //sample mu and update z hat
            gsl_vector_sub (z, z_hat);
            mean_z+=CenterVector(z);
            mean_z+=gsl_ran_gaussian(gsl_r, sqrt(1.0/(double) ni_test) );
            gsl_vector_add_constant (z_hat, mean_z);
        }
        
         //if (t % 10 == 0 && t > w_step) {
         if (t % w_pace == 0 && t > w_step) {
             accept_percent = (double)n_accept/(double)((t+1) * n_mh);
             //cout << "cHyp_old.n_gamma= " << cHyp_old.n_gamma << endl;
             cout << "acceptance percentage = " << setprecision(6) << accept_percent << endl ;
             cout << "# of selected variants per category: " << endl; PrintVector(cHyp_old.m_gamma);
             //cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); cout << endl;
             cout << "loglike: " << loglike_old << endl;
        }
        
        //Save data
        if (t<w_step) {continue;}
        else {
            //save loglikelihood
            gsl_vector_set (LnPost, k_save_sample, loglike_old);
            GV += cHyp_old.pve;
            
            if (cHyp_old.n_gamma > 0){

                region_pip++; //count increase if the model has >0 SNPs
                snps_mcmc_temp="";

                for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
                    // beta_g saved by position
                    pos=SNPrank_vec[rank_old[i]].first;
                    order_i = SNPrank_vec[rank_old[i]].second;

                    betai = gsl_vector_get(beta_old, i);
                    beta_g[pos].first += betai;
                    beta_g[pos].second += 1.0;
                    for (size_t j=0; j < n_type; j++) {
                        if (snp_pos[order_i].indicator_func[j]) {
                            sumbeta2[j] += betai * betai;
                            break;
                        }
                    }
                    snps_mcmc_temp += string(snp_pos[order_i].rs) + string(":") + string(snp_pos[order_i].chr) + string(":") + to_string(snp_pos[order_i].bp) + string(":") + to_string(snp_pos[order_i].a_major) + string(":") + to_string(snp_pos[order_i].a_minor)+ string(";");
                }
                snps_mcmc.push_back(snps_mcmc_temp);

                //if(cHyp_old.m_gamma[0]>0)
                   // sample_sigma0.push_back(make_pair(sumbeta2[0] /(double)cHyp_old.m_gamma[0], cHyp_old.m_gamma[0] ));
                //if(cHyp_old.m_gamma[1]>0)
                   // sample_sigma1.push_back(make_pair(sumbeta2[1] /(double)cHyp_old.m_gamma[1], cHyp_old.m_gamma[1] ));

                for(size_t j=0; j < n_type; j++){
                    if(cHyp_old.m_gamma[j] > 0) 
                        em_gamma[j] += (double)cHyp_old.m_gamma[j];
                    }
                k_save_sample++;
              }
            }
    }
    
    cout<< "MCMC completed ... " << endl << endl;
    cout << "region_pip = " << region_pip << endl;

    accept_percent = (double)n_accept/(double)(total_step * n_mh);
    cout << "gamma acceptance percentage = " << accept_percent << endl ;
    cout << "# of selected variants per category: "; PrintVector(cHyp_old.m_gamma);
    cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); 
    cout << "loglike: " << loglike_old << endl;
    cout << "k_save_sample = " << k_save_sample << endl;

    
    //save all marker information
    if (saveSNP) {
        //cout << "Write genotype txt file after snp_pos is ordered ... \n";
        WriteGenotypeFile(X, snp_pos);
        exit(-1);

        //WriteIniSNP(rank_old, snp_pos);
        //WriteParam(beta_g, snp_pos, pos_loglr, Z_scores, SE_beta, pval_lrt);
        //WriteFGWAS_InputFile(snp_pos, Z_scores, SE_beta);
    }

    //Save temp EM results
    WriteHyptemp(LnPost, em_gamma);
    WriteParam(beta_g, snp_pos, pos_loglr, Z_scores, SE_beta, pval_lrt);
    WriteMCMC(snps_mcmc); // save all active SNPs from MCMC
    
   // gsl_matrix_free(Result_hyp);
   // gsl_matrix_free(Result_gamma);
    // gsl_matrix_free(Sample_m);
    
    gsl_vector_free(sigma_subvec_old);
    gsl_vector_free(sigma_subvec_new);
    gsl_vector_free(LnPost);
    
    gsl_vector_free(z_hat);
    gsl_vector_free(z);
    gsl_vector_free(Xb_new);	
    gsl_vector_free(Xb_old);
    
    //gsl_vector_free(Xb_mcmc);
    
    gsl_matrix_free(Xgamma_old);
    gsl_matrix_free(XtX_old);
    gsl_vector_free(Xtz_old);
    gsl_vector_free(beta_old);
    
    gsl_matrix_free(Xgamma_new);
    gsl_matrix_free(XtX_new);
    gsl_vector_free(Xtz_new);
    gsl_vector_free(beta_new);
    
    delete [] p_gamma;
    beta_g.clear();
    
    return;
}
// end of current version

//calculat LD correlation matrix
void BVSRM::CalcLD(uchar **X, gsl_matrix * LD){

    gsl_vector *xvec_i = gsl_vector_alloc(ni_test);
    gsl_vector *xvec_j = gsl_vector_alloc(ni_test);
    double xtx_ij, cor_ij;

    cout << "set SNPsd as : \n"; 
    if(SNPsd.size() != ns_test){
        SNPsd.clear();
        for (size_t i=0; i<ns_test; ++i) {
            SNPsd.push_back(sqrt(XtX_diagvec[i]));
        }
    }
    PrintVector(SNPsd, 20);

    for (size_t i=0; i<ns_test; ++i) {

        gsl_matrix_set(LD, i, i, 1);
        getGTgslVec(X, xvec_i, i, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);

        for(size_t j=(i+1); j < ns_test; ++j){
            getGTgslVec(X, xvec_j, j, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
            gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
            cor_ij = xtx_ij / (SNPsd[i] * SNPsd[j]);
            gsl_matrix_set(LD, i, j, cor_ij);
            gsl_matrix_set(LD, j, i, cor_ij);
        }
    }

    gsl_vector_free(xvec_i);
    gsl_vector_free(xvec_j);
    SNPsd.clear();

    return;
}




void BVSRM::SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z)
{
	double d1, d2, z_rand=0.0;
	for (size_t i=0; i<z->size; ++i) {
		d1=gsl_vector_get (y, i);
		d2=gsl_vector_get (z_hat, i);
		//y is centerred for case control studies
		if (d1<=0.0) {
			//control, right truncated
			do {
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand>0.0);
		}
		else {
			do {
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand<0.0);
		}
		
		gsl_vector_set (z, i, z_rand);
	}

	return;
}







//JY edit start
void BVSRM::CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX, const gsl_vector *Xtz, gsl_vector *z_res, const size_t &s_size, const double &ztz){
    
    gsl_matrix_const_view X_gsub=gsl_matrix_const_submatrix(Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_gsub = gsl_matrix_const_submatrix(XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xtz_gsub = gsl_vector_const_subvector(Xtz, 0, s_size);
    
    gsl_vector *beta_gamma_hat = gsl_vector_alloc(s_size);
    gsl_matrix *XtXtemp = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(XtXtemp, &XtX_gsub.matrix);
    
    double SSR, R2 ;
    
    if(LapackSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat) != 0)
        EigenSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat);
    gsl_blas_ddot(&Xtz_gsub.vector, beta_gamma_hat, &SSR);
    R2 = (SSR / ztz);

    /* double lambda = 0.0;
    int k=0;
    
   if ((R2 < -0.0) || (R2 > 1.0)) {
        cout << "R2 in Calcres = " << R2 << endl;
        PrintMatrix(XtXtemp, s_size, s_size);
        WriteMatrix(&X_gsub.matrix, "_Xres");
        WriteMatrix(XtXtemp, "_XtXres");
        WriteVector(beta_gamma_hat, "_bres");
        WriteVector(z_res, "_zres");
        //exit(-1);
    }*/
    
if ((R2 < -0.1) || (R2 > 1.1)) {
    /*
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(XtX, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.01;
    // cout << "labmda = " << lambda << endl;
    gsl_vector_view XtXtemp_diag = gsl_matrix_diagonal(XtXtemp);
    
    while ((R2 < -0.1) || (R2 > 1.1)){
        
        cout << "add lambda: negative R2 in calcresidual ...  "<< R2 << "; k = " << k << endl;
        
        gsl_vector_add_constant(&XtXtemp_diag.vector, lambda);
        if(LapackSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat) != 0)
            EigenSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat);
        gsl_blas_ddot(&Xtz_gsub.vector, beta_gamma_hat, &SSR);
        R2 = (SSR / ztz);
        k++;
        cout << "now R2 = " << R2 << "; k = " << k << endl;
        if (k > 9) {
            break;
        }
    }
    
    if ((R2 < 0.0) || (R2 > 1.0)) {
        gsl_vector_memcpy(z_res, z);
    }*/
    gsl_vector_memcpy(z_res, z);
    
}
else if ( (R2 < 0.0) || (R2 > 1.0) ){
    gsl_vector_memcpy(z_res, z);
}
else{
    gsl_blas_dgemv(CblasNoTrans, 1.0, &X_gsub.matrix, beta_gamma_hat, 0.0, z_res);
    gsl_vector_scale(z_res, -1.0);
    gsl_vector_add(z_res, z);
}
    
    gsl_matrix_free(XtXtemp);
    gsl_vector_free(beta_gamma_hat);
    return;
}





//calculate likelihood ratio statistic
double BVSRM::CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec, size_t posj){
    double LR;
    double xtz_res, ztz_res, xtx_vec = XtX_diagvec[posj];
    
    gsl_blas_ddot(z_res, z_res, &ztz_res);
    //gsl_blas_ddot(x_vec, x_vec, &xtx_vec);
    gsl_blas_ddot(x_vec, z_res, &xtz_res);
    //cout << "ztz_res = " << ztz_res << "; xtx_vec = " << xtx_vec << "; xtz_res = " << xtz_res << endl;
    //double ixtx = 1.0 / xtx_vec;
    //double bhat = ixtx * xtz_res;
    //double V = (ixtx * ztz_res - bhat * bhat) / (ni_test);
    //double Z2 = bhat * bhat / V;
    //double VpW = V + Wvar;
    LR = (ni_test)*(log(ztz_res)-log(ztz_res-xtz_res*xtz_res/xtx_vec));
    //sqrt(VpW / V) * exp(-0.5 * Z2 * Wvar / VpW);
    //cout << "log LR = " << BF << ", ";
    return (LR);
}

gsl_ran_discrete_t * BVSRM::MakeProposal(const size_t &o, double *p_BF, uchar **X, const gsl_vector *z_res, const map<size_t, int> &mapRank2in)
{
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
    
    long int orderj;
    size_t posj, rank_j;
    vector<int> j_ind;
    double pmax, psum=0.0, countj = 0.0;
    
    for (size_t j=0; j < ns_neib; ++j){
        orderj = (o - win) + j;
        if((orderj >= 0) && (orderj < (long int)ns_test))
            rank_j = SNPorder_vec[(size_t)orderj].second;
        if((orderj >= 0) && (orderj < (long int)ns_test) && (j != win) && (mapRank2in.count(rank_j) == 0)){
            posj = SNPorder_vec[(size_t)orderj].first;
            getGTgslVec(X, xvec, posj, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
            p_BF[j]=CalcLR(z_res, xvec, posj); //calc loglr
            j_ind.push_back(1);
            countj += 1.0;
        }
        else{
            p_BF[j] = -std::numeric_limits<double>::max();
            j_ind.push_back(0);
        }
    }
    pmax = *std::max_element(p_BF, p_BF+(ns_neib));
    
    for (size_t j=0; j < ns_neib; ++j){
        if(j_ind[j]==1){
            p_BF[j]=exp(p_BF[j]- pmax);
            psum += p_BF[j];
        }
        else{p_BF[j] = 0.0;}
    }
    //
    psum = 1.0/psum;
    for(size_t j=0; j < ns_neib; ++j){
        p_BF[j] *= psum;
    }
    
    gsl_vector_free(xvec);
    
    return (gsl_ran_discrete_preproc(ns_neib, p_BF));
}


bool BVSRM::ColinearTest(uchar ** X, const gsl_matrix * Xtemp, const gsl_matrix * XtX_temp, size_t r_add, size_t s_size)
{
    bool colinear = 0;
    double vreg;
    size_t pos = SNPrank_vec[r_add].first;
    double xtx = XtX_diagvec[pos];
    
    
    gsl_vector *beta_temp = gsl_vector_alloc(s_size);
    gsl_vector *Xtx_temp = gsl_vector_alloc(s_size);
    gsl_vector *xvec_temp = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec_temp, pos, ni_test, ns_test, SNPsd, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
    //gsl_blas_ddot(xvec_temp, xvec_temp, &xtx);
    //cout << "XtX_diagvec[pos] = "<< XtX_diagvec[pos] << "; xtx = " << xtx << endl;
    
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xtemp, 0, 0, Xtemp->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX_temp, 0, 0, s_size, s_size);
    //cout << "XtX in colinear test: \n "; PrintMatrix(&XtX_sub.matrix, s_size, s_size);
    gsl_blas_dgemv(CblasTrans, 1.0, &Xgamma_sub.matrix, xvec_temp, 0.0, Xtx_temp);

    if (LapackSolve(&XtX_sub.matrix, Xtx_temp, beta_temp) !=0)
        EigenSolve(&XtX_sub.matrix, Xtx_temp, beta_temp);
    
    gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
    //cout << "vreg = " << vreg << endl;
    
    double tR2 = 0.95;
    double R2 = (vreg / xtx);
//    int k=0;
//    double lambda = 0.0;
//    gsl_matrix *XtXlu = gsl_matrix_alloc(s_size, s_size);
    
    if ( (R2 >= tR2) && (R2 <= 1.1) ) {
        colinear = 1;
       // cout << "R2 in ColinearTest = " << R2 << endl;
    }
    else if ((R2 < -0.0) || (R2 > 1.0)){
       colinear = 1;
        //cout << "R2 in ColinearTest = " << R2 << "; k = " << k << endl;
        /*
        //PrintMatrix(&XtX_sub.matrix, s_size, s_size);
        WriteMatrix(&Xgamma_sub.matrix, "_Xct");
        WriteMatrix(&XtX_sub.matrix, "_XtXct");
        WriteVector(beta_temp, "_bct");
        WriteVector(xvec_temp, "_xvec");
        
        gsl_matrix_memcpy(XtXlu, &XtX_sub.matrix);
        for (size_t i=0; i<s_size; ++i) {
            lambda += gsl_matrix_get(XtXlu, i, i);
        }
        lambda /= (double)s_size;
        lambda *= 0.01;
        gsl_vector_view XtXlu_diag = gsl_matrix_diagonal(XtXlu);
    
        while ((R2 < -0.1) || (R2 > 1.1)) {
        
            gsl_vector_add_constant(&XtXlu_diag.vector, lambda);
            if (LapackSolve(XtXlu, Xtx_temp, beta_temp) !=0)
                    EigenSolve(XtXlu, Xtx_temp, beta_temp);
            gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
            R2 = (vreg / xtx);
            k++;
            cout << "now R2 = " << R2 << "; k = " << k << endl;
        
            if (k>9) {
                cout << "reached k = " << k << endl;
                WriteMatrix(&XtX_sub.matrix, "_XtX_ct");
                WriteVector(Xtx_temp, "_Xtx_ct");
                WriteVector(beta_temp, "_bct");
                break;
            }
        }
        
        if ((R2 >= tR2) && (R2 <= 1.1)) {
            colinear = 1;
        } else {colinear = 0;} */
    }
    else {
        colinear = 0;
    }
    
 //   gsl_matrix_free(XtXlu);
    gsl_vector_free(xvec_temp);
    gsl_vector_free(beta_temp);
    gsl_vector_free(Xtx_temp);
    
    return colinear;
}

//below fits MCMC for rho=1
void BVSRM::CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty)
{
  time_t time_start=clock();	
  gsl_matrix_const_view X_sub=gsl_matrix_const_submatrix(X, 0, 0, X->size1, s_size);
  gsl_matrix_view XtX_sub=gsl_matrix_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_view Xty_sub=gsl_vector_subvector(Xty, 0, s_size);

  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &X_sub.matrix, &X_sub.matrix, 0.0, &XtX_sub.matrix);
  gsl_blas_dgemv(CblasTrans, 1.0, &X_sub.matrix, y, 0.0, &Xty_sub.vector);

  time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

  return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BVSRM::CalcCC_PVEnZ (gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
  gsl_vector_set_zero(z_hat);
  cHyp.pve=0.0;
  //cHyp.pge=1.0;
  return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BVSRM::CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	
	gsl_blas_ddot (Xb, Xb, &d);
	cHyp.pve=d/(double)ni_test;
	//cHyp.pve/=cHyp.pve+1.0;
	//cHyp.pge=1.0;
	
	gsl_vector_memcpy (z_hat, Xb);

	return;
}


void BVSRM::CreateSnpPosVec(vector<SNPPOS> &snp_pos)
{
    size_t pos;
    string rs;
    string chr;
    long int bp;
    size_t tt=0;
    vector<bool> indicator_func;
    //vector<double> weight;
    //double weight_i;
    double maf;
    string a_minor;
    string a_major;

    //SNPsd.clear();
    
    for (size_t i=0; i < ns_total; ++i){
        if(!indicator_snp[i]) {continue;}
        
        pos=tt;
       // if(tt == pos_loglr[tt].first ) pos = tt;
        //else cout << "error assigning position to snp_pos vector"<< endl;
        
        rs = snpInfo[i].rs_number;
        chr = snpInfo[i].chr;
        bp = snpInfo[i].base_position;
        maf = snpInfo[i].maf;
        a_minor = snpInfo[i].a_minor;
        a_major = snpInfo[i].a_major;

        indicator_func = snpInfo[i].indicator_func;

        //weight = snpInfo[i].weight;
        //weight_i = snpInfo[i].weight_i;
        // cout << weight_i << ", ";
        SNPPOS snp_temp={pos, rs, chr, bp, a_minor, a_major, maf, indicator_func};
        snp_pos.push_back(snp_temp);
        
        //SNPsd.push_back(1.0 / sqrt(2.0 * maf * (1.0 - maf)));
        
        tt++;
    }
    snpInfo.clear();
    
}







