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
#include <stdlib.h> 
#include <bitset>
#include <cstring>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"


#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_integration.h"

#include "gzstream.h"
#include "lapack.h"
#include "lm.h"
#include "param.h"



using namespace std;





void LM::CopyFromParam (PARAM &cPar) 
{
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	file_vcf=cPar.file_vcf;
	
	time_opt=0.0;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	indicator_idv=cPar.indicator_idv;	
	indicator_snp=cPar.indicator_snp;	
	snpInfo=cPar.snpInfo;
	
	return;
}


void LM::CopyToParam (PARAM &cPar) 
{
	cPar.time_opt=time_opt;	
	cPar.pheno_var = pheno_var;	
	return;
}



void LM::WriteFiles () 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".assoc.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
		/*outfile<<"CHR"<<"\t"<<"ID"<<"\t"<<"POS"<<"\t"<<"n_miss"<<"\t"<<"ALT"<<"\t"<<"REF"<<"\t"<<"maf"<<"\t";
		
		if (a_mode==51) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_wald"<<endl;
		} else if (a_mode==52) {
			outfile<<"p_lrt"<<endl;
		} else if (a_mode==53) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_score"<<"\t"<<"score_u"<<"\t"<<"score_v"<<endl;
		} else if (a_mode==54) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_wald"<<"\t"<<"p_lrt"<<"\t"<<"p_score"<<"\t"<<"score_u"<<"\t"<<"score_v"<<endl;
		} else {}
		*/
		
		size_t t=0;
		for (size_t i=0; i<snpInfo.size(); ++i) {
			if (indicator_snp[i]==0) {continue;}
			
			outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t"<<snpInfo[i].a_minor<<"\t"<<snpInfo[i].a_major<<"\t"<<fixed<<setprecision(3)<<snpInfo[i].maf<<"\t";
			
			if (a_mode==51) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_wald <<endl;
			} else if (a_mode==52) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].p_lrt<<endl;
			} else if (a_mode==53) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_score<<"\t"<<sumStat[t].score_u<<"\t"<<sumStat[t].score_v<<endl;
			} else if (a_mode==54) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_wald <<"\t"<<sumStat[t].p_lrt<<"\t"<<sumStat[t].p_score<<"\t"<<sumStat[t].score_u<<"\t"<<sumStat[t].score_v<<endl;
			} else {}
			t++;
		}
	
		
	outfile.close();
	outfile.clear();
	return;
}



void CalcvPv(const gsl_matrix *WtWi, const gsl_vector *Wty, const gsl_vector *Wtx, const gsl_vector *y, const gsl_vector *x,  double &xPwy, double &xPwx)
{
	size_t c_size=Wty->size;
	double d;
	
	gsl_vector *WtWiWtx=gsl_vector_alloc (c_size);
	
	gsl_blas_ddot (x, x, &xPwx);
	gsl_blas_ddot (x, y, &xPwy);
	gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);	
	
	gsl_blas_ddot (WtWiWtx, Wtx, &d);	
	xPwx-=d;
	
	gsl_blas_ddot (WtWiWtx, Wty, &d);	
	xPwy-=d;
	
	gsl_vector_free (WtWiWtx);
	
	return;
}


void CalcvPv(const gsl_matrix *WtWi, const gsl_vector *Wty, const gsl_vector *y, double &yPwy)
{
	size_t c_size=Wty->size;
	double d;
	
	gsl_vector *WtWiWty=gsl_vector_alloc (c_size);
	
	gsl_blas_ddot (y, y, &yPwy);
	gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wty, 0.0, WtWiWty);	
	
	gsl_blas_ddot (WtWiWty, Wty, &d);	
	yPwy-=d;
	
	gsl_vector_free (WtWiWty);
	
	return;
}


//calculate p values and beta/se in a linear model
void LmCalcP (const size_t test_mode, const double yPwy, const double xPwy, const double xPwx, const double df, const size_t n_size, double &beta, double &se, double &p_wald, double &p_lrt, double &p_score)
{
	double yPxy=yPwy-xPwy*xPwy/xPwx;
	double se_wald, se_score;
	
	beta=xPwy/xPwx;
	se_wald=sqrt(yPxy/(df*xPwx) );
	se_score=sqrt(yPwy/((double)n_size*xPwx) );
	
	p_wald=gsl_cdf_fdist_Q (beta*beta/(se_wald*se_wald), 1.0, df);
	p_score=gsl_cdf_fdist_Q (beta*beta/(se_score*se_score), 1.0, df);
	p_lrt=gsl_cdf_chisq_Q ((double)n_size*(log(yPwy)-log(yPxy)), 1);
	
	if (test_mode==3) {se=se_score;} else {se=se_wald;}
	
	return;
}


//Begin of AnalyzeVCF 
void LM::AnalyzeVCF (const gsl_matrix *W, const gsl_vector *y, string &GTfield, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID)
{
	if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string

    igzstream infile (file_vcf.c_str(), igzstream::in);
    if (!infile) {cout<<"error reading vcf genotype file:"<<file_vcf<<endl; exit(-1);}
	
	clock_t time_start=clock();
	
	string line, pheno_id;
	char *pch, *p, *nch=NULL, *n;
	
	double beta=0, se=0, p_wald=0, p_lrt=0, p_score=0, score_u=0, score_v=0;
	size_t n_miss, c_idv=0, c_snp=0, ctest_idv = 0, tab_count, pheno_index, ns_test=0;
	double geno, x_mean;
	int GTpos=0, k=0;
	
	//calculate some basic quantities
	double yPwy, xPwy, xPwx;
	double df=(double)ni_test-(double)W->size2-1.0;

	gsl_vector *x=gsl_vector_alloc (ni_test);
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);		
	gsl_vector *Wty=gsl_vector_alloc (W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);

	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);

	gsl_blas_dgemv (CblasTrans, 1.0, W, y, 0.0, Wty);
	CalcvPv(WtWi, Wty, y, yPwy);
	pheno_var = yPwy / ( (double)ni_test - 1.0 )  ;
	
	//start reading genotypes and analyze	
	while(!safeGetline(infile, line).eof())
    {
        if (c_snp%d_pace==0 || c_snp==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", c_snp, indicator_snp.size()-1);}

        if (line[0] == '#') {
            continue; //skip header
        }
        else {
            if (!indicator_snp[c_snp]) {c_snp++; continue;}
            c_idv=0; //increase to the total individuals ni_total
            x_mean=0.0; n_miss=0; ctest_idv = 0;
        
            pch= (char *)line.c_str();
            for (tab_count=0; pch != NULL; tab_count++) {
                nch=strchr(pch, '\t'); //point to the position of next '\t'

                if ((tab_count == 8) && (c_idv == 0))
                {
                    // parse FORMAT field
                    if (pch[0] == GTfield[0] && pch[1] == GTfield[1] && ((nch==pch+2)||pch[2]==':') ) {
                        GTpos=0; //GT start in the first position
                    }
                    else if (nch == NULL){ cerr << "VCF has FORMAT field but dose not have any genotype\n";}
                    else{
                        k=0; //index of key characters
                        GTpos=0;
                        p=pch;
                        while (p<nch) {
                            if (*p == ':') {
                                if (k >= lkey) {
                                    break;
                                }
                                else {
                                    ++GTpos;
                                    k=0;
                                }
                            }
                            else {
                                if (GTfield[k] == *p) {
                                    ++k;
                                }
                                else { k=0; }
                            }
                            ++p;
                        }
                        if ((p==nch) && (k != lkey)) {
                            cerr << "Cannot find" << GTfield << endl;
                            exit(-1);
                        }
                    }
                }
                else if ( tab_count == SampleVcfPos[ctest_idv] )
                {                   
                    pheno_id = VcfSampleID[c_idv];
                    if (PhenoID2Ind.count(pheno_id) > 0){
                            pheno_index = PhenoID2Ind.at(pheno_id); 
                    }
                    else {
                        cerr << "error: pheno ID matched error ... "<< endl;
                        exit(-1);
                    }
                    if ( !indicator_idv[pheno_index] ) {
                        cerr << "error: pheno is not in sample ... "<< endl;
                        exit(-1);
                        //continue;
                    }
                    else{
                        p = pch; // make p reach to the key index
                        if (GTpos>0) {
                            for (int i=0; (i<GTpos) && (p!=NULL); ++i) {
                                n = strchr(p, ':');
                                p = (n == NULL) ? NULL : n+1;
                            }
                        }
                        if (p==NULL) {
                            geno = -9;//missing
                            gsl_vector_set (x, ctest_idv, geno);
                            n_miss++; c_idv++; ctest_idv++; 
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;
                        }
                        else if ( (p[1] == '/') || (p[1] == '|') ) {
                        //read bi-allelic GT
                            if( (p[0]=='.') && (p[2]=='.')){
                                geno = -9;//missing
                                gsl_vector_set (x, ctest_idv, geno);
                                n_miss++; c_idv++; ctest_idv++;
                                pch = (nch == NULL) ? NULL : nch+1;
                                continue;
                            }
                            else if ( (p[0]=='.') && (p[2]!='.')) {
                                geno = (double)(p[2] -'0');
                            }
                            else if ((p[0]!='.') && p[2]=='.') {
                                geno = (double)(p[0] -'0');
                            }
                            else geno = (double)((p[0] - '0') + (p[2]- '0'));
                        }
                        else {
                            //read dosage data
                            if( (p[0]=='.') && ( (p[1] == '\t') || (p[1] == ':') ) ){
                                geno = -9;
                                gsl_vector_set (x, ctest_idv, geno);
                                n_miss++; c_idv++; ctest_idv++;
                                pch = (nch == NULL) ? NULL : nch+1;
                                continue;                               
                            }else if (isdigit(p[0])){
                                geno = strtod(p, NULL);
                            }else{
                                cerr << "dosage data is not a digit ... " << endl;
                                exit(-1);
                            }                        
                        }
                        gsl_vector_set (x, ctest_idv, geno);
                        if( (geno >= 0.0) && (geno <= 2.0)) {x_mean += geno;}
                        ctest_idv++; // increase analyzed phenotype #
                        c_idv++;
                    }
                }
                else if ( tab_count >= 9 ) { c_idv++; }
                pch = (nch == NULL) ? NULL : nch+1;
            } 	
		ns_test++; c_snp++;

		x_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<ni_test; ++i) {
			if ( gsl_vector_get (x, i) == -9.0 ) 
            {
                gsl_vector_set(x, i, x_mean);
            }
			geno=gsl_vector_get(x, i);
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}		
		
		//calculate statistics		
		time_start=clock();		
		gsl_blas_dgemv(CblasTrans, 1.0, W, x, 0.0, Wtx);		
		CalcvPv(WtWi, Wty, Wtx, y, x, xPwy, xPwx);
		LmCalcP (a_mode-50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald, p_lrt, p_score);	

		// obtain score summary statistics
		if(a_mode == 53){
			score_u = xPwy;
			score_v = pheno_var * xPwx ;
		}
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, score_u, score_v};
		sumStat.push_back(SNPs);
		}
	}	
	cout<<endl;

	gsl_vector_free(x);
	gsl_matrix_free(WtW);
	gsl_matrix_free(WtWi);
	gsl_vector_free(Wty);
	gsl_vector_free(Wtx);
	gsl_permutation_free(pmt);
	
	infile.close();
	infile.clear();
	
	return;
}
// end of AnalyzeVCF


void LM::AnalyzeGeno (const gsl_matrix *W, const gsl_vector *y, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return;}
	
	clock_t time_start=clock();
	
	string line, pheno_id, s;
	char *pch, *nch=NULL;
	
	double beta=0, se=0, p_wald=0, p_lrt=0, p_score=0, score_u=0.0, score_v=0.0;
	int n_miss;
	double geno, x_mean;
	size_t c_idv, ctest_idv, tab_count, pheno_index, c_snp=0;
	
	//calculate some basic quantities
	double yPwy, xPwy, xPwx;
	double df=(double)W->size1-(double)W->size2-1.0;

	gsl_vector *x=gsl_vector_alloc (W->size1); // save vector of genotypes
	gsl_vector *x_miss=gsl_vector_alloc (W->size1);

	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);		
	gsl_vector *Wty=gsl_vector_alloc (W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);

	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);

	gsl_blas_dgemv (CblasTrans, 1.0, W, y, 0.0, Wty);
	CalcvPv(WtWi, Wty, y, yPwy);
	pheno_var = yPwy / ((double)ni_test - 1.0);

	while(!safeGetline(infile, line).eof()){

        if (c_snp%d_pace==0 || c_snp==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", c_snp, indicator_snp.size()-1);}

        pch= (char *)line.c_str();

        if ( (strncmp(line.c_str(), "ID", 2) == 0) ) {continue;} // skip header 
        else{
            if (indicator_snp[c_snp]==0) {c_snp++; continue;} // skip unanalyzed snp            
            c_idv=0; ctest_idv = 0; x_mean = 0.0; n_miss = 0; 
            for (tab_count=0; pch != NULL; tab_count++) {
                nch=strchr(pch, '\t'); //point to the position of next '\t'           
                if(tab_count == SampleVcfPos[ctest_idv] ) 
                {
                    pheno_id = VcfSampleID[c_idv];
                    pheno_index = PhenoID2Ind.at(pheno_id);
                  if ( !indicator_idv[pheno_index] ) {
                       cout << "phenotype of "<< pheno_id<<" is not analyzed."<< endl;
                       pch = (nch == NULL) ? NULL : nch+1;
                       c_idv++;  
                       continue;
                  } 
                  else{
                    // read genotype value
                    if (pch == NULL) {
                        //missing
                        gsl_vector_set (x, ctest_idv, -9.0);
                        n_miss++; c_idv++; ctest_idv++; 
                        pch = (nch == NULL) ? NULL : nch+1;
                        continue;
                    }
                    else {
                        //read dosage data 
                        if( ((pch[0]=='N') && (pch[1] == 'A')) || ((pch[0]=='.') && (pch[1] == '\t'))){
                            gsl_vector_set (x, ctest_idv, -9.0);                       
                            n_miss++; c_idv++; ctest_idv++;
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;                           
                        }else{
                            if (nch == NULL) { s.assign( pch );}
                            else s.assign( pch, nch-pch ); // field string s
                            geno = atof(s.c_str());
                        }  
                    }
                    if( (geno >= 0.0) && (geno <= 2.0)) {
                        gsl_vector_set (x, ctest_idv, geno);
                        x_mean += geno;
                    }else{
                        gsl_vector_set (x, ctest_idv, -9.0);
                        n_miss++; c_idv++; ctest_idv++;
                        pch = (nch == NULL) ? NULL : nch+1;
                        continue;
                    }
                    ctest_idv++;
                    c_idv++;
                  }
                }
                else if(tab_count >= 5){ c_idv++; }
                pch = (nch == NULL) ? NULL : nch+1;
            }
        }
		x_mean/=(double)(ni_test-n_miss);
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (x, i)==-9.0) {gsl_vector_set(x, i, x_mean);}
			geno=gsl_vector_get(x, i);
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}		
		
		//calculate statistics		
		time_start=clock();		
		gsl_blas_dgemv(CblasTrans, 1.0, W, x, 0.0, Wtx);		
		CalcvPv(WtWi, Wty, Wtx, y, x, xPwy, xPwx);
		LmCalcP (a_mode-50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald, p_lrt, p_score);
		if(a_mode == 53){
			score_u = xPwy;
			score_v = pheno_var * xPwx ;
		}
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, score_u, score_v};
		sumStat.push_back(SNPs);
	}	
	cout<<endl;

	gsl_vector_free(x);
	gsl_vector_free(x_miss);

	gsl_matrix_free(WtW);
	gsl_matrix_free(WtWi);
	gsl_vector_free(Wty);
	gsl_vector_free(Wtx);
	gsl_permutation_free(pmt);
	
	infile.close();
	infile.clear();
	
	return;
}


void LM::AnalyzePlink (const gsl_matrix *W, const gsl_vector *y) 
{
	string file_bed=file_bfile+".bed";
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return;}
	
	clock_t time_start=clock();
	
	char ch[1];
	bitset<8> b;	
	
	double beta=0, se=0, p_wald=0, p_lrt=0, p_score=0, score_u = 0.0, score_v = 0.0;
	int n_bit, n_miss, ci_total, ci_test;
	double geno, x_mean;
		
	//calculate some basic quantities
	double yPwy, xPwy, xPwx;
	double df=(double)W->size1-(double)W->size2-1.0;

	gsl_vector *x=gsl_vector_alloc (W->size1);

	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);	
	gsl_vector *Wty=gsl_vector_alloc (W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);

	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);

	gsl_blas_dgemv (CblasTrans, 1.0, W, y, 0.0, Wty);
	CalcvPv(WtWi, Wty, y, yPwy);
	pheno_var = yPwy / ((double)ni_test - 1.0);
		
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }
	
	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	
	for (vector<SNPINFO>::size_type t=0; t<snpInfo.size(); ++t) {
		if (t%d_pace==0 || t==snpInfo.size()-1) {ProgressBar ("Reading SNPs  ", t, snpInfo.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		x_mean=0.0;	n_miss=0; ci_total=0; ci_test=0; 
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==(int)ni_total) {break;}
				if (indicator_idv[ci_total]==0) {ci_total++; continue;}
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(x, ci_test, 2); x_mean+=2.0; }
					else {gsl_vector_set(x, ci_test, 1); x_mean+=1.0; }
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(x, ci_test, 0); }                                  
					else {gsl_vector_set(x, ci_test, -9); n_miss++; }
				}
				
				ci_total++;
				ci_test++;
			}
		}
		
		x_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<ni_test; ++i) {			
			geno=gsl_vector_get(x,i);
			if (geno==-9) {gsl_vector_set(x, i, x_mean); geno=x_mean;}
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		//calculate statistics		
		time_start=clock();	
		gsl_blas_dgemv (CblasTrans, 1.0, W, x, 0.0, Wtx);
		CalcvPv(WtWi, Wty, Wtx, y, x, xPwy, xPwx);		
		LmCalcP (a_mode-50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald, p_lrt, p_score);  

		if(a_mode == 53){
			score_u = xPwy;
			score_v = pheno_var * xPwx ;
		}
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, score_u, score_v};
		sumStat.push_back(SNPs);
	}	
	cout<<endl;
	
	gsl_vector_free(x);

	gsl_matrix_free(WtW);
	gsl_matrix_free(WtWi);	
	gsl_vector_free(Wty);
	gsl_vector_free(Wtx);
	gsl_permutation_free(pmt);
	
	infile.close();
	infile.clear();	
	
	return;
}



//make sure that both y and X are centered already
void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, double &trace_G, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag)
{
    trace_G=0.0;
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
	double yty, xty, xtx, log_lr;
	gsl_blas_ddot(y, y, &yty);

	for (size_t i=0; i<ns_test; ++i) {
        
      //getGTgslVec(X, xvec, i, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
      gsl_blas_ddot(xvec, xvec, &xtx);
	  gsl_blas_ddot(xvec, y, &xty);

	  log_lr=((double)y->size)*(log(yty)-log(yty-xty*xty/xtx));
	  pos_loglr.push_back(make_pair(i,log_lr) );
        trace_G += (xtx );
        //trace_G += (xtx / (double)ni_test);
	}
    //trace_G /= (double)ns_test;
	gsl_vector_free(xvec);
	return;
}


//
void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, const vector<double> &SNPmean, vector<double> &Gvec, vector<double> &XtX_diagvec, const vector<SNPPOS> &snp_pos, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag, const vector<pair<int, double> > &UcharTable)
{
    size_t n_type = snp_pos[0].indicator_func.size();
    Gvec.assign(n_type, 0.0);
    XtX_diagvec.clear();
    
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
	double yty, xty, xtx, log_lr;
	gsl_blas_ddot(y, y, &yty);
    cout << "calcLR: yty = " << yty <<"ns_test = " << ns_test << "; ni_test" << ni_test << endl;
    
	for (size_t i=0; i<ns_test; ++i) {
        
        getGTgslVec(X, xvec, i, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        gsl_blas_ddot(xvec, xvec, &xtx);
        
        XtX_diagvec.push_back(xtx);
        gsl_blas_ddot(xvec, y, &xty);
        
        log_lr=((double)ni_test)*(log(yty)-log(yty-xty*xty/xtx));
        pos_loglr.push_back(make_pair(i,log_lr) );
        
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                Gvec[j] += (xtx / (double)ni_test);
                //Gvec[j] += (xtx );
                continue;
            }
        }
	}
	gsl_vector_free(xvec);
    //cout << "trace_G : " << Gvec[0] << ", " << Gvec[1] << endl;
	return;
}

//calculate Z-Scores and SE(beta); used in EM_block
void MatrixCalcLmLR (uchar **X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr, const size_t &ns_test, const size_t &ni_test, const vector<double> &SNPsd, const vector<double> &SNPmean, vector<double> &Gvec, vector<double> &XtX_diagvec, vector<double> &Z_scores, vector<double> &SE_beta, vector<double> &pval_lrt, const vector<SNPPOS> &snp_pos, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag, const vector<pair<int, double> > &UcharTable)
{
    size_t n_type = snp_pos[0].indicator_func.size();
    Gvec.assign(n_type, 0.0);
    XtX_diagvec.clear();
    Z_scores.clear(); //saved by position
    SE_beta.clear();
    pval_lrt.clear();
    
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
	double yty, xty, xtx, log_lr, beta_est, se_beta;
	gsl_blas_ddot(y, y, &yty);
    //cout << "calcLR: yty = " << yty << endl;
    
	for (size_t i=0; i<ns_test; ++i) {
        
        getGTgslVec(X, xvec, i, ni_test, ns_test, SNPsd, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        gsl_blas_ddot(xvec, xvec, &xtx);
        //if(i < 10) cout << xtx << ";";

        XtX_diagvec.push_back(xtx);
        gsl_blas_ddot(xvec, y, &xty);

        //calculate Zscores and se_beta
        beta_est = xty / xtx;
        se_beta = sqrt(yty / ( ((double)ni_test) * xtx ) );
        Z_scores.push_back(beta_est / se_beta);
        SE_beta.push_back(se_beta);

        // calculate likelihood ratio test statistic values        
        log_lr=((double)ni_test)*(log(yty)-log(yty-xty*xty/xtx));
        pos_loglr.push_back(make_pair(i,log_lr) );
        pval_lrt.push_back(gsl_cdf_chisq_Q (log_lr, 1));
        
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                Gvec[j] += (xtx / (double)ni_test);
                //Gvec[j] += (xtx );
                continue;
            }
        }
	}
	gsl_vector_free(xvec);
    //cout << "trace_G : " << Gvec[0] << ", " << Gvec[1] << endl;
	return;
}
