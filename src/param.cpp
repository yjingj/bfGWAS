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
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"
#include "gzstream.h"
#include "mathfunc.h"
#include "ReadVCF.h"
#include "bvsrm.h"

#include "param.h"
#include "io.h"


using namespace std;


void genMarker::iniRecord(VcfRecord& record){
    
    rs = record.getIDStr();
    chr = record.getChromStr();
    bp = record.get1BasedPosition();
    Ref = record.getRefStr();
    Alt = record.getAltStr();
    
    return;
}

void genMarker::printMarker(){
    
    std::cout << "ID : " << rs <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << bp <<"; ";
    std::cout << "Ref : " << Ref <<"; ";
    std::cout << "Alt : " << Alt <<"\n";
    
    return;
}

void SNPPOS::printMarker(){
    std::cout << "position : " << pos << "; ";
    std::cout << "ID : " << rs <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << bp <<"; ";
    std::cout << "minor allel : " << a_minor <<"; ";
    std::cout << "major allel : " << a_major <<"; \n";
}

void SNPINFO::printMarker(){
    
    std::cout << "ID : " << rs_number <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << base_position <<"; ";
    std::cout << "Ref : " << a_major <<"; ";
    std::cout << "Alt : " << a_minor <<"\n";
    std::cout << "missingness = " << missingness<< "; maf = " << maf << "\n";
}

void printSNPInfo(vector<SNPPOS> &snp_pos, int numMarker)
{
    for (int i=0; i<numMarker; i++) {
        snp_pos[i].printMarker();
    }
}

void CalcWeight(const vector<bool> &indicator_func, vector<double> &weight, const double weight_i)
{
    weight.clear();
    for (size_t i=0; i < indicator_func.size(); i++) {
            if (indicator_func[i]) weight.push_back(weight_i);
            else weight.push_back(0.0);
        }
    if (weight.size() != indicator_func.size()) {
        cout << "Error weight size.\n";
    }
}


PARAM::PARAM(void):
vscale(0.0), iniType(3), saveSNP(0), saveGeno(0), saveLD(0), rv(1.0), 
Compress_Flag(0),
mode_silence (false), a_mode (0), k_mode(1), d_pace (100000),
GTfield("GT"), file_out("result"), 
miss_level(0.05), maf_level(0.005), hwe_level(0.001), 
win(100),nadd_accept(0), ndel_accept(0), nswitch_accept(0),
nother_accept(0), nadd(0), ndel(0),
nswitch(0), nother(0),
h_min(-1), h_max(1.0), h_scale(-1),
rho_min(1.0), rho_max(1.0),	rho_scale(-1),
logp_min(0.0), logp_max(0.0), logp_scale(-1),
s_min(0), s_max(10), 
w_step(50000),	s_step(500000), n_accept(0),
n_mh(10), randseed(2016), error(false),
time_total(0.0), time_G(0.0), time_Omega(0.0)
{}


//read files
//obtain ns_total, ng_total, ns_test, ni_test, n_type
void PARAM::ReadFiles (void) 
{
	string file_str;
	
	// read the set of to be analyzed SNP ids if given file_snps
	if (!file_snps.empty()) {
		if (ReadFile_snps (file_snps, setSnps)==false) {error=true;}
	} else {
		setSnps.clear();
	}

	//read bed/bim/fam files of plink format
	if (!file_bfile.empty()) {
		cout << "Start reading plink bim/fam files ...\n";
		file_str=file_bfile+".bim";
		if (ReadFile_bim (file_str, snpInfo)==false) {error=true;}
        
		file_str=file_bfile+".fam";
		if (ReadFile_fam (file_str, indicator_idv, pheno, InputSampleID, ni_total)==false) {error=true;}
		
		// obtain ni_test, ni_total, PhenoID2Ind before reading genotypes
		ProcessPheno();
		
		file_str=file_bfile+".bed";
		cout << "First time reading Plink bed file: " << file_str << "\n";
		if (ReadFile_bed (file_str, setSnps, indicator_idv, indicator_snp, snpInfo, PhenoID2Ind, ni_test, ni_total, maf_level, miss_level, hwe_level, ns_test, ns_total)==false) {error=true;}
    }else{
    	if (!file_pheno.empty()){
    		cout << "Start reading pheno file ...\n";
        	if (ReadFile_pheno (file_pheno, indicator_idv, pheno, InputSampleID, ni_total)==false)
            	{error=true;}
        	ProcessPheno(); 
        	// obtain ni_test, ni_total, PhenoID2Ind before reading genotypes
    	}else{
    		cout << "No phenotype input file, extracting sample information from the vcf/genotype files.\n";
    		if (!file_vcf.empty()) {
        		getIDVvcf(file_vcf, indicator_idv,  ni_total, GTfield);
      		}else if (!file_geno.empty()) {
				getIDVgeno(file_geno, indicator_idv, ni_total) ;
	  		}else{
	  			cerr << "Unable to get sample information!" << endl;
	  			exit(-1);
	  		}
    	}
    
      	//read vcf file for genotypes
      	if (!file_vcf.empty()) {
        	cout << "Start reading vcf file first time ...\n";
        	indicator_snp.clear();
        	snpInfo.clear();
        	if (ReadFile_vcf(file_vcf, setSnps, indicator_idv, indicator_snp, maf_level, miss_level, hwe_level, snpInfo, ns_test, ns_total, ni_test, GTfield, PhenoID2Ind, VcfSampleID, SampleVcfPos) == false )
            	{error=true;}
      	}else if (!file_geno.empty()) {
	  		//read genotype file 
	  		cout << "Start reading genotype file first time ...\n";
			if (ReadFile_geno (file_geno, setSnps, indicator_idv, indicator_snp, PhenoID2Ind, snpInfo, VcfSampleID, SampleVcfPos, maf_level, miss_level, hwe_level, ns_test, ns_total, ni_test, ni_total)==false) {error=true;}
	  	}
	}

    if ( (!file_anno.empty()) && (!file_func_code.empty()) ) {
    	cout << "Start reading annotation files: " << file_anno << " \nwith code file " << file_func_code << "\n";
        if (ReadFile_anno (file_anno, file_func_code, mapFunc2Code, indicator_snp, snpInfo, n_type, mFunc)==false) {error=true;}
    }
    else {
    	if (Empty_anno (indicator_snp, snpInfo, n_type, mFunc)==false) {error=true;}
    } 

	return;
}



void PARAM::CheckParam (void) 
{	
	struct stat fileInfo;
	string str;
	
	//check parameters
	if (k_mode!=1 && k_mode!=2) {cout<<"error! unknown kinship/relatedness input mode: "<<k_mode<<endl; error=true;}

	if (a_mode!=11 && a_mode!=12 && a_mode!=21 && a_mode!=22 && a_mode!=43 && a_mode!=51 && a_mode!=52 && a_mode!=53 && a_mode!=54 && !saveGeno)   
	{cout<<"error! unknown analysis mode: "<<a_mode<<". make sure -saveGenoe -gk or -lm or -bvsrm or -predict is sepcified correctly."<<endl; error=true;}

	if (miss_level>1) {cout<<"error! missing level needs to be between 0 and 1. current value = "<<miss_level<<endl; error=true;}
	if (maf_level>0.5) {cout<<"error! maf level needs to be between 0 and 0.5. current value = "<<maf_level<<endl; error=true;}
	if (hwe_level>1) {cout<<"error! hwe level needs to be between 0 and 1. current value = "<<hwe_level<<endl; error=true;}
	
	if (h_max<h_min) {cout<<"error! maximum h value must be larger than the minimal value. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	if (rho_max<rho_min) {cout<<"error! maximum rho value must be larger than the minimal value. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max<logp_min) {cout<<"error! maximum logp value must be larger than the minimal value. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}
	
	if (h_max>1) {cout<<"error! h values must be bewtween 0 and 1. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (rho_max>1) {cout<<"error! rho values must be between 0 and 1. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max>0) {cout<<"error! maximum logp value must be smaller than 0. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}
		
	if (h_scale>1.0) {cout<<"error! hscale value must be between 0 and 1. current value = "<<h_scale<<endl; error=true;}
	if (rho_scale>1.0) {cout<<"error! rscale value must be between 0 and 1. current value = "<<rho_scale<<endl; error=true;}
	if (logp_scale>1.0) {cout<<"error! pscale value must be between 0 and 1. current value = "<<logp_scale<<endl; error=true;}

	if (rho_max==1 && rho_min==1 && a_mode==12) {cout<<"error! ridge regression does not support a rho parameter. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
		
	
	//check if files are compatible with each other, and if files exist
	if (!file_bfile.empty()) {
		str=file_bfile+".bim";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bim file: "<<str<<endl; error=true;}
		str=file_bfile+".bed";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bed file: "<<str<<endl; error=true;}
		str=file_bfile+".fam";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .fam file: "<<str<<endl; error=true;}			
	}
	

	str=file_geno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open mean genotype file: "<<str<<endl; error=true;}

	size_t flag=0;
	if (!file_bfile.empty()) {flag++;}
	if (!file_geno.empty()) {flag++;}
	
	if (file_pheno.empty() && (a_mode==43) ) {
		cout<<"error! phenotype file is required."<<endl; error=true;
	}
	
	str=file_snps;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open snps file: "<<str<<endl; error=true;}
	
	
	str=file_anno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open annotation file: "<<str<<endl; error=true;}

	str=file_kin;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open relatedness matrix file: "<<str<<endl; error=true;}
	
	//check if files are compatible with analysis mode

	if ((a_mode==43) && file_kin.empty())  {cout<<"error! missing relatedness file. -predict option requires -k option to provide a relatedness file."<<endl;  error=true;}
	
	return;
}


		

void PARAM::CheckData (void) {


	//calculate ni_total and ni_test, and set indicator_idv to 0 whenever indicator_cvt=0
	//and calculate np_obs and np_miss

	ni_total=(indicator_idv).size();
	
	ni_test=0; 
	for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
		if (indicator_idv[i]==0) {continue;}
		ni_test++;
	}

	np_obs=0; np_miss=0;
	for (size_t j=0; j<indicator_pheno.size(); j++) {					
			if (indicator_pheno[j]==0) {
				np_miss++;
			} else {
				np_obs++;
			}
	}

	if (ni_test==0) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	//output some information
	cout<<"## number of total individuals = "<<ni_total<<endl;
	cout<<"## number of individuals with full phenotypes = "<<ni_test<<endl;

	if (a_mode==43) {
		cout<<"## number of observed data = "<<np_obs<<endl;
		cout<<"## number of missing data = "<<np_miss<<endl;
	}
	
	cout<<"## number of total SNPs = "<<ns_total<<endl;	
	cout<<"## number of analyzed SNPs = "<<ns_test<<endl;

	
	//set parameters for BSLMM
	//and check for predict
	if (a_mode==11 || a_mode==12 || a_mode==13) {
		if (a_mode==11) {n_mh=1;}	
		if (logp_min==0) {logp_min=-1.0*log((double)ns_test);}
	
		if (h_scale==-1) {h_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (rho_scale==-1) {rho_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (logp_scale==-1) {logp_scale=min(1.0, 5.0/sqrt((double)ni_test) );}
        //cout << "h_scale = " << h_scale << "; rho_scale = " << rho_scale<< "; logtheta_scale = " << logp_scale << endl;
        if (vscale <= 0.0) { vscale = min(0.5, 10.0/sqrt((double)ni_test));}
		
		if (h_min==-1) {h_min=0.00000001;}
		if (h_max==-1) {h_max=1.0;}
		
		if (s_max>ns_test) {s_max=ns_test; cout<<"s_max is re-set to the number of analyzed SNPs."<<endl;}
		if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	} 
	
	return;
}


void PARAM::PrintSummary () 
{
	cout<<"pve estimate ="<<endl;
	cout<<"se(pve) ="<<endl;
		
	return;
}


void PARAM::ReadGenotypes (uchar **X, gsl_matrix *K, const bool calc_K) {
 
    string file_str;
    UnCompBufferSize = (ni_test) * sizeof(uchar);
   // cout << "UnCompBufferSize = " << UnCompBufferSize << endl;
    
	if (!file_bfile.empty()) {
		file_str=file_bfile+".bed";
		if (ReadFile_bed (file_str, indicator_idv, indicator_snp, X, K, calc_K, ni_test, ns_test, ni_total, ns_total, SNPmean, CompBuffSizeVec, Compress_Flag)==false) {error=true;}
        //revised
	}
	
    else if(!file_vcf.empty()){
        if ( ReadFile_vcf (file_vcf, indicator_idv, indicator_snp, X, ni_test, ns_test, K, calc_K, GTfield, SNPmean, CompBuffSizeVec, SampleVcfPos, PhenoID2Ind, VcfSampleID, Compress_Flag)==false )
        {error=true;} // revised
    }
    
    else if(!file_geno.empty()){
        if (ReadFile_geno (file_geno, indicator_idv, indicator_snp, X, K, calc_K, ni_test, SNPmean, CompBuffSizeVec, SampleVcfPos, PhenoID2Ind, VcfSampleID, Compress_Flag)==false) {error=true;} //to be revised
    }else{
    	cerr << "one of the genotype files has to be specified." << endl;
    	exit(-1);
    }
    
    return;
    
}

void PARAM::WriteGenotypes(uchar **X){

	string file_str;
    file_str="./output/"+file_out;
    file_str+=".geno";

    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    //write header with VcfSampleID_test 
    //cout << "write header row."<< endl;
    //cout << "ni_test = "<< ni_test << endl;

    outfile<<"ID"<<"\t"<<"CHROM"<<"\t" <<"POS"<<"\t" << "REF"<< "\t" << "ALT"  << "\t";

    for (size_t i=0; i<ni_test; i++) {
    	//if(i < 10){ cout << VcfSampleID_test[i] <<endl; }
        if (i ==(ni_test-1)) {
            outfile << VcfSampleID_test[i] << endl;
        }
        else outfile << VcfSampleID_test[i] << "\t";
    } 

    size_t pos=0;
    double geno_j;
    uchar c;
    
    //cout << "write variant information."<<endl;
    for (size_t i=0; i<ns_total; ++i) {

    	if(!indicator_snp[i]){continue;}
        
    	// save the data 
        outfile<< snpInfo[i].rs_number <<"\t"<< snpInfo[i].chr<<"\t" <<snpInfo[i].base_position << "\t" << snpInfo[i].a_major << "\t" << snpInfo[i].a_minor << "\t";
 
        for (size_t j=0; j < ni_test; j++) {
        	c = X[pos][j];
            geno_j = UcharTable[(int)c].second;
            if(geno_j < 0.0 || geno_j > 2){
            	cout << "ERROR: genotype = " << geno_j <<" for " << snpInfo[i].rs_number <<":"<< snpInfo[i].chr<<":" <<snpInfo[i].base_position << ":" << snpInfo[i].a_major << ":" << snpInfo[i].a_minor << endl;
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
        pos++;
    }
    
    outfile.clear();
    outfile.close();

}


// Calculate kinshiip matrix
void PARAM::CalcKin (gsl_matrix *matrix_kin)  {
	string file_str;
	
	gsl_matrix_set_zero (matrix_kin);
	
	if ( !file_bfile.empty() ) {		
		file_str=file_bfile+".bed";
		if (PlinkKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin)==false) {error=true;}
	}
	else if( !file_geno.empty() ) {
		file_str=file_geno;
		if (GenoKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin, SampleVcfPos, PhenoID2Ind, VcfSampleID)==false) {error=true;}
	}
	else if( !file_vcf.empty() ) {
		file_str=file_vcf;
		if (VCFKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin, GTfield, SampleVcfPos, PhenoID2Ind, VcfSampleID)==false) {error=true;}
	}
	else{
		cerr << "Need input genotype file: plink, bimbam, or VCF!" <<endl;
	}
	
	return;
}
		




/*
void PARAM::CheckCvt () 
{
	if (indicator_cvt.size()==0) {return;}
		
	size_t ci_test=0;
	
	gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);
	
	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}
		ci_test++;
	}

	size_t flag_ipt=0;
	double v_min, v_max;
	set<size_t> set_remove;
	
	//check if any columns is an intercept
	for (size_t i=0; i<W->size2; i++) {
		gsl_vector_view w_col=gsl_matrix_column (W, i);
		gsl_vector_minmax (&w_col.vector, &v_min, &v_max);
		if (v_min==v_max) {flag_ipt=1; set_remove.insert (i);}
	}
	
	//add an intecept term if needed
	if (n_cvt==set_remove.size()) {
		indicator_cvt.clear();
		n_cvt=1;
	} else if (flag_ipt==0) {
		cout<<"no intecept term is found in the cvt file. a column of 1s is added."<<endl;
		for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
			if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
			cvt[i].push_back(1.0);
		}
		
		n_cvt++;
	} else {}	
	
	gsl_matrix_free(W);
	
	return;
}
*/

//reorder phenotypes
void PARAM::ReorderPheno(gsl_vector *y)
{
    if (VcfSampleID.size() < ni_test) {
        cerr << "Sample size in genotype file" <<  VcfSampleID.size() << "< ni_test: " << ni_test << endl;
        exit(-1);
    }
    
    double pheno_i;
    string id;
    size_t c_ind=0, pheno_idx;
    gsl_vector *ytemp=gsl_vector_alloc (ni_test);
    VcfSampleID_test.clear(); // set VCFSampleID_test
    
    //cout << "Number of samples in the genotype file: " << VcfSampleID.size() << endl;
	//cout << "PhenoID2Ind.size() = " << PhenoID2Ind.size() << " within ReorderPheno() function ... " << endl;

	// indicator_idv is of the same order as in the phenotype file
    for (size_t i=0; i < VcfSampleID.size(); i++) {
        id = VcfSampleID[i];

        if(PhenoID2Ind.count(id) > 0){
        	//if (i < 10) cout << id << ", count ="<< PhenoID2Ind.count(id) << ";    "; 
        	pheno_idx = PhenoID2Ind.at(id);

        	if (indicator_idv[pheno_idx] ) {
            	pheno_i = gsl_vector_get(y, pheno_idx);
            	gsl_vector_set(ytemp, c_ind, pheno_i);
            	VcfSampleID_test.push_back(id);
            	c_ind++;
        	}else{
        		indicator_idv[pheno_idx] = 0;
        	}
        }
    }
    ni_test = c_ind;
    cout << "\n Reorder phenotypes, final analyzed sample size ni_test = " << ni_test << endl;
    gsl_vector_memcpy(y, ytemp);
    gsl_vector_free(ytemp);
}


//post-process phentoypes, obtain ni_total, ni_test, PhenoIDInd
void PARAM::ProcessPheno()
{	
	//obtain ni_test
	ni_test=0;

	//obtain PhenoID2Ind map
    PhenoID2Ind.clear();

	for (size_t i=0; i<ni_total; ++i) {
        if (indicator_idv[i]){
            PhenoID2Ind[InputSampleID[i]]= ni_test;
            ni_test++;
        }
        else PhenoID2Ind[InputSampleID[i]] = ULONG_MAX;
    }
    //cout << "Create PhenoID2Ind map; total number of analyzed individual ni_test = " << ni_test << "\n";
    //cout << "PhenoID2Ind map length = " << PhenoID2Ind.size() << "\n";
	
	if (ni_test==0) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	return;
}



//else, indicator_idv to load phenotype
void PARAM::CopyPheno (gsl_vector *y) 
{
	size_t ci_test=0;
	pheno_mean = 0.0;
	
	for (size_t i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0) {continue;}
		gsl_vector_set (y, ci_test, pheno[i]);
		pheno_mean += pheno[i];
		ci_test++;
	}
	pheno_mean /= (double) ci_test;
	
	return;
}





		

