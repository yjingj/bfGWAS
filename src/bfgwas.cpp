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
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <ctime>
#include <cmath>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"  //for functions EigenDecomp
#include "io.h"
#include "bfgwas.h"
#include "lm.h"
#include "bvsrm.h"
#include "mathfunc.h"

using namespace std;


BFGWAS::BFGWAS(void):	
version("bfGWAS_MCMC"), date("08/15/2016"), year("2016")
{}

void BFGWAS::PrintHeader (void)
{
	cout<<endl;
	cout<<"*********************************************************"<<endl;
	cout<<"  Bayesian Functional GWAS --- MCMC (bfGWAS:MCMC) "<<endl;
	cout<<"  Version "<<version<<", "<<date<<"                              "<<endl;
	cout<<"  Visit                                                 "<<endl;
	cout<<"     https://github.com/yjingj/bfGWAS      "<<endl;
	cout<<"  For Possible Updates                                  "<<endl;
	cout<<"  (C) "<<year<<" Jingjing Yang              "<<endl;
	cout<<"  GNU General Public License                            "<<endl;
	cout<<"  For Help, Type ./Estep_mcmc -h                             "<<endl;
	cout<<"*********************************************************"<<endl;
	cout<<endl;
	
	return;
}


void BFGWAS::PrintLicense (void)
{
	cout<<endl;
	cout<<"The Software Is Distributed Under GNU General Public License, But May Also Require The Following Notifications."<<endl;
	cout<<endl;
	
	cout<<"Including Lapack Routines In The Software May Require The Following Notification:"<<endl;
	cout<<"Copyright (c) 1992-2010 The University of Tennessee and The University of Tennessee Research Foundation.  All rights reserved."<<endl;
	cout<<"Copyright (c) 2000-2010 The University of California Berkeley. All rights reserved."<<endl;
	cout<<"Copyright (c) 2006-2010 The University of Colorado Denver.  All rights reserved."<<endl;	
	cout<<endl;
	
	cout<<"$COPYRIGHT$"<<endl;
	cout<<"Additional copyrights may follow"<<endl;
	cout<<"$HEADER$"<<endl;
	cout<<"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:"<<endl;
	cout<<"- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer."<<endl;
	cout<<"- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution."<<endl;
	cout<<"- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission."<<endl;
	cout<<"The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other "
		<<"intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against "
		<<"recipient by any third party for infringement of that parties intellectual property rights. "<<endl;
	cout<<"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT "
		<<"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT "
		<<"OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT "
		<<"LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY "
		<<"THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE "
		<<"OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."<<endl;
	cout<<endl;

	return;
}


void BFGWAS::PrintHelp(size_t option)
{
	if (option==0) {
		cout<<endl; 
		cout<<" bfGWAS_MCMC version "<<version<<", released on "<<date<<endl;
		cout<<" implemented by Jingjing Yang"<<endl; 
		cout<<endl;
		cout<<" type ./Estep_mcmc -h [num] for detailed helps"<<endl;
		cout<<" options: " << endl;
		cout<<" 1: quick guide"<<endl;
		cout<<" 2: file I/O related"<<endl;
		cout<<" 3: variant QC"<<endl;
		cout<<" 4: calculate relatedness matrix"<<endl;
		cout<<" 5: apply single variate linear model"<<endl;
		cout<<" 6: apply BVSR model"<<endl;
		cout<<" 7: save genotype text file or LD matrix"<<endl;
		cout<<" 8: note"<<endl;
		cout<<endl;
	}	
	
	if (option==1) {
		cout<<" QUICK GUIDE " << endl;
		cout<<" Generate a kinship matrix: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -gk [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -vcf [filename] -p [filename] -gk [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -gk [num] -o [prefix]"<<endl<<endl;

		cout<<" Fit a linear model for single variant test: "<<endl<<endl;
		cout<<"         ./gemma -bfile [prefix] -lm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -vcf [filename] -p [filename] -lm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -lm [num] -o [prefix]"<<endl<<endl;

		cout<<" Fit a Bayesian variable selection regression model (BVSRM): "<<endl;
		cout<<"         ./gemma -bfile [prefix] -a [filename] -fcode [filename] -hfile [filename] -bvsrm -o [prefix]"<<endl;
		cout<<"         ./gemma -vcf [filename] -p [filename] -a [filename] -fcode [filename] -hfile [filename] -GTfield [keyword] -bvsrm -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -fcode [filename] -hfile [filename] -bvsrm -o [prefix]"<<endl;
		cout<<endl;
	}
	
	if (option==2) {
		cout<<" FILE I/O RELATED OPTIONS " << endl;
		cout<<" -bfile    [prefix]       "<<" Specify input PLINK binary ped file prefix."<<endl;	
		cout<<"          requires: *.fam, *.bim and *.bed files"<<endl;	
		cout<<"          missing value: -9"<<endl<<endl;

		cout<<" -g        [filename]     "<<" Specify input BIMBAM mean genotype file"<<endl;
		cout<<"          format: rs#1, allele0, allele1, genotype for individual 1, genotype for individual 2, ..."<<endl;	
		cout<<"                  rs#2, allele0, allele1, genotype for individual 1, genotype for individual 2, ..."<<endl;	
		cout<<"                  ..."<<endl;	
		cout<<"          missing value: NA"<<endl<<endl;	

		cout<<" -p        [filename]     "<<" Specify input phenotype file associated with -g and -vcf"<<endl;
		cout<<"          format: individual_1_id    phenotype1"<<endl;	
		cout<<"                  individual_2_id    phenotype2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl<<endl;

		cout<<" -vcf        [filename]     "<<" Specify input VCF genotype file"<<endl<<endl;
		cout<<" -GTfield        [keyword]     "<<" Specify whether to read for genotypes with keyword=GT, or read for doseage data keyword=EC from the VCF genotype file"<<endl<<endl;

		cout<<" -a        [filename]     "<<" Specify input annotation file name (optional)"<<endl;	
		cout<<"          format: rs#1 chr_number base_position anno1"<<endl;	
		cout<<"                  rs#2 chr_number base_position anno2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          in the same order as variants in the genotype file"<<endl;
		cout<<"          missing value: NA or blank"<<endl<<endl;

		cout<<" -fcode        [filename]     "<<" Specify annotation code file"<<endl<<endl;
		cout<<" -hfile        [filename]     "<<" Specify hyper parameter values"<<endl<<endl;

		//cout<<" -k        [filename]     "<<" Specify input kinship/relatedness matrix file name"<<endl<<endl;	

		cout<<" -snps     [filename]     "<<" Specify the variant ids that are to be analyzed"<<endl;
		cout<<"          format: rs#1"<<endl;	
		cout<<"                  rs#2"<<endl;	
		cout<<"                  ..."<<endl<<endl;

		cout<<" -silence                 "<<" Silent terminal display"<<endl<<endl;

		//cout<<" -km       [num]          "<<" Specify input kinship/relatedness file type (default 1)."<<endl;
		//cout<<"          options: 1: \"n by n matrix\" format"<<endl;
		//cout<<"                   2: \"id  id  value\" format"<<endl<<endl;

		//cout<<" -pace     [num]          "<<" Specify terminal display update pace (default 100000 variants or 100000 iterations)."<<endl<<endl;

		cout<<" -o        [prefix]       "<<" Specify output file prefix (default \"result\")"<<endl;  
		cout<<"          output: prefix.cXX.txt or prefix.sXX.txt from kinship/relatedness matrix estimation"<<endl;	
		cout<<"          output: prefix.assoc.txt and prefix.log.txt form association tests"<<endl;	
		cout<<endl;
	}
	
	if (option==3) {
		cout<<" Variant QC OPTIONS " << endl;
		cout<<" -miss     [num]          "<<" specify missingness threshold (default 0.05)" << endl; 
		cout<<" -maf      [num]          "<<" specify minor allele frequency threshold (default 0.005)" << endl; 
		cout<<" -hwe      [num]          "<<" specify HWE test p value threshold (default 0; no test)" << endl; 
		cout<<endl;
	}
	
	if (option==4) {
		cout<<" RELATEDNESS MATRIX CALCULATION OPTIONS " << endl;
		cout<<" -gk       [num]          "<<" specify which type of kinship/relatedness matrix to generate (default 1)" << endl; 
		cout<<"          options: 1: centered XX^T/p"<<endl;
		cout<<"                   2: standardized XX^T/p"<<endl;
		cout<<"          note: non-polymorphic variants are excluded "<<endl;
		cout<<endl;
	}
	
	if (option==5) {
		cout<<" LINEAR MODEL OPTIONS " << endl;		
		cout<<" -lm       [num]         "<<" specify analysis options (default 1)."<<endl;
		cout<<"          options: 1: Wald test"<<endl;
		cout<<"                   2: Likelihood ratio test"<<endl;
		cout<<"                   3: Score test"<<endl;
		cout<<"                   4: 1-3"<<endl;
		cout<<endl;
	}
	
	if (option==6) {
		cout<<" BVSRM OPTIONS " << endl;	
		cout<<" -bvsrm	   "<<" apply BVSR model "<<endl;
		
		//cout<<"   MCMC OPTIONS" << endl;
		//cout<<"   Prior" << endl;
		//cout<<" -hmin     [num]          "<<" specify minimum value for h (default 0)" << endl; 
		//cout<<" -hmax     [num]          "<<" specify maximum value for h (default 1)" << endl; 
		//cout<<" -rmin     [num]          "<<" specify minimum value for rho (default 0)" << endl; 
		//cout<<" -rmax     [num]          "<<" specify maximum value for rho (default 1)" << endl; 
		//cout<<" -pmin     [num]          "<<" specify minimum value for log10(pi) (default log10(1/p), where p is the number of analyzed SNPs )" << endl; 
		//cout<<" -pmax     [num]          "<<" specify maximum value for log10(pi) (default log10(1) )" << endl; 			

		cout<<" -smin     [num]          "<<" specify minimum value for |gamma| (default 0)" << endl; 
		cout<<" -smax     [num]          "<<" specify maximum value for |gamma| (default 300)" << endl;
		cout<<" -rv     [num]          "<<" specify the residual variance value (default 1 for standardized phenotypes)" << endl; 		
		cout<<" -win   [num]          "<<" specify the neighbourhood window size (default 100)" << endl; 
		cout<<" -w        [num]          "<<" specify burn-in steps (default 100,000)" << endl; 
		cout<<" -s        [num]          "<<" specify sampling steps (default 1,000,000)" << endl; 
		cout<<" -comp                    "<<" specify whether to impement in-memory compression (default no comp flag)" << endl;
		
		cout<<" -initype     [num]          "<<" specify the initial model for MCMC: " << endl;
		cout<<"          option 1: start with the most significant variant"<<endl;
		cout<<"          option 2: start with the significant variants that acheive genome-wide significance"<<endl;
		cout<<"          option 3: start with a model stepwise selection"<<endl;

		cout<<" -seed     [num]          "<<" specify random seed (a random seed is generated by default)" << endl; 	
		cout<<" -mh       [num]          "<<" specify number of MH steps in each iteration (default 10)" << endl; 
		cout<<"          requires quantitative phenotypes (0/1 for case/control)"<<endl;	
		cout<<endl<<endl;
	}

	if (option==7) {
		cout<<" -saveGeno  "<<" specify whether to save a tab delimited genotype text file from the inpute genotype file" << endl;
		cout<<" -saveLD  "<<" specify whether to save LD matrix from the analyzed variants" << endl;
		cout<<" -saveSNP  "<<" specify whether to save a tab delimited genotype text file for the analyzed variants in bvsrm" << endl;
		cout<<endl;
	}
	
	if (option==8) {
		cout<<" NOTE"<<endl;
		cout<<" 1. Only individuals with quantitative non-NA phenotoypes will be analyzed."<<endl;
		cout<<" 2. Missing genotoypes will be repalced by the sample mean genotype."<<endl;
		cout<<endl;
	}
	
	return;
}



void BFGWAS::Assign(int argc, char ** argv, PARAM &cPar)
{
	string str;

	for(int i = 1; i < argc; i++) {		
		if (strcmp(argv[i], "-bfile")==0 || strcmp(argv[i], "-b")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_bfile=str;
		}
		else if (strcmp(argv[i], "-silence")==0) {
			cPar.mode_silence=true;
		}
		else if (strcmp(argv[i], "-g")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_geno=str;
		}
        else if (strcmp(argv[i], "-vcf")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_vcf=str;
		}
        else if (strcmp(argv[i], "-iniSNP")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.iniSNPfile=str;
        }
        else if (strcmp(argv[i], "-hfile")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.hypfile=str;
        }
        else if (strcmp(argv[i], "-GTfield")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.GTfield=str;
        }
        else if (strcmp(argv[i], "-fsample")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.file_sample=str;
        }
		else if (strcmp(argv[i], "-p")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_pheno=str;
		}
		else if (strcmp(argv[i], "-a")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_anno=str;
		}
        else if (strcmp(argv[i], "-fcode")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_func_code=str;
		}
		else if (strcmp(argv[i], "-k")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_kin=str;
		}
		else if (strcmp(argv[i], "-snps")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_snps=str;
		}
		else if (strcmp(argv[i], "-km")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.k_mode=atoi(str.c_str());
		}		
		else if (strcmp(argv[i], "-o")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_out=str;
		}		
		else if (strcmp(argv[i], "-miss")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.miss_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-maf")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			if (cPar.maf_level!=-1) {cPar.maf_level=atof(str.c_str());}
		}
		else if (strcmp(argv[i], "-hwe")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.hwe_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-notsnp")==0) {
			cPar.maf_level=-1;
		}
		else if (strcmp(argv[i], "-gk")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -lm -bvsrm options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=21; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=20+atoi(str.c_str());
		}	
		else if (strcmp(argv[i], "-lm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -lm -bvsrm options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=51; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=50+atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-bvsrm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -lm -bvsrm options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=11; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=10+atoi(str.c_str());
		}
        else if (strcmp(argv[i], "-vscale")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.vscale=atof(str.c_str());
        }
        else if (strcmp(argv[i], "-rv")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.rv=atof(str.c_str());
        }
		else if (strcmp(argv[i], "-w")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.w_step=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-s")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_step=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-smin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_min=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-smax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_max=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-hmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-hmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-pmin")==0) {
			if(argv[i+1] == NULL) {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.logp_min=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-pmax")==0) {
			if(argv[i+1] == NULL) {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.logp_max=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-hscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_scale=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_scale=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-pscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
            cPar.logp_scale=atof(str.c_str());
			//cPar.logp_scale=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-seed")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.randseed=atol(str.c_str());
		}
		else if (strcmp(argv[i], "-mh")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n_mh=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-pace")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.d_pace=atoi(str.c_str());
		}
        else if (strcmp(argv[i], "-win")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.win=100;continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.win=atol(str.c_str());
		}
        else if (strcmp(argv[i], "-initype")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.iniType=3;continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.iniType=atoi(str.c_str());
        }
        else if (strcmp(argv[i], "-fixhyp")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.FIXHYP=0;continue;}
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.FIXHYP=atoi(str.c_str());
        }
        else if (strcmp(argv[i], "-saveSNP")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.saveSNP=1; }
        }
        else if (strcmp(argv[i], "-saveGeno")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.saveGeno=1;}
        }
        else if (strcmp(argv[i], "-saveLD")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.saveLD=1; }
        }
        else if (strcmp(argv[i], "-comp")==0) {
            if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.Compress_Flag=1;}
        }
		else {cout<<"error! unrecognized option: "<<argv[i]<<endl; cPar.error=true; continue;}
	}
	
	return;
}



void BFGWAS::BatchRun (PARAM &cPar) 
{
	clock_t time_begin=clock(), time_start=clock();
    
    //cout << "create UcharTable ...\n";
    //vector<pair<long long int, double> > UcharTable;
    //CreateUcharTable(UcharTable);
    	
	//Read Files for the first time and filt variants
	cout<<"Reading Files ... " << endl;
	cPar.ReadFiles();
	if (cPar.error==true) {cout<<"error! fail to read files. "<<endl; return;}
    cout << "Readin files first time cost " << (clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0) << "mints \n";

    //Save Genotype file 
    if(cPar.saveGeno){

		//Save all genotypes
		cPar.indicator_snp.assign(cPar.ns_total, 1);
		cPar.indicator_idv.assign(cPar.ni_total, 1);
		cPar.ni_test = cPar.ni_total;
		cPar.ns_test = cPar.ns_total;

		cPar.CheckData();
		if (cPar.error==true) {cout<<"error! fail to check data. "<<endl; return;}
    	cout << "Pass data check." << endl;

		gsl_matrix *G=gsl_matrix_alloc (cPar.ni_test, cPar.ni_test);
		gsl_vector *y=gsl_vector_alloc (cPar.ni_test); // phenotype
		
		//set phenotype vector y		
		cout << "copy phenotype success ... "<< endl;
		cPar.CopyPheno (y);

		if ( (!cPar.file_vcf.empty()) || (!cPar.file_geno.empty()) ) {
        	// reorder y for reading vcf/genotype files
        	cout << "Reorder phenotype for reading vcf files ... "<< endl;
        	cPar.ReorderPheno(y);
    	} // reorder y for reading vcf files

        //read genotypes X 
        clock_t time_readfile = clock();
        uchar ** X_Genotype = new uchar*[cPar.ns_test];
        cPar.ReadGenotypes (X_Genotype, G, false); 
            
        cout << "load genotype data cost " << (clock()-time_readfile)/(double(CLOCKS_PER_SEC)*60.0) << "mints\n";
        
        cPar.WriteGenotypes(X_Genotype);

        gsl_matrix_free(G);
        gsl_vector_free(y);

        cout << "writting genotype file success ... "<< endl; 
        //exit(EXIT_SUCCESS);
    }else{
    	cPar.CheckData();
		if (cPar.error==true) {cout<<"error! fail to check data. "<<endl; return;}
    	cout << "Pass data check." << endl;
    }

    
	
	//Generate Kinship matrix
	if (cPar.a_mode==21 || cPar.a_mode==22) {  
		cout<<"Calculating Relatedness Matrix ... "<<endl;
		
		gsl_matrix *G=gsl_matrix_alloc (cPar.ni_test, cPar.ni_test);
		
		time_start=clock();
		cPar.CalcKin (G);
		cPar.time_G=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		if (cPar.error==true) {cout<<"error! fail to calculate relatedness matrix. "<<endl; return;}
		
		if (cPar.a_mode==21) {
			string file_str = "./output/" + cPar.file_out;
			file_str += ".cXX.txt";
			WriteMatrix (G, file_str);
		} else {
			string file_str = "./output/" + cPar.file_out;
			file_str += ".sXX.txt";
			WriteMatrix (G, file_str);
		}
		
		gsl_matrix_free (G);
	}
	
	
	//LM
	if (cPar.a_mode==51 || cPar.a_mode==52 || cPar.a_mode==53 || cPar.a_mode==54) {  

		gsl_vector *Y=gsl_vector_alloc (cPar.ni_test); //set phenotype vector Y
		gsl_matrix *W=gsl_matrix_alloc (cPar.ni_test, 1); // intercept column of 1's; or covariates
		gsl_matrix_set_all(W, 1);

		cPar.CopyPheno (Y);
		if ( (!cPar.file_vcf.empty()) || (!cPar.file_geno.empty()) ) {
        	// reorder y for reading vcf/genotype files
        	cout << "Reorder y for reading vcf or geno files ... "<< endl;
        	cPar.ReorderPheno(Y);
    	}

		//Fit LM 
			LM cLm;
			cLm.CopyFromParam(cPar);
						
			if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, Y);
			} 
			else if( !cPar.file_geno.empty() ) {
				cLm.AnalyzeGeno (W, Y, cPar.SampleVcfPos, cPar.PhenoID2Ind, cPar.VcfSampleID);
			} 
			else if( !cPar.file_vcf.empty() ) {
				cLm.AnalyzeVCF(W, Y, cPar.GTfield, cPar.SampleVcfPos, cPar.PhenoID2Ind, cPar.VcfSampleID);
			}
			
			cLm.WriteFiles();
			cLm.CopyToParam(cPar);

		//release all matrices and vectors
		gsl_vector_free (Y);
		gsl_matrix_free (W);
	} 

	
	//BVSRM
	if (cPar.a_mode==11 ) {
        
		gsl_vector *y=gsl_vector_alloc (cPar.ni_test); // phenotype
		gsl_matrix *W=gsl_matrix_alloc (y->size, 1); // intercept column of 1's
		gsl_matrix_set_all(W, 1); //covariate matx
		gsl_matrix *G=gsl_matrix_alloc (y->size, y->size); // kinship matrix
		
		// set phenotype vector y		
		// cout << "copy phenotype success ... "<< endl;
		cPar.CopyPheno (y);
        
        if ( (!cPar.file_vcf.empty()) || (!cPar.file_geno.empty()) ) {
        	// reorder y for reading vcf/genotype files
        	cout << "Reorder y for reading vcf files ... "<< endl;
        	cPar.ReorderPheno(y);
    	}
        
        //center y, even for case/control data
        cout << "Center phenotype ... "<< endl;
        cPar.pheno_mean=CenterVector(y);
        cout << "pheno_mean = " << cPar.pheno_mean << "\n";
        //PrintVector(y);
        
		//run bvsrm as if rho==1
        //read genotypes X (not UtX)
            clock_t time_readfile = clock();
 
            //cout << "allocate uchar* vector with length : " << cPar.ns_test << endl;
            uchar ** X_Genotype = new uchar*[cPar.ns_test];
            cPar.ReadGenotypes (X_Genotype, G, false); //load genotypes
            
            cout << "load genotype data cost " << (clock()-time_readfile)/(double(CLOCKS_PER_SEC)*60.0) << "mints\n";
            
           // PrintVector(cPar.CompBuffSizeVec);
            //cout << "UnCompBufferSize = " << cPar.UnCompBufferSize << endl;
           //print(X_Genotype, 10, 10, cPar.CompBuffSizeVec, cPar.UnCompBufferSize);
           //print(X_Genotype, 10, 10);

            gsl_matrix_free(G);
            gsl_matrix_free(W);

            //perform BSVRM analysis
            BVSRM cBvsrm;
           // cout << "copy data from param ...\n";
            cBvsrm.CopyFromParam(cPar);
            cBvsrm.ns_neib = 2 * cBvsrm.win + 1;
         
           // cout << "start bslmm.mcmc ...\n";
            time_start=clock();
            cBvsrm.MCMC(X_Genotype, y, 1);
            cPar.time_opt=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            cBvsrm.CopyToParam(cPar);
            
            FreeUCharMatrix(X_Genotype, cPar.ns_test);
                 
		gsl_vector_free (y);
			
    }
		
	cPar.time_total=(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0);
	
	return;
}


void BFGWAS::WriteLog (int argc, char ** argv, PARAM &cPar) 
{
	string file_str;
	file_str="./output/"+cPar.file_out;
	file_str+=".log.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing log file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"##"<<endl;
	outfile<<"## GEMMA Version = "<<version<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Command Line Input = ";
	for(int i = 1; i < argc; i++) {	
		outfile<<argv[i]<<" ";
	}
	outfile<<endl;

	outfile<<"##"<<endl;
	time_t  rawtime; 
	time(&rawtime);
	tm *ptm = localtime (&rawtime);

	outfile<<"## Date = "<<asctime(ptm)<<endl;
	  //ptm->tm_year<<":"<<ptm->tm_month<<":"<<ptm->tm_day":"<<ptm->tm_hour<<":"<<ptm->tm_min<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Summary Statistics:"<<endl;
	outfile<<"## number of total individuals = "<<cPar.ni_total<<endl;	
	outfile<<"## number of analyzed individuals = "<<cPar.ni_test<<endl;
	outfile<<"## number of total SNPs = "<<cPar.ns_total<<endl;	
	outfile<<"## number of analyzed SNPs = "<<cPar.ns_test<<endl;
	
	if (cPar.a_mode==11) {
		outfile<<"## Phenotype mean = "<<cPar.pheno_mean<<endl;	
		outfile<<"## Phenotype var = "<<cPar.pheno_var<<endl;
		outfile<<"##"<<endl;
		outfile<<"## MCMC related:"<<endl;	
		//outfile<<"## initial value of h = "<<cPar.cHyp_initial.h<<endl;
		//outfile<<"## initial value of rho = "<<cPar.cHyp_initial.rho<<endl;
		//outfile<<"## initial value of pi = "<<exp(cPar.cHyp_initial.logp)<<endl;
		outfile<<"## initial value of |gamma| = "<<cPar.cHyp_initial.n_gamma<<endl;
		outfile<<"## random seed = "<<cPar.randseed<<endl;
		outfile<<"## acceptance ratio = "<<(double)cPar.n_accept/(double)((cPar.w_step+cPar.s_step)*cPar.n_mh)<<endl;

		outfile<<"## Region_PIP = "<<(double)cPar.region_pip/(double)(cPar.s_step)<<endl;

	}else if (cPar.a_mode >= 51 && cPar.a_mode <= 54){
		outfile<<"## Phenotype mean = "<<cPar.pheno_mean<<endl;	
		outfile<<"## Phenotype var = "<<cPar.pheno_var<<endl;	
	}
	
	outfile<<"##"<<endl;
	outfile<<"## Computation Time:"<<endl;
	outfile<<"## total computation time = "<<cPar.time_total<<" min "<<endl;
	outfile<<"## computation time break down: "<<endl;


	if (cPar.a_mode==11) {
		outfile<<"##      time on proposal = "<<cPar.time_Proposal<<" min "<<endl;
		outfile<<"##      time on mcmc = "<<cPar.time_opt<<" min "<<endl;
		outfile<<"##      time on Omega = "<<cPar.time_Omega<<" min "<<endl;
        
        outfile << "Accept #add=" << cPar.nadd_accept<< ", total add step = " << cPar.nadd<<endl;
        outfile << "Accept #delete=" << cPar.ndel_accept<< ", total delete step = " << cPar.ndel<<endl;
        outfile << "Accept #switch=" << cPar.nswitch_accept<< ", total switch step = " << cPar.nswitch<<endl;
        outfile << "Accept #other=" << cPar.nother_accept<< ", total other step = " << cPar.nother<<endl;
	}

	outfile<<"##"<<endl;
	
	outfile.close();
	outfile.clear();
	return;
}


