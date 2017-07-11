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
#include <iomanip>
#include <bitset>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <ctype.h>

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
#include "compress.h"
#include "io.h"


using namespace std;

// define to_string function : convert to string
template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}


//Print process bar
void ProgressBar (string str, double p, double total)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <30; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%\r"<<flush;
	
	return;
}


//Print process bar (with acceptance ratio)
void ProgressBar (string str, double p, double total, double ratio)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <30; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%  "<< "& acceptance ratio "<<ratio<<"\r"<<flush;
	
	
	return;
}

// in case files are ended with "\r" or "\r\n"
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

//Read snp file
bool ReadFile_snps (const string &file_snps, set<string> &setSnps)
{
	setSnps.clear();

	ifstream infile (file_snps.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open snps file: "<<file_snps<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		setSnps.insert(ch_ptr); 
	}
	
    infile.clear();
	infile.close();
	
	return true;
}


void SetMAFCode (const double &maf, string &func_type){
	if( maf >= 0.005 && maf < 0.01) func_type+="-maf-range1";
    else if ( maf >= 0.01 && maf < 0.05) func_type+="-maf-range2";
    else if ( maf >= 0.05 && maf < 0.1) func_type+="-maf-range3";
    else if ( maf >= 0.1 && maf < 0.2) func_type+="-maf-range4";
    else if ( maf >= 0.2 && maf < 0.3) func_type+="-maf-range5";
    else if ( maf >= 0.3 && maf < 0.4) func_type+="-maf-range6";
    else if ( maf >= 0.4 ) func_type+="-maf-range7";
}

//Read function annotation file
bool ReadFile_anno (const string &file_anno, const string &file_func_code, map<string, int> &mapFunc2Code, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, size_t &n_type, vector<size_t> &mFunc)
{
    string line;
    char *pch, *nch;

    //load in unique function codes
    string func_type;
    int func_code, snp_nfunc;
    
    // Load function_code file first, create a hash map between func_type and code
    // cout<<"Reading annotation code file: "<<file_func_code<<endl; 
    igzstream infile_code (file_func_code.c_str(), igzstream::in);
    if (!infile_code) {cout<<"error opening annotation file: "<<file_func_code<<endl; return false;}

    while (!safeGetline(infile_code, line).eof()) {
        
        if (line[0] == '#') {
            pch = (char *)line.c_str();
            nch = strchr(pch, '\t');
            n_type = strtol(nch, NULL, 0);
            // cout << "Number of annotation categories" << n_type << endl;            
            mFunc.assign(n_type, 0);
            continue;
        }
        else {
            pch = (char *)line.c_str();
            nch = strchr(pch, '\t');
            func_type.assign(pch, nch-pch);
            func_code = strtol(nch, NULL, 0);
            //cout << func_type << ":" << func_code << endl;
            mapFunc2Code[func_type] = func_code;
        }
    }
    infile_code.close();
    infile_code.clear();
    
    // Load annotation file...
    // cout<<"Reading annotation file: "<<file_anno<<endl; 
    igzstream infile (file_anno.c_str(), igzstream::in);
    if (!infile) {cout<<"error opening annotation file: "<<file_anno<<endl; return false;}
    
    // read function annotation file
    string rs, chr;
    long int b_pos=0;
    size_t snp_i = 0;
    double maf_temp;

    while (!safeGetline(infile, line).eof()) {
        if (line[0] == '#' || (line[0] == 'I' && line[1] == 'D')) {
            continue;
        }
        else {
          if (!indicator_snp[snp_i]) {
          	// SNP is excluded from analysis
            //  pch=(char *)line.c_str();
            //  nch = strchr(pch, '\t');
            //  rs.assign(pch, nch-pch);
              /*if (snpInfo[snp_i].rs_number.compare(rs) != 0) {
                  cerr << "annotation file ID dose not match genotype file ID...\n";
                  exit(-1);
              }*/ // do not check for variant ID
              snp_i++;
              continue;
          }
          else{
            pch=(char *)line.c_str();
            nch = strchr(pch, '\t');
            rs.assign(pch, nch-pch);

            //cout << "snp_i=" << snp_i << ";anno rs=" << rs << "; snpInfo.rs_number=" << snpInfo[snp_i].rs_number << endl;
            /*if (snpInfo[snp_i].rs_number.compare(rs) != 0) {
                cerr << "annotation file ID dose not match genotype file ID...\n";
                exit(-1);
            }*/
            
            pch = (nch == NULL) ? NULL : nch+1;
            nch = strchr(pch, '\t');
            chr.assign(pch, nch-pch);

            pch = (nch == NULL) ? NULL : nch+1;
            nch = strchr(pch, '\t');
            b_pos = strtol(pch, NULL, 0);

            pch = (nch == NULL) ? NULL : nch+1;
            snp_nfunc = 0;
            snpInfo[snp_i].indicator_func.assign(n_type, 0);
            maf_temp = snpInfo[snp_i].maf;
        	if(maf_temp > 0.5) maf_temp = 1.0 - maf_temp;

            //if (snp_i < 5)  cout << rs << ":chr" << chr << ":bp"<< b_pos <<endl;
            if( isalpha(pch[0]) || isdigit(pch[0]) ){
            	//pch[0] is a letter or number
            	while (pch != NULL) {
	                nch = strchr(pch, ',');
	                if (nch == NULL) func_type.assign(pch);
	                else func_type.assign(pch, nch-pch);

	                // consider MAF range
	                // if((func_type.compare("others") == 0) || (func_type.compare("nonsyn") == 0) ) SetMAFCode(maf_temp, func_type);

	                func_code = mapFunc2Code[func_type];
	                //if(snp_i < 10)  cout << func_type << " with code " << func_code << endl;
	                if(!snpInfo[snp_i].indicator_func[func_code])
	                {
	                    snpInfo[snp_i].indicator_func[func_code] = 1;
	                    snp_nfunc++;
	                }
	                pch = (nch == NULL) ? NULL : nch+1;
            	}
        	}
        	else{
        		func_type.assign("NA");
        		// consider MAF range
        		// SetMAFCode(maf_temp, func_type);

        		func_code = mapFunc2Code[func_type];
        		//if(snp_i < 10)  cout << "NA" << " with code " << func_code << endl;
        		if(!snpInfo[snp_i].indicator_func[func_code])
	                {
	                    snpInfo[snp_i].indicator_func[func_code] = 1;
	                    snp_nfunc++;
	                }
        	}
            
            //if ((snp_nfunc > 0) && (snp_nfunc <= n_type))
            if (snp_nfunc == 1)
              {
                  snpInfo[snp_i].weight_i = 1.0 ;// / (double)snp_nfunc;
                  mFunc[func_code]++;
                  // CalcWeight(snpInfo[snp_i].indicator_func, snpInfo[snp_i].weight, snpInfo[snp_i].weight_i);
              }
            else if (snp_nfunc == 0) {
                snpInfo[snp_i].weight_i = 0.0;
                indicator_snp[snp_i] = 0;
                cout << "function annotation is NULL \n ";
            }
            else {cerr << "ERROR: snp_nfunc = " <<snp_nfunc<< " ... \n"; exit(-1);}
            snp_i++;
          }
        }
    }
    cout << "Number of annotation categories: " << n_type << endl;
    cout << "Number of variants per category: "; PrintVector(mFunc);
    //cout << "total snp number = " << snp_i << endl;
    
    infile.close();
    infile.clear();	
    
    return true;
}

//Empty Annotation
bool Empty_anno (vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, size_t &n_type, vector<size_t> &mFunc)
{
    cout << "Empty annotation file, all variants are treated as of one category!" << endl;

    n_type = 1; // all variants are of one annotation
    mFunc.assign(1, 0);

    for(size_t i = 0; i < indicator_snp.size(); i++){
        if(indicator_snp[i] == 0) continue;
        snpInfo[i].indicator_func.assign(n_type, 1);
        snpInfo[i].weight_i = 1.0 ;
        mFunc[0]++;
    }

    cout << "Number of annotation categories: " << n_type << endl;
    cout << "Number of variants per category: "; PrintVector(mFunc);
    
    return true;
}



//Read geno/VCF phenotype file, 
bool ReadFile_pheno (const string &file_pheno, vector<bool> &indicator_idv, vector<double> &pheno, vector<string> &InputSampleID, size_t & ni_total)
{
	indicator_idv.clear();
	pheno.clear();
	
    cout << "open phenotype file ... " << file_pheno << "\n";
    
	igzstream infile (file_pheno.c_str(), igzstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_pheno<<endl; return false;}

	string line;
	char *ch_ptr;
  
	string id;
	
    size_t numPheno=0;

	while (!safeGetline(infile, line).eof()) {

		// first column: sample id
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		id=ch_ptr;
		InputSampleID.push_back(id); //load first column as Sample IDs.
        //if(numPheno < 10) {cout <<"pheno id "<< id << endl;}
		
		// second column: NA or quantitative pheno
		ch_ptr=strtok (NULL, " , \t");
        //if(numPheno < 10) {cout <<"pheno value "<< ch_ptr << endl;}
        
		if (strcmp(ch_ptr, "NA")==0) {
			indicator_idv.push_back(0); 
			pheno.push_back(-9);
		}
        else
        {
            indicator_idv.push_back(1);
            pheno.push_back( atof(ch_ptr) );
        }
					
        numPheno++;
	}
    //cout << "Load numPheno = " << numPheno << "\n";
    ni_total = indicator_idv.size();
 
	infile.close();
	infile.clear();	
	
	return true;
}


//Read .bim file (SNP information)
bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo)
{
	snpInfo.clear();
	
    cout << "Start reading bim file: " << file_bim << "\n";
	ifstream infile (file_bim.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening .bim file: "<<file_bim<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	long int b_pos=0;
	string chr;
	double cM;
	string major;
	string minor;
    
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		chr=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		cM=atof(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		b_pos=atol(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		major=ch_ptr;

        if(rs.compare(".") == 0 || rs.empty()){
                rs = chr + ":" + to_string(b_pos) + ":" + minor + ":" + major;
        }
		
        SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, -9, -9, -9, indicator_func_temp, weight_temp, 0.0};
		snpInfo.push_back(sInfo);
	}
	
	infile.close();
	infile.clear();
    cout << "Success reading bim file.\n";
	return true;
}


//Read geno file to get information of indicator_idv; ni_total
bool getIDVgeno(const string &file_geno, vector<bool> &indicator_idv, size_t & ni_total)
{
    indicator_idv.clear();
    
    igzstream infile (file_geno.c_str(), igzstream::in);
    if (!infile) {cout<<"error! fail to open phenotype file: "<<file_geno<<endl; return false;}

    string line;
    char *pch, *nch=NULL;
    size_t tab_count;

    !safeGetline(infile, line).eof(); // read first header line
        
    for (tab_count=0; pch != NULL; tab_count++) {
        nch=strchr(pch, '\t'); //point to the position of next '\t'
        if (tab_count > 4) {
            indicator_idv.push_back(1);
        }
        pch = (nch == NULL) ? NULL : nch+1;
    }
    
    ni_total = indicator_idv.size();
    cout << "ni_total = " << ni_total << "\n";
 
    infile.close();
    infile.clear(); 
    
    return true;
}

//Read vcf file to get information of indicator_idv; ni_total
bool getIDVvcf(const string &file_vcf, vector<bool> &indicator_idv, size_t & ni_total, string &GTfield)
{
    indicator_idv.clear();
    
    igzstream infile (file_vcf.c_str(), igzstream::in);
    if (!infile) {cout<<"error! fail to open phenotype file: "<<file_vcf<<endl; return false;}

    string line;
    char *pch, *nch=NULL;
    size_t tab_count;

    while(!safeGetline(infile, line).eof()) // read first header line
    {
        if (line[0] == '#') {
            continue;
        }
        else{
            for(tab_count=0; pch != NULL; tab_count++) {
                nch=strchr(pch, '\t'); //point to the position of next '\t'
                if (tab_count >= 9) {
                    indicator_idv.push_back(1);
                }
                pch = (nch == NULL) ? NULL : nch+1;
            }
            break;
        }
    }
    
    ni_total = indicator_idv.size();
    cout << "ni_total = " << ni_total << "\n";
 
    infile.close();
    infile.clear(); 
    
    return true;
}


//Read .fam file
bool ReadFile_fam (const string &file_fam, vector<bool> &indicator_idv, vector<double> &pheno, vector<string> & InputSampleID, size_t &ni_total)
{
	indicator_idv.clear();
	pheno.clear();
	
    cout << "Start reading fam file: " << file_fam <<"\n";
	igzstream infile (file_fam.c_str(), igzstream::in);
	if (!infile) {cout<<"error opening .fam file: "<<file_fam<<endl; return false;}

	string line, id;
	char *ch_ptr;
	double p;
    InputSampleID.clear(); // save sample IDs
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " \t"); // Family ID
		ch_ptr=strtok (NULL, " \t"); //individual or sample ID
		id=ch_ptr;
        InputSampleID.push_back(id);

		ch_ptr=strtok (NULL, " \t"); // paternal id
		ch_ptr=strtok (NULL, " \t"); // maternal id
		ch_ptr=strtok (NULL, " \t"); // sex
		ch_ptr=strtok (NULL, " \t"); // phenotype
		
		if (strcmp(ch_ptr, "NA")==0) {
			indicator_idv.push_back(0); 
			pheno.push_back(-9);
		} else {
			p=atof(ch_ptr);
			if (p==-9) {
				indicator_idv.push_back(0); 
				pheno.push_back(-9);
			}
			else {
				indicator_idv.push_back(1); 
				pheno.push_back(p);
			}
		}
	}
 
    ni_total = indicator_idv.size();
    cout << "ni_total = " << ni_total << "\n";

	infile.close();
	infile.clear();
    cout << "Success reading fam file.\n";
	return true;
}

bool CreatVcfHash(const string &file_vcf, StringIntHash &sampleID2vcfInd, const string &file_sample){
        
    VcfFileReader inFile;
    VcfHeader header;
    
    if(!inFile.open(file_vcf.c_str(), header, file_sample.c_str(), NULL, NULL))
    {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(1);
    }
    
    uint numSample = (uint)header.getNumSamples();
    cout << "numSample = " << numSample << endl;
    String sample_name;
	for (size_t i=0; i<numSample; ++i) {
        sample_name = header.getSampleName(i);
        sampleID2vcfInd.Add(sample_name, i);
	}
    cout << "\n create hash sampleID to vcf index success...\n";
    return true;
}

/* void GetVcfPos(const vector<string> &VcfSampleID, const map<string, size_t> &PhenoID2Ind, vector <size_t> &SampleVcfPos)
{
    size_t yidx;
    string sampleid;
    SampleVcfPos.clear();
    
    for (size_t i=0; i < VcfSampleID.size(); i++) {
        sampleid = VcfSampleID[i];
        if (PhenoID2Ind.count(sampleid) == 0) continue;
        else {
            yidx = PhenoID2Ind.at(sampleid);
            SampleVcfPos.push_back(i);
        }
    }
} */

// Read VCF genotype file, the first time,
bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ns_total, size_t &ni_test, string &GTfield, const map<string, size_t> &PhenoID2Ind, vector<string> &VcfSampleID, vector<size_t> &SampleVcfPos)
{
    if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string
    cout << "Load VCF file genotype field: " << GTfield << endl;
    
    VcfSampleID.clear();
    SampleVcfPos.clear(); // with length = ni_total
    indicator_snp.clear();
    snpInfo.clear();
    ns_test=0; // variable defined in param.h
    
    igzstream infile(file_vcf.c_str(), igzstream::in);
    cout << "open vcf file ...\n";
    if(!infile) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(-1);
    }
        
    long int b_pos = 0; string chr;
    string rs, major, minor, s, pheno_id, line; 
    size_t pheno_index, n_miss, n_0, n_1, n_2, c_idv=0, ctest_idv=0, tab_count;
    double maf, geno, geno_old, cM=-9;
    int flag_poly, GTpos=0, k=0;  // flag polymophysum variant
    char *pch, *p, *nch=NULL, *n;

    gsl_vector *genotype = gsl_vector_alloc(ni_test);
    vector<bool> genotype_miss(ni_test, 0);
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;

    //cout << "PhenoID2Ind.size() = " << PhenoID2Ind.size() << "Before first time load vcf file ... " << endl;

    // cout << "start reading record ... \n";
  while(!safeGetline(infile, line).eof()) {
        if (line[0] == '#') {
           if (strncmp(line.c_str(), "#CHROM", 6) == 0) {
               pch= (char *)line.c_str();
             //parse for individual IDs, save VCFsampleID, create SampleVcfPos
               for (tab_count=0; pch != NULL; tab_count++) {
                   nch=strchr(pch, '\t'); //point to the position of next '\t'
                   if (tab_count>8) {
                       if (nch == NULL) { s.assign( pch );}
                       else s.assign( pch, nch-pch );
                       VcfSampleID.push_back(s);
                       if (PhenoID2Ind.count(s)>0) {
                       		//cout << "id = " << s << "tab_count = " << tab_count << ", ";
                           SampleVcfPos.push_back(tab_count); //record tab_position
                       }
                   }
                   pch = (nch == NULL) ? NULL : nch+1;
               }
               cout << "\n Matched phenotype sample IDs in the VCF file " << SampleVcfPos.size() << "\n";
            }
            continue;
        }
    else{
        c_idv=0; ctest_idv = 0; n_0=0; n_1=0; n_2=0;
        maf=0; n_miss=0; flag_poly=0; geno_old=-9;
        
        pch= (char *)line.c_str();
        
        for (tab_count=0; pch != NULL; tab_count++) {

            nch=strchr(pch, '\t'); //point to the position of next '\t'
            
            if (tab_count<5) {
                if (nch == NULL) { s.assign( pch );}
                else s.assign( pch, nch-pch ); // field string s
                
                switch (tab_count) {
                    case 0:
                        chr = s; break;
                    case 1:
                        b_pos=atol(s.c_str()); break;
                    case 2:
			            rs = s; break;
                    case 3:
                        minor = s; break;
                    case 4:
                        major = s; break;
                    default:
                        break;
                }
                if(rs.compare(".") == 0 || rs.empty()){
                    rs = chr + ":" + to_string(b_pos) + ":" + minor + ":" + major;
                }
                if (setSnps.size()!=0 && setSnps.count(rs)==0) {
                    indicator_snp.push_back(0);
                    pch = (nch == NULL) ? NULL : nch+1;
                    continue;
                }
            }
            
            else if ((tab_count == 6) && (pch[0] == 'F')){
                SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
                snpInfo.push_back(sInfo); //save marker information
                indicator_snp.push_back(0);
                pch = (nch == NULL) ? NULL : nch+1;
                continue;
                //failed filter, continue with next record
            }
            
            else if ((tab_count == 8) && (ns_test == 0))
            {
                // cout << "parse FORMAT field" << endl;
            	// cout << "; ni_test = " <<ni_test << "; n_miss start = " << n_miss << "\n";
                // parse FORMAT field
                if (pch[0] == GTfield[0] && pch[1] == GTfield[1] && ((nch==pch+2)||pch[2]==':') ) {
                    GTpos=0; //GT start in the first position
                    //cout << "GT start in the first position" << endl;
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
                        cerr << "Cannot find" << GTfield << "at marker" << chr << ":" << b_pos << endl;
                        exit(-1);
                    }
                }
            }
            else if ( tab_count == SampleVcfPos[ctest_idv] )
                {
                	//cout << "SampleVcfPos[ctest_idv] = " << SampleVcfPos[ctest_idv] << ", ";
                	//cout << "tab_count = "<< tab_count << ", c_idv = " << c_idv << ", " << "ctest_idv = " << ctest_idv << endl;
                	//cout << pheno_id << " "<< PhenoID2Ind.count(pheno_id) << endl;

                	pheno_id = VcfSampleID[c_idv];            
                	if (PhenoID2Ind.count(pheno_id) > 0){
                		pheno_index = PhenoID2Ind.at(pheno_id);
                	}
                	else {
                		cerr << "phenotype id matched error ... " << endl;
                  	    exit(-1);
                		continue;
                	}

                  if ( !indicator_idv[pheno_index] ) {
                  	   cerr << "error: pheno Ind is 0 ... " << endl;
                  	   exit(-1);
                      //c_idv++; pch = (nch == NULL) ? NULL : nch+1; continue;
                  }
                  else{
                    p = pch; // make p reach to the key index
                    //cout << "GTpos : " << GTpos << endl;
                    //cout << p << endl;

                    if (GTpos>0) {
                        for (int i=0; (i<GTpos) && (p!=NULL); ++i) {
                            n = strchr(p, ':');
                            p = (n == NULL) ? NULL : n+1;
                        }
                    }
                    
                    // pgeno=p;
                    // n = strchr (p, '\t'); // pop out first GC/EC field
                    // cout <<"after parse for EC: " << (n - pgeno) << endl;

                    if (p==NULL) {
                        geno = -9;//missing
                    }
                    //start here
                    else if ( (p[1] == '/') || (p[1] == '|') ) {
                        //read bi-allelic GT
                            if( (p[0]=='.') && (p[2]=='.')){
                                geno = -9;//missing
                            }
                            else if ( (p[0]=='.') && (p[2]!='.')) {
                                geno = (double)(p[2] -'0');
                                if(geno != 1 && geno != 0){ 
                                    geno = -9; // multi-allelic
                                }
                            }
                            else if ((p[0]!='.') && p[2]=='.') {
                                geno = (double)(p[0] -'0');
                                if(geno != 1 && geno != 0){ geno = -9; } // multi-allelic
                            }
                            else {
                                geno = (double)((p[0] - '0') + (p[2]- '0'));
                                if(geno != 1 && geno != 0 && geno != 2){ geno = -9; } // multi-allelic
                            }                   
                        }
                        else if ( GTfield != "GT" ) {
                            //read dosage data
                            if( (p[0]=='.') && ( (p[1] == '\t') || (p[1] == ':') ) ){
                                geno = -9; // missing                           
                            }else if (isdigit(p[0])){
                                geno = strtod(p, NULL);
                                if(geno < 0 || geno > 2) {geno = -9;} // invalid dosage
                            }else{
                                cout << chr + ":" + to_string(b_pos) + ":" + minor + ":" + major << "; Pheno_ID = " << pheno_id << endl;
                                cerr << " has dosage data that is not a digit ... " << endl;
                                exit(-1);
                            }                        
                        }else{
                            geno = -9; // Not in GT format with GTfield="GT"
                        }
                                                // Missing or multi-allelic
                        if(geno == -9){
                            genotype_miss[ctest_idv]=1; 
                            n_miss++; c_idv++; ctest_idv++;
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;
                        }
                        else if( (geno >= 0.0) && (geno <= 2.0)) 
                            {
                                if (geno>=0 && geno<=0.5) {n_0++; maf+=geno;}
                                if (geno>0.5 && geno<1.5) {n_1++; maf+=geno;}
                                if (geno>=1.5 && geno<=2.0) {n_2++; maf+=geno;}
                            }                       
                        else {
                            cout << "ERROR: geno falls outside [0, 2]! " << geno <<";" << pheno_id << endl;
                            exit(-1); 
                        }
                        // end here

                    // cout << "geno = " << geno << endl;
                    // if(n_miss > 600) exit(-1);

                    if (flag_poly==0) {geno_old=geno; flag_poly=2;}
                    if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}

                    gsl_vector_set (genotype, ctest_idv, geno);
                    ctest_idv++;
                    c_idv++;
                  }
                }
                else if(tab_count >= 9){
                  	c_idv++;
                  }
            pch = (nch == NULL) ? NULL : nch+1;
        }

	//cout << "Total sample number " << c_idv << "; analyzed sample number " << ctest_idv << "\n";

        if (ctest_idv != ni_test) {
            cerr << "matched sample size in the genotype file " << ctest_idv << " dose not equal to analyzed sample size " << ni_test << "\n";
            exit(-1);
        }

        if(ni_test > n_miss){
            maf/=2.0*(double)(ni_test-n_miss);
        }else{
            maf = 0.0;
        }
        
        //cout << "maf = " << maf << "; ni_test = " <<ni_test << "; n_miss = " << n_miss << "\n";
        // exit(-1);
        
        SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
        snpInfo.push_back(sInfo); //save marker information
        
        // filter by missing rate
        if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
        // cout << "pass missness criteron...\n";
        
        if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {
            //cout << rs <<":"<<chr<<":"<<b_pos<<":"<<major<<":"<<minor << " filtered by polymorphism \n";
            indicator_snp.push_back(0); 
            continue;
        }
        
        if ( (maf < maf_level || maf > (1.0-maf_level)) && maf_level!=-1 ) {
            //cout << rs <<":"<<chr<<":"<<b_pos<<":"<<major<<":"<<minor << " filtered by MAF cutoff \n";
            indicator_snp.push_back(0); 
            continue;
        }
        //cout << "pass maf criteron...\n";
        
        if (flag_poly!=1) {
            //cout << rs <<":"<<chr<<":"<<b_pos<<":"<<major<<":"<<minor << " filtered by flag_poly \n";
            indicator_snp.push_back(0); 
            continue;}
        // cout << "pass poly criteron...\n";
        
        if (hwe_level!=0) {
            if (CalcHWE(n_0, n_2, n_1)<hwe_level) {
                //cout << rs <<":"<<chr<<":"<<b_pos<<":"<<major<<":"<<minor << " filtered by HWE \n";
                indicator_snp.push_back(0); 
                continue;
            }
        }
        //  cout << "pass hwe criteron...\n";
        
        //filter SNP if it is correlated with W
        for (size_t i=0; i<ni_test; ++i) {
            if (genotype_miss[i]) {
                geno=maf*2.0;
                gsl_vector_set (genotype, i, geno);
            }
        }

        indicator_snp.push_back(1);
        ns_test++;
        //if (ns_test < 5) {sInfo.printMarker(); PrintVector(genotype, 10);}
        }
    }

    ns_total = indicator_snp.size();

    // cout << "genotype vector:\n";
    // PrintVector(genotype, 10);
    //cout << "VCF tab_count = " << tab_count << endl;
    cout << "vcf read first time success ... \n";
    // cout << "analyzed sample size ns_test = " << ns_test << "; loaded sample size ns_total = " << ns_total<<"\n";
     
    gsl_vector_free (genotype);
    infile.clear();
    infile.close();

    //cout << "PhenoID2Ind.size() = " << PhenoID2Ind.size() << "in the end of first vcf file loading... " << endl;
    return true;
}

//Read genotype file, the first time
bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const map<string, size_t> &PhenoID2Ind, vector<SNPINFO> &snpInfo, vector<string> &VcfSampleID, vector<size_t> &SampleVcfPos, const double &maf_level, const double &miss_level, const double &hwe_level, size_t &ns_test, size_t &ns_total, const size_t &ni_test, const size_t &ni_total) 
{
	indicator_snp.clear();
	snpInfo.clear();
    ns_test = 0;
	
	igzstream infile (file_geno.c_str(), igzstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}

	gsl_vector *genotype=gsl_vector_alloc (indicator_idv.size());
	
    char *pch, *nch=NULL;
	long int b_pos=0;
	string s, rs, line, chr, major, minor, pheno_id;
	double cM=-9, maf, geno, geno_old;
	size_t c_idv, ctest_idv, n_miss, n_0, n_1, n_2, tab_count, pheno_index;
	int flag_poly;

    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
    
	while (!safeGetline(infile, line).eof()) {

        pch= (char *)line.c_str();

        // read first header line
        if (strncmp(line.c_str(), "ID", 2) == 0) {
            //parse for individual IDs, save VCFsampleID, create SampleVcfPos
            for (tab_count=0; pch != NULL; tab_count++) {
                nch=strchr(pch, '\t'); //point to the position of next '\t'
                if (tab_count > 4) {
                    if (nch == NULL) { s.assign( pch );}
                    else s.assign( pch, nch-pch );
                    VcfSampleID.push_back(s);
                    if (PhenoID2Ind.count(s)>0) {
                        //cout << "id = " << s << "tab_count = " << tab_count << ", ";
                        SampleVcfPos.push_back(tab_count); //record tab_position
                    }
                }
                pch = (nch == NULL) ? NULL : nch+1;
            }
            cout << "\n Matched phenotype Sample IDs in the genotype file " << SampleVcfPos.size() << "\n";
            continue;
        }else{
            nch=strchr(pch, '\t'); // parse ID first
            if (nch == NULL) { s.assign( pch );}
            else s.assign( pch, nch-pch ); // field string s
            rs = s;

            if (setSnps.size()!=0 && setSnps.count(rs)==0) {
                indicator_snp.push_back(0); continue;
            }

            pch = (nch == NULL) ? NULL : nch+1;

            c_idv=0; ctest_idv = 0; n_0=0; n_1=0; n_2=0;
            maf=0; n_miss=0; flag_poly=0; geno_old=-9;
            vector<bool> genotype_miss(ni_test, 0);
            
            for (tab_count=1; pch != NULL; tab_count++) {

                nch=strchr(pch, '\t'); //point to the position of next '\t'
            
                if (tab_count<5) {
                    if (nch == NULL) { s.assign( pch );}
                    else s.assign( pch, nch-pch ); // field string s

                    switch (tab_count) {
                        case 1:
                            chr = s; 
                            break;
                        case 2:
                            b_pos=atol(s.c_str()); break;
                        case 3:
                            minor = s; break;
                        case 4:
                            major = s; break;
                        default:
                            break;
                    }
                    pch = (nch == NULL) ? NULL : nch+1;  
                }
                if(rs.compare(".") == 0 || rs.empty()){
                    rs = chr + ":" + to_string(b_pos) + ":" + minor + ":" + major;
                }
                else if ( tab_count == SampleVcfPos[ctest_idv] )
                {
                    pheno_id = VcfSampleID[c_idv];            
                    pheno_index = PhenoID2Ind.at(pheno_id); 
                        //should exist in the PhenoID2Ind                    

                  if ( !indicator_idv[pheno_index] ) {
                       cout << "phenotype of "<<rs<<"with id "<<pheno_id<<" is not analyzed."<< endl;
                       pch = (nch == NULL) ? NULL : nch+1;
                       c_idv++;  
                       continue;
                  } 
                  else{
                    // read genotype value
                    if (pch == NULL) {
                        geno = -9;//missing
                        genotype_miss[ctest_idv]=1; n_miss++; c_idv++; ctest_idv++; 
                        pch = (nch == NULL) ? NULL : nch+1;
                        continue;
                    }
                    else {
                        //read dosage data 
                        if( ((pch[0]=='N') && (pch[1] == 'A')) || ((pch[0]=='.') && (pch[1] == '\t'))  ){
                            geno = -9;
                            genotype_miss[ctest_idv]=1; 
                            n_miss++; c_idv++; ctest_idv++;
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;                           
                        }else{
                            if (nch == NULL) { s.assign( pch );}
                            else s.assign( pch, nch-pch ); // field string s
                            geno = atof(s.c_str());
                        }  
                    }
                    // cout << "geno = " << geno << endl;
                    // if(n_miss > 600) exit(-1);

                    if (geno>=0 && geno<=0.5) {n_0++; maf+=geno;}
                    if (geno>0.5 && geno<1.5) {n_1++; maf+=geno;}
                    if (geno>=1.5 && geno<=2.0) {n_2++; maf+=geno;}

                    gsl_vector_set (genotype, ctest_idv, geno);

                    if (flag_poly==0) {geno_old=geno; flag_poly=2;}
                    if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}

                    ctest_idv++;
                    c_idv++;
                  }
                }
            else if(tab_count >= 5){
                    c_idv++;
            }
            pch = (nch == NULL) ? NULL : nch+1;
        }

		maf/=2.0*(double)(ni_test-n_miss);	
		
		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
		snpInfo.push_back(sInfo);
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
		
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
		
        // replace missing genotypes with sample mean.
		for (size_t i=0; i<ni_test; ++i) {
			if ( genotype_miss[i] )
            {
                geno=maf*2.0; 
                gsl_vector_set (genotype, i, geno);
            }
		}
		
		indicator_snp.push_back(1); 
		ns_test++;
      }
	}
	
    ns_total = indicator_snp.size();
	gsl_vector_free (genotype);
	
	infile.close();
	infile.clear();	
	
	return true;
}


//Read bed file, the first time
bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, const map<string, size_t> &PhenoID2Ind, const size_t &ni_test, const size_t &ni_total, const double &maf_level, const double &miss_level, const double &hwe_level, size_t &ns_test, size_t &ns_total) 
{
	indicator_snp.clear();
	ns_total=snpInfo.size();
    ns_test=0;

    gsl_vector *genotype=gsl_vector_alloc (ni_test);
    gsl_vector *genotype_miss=gsl_vector_alloc (ni_test);
	
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
    
	double geno, geno_old, maf;
	size_t c_idv=0, n_miss, n_0, n_1, n_2, c, n_bit;
	char ch[1];
	bitset<8> b;
    int flag_poly;
  	
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}
    
	//ignore the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}

	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
		infile.seekg(t*n_bit+3); //n_bit, and 3 is the number of magic numbers
		
		if (setSnps.size()!=0 && setSnps.count(snpInfo[t].rs_number)==0) {
			snpInfo[t].n_miss=-9;
			snpInfo[t].missingness=-9;
			snpInfo[t].maf=-9;
			indicator_snp.push_back(0);
			continue;
		}
        
		//read genotypes
		c=0; maf=0.0; n_miss=0; n_0=0; n_1=0; n_2=0;
        flag_poly=0; geno_old=-9;
		c_idv=0; 

        gsl_vector_set_zero(genotype_miss);

		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {
                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c==ni_total) {break;}
                c++;
				if (indicator_idv[c]==0) {continue;} //skip the sample
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); maf+=2.0; n_2++;}
					else {gsl_vector_set(genotype, c_idv, 1.0); maf+=1.0; n_1++;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); maf+=0.0; n_0++;}
					else {
                        gsl_vector_set(genotype, c_idv, -9.0);
                        gsl_vector_set(genotype_miss, c_idv, 1.0); n_miss++; 
                    }
				}

                geno = gsl_vector_get(genotype, c_idv);
                if (flag_poly==0) {geno_old=geno; flag_poly=2;}
                if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
                
				c_idv++;
			}
		}
		maf/=2.0*(double)(ni_test-n_miss);
		
		snpInfo[t].n_miss=n_miss;
		snpInfo[t].missingness=(double)n_miss/(double)ni_test;
		snpInfo[t].maf=maf;

		// Apply filter for SNPs
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

        if (flag_poly!=1) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
		
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {
                geno=maf*2.0; 
                gsl_vector_set (genotype, i, geno);
            }
		}
		
		indicator_snp.push_back(1);
		ns_test++;
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
    
	infile.close();
	infile.clear();
	
	return true;
}


void ReadFile_kin (const string &file_kin, vector<bool> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G)
{
	igzstream infile (file_kin.c_str(), igzstream::in);
    //	ifstream infile (file_kin.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open kinship file: "<<file_kin<<endl; error=true; return;}
	
	size_t ni_total=indicator_idv.size();
	
	gsl_matrix_set_zero (G);
	
	string line;
	char *ch_ptr;
	double d;
	
	if (k_mode==1) {
		size_t i_test=0, i_total=0, j_test=0, j_total=0;
		while (getline(infile, line)) {
			if (i_total==ni_total) {cout<<"error! number of rows in the kinship file is larger than the number of phentypes."<<endl; error=true;}
			
			if (indicator_idv[i_total]==0) {i_total++; continue;}
			
			j_total=0; j_test=0;
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			while (ch_ptr!=NULL) {
				if (j_total==ni_total) {cout<<"error! number of columns in the kinship file is larger than the number of phentypes for row = "<<i_total<<endl; error=true;}
				
				d=atof(ch_ptr);
				if (indicator_idv[j_total]==1) {gsl_matrix_set (G, i_test, j_test, d); j_test++;}
				j_total++;
				
				ch_ptr=strtok (NULL, " , \t");
			}
			if (j_total!=ni_total) {cout<<"error! number of columns in the kinship file do not match the number of phentypes for row = "<<i_total<<endl; error=true;}
			i_total++; i_test++;
		}
		if (i_total!=ni_total) {cout<<"error! number of rows in the kinship file do not match the number of phentypes."<<endl; error=true;}
	}
	else {
		map<size_t, size_t> mapID2ID;
		size_t c=0;
		for (size_t i=0; i<indicator_idv.size(); i++) {
			if (indicator_idv[i]==1) {mapID2ID[i]=c; c++;}
		}
		
		string id1, id2;
		double Cov_d;
		size_t n_id1, n_id2;
		
		while (getline(infile, line)) {
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			id1=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			id2=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			d=atof(ch_ptr);
			if (mapID2num.count(id1)==0 || mapID2num.count(id2)==0) {continue;}
			if (indicator_idv[mapID2num[id1]]==0 || indicator_idv[mapID2num[id2]]==0) {continue;}
			
			n_id1=mapID2ID[mapID2num[id1]];
			n_id2=mapID2ID[mapID2num[id2]];
			
			Cov_d=gsl_matrix_get(G, n_id1, n_id2);
			if (Cov_d!=0 && Cov_d!=d) {cout<<"error! redundant and unequal terms in the kinship file, for id1 = "<<id1<<" and id2 = "<<id2<<endl;}
			else {
				gsl_matrix_set(G, n_id1, n_id2, d);
				gsl_matrix_set(G, n_id2, n_id1, d);
			}
		}
	}
	
	infile.close();
	infile.clear();
	
	return;
}


//read genotype text file and calculate kinship matrix
bool GenoKin (const string &file_geno, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line, pheno_id, s;
	char *pch, *nch=NULL;
	size_t n_miss, c_idv, ctest_idv, c_snp=0, ns_test=0, tab_count, pheno_index;
	double d, geno_mean, geno_var, geno;
	size_t ni_test=matrix_kin->size1;
	gsl_vector *geno_vec=gsl_vector_alloc (ni_test);

    while(!safeGetline(infile, line).eof()){

        if (c_snp%display_pace==0 || c_snp==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", c_snp, indicator_snp.size()-1);}

        pch= (char *)line.c_str();

        if ( (strncmp(line.c_str(), "ID", 2) == 0) ) {continue;} // skip header 
        else{
            if (indicator_snp[c_snp]==0) {c_snp++; continue;} // skip unanalyzed snp            
            c_idv=0; ctest_idv = 0; geno_mean = 0.0; n_miss = 0; geno_var = 0.0;
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
                        geno = -9;//missing
                        gsl_vector_set (geno_vec, ctest_idv, -9.0);
                        n_miss++; c_idv++; ctest_idv++; 
                        pch = (nch == NULL) ? NULL : nch+1;
                        continue;
                    }
                    else {
                        //read dosage data 
                        if( ((pch[0]=='N') && (pch[1] == 'A')) || ((pch[0]=='.') && (pch[1] == '\t'))){
                            geno = -9;
                            gsl_vector_set (geno_vec, ctest_idv, -9.0);                       
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
                        gsl_vector_set (geno_vec, ctest_idv, geno);
                        geno_mean += geno;
                    }else{
                        gsl_vector_set (geno_vec, ctest_idv, -9.0);
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
		
		geno_mean/=(double)(ni_test-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_test;
		geno_var-=geno_mean*geno_mean;
        //		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (geno_vec, i)==-9.0) {gsl_vector_set(geno_vec, i, geno_mean);}
		}
		
		gsl_vector_add_constant (geno_vec, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno_vec, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno_vec, matrix_kin);}
			else {
                cout<<"Unknown kinship mode."<<endl;
                exit(-1);
            }
		}
		ns_test++;
        c_snp++;
    }
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_test; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno_vec);
	
	infile.close();
	infile.clear();
	
	return true;
}


bool PlinkKin (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
    
	char ch[1];
	bitset<8> b;
	
	size_t n_miss, ci_total, ci_test;
	double d, geno_mean, geno_var;
	
    size_t ni_total = indicator_idv.size();
	size_t ni_test=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_test);
    
	size_t ns_test=0;
	size_t n_bit;
	
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }
    
	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		geno_mean=0.0;	n_miss=0; ci_total=0; geno_var=0.0; ci_test = 0;
		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==ni_total) {break;}
                if (indicator_idv[ci_total] == 0) {ci_total++; continue;}

				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(geno, ci_test, 2.0); geno_mean+=2.0; geno_var+=4.0; }
					else {gsl_vector_set(geno, ci_test, 1.0); geno_mean+=1.0; geno_var+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(geno, ci_test, 0.0); }
					else {gsl_vector_set(geno, ci_test, -9.0); n_miss++; }
				}
                ci_test++;
				ci_total++;
			}
		}
        
		geno_mean/=(double)(ni_test-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_test;
		geno_var-=geno_mean*geno_mean;
        //		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_test; ++i) {
			d=gsl_vector_get(geno,i);
			if (d==-9.0) {gsl_vector_set(geno, i, geno_mean);}
		}
		
		gsl_vector_add_constant (geno, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);}
			else {
                cout<<"Unknown kinship mode."<<endl;
                exit(1);
            }
		}
		
		ns_test++;
    }
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_test; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno);
	
	infile.close();
	infile.clear();

	return true;
}

//read VCF file for the 2nd time and calculate kinship matrix ** NEED to be rewritten **
bool VCFKin (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin, string &GTfield, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID)
{
    if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string

    igzstream infile (file_vcf.c_str(), igzstream::in);
    if (!infile) {cout<<"error reading vcf genotype file:"<<file_vcf<<endl; exit(-1);}
    
    double geno, geno_mean, geno_var, d;
    size_t n_miss, c_idv=0, c_snp=0, ctest_idv = 0;
    
    char *pch, *p, *nch=NULL, *n;
    size_t tab_count, pheno_index;
    int GTpos=0, k=0;
    string line, pheno_id;
    
    size_t ni_test=matrix_kin->size1;
    size_t ns_test=0;
    gsl_vector *geno_vec=gsl_vector_alloc (ni_test);
    
    
    while(!safeGetline(infile, line).eof())
    {
        if (c_snp%display_pace==0 || c_snp==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", c_snp, indicator_snp.size()-1);}

        if (line[0] == '#') {
            continue; //skip header
        }
        else {
            if (!indicator_snp[c_snp]) {c_snp++; continue;}
            c_idv=0; //increase to the total individuals ni_total
            geno_mean=0.0; n_miss=0; geno_var=0.0; ctest_idv = 0;
        
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
                        }
                        else if ( (p[1] == '/') || (p[1] == '|') ) {
                        //read bi-allelic GT
                            if( (p[0]=='.') && (p[2]=='.')){
                                geno = -9;//missing
                            }
                            else if ( (p[0]=='.') && (p[2]!='.')) {
                                geno = (double)(p[2] -'0');
                            }
                            else if ((p[0]!='.') && p[2]=='.') {
                                geno = (double)(p[0] -'0');
                            }
                            else geno = (double)((p[0] - '0') + (p[2]- '0'));
                        }
                        else if ( GTfield != "GT" ) {
                            //read dosage data
                            if( (p[0]=='.') && ( (p[1] == '\t') || (p[1] == ':') ) ){
                                geno = -9;                              
                            }else if (isdigit(p[0])){
                                geno = strtod(p, NULL);
                            }else{
                                cout << "Pheno_ID = " << pheno_id  << endl;
                                cerr << " has dosage data that is not a digit ... " << endl;
                                exit(-1);
                            }                        
                        }else{
                            geno = -9;//missing
                        }

                        if(geno == -9){ // missing 
                            gsl_vector_set (geno_vec, ctest_idv, geno);
                            n_miss++; c_idv++; ctest_idv++; 
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;
                        }

                        gsl_vector_set (geno_vec, ctest_idv, geno);
                        if( (geno >= 0.0) && (geno <= 2.0)) {geno_mean += geno;}
                        ctest_idv++; // increase analyzed phenotype #
                        c_idv++;
                    }
                }
                else if ( tab_count >= 9 )
                {
                    c_idv++;
                }
                pch = (nch == NULL) ? NULL : nch+1;
            }        
        geno_mean/=(double)(ni_test-n_miss);
        geno_var+=geno_mean*geno_mean*(double)n_miss;
        geno_var/=(double)ni_test;
        geno_var-=geno_mean*geno_mean;
        //      geno_var=geno_mean*(1-geno_mean*0.5);
        
        for (size_t i=0; i<ni_test; ++i) {
            if ( gsl_vector_get (geno_vec, i) == -9.0 ) 
            {
                gsl_vector_set(geno_vec, i, geno_mean);
            }
        }
        
        gsl_vector_add_constant (geno_vec, -1.0*geno_mean);
        
        if (geno_var!=0) {
            if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno_vec, matrix_kin);}
            else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno_vec, matrix_kin);}
            else {cout<<"Unknown kinship mode."<<endl;}
        }
        
        ns_test++;
      }//end of ifelse
    } // end of while
    cout<<endl;
    
    gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
    
    for (size_t i=0; i<ni_test; ++i) {
        for (size_t j=0; j<i; ++j) {
            d=gsl_matrix_get (matrix_kin, j, i);
            gsl_matrix_set (matrix_kin, i, j, d);
        }
    }
    
    gsl_vector_free (geno_vec);
    
    infile.close();
    infile.clear();
    
    return true;
}

//Read VCF genotype file, the second time, recode genotype and calculate K
bool ReadFile_vcf (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** X, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, string &GTfield, vector<double> &SNPmean, vector <size_t> &CompBuffSizeVec, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID, bool Compress_Flag)
{
    if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string

    // Open the VCF file.
    igzstream infile(file_vcf.c_str(), igzstream::in);
    cout << "open vcf file second time ...\n";
    if(!infile) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(-1);
    }
    
    if (calc_K==true) {gsl_matrix_set_zero (K);}
    
    gsl_vector *genotype=gsl_vector_alloc (ni_test);
    uchar *geno_uchar = new uchar[ni_test];
    
    size_t sourceBufferSize = (ni_test) * sizeof(uchar);
    const size_t BufferSize = (size_t)(compressBound(sourceBufferSize));
    uchar * TempCompBuffer = (uchar*)malloc(BufferSize);
    uchar * TempBuffer = (uchar*)malloc(sourceBufferSize);
    size_t compressedBufferSize = BufferSize;
    //cout << "Source Buffer Size = " << sourceBufferSize << "; Comp Buffer Bound = " << BufferSize  << endl;
    CompBuffSizeVec.clear();
    SNPmean.clear();

    double geno, geno_mean, vtx;
    size_t n_miss, c_idv=0, c_snp=0, ctest_snp = 0, ctest_idv=0;
    int result;
    
    char *pch, *p, *nch=NULL, *n;
    size_t tab_count;
    int GTpos=0, k=0;
    string line, pheno_id;
    size_t pheno_index;
    
    //cout << "PhenoID2Ind.size() = " << PhenoID2Ind.size() << " before second vcf file loading... " << endl;

    while(!safeGetline(infile, line).eof())
    {
        if (line[0] == '#') {
            continue; //skip header
        }
        else {
            if (!indicator_snp[c_snp]) {c_snp++; continue;}
            c_idv=0; //increase to the total individuals ni_total
            ctest_idv=0; // increase to the total analyzed individuals
            geno_mean=0.0; n_miss=0;
            vector<bool> genotype_miss(ni_test, 0);
        
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
                        }
                        else if ( (p[1] == '/') || (p[1] == '|') ) {
                        //read bi-allelic GT
                        	if( (p[0]=='.') && (p[2]=='.')){
                        		geno = -9;//missing
                        	}
                        	else if ( (p[0]=='.') && (p[2]!='.')) {
                            	geno = (double)(p[2] -'0');
                                if(geno != 1 && geno != 0){ 
                                    geno = -9; // multi-allelic
                                }
                        	}
                        	else if ((p[0]!='.') && p[2]=='.') {
                           		geno = (double)(p[0] -'0');
                                if(geno != 1 && geno != 0){ geno = -9; } // multi-allelic
                        	}
                        	else {
                                geno = (double)((p[0] - '0') + (p[2]- '0'));
                                if(geno != 1 && geno != 0 && geno != 2){ geno = -9; } // multi-allelic
                            }                   
                    	}
                    	else if (GTfield != "GT") {
                        	//read dosage data
                        	if( (p[0]=='.') && ( (p[1] == '\t') || (p[1] == ':') ) ){
                        		geno = -9; // missing                      		
                        	}else if (isdigit(p[0])){
                        		geno = strtod(p, NULL);
                                if(geno < 0 || geno > 2) {geno = -9;} // invalid dosage
                        	}else{
                                cout << "Pheno_ID = " << pheno_id << endl;
                        		cerr << " has dosage data that is not a digit ... " << endl;
                        		exit(-1);
                        	}                        
                    	}

                        // Missing or multi-allelic
                        if(geno == -9){
                            genotype_miss[ctest_idv]=1; 
                            n_miss++; c_idv++; ctest_idv++;
                            pch = (nch == NULL) ? NULL : nch+1;
                            continue;
                        }
                        else if( (geno >= 0.0) && (geno <= 2.0)) 
                            {geno_mean += geno;}                       
                        else {
                            cout << "ERROR: geno falls outside [0, 2]! " << geno <<";" << pheno_id << endl;
                            exit(-1); 
                        }

                        gsl_vector_set (genotype, ctest_idv, geno);
                        ctest_idv++; // increase analyzed phenotype #
                        c_idv++;
                    }
                }else if (tab_count >= 9 ){
                	c_idv++; 
                }
                pch = (nch == NULL) ? NULL : nch+1;
            } // for tab_count
        if(ni_test <= n_miss){
            geno_mean = 0.0;
        }else{
            geno_mean/=(double)(ni_test-n_miss);
            if(geno_mean < 0.0 || geno_mean > 2.0){
                cout << "ERROR: geno_mean falls outside [0, 2]! \n ";
                exit(-1); 
            }else{
                SNPmean.push_back(geno_mean);
            }
        }
        // cout << "geno_mean = " << geno_mean << endl; 

        
        for (size_t i=0; i < ni_test; ++i) {
                if (genotype_miss[i]) {geno=geno_mean; gsl_vector_set (genotype, i, geno);}
                // do not center genotype data in UCHAR**
                else { geno = gsl_vector_get (genotype, i);}
                geno_uchar[i] = DoubleToUchar(geno);
                //UtX[ctest_snp][i] = DoubleToUchar(geno);
                //if (ctest_snp==0 && i < 10) cout << geno << ":" << (int)geno_uchar[i] << ", ";
            }
        gsl_vector_add_constant(genotype, -geno_mean); // center genotype gsl_vector here
            
            if (Compress_Flag) {
                compressedBufferSize = BufferSize;
                result = compress(TempCompBuffer, &compressedBufferSize, geno_uchar, sourceBufferSize);
                if (result != Z_OK) {
                    zerr(result);
                    exit(-1);
                }
                else {
                    X[ctest_snp] = (uchar*)malloc(compressedBufferSize);
                    memcpy(X[ctest_snp], TempCompBuffer, compressedBufferSize);
                    CompBuffSizeVec.push_back(compressedBufferSize);
                    
                    // UnCompBufferSize=sourceBufferSize;
                    //  result = uncompress(TempBuffer, &UnCompBufferSize, UtX[c_snp],compressedBufferSize);
                    //  if(c_snp < 10)  {
                    //    zerr(result);
                    //cout << "uncompressed buffer size = " << UnCompBufferSize << endl;
                    //  PrintVector(TempBuffer, 10);
                    // }
                    // cout << "compressed Buffer size = " << compressedBufferSize << endl;
                }
            }
            else {
                X[ctest_snp] = (uchar*)malloc(sourceBufferSize);
                memcpy(X[ctest_snp], geno_uchar, ni_test);
            }            
            
            //JY add
            /*gsl_blas_ddot(genotype, genotype, &vtx);
            if(vtx < 0.00000001)
            {cout << "snp has x'x = " << setprecision(9) << vtx << endl;} */
            
            if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
            
            c_snp++;
            ctest_snp++;
        }
    }
    // cout << "ctest_snp = " << c_snp << "; ns_test = " << ns_test << endl;
     //cout << "SNPmean size = " << SNPmean.size() << endl;

    if (calc_K==true) {
        gsl_matrix_scale (K, 1.0/(double)ns_test);
        
        for (size_t i=0; i<genotype->size; ++i) {
            for (size_t j=0; j<i; ++j) {
                geno=gsl_matrix_get (K, j, i);
                gsl_matrix_set (K, i, j, geno);
            }
        }
    }
    
    free(TempBuffer);
    free(TempCompBuffer);
    gsl_vector_free(genotype);
    delete [] geno_uchar;
    infile.clear();
    infile.close();
    
    //cout << "PhenoID2Ind.size() = " << PhenoID2Ind.size() << " after second vcf file loading... " << endl;
    cout << "read vcf file second time success ... \n" ;
    return true;
}


//Read bimbam genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_geno (const string &file_geno, const vector<bool> &indicator_idv, const vector<bool> &indicator_snp, uchar **X, gsl_matrix *K, const bool calc_K, const size_t ni_test, vector<double> &SNPmean, vector <size_t> &CompBuffSizeVec, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID, bool Compress_Flag)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line, pheno_id, s;
	char *pch, *nch=NULL;
    double geno, geno_mean;
    size_t n_miss, c_idv, ctest_idv, c_snp=0, ctest_snp=0, tab_count, pheno_index;
    int result;
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
	gsl_vector *genotype=gsl_vector_alloc (ni_test);
    uchar *geno_uchar = new uchar[ni_test];

    size_t sourceBufferSize = (ni_test) * sizeof(uchar);
    const size_t BufferSize = (size_t)(compressBound(sourceBufferSize));
    uchar * TempCompBuffer = (uchar*)malloc(BufferSize);
    uchar * TempBuffer = (uchar*)malloc(sourceBufferSize);
    size_t compressedBufferSize = BufferSize;
    //cout << "Source Buffer Size = " << sourceBufferSize << "; Comp Buffer Bound = " << BufferSize  << endl;
    CompBuffSizeVec.clear();
    SNPmean.clear();

	
    while(!safeGetline(infile, line).eof()){

        pch= (char *)line.c_str();

        if ( (strncmp(line.c_str(), "ID", 2) == 0) ) {continue;} // skip header 
        else{
		  if (indicator_snp[c_snp]==0) {c_snp++; continue;}

            c_idv=0; ctest_idv = 0; geno_mean = 0.0; n_miss = 0;
            vector<bool> genotype_miss(ni_test, 0);

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
                        geno = -9;//missing
                        genotype_miss[ctest_idv]=1; n_miss++; c_idv++; ctest_idv++; 
                        pch = (nch == NULL) ? NULL : nch+1;
                        continue;
                    }
                    else {
                        //read dosage data 
                        if( ((pch[0]=='N') && (pch[1] == 'A')) || ((pch[0]=='.') && (pch[1] == '\t'))){
                            geno = -9;
                            genotype_miss[ctest_idv]=1; 
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
                        gsl_vector_set (genotype, ctest_idv, geno);
                        geno_mean += geno;
                    }else{
                        gsl_vector_set (genotype, ctest_idv, -9.0);
                        genotype_miss[ctest_idv]=1; 
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

        geno_mean/=(double)(ni_test-n_miss);
        SNPmean.push_back(geno_mean);

        for (size_t i=0; i < ni_test; ++i) {
                if (genotype_miss[i]) {
                    gsl_vector_set (genotype, i, geno_mean);
                    geno=geno_mean; 
                }
                else { geno = gsl_vector_get (genotype, i);}
                geno_uchar[i] = DoubleToUchar(geno);
            }
        gsl_vector_add_constant(genotype, -geno_mean); // center genotype gsl_vector here

        if (Compress_Flag) {
                compressedBufferSize = BufferSize;
                result = compress(TempCompBuffer, &compressedBufferSize, geno_uchar, sourceBufferSize);
                if (result != Z_OK) {
                    zerr(result);
                    exit(-1);
                }
                else {
                    X[ctest_snp] = (uchar*)malloc(compressedBufferSize);
                    memcpy(X[ctest_snp], TempCompBuffer, compressedBufferSize);
                    CompBuffSizeVec.push_back(compressedBufferSize);
                }
            }
        else {
                X[ctest_snp] = (uchar*)malloc(sourceBufferSize);
                memcpy(X[ctest_snp], geno_uchar, ni_test);
            }

		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
        ctest_snp++;
      }
	}

    if(c_snp != indicator_snp.size() ){cerr << "compressed variant number dose not equal to analyzed variant number!" << endl; exit(-1);}
	
    if (calc_K==true) {
        gsl_matrix_scale (K, 1.0/(double)ctest_snp);
        
        for (size_t i=0; i<genotype->size; ++i) {
            for (size_t j=0; j<i; ++j) {
                geno=gsl_matrix_get (K, j, i);
                gsl_matrix_set (K, i, j, geno);
            }
        }
    }

    free(TempBuffer);
    free(TempCompBuffer);
	gsl_vector_free (genotype);	
    delete [] geno_uchar;
    
	infile.clear();
	infile.close();
	
	return true;
}



//Read BED genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_bed (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar **X, gsl_matrix *K, const bool calc_K, const size_t ni_test, const size_t ns_test, const size_t ni_total, const size_t ns_total, vector<double> &SNPmean, vector <size_t> &CompBuffSizeVec, bool Compress_Flag)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
	
	char ch[1];
	bitset<8> b;
    double geno, geno_mean;
	size_t n_bit, n_miss, c_idv=0, c_snp=0, c=0;
    int result;

    CompBuffSizeVec.clear();
    SNPmean.clear();

	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}
	
	//first three magic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
    //cout << "ni_test = " << ni_test << endl;
	gsl_vector *genotype = gsl_vector_alloc (ni_test);
    uchar *geno_uchar = new uchar[ni_test];
    size_t sourceBufferSize = (ni_test) * sizeof(uchar);
   // size_t UnCompBufferSize=sourceBufferSize;
    
    const size_t BufferSize = (size_t)(compressBound(sourceBufferSize));
    uchar * TempCompBuffer = (uchar*)malloc(BufferSize);
    uchar * TempBuffer = (uchar*)malloc(sourceBufferSize);

    size_t compressedBufferSize = BufferSize;
    //cout << "Source Buffer Size = " << sourceBufferSize << "; Comp Buffer Bound = " << BufferSize  << endl;

	//start reading genotypes
	for (size_t t=0; t<ns_total; ++t) {
		if (indicator_snp[t]==0) {continue;}
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes for the t_th snp
		c_idv=0; geno_mean=0.0; n_miss=0; c=0;
		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c == (size_t)ni_total) {break;}
				if (indicator_idv[c]==0) {c++; continue;}
				c++;
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); geno_mean+=2.0;}
					else {gsl_vector_set(genotype, c_idv, 1.0); geno_mean+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); geno_mean+=0.0;}
					else {gsl_vector_set(genotype, c_idv, -9.0); n_miss++;}
				}
				c_idv++;
			}
		}
        if (n_miss > 0) cout << "n_miss = " << n_miss << endl;
		if(c_idv != (size_t)ni_test) cout << "# of readed individuals not equal to ni_test \n";
        
		geno_mean/=(double)(ni_test-n_miss);
        SNPmean.push_back(geno_mean);

        //if(geno_mean < 0.00000001) cout << "SNP_" << c_snp << "has geno_mean =" << geno_mean << endl;
        
		for (size_t i=0; i<ni_test; ++i) {
			geno=gsl_vector_get (genotype, i);
			if (geno==-9.0) {geno=geno_mean; gsl_vector_set (genotype, i, geno);}
            geno_uchar[i] = DoubleToUchar(geno);
		}
        gsl_vector_add_constant(genotype, -geno_mean); // center genotypes
        
        if (Compress_Flag) {
            compressedBufferSize = BufferSize;
            result = compress(TempCompBuffer, &compressedBufferSize, geno_uchar, sourceBufferSize);
            if (result != Z_OK) {
                zerr(result);
                exit(-1);
            }
            else {
                X[c_snp] = (uchar*)malloc(compressedBufferSize);
                memcpy(X[c_snp], TempCompBuffer, compressedBufferSize);
                CompBuffSizeVec.push_back(compressedBufferSize);
            }
        }
        else{
            X[c_snp] = (uchar*)malloc(sourceBufferSize);
            memcpy(X[c_snp], geno_uchar, sourceBufferSize);
        }
		
		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
	}
    //cout << "compressed Buffer size = " << compressedBufferSize << endl;
    //cout << "CompBuffSizeVec length = " << CompBuffSizeVec.size() << endl;
    
	if(c_snp != ns_test) cout <<"# of readed SNP not equal to ns_test \n";
	
	if (calc_K==true) {
		gsl_matrix_scale (K, 1.0/(double)ns_test);
		
		for (size_t i=0; i<genotype->size; ++i) {
			for (size_t j=0; j<i; ++j) {
				geno=gsl_matrix_get (K, j, i);
				gsl_matrix_set (K, i, j, geno);
			}
		}
	}
	
    free(TempBuffer);
    free(TempCompBuffer);
	gsl_vector_free (genotype);
    delete [] geno_uchar;
    
	infile.clear();
	infile.close();
	
	return true;
}


bool CountFileLines (const string &file_input, size_t &n_lines)
{
	igzstream infile (file_input.c_str(), igzstream::in);
	//ifstream infile (file_input.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open file: "<<file_input<<endl; return false;}
    
	n_lines=count(istreambuf_iterator<char>(infile), istreambuf_iterator<char>(), '\n');
	infile.seekg (0, ios::beg);
	
	return true;
}

void WriteMatrix(const gsl_matrix * X, const string file_str){
    //string file_str = "./output/"+file_out;
    //file_str += filename;

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for(size_t i=0; i<X->size1; ++i){
        for(size_t j = 0; j < X->size2; ++j){
            outfile << scientific << setprecision(6) << gsl_matrix_get(X, i, j) << " " ;
        }
        outfile << endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
} //write gsl_matrix X with filename = ***.txt

void WriteVector(const gsl_vector * X, const string file_str){
    // string file_str = "./output/"+file_out;
    // file_str += filename;

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for(size_t i=0; i<X->size; ++i){
        outfile << scientific << setprecision(6) << gsl_vector_get(X, i)<< endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
}

//write gsl_vector X with filename = ***.txt


