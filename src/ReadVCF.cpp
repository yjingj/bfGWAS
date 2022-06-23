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


#include "ReadVCF.h"
#include "bvsrm.h"

//parse for individual IDs and number
/*int parseInds(char* line, std::set<std::string>& idset, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;
    
    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
        nch = strchr(pch, '\t');
        if ( i >= startIdx ) {
            std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
            if ( idset.empty() || (idset.find(id) != idset.end()) ) {
                icols.push_back(i - startIdx);
                inds.push_back(id);
                ++j;
            }
        }
        pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    //notice("icols.size() = %d",(int)icols.size());
    return j;
}*/

void CreateUcharTable(vector<pair<int, double> > &UcharTable)
{
    UcharTable.clear();
    double second_val;
    for (int i=0; i <= UCHAR_MAX; i++) {
        second_val = ((double)i) * 0.01;
        UcharTable.push_back(make_pair(i, second_val));
        //cout << (uchar)i << ", " << i << ", " << second_val << endl;
    }
}


double getDoubleDosageFromRecord(VcfRecord& record, const uint smNum)
{
    //Read EC, dosage data from vcf in hex string
    double dosage = -9.0;
    
    const std::string* ecStrPtr = record.getGenotypeInfo().getString("EC", smNum);
    if(ecStrPtr != NULL)
    {
        //cout << *ecStrPtr << " ";
        dosage = strtod(ecStrPtr->c_str(), NULL);
    }
    
    //read GT
    else{
        int alleleIndex;
        dosage = 0.0;
        for(int i = 0; i < 2; i++)
        {
            alleleIndex = record.getGT(smNum, i);
            if(alleleIndex == VcfGenotypeSample::MISSING_GT)
            { dosage = -9.0; break;}
            else if(alleleIndex > 0)
            {
                dosage += 1.0;
            }
        }
    }

    if (dosage < 0.0 || dosage > 2.0) {
        dosage = -9.0;
    }
    
    return dosage;
    
}


uchar getUcharDosageFromRecord(VcfRecord& record, const uint smNum)
{
    //Read EC, dosage data from vcf in hex string
    uchar c;
    int intc;
    float dosage;
    
    const std::string* ecStrPtr = record.getGenotypeInfo().getString("EC", smNum);
    if(ecStrPtr != NULL)
    {
        dosage = strtof(ecStrPtr->c_str(), NULL);
        c=FloatToUchar(dosage);
    }

//read GT
else{
    intc = 0;
    int alleleIndex;
    for(int i = 0; i < 2; i++)
    {
        alleleIndex = record.getGT(smNum, i);
        if(alleleIndex == VcfGenotypeSample::MISSING_GT)
        { intc = -2; break;}
        else if(alleleIndex > 0)
        { intc += 1; }
    }
    c = IntToUchar(intc);
    }
    return c;
}

double StringToDouble(const char* s) {
    double f = strtod(s, NULL);
    return f;
}

float StringToFloat(const char* s) {
    float f = strtof(s, NULL);
    return f;
}

float UcharToFloat(const uchar c){
    if (c != UCHAR_MAX) {
        int intc = c;
        return (((float)intc) * 0.01);
    }
    else return -9.0;
}

double UcharToDouble(const uchar c){
    if (c != UCHAR_MAX) {
        int intc = c;
        return (((double)intc) * 0.01);
    }
    else return -9.0;
}

uchar FloatToUchar(const float dosage){
    if (dosage >= 0.0  && dosage <= 2.0) {
        return  ((int)(dosage*100.0));
    }
    else return UCHAR_MAX;
}

uchar DoubleToUchar(const double dosage){
    if (dosage >= 0.0  && dosage <= 2.0) {
        return ((int)(dosage*100.0));
    }
    else return UCHAR_MAX;
}

uchar IntToUchar(const int intc){
    if (intc >= 0  && intc <= 2) {
        return (intc * 100);
    }
    else return UCHAR_MAX;
}

// get genotype vector for given column
void getGTgslVec(uchar ** X, gsl_vector *xvec, size_t marker_i, const size_t ni_test, const size_t ns_test, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize){
    
    if (marker_i < ns_test ) {
        
        double geno, geno_mean = 0.0;
        
            size_t compressedBufferSize = CompBuffSizeVec[marker_i];
            uchar * UnCompBuffer = (uchar*)malloc(UnCompBufferSize);
            
            int result = uncompress(UnCompBuffer, &UnCompBufferSize, X[marker_i],compressedBufferSize);
            
            if(result != Z_OK)
            {
                zerr(result);
                exit(-1);
            }
            else{
                for (size_t j=0; j<ni_test; j++) {
                    geno = UcharToDouble(UnCompBuffer[j]);
                    //cout << geno << ", ";
                    if (geno < 0.0 || geno > 2.0) {
                        cerr << "wrong genotype value = " << geno << endl;
                        exit(-1);
                    }
                    gsl_vector_set(xvec, j, geno);
                    geno_mean += geno;
                }
                geno_mean /= (double)ni_test;
                //geno_mean = -geno_mean;
                gsl_vector_add_constant(xvec, -geno_mean); // center genotypes here
            }
            free(UnCompBuffer);
    }
    else {
        std::cerr << "Marker index out of range \n";
        exit(-1);
    }
}

//used in MCMC()
void getGTgslVec(uchar ** X, gsl_vector *xvec, size_t marker_i, const size_t ni_test, const size_t ns_test, const vector<double> &SNPsd, const vector<double> &SNPmean, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize, bool Compress_Flag, const vector<pair<int, double> > &UcharTable){
    
    if (marker_i < ns_test ) {
        
        int c;
        double geno;
        double geno_mean = SNPmean[marker_i];
        /*if (geno_mean > 200 || geno_mean < 0)
        {
            cout << "\n geno_mean=" << geno_mean << endl;
            exit(-1);
         }*/
        
        if (Compress_Flag) {
            size_t compressedBufferSize = CompBuffSizeVec[marker_i];
            uchar * UnCompBuffer = (uchar*)malloc(UnCompBufferSize);
            
            int result = uncompress(UnCompBuffer, &UnCompBufferSize, X[marker_i],compressedBufferSize);
            
            if(result != Z_OK)
            {
                zerr(result);
                exit(-1);
            }
            else{
                for (size_t j=0; j<ni_test; j++) {

                    c = (int)UnCompBuffer[j];
                    geno = UcharTable[c].second;
                    //if(i < 20)  cout << geno << ",";
                    gsl_vector_set(xvec, j, geno);
                }
            }
            free(UnCompBuffer);
        }
        else{
            for (size_t j=0; j<ni_test; j++) {

               c = (int)X[marker_i][j];
                geno = UcharTable[c].second;
                //if(j < 20)  cout << geno << ",";
                gsl_vector_set(xvec, j, geno);
            }
        }
        gsl_vector_add_constant(xvec, -geno_mean);
    }
    else {
        std::cerr << "Marker index out of range \n";
        exit(-1);
    }
}


void getGTgslVec(uchar ** X, gsl_vector *xvec, size_t marker_i, const size_t ni_test, const size_t ns_test){

    double geno, geno_mean = 0.0;
    if (marker_i < ns_test ) {
        for (size_t j=0; j<ni_test; j++) {
            geno = UcharToDouble(X[marker_i][j]);
            //cout << geno << ", ";
            if (geno < 0.0 || geno > 2.0) {
                cerr << "wrong genotype value = " << geno << endl;
                exit(1);
            }
            gsl_vector_set(xvec, j, geno);
            geno_mean += geno;
        }
        geno_mean /= (double)ni_test;
        //geno_mean = -geno_mean;
        gsl_vector_add_constant(xvec, -geno_mean); // center genotypes here
    }
    else {
        std::cerr << "Error return genotype vector...\n";
        exit(1);
    }
    //cout << endl;
}

//get genotype matrix for given column vector
bool getGTgslMat(uchar ** X, gsl_matrix *Xgsl, std::vector<size_t> marker_idx, const size_t ni_test, const size_t ns_test){
    
    size_t marker_i;
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
    
    for (size_t i=0; i < marker_idx.size(); i++) {
        marker_i = marker_idx[i];
        getGTgslVec(X, xvec, marker_i, ni_test, ns_test);
        gsl_matrix_set_col(Xgsl, i, xvec);
    }
    
    gsl_vector_free(xvec);
    return 1;
}

// print sub genotype uchar array
bool print(const char* description, uchar **genotypes, uint numMarkers, uint numSamples, std::vector<String> &sampleIDs)
{
    std::cout << "\n\t\t\t" << description << "\n\t\t";
    std::cout << "Sample ID : ";
    std::cout << "\n";
    for(uint i = 0; i < numSamples; i++)
    {
        std::cout << "\t" << sampleIDs[i];
    }
    std::cout << "\n";
    
    uchar c;
    float fc;
    std::cout << "\n" << "Dosage GT: \n";
    for(uint i = 0; i < numMarkers; i++)
    {
      for(uint j = 0; j < numSamples; j++)
        {
            std::cout << "\t";
            if(genotypes[i][j] != UCHAR_MAX)
            {
                c = genotypes[i][j];
                fc = (float)c;
                std::cout << fc * 0.01;
            }
            else
                std::cout << "GT missing";
        }
        std::cout << "\n";
    }
    
    std::cout << "Succesfully output sub-genotype\n";
    return 1;
}

void print(uchar **UtX, uint numMarkers, uint numSamples, std::vector <size_t> &CompBuffSizeVec, size_t UnCompBufferSize)
{
    int intc;
    float fc;
    std::cout << "\n" << "Dosage GT: \n";
    size_t compressedBufferSize;
    uchar * TempBuffer = (uchar *)malloc(UnCompBufferSize);
    
    for(uint i = 0; i < numMarkers; i++)
    {
      // cout << "compressedBufferSize = " << CompBuffSizeVec[i] << endl;
        compressedBufferSize = CompBuffSizeVec[i];
        
        int result = uncompress(TempBuffer, &UnCompBufferSize, UtX[i],compressedBufferSize);
        if (result != Z_OK) {
            zerr(result);
            exit(-1);
        }
        else{
            //cout << "Adjusted UnCompBuffer with size : " << UnCompBufferSize << endl;
            for(uint j = 0; j < numSamples; j++)
            {
             intc = (int)TempBuffer[j];
            // std::cout << intc << ":" ;
             if(TempBuffer[j] != UCHAR_MAX)
             {
                 fc = (float)(intc);
                 std::cout << fc * 0.01 << ", ";
             }
             else
                 std::cout << "GT missing";
            }
            std::cout << "\n";
        }
    }
    free(TempBuffer);
    std::cout << "Succesfully output sub-genotype\n";
}


bool print(uchar **genotypes, uint numMarkers, uint numSamples)
{
    int intc;
    float fc;
    std::cout << "\n" << "Dosage GT: \n";
    for(uint i = 0; i < numMarkers; i++)
    {
        for(uint j = 0; j < numSamples; j++)
        {
            std::cout << "\t";
            intc = (int)genotypes[i][j];
            std::cout << intc << ":" ;
            
            if(genotypes[i][j] != UCHAR_MAX)
            {
                fc = (float)(intc);
                std::cout << fc * 0.01;
            }
            else
                std::cout << "GT missing";
        }
        std::cout << "\n";
    }
    
    std::cout << "Succesfully output sub-genotype\n";
    return 1;
}




