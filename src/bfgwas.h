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

#ifndef __BFGWAS_H__                
#define __BFGWAS_H__

#include "param.h"
#include "VcfFileReader.h"
#include "ReadVCF.h"
#include "mathfunc.h" 

using namespace std;

class BFGWAS {

public:			
	//parameters
	string version;
	string date;
	string year;
	
	//constructor
	BFGWAS(void);
	
	//functions
	void PrintHeader (void);
	void PrintHelp (size_t option);
	void PrintLicense (void);
	void Assign (int argc, char **argv, PARAM &cPar);
	void BatchRun (PARAM &cPar);
	void WriteLog (int argc, char **argv, PARAM &cPar);
};


#endif

