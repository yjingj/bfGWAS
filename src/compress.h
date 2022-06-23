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

    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "zlib.h"
#include "zconf.h"

//#include "stdafx.h"   // not actually needed
//#define ZLIB_WINAPI   // actually actually needed (for linkage)

//#include "windows.h"  // get BYTE et al.
//#pragma comment(lib, "zlibwapi.lib") // for access to the DLL


#ifndef uchar
typedef unsigned char uchar;
#endif
typedef unsigned char BYTE;



void zerr(int ret);

int GetMaxCompressedLen(int nLenSrc );

int CompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );

int UncompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );
