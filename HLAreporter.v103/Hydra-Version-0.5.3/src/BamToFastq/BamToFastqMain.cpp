/*
  ***************************************************************************
   bamToFastqMain.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.
 ***************************************************************************
*/

#define PROGRAM_NAME "bamToFastq"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "BamToFastq.h"
#include "Version.h"
using namespace std;

// BamTools includes
#include "BamReader.h"
using namespace BamTools;



// function declarations
void ShowHelp(void);
	

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	bool haveInBam     = false;
	bool haveFastq1    = false;
	bool haveFastq2    = false;
	bool useMateTags   = false;
		
	// input files
	string inBamFile;
	
	//output files
	string fastq1, fastq2;

	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();


	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

 		if (PARAMETER_CHECK("-bam", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveInBam = true;
				inBamFile = argv[i + 1];
				i++;
			}
		}
		else if (PARAMETER_CHECK("-fq1", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveFastq1 = true;
				fastq1 = argv[i + 1];
				i++;
			}
		}
		else if (PARAMETER_CHECK("-fq2", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveFastq2 = true;
				fastq2 = argv[i + 1];
				i++;
			}
		}
		else if (PARAMETER_CHECK("-useTags", 8, parameterLength)) {
			useMateTags = true;
		}					
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have a BAM file
	if (!haveInBam) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -bam. " << endl << "*****" << endl;
	  showHelp = true;
	}
	// make sure we have an end1 FASTQ file
	if (!haveFastq1) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -fq1. " << endl << "*****" << endl;
	  showHelp = true;
	}
	// make sure we have an end2 FASTQ file
	if (!haveFastq2) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -fq2. " << endl << "*****" << endl;
	  showHelp = true;
	}

	
	// let 'er rip.
	if (!showHelp) {
		
		BamToFastq b2fq(inBamFile, fastq1, fastq2, useMateTags);
		
		// convert the alignments to fastq
		b2fq.ConvertBamToFastq();
	}
	else {
		ShowHelp();
	}
}


void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Creates paired-end FASTQ files from a BAM file" << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <BAM> -fq1 <FQ> -fq2 <FQ>" << endl << endl;

	cerr << "Options:" << endl;
	cerr << "  -useTags\tCreate FASTQ based on the mate info" << endl;
	cerr << "  \tin the BAM R2 and Q2 tags." << endl << endl;
	
	exit(1);	
}

