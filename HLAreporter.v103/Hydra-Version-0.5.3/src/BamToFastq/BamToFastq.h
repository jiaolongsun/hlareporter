/*
***************************************************************************
bamToFastq.h (c) 2009 Aaron Quinlan

Hall Lab
Department of Biochemistry and Molecular Genetics
University of Virginia

All rights reserved.

Filters BAM alignments based upon user-defined criteria.
***************************************************************************
*/
#include "BamAux.h"
#include "BamReader.h"
using namespace BamTools;

#include "SequenceUtilities.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <map>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BamToFastq {

public:

	// constructor 
	BamToFastq(string bamFile, string fastq1, string fastq2, bool useMateTags);

	// destructor
	~BamToFastq(void);

	void ConvertBamToFastq();
		
private:

	// returns true if valid ends of same pair.
	bool CheckPairEnds(BamAlignment &end1, BamAlignment &end2);

	void PrintPairAsFastq();
		
	string _bamFile;

	BamAlignment _end1;
	BamAlignment _end2;

	string _fastq1, _fastq2;	// the names of the fastq output files
	bool _useMateTags;			// whether or not the mate sequence should be 
								// extracted from the R2 BAM tag.
};
