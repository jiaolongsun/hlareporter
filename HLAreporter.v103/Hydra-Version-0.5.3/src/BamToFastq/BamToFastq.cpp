/*
  ***************************************************************************
   bamFilterPaired.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.

   Filters BAM alignments based upon user-defined criteria.
 ***************************************************************************
*/

#include "BamToFastq.h"

// constructor
BamToFastq::BamToFastq(string bamFile, string fastq1, string fastq2, bool useMateTags)
:_bamFile(bamFile)
,_fastq1(fastq1)
,_fastq2(fastq2)
,_useMateTags(useMateTags)
{ }


// destructor
BamToFastq::~BamToFastq(void) {}


// set the alignments for each end of this BAM pair. 
bool BamToFastq::CheckPairEnds(BamAlignment &end1, BamAlignment &end2) {
	_end1 = end1;
	_end2 = end2;
	
	// make sure the ends are from the same read-pair.  else, gripe.
	if (_end1.Name == _end2.Name) return true;
	else return false;
}


void BamToFastq::ConvertBamToFastq() {
	
	// open the 1st fastq file for writing
	ofstream fq1(_fastq1.c_str(), ios::out);
	if ( !fq1 ) {
		cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
		exit (1);
	}
	// open the 2nd fastq file for writing
	ofstream fq2(_fastq2.c_str(), ios::out);
	if ( !fq2 ) {
		cerr << "Error: The second fastq file (" << _fastq2 << ") could not be opened.  Exiting!" << endl;
		exit (1);
	}
	
	// maps to store the whether a read has been reported or not.
	map<string, bool, less<string> > fq1Reads, fq2Reads;
	
	
	// open the BAM file
	BamReader reader;	
	reader.Open(_bamFile);
	
	BamAlignment alignment;
	while (reader.GetNextAlignment(alignment)) {
		
		// extract the sequence and qualities for the BAM "query"
		string sequence  = alignment.QueryBases;
		string qualities = alignment.Qualities;
		if (alignment.IsReverseStrand() == true) {
			CSequenceUtilities::GetReverseComplement(sequence);
			CSequenceUtilities::ReverseSequence(qualities);
		}
		
		if (_useMateTags == false) {
			
			// report the read/end if it has not been reported already
			if (alignment.IsFirstMate() == true) {
				if (fq1Reads[alignment.Name] != true) {
					fq1 << "@" << alignment.Name << "/1" << endl;
					fq1 << sequence << endl;
					fq1 << "+" << endl;
					fq1 << qualities << endl;
				
					// note that this read has been processed
					fq1Reads[alignment.Name] = true;				
				}			
			}
			else {
				if (fq2Reads[alignment.Name] != true) {
					fq2 << "@" << alignment.Name << "/2" <<endl;
					fq2 << sequence << endl;
					fq2 << "+" << endl;
					fq2 << qualities << endl;
					
					// note that this read has been processed
					fq2Reads[alignment.Name] = true;							
				}
			}
		}
		else {  // construct FASTQ entries from the query info and the R2 and Q2 tags.
			
			// assume the R2 and Q2 tags are on the + strand.
			string mateSequence, mateQualities;
			alignment.GetMateSequence(mateSequence);
			alignment.GetMateQualities(mateQualities);
			
			// report the read/end if it has not been reported already
			if (fq1Reads[alignment.Name] == false && fq2Reads[alignment.Name] == false) {
				if (alignment.IsFirstMate() == true) {
					fq1 << "@" << alignment.Name << "/1" << endl;
					fq1 << sequence << endl;
					fq1 << "+" << endl;
					fq1 << qualities << endl;
		
					// note that this read has been processed
					fq1Reads[alignment.Name] = true;
			
					fq2 << "@" << alignment.Name << "/2" <<endl;
					fq2 << mateSequence << endl;
					fq2 << "+" << endl;
					fq2 << mateQualities << endl;
			
					// note that this read has been processed
					fq2Reads[alignment.Name] = true;	
				}
				else {
					fq2 << "@" << alignment.Name << "/2" <<endl;
					fq2 << sequence << endl;
					fq2 << "+" << endl;
					fq2 << qualities << endl;
				
					// note that this read has been processed
					fq2Reads[alignment.Name] = true;
					
					fq1 << "@" << alignment.Name << "/1" <<endl;
					fq1 << mateSequence << endl;
					fq1 << "+" << endl;
					fq1 << mateQualities << endl;
			
					// note that this read has been processed
					fq1Reads[alignment.Name] = true;						
				}
			}
		}
	}
	// close up shop
	reader.Close();
}
