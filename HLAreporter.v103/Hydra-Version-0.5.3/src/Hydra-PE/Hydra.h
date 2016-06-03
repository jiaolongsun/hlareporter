#include <algorithm>
#include <numeric>    // for the accumulate algorithm  
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <map>
#include <vector>

using namespace std;


// prevent duplicate includes of this header.
#pragma once


// Constants
const double UNIQ_WEIGHT = 1.0;
const double ANCH_WEIGHT = 0.5;
const double MULT_WEIGHT = 0.1;

const char UNIQ_TYPE = 1;
const char ANCH_TYPE = 2;
const char MULT_TYPE = 3;


// Main structure to store a specific mapping for a paired-end alignment
struct PAIR {

	string readId;                		// the read id.

	string chrom1, chrom2;     			// the mapped chrom for read1, read2
	short mate1, mate2;                 // which end of the mate-pair
	int start1, start2;           		// the mapped start position for read1, read2
	int end1, end2;           			// the mapped end position for read1, read2
	string strand1, strand2;   			// the mapping orientation for read1, read2

	unsigned short edit1, edit2;        // the number of mismatches for read1, read2
	unsigned int mappings1, mappings2;  // the number of mappings observed on each end.
	char mappingType;        		    // the mapping type: 1 = dis_unique, 2 = anchor, 3 = multiple, 4 = orphan

	unsigned int nonOverlap;            // this is a temporary variable used to compute the 
									    // number of bases in this pair the DO NOT overlap
									    // with the current contig seed that is being evaluated.
									
	unsigned int support;               // the number of reads in support of this read
	bool used;
	bool include;
};


// structure for the coordinates of an alignment
// Using this struct as a key to a map.  Thus, need the "<" overload so that
// the map can compare BED3.
struct BED3 {
	
	// data
	string chrom;
	int start;
	int end;
	
	// constructor
	BED3 (const string c, const int s, const int e) :
		chrom(c), 
		start(s), 
		end(e) 
	{}
	
	// comparison operator for maps keyed on this structure
	bool operator < (const BED3 &a) const
	{ 
		// first, sort by chrom (asc)
		if (chrom < a.chrom)      return true;
		else if (chrom > a.chrom) return false;

		// second, within a chrom sort by start
		if (start < a.start)      return true;
		else if (start > a.start) return false;
		
		if (end < a.end)      return true;
		else return false;
	}
};


// structure to track the "contigs" that a given mapping 
// has been assigned to. 
struct IN_CONTIG {
	int contig;
	bool include;
	double weightedSupport;
	int totalMM;
	int contigSize;
	bool unlinked;
};


// Using this struct as a key to a map.  Thus, need the "<" overload so that
// the map can compare CHROMS_AND_STRANDS.
struct CHROMS_AND_STRANDS {
	// data
	string chrom1;
	string chrom2;
	string strand1;
	string strand2;
	
	// constructor
	CHROMS_AND_STRANDS (const string chrom1, const string chrom2,
						const string strand1, const string strand2) :
		chrom1(chrom1),
		chrom2(chrom2),
		strand1(strand1),
		strand2(strand2) 
	{}
	
	// comparison operator for maps keyed on this structure
	bool operator < (const CHROMS_AND_STRANDS &a) const
	{ 
		// first, sort by chrom1 (asc)
		if (chrom1 < a.chrom1)      return true;
		else if (chrom1 > a.chrom1) return false;

		// second, sort by chrom2 (asc)
		if (chrom2 < a.chrom2)      return true;
		else if (chrom2 > a.chrom2) return false;
		
		// third, sort by strand1 (asc)
		if (strand1 < a.strand1)      return true;
		else if (strand1 > a.strand1) return false;
		
		// lastly, sort by strand2 (asc)
		if (strand2 < a.strand2) return true;
		else return false;
	}
};



// TYPEDEFS for common data structures


typedef map<CHROMS_AND_STRANDS, vector< vector<PAIR> >, std::less<CHROMS_AND_STRANDS> > putativeContigMap;

typedef multimap<int, vector<PAIR>, std::less<int> > circosContigMap;

typedef map<string, vector<IN_CONTIG>, std::less<string> > read2ContigsMap;

typedef map<string, vector<PAIR>, std::less<string> > read2MappingsMap; 

// vector of paired-end mappings
typedef vector<PAIR> pairVector;

// map of paired-end mappings key by chrom-chrom-strand-strand
typedef map<CHROMS_AND_STRANDS, vector<PAIR>, std::less<CHROMS_AND_STRANDS> > pairMap;

// vector of contigs containing all mappings in each contig
typedef vector< vector<PAIR> > contigVector;


//*************************************************
// Template Functions
//*************************************************

// templated function to convert objects to strings
template <typename T>
std::string ToString(const T & value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}



class HydraPE {

public:

	//************************************************
	// Public methods and elements
	//************************************************	
	
	// constructor 
	HydraPE(int lengthDev, int spanDev, int minSupport, int maxLinkedDistance, 
				bool ignoreSize, bool lumpInversions, string mappingUsage, int editBeyondBest);

	// destructor
	~HydraPE(void);

	void Tokenize(const string &, vector<string> &, const string &);  	// tokenize a line into a vector of strings
	void ReadDiscordantFile(const string &discordantFile);				// process a discordant file in search of SVs
	void SortFragments();												// find clusters of discordant mappings
	void AssembleContigs();												// assemble clusters of discordant mappings
	void AllowOneContigPerRead();										// find the most likely contig for each pair
	void ReportSVCalls(string &);										// report the final SV calls to output files.


private:

	//************************************************
	// Private methods and elements
	//************************************************
	
	// private processing variables
	int lengthDev;
	int spanDev;
	int minSupport;
	int maxLinkedDistance;
	bool ignoreSize;
	bool lumpInversions;
	string mappingUsage;
	int editBeyondBest;
	
	// private data structures for main processing.
	pairMap mappingsByChromAndStrand;
	putativeContigMap putativeContigs;
	contigVector contigs;
	read2ContigsMap read2Contigs;
	read2MappingsMap read2Mappings;
	
	// private methods
	void CullMappingsByMisMatches(pairVector &mappings);	
	void CorrectMateOrder(PAIR &pair);
	void SwapEnds(PAIR &pair);
	void AddMappingsToMasterMap(const pairVector &mappings);
};
