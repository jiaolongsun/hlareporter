/******************************************************************************
Program: 		hydra pem
Description: 	A comprehensize algorithm for detecting structural 
				variation in all genomic regions using paired-end mapping.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*******************************************************************************/
#include "Hydra.h"
#include "Sort.h"
#include "Ancillary.h"


// make
HydraPE::HydraPE(int lengthDev, int spanDev, int minSupport, int maxLinkedDistance, 
						 bool ignoreSize, bool lumpInversions, string mappingUsage, int editBeyondBest) 
: lengthDev(lengthDev)
, spanDev(spanDev)
, minSupport(minSupport)
, maxLinkedDistance(maxLinkedDistance)
, ignoreSize(ignoreSize)
, lumpInversions(lumpInversions)
, mappingUsage(mappingUsage)
, editBeyondBest(editBeyondBest)
{}


// and break
HydraPE::~HydraPE(void) {}


/******************************************************************************
* ReadDiscordantFile
*
* Loads the discordant mapping file into a map that can be used to 
* assemble contigs from "similarly-dsicordant" mappings.
******************************************************************************/
void HydraPE::ReadDiscordantFile(const string &discordantFile) {

	// open the mapping files for reading
	ifstream mappings(discordantFile.c_str(), ios::in);

	if ( !mappings ) {
		cerr << "Error: The requested mappings file (" << discordantFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	else if (mappings.is_open()) {

		string mappingLine;
		string read2; 			// placeholder for reading from file
		string currRead = "";
		string prevRead = "";
		PAIR line;
		int whichMateIsBlock1;
		
		// vector of mappings for the current read.
		vector<PAIR> currReadMappings;

		// reserve some memory so we don't have to reallocate (as much)
		currReadMappings.reserve(1000);      

		// loop through the current mapping file and store each line in a PAIR struct.
		while (mappings >> line.chrom1  >> line.start1  >> line.end1
						>> line.chrom2  >> line.start2  >> line.end2
						>> line.readId  >> whichMateIsBlock1
						>> line.strand1 >> line.strand2
						>> line.edit1   >> line.edit2) {

			// get the id for the current mapping							
			currRead = line.readId;	
			
			// set the ends of the pair
			if (whichMateIsBlock1 == 1) {
				line.mate1 = 1;
				line.mate2 = 2;
			}
			else {
				line.mate1 = 2;
				line.mate2 = 1;				
			}	 
			
			// initialize the other fields for the PAIR structure
			line.support   = 0;
			line.used      = false;
			line.include   = true;

			// ensure that the mappings for each end
			// of the pair are in the correct order.
			CorrectMateOrder(line);

			/*--------------------------------------------------------------------
			  Add this mapping to the list of mappings for this READ.
			  We will use all of the mappings for a read to reduce the number of 
			  mappings to those that have the least number of mismatches.	  
			  --------------------------------------------------------------------*/
			if ( (currRead != prevRead) && (prevRead != "") ) {
				
				// reduce the mappings for this pair to those with the least
				// edit distance and add the remaning mappings to
				CullMappingsByMisMatches(currReadMappings);
				AddMappingsToMasterMap(currReadMappings);

				// reset the list of mappings for the current read
				// and add the first mapping for the current read
				currReadMappings.clear();
				currReadMappings.push_back(line);
			}
			else {
				currReadMappings.push_back(line);
			}
			prevRead = currRead;
		}
		// add the best mappings for the last pair in the file.
		CullMappingsByMisMatches(currReadMappings);
		AddMappingsToMasterMap(currReadMappings);
	}
}


/******************************************************************************
* CorrectMateOrder
*
* Swaps the ends of a pair so that:
*    1. the leftmost end is first when an intra-chromosomal.
*	 2. the lexigraphically lower chrom is first when an inter-chromosomal.
******************************************************************************/
void HydraPE::CorrectMateOrder(PAIR &pair) {
	if (pair.chrom1 == pair.chrom2) {
		if (pair.start2 < pair.start1) {
			SwapEnds(pair);
		}
	}
	else {
		if (pair.chrom2 < pair.chrom1) {
			SwapEnds(pair);
		}	
	}
}


/******************************************************************************
* SwapEnds
*
* Swaps the ends of a pair.
******************************************************************************/
void HydraPE::SwapEnds(PAIR &pair) {

	// hold end1 in temp
	string tmpChrom        = pair.chrom1;
	short tmpMate          = pair.mate1;
	int tmpStart           = pair.start1;
	int tmpEnd             = pair.end1;
	string tmpStrand       = pair.strand1;
	unsigned short tmpEdit = pair.edit1;
	
	// set end1 to end2
	pair.chrom1  = pair.chrom2;
	pair.mate1   = pair.mate2;
	pair.start1  = pair.start2;
	pair.end1    = pair.end2;
	pair.strand1 = pair.strand2;
	pair.edit1   = pair.edit2;
	
	// set end2 to end1
	pair.chrom2  = tmpChrom;
	pair.mate2   = tmpMate;
	pair.start2  = tmpStart;
	pair.end2    = tmpEnd;
	pair.strand2 = tmpStrand;
	pair.edit2   = tmpEdit;
}


/******************************************************************************
* AddMappingsToMasterMap
*
* Creates a map where the key is chrom1,chrom2,strand1,strand2 and the value
* is a vector of all mappings with the same chroms and strands.  In effect, 
* constructing this map does much of the sorting work for us ahead of time.
******************************************************************************/
void HydraPE::AddMappingsToMasterMap(const pairVector &mappings) {

	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		// force both types of inversion mappings to be the same to increase sensitivity.
		if (this->lumpInversions == true) {
			if ( ((mapIter->strand1 == "+") && (mapIter->strand2 == "+")) || ((mapIter->strand1 == "-") && (mapIter->strand2 == "-")) ) {
				CHROMS_AND_STRANDS key(mapIter->chrom1, mapIter->chrom2, "+", "+");
				this->mappingsByChromAndStrand[key].push_back(*mapIter);
			}
			else {
				CHROMS_AND_STRANDS key(mapIter->chrom1, mapIter->chrom2, mapIter->strand1, mapIter->strand2);
				this->mappingsByChromAndStrand[key].push_back(*mapIter);
			}	
		}
		else {
			CHROMS_AND_STRANDS key(mapIter->chrom1, mapIter->chrom2, mapIter->strand1, mapIter->strand2);
			this->mappingsByChromAndStrand[key].push_back(*mapIter);
		}
	}
}


void HydraPE::SortFragments() {

	// vector to store the putative contig blocks
	vector< vector<PAIR> > lengthAndPosClusters;

	// Loop through the reads in each set
	pairMap::iterator iter = this->mappingsByChromAndStrand.begin();

	while (iter != this->mappingsByChromAndStrand.end() ) {

		CHROMS_AND_STRANDS key = iter->first;
		vector<PAIR> reads = iter->second;

		//****************************************************************
		// Phase 1: Sort the fragments from this chrom-chrom-strand-strand
		//          set by length.  This will be the first cut at
		//          indentifying contigs
		//****************************************************************
		sort(reads.begin(), reads.end(), byFragmentSize);

		// Create a vector to see if the adjacent (sorted) pair lengths
		// are within our expected length deviation.
		// Those that are form putative contig blocks.
		vector<bool> lengthDiffs(reads.size());

		// for each adjacent mapping (sorted by length), flag whether or 
		// not the lengths are within our expected length deviation.
		for (int i = 0; i < reads.size() - 1; i++) {
			lengthDiffs[i] = doLengthsSupportOneAnother(reads[i+1], reads[i], this->lengthDev);
		}
		
		// vector to store the putative contig blocks
		vector< vector<PAIR> > lengthClusters;

		// loop through the lengthDiff flags from above and find 
		// putative contig blocks
		for (int i = 0; i < lengthDiffs.size(); i++) {

			// have we found the start of a new putative block?
			if (lengthDiffs[i]) {

				vector<PAIR> cluster;
				// add the first pair
				cluster.push_back(reads[i]);
				i++;

				// keep adding pairs until the block ends
				while (lengthDiffs[i]) {
					cluster.push_back(reads[i]);
					i++;			
				}
				// add the last pair in the block (the last pair is always false per the above loop)
				cluster.push_back(reads[i]);
				lengthClusters.push_back(cluster);
			}
		}

		//****************************************************************
		// Phase 2: Now sort the fragments by the first pos1.  The "blocks"
		//          that remain are blocks of true contigs.  We will send
		//          each pair from the blocks to a "brute force", N x N
		//          contig maker.
		//****************************************************************

		for (int l = 0; l < lengthClusters.size(); l++) {

			vector<PAIR> currLengthCluster = lengthClusters[l];
			sort(currLengthCluster.begin(), currLengthCluster.end(), byStart1);	    

			// Create a vector to see if the adjacent (sorted) pair lengths
			// are within our expected length deviation.
			// Those that are form putative contig blocks.
			vector<bool> nonOverlapDiffs(currLengthCluster.size());

			// for each adjacent mapping (sorted by length), flag whether or
			// not the lengths are within our expected span deviation.
			for (int lc = 0; lc < currLengthCluster.size() - 1; lc++) {		
				nonOverlapDiffs[lc] = doSpansSupportOneAnother(currLengthCluster[lc], currLengthCluster[lc+1], this->spanDev);				
			}

			// loop through the nnOverlap flags from above and find
			// putative contig blocks
			for (int n = 0; n < nonOverlapDiffs.size(); n++) {

				// have we found the start of a new putative block?
				if (nonOverlapDiffs[n]) {

					vector<PAIR> posCluster;

					// add the first pair
					posCluster.push_back(currLengthCluster[n]);
					n++;

					// keep adding pairs until the block ends
					while (nonOverlapDiffs[n]) {
						posCluster.push_back(currLengthCluster[n]);
						n++;
					}
					// add the last pair in the block (the last pair is always false per the above loop)
					posCluster.push_back(currLengthCluster[n]);
					lengthAndPosClusters.push_back(posCluster);
				}
			}	  
		}
		this->putativeContigs[key] = lengthAndPosClusters;

		// move onto the next set of 
		// chrom-chrom-strand-strand set of mappings
		lengthAndPosClusters.clear();
		iter++;
	}
}



void HydraPE::AssembleContigs() {

	putativeContigMap::iterator iter = this->putativeContigs.begin();

	while (iter != this->putativeContigs.end() ) {

		// the "key" is the chrom/chrom/strand/strand
		CHROMS_AND_STRANDS key = iter->first;

		// the list of all the putative contigs for this key
		vector< vector<PAIR> > putativeContigs = iter->second;
		
		// Loop through all of the putative contigs to test whether or 
		// not they hold up.
		for (int c = 0; c < putativeContigs.size(); c++) {

			vector<PAIR> currContig = putativeContigs[c];
			int  nonOverlap;
			int  lengthDiff;
			int  support;
			long totalSupport;

			//***************************************************************
			// Compute the support for each read in this class,
			// regardless of the mapping class.
			//***************************************************************
			for (int i = 0; i < currContig.size(); i++) {

				// NEW: Requires 1bp of overlap
				support = 0;
				for (int j = 0; j < currContig.size(); j++) {
					if (j != i) {
						bool haveLengthSupport = doLengthsSupportOneAnother(currContig[i], currContig[j], this->lengthDev);
						bool haveSpanSupport   = doSpansSupportOneAnother(currContig[i], currContig[j], this->spanDev);
						if ( haveLengthSupport && haveSpanSupport ) support++;	
					}
				}
				// update the support for this mapping
				currContig[i].support = support;
				totalSupport += support;
			}


			// sort the reads first by mappingType (asc), then by support (desc). 
			sort(currContig.begin(), currContig.end(), byTypeAndSupport);

			//*******************************************************************
			// Continue as long as:
			// 1. There are reads left to assemble
			// 2. The current seed (reads[0]) has support. 
			// NOTE: after each cycle, reads[0] will be the potential seed.
			//*******************************************************************
			std::map<string, bool> readAlreadyUsed;
			int pass = 0;

			while(currContig.size() > 1) {
				pass++;
				vector<PAIR> contig;

				if (currContig[0].support > 0) {

					// the list of pairs that support the seed...aka a contig
					contig.push_back(currContig[0]);

					// remove the seed from the pool of pairs
					currContig[0].used = true;
					readAlreadyUsed[currContig[0].readId] = true;

					// compute the non-overlap between all the other reads and the seed
					for (int i = 0; i < currContig.size(); i++) {
						currContig[i].nonOverlap = getNonOverlap(currContig[i], contig[0]);
					}

					// sort the reads by non-overlap relative to the seed.
					sort(currContig.begin(),currContig.end(),byNonOverlap);

					// Now, loop through the pairs in order of how well they overlap
					// with the seed.  That is, in ascending order of "non-overlap".
					for (int currPair = 1; currPair < currContig.size(); currPair++) {

						bool supportsAll = true;

						// Make sure the currPair is in support with ALL of the other pairs in the contig
						int contigSize = contig.size();
						
						for(int contigPair = 0; contigPair < contigSize; contigPair++) {
							bool haveLengthSupport = doLengthsSupportOneAnother(contig[contigPair], currContig[currPair], this->lengthDev);
							bool haveSpanSupport   = doSpansSupportOneAnother(contig[contigPair], currContig[currPair], this->spanDev);
							if ( haveLengthSupport == false || haveSpanSupport == false) {
								supportsAll = false;
								break;  // no need to continue testing if it doesn't support >=1 other mapping
							}
						}
						
						// If the current pair is in support with ALL the other
						// pairs in the contig, then add it to the contig and remove it from
						// the working list of pairs.
						if (supportsAll) {
							contig.push_back(currContig[currPair]);
						}
					}

					// add this contig to the list of all contigs
					// IF there is more than just the seed
					if (contig.size() >= this->minSupport) {
						
						this->contigs.push_back(contig);

						for (int c = 0; c < contig.size(); c++) {
							// remove each read from the "pool" of mappings
							currContig[c].used = true;
							totalSupport       -= currContig[c].support;
						}
						contig.clear();
					}
					else contig.clear();

					// sort the reads so the unused mappings are at the "top"
					sort(currContig.begin(),currContig.end(),byUsed);

					// delete the used mappings
					vector<PAIR>::iterator it = currContig.begin();
					while ((it != currContig.end()) && !(it->used)) {
						it++;
					}
					currContig.erase(it, currContig.end());

					// resort the reads by type and support for the next seed search.
					sort(currContig.begin(),currContig.end(),byTypeAndSupport);	    
				}
				else {
					currContig.erase(currContig.begin());
				}	
			}
			currContig.clear();
		}

		// move on to the next set of putative chr-chr-str-str contigs
		iter++;
		this->putativeContigs.erase(key);
	}
}



//***************************************************
// Method: AllowOnlyOneContigPerRead
// Purpose:
//
// Outline:
//***************************************************
void HydraPE::AllowOneContigPerRead() {

	int contigId = 0;

	// loop through all the contigs that were found
	for (int j = 0; j < this->contigs.size(); j++) {
		
		// used to track the "footprint" of contig
		int minStart1 = INT_MAX;					
		int minStart2 = INT_MAX;			
		int maxEnd1   = 0;
		int maxEnd2   = 0;
		
		// get the mappings that comprise this contig
		vector<PAIR> mappings = contigs[j];

		// get the size of the contig
		int contigSize         = getClusterSizeAll(mappings, minStart1, minStart2, maxEnd1, maxEnd2);

		// get the total mismatches among all of the mappings for this contig.
		int totalMM            = getTotalMMAmongAllMappings(mappings);
		
		// compute the weighted support for this contig
		double weightedSupport = getTotalWeightedSupportAmongAllMappings(mappings);

		// check to see if this contig is unlinked or not.
		bool unlinked          = isVariantUnlinked(mappings, this->maxLinkedDistance);
		 
		// for each read in the contig, add the contig summary information 
		// so that we can decide which contig is best for this read.
		for (int m = 0; m < mappings.size(); m++) {
			IN_CONTIG contigInstance;
			contigInstance.contig          = j;
			contigInstance.include         = true;
			contigInstance.weightedSupport = weightedSupport;
			contigInstance.totalMM         = totalMM;
			contigInstance.contigSize      = contigSize;
			contigInstance.unlinked        = unlinked;
			
			this->read2Contigs[mappings[m].readId].push_back(contigInstance);
		}
	}    

	// Now, find the best contig for each read and then exlude the read
	// from the inferior contigs
	read2ContigsMap::const_iterator iter = this->read2Contigs.begin();
	while (iter != this->read2Contigs.end() ) {

		string read = iter->first;
		vector<IN_CONTIG> readContigs = iter->second;

		// only bother choosing a best contig if there's more than one choice
		if (readContigs.size() > 1) {

			bool anyUnlinked = areAnyContigsUnlinked(readContigs);
			
			if ( (anyUnlinked = true) || (this->ignoreSize == true) ) {
				// sort the contig with the most support and least mismatches to the "top"
				sort(readContigs.begin(),readContigs.end(),byWeightedSupportAndTotalMM);
			}
			else {
				// sort the contig with the most support and __SMALLEST SIZE__ to the "top"
				sort(readContigs.begin(),readContigs.end(),byWeightedSupportAndSize);
			}

			// exclude this read from inferior contigs
			for (int j = 1; j < readContigs.size(); j++) {
				// get the contig id to which the read DOES NOT belong.
				int contig            = readContigs[j].contig;
				// get the mappings that comprise this contig
				vector<PAIR> mappings = contigs[contig];
				// exclude this read from this contig
				for (int m = 0; m < mappings.size(); m++) {
					if (mappings[m].readId == read) {
						this->contigs[contig][m].include = false;
					}
				}
			}
		}
		iter++;
	}
} 



//***************************************************
// Method: SummarizeContigs
// Purpose:
//
// Outline:
//***************************************************
void HydraPE::ReportSVCalls(string &outStub) {

	string allFile    = outStub;
	string finalFile  = outStub;
	string detailFile = outStub;
	
	detailFile.append(".detail");
	finalFile.append(".final");
	allFile.append(".all");
	
	// open the all contig file for writing
	ofstream all(allFile.c_str(), ios::out);
	if ( !all ) {
		cerr << "Error: The file of _all_ breakpoints (" << allFile << ") could not be opened.  Exiting!" << endl;
		exit (1);
	}
	
	// open the finall contig file for writing
	ofstream final(finalFile.c_str(), ios::out);
	if ( !final ) {
		cerr << "Error: The file of _final_ breakpoints (" << finalFile << ") could not be opened.  Exiting!" << endl;
		exit (1);
	}
	
	// open the detail contig file for writing
	ofstream detail(detailFile.c_str(), ios::out);
	if ( !detail ) {
		cerr << "Error: The file of breakpoints details (" << detailFile << ") could not be opened.  Exiting!" << endl;
		exit (1);
	}

	int contigId = 0;						// incremental id for each event that has the required support.
	int numMappings;						// number of mappings in contig
	
	int minStart1, minStart2;				// used to track the "footprint" of contig
	int maxEnd1, maxEnd2;
	
	int totalMM1, totalMM2;					// summary info
	int numUniqueMappers, numAnchoredMappers, numMultipleMappers;
	int totalMappings1, totalMappings2;
	int totalQual1, totalQual2;
	int meanMappings1, meanMappings2;
	int meanQual1, meanQual2;
	float meanMM1, meanMM2;
			
	double allWeightedSupport;
	int finalSupport;
	double finalWeightedSupport;
			
	// loop through all the contigs that were found 
	for (int j = 0; j < this->contigs.size(); j++) {

		vector<PAIR> contigMappings  = contigs[j];			    // get the mappings for this contig
		numMappings                  = contigMappings.size();	// get the number of mappings for this contig
		
		// The "footprint" of the breakpoint
		minStart1 = minStart2 = INT_MAX;	
		maxEnd1   = maxEnd2   = 0;
		
		totalMM1         = totalMM2 = 0;
		numUniqueMappers = numAnchoredMappers = numMultipleMappers = 0;
		totalMappings1   = totalMappings2 = 0;
		totalQual1       = totalQual2 = 0;
		meanMappings1    = meanMappings2 = 0;
		meanQual1        = meanQual2 = 0;
		meanMM1          = meanMM2 = 0.0;
		
		allWeightedSupport   = 0.0;
		finalSupport         = 0;
		finalWeightedSupport = 0.0;
		
		int numAllUniques    = 0;
		int numFinalUniques  = 0;
		int contigSize;  // what is the "span" of the breakpoint
		
		// get the contig size and "footprints" based on the min and max mapping coordinates
		// if there is final support, only use the mappings that were _included_ in the cluster
		// otherwise, define the footprint based on all the mappings that were once in the cluster
		if (hasFinalSupport(contigMappings) == true) {
			contigSize = getClusterSizeFinal(contigMappings, minStart1, minStart2, maxEnd1, maxEnd2);
		}
		else {
			contigSize = getClusterSizeAll(contigMappings, minStart1, minStart2, maxEnd1, maxEnd2);			
		} 

		// (1)  Get _all_ of the unique ids in this cluster
		getNumUniquePairs(contigMappings, numFinalUniques, numAllUniques);
		
		// (2)  Get the total edit distance on each end of this cluster
		getTotalEditDistance(contigMappings, totalMM1, totalMM2);
		
		// (3)  Get the total number of mappings on each end of this cluster
		getTotalNumMappings(contigMappings, totalMappings1, totalMappings2);
		
		// (4)  Tally the final, finalWeighted and allWeightedSupport for this cluster
		computeSupport(contigMappings, finalSupport, finalWeightedSupport, allWeightedSupport, 
						numUniqueMappers, numAnchoredMappers, numMultipleMappers);

		// calculate quality statistics for each contig.
		if (finalSupport) {
			meanMappings1 = totalMappings1 / finalSupport;
			meanMappings2 = totalMappings2 / finalSupport;
			meanQual1     = totalQual1 / finalSupport;
			meanQual2     = totalQual2 / finalSupport; 
			meanMM1       = float(totalMM1) / float(finalSupport);
			meanMM2       = float(totalMM2) / float(finalSupport);
		}

		// report this contigs as longs as there is at least one uniq
		// pair in the contig.
		if (allWeightedSupport > 0) {
			
			contigId++;
						
			// BEDPE-style
			all << contigMappings[0].chrom1 << "\t" << minStart1 << "\t" << maxEnd1 << "\t" <<
				   contigMappings[0].chrom2 << "\t" << minStart2 << "\t" << maxEnd2 << "\t" <<
				   contigId << "\t" << numFinalUniques << "\t" << contigMappings[0].strand1 << "\t" << 
				   contigMappings[0].strand2 << "\t" << meanMM1 << "\t" << meanMM2 << "\t" << 
				   meanMappings1 << "\t" << meanMappings2 << "\t" << contigSize << "\t" << numMappings << "\t" << 
				   allWeightedSupport << "\t" << finalSupport << "\t" << finalWeightedSupport << "\t" << 
				   numUniqueMappers << "\t" << numAnchoredMappers << "\t" << numMultipleMappers << "\t" << endl;
					
			if (numFinalUniques >= this->minSupport) {
				// BEDPE-style
				final << contigMappings[0].chrom1 << "\t" << minStart1 << "\t" << maxEnd1 << "\t" <<
					   contigMappings[0].chrom2 << "\t" << minStart2 << "\t" << maxEnd2 << "\t" <<
					   contigId << "\t" << numFinalUniques << "\t" << contigMappings[0].strand1 << "\t" << 
					   contigMappings[0].strand2 << "\t" << meanMM1 << "\t" << meanMM2 << "\t" << 
					   meanMappings1 << "\t" << meanMappings2 << "\t" << contigSize << "\t" << numMappings << "\t" << 
					   allWeightedSupport << "\t" << finalSupport << "\t" << finalWeightedSupport << "\t" << 
					   numUniqueMappers << "\t" << numAnchoredMappers << "\t" << numMultipleMappers << "\t" << endl;
			}	

			// write the read/mapping information for this contig to the contig detail file.
			for (int i = 0; i < numMappings; i++) {
				if (contigMappings[i].include) {
					detail << contigMappings[i].chrom1 << "\t" << contigMappings[i].start1 << "\t" << contigMappings[i].end1 << "\t"
						   << contigMappings[i].chrom2 << "\t" << contigMappings[i].start2 << "\t" << contigMappings[i].end2 << "\t"
						   << contigMappings[i].readId << "\t" << contigMappings[i].mate1 << "\t" << contigMappings[i].strand1 << "\t"
						   << contigMappings[i].strand2 << "\t" << contigMappings[i].edit1 << "\t" << contigMappings[i].edit2 << "\t"
						   << contigMappings[i].mappings1 << "\t" << contigMappings[i].mappings2 << "\t" 
						   << int(contigMappings[i].mappingType) << "\t" << "Y" << "\t" << contigId << endl;  										
				}
				else {
					detail << contigMappings[i].chrom1 << "\t" << contigMappings[i].start1 << "\t" << contigMappings[i].end1 << "\t"
						   << contigMappings[i].chrom2 << "\t" << contigMappings[i].start2 << "\t" << contigMappings[i].end2 << "\t"
						   << contigMappings[i].readId << "\t" << contigMappings[i].mate1 << "\t" << contigMappings[i].strand1 << "\t"
						   << contigMappings[i].strand2 << "\t" << contigMappings[i].edit1 << "\t" << contigMappings[i].edit2 << "\t"
						   << contigMappings[i].mappings1 << "\t" << contigMappings[i].mappings2 << "\t" 
						   << int(contigMappings[i].mappingType) << "\t" << "N" << "\t" << contigId << endl;
				}
			}
		}
	}
	all.close();
	final.close();
	detail.close();
}



//***********************************************************
// Method: CullMappingsByMisMatches
// Purpose: To reduce the mappings for a given pair
//          to those with the least mismatches.
//
// Outline: 
//	1.  Find the mapping(s) with the least mismatches.
//	2.  Sort the mappings so the least MM set is on "top".	 
//	3.  Ditch all but the least MM set.
//***********************************************************
void HydraPE::CullMappingsByMisMatches(pairVector &pairMappings) {

	// set of unique mappings for each mate in a pair
	std::map<BED3, bool> end1Mappings, end2Mappings;
		
	unsigned short minMM = USHRT_MAX;
	int totalMM;
	
	// loop through the mappings and track the unique 
	// mappings on each end as well as the minimum edit
	// distance observed for any one mapping.  this minimum
	// is the standard for all other mappings.	
	pairVector::const_iterator mapIter = pairMappings.begin();
	pairVector::const_iterator mapEnd = pairMappings.end();
	for (; mapIter != mapEnd; ++mapIter) {		
		if (mapIter->mate1 == 1) {
			BED3 end1Mapping(mapIter->chrom1, mapIter->start1, mapIter->end1);
			BED3 end2Mapping(mapIter->chrom2, mapIter->start2, mapIter->end2);
			end1Mappings[end1Mapping] = true;
			end2Mappings[end2Mapping] = true;			
		}
		else {
			BED3 end2Mapping(mapIter->chrom1, mapIter->start1, mapIter->end1);
			BED3 end1Mapping(mapIter->chrom2, mapIter->start2, mapIter->end2);
			end1Mappings[end1Mapping] = true;
			end2Mappings[end2Mapping] = true;		
		}
		
		// does this mapping have the least mismatches?
		totalMM = getTotalMM(*mapIter);
		if (totalMM < minMM)
			minMM = totalMM;
	}   

	// Figure out what type of pair this is based on the alignments.
	char mapType;
	int end1Size = end1Mappings.size();
	int end2Size = end2Mappings.size();
	
	if ((end1Size == 1) && (end2Size == 1))
		mapType = UNIQ_TYPE;
	else if ((end1Size == 1) || (end2Size == 1))
		mapType = ANCH_TYPE;
	else
		mapType = MULT_TYPE;

	// sort the reads so the mappings with the least mismatches are at the "top"                                                                                                                                                    
	sort(pairMappings.begin(), pairMappings.end(), byTotalMM);
                                                                                                                                        
	// default to assuming we want to use just the best mappings
	// "withinBest" = allow up to editBeyondBest edits worse than minMM
	// "all" = we don't want to delete anything.
	int editDistanceCutoff = minMM;
	if (this->mappingUsage == "withinBest") 
		editDistanceCutoff = minMM + editBeyondBest;
	else if (this->mappingUsage == "all")
		editDistanceCutoff = INT_MAX;
	
	vector<PAIR>::iterator mappingsIter = pairMappings.begin();
	while ((mappingsIter != pairMappings.end()) && (getTotalMM(*mappingsIter) <= editDistanceCutoff)) {
		mappingsIter->mappingType = mapType; // set the mapping type for all of the mappings within the edit distance cutoff
		mappingsIter++;
	}

	// ditch all but the least mm set.
	pairMappings.erase(mappingsIter, pairMappings.end());
	
	// Now, store the number of mappings observed on each end of this pair
	pairVector::iterator pairIter = pairMappings.begin();
	pairVector::iterator pairEnd  = pairMappings.end();
	for (; pairIter != pairEnd; ++pairIter) {
		pairIter->mappings1 = end1Size;
		pairIter->mappings2 = end2Size;
	}		
}
