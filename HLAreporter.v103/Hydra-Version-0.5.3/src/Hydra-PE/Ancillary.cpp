/********************************************************************************
Program: 		Ancillary.cpp
Description: 	Ancillary helper functions for Hydra data structures.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*********************************************************************************/
#include "Ancillary.h"


/******************************************************************************
   Utility functions.
*******************************************************************************/
int getFragSize (const PAIR &pair) {
	/* old, Contigger */
	//return (a.start2 - a.start1 + 1);
	return (pair.start2 - pair.end1 + 1);	
}


bool doSpansSupportCommonBreakpoint (const PAIR &a, const PAIR &b) {
	/***************************************************************
	 We want to enforce that two mappings have acceptable non-overlap
	 AND that they support the same putative breakpoint region.
	
	For example, we want this:
	
	     s1   e1                     s2    e2
	(A)   ======......................======
	(B)                   ======...............................=====
	                      s1   e1                             s2   e2
	
	But NOT this:
	     s1   e1                     s2    e2
	(A)   ======......................======
	(B)                          ======...............................=====
	                             s1   e1                             s2   e2	
	****************************************************************/
	
	if ((a.chrom1 == a.chrom2) && (b.chrom1 == b.chrom2)) {
		int overlap = (min(a.end2, b.end2) - max(a.start1, b.start1));
		if ( (overlap > 0) && 
		     (a.end1 < b.start2) && 
		     (b.end1 < a.start2))
			return true;
		else
			return false;
	}
	else return true;  // can't enforce this for inter-chromosomals, so just return true
}


int getNonOverlap (const PAIR &a, const PAIR &b) {
	/* old, Contigger */
	// return ( abs(b.start1 - a.start1) + abs(b.start2 - a.start2) );
	return ( abs(b.start1 - a.start1) + abs(b.end2 - a.end2) );
}


bool doLengthsSupportOneAnother (const PAIR &a, const PAIR &b, int lengthDev) {
	return abs( getFragSize(a) - getFragSize(b) ) <= lengthDev;
}


bool doSpansSupportOneAnother (const PAIR &a, const PAIR &b, int spanDev) {
	/********************************************************************************
	  Two mapping spans support one another so long as there is at least 1bp of overlap
	  __and__ the non-overlap between the spans is within the tolerance of the 
	  library(ies).  The latter constraint prevents situations where you have one
	  base pair of overlap but the mappings are say, 100kb.  Exampli gratis:
	 
	  +...................................-
	                                    +.........................................-
	  
	   [<-------------100kb----------------]
	                                    [<-------------100kb----------------------]
	
	   Yet, the two mappings come from a library where the median + 10(mad) is 2000bp.
	   In such a case, the two mappings clearly do NOT support the same breakpoint.
	***********************************************************************************/
	if ( (doSpansSupportCommonBreakpoint(a, b) == true) && (getNonOverlap(a, b) <= spanDev) ) {
		return true;
	}
	return false;
}


int getTotalMM (const PAIR &pair) {
	return (pair.edit1 + pair.edit2);
}


int getTotalMMAmongAllMappings(const pairVector &mappings) {
	
	int totalMM = 0;
	
	// compute the total number of mismatches among all of the mappings
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) totalMM = getTotalMM(*mapIter);

	return totalMM;
}


double getTotalWeightedSupportAmongAllMappings(const pairVector &mappings) {
	double totalWeightedSupport = 0.0;
	
	// compute the total weighted support among all of the mappings
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();	
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->mappingType == UNIQ_TYPE) {
			totalWeightedSupport += UNIQ_WEIGHT;
		}
		else if (mapIter->mappingType == ANCH_TYPE) {
			totalWeightedSupport += ANCH_WEIGHT;
		}
		else if (mapIter->mappingType == MULT_TYPE) {
			totalWeightedSupport += MULT_WEIGHT;
		}
	}
	return totalWeightedSupport;
}


// Return true if any of the mappings have include == true.
// That is, there is at least one mapping left in this putative cluster.
bool hasFinalSupport(const pairVector &mappings) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->include == true) return true;
	}
	return false;
}


int getClusterSizeAll(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->start1 < minStart1)    minStart1 = mapIter->start1;
		if (mapIter->end1 > maxEnd1)          maxEnd1 = mapIter->end1;
		if (mapIter->start2 < minStart2)    minStart2 = mapIter->start2;
		if (mapIter->end2 > maxEnd2)          maxEnd2 = mapIter->end2;
	}
	return (maxEnd2 - minStart1) + 1;
}


int getClusterSizeFinal(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		// only compute add to footprint if included in final call.
		if (mapIter->include == true) {
			if (mapIter->start1 < minStart1)    minStart1 = mapIter->start1;
			if (mapIter->end1 > maxEnd1)          maxEnd1 = mapIter->end1;
			if (mapIter->start2 < minStart2)    minStart2 = mapIter->start2;
			if (mapIter->end2 > maxEnd2)          maxEnd2 = mapIter->end2;
		}
	}
	return (maxEnd2 - minStart1) + 1;
}


bool isVariantUnlinked(const pairVector &mappings, int maxDistance) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if ( mapIter->chrom1 != mapIter->chrom2 ) return true;
		else if ( (mapIter->end2 - mapIter->start1) > maxDistance ) return true;
	}
	return false;
}


bool areAnyContigsUnlinked(const vector<IN_CONTIG> &contigs) {
	vector<IN_CONTIG>::const_iterator contigIter = contigs.begin();
	vector<IN_CONTIG>::const_iterator contigEnd  = contigs.end();
	for (; contigIter != contigEnd; ++contigIter) {
		if ( contigIter->unlinked ) return true;
	}
	return false;
}


void getNumUniquePairs(const pairVector &mappings, int &numFinalUniques, int &numAllUniques) {

	// maps to track the unique pair ids that are clustered.
	std::map<string, short, std::less<string> > allUniqIds;
	std::map<string, short, std::less<string> > finalUniqIds;
	
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->include == true) finalUniqIds[mapIter->readId]++;
		allUniqIds[mapIter->readId]++;		
	}
	numFinalUniques = finalUniqIds.size();
	numAllUniques   = allUniqIds.size();
}


void getTotalEditDistance(const pairVector &mappings, int &totalMM1, int &totalMM2) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->include == true) {
			totalMM1 += mapIter->edit1;
			totalMM2 += mapIter->edit2;
		}
	}
}


void getTotalNumMappings(const pairVector &mappings, int &totalMappings1, int &totalMappings2) {
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		if (mapIter->include == true) {
			if (mapIter->mate1 == 1) {
				totalMappings1 += mapIter->mappings1;
				totalMappings2 += mapIter->mappings2;
			}
			else {
				totalMappings1 += mapIter->mappings2;
				totalMappings2 += mapIter->mappings1;				
			}
		}
	}
}


void computeSupport(const pairVector &mappings, int &finalSupport, double &finalWeightedSupport, double &allWeightedSupport, 
				    int &numUniqueMappers, int &numAnchoredMappers, int &numMultipleMappers) {
					
	pairVector::const_iterator mapIter = mappings.begin();
	pairVector::const_iterator mapEnd  = mappings.end();
	for (; mapIter != mapEnd; ++mapIter) {
		
		if (mapIter->mappingType == UNIQ_TYPE) {		
			if (mapIter->include) {
				finalSupport++;
				finalWeightedSupport += UNIQ_WEIGHT;
				numUniqueMappers++;
			}
			allWeightedSupport += UNIQ_WEIGHT;
		}
		else if (mapIter->mappingType == ANCH_TYPE) {
			if (mapIter->include) {
				finalSupport++;
				finalWeightedSupport += ANCH_WEIGHT;
				numAnchoredMappers++;
			}
			allWeightedSupport += ANCH_WEIGHT;
		}
		else if (mapIter->mappingType == MULT_TYPE) {
			if (mapIter->include) {
				finalSupport++;
				finalWeightedSupport += MULT_WEIGHT;
				numMultipleMappers++;
			}
			allWeightedSupport += MULT_WEIGHT;
		}
	}				
}
