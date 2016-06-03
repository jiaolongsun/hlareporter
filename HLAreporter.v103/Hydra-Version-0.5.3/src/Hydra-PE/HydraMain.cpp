#include <iostream>
#include <ctime>
#include "Hydra.h"
#include "Version.h"

using namespace std;

// define our program name
const string PROGRAM_NAME  = "hydra";

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
bool IsNumber(const std::string& s);

int main(int argc, char* argv[]) {

	// initialization
	ostringstream sb;
	char* end_ptr = NULL;

	// our configuration variables
	bool showHelp          = false;

	int lengthDev         = 0; 
	int spanDev           = 0;
	int minSupport        = 2;
	int maxLinkedDistance = 1000000; // 1Mb by default
	
	
	// input files
	string discordantsFile;

	// output files
	string outFile;

	// checks for existence of parameters
	bool haveDiscordants = false;
	bool haveOutFile = false;
	bool haveLD = false;
	bool haveSD = false;
	bool haveMinSupport = false;
	bool haveMaxLinkedDistance = false;
	bool lumpInversions = false;
	bool ignoreSize = false;
	string mappingUsage = "best";
	unsigned int editBeyondBest = 0;
	

	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-h", 2, parameterLength) || 
		PARAMETER_CHECK("--help", 5, parameterLength)) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();

	// do some parsing (all of these parameters require 2 strings)
	cerr << endl << "Parameters:" << endl;
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-mld", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveLD = true;
				lengthDev = atoi(argv[i + 1]);
				cerr << "  Using maximum length difference (-mld) of: " << lengthDev << endl; 
				i++;
			}
		} 
		else if(PARAMETER_CHECK("-mno", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveSD = true;
				spanDev = atoi(argv[i + 1]);
				cerr << "  Using maximum non-overlap (-mno) of: " << spanDev << endl;
				i++;
			}
		} 
		else if(PARAMETER_CHECK("-in", 3, parameterLength)) {
			if ((i+1) < argc) {
				haveDiscordants = true;
				discordantsFile = argv[i + 1];
				cerr << "  Discordant mappings input file (-in): " << discordantsFile << endl;
				i++;
			}
		} 
		else if(PARAMETER_CHECK("-out", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveOutFile = true;
				outFile = argv[i + 1];
				cerr << "  Summary breakpoint output file (-out): " << outFile << endl;
				cerr << "  Detailed breakpoint output file: " << outFile << ".detail" << endl;
				i++;
			}
		}
		else if(PARAMETER_CHECK("-ms", 3, parameterLength)) {
			if ((i+1) < argc) {
				haveMinSupport = true;
				minSupport = atoi(argv[i + 1]);
				cerr << "  Minimum number of supporting pairs (-ms): " << minSupport << endl;
				i++;
			}
		}
		else if(PARAMETER_CHECK("-lnk", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveMaxLinkedDistance = true;
				maxLinkedDistance = atoi(argv[i + 1]);
				cerr << "  Maximum intrachromosomal distance before unlinked (-maxIntra): " << maxLinkedDistance << endl;
				i++;
			}
		}
		else if(PARAMETER_CHECK("-is", 3, parameterLength)) {
			ignoreSize = true;
			cerr << "  Break cluster ties based on edit distance instead of size. " << endl;
		}
		else if(PARAMETER_CHECK("-li", 3, parameterLength)) {
			lumpInversions = true;
			cerr << "  Combine +/+ and -/- mappings when finding inversions. " << endl;
		}
		else if(PARAMETER_CHECK("-use", 4, parameterLength)) {
			
			if ((i+1) < argc) {
				string useValue = argv[i + 1];
				if (useValue == "best") {
					mappingUsage = "best";
					editBeyondBest = 0;
					cerr << "  Using *best* mappings." << endl;
				}
				else if (useValue == "all") {
					mappingUsage = "all";
					editBeyondBest = INT_MAX;
					cerr << "  Using *all* mappings." << endl;
				}
				else if ( IsNumber(useValue) ) {
					mappingUsage = "withinBest";
					editBeyondBest = atoi(useValue.c_str());
					cerr << "  Using all mappings within " << editBeyondBest << "edits of best." << endl;					
				}
				i++;
			}
		}
		else {
			cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	if (!haveDiscordants) {
		cerr << "*****ERROR: You must specify an input mapping file.*****" << endl << endl;
		showHelp = true;
	}
	if (!haveOutFile) {
		cerr << "*****ERROR: You must specify an output file.*****" << endl << endl;
		showHelp = true;
	}

	if (!showHelp) {

		cerr << endl << "Processing: " << endl;;
		HydraPE *events = new HydraPE(lengthDev, spanDev, minSupport,
												maxLinkedDistance, ignoreSize, lumpInversions,
												mappingUsage, editBeyondBest);

		if (haveDiscordants) {
			cerr << "  Loading discordant mappings." << endl;  
			events->ReadDiscordantFile(discordantsFile);
		} 
		events->SortFragments();

		// Assemble logical contigs from the mappings
		cerr << "  Refining breakpoint clusters." << endl;
		events->AssembleContigs();

		// Ensure that each read is only in one contig
		cerr << "  Resolving mappings belonging to multiple breakpoint clusters." << endl;
		events->AllowOneContigPerRead();

		// Summarize the contigs
		cerr << "  Writing breakpoint files." << endl;
		events->ReportSVCalls(outFile);

		cerr << "\n" << PROGRAM_NAME << " has completed." << endl;
		return 0;
	}
	else {
		ShowHelp();
	}
}


void ShowHelp(void) {
	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Calls SV breakpoints from discordant paired-end mappings." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -in <mappings> -out <breakpoints> -mld <bp> -mno <bp>" << endl << endl;

	cerr << "Options:" << endl;
	cerr << "  -in\tFile of discordant mappings." << endl;

	cerr << "  -out\tStub for output files." << endl;

	cerr << "  -mld\tMaximum allowable length difference b/w mappings." << endl;
	cerr << "  \tTypically set to 10 * m.a.d. of the DNA fragment libraries." << endl;
	cerr << "  \tsee: http://en.wikipedia.org/wiki/Median_absolute_deviation" << endl << endl;
	
	cerr << "  -mno\tMaximum allowable non-overlap b/w mappings." << endl;
	cerr << "  \tTypically set to median + (20 * m.a.d.) of the DNA fragment libraries." << endl << endl;;
	
	cerr << "  -ms\tMinimum number of pairs required for variant to be called." << endl;
	cerr << "\tDefault: 2" << endl << endl;

	cerr << "  -lnk\tMaximum intrachromosomal distance allowed before a" << endl;
	cerr << "  \tvariant is considered to be between unlinked DNA segments." << endl;
	cerr << "  \tDefault: 1000000 (i.e., 1Mb)" << endl << endl;
	
	cerr << "  -is\tChoose most likely variant (when a tie exists) based on" << endl;
	cerr << "  \tleast edit distance rather than size." << endl << endl;

	cerr << "  -li\tCombine +/+ and -/- mappings when screening for inversions." << endl;
	cerr << "  \tThis increases sensitivity in low coverage." << endl << endl;

	cerr << "  -use\tWhich mappings should be used for each pair?" << endl;
	cerr << "  \t\"best\"\tUse the mappings with the least edit distance for each pair." << endl;
	cerr << "  \t\tDefault." << endl;
	cerr << "  \t\"all\"\tUse all mappings for each pair." << endl;
	cerr << "  \t<INT>\tUse the best plus those within <INT> edit distance of best." << endl;
	
	// end the program here
	exit(1);
}


bool IsNumber(const std::string& s) {
   for (int i = 0; i < s.length(); i++) {
       if (!std::isdigit(s[i]))
           return false;
   }
   return true;
}

