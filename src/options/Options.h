#ifndef OPTIONS_H_
#define OPTIONS_H_

/* DEFAULT OPTIONS */

#define HELP false
#define K 12
#define BLOCK_LENGTH 15
#define LAYOUT false
#define OUTPUT_PREFIX "GapFiller_output"
#define GZIP false
#define BZIP2 false
#define SEED_FASTA false
#define SEED_SAM false
#define QUERY_FASTA false
#define QUERY_SAM false
#define SEED_AS_QUERY true
#define OVERLAP 30
#define GLOBAL_MISMATCH 5
#define EXT_THRESHOLD 2
#define LIMITED_EXTENSION false
#define NO_CYCLE false
#define REVERSE_SEED false
#define REVERSE_MATE true
#define VERBOSE false

/* INCLUDES */

#include <vector>
using namespace std;

#include <common.h>
#include "options/Menu.h"

class Options {

public:

	Options (); // set to default values
	void set_values (Menu & values);

	bool 			help;
	string 			outputPrefix;
	bool			gzip;
	bool			bzip2;
	unsigned int 	k;
	unsigned int	blockLength;
	bool 			layout;

	vector<string> 	seed_fw_files; // conditionally required
	vector<string> 	seed_rv_files; // conditionally required
	vector<string> 	seed_sam_files; // conditionally required

	unsigned int 	seed_ins; // required
	unsigned int	seed_var; // required
	bool 			seed_fasta;
	bool			seed_sam;

	vector<string> 	query_files;
	vector<string> 	query_sam_files;

	bool 			query_fasta;
	bool 			query_sam;
	bool 			seed_as_query;

	unsigned int 	overlap;
	unsigned int 	globalMismatch;
	unsigned int 	extThreshold;
	bool 			limitedExtension;
	unsigned int	maxReadsToExtend;
	bool 			no_cycle;
	bool 			reverse_seed;
	bool 			reverse_mate;
	bool 			verbose;

};

#endif /* OPTIONS_H_ */
