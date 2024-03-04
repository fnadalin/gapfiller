#ifndef HASH_H_
#define HASH_H_

#define _HASH_VERSION 0

#include <iostream>
#include <vector>
using namespace std;

#include "errors/Generic_Exception.h"
#include "errors/Incorrect_Format.h"
using namespace errors;

#include "options/Options.h"
#include "io/Auto_Unzip.h"
#include "io/Fastq.h"
#include "io/Sam.h"
using namespace fastq;

#include "data_structures/include.h"
#include "data_structures/Reads.h"

// typedef pair <unsigned int, bool> readOriented;

void process_mem_usage(double& vm_usage, double& resident_set);

class Hash {

public:

	Hash (unsigned int k_hash, unsigned int b, unsigned int overlap);
	~Hash();

	void					store_reads (const Options & Opt);
	void 					fill ();
	void					delete_HASHcounter ();

	unsigned long int 		ComputeFingerprint (unsigned int pos, bool orientation, unsigned int start);

	unsigned int 			k;
	unsigned long int 		q;
	unsigned int 			blockLength;
	unsigned int			overlap;
	unsigned long int 		h;
	unsigned int 			numReads;
	unsigned int			maxReadLength;

	vector <readOriented> 	HASHvalues;
	vector <Reads>			readsMulti;
	unsigned int * 			HASHcounter;

protected:

	size_t 					limit;

	void					read_fastq (const vector <string> &, const vector <string> &);
	void					read_fastq (const vector <string> & filename);
	void 					read_sam (const vector <string> & filename) throw (Data_Exception);

};


#endif /* HASH_H_ */
