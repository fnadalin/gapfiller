/*
 ============================================================================
 Name        : GapFiller.cpp
 Author      : Francesca Nadalin, Francesco Vezzi
 Version     : 2.1.1
 Description : A de novo local assembler for paired reads
 ============================================================================
 */

#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
using namespace std;

#include <sys/time.h>
#include <sys/resource.h>

#include "modules/Assembly_only_sequence.h"
#include "modules/Assembly_store_layout.h"

int main(int argc, char * argv[]) {
	clock_t start = clock();

	// Options
	Menu M;
	M.add_options();
	M.parse_options(argc, argv);
	Options Opt;
	Opt.set_values(M);

	// Build Hash table
	Hash * HashHead = new Hash(Opt.k, Opt.blockLength, Opt.overlap);
	HashHead->store_reads(Opt); // store reads into readsMulti
	HashHead->fill(); // compute fingerprints

	// Assembly
	if (Opt.layout) {
		Assembly_store_layout A(Opt, HashHead);
		A.execute(HashHead);
		// A.execute_DEBUG (HashHead);
	} else {
		Assembly_only_sequence A(Opt, HashHead);
		if (Opt.seed_as_query) {
			A.execute(HashHead);
		} else if (Opt.seed_fasta) {
			A.execute_fasta_query(HashHead);
		} else { // Opt.seed_sam
			A.execute_sam_query(HashHead);
		}
	}

	double time2 = clock();
	cout << "Wall time: " << (double) (time2 - start) / CLOCKS_PER_SEC << " s"
			<< endl;

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	cout << "RSS: " << usage.ru_maxrss << endl;

	return 0;
}

